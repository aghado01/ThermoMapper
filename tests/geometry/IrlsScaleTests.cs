using System;
using System.Linq;
using Maths.Geometry.Estimators.Ambient;
using Maths.Geometry.Estimators.Intrinsic;
using Maths.Geometry.Solver;
using Maths.LinAlg;
using Xunit;

namespace Maths.Geometry.Tests;

/// <summary>
/// The IRLS solver's spread unit (<see cref="IrlsScale.Spread"/>, the default): the median and mean of c + s·u equal
/// c + s·(the estimate for u), across spreads 1e-12 to 1e6 on flat space and down to the manifold's resolution on the
/// sphere and the Grassmannian. Truth comes from outside the solver: a plain weighted Weiszfeld at scale 1 for the flat
/// median, and the tangent-space construction exp_{p₀}(s·m_u), exact to O(s³), on curved manifolds.
/// </summary>
public sealed class IrlsScaleTests
{
    private const int MaxIterations = 200;   // IrlsOptions.Default
    private const double MachineEpsilon = 2.220446049250313e-16;
    private static readonly double[] Spreads = Enumerable.Range(-12, 19).Select(k => Math.Pow(10, k)).ToArray();

    /// <summary>
    /// A seeded sweep over dimension, n, equal and unequal weights and an occasional far outlier, centred at the
    /// origin and at 100·s from it. The median (warm-started at the medoid, as the callers do, and at the weighted
    /// mean), the mean and <see cref="EuclideanMedian"/> (warm-started at the mean, as GMM does) keep a relative
    /// accuracy that does not depend on s, and no s costs the median more iterations than s = 1. Most draws converge in
    /// 10 to 70 iterations; one (d = 3, n = 5) has a data point that barely fails the optimality condition
    /// (‖G‖/w = 1.07), so Weiszfeld contracts at about 0.93 and reaches MaxIterations at every s alike, within the
    /// accuracy. Some draws have a data point for their median, which the optimality check returns.
    /// </summary>
    [Fact]
    public void FlatMedianAndMean_AreScaleEquivariant_FromSpread1e12To1e6()
    {
        var rng = new Random(20261108);
        for (int draw = 0; draw < 12; draw++)
        {
            FlatProblem problem = FlatProblem.Draw(rng);
            int[] reference = CheckFlat(problem, new double[problem.D], 1.0, roundingUlps: 0, $"draw {draw}, scale 1");
            foreach (bool covariant in new[] { false, true })
                foreach (double s in Spreads)
                {
                    double[] c = covariant ? Scaled(problem.Direction, 100.0 * s) : new double[problem.D];
                    CheckFlat(problem, c, s, roundingUlps: 0, $"draw {draw}, centre {(covariant ? "100s" : "0")}, s {s:G3}", reference);
                }
        }
    }

    /// <summary>
    /// Seven collinear points in R³ with unequal weights: the median is the weighted 1-D median, a data point. From the
    /// mean and from the leftmost point the solver returns that point exactly, at every spread, in a few iterations.
    /// </summary>
    [Fact]
    public void CollinearMedian_IsTheDataPoint_ReturnedByTheOptimalityCheck()
    {
        double[] t = [-3, -1, 0, 0.5, 2.5, 4, 7];
        double[] w = [1.0, 0.5, 2.0, 1.5, 0.7, 1.2, 0.9];
        const int medianIndex = 3;   // cumulative weight reaches half the total (3.9) at t = 0.5
        double[] direction = [1.0 / 3, 2.0 / 3, -2.0 / 3];
        var euclid = new EuclideanVectorManifold(3);
        foreach (double s in Spreads)
        {
            double[] c = [3 * s, -s, 2 * s];
            double[][] x = t.Select(ti => c.Select((cj, j) => cj + s * ti * direction[j]).ToArray()).ToArray();
            foreach (double[] warm in new[] { WeightedMean(x, w), x[0] })
            {
                var counter = new Counter();
                double[] median = (double[])warm.Clone();
                GeometricMedian.Compute(new Counting<EuclideanVectorManifold>(euclid, counter), x, w, median);
                Assert.Equal(x[medianIndex], median);
                Assert.True(counter.Steps < 30, $"s {s:G3}: {counter.Steps} iterations");
            }
        }
    }

    /// <summary>
    /// A point given twice with weight 1 each, and two more points: the duplicated point holds half the weight and is
    /// the median, though each copy alone fails the optimality condition (‖G‖ = √2 > 1). The check sums the weight
    /// sitting on the point, so the solver returns it, from the mean and from the medoid, at every spread.
    /// </summary>
    [Fact]
    public void DuplicatedMedian_HoldingHalfTheWeight_IsReturned()
    {
        var euclid = new EuclideanVectorManifold(2);
        double[] w = [1.0, 1.0, 1.0, 1.0];
        foreach (double s in Spreads)
        {
            double[] a = [5 * s, -2 * s];
            double[][] x = [a, (double[])a.Clone(), [a[0] + s, a[1]], [a[0], a[1] + s]];
            foreach (double[] warm in new[] { WeightedMean(x, w), x[GeometricMedian.MedoidIndex(euclid, x)] })
            {
                var counter = new Counter();
                double[] median = (double[])warm.Clone();
                GeometricMedian.Compute(new Counting<EuclideanVectorManifold>(euclid, counter), x, w, median);
                Assert.Equal(a, median);
                Assert.True(counter.Steps < 30, $"s {s:G3}: {counter.Steps} iterations");
            }
        }
    }

    /// <summary>
    /// A seeded sweep over dimension, n, equal and unequal weights and an occasional far outlier, centred 10 from the
    /// origin at every s and warm-started at the medoid and at the mean. The data are then representable only to
    /// about ε·‖c‖: the median, the mean and <see cref="EuclideanMedian"/> (warm-started at the mean, as GMM does) stay
    /// within a fixed number of those ulps of the truth for the represented data, and no s costs the median more
    /// iterations than s = 1 at the origin.
    /// </summary>
    [Fact]
    public void FlatMedianAndMean_FarFromTheOrigin_AreLimitedOnlyByTheCentresRounding()
    {
        var rng = new Random(81102026);
        for (int draw = 0; draw < 6; draw++)
        {
            FlatProblem problem = FlatProblem.Draw(rng);
            int[] reference = CheckFlat(problem, new double[problem.D], 1.0, roundingUlps: 0, $"draw {draw}, scale 1");
            double[] c = Scaled(problem.Direction, 10.0);
            foreach (double s in Spreads)
                CheckFlat(problem, c, s, roundingUlps: 1000, $"draw {draw}, centre 10, s {s:G3}", reference);
        }
    }

    /// <summary>
    /// On S⁹ and the Grassmannians Gr(2, 10) and Gr(3, 20), 21 points exp_{p₀}(s·uᵢ) with spreads from 1e-1 down to
    /// 1e-12 on the Grassmannians, whose angles are atan2-based, and 1e-6 on the sphere, whose are still arccos-based.
    /// The allowance holds the truth's O(s²) relative error, the 1e-6 accuracy of the flat case, and floor/s: 1e-14/s on
    /// the Grassmannians for the rounding of the unit-length points themselves (measured about 5e-16/s), 1e-10/s on the
    /// sphere, whose Exp drops tangent steps below 1e-12. No median takes more than 50 iterations.
    /// </summary>
    [Theory]
    [InlineData(10, 0, 6, 1e-10, 7001)]    // S⁹
    [InlineData(10, 2, 12, 1e-14, 7002)]   // Gr(2, 10)
    [InlineData(20, 3, 12, 1e-14, 7003)]   // Gr(3, 20)
    public void CurvedMedianAndMean_AreScaleEquivariant_DownToTheManifoldResolution(int n, int r, int deepest, double floor, int seed)
    {
        var rng = new Random(seed);
        if (r == 0) CheckCurved(new SphericalManifold(n), n, 1, sphere: true, deepest, floor, rng);
        else CheckCurved(new GrassmannManifold(n, r), n, r, sphere: false, deepest, floor, rng);
    }

    /// <summary>
    /// Frames on the cut locus: xz meets the three xy-like frames at a principal angle of exactly π/2, where the
    /// Grassmann Log map drops that direction and returns a vector short of the distance (here zero). Read as the
    /// distance, xz's radius is π/2, not 0: it is not coincident with xy, and the optimality check does not certify
    /// xz. xy, which holds half the weight with its rotated copy, is the median and is returned exactly from the
    /// medoid. (From xz itself every Log is zero and the iteration cannot leave it, which is why callers warm-start
    /// at the medoid.)
    /// </summary>
    [Fact]
    public void GrassmannMedian_WithFramesOnTheCutLocus_ReadsTheirRadiusAsTheDistance()
    {
        double[] xz = [1, 0, 0, 0, 0, 1];          // 3×2, column-major
        double[] xy = [1, 0, 0, 0, 1, 0];
        double[] xyRotated = [0, 1, 0, -1, 0, 0];
        double[] xyTilted = [Math.Cos(0.04), 0, Math.Sin(0.04), 0, 1, 0];
        double[][] frames = [xz, xy, xyRotated, xyTilted];
        var grass = new GrassmannManifold(3, 2);

        double[] median = (double[])frames[GeometricMedian.MedoidIndex(grass, frames)].Clone();
        GeometricMedian.Compute(grass, frames, [1.0, 1.0, 1.0, 1.0], median);
        Assert.Equal(xy, median);
    }

    /// <summary>
    /// When every point sits on the warm start the spread is 0: both medians return that point at once, and
    /// <see cref="EuclideanMedian"/>'s scatter is zero rather than the 0/0 of an unguarded floor ε·σ.
    /// </summary>
    [Fact]
    public void IdenticalPoints_GiveThatPoint_AndAZeroScatter()
    {
        double[][] x = Enumerable.Repeat(new[] { 1e-9, -3e-9, 2e-9 }, 5).ToArray();
        double[] w = [1.0, 2.0, 0.5, 1.0, 1.5];

        var counter = new Counter();
        double[] median = (double[])x[0].Clone();
        GeometricMedian.Compute(new Counting<EuclideanVectorManifold>(new EuclideanVectorManifold(3), counter), x, w, median);
        Assert.Equal(x[0], median);
        Assert.Equal(0, counter.Steps);

        double[] ambient = (double[])x[0].Clone();
        double[] scatter = new double[9];
        EuclideanMedian.ComputeWithScatter(x, Enumerable.Range(0, 5).ToArray(), w, ambient, scatter);
        Assert.Equal(x[0], ambient);
        Assert.All(scatter, v => Assert.Equal(0.0, v));
    }

    /// <summary>
    /// <see cref="IrlsScale.Data"/> reads the constants in the data's units and reproduces the solver as it was before
    /// the spread unit (afc814d) bit for bit on a Euclidean median and on <see cref="EuclideanMedian"/> with its scatter.
    /// The Grassmann median and mean are pinned at the Grassmannian's atan2 angles: they moved by up to 7e-11 from
    /// afc814d's arccos angles, under the same 1e-8 absolute tolerance. One more Grassmann median, of seven frames, has a
    /// data point for its median: frame 0, where ‖G‖ ≈ 0.90 ≤ 1. The optimality check at the nearest data point returns
    /// it exactly; before that check the solver stopped about 2e-4 from it.
    /// </summary>
    [Fact]
    [Trait("Tier", "Pinned")]
    public void DataScale_ReproducesTheSolverInDataUnits_BitForBit()
    {
        IrlsOptions data = IrlsOptions.Default with { Scale = IrlsScale.Data };

        double[][] e =
        [
            [0.3, -1.2, 2.0], [1.1, 0.4, -0.7], [-0.8, 0.9, 0.5], [2.4, -0.3, 1.1], [0.0, 1.7, -1.4], [-1.5, -0.6, 0.2],
        ];
        double[] we = [1.0, 0.5, 2.0, 1.5, 0.8, 1.2];
        double[] mean = new double[3];
        for (int i = 0; i < 6; i++) for (int j = 0; j < 3; j++) mean[j] += e[i][j] / 6;

        double[] median = (double[])mean.Clone();
        GeometricMedian.Compute(new EuclideanVectorManifold(3), e, we, median, data);
        Assert.Equal(new[] { -0.3771572745186732, 0.4475843601669949, 0.4489833045052648 }, median);

        double[][] offsets =
        [
            [0.0, 0, 0, 0, 0, 0], [0.05, -0.02, 0.03, -0.04, 0.02, 0.01], [-0.04, 0.03, -0.02, 0.03, -0.01, 0.02],
            [0.03, 0.02, -0.04, 0.01, 0.04, -0.02], [-0.02, -0.04, 0.02, 0.04, -0.03, 0.03],
            [0.08, -0.05, 0.05, -0.08, 0.03, 0.04], [-0.07, 0.04, -0.06, 0.06, -0.04, -0.03],
        ];
        double[][] frames = offsets.Select(o =>
        {
            var f = new double[10];
            f[0] = 1; f[6] = 1;
            int t = 0;
            for (int col = 0; col < 2; col++) for (int row = 2; row < 5; row++) f[col * 5 + row] = o[t++];
            MatrixOps.Orthonormalize(f, 5, 2);
            return f;
        }).ToArray();
        double[] init = new double[10];
        foreach (double[] f in frames) for (int i = 0; i < 10; i++) init[i] += f[i] / frames.Length;
        MatrixOps.Orthonormalize(init, 5, 2);
        double[] gw = Enumerable.Repeat(1.0, frames.Length).ToArray();

        double[] gMedian = (double[])init.Clone();
        GeometricMedian.Compute(new GrassmannManifold(5, 2), frames, gw, gMedian, data);
        Assert.Equal(frames[0], gMedian);

        double[][] inner = frames.Skip(1).ToArray();
        double[] innerInit = new double[10];
        foreach (double[] f in inner) for (int i = 0; i < 10; i++) innerInit[i] += f[i] / inner.Length;
        MatrixOps.Orthonormalize(innerInit, 5, 2);
        double[] innerMedian = (double[])innerInit.Clone();
        GeometricMedian.Compute(
            new GrassmannManifold(5, 2), inner, Enumerable.Repeat(1.0, inner.Length).ToArray(), innerMedian, data);
        Assert.Equal(
            new[]
            {
                0.9999826582733471, 4.4339852850683335E-05, 0.0008525723808938266, -0.00030146060784903004,
                -0.005819229191660574, -1.9387987999012176E-06, 0.9999190835745473, 0.009274690978703464,
                0.0015753784497430084, 0.008562977850508,
            },
            innerMedian);

        double[] gMean = (double[])init.Clone();
        GeometricMean.Compute(new GrassmannManifold(5, 2), frames, gw, gMean, data);
        Assert.Equal(
            new[]
            {
                0.9999827000029012, 2.3985892716460462E-08, 0.004271230044846536, -0.0028580817238377633,
                -0.002861408336970325, 1.2148160770930532E-05, 0.999969224670312, 0.00289135840541956,
                0.001417337975397572, 0.007154073267174528,
            },
            gMean);

        double[] ambient = (double[])mean.Clone();
        double[] scatter = new double[9];
        EuclideanMedian.ComputeWithScatter(e, Enumerable.Range(0, 6).ToArray(), we, ambient, scatter, scale: IrlsScale.Data);
        Assert.Equal(new[] { -0.3771572680082368, 0.44758435399766566, 0.44898330453321866 }, ambient);
        Assert.Equal(
            new[]
            {
                1.1363906519075286, -0.1964486633971708, 0.14631067513452747, -0.1964486633971708, 0.6303202155591183,
                -0.3351119768184708, 0.14631067513452747, -0.3351119768184708, 0.5142937573518758,
            },
            scatter);
    }

    // ── Flat checks ──────────────────────────────────────────────────────────

    // Allowances are in the u-units of the truth: 1e-6 of the data's spread for the medians, 1e-12 for the mean, plus
    // roundingUlps·ε·‖c‖/s for a centre far from the origin. Returns the medians' iteration counts (medoid and mean
    // warm starts); with stepBound, each is at most its bound plus 2.
    private static int[] CheckFlat(FlatProblem problem, double[] c, double s, int roundingUlps, string label, int[]? stepBound = null)
    {
        int d = problem.D, n = problem.U.Length;
        double[][] x = problem.U.Select(u => u.Select((v, j) => c[j] + s * v).ToArray()).ToArray();
        // u' = (x − c)/s: the data as represented, which the truth must use when c is far from the origin.
        double[][] represented = x.Select(xi => xi.Select((v, j) => (v - c[j]) / s).ToArray()).ToArray();
        double[] medianTruth = WeightedWeiszfeld(represented, problem.W);
        double[] meanTruth = WeightedMean(represented, problem.W);
        double spread = MedianDistance(represented, medianTruth);
        double rounding = roundingUlps * MachineEpsilon * Norm(c) / s;

        var euclid = new EuclideanVectorManifold(d);
        double[] warm = x[GeometricMedian.MedoidIndex(euclid, x)];

        double[][] starts = [warm, WeightedMean(x, problem.W)];
        var steps = new int[starts.Length];
        for (int k = 0; k < starts.Length; k++)
        {
            var counter = new Counter();
            double[] median = (double[])starts[k].Clone();
            GeometricMedian.Compute(new Counting<EuclideanVectorManifold>(euclid, counter), x, problem.W, median);
            steps[k] = counter.Steps;
            if (stepBound is not null)
                Assert.True(steps[k] <= stepBound[k] + 2, $"{label}: {steps[k]} iterations against {stepBound[k]} at scale 1");
            double medianError = Gap(median, c, s, medianTruth);
            Assert.True(medianError <= 1e-6 * spread + rounding, $"{label}: median error {medianError:G3} (spread {spread:G3})");
        }

        double[] mean = (double[])warm.Clone();
        GeometricMean.Compute(euclid, x, problem.W, mean);
        double meanError = Gap(mean, c, s, meanTruth);
        Assert.True(meanError <= 1e-12 * spread + rounding, $"{label}: mean error {meanError:G3}");

        double[] ambient = WeightedMean(x, problem.W);
        EuclideanMedian.Compute(x, Enumerable.Range(0, n).ToArray(), problem.W, ambient);
        double ambientError = Gap(ambient, c, s, medianTruth);
        Assert.True(ambientError <= 1e-6 * spread + rounding, $"{label}: EuclideanMedian error {ambientError:G3}");
        return steps;
    }

    private sealed record FlatProblem(int D, double[][] U, double[] W, double[] Direction)
    {
        public static FlatProblem Draw(Random rng)
        {
            int[] dims = [2, 3, 5, 8, 20];
            int d = dims[rng.Next(dims.Length)];
            int n = 5 + rng.Next(36);
            var u = new double[n][];
            for (int i = 0; i < n; i++) u[i] = Gaussian(rng, d, 1.0);
            if (rng.NextDouble() < 0.4) u[n - 1] = Scaled(u[n - 1], 50.0);   // a far outlier
            bool equal = rng.NextDouble() < 0.5;
            var w = new double[n];
            for (int i = 0; i < n; i++) w[i] = equal ? 1.0 : 0.2 + 1.8 * rng.NextDouble();
            double[] direction = Gaussian(rng, d, 1.0);
            return new FlatProblem(d, u, w, Scaled(direction, 1.0 / Norm(direction)));
        }
    }

    // ── Curved checks ────────────────────────────────────────────────────────

    private static void CheckCurved<TM>(TM manifold, int n, int r, bool sphere, int deepest, double floor, Random rng)
        where TM : struct, IRiemannianManifold
    {
        const int count = 21;
        int dim = n * r;
        double[] w = Enumerable.Repeat(1.0, count).ToArray();
        for (int draw = 0; draw < 3; draw++)
        {
            double[] p0 = Gaussian(rng, dim, 1.0);
            if (sphere) p0 = Scaled(p0, 1.0 / Norm(p0));
            else MatrixOps.Orthonormalize(p0, n, r);
            var u = new double[count][];
            for (int i = 0; i < count; i++) u[i] = Horizontal(p0, Gaussian(rng, dim, 1.0 / Math.Sqrt(dim)), n, r);
            double[] mu = WeightedWeiszfeld(u, w);
            double[] ubar = WeightedMean(u, w);

            for (int k = 1; k <= deepest; k++)
            {
                double s = Math.Pow(10, -k);
                string label = $"{(sphere ? $"S^{n - 1}" : $"Gr({r}, {n})")}, draw {draw}, s {s:G3}";
                double[][] x = u.Select(ui => Exp(manifold, p0, Scaled(ui, s))).ToArray();
                double[] medianTruth = Exp(manifold, p0, Scaled(mu, s));
                double[] meanTruth = Exp(manifold, p0, Scaled(ubar, s));
                double[] warm = x[GeometricMedian.MedoidIndex(manifold, x)];

                var counter = new Counter();
                double[] median = (double[])warm.Clone();
                GeometricMedian.Compute(new Counting<TM>(manifold, counter), x, w, median);
                Assert.True(counter.Steps < 50, $"{label}: median took {counter.Steps} iterations");
                double medianError = CurvedGap(median, medianTruth, n, r, sphere) / s;
                Assert.True(medianError <= 10 * s * s + 1e-6 + floor / s, $"{label}: median error {medianError:G3}");

                double[] mean = (double[])warm.Clone();
                GeometricMean.Compute(manifold, x, w, mean);
                double meanError = CurvedGap(mean, meanTruth, n, r, sphere) / s;
                Assert.True(meanError <= 10 * s * s + 1e-6 + floor / s, $"{label}: mean error {meanError:G3}");
            }
        }
    }

    // Sphere: chordal gap of the normalized estimate. Grassmann: ‖(I − TTᵀ)E‖_F, stable at small angles.
    private static double CurvedGap(double[] estimate, double[] truth, int n, int r, bool sphere)
    {
        if (sphere) return Norm(Subtract(Scaled(estimate, 1.0 / Norm(estimate)), truth));
        return Norm(Horizontal(truth, estimate, n, r));
    }

    // The component of z (n×r, column-major) orthogonal to span(y); for r = 1 the projection off the unit vector y.
    private static double[] Horizontal(double[] y, double[] z, int n, int r)
    {
        var result = (double[])z.Clone();
        for (int j = 0; j < r; j++)
            for (int k = 0; k < r; k++)
            {
                double m = 0;
                for (int i = 0; i < n; i++) m += y[i + k * n] * z[i + j * n];
                for (int i = 0; i < n; i++) result[i + j * n] -= y[i + k * n] * m;
            }
        return result;
    }

    private static double[] Exp<TM>(TM manifold, double[] p, double[] v) where TM : struct, IRiemannianManifold
    {
        var q = new double[p.Length];
        manifold.ExpMap(p, v, q);
        return q;
    }

    // ── Truth from outside the solver ────────────────────────────────────────

    // Plain weighted Weiszfeld at scale 1, from the weighted mean, until the step is below 1e-15 of the data's spread
    // about that mean. On landing on data points it stops if the weight there outweighs the pull of the rest (‖G‖),
    // and otherwise steps with those points weighted 0.
    private static double[] WeightedWeiszfeld(double[][] x, double[] w)
    {
        int d = x[0].Length;
        double[] m = WeightedMean(x, w);
        double stop = 1e-15 * MedianDistance(x, m);
        var num = new double[d];   // Σ wᵢ(xᵢ − m)/rᵢ: the step's numerator and the pull G of the points off m
        for (int it = 0; it < 100_000; it++)
        {
            Array.Clear(num);
            double den = 0, weightOn = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double r = 0;
                for (int j = 0; j < d; j++) { double e = x[i][j] - m[j]; r += e * e; }
                r = Math.Sqrt(r);
                if (r < 1e-300) { weightOn += w[i]; continue; }
                for (int j = 0; j < d; j++) num[j] += w[i] * (x[i][j] - m[j]) / r;
                den += w[i] / r;
            }
            if (weightOn > 0 && Norm(num) <= weightOn) break;
            double shift = 0;
            for (int j = 0; j < d; j++) { double step = num[j] / den; shift += step * step; m[j] += step; }
            if (Math.Sqrt(shift) < stop) break;
        }
        return m;
    }

    private static double[] WeightedMean(double[][] x, double[] w)
    {
        double total = w.Sum();
        var m = new double[x[0].Length];
        for (int i = 0; i < x.Length; i++) for (int j = 0; j < m.Length; j++) m[j] += w[i] * x[i][j] / total;
        return m;
    }

    private static double MedianDistance(double[][] x, double[] centre)
    {
        double[] r = x.Select(xi => Norm(Subtract(xi, centre))).OrderBy(v => v).ToArray();
        return r[r.Length / 2];
    }

    // ‖(estimate − c)/s − truth‖, in the truth's u-units.
    private static double Gap(double[] estimate, double[] c, double s, double[] truth)
    {
        double sum = 0;
        for (int j = 0; j < truth.Length; j++) { double e = (estimate[j] - c[j]) / s - truth[j]; sum += e * e; }
        return Math.Sqrt(sum);
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    private static double[] Gaussian(Random rng, int d, double sd)
    {
        var v = new double[d];
        for (int j = 0; j < d; j++)
        {
            double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
            v[j] = sd * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
        return v;
    }

    private static double[] Scaled(double[] v, double a) => v.Select(x => a * x).ToArray();
    private static double[] Subtract(double[] a, double[] b) => a.Select((x, i) => x - b[i]).ToArray();
    private static double Norm(double[] v) => Math.Sqrt(v.Sum(x => x * x));

    private sealed class Counter { public int Steps; }

    // Counts Exp calls: under the spread unit each iteration takes exactly one step.
    private readonly struct Counting<T> : IRiemannianManifold where T : struct, IRiemannianManifold
    {
        private readonly T _inner;
        private readonly Counter _counter;
        public Counting(T inner, Counter counter) { _inner = inner; _counter = counter; }
        public static bool IsFlat => T.IsFlat;
        public int Dimension => _inner.Dimension;
        public double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q) => _inner.Distance(p, q);
        public void LogMap(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst) => _inner.LogMap(p, q, dst);
        public void ExpMap(ReadOnlySpan<double> p, ReadOnlySpan<double> v, Span<double> dst) { _counter.Steps++; _inner.ExpMap(p, v, dst); }
        public void AddScaled(Span<double> dst, ReadOnlySpan<double> v, double s) => _inner.AddScaled(dst, v, s);
        public double Norm(ReadOnlySpan<double> p, ReadOnlySpan<double> v) => _inner.Norm(p, v);
    }
}
