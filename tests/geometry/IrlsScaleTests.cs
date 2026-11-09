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
    /// A seeded sweep over dimension, n, equal and unequal weights and an occasional far outlier, centred 10 from the
    /// origin at every s and warm-started at the medoid as the callers do. The data are then representable only to
    /// about ε·‖c‖: the median, the mean and <see cref="EuclideanMedian"/> (warm-started at the mean, as GMM does) stay
    /// within a fixed number of those ulps of the truth for the represented data, and no median runs to MaxIterations.
    /// </summary>
    [Fact]
    public void FlatMedianAndMean_FarFromTheOrigin_AreLimitedOnlyByTheCentresRounding()
    {
        var rng = new Random(81102026);
        for (int draw = 0; draw < 6; draw++)
        {
            FlatProblem problem = FlatProblem.Draw(rng);
            double[] c = Scaled(problem.Direction, 10.0);
            foreach (double s in Spreads)
                CheckFlat(problem, c, s, roundingUlps: 1000, $"draw {draw}, centre 10, s {s:G3}");
        }
    }

    /// <summary>
    /// On S⁹ and the Grassmannians Gr(2, 10) and Gr(3, 20), 21 points exp_{p₀}(s·uᵢ) with spreads 1e-1 to 1e-6. The
    /// allowance holds the truth's O(s²) relative error, the 1e-6 accuracy of the flat case, and 1e-10/s for the
    /// manifolds' Exp, which drops tangent steps below 1e-12. No median takes more than 50 iterations.
    /// </summary>
    [Theory]
    [InlineData(10, 0, 7001)]   // S⁹
    [InlineData(10, 2, 7002)]   // Gr(2, 10)
    [InlineData(20, 3, 7003)]   // Gr(3, 20)
    public void CurvedMedianAndMean_AreScaleEquivariant_DownToTheManifoldResolution(int n, int r, int seed)
    {
        var rng = new Random(seed);
        if (r == 0) CheckCurved(new SphericalManifold(n), n, 1, sphere: true, rng);
        else CheckCurved(new GrassmannManifold(n, r), n, r, sphere: false, rng);
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
    /// the spread unit (afc814d) bit for bit: a Euclidean median, a Grassmann median and mean, and
    /// <see cref="EuclideanMedian"/> with its scatter.
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
        Assert.Equal(
            new[]
            {
                0.9999999772287292, 2.6123219160401057E-05, 0.00010966421803994649, -4.216689096171608E-05,
                -0.0001762266467619848, -2.6094562911371732E-05, 0.9999999120481814, 0.00026854019792198436,
                9.108648728956332E-05, 0.00030791576352311054,
            },
            gMedian);

        double[] gMean = (double[])init.Clone();
        GeometricMean.Compute(new GrassmannManifold(5, 2), frames, gw, gMean, data);
        Assert.Equal(
            new[]
            {
                0.9999827000026962, 2.3985945202253938E-08, 0.004271230080841782, -0.0028580817189473025,
                -0.0028614083598075533, 1.2148160713673067E-05, 0.9999692246700764, 0.002891358423609038,
                0.0014173379524098605, 0.007154073297315209,
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
    // roundingUlps·ε·‖c‖/s for a centre far from the origin.
    private static void CheckFlat(FlatProblem problem, double[] c, double s, int roundingUlps, string label)
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

        var counter = new Counter();
        double[] median = (double[])warm.Clone();
        GeometricMedian.Compute(new Counting<EuclideanVectorManifold>(euclid, counter), x, problem.W, median);
        Assert.True(counter.Steps < MaxIterations, $"{label}: median ran to MaxIterations");
        double medianError = Gap(median, c, s, medianTruth);
        Assert.True(medianError <= 1e-6 * spread + rounding, $"{label}: median error {medianError:G3} (spread {spread:G3})");

        double[] mean = (double[])warm.Clone();
        GeometricMean.Compute(euclid, x, problem.W, mean);
        double meanError = Gap(mean, c, s, meanTruth);
        Assert.True(meanError <= 1e-12 * spread + rounding, $"{label}: mean error {meanError:G3}");

        double[] ambient = WeightedMean(x, problem.W);
        EuclideanMedian.Compute(x, Enumerable.Range(0, n).ToArray(), problem.W, ambient);
        double ambientError = Gap(ambient, c, s, medianTruth);
        Assert.True(ambientError <= 1e-6 * spread + rounding, $"{label}: EuclideanMedian error {ambientError:G3}");
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

    private static void CheckCurved<TM>(TM manifold, int n, int r, bool sphere, Random rng)
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

            // Below 1e-6 the arccos in Grassmann's Log drops principal angles under about 1e-8 (out of this test's scope).
            for (int k = 1; k <= 6; k++)
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
                Assert.True(medianError <= 10 * s * s + 1e-6 + 1e-10 / s, $"{label}: median error {medianError:G3}");

                double[] mean = (double[])warm.Clone();
                GeometricMean.Compute(manifold, x, w, mean);
                double meanError = CurvedGap(mean, meanTruth, n, r, sphere) / s;
                Assert.True(meanError <= 10 * s * s + 1e-6, $"{label}: mean error {meanError:G3}");
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

    // Plain weighted Weiszfeld at scale 1, from the weighted mean, to machine precision; it stops on a data point.
    private static double[] WeightedWeiszfeld(double[][] x, double[] w)
    {
        int d = x[0].Length;
        double[] m = WeightedMean(x, w);
        for (int it = 0; it < 100_000; it++)
        {
            var num = new double[d];
            double den = 0;
            bool onPoint = false;
            for (int i = 0; i < x.Length; i++)
            {
                double r = Norm(Subtract(x[i], m));
                if (r < 1e-300) { onPoint = true; continue; }
                for (int j = 0; j < d; j++) num[j] += w[i] * (x[i][j] - m[j]) / r;
                den += w[i] / r;
            }
            double shift = 0;
            for (int j = 0; j < d; j++) { double step = num[j] / den; shift += step * step; m[j] += step; }
            if (Math.Sqrt(shift) < 1e-17 || (onPoint && Math.Sqrt(shift) < 1e-15)) break;
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
