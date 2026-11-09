using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 4 gates (the active-subspace brief) on dependent inputs: a law given as a joint sample (<see cref="JointSample"/>), the
/// sampled route at its points, and whitening by its full covariance (<see cref="Whitening"/>):
/// <list type="bullet">
/// <item>against a product law the exact route already gives: a sample drawn from it, for a scalar fit with a categorical input, a
/// direct functional fit, a basis fit (whose sampled route is new here, product law or sample) and co-active pairs, every entry within
/// Monte Carlo error;</item>
/// <item>whitening under a sample: the readout's C is L̂ᵀC_boxL̂ with L̂ the sample covariance's symmetric root, computed here by
/// Denman–Beavers iteration; standardizing its diagonal; weighted points as replicated ones; the whitened spectrum the same in the
/// inputs' own units as in the box's;</item>
/// <item>equivariance under affine input maps: on reference gradient fields from outside any model, under random invertible maps
/// across scales, a seeded sweep, whitened eigenvalues invariant and the maps between whitened frames orthogonal, where the box
/// metric's spectrum moves;</item>
/// <item>dependent examples: 2307.11241's triangle law (§4.1.1), whose C the paper gives in closed form and reaches only through a
/// mixture of products, here from a BASS fit and a joint sample; and a band law where the box metric and whitening lead in
/// perpendicular directions, both reported;</item>
/// <item>the refusals.</item>
/// </list>
/// </summary>
public sealed class JointSampleActiveSubspaceTests(ITestOutputHelper output)
{
    // ── Samples ──────────────────────────────────────────────────────────────────────────────────────────────────────

    // A joint sample of named continuous columns, a categorical one drawn uniformly over its levels when given.
    private static JointSample SampleOf(string[] names, double[][] columns, double[]? weights = null, (string Name, string[] Levels)? categorical = null, int seed = 1)
    {
        var cols = new List<PredictionColumn>();
        for (int j = 0; j < names.Length; j++) cols.Add(new ContinuousColumn(names[j], columns[j]));
        if (categorical is { } c)
        {
            var rng = new Xoshiro256PlusPlus(seed);
            cols.Add(new CategoricalColumn(c.Name, [.. columns[0].Select(_ => c.Levels[rng.NextInt(c.Levels.Length)])]));
        }
        return new JointSample(new PredictionCases(cols), weights);
    }

    private static double[][] Uniform(Xoshiro256PlusPlus rng, int n, int p)
    {
        var x = new double[p][];
        for (int j = 0; j < p; j++) x[j] = new double[n];
        for (int r = 0; r < n; r++)
            for (int j = 0; j < p; j++) x[j][r] = rng.NextDouble();
        return x;
    }

    // Uniform on the triangle 0 < x₂ < x₁ < 1: the larger and smaller of two uniforms.
    internal static double[][] Triangle(Xoshiro256PlusPlus rng, int n)
    {
        double[] x1 = new double[n], x2 = new double[n];
        for (int r = 0; r < n; r++)
        {
            double u = rng.NextDouble(), v = rng.NextDouble();
            x1[r] = Math.Max(u, v);
            x2[r] = Math.Min(u, v);
        }
        return [x1, x2];
    }

    // Uniform on the band |x₁ − x₂| < h in the unit square, by rejection.
    internal static double[][] Band(Xoshiro256PlusPlus rng, int n, double h)
    {
        double[] x1 = new double[n], x2 = new double[n];
        for (int r = 0; r < n;)
        {
            double u = rng.NextDouble(), v = rng.NextDouble();
            if (Math.Abs(u - v) >= h) continue;
            x1[r] = u;
            x2[r++] = v;
        }
        return [x1, x2];
    }

    private static double TruncatedNormal(Xoshiro256PlusPlus rng, double mean, double sd)
    {
        while (true)
        {
            double x = mean + sd * Normal.SampleStandard(rng);
            if (x >= 0.0 && x <= 1.0) return x;
        }
    }

    private static int[] Spread(int count, int k) => [.. Enumerable.Range(0, k).Select(i => (int)((long)i * (count - 1) / (k - 1)))];

    // ── Against a product law ────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// A scalar fit with a categorical input: a joint sample drawn from the product law of a truncated normal on x0, the uniform law on x1
    /// and uniform levels of z, against the exact route under that law; every entry of every draw's C within 4.5 of the sampled route's
    /// standard errors. Under a sample the route is sampled, the product measure and the checks absent, and N the sample's size.
    /// </summary>
    [Fact]
    public void ScalarFit_SampleOfAProductLaw_MeetsTheExactRoute()
    {
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = SensitivityFixture.FitScalar(1);
        const int n = 16384;
        var rng = new Xoshiro256PlusPlus(41);
        double[] x0 = new double[n], x1 = new double[n];
        for (int r = 0; r < n; r++)
        {
            x0[r] = TruncatedNormal(rng, 0.4, 0.25);
            x1[r] = rng.NextDouble();
        }
        JointSample sample = SampleOf(["x0", "x1"], [x0, x1], categorical: ("z", ["a", "b", "c"]), seed: 43);
        int[] draws = Spread(fit.Draws.Sum(c => c.Count), 12);
        var spec = new ActiveSubspaceSpec { Draws = draws, Dimension = 2 };
        ActiveSubspaceAnalysis exact = fit.ActiveSubspace(data, spec with
        {
            Measures = new Dictionary<string, InputMeasure> { ["x0"] = new MixtureMeasure([new NormalComponent(1.0, 0.4, 0.25)]) },
        });
        ActiveSubspaceAnalysis sampled = fit.ActiveSubspace(data, spec with { Sample = sample });
        double worst = Worst(exact, sampled);
        output.WriteLine($"{draws.Length} draws, {n} points: worst |exact − sampled| / se {worst:F2}");
        Assert.True(worst <= 4.5, $"{worst:F2}");
        Assert.Equal(ActiveSubspaceRoute.Sampled, sampled.Route);
        Assert.Null(sampled.Measure);
        Assert.Null(sampled.Poincare);
        Assert.Equal(n, sampled.Points);
        Assert.Equal(["z"], sampled.Excluded);
    }

    /// <summary>
    /// A direct functional fit under a full W_y, and a basis fit, against their exact routes: the direct fit under a uniform joint
    /// sample; the basis fit on its new sampled route, by the product law's own points and by a joint sample.
    /// </summary>
    [Fact]
    public void FunctionalAndBasisFits_SampleOfAProductLaw_MeetTheExactRoute()
    {
        var rng = new Xoshiro256PlusPlus(47);
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        var metric = new OutputMetric(FunctionalActiveSubspaceTests.RandomMetric(data.Response.Locations, 53));
        var spec = new ActiveSubspaceSpec { Draws = Spread(fit.Draws.Sum(c => c.Count), 8), OutputMetric = metric };
        double direct = Worst(fit.ActiveSubspace(data, spec), fit.ActiveSubspace(data, spec with { Sample = SampleOf(["x0", "x1"], Uniform(rng, 8192, 2)) }));

        BassBasisFit basis = BassOutputComponentsTests.FitOf();
        var basisSpec = new ActiveSubspaceSpec { Draws = Spread(basis.JointDraws.Count, 6), OutputMetric = new OutputMetric(FunctionalActiveSubspaceTests.RandomMetric(basis.Decomposition.Locations, 59)) };
        ActiveSubspaceAnalysis exact = basis.ActiveSubspace(basisSpec);
        ActiveSubspaceAnalysis drawn = basis.ActiveSubspace(basisSpec with { Route = ActiveSubspaceRoute.Sampled, Points = 8192, Seed = 7 });
        ActiveSubspaceAnalysis given = basis.ActiveSubspace(basisSpec with { Sample = SampleOf(["z0", "z1"], Uniform(rng, 8192, 2)) });
        double byLaw = Worst(exact, drawn), bySample = Worst(exact, given);
        output.WriteLine($"direct functional, uniform sample: worst {direct:F2} se; basis, the law's points: {byLaw:F2} se, a uniform sample: {bySample:F2} se; " +
                         $"basis sampled checks: {drawn.Poincare!.Violations} violations");
        Assert.True(direct <= 4.5, $"direct {direct:F2}");
        Assert.True(byLaw <= 4.5, $"basis, law {byLaw:F2}");
        Assert.True(bySample <= 4.5, $"basis, sample {bySample:F2}");
        Assert.Equal(0, drawn.Poincare!.Violations);
        Assert.NotNull(drawn.MonteCarloError);
    }

    /// <summary>Co-active pairs under a uniform joint sample against the exact route: two scalar fits, and two outputs of a basis fit by
    /// a sample and by its new sampled route.</summary>
    [Fact]
    public void CoActivePairs_SampleOfAProductLaw_MeetTheExactRoute()
    {
        var rng = new Xoshiro256PlusPlus(61);
        (BassData df, BassData dg) = CoActiveFixture.Data();
        (BassFit f, BassFit g) = CoActiveFixture.Fits();
        int[] sel = Spread(f.Draws.Sum(c => c.Count), 8);
        var spec = new CoActiveSpec { Draws = sel, OtherDraws = sel };
        double fits = Worst(f.CoActiveSubspace(df, g, dg, spec), f.CoActiveSubspace(df, g, dg, spec with { Sample = SampleOf(["x0", "x1"], Uniform(rng, 8192, 2)) }));

        BassBasisFit basis = BassOutputComponentsTests.FitOf();
        int m = basis.Decomposition.Locations;
        double[] a = [.. Enumerable.Range(0, m).Select(l => 0.2 + 0.3 * Math.Sin(1.3 * l))], b = [.. Enumerable.Range(0, m).Select(l => l % 3 == 0 ? 1.0 : -0.4)];
        var basisSpec = new CoActiveSpec { Draws = Spread(basis.JointDraws.Count, 5), Output = a, OtherOutput = b };
        CoActiveAnalysis exact = basis.CoActiveSubspace(basisSpec);
        double byLaw = Worst(exact, basis.CoActiveSubspace(basisSpec with { Route = ActiveSubspaceRoute.Sampled, Points = 8192, Seed = 3 }));
        double bySample = Worst(exact, basis.CoActiveSubspace(basisSpec with { Sample = SampleOf(["z0", "z1"], Uniform(rng, 8192, 2)) }));
        output.WriteLine($"two fits, uniform sample: worst {fits:F2} se; basis outputs, the law's points {byLaw:F2} se, a sample {bySample:F2} se");
        Assert.True(fits <= 4.5, $"{fits:F2}");
        Assert.True(byLaw <= 4.5, $"{byLaw:F2}");
        Assert.True(bySample <= 4.5, $"{bySample:F2}");
    }

    // ── Whitening under a sample ─────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Under a weighted joint sample from the triangle law: the whitened readout's map is the sample covariance's symmetric root L̂,
    /// computed here by Denman–Beavers iteration, and each draw's C is L̂ᵀC_boxL̂; standardizing takes L̂'s place with diag(σ̂); and the
    /// sample with integer weights is the sample with each point repeated.
    /// </summary>
    [Fact]
    public void Whitened_UnderASample_IsTheCovariancesRootsCongruence()
    {
        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        var rng = new Xoshiro256PlusPlus(67);
        const int n = 3000;
        double[][] x = Triangle(rng, n);
        double[] weights = [.. Enumerable.Range(0, n).Select(_ => (double)(1 + rng.NextInt(3)))];
        JointSample sample = SampleOf(["x0", "x1"], x, weights);
        var spec = new ActiveSubspaceSpec { Draws = Spread(f.Draws.Sum(c => c.Count), 6), Sample = sample, Dimension = 2 };
        ActiveSubspaceAnalysis box = f.ActiveSubspace(df, spec);
        ActiveSubspaceAnalysis white = f.ActiveSubspace(df, spec with { Metric = InputMetric.Whitened });
        ActiveSubspaceAnalysis std = f.ActiveSubspace(df, spec with { Metric = InputMetric.Standardized });

        double[,] sigma = WeightedCovariance(x, weights);
        double[,] root = DenmanBeavers(sigma);
        double[,] diag = { { Math.Sqrt(sigma[0, 0]), 0.0 }, { 0.0, Math.Sqrt(sigma[1, 1]) } };
        double mapGap = FunctionalActiveSubspaceTests.RelativeGap(white.Map, root), worst = 0.0;
        for (int d = 0; d < box.Posterior.Draws.Count; d++)
        {
            double[,] c = box.Posterior.Draws[d].C;
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(white.Posterior.Draws[d].C, Congruence(c, root)));
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(std.Posterior.Draws[d].C, Congruence(c, diag)));
        }

        // Integer weights as repeats.
        var repeated = new double[2][];
        for (int j = 0; j < 2; j++) repeated[j] = [.. Enumerable.Range(0, n).SelectMany(r => Enumerable.Repeat(x[j][r], (int)weights[r]))];
        ActiveSubspaceAnalysis replicated = f.ActiveSubspace(df, spec with { Metric = InputMetric.Whitened, Sample = SampleOf(["x0", "x1"], repeated) });
        double replication = FunctionalActiveSubspaceTests.RelativeGap(replicated.Map, white.Map);
        for (int d = 0; d < box.Posterior.Draws.Count; d++)
            replication = Math.Max(replication, FunctionalActiveSubspaceTests.RelativeGap(replicated.Posterior.Draws[d].C, white.Posterior.Draws[d].C));
        output.WriteLine($"map against Denman–Beavers {mapGap:E2}; whitened and standardized against the congruences {worst:E2}; weights against repeats {replication:E2}; " +
                         $"the sample's correlation {sigma[0, 1] / Math.Sqrt(sigma[0, 0] * sigma[1, 1]):F3}");
        Assert.True(mapGap <= 1e-12, $"map {mapGap:E2}");
        Assert.True(worst <= 1e-12, $"congruence {worst:E2}");
        Assert.True(replication <= 1e-12, $"repeats {replication:E2}");
        Assert.Equal(InputMetric.Whitened, white.Metric);
        Assert.NotEqual(0.0, white.Map[0, 1]);
    }

    /// <summary>
    /// The whitened spectrum does not depend on the inputs' units: on a fit whose boxes are [−3, 5] and [10, 12], under a correlated
    /// sample, each draw's whitened eigenvalues against those of L_rawᵀC_rawL_raw with C_raw the raw-unit readout and L_raw the root of
    /// the sample's covariance in raw units; the box and raw spectra differ.
    /// </summary>
    [Fact]
    public void Whitened_IsTheSameInTheInputsOwnUnits()
    {
        (BassData data, BassFit fit) = WideFixture.Value;
        var rng = new Xoshiro256PlusPlus(71);
        const int n = 4096;
        double[] x0 = new double[n], x1 = new double[n];
        for (int r = 0; r < n; r++)
        {
            double u = rng.NextDouble(), v = rng.NextDouble();
            x0[r] = -3.0 + 8.0 * u;
            x1[r] = 10.0 + 2.0 * (0.7 * u + 0.3 * v);
        }
        var spec = new ActiveSubspaceSpec { Draws = Spread(fit.Draws.Sum(c => c.Count), 6), Sample = SampleOf(["x0", "x1"], [x0, x1]), Dimension = 2 };
        ActiveSubspaceAnalysis raw = fit.ActiveSubspace(data, spec with { Metric = InputMetric.Raw });
        ActiveSubspaceAnalysis box = fit.ActiveSubspace(data, spec);
        ActiveSubspaceAnalysis white = fit.ActiveSubspace(data, spec with { Metric = InputMetric.Whitened });
        double[,] rawRoot = DenmanBeavers(WeightedCovariance([x0, x1], null));
        double worst = 0.0, boxRaw = 0.0;
        for (int d = 0; d < white.Posterior.Draws.Count; d++)
        {
            double[] expected = Eigenvalues(Congruence(raw.Posterior.Draws[d].C, rawRoot));
            double[] actual = white.Posterior.Draws[d].Eigenvalues;
            for (int i = 0; i < 2; i++) worst = Math.Max(worst, Math.Abs(actual[i] - expected[i]) / expected[0]);
            boxRaw = Math.Max(boxRaw, Math.Abs(box.Posterior.Draws[d].Eigenvalues[0] / raw.Posterior.Draws[d].Eigenvalues[0] - 1.0));
        }
        output.WriteLine($"whitened eigenvalues against raw units' {worst:E2}; the box's leading eigenvalue against raw units' differs by up to {boxRaw:P0}");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
        Assert.True(boxRaw > 1.0, $"{boxRaw:F3}");
    }

    private static readonly Lazy<(BassData, BassFit)> WideFixture = new(() =>
    {
        var rng = new Xoshiro256PlusPlus(73);
        const int n = 80;
        double[] x0 = new double[n], x1 = new double[n], y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = -3.0 + 8.0 * (i % 10) / 9.0;
            x1[i] = 10.0 + 2.0 * (i * 7 % 11) / 10.0;
            y[i] = Math.Sin(0.5 * x0[i]) + 2.0 * (x1[i] - 11.0) * (x1[i] - 11.0) + 0.3 * x0[i] * (x1[i] - 11.0) + 0.05 * (rng.NextDouble() - 0.5);
        }
        var data = new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
        return (data, Bass.Fit(SensitivityFixture.Config(), data));
    });

    // ── Affine equivariance, fit-agnostic ────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// A seeded sweep of 200 cases: p from 2 to 5, a correlated sample of 50 to 400 points (weighted in half the cases) and a reference
    /// gradient field (a quadratic plus sine ridges, from outside any model), mapped by a random invertible z = Ax + c with scales from
    /// 10⁻³ to 10³. Whitening by <see cref="Whitening"/>: C_u = L_xᵀC_xL_x against C_u′ = L_zᵀC_zL_z with C_z = A⁻ᵀC_xA⁻¹ the same
    /// field's in z. The eigenvalues agree, and Q = L_x⁻¹A⁻¹L_z is orthogonal with C_u′ = QᵀC_uQ: whitened subspaces, read as functions
    /// of the inputs, are the same in both frames. The box metric's spectrum, C_x against C_z, moves.
    /// </summary>
    [Fact]
    public void Whitening_IsAffineEquivariant_SeededSweep()
    {
        var rng = new Xoshiro256PlusPlus(20261108);
        double worstEigen = 0.0, worstOrthogonal = 0.0, worstCongruence = 0.0, boxMove = 0.0;
        for (int s = 0; s < 200; s++)
        {
            int p = 2 + rng.NextInt(4), n = 50 + rng.NextInt(351);
            double[,] mix = RandomMatrix(rng, p, 1.0);
            var x = new double[p][];
            for (int j = 0; j < p; j++) x[j] = new double[n];
            for (int r = 0; r < n; r++)
            {
                double[] u = [.. Enumerable.Range(0, p).Select(_ => rng.NextDouble() < 0.5 ? Normal.SampleStandard(rng) : 2.0 * rng.NextDouble() - 1.0)];
                for (int j = 0; j < p; j++)
                    for (int k = 0; k < p; k++) x[j][r] += mix[j, k] * u[k];
            }
            double[]? w = s % 2 == 0 ? null : [.. Enumerable.Range(0, n).Select(_ => 0.1 + rng.NextDouble())];
            Func<double[], double[]> grad = Field(rng, p);
            double scale = Math.Pow(10.0, -3.0 + 6.0 * rng.NextDouble());
            double[,] a = RandomMatrix(rng, p, scale);
            double[] c = [.. Enumerable.Range(0, p).Select(_ => scale * (rng.NextDouble() - 0.5) * 10.0)];
            double[,] aInv = Inverse(a);

            double[,] cx = new double[p, p], cz = new double[p, p];
            var z = new double[p][];
            for (int j = 0; j < p; j++) z[j] = new double[n];
            double[] wn = w is null ? [.. Enumerable.Repeat(1.0 / n, n)] : [.. w.Select(v => v / w.Sum())];
            for (int r = 0; r < n; r++)
            {
                double[] point = [.. Enumerable.Range(0, p).Select(j => x[j][r])];
                double[] g = grad(point);
                double[] gz = new double[p];   // ∇_z = A⁻ᵀ∇_x
                for (int i = 0; i < p; i++)
                    for (int k = 0; k < p; k++) gz[i] += aInv[k, i] * g[k];
                for (int j = 0; j < p; j++)
                {
                    z[j][r] = c[j];
                    for (int k = 0; k < p; k++) z[j][r] += a[j, k] * point[k];
                }
                for (int i = 0; i < p; i++)
                    for (int j = 0; j < p; j++)
                    {
                        cx[i, j] += wn[r] * g[i] * g[j];
                        cz[i, j] += wn[r] * gz[i] * gz[j];
                    }
            }
            double[,] lx = Whitening.SymmetricRoot(Whitening.Moments(x, w).Covariance), lz = Whitening.SymmetricRoot(Whitening.Moments(z, w).Covariance);
            double[,] cu = Congruence(cx, lx), cuz = Congruence(cz, lz);
            double[] ex = Eigenvalues(cu), ez = Eigenvalues(cuz);
            for (int i = 0; i < p; i++) worstEigen = Math.Max(worstEigen, Math.Abs(ex[i] - ez[i]) / ex[0]);
            double[,] q = Multiply(Multiply(Inverse(lx), aInv), lz);
            double[,] qtq = Multiply(Transpose(q), q);
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) worstOrthogonal = Math.Max(worstOrthogonal, Math.Abs(qtq[i, j] - (i == j ? 1.0 : 0.0)));
            worstCongruence = Math.Max(worstCongruence, FunctionalActiveSubspaceTests.RelativeGap(cuz, Multiply(Multiply(Transpose(q), cu), q)));
            double[] bx = Eigenvalues(cx), bz = Eigenvalues(cz);
            boxMove = Math.Max(boxMove, Math.Abs(Math.Log(bz[0] / bx[0])));
        }
        output.WriteLine($"200 cases: whitened eigenvalues {worstEigen:E2}, QᵀQ − I {worstOrthogonal:E2}, C_u′ against QᵀC_uQ {worstCongruence:E2}; " +
                         $"the box spectrum's leading eigenvalue moves by a factor up to {Math.Exp(boxMove):E1}");
        Assert.True(worstEigen <= 1e-9, $"eigenvalues {worstEigen:E2}");
        Assert.True(worstOrthogonal <= 1e-9, $"orthogonality {worstOrthogonal:E2}");
        Assert.True(worstCongruence <= 1e-9, $"congruence {worstCongruence:E2}");
        Assert.True(boxMove > 1.0);
    }

    // A gradient field from outside any model: ½xᵀQx + bᵀx plus two sine ridges.
    private static Func<double[], double[]> Field(Xoshiro256PlusPlus rng, int p)
    {
        double[,] q = RandomMatrix(rng, p, 1.0);
        double[] b = [.. Enumerable.Range(0, p).Select(_ => Normal.SampleStandard(rng))];
        double[][] w = [[.. Enumerable.Range(0, p).Select(_ => Normal.SampleStandard(rng))], [.. Enumerable.Range(0, p).Select(_ => Normal.SampleStandard(rng))]];
        return x =>
        {
            var g = new double[p];
            for (int i = 0; i < p; i++)
            {
                g[i] = b[i];
                for (int k = 0; k < p; k++) g[i] += 0.5 * (q[i, k] + q[k, i]) * x[k];
            }
            foreach (double[] wk in w)
            {
                double t = 0.0;
                for (int k = 0; k < p; k++) t += wk[k] * x[k];
                double slope = Math.Cos(1.5 * t) * 1.5;
                for (int i = 0; i < p; i++) g[i] += slope * wk[i];
            }
            return g;
        };
    }

    /// <summary>
    /// The primitives: the symmetric root squares to Σ and is symmetric, on chosen and random covariances with condition numbers to 10⁸
    /// and scales from 10⁻⁶ to 10⁶; weighted moments are replicated ones; and the refusals: a singular covariance, an asymmetric one, a
    /// non-finite entry, negative or zero-sum weights.
    /// </summary>
    [Fact]
    public void Whitening_Primitives()
    {
        var rng = new Xoshiro256PlusPlus(79);
        double worstSquare = 0.0, worstSymmetry = 0.0;
        for (int s = 0; s < 120; s++)
        {
            int p = 1 + rng.NextInt(6);
            double scale = Math.Pow(10.0, -6.0 + 12.0 * rng.NextDouble()), condition = Math.Pow(10.0, 8.0 * rng.NextDouble());
            double[,] v = Orthogonal(rng, p);
            var sigma = new double[p, p];
            for (int k = 0; k < p; k++)
            {
                double lambda = scale * Math.Pow(condition, -(double)k / Math.Max(1, p - 1));
                for (int i = 0; i < p; i++)
                    for (int j = 0; j < p; j++) sigma[i, j] += lambda * v[i, k] * v[j, k];
            }
            for (int i = 0; i < p; i++)
                for (int j = i + 1; j < p; j++) sigma[j, i] = sigma[i, j];
            double[,] root = Whitening.SymmetricRoot(sigma);
            worstSquare = Math.Max(worstSquare, FunctionalActiveSubspaceTests.RelativeGap(Multiply(root, root), sigma));
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) worstSymmetry = Math.Max(worstSymmetry, Math.Abs(root[i, j] - root[j, i]));
        }
        double[][] cols = [[0.1, 0.5, 0.9, 0.3], [2.0, 1.0, 4.0, 3.0]];
        (double[] mw, double[,] cw) = Whitening.Moments(cols, [1.0, 2.0, 1.0, 3.0]);
        (double[] mr, double[,] cr) = Whitening.Moments([[0.1, 0.5, 0.5, 0.9, 0.3, 0.3, 0.3], [2.0, 1.0, 1.0, 4.0, 3.0, 3.0, 3.0]]);
        output.WriteLine($"root squared against Σ {worstSquare:E2}, symmetry {worstSymmetry:E2}; weighted against repeated {FunctionalActiveSubspaceTests.RelativeGap(cw, cr):E2}");
        Assert.True(worstSquare <= 1e-8, $"{worstSquare:E2}");
        Assert.Equal(0.0, worstSymmetry);
        Assert.True(FunctionalActiveSubspaceTests.RelativeGap(cw, cr) <= 1e-15);
        Assert.Equal(mr[0], mw[0], 15);
        Assert.Throws<ArgumentException>(() => Whitening.SymmetricRoot(new[,] { { 1.0, 1.0 }, { 1.0, 1.0 } }));
        Assert.Throws<ArgumentException>(() => Whitening.SymmetricRoot(new[,] { { 1.0, 0.2 }, { 0.3, 1.0 } }));
        Assert.Throws<ArgumentException>(() => Whitening.SymmetricRoot(new[,] { { 1.0, double.NaN }, { double.NaN, 1.0 } }));
        Assert.Throws<ArgumentException>(() => Whitening.Moments(cols, [1.0, -1.0, 1.0, 1.0]));
        Assert.Throws<ArgumentException>(() => Whitening.Moments(cols, [0.0, 0.0, 0.0, 0.0]));
        Assert.Throws<ArgumentException>(() => Whitening.Moments(cols, [1.0, 1.0]));
    }

    // ── Dependent examples ───────────────────────────────────────────────────────────────────────────────────────────

    // A region's uniform law by quadrature: x₁ on panels split at the breaks, x₂ between the region's limits at x₁; exact for
    // polynomials of the rule's degree. Returns points and weights summing to 1.
    private static (double[][] X, double[] W) RegionRule(Func<double, (double Lower, double Upper)> limits, double[] breaks, int nodes)
    {
        GaussLegendre rule = GaussLegendre.Of(nodes);
        var x1 = new List<double>();
        var x2 = new List<double>();
        var w = new List<double>();
        for (int k = 0; k + 1 < breaks.Length; k++)
        {
            double a = breaks[k], b = breaks[k + 1];
            for (int i = 0; i < nodes; i++)
            {
                double u = 0.5 * (a + b) + 0.5 * (b - a) * rule.Nodes[i], wu = 0.5 * (b - a) * rule.Weights[i];
                (double lo, double hi) = limits(u);
                for (int j = 0; j < nodes; j++)
                {
                    x1.Add(u);
                    x2.Add(0.5 * (lo + hi) + 0.5 * (hi - lo) * rule.Nodes[j]);
                    w.Add(wu * 0.5 * (hi - lo) * rule.Weights[j]);
                }
            }
        }
        double total = w.Sum();
        return ([[.. x1], [.. x2]], [.. w.Select(v => v / total)]);
    }

    private static double[,] RuleC(PolynomialReference f, (double[][] X, double[] W) rule)
    {
        var c = new double[2, 2];
        var g = new double[2];
        for (int r = 0; r < rule.W.Length; r++)
        {
            f.Gradient([rule.X[0][r], rule.X[1][r]], g);
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++) c[i, j] += rule.W[r] * g[i] * g[j];
        }
        return c;
    }

    // A BASS fit of f on 500 Latin hypercube points of the unit square, on the declared domain [0, 1], noise free.
    private static (BassData Data, BassFit Fit) FitOnTheSquare(PolynomialReference f, int dataSeed, int chainSeed)
    {
        var rng = new Xoshiro256PlusPlus(dataSeed);
        const int n = 500;
        double[][] x = [Latin(rng, n), Latin(rng, n)];
        double[] y = [.. Enumerable.Range(0, n).Select(r => f.Value([x[0][r], x[1][r]]))];
        var domain = new InputDomain(0.0, 1.0);
        var data = new BassData([new ContinuousPredictor("x1", x[0], domain), new ContinuousPredictor("x2", x[1], domain)], new ScalarResponse(y));
        BassConfig config = SensitivityFixture.Config(maxBasis: 40) with
        {
            Sampling = new BassSampling { Chains = 2, MasterSeed = chainSeed, BurnIn = 2_000, Samples = 30, Thin = 10 },
        };
        return (data, Bass.Fit(config, data));
    }

    private static double[] Latin(Xoshiro256PlusPlus rng, int n)
    {
        int[] order = [.. Enumerable.Range(0, n)];
        for (int i = n - 1; i > 0; i--)
        {
            int j = rng.NextInt(i + 1);
            (order[i], order[j]) = (order[j], order[i]);
        }
        return [.. order.Select(k => (k + rng.NextDouble()) / n)];
    }

    /// <summary>
    /// 2307.11241's linear constraint (§4.1.1): f = x₁² + x₁x₂ + x₂³/9 under ρ₂ = 2·𝟙(0 &lt; x₂ &lt; x₁ &lt; 1), whose C the paper writes as
    /// [[1710, 741], [741, 322]]/540 and reaches only through a mixture of L products. Here: the closed form against quadrature of the
    /// triangle; a BASS fit on the square read under a joint sample of the law, E[C] in the box metric against it and, whitened, against
    /// L̂ᵀCL̂ with L̂ the sample covariance's root. C is nearly rank one (λ₁ carries 99.96% of tr C), so the leading direction, read as a
    /// function of the inputs, is the same in both metrics: an exactly rank-k C's subspace is the same set of input functionals in every
    /// metric, which then decides only among the directions a C of higher rank leaves.
    /// </summary>
    [Fact]
    public void TriangleLaw_ThePapersLinearConstraint()
    {
        var f = new PolynomialReference("2307.11241 eq. simple", [new Monomial(1.0, [2, 0]), new Monomial(1.0, [1, 1]), new Monomial(1.0 / 9.0, [0, 3])]);
        (double[][] X, double[] W) rule = RegionRule(u => (0.0, u), [0.0, 1.0], 6);
        double[,] truth = RuleC(f, rule);
        double[,] paper = { { 1710.0 / 540, 741.0 / 540 }, { 741.0 / 540, 322.0 / 540 } };
        Assert.True(FunctionalActiveSubspaceTests.RelativeGap(truth, paper) <= 1e-14, "the paper's closed form");
        double[,] sigmaTruth = WeightedCovariance(rule.X, rule.W);
        Assert.True(FunctionalActiveSubspaceTests.RelativeGap(sigmaTruth, new[,] { { 1.0 / 18, 1.0 / 36 }, { 1.0 / 36, 1.0 / 18 } }) <= 1e-14, "Σ");

        (BassData data, BassFit fit) = FitOnTheSquare(f, 20231117, 29);
        double[][] x = Triangle(new Xoshiro256PlusPlus(83), 20000);
        var spec = new ActiveSubspaceSpec { Sample = SampleOf(["x1", "x2"], x) };
        ActiveSubspaceAnalysis box = fit.ActiveSubspace(data, spec);
        ActiveSubspaceAnalysis white = fit.ActiveSubspace(data, spec with { Metric = InputMetric.Whitened });
        double[,] l = DenmanBeavers(WeightedCovariance(x, null));
        double boxError = Frobenius(box.Posterior.MeanC, truth) / Frobenius(truth, new double[2, 2]);
        double[,] whiteTruth = Congruence(truth, l);
        double whiteError = Frobenius(white.Posterior.MeanC, whiteTruth) / Frobenius(whiteTruth, new double[2, 2]);
        double[] boxLead = Lead(box.Posterior.MeanCVectors), whiteLead = Functional(white.Map, Lead(white.Posterior.MeanCVectors));
        double between = Angle(boxLead, whiteLead), toTruth = Angle(boxLead, Lead(Vectors(truth)));
        output.WriteLine($"box: E[C] = [[{box.Posterior.MeanC[0, 0]:F4}, {box.Posterior.MeanC[0, 1]:F4}], [·, {box.Posterior.MeanC[1, 1]:F4}]] against " +
                         $"[[{truth[0, 0]:F4}, {truth[0, 1]:F4}], [·, {truth[1, 1]:F4}]], relative error {boxError:P2}; leading direction ({boxLead[0]:F4}, {boxLead[1]:F4}), " +
                         $"{toTruth:E2} rad from the truth's; whitened: relative error {whiteError:P2}, leading direction as an input functional " +
                         $"({whiteLead[0]:F4}, {whiteLead[1]:F4}), {between:E2} rad from the box's; whitened eigenvalues {white.Posterior.MeanCEigenvalues[0]:F4}, " +
                         $"{white.Posterior.MeanCEigenvalues[1]:F5}, box {box.Posterior.MeanCEigenvalues[0]:F4}, {box.Posterior.MeanCEigenvalues[1]:F5}");
        Assert.True(boxError < 0.05, $"box {boxError:P2}");
        Assert.True(whiteError < 0.05, $"whitened {whiteError:P2}");
        Assert.True(toTruth < 0.02, $"{toTruth:E2}");
        Assert.True(between < 0.02, $"{between:E2}");
    }

    /// <summary>
    /// A law and a function where the box metric and whitening lead in perpendicular directions: f = 16(x₁ − x₂)² + (x₁ + x₂)² under the
    /// uniform law on the band |x₁ − x₂| &lt; 0.2 (correlation 0.93). f changes fastest across the band, so the box metric leads with
    /// x₁ − x₂; the law barely moves that way, so whitening, which weighs each direction by how far the law moves along it, leads with
    /// x₁ + x₂. Both truths by quadrature of the band; a BASS fit on the square read under a joint sample; each readout's leading
    /// direction, as an input functional, against its own truth; both reported.
    /// </summary>
    [Fact]
    public void BandLaw_BoxAndWhitenedLeadPerpendicularly()
    {
        const double h = 0.2;
        var f = new PolynomialReference("band", [new Monomial(17.0, [2, 0]), new Monomial(-30.0, [1, 1]), new Monomial(17.0, [0, 2])]);   // 16(x₁−x₂)² + (x₁+x₂)²
        (double[][] X, double[] W) rule = RegionRule(u => (Math.Max(0.0, u - h), Math.Min(1.0, u + h)), [0.0, h, 1.0 - h, 1.0], 6);
        double[,] truth = RuleC(f, rule), sigma = WeightedCovariance(rule.X, rule.W);
        double[,] lTruth = DenmanBeavers(sigma);
        double[] boxTruth = Lead(Vectors(truth)), whiteTruth = Functional(lTruth, Lead(Vectors(Congruence(truth, lTruth))));

        (BassData data, BassFit fit) = FitOnTheSquare(f, 20231119, 31);
        double[][] x = Band(new Xoshiro256PlusPlus(89), 20000, h);
        var spec = new ActiveSubspaceSpec { Sample = SampleOf(["x1", "x2"], x) };
        ActiveSubspaceAnalysis box = fit.ActiveSubspace(data, spec);
        ActiveSubspaceAnalysis white = fit.ActiveSubspace(data, spec with { Metric = InputMetric.Whitened });
        double[] boxLead = Lead(box.Posterior.MeanCVectors), whiteLead = Functional(white.Map, Lead(white.Posterior.MeanCVectors));
        double boxAngle = Angle(boxLead, boxTruth), whiteAngle = Angle(whiteLead, whiteTruth), between = Angle(boxLead, whiteLead);
        double gap = box.Posterior.MeanCEigenvalues[1] / box.Posterior.MeanCEigenvalues[0], whiteGap = white.Posterior.MeanCEigenvalues[1] / white.Posterior.MeanCEigenvalues[0];
        output.WriteLine($"truths: box leads ({boxTruth[0]:F4}, {boxTruth[1]:F4}), whitened ({whiteTruth[0]:F4}, {whiteTruth[1]:F4}); the sample's correlation " +
                         $"{sigma[0, 1] / Math.Sqrt(sigma[0, 0] * sigma[1, 1]):F3}; the readouts: box ({boxLead[0]:F4}, {boxLead[1]:F4}) {boxAngle:F4} rad from its truth, " +
                         $"λ₂/λ₁ {gap:F3}; whitened ({whiteLead[0]:F4}, {whiteLead[1]:F4}) {whiteAngle:F4} rad from its truth, λ₂/λ₁ {whiteGap:F4}; " +
                         $"{between:F4} rad apart");
        Assert.True(Angle(boxTruth, whiteTruth) > 1.5, "the truths lead perpendicularly");
        Assert.True(boxAngle < 0.1, $"box {boxAngle:F4}");
        Assert.True(whiteAngle < 0.1, $"whitened {whiteAngle:F4}");
        Assert.True(between > 1.4, $"{between:F4}");
    }

    // ── Refusals ─────────────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The readout refuses a sample beside Measures; the exact route under a sample; points outside the boxes; weights of the wrong count,
    /// negative, of zero sum, or with one positive; a missing column; whitening a degenerate sample, and standardizing an input that does
    /// not vary over it. The box metric reads the degenerate sample.
    /// </summary>
    [Fact]
    public void JointSample_RefusesWhatItDoesNotRead()
    {
        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        var rng = new Xoshiro256PlusPlus(97);
        double[][] x = Uniform(rng, 200, 2);
        var spec = new ActiveSubspaceSpec { Draws = [0, 1], Sample = SampleOf(["x0", "x1"], x) };
        Assert.Throws<ArgumentException>(() => f.ActiveSubspace(df, spec with { Measures = new Dictionary<string, InputMeasure> { ["x0"] = new UniformMeasure() } }));
        Assert.Throws<ArgumentException>(() => f.ActiveSubspace(df, spec with { Route = ActiveSubspaceRoute.Exact }));
        double[][] outside = [[.. x[0]], [.. x[1]]];
        outside[0][5] = 1.5;
        Assert.Throws<ArgumentException>(() => f.ActiveSubspace(df, spec with { Sample = SampleOf(["x0", "x1"], outside) }));
        Assert.Throws<ArgumentException>(() => f.ActiveSubspace(df, spec with { Sample = SampleOf(["x0", "x1"], x, [1.0, 2.0]) }));
        Assert.Throws<ArgumentException>(() => f.ActiveSubspace(df, spec with { Sample = SampleOf(["x0", "x1"], x, [.. Enumerable.Range(0, 200).Select(r => r == 3 ? -1.0 : 1.0)]) }));
        Assert.Throws<ArgumentException>(() => f.ActiveSubspace(df, spec with { Sample = SampleOf(["x0", "x1"], x, new double[200]) }));
        Assert.Throws<ArgumentException>(() => f.ActiveSubspace(df, spec with { Sample = SampleOf(["x0", "x1"], x, [.. Enumerable.Range(0, 200).Select(r => r == 7 ? 1.0 : 0.0)]) }));
        Assert.Throws<ArgumentException>(() => f.ActiveSubspace(df, spec with { Sample = SampleOf(["x0"], [x[0]]) }));
        double[][] diagonal = [[.. x[0]], [.. x[0]]];
        Assert.Throws<ArgumentException>(() => f.ActiveSubspace(df, spec with { Sample = SampleOf(["x0", "x1"], diagonal), Metric = InputMetric.Whitened }));
        double[][] flat = [[.. x[0]], [.. Enumerable.Repeat(0.5, 200)]];
        Assert.Throws<ArgumentException>(() => f.ActiveSubspace(df, spec with { Sample = SampleOf(["x0", "x1"], flat), Metric = InputMetric.Standardized }));
        Assert.Equal(ActiveSubspaceRoute.Sampled, f.ActiveSubspace(df, spec with { Sample = SampleOf(["x0", "x1"], diagonal) }).Route);
        Assert.Throws<ArgumentException>(() => SensitivityFixture.FitScalar(1).ActiveSubspace(SensitivityFixture.Scalar(), spec));   // no z column
    }

    // ── Helpers ──────────────────────────────────────────────────────────────────────────────────────────────────────

    // The worst |exact − sampled|/se over each draw's C (an analysis) or each pair's three matrices.
    private static double Worst(ActiveSubspaceAnalysis exact, ActiveSubspaceAnalysis sampled)
    {
        double worst = 0.0;
        for (int d = 0; d < exact.Posterior.Draws.Count; d++)
            worst = Math.Max(worst, Worst(exact.Posterior.Draws[d].C, sampled.Posterior.Draws[d].C, sampled.MonteCarloError![d]));
        return worst;
    }

    private static double Worst(CoActiveAnalysis exact, CoActiveAnalysis sampled)
    {
        double worst = 0.0;
        for (int s = 0; s < exact.Posterior.Pairs.Count; s++)
        {
            CoActivePair e = exact.Posterior.Pairs[s], q = sampled.Posterior.Pairs[s];
            CoActiveError se = sampled.MonteCarloError![s];
            worst = Math.Max(worst, Math.Max(Worst(e.Cf, q.Cf, se.Cf), Math.Max(Worst(e.Cg, q.Cg, se.Cg), Worst(e.Cfg, q.Cfg, se.Cfg))));
        }
        return worst;
    }

    private static double Worst(double[,] exact, double[,] sampled, double[,] se)
    {
        double worst = 0.0;
        for (int i = 0; i < exact.GetLength(0); i++)
            for (int j = 0; j < exact.GetLength(1); j++)
            {
                double gap = Math.Abs(exact[i, j] - sampled[i, j]);
                if (se[i, j] > 0.0) worst = Math.Max(worst, gap / se[i, j]);
                else Assert.True(gap <= 1e-10 * Math.Max(1.0, Math.Abs(exact[i, j])), $"a zero error with gap {gap:E2}");
            }
        return worst;
    }

    // The weighted covariance about the mean, the sample as its own law, two passes.
    internal static double[,] WeightedCovariance(double[][] x, double[]? w)
    {
        int p = x.Length, n = x[0].Length;
        double total = w?.Sum() ?? n;
        var mean = new double[p];
        for (int j = 0; j < p; j++)
            for (int r = 0; r < n; r++) mean[j] += (w?[r] ?? 1.0) / total * x[j][r];
        var c = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
                for (int r = 0; r < n; r++) c[i, j] += (w?[r] ?? 1.0) / total * (x[i][r] - mean[i]) * (x[j][r] - mean[j]);
        return c;
    }

    // Σ^{1/2} by Denman–Beavers iteration, independent of an eigendecomposition: Y → Σ^{1/2}, Z → Σ^{−1/2}.
    internal static double[,] DenmanBeavers(double[,] sigma)
    {
        int p = sigma.GetLength(0);
        double[,] y = (double[,])sigma.Clone(), z = Identity(p);
        for (int k = 0; k < 60; k++)
        {
            double[,] yi = Inverse(y), zi = Inverse(z);
            var ny = new double[p, p];
            var nz = new double[p, p];
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                {
                    ny[i, j] = 0.5 * (y[i, j] + zi[i, j]);
                    nz[i, j] = 0.5 * (z[i, j] + yi[i, j]);
                }
            y = ny;
            z = nz;
        }
        for (int i = 0; i < p; i++)
            for (int j = i + 1; j < p; j++) y[i, j] = y[j, i] = 0.5 * (y[i, j] + y[j, i]);
        return y;
    }

    private static double[,] Identity(int p)
    {
        var m = new double[p, p];
        for (int i = 0; i < p; i++) m[i, i] = 1.0;
        return m;
    }

    // Gauss–Jordan with partial pivoting.
    private static double[,] Inverse(double[,] a)
    {
        int p = a.GetLength(0);
        var m = new double[p, 2 * p];
        for (int i = 0; i < p; i++)
        {
            for (int j = 0; j < p; j++) m[i, j] = a[i, j];
            m[i, p + i] = 1.0;
        }
        for (int c = 0; c < p; c++)
        {
            int pivot = c;
            for (int r = c + 1; r < p; r++)
                if (Math.Abs(m[r, c]) > Math.Abs(m[pivot, c])) pivot = r;
            for (int j = 0; j < 2 * p; j++) (m[c, j], m[pivot, j]) = (m[pivot, j], m[c, j]);
            double d = m[c, c];
            for (int j = 0; j < 2 * p; j++) m[c, j] /= d;
            for (int r = 0; r < p; r++)
            {
                if (r == c) continue;
                double factor = m[r, c];
                if (factor == 0.0) continue;
                for (int j = 0; j < 2 * p; j++) m[r, j] -= factor * m[c, j];
            }
        }
        var inv = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) inv[i, j] = m[i, p + j];
        return inv;
    }

    private static double[,] Multiply(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), k = a.GetLength(1), p = b.GetLength(1);
        var c = new double[n, p];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++)
            {
                double s = 0.0;
                for (int q = 0; q < k; q++) s += a[i, q] * b[q, j];
                c[i, j] = s;
            }
        return c;
    }

    private static double[,] Transpose(double[,] a)
    {
        var t = new double[a.GetLength(1), a.GetLength(0)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) t[j, i] = a[i, j];
        return t;
    }

    private static double[,] Congruence(double[,] c, double[,] l) => Multiply(Multiply(Transpose(l), c), l);

    // A random matrix of the given scale with singular values within a factor of about 30.
    private static double[,] RandomMatrix(Xoshiro256PlusPlus rng, int p, double scale)
    {
        double[,] u = Orthogonal(rng, p), v = Orthogonal(rng, p);
        var m = new double[p, p];
        for (int k = 0; k < p; k++)
        {
            double s = scale * Math.Pow(30.0, -rng.NextDouble());
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) m[i, j] += s * u[i, k] * v[j, k];
        }
        return m;
    }

    // A random orthogonal matrix, Gram–Schmidt on normal columns.
    private static double[,] Orthogonal(Xoshiro256PlusPlus rng, int p)
    {
        var q = new double[p, p];
        for (int k = 0; k < p; k++)
        {
            var v = new double[p];
            for (int i = 0; i < p; i++) v[i] = Normal.SampleStandard(rng);
            for (int pass = 0; pass < 2; pass++)
                for (int j = 0; j < k; j++)
                {
                    double d = 0.0;
                    for (int i = 0; i < p; i++) d += v[i] * q[i, j];
                    for (int i = 0; i < p; i++) v[i] -= d * q[i, j];
                }
            double n = Math.Sqrt(v.Sum(x => x * x));
            for (int i = 0; i < p; i++) q[i, k] = v[i] / n;
        }
        return q;
    }

    // Eigenvalues, descending.
    private static double[] Eigenvalues(double[,] c) => [.. DenseEigen.DecomposeSymmetric(c).Eigenvalues.OrderByDescending(x => x)];

    // Eigenvectors as columns, by descending eigenvalue.
    private static double[,] Vectors(double[,] c)
    {
        EigenResult e = DenseEigen.DecomposeSymmetric(c);
        int p = c.GetLength(0);
        int[] order = [.. Enumerable.Range(0, p).OrderByDescending(i => e.Eigenvalues[i])];
        var v = new double[p, p];
        for (int k = 0; k < p; k++)
            for (int i = 0; i < p; i++) v[i, k] = e.Eigenvectors[order[k]][i];
        return v;
    }

    private static double[] Lead(double[,] vectors) => [vectors[0, 0], vectors[1, 0]];

    // A metric direction v read as a function of the box coordinates: u = L⁻¹(x − μ), so vᵀu = (L⁻ᵀv)ᵀ(x − μ); unit length.
    private static double[] Functional(double[,] l, double[] v)
    {
        double[,] li = Inverse(l);
        double[] g = [li[0, 0] * v[0] + li[1, 0] * v[1], li[0, 1] * v[0] + li[1, 1] * v[1]];
        double n = Math.Sqrt(g[0] * g[0] + g[1] * g[1]);
        return [g[0] / n, g[1] / n];
    }

    private static double Angle(double[] a, double[] b)
        => Math.Acos(Math.Min(1.0, Math.Abs(a[0] * b[0] + a[1] * b[1]) / Math.Sqrt((a[0] * a[0] + a[1] * a[1]) * (b[0] * b[0] + b[1] * b[1]))));

    private static double Frobenius(double[,] a, double[,] b)
    {
        double s = 0.0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) s += (a[i, j] - b[i, j]) * (a[i, j] - b[i, j]);
        return Math.Sqrt(s);
    }
}
