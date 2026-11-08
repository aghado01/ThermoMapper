using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// KHAOS's readouts past the paper (slice 4), on Legendre fits to Ishigami's function under KHAOS's defaults (the modified
/// g-prior with g₀² drawn, the zero centre, Jeffreys σ²):
/// <list type="bullet">
/// <item><b>Sobol from the coefficients</b> (<see cref="BassSensitivity.FromCoefficients"/>) against the moment route
/// (<see cref="BassSensitivity.Of"/>), draw for draw, under the product uniform law on the boxes.</item>
/// <item><b>Lean SG-4(b)'s identities</b> against tensor Gauss–Legendre integration of each draw's own surface and gradient, the
/// rule exact for its polynomial degrees: Var f, each main variance, each total as E[(f − E_v f)²], the Jacobi measure
/// E[x_v(1 − x_v)(∂_v f)²], the bound total ≤ ½·Jacobi, and its equality on a fit of degree 1.</item>
/// <item><b>The gradient readout</b> (<see cref="BassFit.Gradient"/>) against central differences of the predicted surface, on a
/// Legendre fit whose declared domain is not [0, 1], so the chain rule through the box shows, and on a hinge fit.</item>
/// <item><b>Moments under other measures:</b> the moment route's decomposition of a Legendre fit under a truncated mixture
/// (normal and Student components, both weightings) and under the empirical marginals, against brute-force tensor sums over each
/// measure with its density written out here; and the coefficient route's refusals.</item>
/// </list>
/// </summary>
public sealed class KhaosReadoutTests(ITestOutputHelper output)
{
    // Ishigami's function on z = π(2x − 1), x ∈ [0, 1]: sin z₁ + 7 sin² z₂ + 0.1 z₃⁴ sin z₁.
    private static double Ishigami(double x0, double x1, double x2)
    {
        double z1 = Math.PI * (2 * x0 - 1), z2 = Math.PI * (2 * x1 - 1), z3 = Math.PI * (2 * x2 - 1);
        return Math.Sin(z1) + 7 * Math.Pow(Math.Sin(z2), 2) + 0.1 * Math.Pow(z3, 4) * Math.Sin(z1);
    }

    // n cases on [lo, hi]³ with that domain declared, y Ishigami's on the box mapped to [0, 1], plus 0.3 noise.
    internal static BassData Data(int n, int seed, double lo = 0.0, double hi = 1.0)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[3][];
        for (int j = 0; j < 3; j++) x[j] = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < 3; j++) x[j][i] = lo + (hi - lo) * rng.NextDouble();
            y[i] = Ishigami((x[0][i] - lo) / (hi - lo), (x[1][i] - lo) / (hi - lo), (x[2][i] - lo) / (hi - lo)) + 0.3 * Normal.SampleStandard(rng);
        }
        var domain = new InputDomain(lo, hi);
        return new BassData([.. Enumerable.Range(0, 3).Select(j => new ContinuousPredictor($"x{j}", x[j], domain))], new ScalarResponse(y));
    }

    internal static BassConfig Config(int degree, int order = 2, int seed = 5) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = order, Legendre = new LegendreBasis { MaxDegree = degree } },
        Count = new BassCountPrior(4, 0.05),
        Coefficients = BassCoefficientPrior.ModifiedG(TauPrior.Gamma(0.001, 1000.0)),
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.Jeffreys,
        Sampling = new BassSampling { Chains = 2, MasterSeed = seed, BurnIn = 1_000, Samples = 40, Thin = 10 },
    };

    // The n-point Gauss–Legendre rule mapped to [0, 1], its weights summing to 1.
    private static (double[] X, double[] W) Rule(int n)
    {
        GaussLegendre g = GaussLegendre.Of(n);
        double[] t = g.Nodes.ToArray(), w = g.Weights.ToArray();
        bool symmetric = t.Min() < 0.0;
        double sum = w.Sum();
        return ([.. t.Select(v => symmetric ? 0.5 * (v + 1.0) : v)], [.. w.Select(v => v / sum)]);
    }

    // Each selected draw's surface at the cases: draw × case.
    private static double[,] Surface(BassFit fit, BassData data, double[][] x, int[] draws)
    {
        BassPrediction p = fit.Predict(data, PredictionCases.Of([.. Enumerable.Range(0, x.Length).Select(j => new ContinuousPredictor($"x{j}", x[j]))]),
                                       new BassPredictOptions { Draws = draws, Predictive = false });
        var f = new double[draws.Length, x[0].Length];
        for (int d = 0; d < draws.Length; d++)
            for (int i = 0; i < x[0].Length; i++) f[d, i] = p.Latent[d, i, 0];
        return f;
    }

    private static BassGradient Gradient(BassFit fit, BassData data, double[][] x, int[] draws)
        => fit.Gradient(data, PredictionCases.Of([.. Enumerable.Range(0, x.Length).Select(j => new ContinuousPredictor($"x{j}", x[j]))]), draws);

    private static void Close(double expected, double actual, double tolerance, string what)
        => Assert.True(Math.Abs(expected - actual) <= tolerance * Math.Max(1.0, Math.Abs(expected)), $"{what}: {actual:R} against {expected:R}");

    // ---- Sobol from the coefficients against the moment route ----------------------------------------------------------------------

    [Fact]
    public void Coefficients_MeetTheMomentRoute_DrawForDraw()
    {
        BassData data = Data(120, 20261108);
        BassFit fit = Bass.Fit(Config(8, order: 3), data);
        long t0 = System.Diagnostics.Stopwatch.GetTimestamp();
        SobolAnalysis moments = fit.Sensitivity(data);
        long t1 = System.Diagnostics.Stopwatch.GetTimestamp();
        PolynomialChaosSensitivity chaos = BassSensitivity.FromCoefficients(fit, data);
        long t2 = System.Diagnostics.Stopwatch.GetTimestamp();
        double Ms(long a0, long a1) => 1000.0 * (a1 - a0) / System.Diagnostics.Stopwatch.Frequency;
        SobolDecomposition a = chaos.Sobol.Decompositions[0], b = moments.Decompositions[0];
        Assert.Equal(b.Subsets.Select(u => string.Join(",", u)), a.Subsets.Select(u => string.Join(",", u)));
        Assert.Equal(moments.MaxOrder, chaos.Sobol.MaxOrder);
        Assert.Equal(b.Weights, a.Weights);
        Assert.All(chaos.Sobol.Measure.Declared, m => Assert.IsType<UniformMeasure>(m));
        double worst = 0.0;
        for (int d = 0; d < a.Variance.Length; d++)
        {
            worst = Math.Max(worst, Math.Abs(a.Variance[d] - b.Variance[d]) / b.Variance[d]);
            for (int u = 0; u < a.Subsets.Count; u++) worst = Math.Max(worst, Math.Abs(a.MainVariance[d, u] - b.MainVariance[d, u]) / b.Variance[d]);
            for (int v = 0; v < 3; v++) worst = Math.Max(worst, Math.Abs(a.TotalVariance[d, v] - b.TotalVariance[d, v]) / b.Variance[d]);
            Assert.Equal(b.Status[d], a.Status[d]);
        }
        output.WriteLine($"{a.Variance.Length} draws, M up to {fit.Draws.SelectMany(c => c).Max(x => x.Structure.TermIds.Length)}, {a.Subsets.Count} subsets: " +
                         $"largest gap to the moment route {worst:E2} of Var f; wall time, indicative: the moment route {Ms(t0, t1):F1} ms, the coefficients {Ms(t1, t2):F1} ms");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
    }

    // ---- SG-4(b)'s identities against numeric integration -----------------------------------------------------------------------

    // At degree 2 and order 2 many draws hold an input only to degree 1 in every term, where its total is half its Jacobi measure
    // (tensorTotal_eq_half); those are checked wherever they occur.
    [Theory]
    [InlineData(8, 3)]
    [InlineData(2, 2)]
    public void Identities_MeetGaussLegendreIntegration(int degree, int order)
    {
        BassData data = Data(120, 20261109);
        BassFit fit = Bass.Fit(Config(degree, order), data);
        PolynomialChaosSensitivity chaos = BassSensitivity.FromCoefficients(fit, data);
        SobolDecomposition s = chaos.Sobol.Decompositions[0];
        int[] draws = [.. Enumerable.Range(0, s.Variance.Length).Where(d => d % 10 == 0)];
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        int equalities = 0;

        // The tensor rule on [0, 1]³: 10 nodes are exact to degree 19 per variable, and f², x(1 − x)(∂f)² reach 16.
        (double[] r, double[] w) = Rule(10);
        int k = r.Length, points = k * k * k;
        var x = new double[3][];
        for (int j = 0; j < 3; j++) x[j] = new double[points];
        for (int i = 0; i < points; i++) (x[0][i], x[1][i], x[2][i]) = (r[i / (k * k)], r[i / k % k], r[i % k]);
        double[,] f = Surface(fit, data, x, draws);
        BassGradient g = Gradient(fit, data, x, draws);
        int At(int a, int b, int c) => (a * k + b) * k + c;

        double worst = 0.0, slack = double.PositiveInfinity;
        for (int q = 0; q < draws.Length; q++)
        {
            int d = draws[q];
            double mean = 0.0, variance = 0.0, second = 0.0;
            for (int a = 0; a < k; a++)
                for (int b = 0; b < k; b++)
                    for (int c = 0; c < k; c++) mean += w[a] * w[b] * w[c] * f[q, At(a, b, c)];
            for (int a = 0; a < k; a++)
                for (int b = 0; b < k; b++)
                    for (int c = 0; c < k; c++)
                    {
                        double ww = w[a] * w[b] * w[c], v = f[q, At(a, b, c)];
                        variance += ww * (v - mean) * (v - mean);   // about the mean, two passes, so a near-constant draw keeps its digits
                        second += ww * v * v;
                    }
            double scale = Math.Max(Math.Max(variance, 1e-12 * second), 1e-300);
            Close(pooled[d].Coefficients.Skip(1).Sum(c => c * c), s.Variance[d], 1e-12, $"draw {d}: the decomposition's draw order is the pooled one");
            worst = Math.Max(worst, Math.Abs(variance - s.Variance[d]) / scale);
            for (int v = 0; v < 3; v++)
            {
                // E[f | x_v] over the other two, its variance the main one's; the total, E over the others of Var(f | others) in x_v.
                double main = 0.0, total = 0.0, jacobi = 0.0;
                for (int i = 0; i < k; i++)
                {
                    double conditional = 0.0;
                    for (int j1 = 0; j1 < k; j1++)
                        for (int j2 = 0; j2 < k; j2++)
                            conditional += w[j1] * w[j2] * f[q, Index(v, i, j1, j2, k)];
                    main += w[i] * (conditional - mean) * (conditional - mean);
                }
                for (int j1 = 0; j1 < k; j1++)
                    for (int j2 = 0; j2 < k; j2++)
                    {
                        double m1 = 0.0, spread = 0.0;
                        for (int i = 0; i < k; i++)
                        {
                            m1 += w[i] * f[q, Index(v, i, j1, j2, k)];
                            double dv = g.Partials[v][q, Index(v, i, j1, j2, k), 0];
                            jacobi += w[i] * w[j1] * w[j2] * r[i] * (1.0 - r[i]) * dv * dv;
                        }
                        for (int i = 0; i < k; i++) spread += w[i] * Math.Pow(f[q, Index(v, i, j1, j2, k)] - m1, 2);
                        total += w[j1] * w[j2] * spread;
                    }
                int single = Enumerable.Range(0, s.Subsets.Count).First(u => s.Subsets[u].Length == 1 && s.Subsets[u][0] == v);
                worst = Math.Max(worst, Math.Abs(main - s.MainVariance[d, single]) / scale);
                worst = Math.Max(worst, Math.Abs(total - s.TotalVariance[d, v]) / scale);
                worst = Math.Max(worst, Math.Abs(jacobi - chaos.Jacobi[d, v]) / Math.Max(chaos.Jacobi[d, v], scale));
                // SG-4: the total is at most half the Jacobi measure, with equality at degree 1.
                Assert.True(s.TotalVariance[d, v] <= 0.5 * chaos.Jacobi[d, v] * (1 + 1e-12) + 1e-300, $"draw {d}, x{v}: total above half the Jacobi measure");
                BassTerm[] holding = [.. pooled[d].Structure.TermIds.Select(id => fit.Terms[id]).Where(t => t.Legendre.Any(l => l.Variable == v))];
                if (holding.Length > 0 && holding.All(t => t.Legendre.Single(l => l.Variable == v).Degree == 1))
                {
                    Close(0.5 * chaos.Jacobi[d, v], s.TotalVariance[d, v], 1e-12, $"draw {d}, x{v}: degree 1's equality");
                    equalities++;
                }
                else if (s.TotalVariance[d, v] > 0) slack = Math.Min(slack, 0.5 * chaos.Jacobi[d, v] / s.TotalVariance[d, v]);
            }
        }
        output.WriteLine($"degree {degree}, order {order}: {draws.Length} draws on a {k}³ Gauss–Legendre rule; largest gap of Var f, the main and total " +
                         $"variances and the Jacobi measure to the coefficients' {worst:E2}; {equalities} equalities at degree 1, elsewhere the least " +
                         $"½·Jacobi over total {slack:F3}");
        Assert.True(worst <= 1e-9, $"{worst:E2}");
        if (degree == 2) Assert.True(equalities > 0, "no draw held an input to degree 1 alone");
    }

    // The flat index of the point whose variable v sits at node i and the other two at j1, j2 (in ascending order).
    private static int Index(int v, int i, int j1, int j2, int k) => v switch
    {
        0 => (i * k + j1) * k + j2,
        1 => (j1 * k + i) * k + j2,
        _ => (j1 * k + j2) * k + i,
    };

    // ---- the gradient readout against central differences -------------------------------------------------------------------------

    [Theory]
    [InlineData("legendre")]
    [InlineData("hinge")]
    public void Gradient_MeetsCentralDifferences(string basis)
    {
        // A Legendre fit on [−π, π]³ with that domain declared, so ∂/∂z = ∂/∂x / 2π; a hinge fit of degree 1 on BassSweepTests' data.
        bool legendre = basis == "legendre";
        BassData data = legendre ? Data(120, 20261110, -Math.PI, Math.PI) : BassSweepTests.Data(large: true);
        BassConfig config = legendre ? Config(8, 3) : BassSweepTests.Config(maxBasis: 6) with { Sampling = new BassSampling { Chains = 2, MasterSeed = 3, BurnIn = 300, Samples = 20, Thin = 5 } };
        BassFit fit = Bass.Fit(config, data);
        string[] names = [.. data.Predictors.OfType<ContinuousPredictor>().Select(p => p.Name)];
        double[] lo = [.. data.Predictors.OfType<ContinuousPredictor>().Select(p => p.Values.ToArray().Min())];
        double[] hi = [.. data.Predictors.OfType<ContinuousPredictor>().Select(p => p.Values.ToArray().Max())];
        var rng = new Xoshiro256PlusPlus(17);
        const int cases = 40;
        double[][] x = [.. Enumerable.Range(0, names.Length).Select(j => Enumerable.Range(0, cases).Select(_ => lo[j] + (hi[j] - lo[j]) * (0.05 + 0.9 * rng.NextDouble())).ToArray())];
        int[] draws = [.. Enumerable.Range(0, fit.Draws.Sum(c => c.Count)).Where(d => d % 5 == 0)];
        PredictionCases Cases(double[][] at) => PredictionCases.Of([.. names.Select((name, j) => new ContinuousPredictor(name, at[j]))]);
        BassGradient g = fit.Gradient(data, Cases(x), draws);
        Assert.Equal(names, g.Inputs);
        double worst = 0.0;
        for (int v = 0; v < names.Length; v++)
        {
            double h = 1e-6 * (hi[v] - lo[v]);
            double[][] up = [.. x.Select((col, j) => j == v ? col.Select(t => t + h).ToArray() : col)];
            double[][] down = [.. x.Select((col, j) => j == v ? col.Select(t => t - h).ToArray() : col)];
            var options = new BassPredictOptions { Draws = draws, Predictive = false };
            BassPrediction pu = fit.Predict(data, Cases(up), options), pd = fit.Predict(data, Cases(down), options);
            for (int d = 0; d < draws.Length; d++)
                for (int i = 0; i < cases; i++)
                {
                    double difference = (pu.Latent[d, i, 0] - pd.Latent[d, i, 0]) / (2 * h), partial = g.Partials[v][d, i, 0];
                    worst = Math.Max(worst, Math.Abs(partial - difference) / Math.Max(1.0, Math.Abs(partial)));
                }
            Assert.Same(g.Partials[v], g.Of(names[v]));
        }
        output.WriteLine($"{basis}: {draws.Length} draws × {cases} cases × {names.Length} inputs, largest relative gap to central differences {worst:E2}");
        Assert.True(worst <= 1e-5, $"{worst:E2}");
    }

    [Fact]
    public void Coefficients_RefuseAHingeFit_WhoseCoefficientsAreNotOrthonormal()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassFit fit = Bass.Fit(BassSweepTests.Config(maxBasis: 4) with { Sampling = new BassSampling { Chains = 1, MasterSeed = 3, BurnIn = 100, Samples = 10 } }, data);
        Assert.Contains("not a Legendre product", Assert.Throws<ArgumentException>(() => BassSensitivity.FromCoefficients(fit, data)).Message);
    }

    // ---- moments under other measures ----------------------------------------------------------------------------------------------

    // A truncated mixture's density on [0, 1] in the scaled units, written from the components' own densities.
    private static Func<double, double> MixtureDensity(MixtureMeasure m)
    {
        (double[] r, double[] w) = Rule(200);
        Func<double, double>[] parts = [.. m.Components.Select<MixtureComponent, Func<double, double>>(c => c switch
        {
            NormalComponent n => t => Math.Exp(-0.5 * Math.Pow((t - n.Mean) / n.Sd, 2)) / (n.Sd * Math.Sqrt(2 * Math.PI)),
            StudentComponent s => t =>
            {
                double z = (t - s.Location) / s.Scale, nu = s.DegreesOfFreedom;
                return Math.Exp(SpecialFunctions.LogGamma(0.5 * (nu + 1)) - SpecialFunctions.LogGamma(0.5 * nu) - 0.5 * Math.Log(nu * Math.PI)
                                - 0.5 * (nu + 1) * Math.Log(1 + z * z / nu)) / s.Scale;
            },
            _ => throw new ArgumentException("A normal or Student component."),
        })];
        double[] masses = [.. parts.Select(p => r.Select((t, i) => w[i] * p(t)).Sum())];
        double[] weights = [.. m.Components.Select(c => c.Weight)];
        if (m.Weighting == MixtureWeighting.BeforeTruncation)
        {
            double total = weights.Select((wk, j) => wk * masses[j]).Sum();
            return t => parts.Select((p, j) => weights[j] * p(t)).Sum() / total;
        }
        double sum = weights.Sum();
        return t => parts.Select((p, j) => weights[j] * p(t) / masses[j]).Sum() / sum;
    }

    [Theory]
    [InlineData("mixture")]
    [InlineData("empirical")]
    public void MomentRoute_UnderAnotherMeasure_MeetsBruteForce(string measure)
    {
        // A small Legendre fit on [0, 1]³ (degree 4, order 2), its inputs' box [0, 1], so the scaled units are the raw ones.
        BassData data = Data(30, 20261111);
        BassFit fit = Bass.Fit(Config(4, 2), data);
        IReadOnlyDictionary<string, InputMeasure> measures = measure == "mixture"
            ? new Dictionary<string, InputMeasure>
            {
                ["x0"] = new MixtureMeasure([new NormalComponent(0.6, 0.3, 0.1), new StudentComponent(0.4, 4, 0.7, 0.15)]),
                ["x1"] = new MixtureMeasure([new NormalComponent(1.0, 0.5, 0.2)]) { Weighting = MixtureWeighting.ConditionalOnTruncation },
                ["x2"] = new UniformMeasure(),
            }
            : new Dictionary<string, InputMeasure> { ["x0"] = new EmpiricalMeasure(), ["x1"] = new EmpiricalMeasure(), ["x2"] = new EmpiricalMeasure() };
        var spec = new SensitivitySpec { Measures = measures };
        SobolAnalysis analysis = BassSensitivity.Of(fit, data, spec);
        SobolDecomposition s = analysis.Decompositions[0];
        for (int v = 0; v < 3; v++)   // the readout says the measure it holds under
            Assert.IsType(measure == "mixture" ? (v < 2 ? typeof(MixtureMeasure) : typeof(UniformMeasure)) : typeof(EmpiricalMeasure), analysis.Measure.Declared[v]);

        // Each input's nodes and weights under its measure: 40 Gauss–Legendre nodes weighted by the density, or the training values.
        (double[] Nodes, double[] Weights)[] rules = new (double[], double[])[3];
        for (int v = 0; v < 3; v++)
        {
            if (measure == "empirical")
            {
                double[] values = data.Predictors[v] is ContinuousPredictor c ? c.Values.ToArray() : [];
                rules[v] = (values, [.. values.Select(_ => 1.0 / values.Length)]);
                continue;
            }
            (double[] r, double[] w) = Rule(40);
            Func<double, double> density = measures[$"x{v}"] is MixtureMeasure mix ? MixtureDensity(mix) : _ => 1.0;
            double[] dw = [.. r.Select((t, i) => w[i] * density(t))];
            double sum = dw.Sum();
            rules[v] = (r, [.. dw.Select(t => t / sum)]);
        }
        int k0 = rules[0].Nodes.Length, k1 = rules[1].Nodes.Length, k2 = rules[2].Nodes.Length, points = k0 * k1 * k2;
        var x = new double[3][];
        for (int j = 0; j < 3; j++) x[j] = new double[points];
        for (int i = 0; i < points; i++)
            (x[0][i], x[1][i], x[2][i]) = (rules[0].Nodes[i / (k1 * k2)], rules[1].Nodes[i / k2 % k1], rules[2].Nodes[i % k2]);
        int[] draws = [.. Enumerable.Range(0, s.Variance.Length).Where(d => d % 16 == 0)];
        double[,] f = Surface(fit, data, x, draws);

        double worst = 0.0;
        for (int q = 0; q < draws.Length; q++)
        {
            int d = draws[q];
            double Weight(int i) => rules[0].Weights[i / (k1 * k2)] * rules[1].Weights[i / k2 % k1] * rules[2].Weights[i % k2];
            double mean = 0.0, second = 0.0;
            for (int i = 0; i < points; i++) { mean += Weight(i) * f[q, i]; second += Weight(i) * f[q, i] * f[q, i]; }
            double variance = second - mean * mean;
            worst = Math.Max(worst, Math.Abs(variance - s.Variance[d]) / variance);
            for (int v = 0; v < 3; v++)
            {
                // Var(E[f | x_v]) and E[Var(f | x_−v)], by sums over the tensor grid.
                int kv = rules[v].Nodes.Length;
                var conditional = new double[kv];
                var byOthers = new Dictionary<int, (double M1, double M2, double W)>();
                for (int i = 0; i < points; i++)
                {
                    int iv = v == 0 ? i / (k1 * k2) : v == 1 ? i / k2 % k1 : i % k2;
                    double wo = Weight(i) / rules[v].Weights[iv];
                    conditional[iv] += wo * f[q, i];
                    int others = v == 0 ? i % (k1 * k2) : v == 1 ? (i / (k1 * k2)) * k2 + i % k2 : i / k2;
                    (double m1, double m2, double wt) = byOthers.GetValueOrDefault(others);
                    byOthers[others] = (m1 + rules[v].Weights[iv] * f[q, i], m2 + rules[v].Weights[iv] * f[q, i] * f[q, i], wo);
                }
                double main = conditional.Select((c, iv) => rules[v].Weights[iv] * (c - mean) * (c - mean)).Sum();
                double total = byOthers.Values.Sum(e => e.W * (e.M2 - e.M1 * e.M1));
                int single = Enumerable.Range(0, s.Subsets.Count).First(u => s.Subsets[u].Length == 1 && s.Subsets[u][0] == v);
                worst = Math.Max(worst, Math.Abs(main - s.MainVariance[d, single]) / variance);
                worst = Math.Max(worst, Math.Abs(total - s.TotalVariance[d, v]) / variance);
            }
        }
        output.WriteLine($"{measure}: {draws.Length} draws over {points:N0} grid points; largest gap of Var f, the main and total variances to the " +
                         $"moment route {worst:E2} of Var f; the moments' own quadrature error {analysis.MomentError:E2}");
        Assert.True(worst <= 1e-8, $"{worst:E2}");
    }
}
