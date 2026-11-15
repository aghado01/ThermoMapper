using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Sensitivity;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The sampled-sensitivity brief's slice 1 gates: Sobol indices on the sampled route for a scalar fit (<see cref="BassSensitivity.Of"/>
/// with <see cref="SensitivityRoute"/>), under both main-variance forms while the form's default waits on the user.
/// <list type="bullet">
/// <item>Separable fits read on both routes: every draw's sampled pieces against its exact ones, replicated over seeds, within their
/// stated errors, on the sensitivity fixtures under every kind of product law and on a Legendre fit.</item>
/// <item>A hand-built directional fit against tensor quadrature of its surface, a method apart from the estimator.</item>
/// <item>The route's resolution and refusals; the exact route bit for bit under Auto; the points shared by every draw and drawn from the
/// seed; the summary's error rows; and the Poincaré check on a directional fit, its totals sampled here and its C E8's.</item>
/// </list>
/// </summary>
public sealed class SampledSensitivityTests(ITestOutputHelper output)
{
    private static readonly Dictionary<string, InputMeasure> MixedLaws = new()
    {
        ["x0"] = new MixtureMeasure([new NormalComponent(0.5, 0.35, 0.15), new NormalComponent(0.5, 0.75, 0.25)]),
        ["x1"] = new UniformMeasure(0.1, 0.9),
        ["z"] = new CategoricalMeasure(new Dictionary<string, double> { ["a"] = 0.2, ["b"] = 0.5, ["c"] = 0.3 }),
    };

    private static readonly Dictionary<string, InputMeasure> EmpiricalAndStudent = new()
    {
        ["x0"] = new EmpiricalMeasure(),
        ["x1"] = new MixtureMeasure([new StudentComponent(1.0, 4.0, 0.5, 0.2)]),
    };

    private static readonly Lazy<BassFit> Legendre = new(() => Bass.Fit(KhaosReadoutTests.Config(5), KhaosReadoutTests.Data(120, 20261110)));

    // The separable references: the fit, its data and the spec's laws.
    private static (BassFit Fit, BassData Data, Dictionary<string, InputMeasure>? Laws) Separable(string name) => name switch
    {
        "hinges, mixed laws" => (SensitivityFixture.FitScalar(), SensitivityFixture.Scalar(), MixedLaws),
        "hinges, the boxes" => (SensitivityFixture.FitScalar(), SensitivityFixture.Scalar(), null),
        "degree 2, empirical and Student" => (SensitivityFixture.FitScalar(2), SensitivityFixture.Scalar(), EmpiricalAndStudent),
        "Legendre" => (Legendre.Value, KhaosReadoutTests.Data(120, 20261110), null),
        _ => throw new ArgumentException(name),
    };

    public static TheoryData<string, MainVarianceForm> SeparableCases()
    {
        var cases = new TheoryData<string, MainVarianceForm>();
        foreach (string name in new[] { "hinges, mixed laws", "hinges, the boxes", "degree 2, empirical and Student", "Legendre" })
            foreach (MainVarianceForm form in new[] { MainVarianceForm.Saltelli, MainVarianceForm.TwoDifference }) cases.Add(name, form);
        return cases;
    }

    internal sealed record Piece(string Name, double Truth, double Estimate, double Error);

    // Every piece of draw d: the variance, each subset's main variance and index, each input's total variance and index, against the
    // truth's decomposition (the exact route's, or quadrature's).
    internal static List<Piece> Pieces(SobolDecomposition sampled, int d, Func<string, int[]?, double> truth)
    {
        var pieces = new List<Piece> { new("V", truth("V", null), sampled.Variance[d], sampled.VarianceError![d]) };
        for (int k = 0; k < sampled.Subsets.Count; k++)
        {
            int[] u = sampled.Subsets[k];
            string label = string.Join(",", u);
            pieces.Add(new($"V_{label}", truth("V_u", u), sampled.MainVariance[d, k], sampled.MainVarianceError![d, k]));
            if (sampled.IsDefined(d)) pieces.Add(new($"S_{label}", truth("S_u", u), sampled.MainIndex(d, k), sampled.MainIndexError![d, k]));
        }
        for (int i = 0; i < sampled.Inputs.Count; i++)
        {
            pieces.Add(new($"T_{i}", truth("T", [i]), sampled.TotalVariance[d, i], sampled.TotalVarianceError![d, i]));
            if (sampled.IsDefined(d)) pieces.Add(new($"ST_{i}", truth("ST", [i]), sampled.TotalIndex(d, i), sampled.TotalIndexError![d, i]));
        }
        return pieces;
    }

    internal sealed record Tally(int Cells, int Exact, double WorstBias, string WorstBiasCell, double PooledSd, double LowestSd, double HighestSd,
                                 double Beyond3, double WorstZ, string WorstZCell);

    // Over seeds, per cell (draw and piece): the mean z in units of its own error, √R·mean z, and the spread of z; pooled over cells, the
    // spread and the tail. A cell whose truth and every estimate are zero within rounding is an exact zero.
    internal static Tally Summarize(List<Dictionary<string, Piece>> runs, Func<string, double> scale)
    {
        int r = runs.Count, exact = 0;
        double worstBias = 0.0, low = double.PositiveInfinity, high = 0.0, worstZ = 0.0;
        string worstCell = "", worstZCell = "";
        var all = new List<double>();
        foreach (string cell in runs[0].Keys)
        {
            Piece[] seen = [.. runs.Select(run => run[cell])];
            double s = scale(cell);
            if (Math.Abs(seen[0].Truth) <= 1e-12 * s && seen.All(p => Math.Abs(p.Estimate) <= 1e-10 * s))
            {
                exact++;
                continue;
            }
            Assert.All(seen, p => Assert.True(p.Error > 0.0, $"{cell}: estimate {p.Estimate:R}, truth {p.Truth:R}, no error"));
            double[] z = [.. seen.Select(p => (p.Estimate - p.Truth) / p.Error)];
            all.AddRange(z);
            double mean = z.Average(), sd = Math.Sqrt(z.Sum(x => (x - mean) * (x - mean)) / (r - 1));
            double bias = Math.Abs(mean) * Math.Sqrt(r);
            if (bias > worstBias) (worstBias, worstCell) = (bias, cell);
            low = Math.Min(low, sd);
            high = Math.Max(high, sd);
            double peak = z.MaxBy(Math.Abs);
            if (Math.Abs(peak) > Math.Abs(worstZ)) (worstZ, worstZCell) = (peak, cell);
        }
        double m = all.Average(), pooled = Math.Sqrt(all.Sum(x => (x - m) * (x - m)) / (all.Count - 1));
        return new Tally(runs[0].Count - exact, exact, worstBias, worstCell, pooled, low, high, all.Count(x => Math.Abs(x) > 3.0) / (double)all.Count,
                         worstZ, worstZCell);
    }

    // The cut for √R·|mean z| over `cells` cells, each near a standard normal when the estimate is unbiased: Bonferroni's at a
    // family-wise level of 10⁻³, so the cut grows with the cells read (4.2 at 42 cells, 5.2 at 6,000).
    internal static double BiasCut(int cells) => Maths.Distributions.Normal.StandardQuantile(1.0 - 0.0005 / cells);

    // The pooled spread of z over a fit's cells: the cells of one seed share its points, so the pooled spread has about the seeds' few
    // degrees of freedom. Slices 1 and 2's eight separable tables read 0.909 to 1.082, which makes ±0.1 about a two-sigma band; the
    // band here is ±0.15. The primitives' gates (slice 0, 300 replicates) read calibration sharply, and these read it against gross error.
    internal const double PooledLow = 0.85, PooledHigh = 1.15;

    private void Report(string name, Tally t) =>
        output.WriteLine($"{name}: {t.Cells} cells, {t.Exact} exact zeros; worst √R·|mean z| {t.WorstBias:F2} ({t.WorstBiasCell}); pooled sd of z "
                         + $"{t.PooledSd:F3}, per cell {t.LowestSd:F2} to {t.HighestSd:F2}; beyond 3: {t.Beyond3:P2}; worst z {t.WorstZ:F2} ({t.WorstZCell})");

    [Theory]
    [MemberData(nameof(SeparableCases))]
    public void SeparableFits_SampledMeetsExact_WithinTheStatedError(string name, MainVarianceForm form)
    {
        const int Seeds = 24, N = 2048;
        (BassFit fit, BassData data, Dictionary<string, InputMeasure>? laws) = Separable(name);
        SobolDecomposition exact = BassSensitivity.Of(fit, data, new SensitivitySpec { Measures = laws, Route = SensitivityRoute.Exact }).Decompositions[0];
        var runs = new List<Dictionary<string, Piece>>();
        foreach (int seed in SeedTree.Derive(20261115, Seeds))
        {
            SobolAnalysis analysis = BassSensitivity.OfInForm(fit, data, new SensitivitySpec { Measures = laws, Route = SensitivityRoute.Sampled, Points = N, Seed = seed }, form);
            Assert.Equal(SensitivityRoute.Sampled, analysis.Route);
            Assert.Equal(N, analysis.Points);
            SobolDecomposition sampled = analysis.Decompositions[0];
            Assert.Equal(exact.Subsets.Select(u => string.Join(",", u)), sampled.Subsets.Select(u => string.Join(",", u)));
            var cells = new Dictionary<string, Piece>();
            for (int d = 0; d < exact.Count; d++)
            {
                Assert.Equal(exact.IsDefined(d), sampled.IsDefined(d));
                int draw = d;
                double Truth(string piece, int[]? u) => piece switch
                {
                    "V" => exact.Variance[draw],
                    "V_u" => exact.MainVariance[draw, IndexOf(exact, u!)],
                    "S_u" => exact.MainIndex(draw, IndexOf(exact, u!)),
                    "T" => exact.TotalVariance[draw, u![0]],
                    _ => exact.TotalIndex(draw, u![0]),
                };
                foreach (Piece p in Pieces(sampled, d, Truth)) cells[$"draw {d} {p.Name}"] = p;
            }
            runs.Add(cells);
        }
        Tally t = Summarize(runs, cell => cell.Split(' ')[2].StartsWith('S') ? 1.0 : exact.Variance[int.Parse(cell.Split(' ')[1])]);
        Report($"{name}, {form}", t);
        Assert.True(t.WorstBias <= BiasCut(t.Cells), $"{t.WorstBiasCell}: √R·|mean z| {t.WorstBias:F2} against {BiasCut(t.Cells):F2}");
        Assert.True(t.PooledSd >= PooledLow && t.PooledSd <= PooledHigh, $"pooled sd of z {t.PooledSd:F3}");
        Assert.True(t.Beyond3 <= 0.01, $"{t.Beyond3:P2} beyond 3");
    }

    // ── A directional fit by hand ───────────────────────────────────────────────────────────────────────────────────────────

    private static readonly ProjectionDirection Dense = new([0, 1, 2], [0.6, -0.48, 0.64]);
    private static readonly ProjectionDirection Single = new([1], [1.0]);
    private static readonly ProjectionDirection Pair = new([2, 0], [-0.6, 0.8]);

    // DirectionalHingeTests' data: x0 on [0, 1], x1 on [−2, 5] skewed toward −2, x2 on [10, 13] inside a declared domain [9, 14].
    internal static BassData DirectionalData()
    {
        const int n = 80;
        var x0 = new double[n];
        var x1 = new double[n];
        var x2 = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            double u = i * 0.6180339887498949 % 1.0, v = i * 0.4142135623730950 % 1.0, w = i * 0.7320508075688772 % 1.0;
            x0[i] = u;
            x1[i] = -2.0 + 7.0 * v * v;
            x2[i] = 10.0 + 3.0 * w;
            y[i] = Math.Sin(2.0 * (0.6 * x0[i] + 0.1 * x1[i])) + 0.3 * x2[i] + 0.01 * Math.Cos(17.0 * i);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1), new ContinuousPredictor("x2", x2, new InputDomain(9.0, 14.0))],
                            new ScalarResponse(y));
    }

    // Two hand-built draws, their term ids ascending: (Dense, −0.5) and (Pair, 0.3); then (Pair, 0.3), (Single, 0.0) and (Dense, 0.7).
    private static readonly (ProjectionDirection Direction, double Knot)[] Hinges = [(Dense, -0.5), (Pair, 0.3), (Single, 0.0), (Dense, 0.7)];
    private static readonly (int[] Hinges, double[] Beta)[] HandDraws = [([0, 1], [0.4, 1.5, -0.7]), ([1, 2, 3], [-0.2, 1.1, 0.9, -2.0])];

    internal static BassFit DirectionalFit(BassData data, int degree)
    {
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 6, Degree = degree },
            Sampling = new BassSampling { Chains = 1, MasterSeed = 36, BurnIn = 20, Samples = 2 },
        };
        BassFit run = Bass.Fit(config, data);
        var table = new BassTermTable();
        int[] ids = [.. Hinges.Select(h => table.Intern(new BassTerm(directional: new DirectionalHinge(h.Direction, h.Knot))))];
        var store = new DrawStore<BassDraw>();
        foreach ((int[] hinges, double[] beta) in HandDraws) store.Add(new BassDraw(new BassStructure([.. hinges.Select(h => ids[h])]), 1.0, 1.0, beta, 1.0));
        return run with { Terms = table, Draws = [store], DrawWeights = null };
    }

    // The hand surface of a draw in raw units, DirectionalHingeTests.Hand's arithmetic: z by the training mean and sample deviation, each
    // hinge over its normalizer c.
    internal static Func<double[], double> HandSurface(BassData data, int degree, int draw)
    {
        (double Mean, double Deviation)[] moments = DirectionalHingeTests.Moments(data);
        double[][] probe = [[0.5], [0.0], [11.0]];
        (int[] hinges, double[] beta) = HandDraws[draw];
        double[] c = [.. hinges.Select(h => DirectionalHingeTests.Hand(data, probe, Hinges[h].Direction, Hinges[h].Knot, degree, HingeNormalization.Box).C)];
        return x =>
        {
            double f = beta[0];
            for (int m = 0; m < hinges.Length; m++)
            {
                (ProjectionDirection d, double knot) = Hinges[hinges[m]];
                double p = 0.0;
                for (int k = 0; k < d.Count; k++) p += d.Weights[k] * (x[d.Variables[k]] - moments[d.Variables[k]].Mean) / moments[d.Variables[k]].Deviation;
                f += beta[m + 1] * Math.Pow(Math.Max(p - knot, 0.0), degree) / c[m];
            }
            return f;
        };
    }

    // The boxes the default law is uniform on: the declared domain, else the training range.
    private static ReferenceBox Boxes(BassFit fit) => new([.. fit.Transform.Predictors.Select(p => ((RangeScale)p).Min)], [.. fit.Transform.Predictors.Select(p => ((RangeScale)p).Max)]);

    public static TheoryData<int, MainVarianceForm> DirectionalCases()
    {
        var cases = new TheoryData<int, MainVarianceForm>();
        foreach (int degree in new[] { 1, 2 })
            foreach (MainVarianceForm form in new[] { MainVarianceForm.Saltelli, MainVarianceForm.TwoDifference }) cases.Add(degree, form);
        return cases;
    }

    [Theory]
    [MemberData(nameof(DirectionalCases))]
    public void DirectionalFit_SampledMeetsQuadrature_WithinTheStatedError(int degree, MainVarianceForm form)
    {
        const int Seeds = 60, N = 4096;
        BassData data = DirectionalData();
        BassFit fit = DirectionalFit(data, degree);
        ReferenceBox box = Boxes(fit);
        QuadratureReference[] truth = [.. Enumerable.Range(0, HandDraws.Length).Select(d => new QuadratureReference($"draw {d}", HandSurface(data, degree, d), box, 48))];

        // The truth's own error: the same quadrature on half the panels.
        double coarse = 0.0;
        for (int d = 0; d < truth.Length; d++)
        {
            var half = new QuadratureReference($"draw {d}", HandSurface(data, degree, d), box, 24);
            foreach (int[] u in AllSubsets(3)) coarse = Math.Max(coarse, Math.Abs(half.MainVariance(u) - truth[d].MainVariance(u)) / truth[d].Variance);
        }

        var runs = new List<Dictionary<string, Piece>>();
        foreach (int seed in SeedTree.Derive(36, Seeds))
        {
            SobolAnalysis analysis = BassSensitivity.OfInForm(fit, data, new SensitivitySpec { Points = N, Seed = seed }, form);
            Assert.Equal(SensitivityRoute.Sampled, analysis.Route);
            SobolDecomposition sampled = analysis.Decompositions[0];
            Assert.Equal(7, sampled.Subsets.Count);
            var cells = new Dictionary<string, Piece>();
            for (int d = 0; d < sampled.Count; d++)
            {
                QuadratureReference q = truth[d];
                double Truth(string piece, int[]? u) => piece switch
                {
                    "V" => q.Variance,
                    "V_u" => q.MainVariance(u!),
                    "S_u" => q.MainVariance(u!) / q.Variance,
                    "T" => q.TotalVariance(u![0]),
                    _ => q.TotalVariance(u![0]) / q.Variance,
                };
                foreach (Piece p in Pieces(sampled, d, Truth)) cells[$"draw {d} {p.Name}"] = p;
            }
            runs.Add(cells);
        }
        Tally t = Summarize(runs, cell => cell.Split(' ')[2].StartsWith('S') ? 1.0 : truth[int.Parse(cell.Split(' ')[1])].Variance);
        Report($"degree {degree}, {form}", t);
        output.WriteLine($"degree {degree}: the quadrature's main variances move by at most {coarse:E2} of V on halving its panels");
        Assert.True(coarse <= 1e-3, $"quadrature at {coarse:E2}");
        Assert.True(t.WorstBias <= BiasCut(t.Cells), $"{t.WorstBiasCell}: √R·|mean z| {t.WorstBias:F2} against {BiasCut(t.Cells):F2}");
        Assert.True(t.PooledSd >= PooledLow && t.PooledSd <= PooledHigh, $"pooled sd of z {t.PooledSd:F3}");
        Assert.True(t.LowestSd >= 0.7 && t.HighestSd <= 1.35, $"per-cell sd of z from {t.LowestSd:F2} to {t.HighestSd:F2}");
        Assert.True(t.Beyond3 <= 0.01, $"{t.Beyond3:P2} beyond 3");
    }

    // ── The route, the points and the summary ────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Auto_IsTheExactRoute_BitForBit_OnASeparableFit()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        SobolAnalysis auto = BassSensitivity.Of(fit, data, new SensitivitySpec { Measures = MixedLaws });
        SobolAnalysis exact = BassSensitivity.Of(fit, data, new SensitivitySpec { Measures = MixedLaws, Route = SensitivityRoute.Exact });
        Assert.Equal(SensitivityRoute.Exact, auto.Route);
        Assert.Equal(0, auto.Points);
        SobolDecomposition a = auto.Decompositions[0], e = exact.Decompositions[0];
        Assert.Equal(e.MainVariance, a.MainVariance);
        Assert.Equal(e.TotalVariance, a.TotalVariance);
        Assert.Equal(e.Variance, a.Variance);
        Assert.Equal(e.Status, a.Status);
        Assert.Null(a.VarianceError);
        Assert.Null(a.MainIndexError);
    }

    [Fact]
    public void Route_RefusesWithItsReason()
    {
        BassData data = DirectionalData();
        BassFit fit = DirectionalFit(data, 1);
        ArgumentException exact = Assert.Throws<ArgumentException>(() => fit.Sensitivity(data, new SensitivitySpec { Route = SensitivityRoute.Exact }));
        Assert.Contains("not a product of one-variable factors", exact.Message);
        Assert.Contains(fit.Terms[fit.Draws[0][0].Structure.TermIds[0]].ToString(), exact.Message);   // it names the term it meets
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.Sensitivity(data, new SensitivitySpec { Points = 1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.Sensitivity(data, new SensitivitySpec { Route = (SensitivityRoute)7 }));

        BassFit functional = SensitivityFixture.FitFunctional();
        ArgumentException later = Assert.Throws<ArgumentException>(() => BassSensitivity.Of(functional, SensitivityFixture.Functional(), new SensitivitySpec { Route = SensitivityRoute.Sampled }));
        Assert.Contains("slice 3", later.Message);
    }

    [Fact]
    public void Points_AreDrawnFromTheSeed_AndSharedByEveryDraw()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        var spec = new SensitivitySpec { Measures = MixedLaws, Route = SensitivityRoute.Sampled, Points = 512, Seed = 3 };
        SobolDecomposition all = BassSensitivity.Of(fit, data, spec).Decompositions[0];
        SobolDecomposition again = BassSensitivity.Of(fit, data, spec).Decompositions[0];
        Assert.Equal(all.MainVariance, again.MainVariance);
        Assert.Equal(all.TotalIndexError, again.TotalIndexError);
        SobolDecomposition other = BassSensitivity.Of(fit, data, spec with { Seed = 4 }).Decompositions[0];
        Assert.NotEqual(all.Variance, other.Variance);

        // A draw read alone reads the same points: its row, bit for bit, whatever else is read beside it.
        foreach (int d in new[] { 0, 17, all.Count - 1 })
        {
            SobolDecomposition one = BassSensitivity.Of(fit, data, spec with { Draws = [d] }).Decompositions[0];
            Assert.Equal(all.Variance[d], one.Variance[0]);
            Assert.Equal(all.VarianceError![d], one.VarianceError![0]);
            for (int i = 0; i < all.Inputs.Count; i++)
            {
                Assert.Equal(all.TotalVariance[d, i], one.TotalVariance[0, i]);
                Assert.Equal(all.TotalIndexError![d, i], one.TotalIndexError![0, i]);
            }
        }
    }

    [Fact]
    public void Summary_CarriesTheErrorRows()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        var spec = new SensitivitySpec { Measures = MixedLaws, Route = SensitivityRoute.Sampled, Points = 512 };
        SobolAnalysisSummary summary = fit.SensitivitySummary(data, spec);
        SobolDecomposition s = fit.Sensitivity(data, spec).Decompositions[0];
        Assert.Equal(SensitivityRoute.Sampled, summary.Route);
        Assert.Equal(512, summary.Points);
        SobolSummary one = summary.Decompositions[0];
        int[] defined = [.. Enumerable.Range(0, s.Count).Where(s.IsDefined)];
        for (int k = 0; k < s.Subsets.Count; k++)
            Assert.Equal(defined.Average(d => s.MainIndexError![d, k]), one.MainIndexError![k], 12);
        for (int i = 0; i < s.Inputs.Count; i++)
            Assert.Equal(Enumerable.Range(0, s.Count).Average(d => s.TotalVarianceError![d, i]), one.TotalVarianceError![i], 12);
        Assert.Equal(s.VarianceError!.Average(), one.VarianceError!.Value, 12);
        Assert.Equal(s.Status.Count(x => (x & DecompositionStatus.BelowZeroBeyondError) != 0), one.BelowZeroBeyondError);

        SobolAnalysisSummary exact = fit.SensitivitySummary(data, spec with { Route = SensitivityRoute.Exact });
        Assert.Equal(SensitivityRoute.Exact, exact.Route);
        Assert.Null(exact.Decompositions[0].MainIndexError);
        Assert.Null(exact.Decompositions[0].VarianceError);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void Poincare_OnADirectionalFit_TotalsSampledAgainstEightsC(int degree)
    {
        // Decision 7, read here before E8's readout takes it: each draw's sampled total(v) against C_P(v)·C_vv, C from E8's sampled route in
        // box units and C_P the input's constant, within three standard errors of the two independent estimates.
        BassData data = DirectionalData();
        BassFit fit = DirectionalFit(data, degree);
        SobolAnalysis sobol = fit.Sensitivity(data, new SensitivitySpec { Points = 8192, Seed = 5 });
        ActiveSubspaceAnalysis c = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Points = 8192, Seed = 6 });
        SobolDecomposition s = sobol.Decompositions[0];
        Assert.Equal(s.Inputs, c.Inputs);
        double worst = 0.0, worstZ = double.NegativeInfinity;
        for (int i = 0; i < s.Inputs.Count; i++)
        {
            var scale = (RangeScale)fit.Transform.Predictors[i];
            InputPoincareConstant constant = InputPoincareConstant.Of(sobol.Measure, i, scale.Max - scale.Min);
            Assert.Equal(PoincareConstantKind.Exact, constant.Kind);
            for (int d = 0; d < s.Count; d++)
            {
                double bound = constant.Value * c.Posterior.Draws[d].C[i, i];
                double error = Math.Sqrt(Math.Pow(s.TotalVarianceError![d, i], 2) + Math.Pow(constant.Value * c.MonteCarloError![d][i, i], 2));
                double total = s.TotalVariance[d, i];
                if (bound > 0.0) worst = Math.Max(worst, total / bound);
                worstZ = Math.Max(worstZ, (total - bound) / Math.Max(error, double.Epsilon));
                Assert.True(total <= bound + 3.0 * error, $"draw {d}, {s.Inputs[i]}: total {total:G6} against the bound {bound:G6} ± {error:G3}");
            }
        }
        output.WriteLine($"degree {degree}: largest total / bound {worst:F3}; largest (total − bound) / error {worstZ:F2}");
    }

    private static int IndexOf(SobolDecomposition s, int[] u)
    {
        for (int k = 0; k < s.Subsets.Count; k++)
            if (s.Subsets[k].SequenceEqual(u)) return k;
        throw new ArgumentException($"No subset {string.Join(",", u)}.");
    }

    private static IEnumerable<int[]> AllSubsets(int p)
    {
        for (int mask = 1; mask < 1 << p; mask++) yield return [.. Enumerable.Range(0, p).Where(k => (mask >> k & 1) != 0)];
    }
}

/// <summary>
/// A function on a box whose ANOVA pieces come from tensor Gauss–Legendre quadrature (<see cref="SobolQuadrature"/>) on equal panels,
/// a method apart from the estimators: for a surface whose kinks no panel follows, such as a directional hinge's, the error falls with the
/// panels' width and is read by halving them.
/// </summary>
internal sealed class QuadratureReference : SobolReference
{
    private readonly Func<double[], double> _f;
    private readonly int _panels;
    private readonly Dictionary<string, double> _closed = [];
    private readonly IReadOnlyList<QuadratureRule> _rules;
    private readonly double[] _x;

    public QuadratureReference(string name, Func<double[], double> f, ReferenceBox law, int panels) : base(name, law)
    {
        _f = f;
        _panels = panels;
        _x = new double[law.Dimension];
        _rules = SobolQuadrature.Rules(this);
    }

    public override double Value(ReadOnlySpan<double> x)
    {
        x.CopyTo(_x);
        return _f(_x);
    }

    public override (int Nodes, int Panels) Quadrature(int v) => (3, _panels);

    private double? _mean;

    public override double Mean => _mean ??= SobolQuadrature.Mean(this, _rules);

    public override double ClosedVariance(int[] u)
    {
        string key = string.Join(",", u);
        if (!_closed.TryGetValue(key, out double tau)) _closed[key] = tau = SobolQuadrature.ClosedVariance(this, _rules, u);
        return tau;
    }

    public override double ConditionalMean(int[] u, ReadOnlySpan<double> values) => SobolQuadrature.ConditionalMean(this, _rules, u, values);
}
