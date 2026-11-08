using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// The engine's Legendre fits under KHAOS's raw ridge (the KHAOS brief's slice 1) against slice 0's corrected R references:
/// <list type="bullet">
/// <item><b>The exact law.</b> On slice 0's enumerable space (p = 4, d_max = 3, q_max = 2, M ≤ 2, khaos's τ² = 10, h₁ = 4,
/// h₂ = 40/n, σ² Jeffreys), the engine's posterior over the 466 structures, from its target, against
/// <c>khaos_exactlaw_oracle.R</c>'s, computed in R from <c>khaos</c>'s own <c>make_basis</c> and the ridge's closed form: P(M) and each
/// term's inclusion.</item>
/// <item><b>The readouts.</b> On Ishigami's function with one inert input (slice 0's readout fixture), the engine's fit against
/// <c>adaptive_khaos_ridge</c> with slice 0's corrections (<c>corrected=1</c>): M, σ², the surface on a grid and the total Sobol
/// indices. Two chains are not one, so the comparison is of posterior summaries within Monte Carlo reach, each difference stated.</item>
/// </list>
/// Each runs twice: on the raw inputs, in (0, 1), with the domain [0, 1] declared (<see cref="InputDomain"/>), so the engine reads
/// the values <c>khaos</c> reads; and on the inputs rescaled to [0, 1] exactly, the box the engine's range scaling reads without a
/// domain. The two are different posteriors: on Ishigami's raw inputs khaos reads M about 10.2 and σ² about 0.23 against the noise's
/// 0.25, and on the rescaled ones about 6.2 and 0.95. A sparse Legendre posterior is not invariant to an affine map of the inputs
/// (likely Ishigami's exact symmetry about ½, which keeps its raw expansion sparse).
/// </summary>
public sealed class KhaosEngineFitParityTests(ITestOutputHelper output)
{
    // Inputs rescaled to [0, 1] exactly, so the engine's range scaling is the identity and R reads the same values.
    private static double[][] Rescaled(double[][] rows, int inputs)
    {
        double[][] copy = [.. rows.Select(r => (double[])r.Clone())];
        for (int j = 0; j < inputs; j++)
        {
            double min = copy.Min(r => r[j]), max = copy.Max(r => r[j]);
            foreach (double[] r in copy) r[j] = (r[j] - min) / (max - min);
        }
        return copy;
    }

    // The rows as data, each input declaring `domain` (null: none).
    private static BassData ToData(double[][] rows, int inputs, InputDomain? domain = null)
        => new([.. Enumerable.Range(0, inputs).Select(j => new ContinuousPredictor($"x{j}", [.. rows.Select(r => r[j])], domain))],
               new ScalarResponse([.. rows.Select(r => r[inputs])]));

    // Raw rows with [0, 1] declared, or the rows rescaled onto it.
    private static (double[][] Rows, InputDomain? Domain) Case(double[][] raw, int inputs, bool declared)
        => declared ? (raw, new InputDomain(0.0, 1.0)) : (Rescaled(raw, inputs), null);

    private static BassConfig Config(int maxDegree, int maxOrder, double tau2, int n, int? maxBasis, BassSampling sampling) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = maxOrder, MaxBasis = maxBasis, Legendre = new LegendreBasis { MaxDegree = maxDegree } },
        Count = new BassCountPrior(4, 40.0 / n),
        Coefficients = BassCoefficientPrior.RawRidge(TauPrior.Held(1.0 / tau2)),
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.Jeffreys,
        Sampling = sampling,
    };

    // Slice 0's exact-law fixture (n = 40, its seed), its inputs in (0, 1).
    private static double[][] LawRows()
    {
        var rng = new Xoshiro256PlusPlus(20261109);
        double[][] rows = [.. Enumerable.Range(0, 40).Select(_ =>
        {
            double[] x = [.. Enumerable.Range(0, 4).Select(_ => rng.NextDouble())];
            double y = 0.5 * ShiftedLegendre.Value(1, x[0]) + 0.35 * ShiftedLegendre.Value(1, x[1]) * ShiftedLegendre.Value(1, x[2]) + Normal.SampleStandard(rng);
            return x.Append(y).ToArray();
        })];
        return rows;
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void ExactLaw_IsTheRReferencesUnderTheRawRidge(bool declared)
    {
        if (!ROracle.IsAvailable) return;
        (double[][] rows, InputDomain? domain) = Case(LawRows(), 4, declared);
        if (GbassOracle.Run("khaos_exactlaw_oracle.R", rows, "ridge", "corrected", "200", "3") is not { } r) return;
        BassModel model = BassModel.Create(Config(3, 2, 10.0, rows.Length, 2, new BassSampling { Chains = 1, BurnIn = 0, Samples = 1 }), ToData(rows, 4, domain));
        List<BassTerm> terms = [.. model.Prior.Legendre!.Terms()];
        double logTau = Math.Log(model.Target.HeldTau!.Value);
        var sets = new List<BassTerm[]> { Array.Empty<BassTerm>() };
        for (int a = 0; a < terms.Count; a++)
        {
            sets.Add([terms[a]]);
            for (int b = a + 1; b < terms.Count; b++) sets.Add([terms[a], terms[b]]);
        }
        double[] log = [.. sets.Select(s =>
        {
            BassState state = new BassStructureMoves(model, s, logTau).Start;
            return model.Prior.LogCollectionWeight(s.Length) + s.Sum(model.Prior.LogPrior) + model.Target.Evaluate(state).LogLikelihood(1.0);
        })];
        double top = log.Max();
        double[] post = [.. log.Select(v => Math.Exp(v - top))];
        double total = post.Sum();
        for (int k = 0; k < post.Length; k++) post[k] /= total;

        JsonElement exact = r.GetProperty("exact");
        Assert.Equal(sets.Count, exact.GetProperty("structures").GetInt32());
        double[] rM = GbassOracle.Numbers(exact.GetProperty("M"));
        double[] eM = [.. Enumerable.Range(0, 3).Select(m => Enumerable.Range(0, sets.Count).Where(k => sets[k].Length == m).Sum(k => post[k]))];
        double worst = rM.Zip(eM, (x, y) => Math.Abs(x - y)).Max();

        // Each term's inclusion, R's terms keyed by their multi-index.
        double[] rInclusion = GbassOracle.Numbers(exact.GetProperty("inclusion"));
        int[][] rTerms = [.. r.GetProperty("terms").EnumerateArray().Select(t => GbassOracle.Numbers(t).Select(v => (int)v).ToArray())];
        Assert.Equal(terms.Count, rTerms.Length);
        for (int i = 0; i < rTerms.Length; i++)
        {
            var term = new BassTerm(legendre: rTerms[i].Select((d, v) => (d, v)).Where(t => t.d > 0).Select(t => new LegendreFactor(t.v, t.d)));
            double engine = Enumerable.Range(0, sets.Count).Where(k => sets[k].Contains(term)).Sum(k => post[k]);
            worst = Math.Max(worst, Math.Abs(engine - rInclusion[i]));
        }
        output.WriteLine($"{(declared ? "raw inputs, [0, 1] declared" : "inputs rescaled")}; raw ridge on 466 structures: P(M = 0, 1, 2) engine {string.Join(" / ", eM.Select(v => v.ToString("F6")))}, R {string.Join(" / ", rM.Select(v => v.ToString("F6")))}; " +
                         $"worst |Δ| over P(M) and the 30 inclusions {worst:E2}");
        Assert.True(worst < 1e-10, $"the engine's exact law differs from R's by {worst:E2}");
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void Readouts_MeetTheCorrectedRidge_OnIshigami(bool declared)
    {
        if (!ROracle.IsAvailable) return;
        const int inputs = KhaosReadoutParityTests.Inputs, n = KhaosReadoutParityTests.Observations;
        (double[][] rows, InputDomain? domain) = Case(KhaosReadoutParityTests.Data(n, 20261107), inputs, declared);
        double[][] grid = KhaosReadoutParityTests.Grid();
        string gridFile = GbassOracle.WriteCsv(grid);
        JsonElement? reference;
        try { reference = GbassOracle.Run("khaos_readouts_oracle.R", rows, "ridge", gridFile, "6000", "3000", "11", "degree=8", "order=3", "corrected=1"); }
        finally { System.IO.File.Delete(gridFile); }
        if (reference is not { } r) return;

        BassData data = ToData(rows, inputs, domain);
        // Births are rare here (about 2% accepted on the raw inputs), so a chain started at M = 0 needs its 10,000 steps to grow: at 3,000
        // the raw inputs' chains still sat at M 9.8 against the law's 10.1 (eight engine chains of 30,000 and two khaos runs agree on it).
        var sampling = new BassSampling { Chains = 2, MasterSeed = 11, BurnIn = 10000, Samples = 10000, ReportEvidence = false };
        var clock = System.Diagnostics.Stopwatch.StartNew();
        BassFit fit = Bass.Fit(Config(8, 3, 1e5, n, null, sampling), data);
        double seconds = clock.Elapsed.TotalSeconds;
        BassDraw[] draws = [.. fit.Draws.SelectMany(s => s)];
        double engineM = draws.Average(d => d.Structure.TermIds.Length), engineS2 = draws.Average(d => d.Sigma2);
        double rM = r.GetProperty("M").GetProperty("mean").GetDouble(), rS2 = r.GetProperty("s2").GetProperty("mean").GetDouble();

        var cases = new PredictionCases([.. Enumerable.Range(0, inputs).Select(j => new ContinuousColumn($"x{j}", [.. grid.Select(g => g[j])]))]);
        double[,] mean = fit.Predict(data, cases, new BassPredictOptions { Predictive = false }).Mean();
        double[] rSurface = GbassOracle.Numbers(r.GetProperty("surface").GetProperty("mean"));
        double surfaceGap = Math.Sqrt(Enumerable.Range(0, grid.Length).Average(i => (mean[i, 0] - rSurface[i]) * (mean[i, 0] - rSurface[i])));
        double truthGap = Math.Sqrt(Enumerable.Range(0, grid.Length).Average(i => Math.Pow(mean[i, 0] - KhaosReadoutParityTests.Ishigami(grid[i]), 2)));

        double[] engineT = fit.Sensitivity(data).Decompositions[0].MeanTotalIndices();
        double[] rT = GbassOracle.Numbers(r.GetProperty("sobol").GetProperty("T_f_mean"));
        double sobolGap = engineT.Zip(rT, (x, y) => Math.Abs(x - y)).Max();
        int repeated = draws.Sum(d => d.Structure.TermIds.Count(id => fit.Terms[id].Legendre.Select(f => f.Variable).Distinct().Count() != fit.Terms[id].Legendre.Count));

        output.WriteLine($"{(declared ? "raw inputs, [0, 1] declared" : "inputs rescaled to [0, 1]")}: engine {seconds:F1} s, R {r.GetProperty("seconds").GetDouble():F1} s; M {engineM:F2} against {rM:F2}; σ² {engineS2:G4} against {rS2:G4} (noise 0.25)");
        output.WriteLine($"the surface on 16 points: engine against R's mean, root mean square {surfaceGap:F3}; engine against the function {truthGap:F3} (its sd 3.72)");
        output.WriteLine($"total indices over Var f: engine {string.Join(", ", engineT.Select((t, j) => $"x{j + 1} {t:F3}"))}; R {string.Join(", ", rT.Select((t, j) => $"x{j + 1} {t:F3}"))}; worst |Δ| {sobolGap:F3}");
        output.WriteLine($"kept terms repeating a variable: {repeated}; moves {string.Join(", ", fit.Diagnostics.MoveKinds.Entries.Where(e => e.Attempts > 0).Select(e => $"{e.Kind} {e.Accepted}/{e.Attempts}"))}");
        Assert.Equal(0, repeated);
        Assert.True(Math.Abs(engineM - rM) < 2.0, $"M {engineM:F2} against {rM:F2}");
        Assert.True(Math.Abs(engineS2 / rS2 - 1.0) < 0.15, $"σ² {engineS2:G4} against {rS2:G4}");
        Assert.True(surfaceGap < 0.3, $"the surfaces differ by {surfaceGap:F3}");
        Assert.True(sobolGap < 0.04, $"the total indices differ by up to {sobolGap:F3}");
    }

    // ---- KHAOS's modified g-prior (slice 2) ----------------------------------------------------------------------------------------

    // adaptive_khaos_gprior's defaults: λ ~ Gamma(4, rate 4/n), g₀² ~ IG(0.001, 1,000) as τ = 1/g₀² ~ Gamma(0.001, rate 1,000), ζ = 1,
    // σ² Jeffreys, the intercept inside under the zero centre.
    private static BassConfig GConfig(int maxDegree, int maxOrder, int n, int? maxBasis, BassSampling sampling) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = maxOrder, MaxBasis = maxBasis, Legendre = new LegendreBasis { MaxDegree = maxDegree } },
        Count = new BassCountPrior(4, 4.0 / n),
        Coefficients = new BassCoefficientPrior(new ModifiedGForm(1.0), TauPrior.Gamma(0.001, 1000)),
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.Jeffreys,
        Sampling = sampling,
    };

    private static double Simpson(Func<double, double> h, double a, double b, double top)
    {
        double step = (b - a) / 2000.0, sum = 0.0;
        for (int i = 0; i <= 2000; i++)
        {
            double f = Math.Exp(h(a + i * step) - top);
            sum += i == 0 || i == 2000 ? f : i % 2 == 1 ? 4.0 * f : 2.0 * f;
        }
        return step / 3.0 * sum;
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void ExactLaw_IsTheRReferencesUnderTheModifiedGPrior(bool declared)
    {
        if (!ROracle.IsAvailable) return;
        (double[][] rows, InputDomain? domain) = Case(LawRows(), 4, declared);
        if (GbassOracle.Run("khaos_exactlaw_oracle.R", rows, "gprior", "corrected", "200", "3") is not { } r) return;
        BassModel model = BassModel.Create(GConfig(3, 2, rows.Length, 2, new BassSampling { Chains = 1, BurnIn = 0, Samples = 1 }), ToData(rows, 4, domain));
        List<BassTerm> terms = [.. model.Prior.Legendre!.Terms()];
        var sets = new List<BassTerm[]> { Array.Empty<BassTerm>() };
        for (int a = 0; a < terms.Count; a++)
        {
            sets.Add([terms[a]]);
            for (int b = a + 1; b < terms.Count; b++) sets.Add([terms[a], terms[b]]);
        }
        // Each structure's weight with g₀² integrated over u = ln τ by Simpson's rule on the oracle's nodes, reflected: 2,000 panels
        // between g₀² = 10¹⁶, the cuts and 10⁻⁶.
        double[] cuts = [300.0, 1000.0, 10000.0];
        double[] bounds = [.. new[] { 1e16 }.Concat(cuts.Reverse()).Append(1e-6).Select(g => -Math.Log(g))];
        var log = new double[sets.Count];
        var below = new double[sets.Count][];
        for (int k = 0; k < sets.Count; k++)
        {
            BassState state = new BassStructureMoves(model, sets[k], 0.0).Start;
            double F(double u)
            {
                BassEvaluation e = model.Target.Evaluate(state.WithLogTau(u));
                return e.LogPrior + e.LogLikelihood(1.0);
            }
            double top = Enumerable.Range(0, 2001).Max(i => F(bounds[0] + (bounds[^1] - bounds[0]) * i / 2000.0));
            double[] parts = [.. Enumerable.Range(0, bounds.Length - 1).Select(j => Simpson(F, bounds[j], bounds[j + 1], top))];
            log[k] = top + Math.Log(parts.Sum());
            below[k] = [.. Enumerable.Range(0, cuts.Length).Select(c => parts.Skip(cuts.Length - c).Sum() / parts.Sum())];
        }
        double max = log.Max();
        double[] post = [.. log.Select(v => Math.Exp(v - max))];
        double total = post.Sum();
        for (int k = 0; k < post.Length; k++) post[k] /= total;

        JsonElement exact = r.GetProperty("exact");
        Assert.Equal(sets.Count, exact.GetProperty("structures").GetInt32());
        double[] rM = GbassOracle.Numbers(exact.GetProperty("M"));
        double[] eM = [.. Enumerable.Range(0, 3).Select(m => Enumerable.Range(0, sets.Count).Where(k => sets[k].Length == m).Sum(k => post[k]))];
        double worst = rM.Zip(eM, (x, y) => Math.Abs(x - y)).Max();
        double[] rInclusion = GbassOracle.Numbers(exact.GetProperty("inclusion"));
        int[][] rTerms = [.. r.GetProperty("terms").EnumerateArray().Select(t => GbassOracle.Numbers(t).Select(v => (int)v).ToArray())];
        for (int i = 0; i < rTerms.Length; i++)
        {
            var term = new BassTerm(legendre: rTerms[i].Select((d, v) => (d, v)).Where(t => t.d > 0).Select(t => new LegendreFactor(t.v, t.d)));
            worst = Math.Max(worst, Math.Abs(Enumerable.Range(0, sets.Count).Where(k => sets[k].Contains(term)).Sum(k => post[k]) - rInclusion[i]));
        }
        double[] rCdf = GbassOracle.Numbers(exact.GetProperty("g0sq_cdf"));
        double[] eCdf = [.. Enumerable.Range(0, cuts.Length).Select(c => Enumerable.Range(0, sets.Count).Sum(k => post[k] * below[k][c]))];
        double cdfGap = rCdf.Zip(eCdf, (x, y) => Math.Abs(x - y)).Max();
        output.WriteLine($"{(declared ? "raw inputs, [0, 1] declared" : "inputs rescaled")}; modified g-prior, g₀² integrated, on 466 structures: P(M = 0, 1, 2) engine " +
                         $"{string.Join(" / ", eM.Select(v => v.ToString("F6")))}, R {string.Join(" / ", rM.Select(v => v.ToString("F6")))}; worst |Δ| over P(M) and the 30 inclusions {worst:E2}; " +
                         $"P(g₀² < 300, 1,000, 10⁴) engine {string.Join(" / ", eCdf.Select(v => v.ToString("F6")))}, worst |Δ| {cdfGap:E2}");
        Assert.True(worst < 1e-9, $"the engine's exact law differs from R's by {worst:E2}");
        Assert.True(cdfGap < 1e-9, $"g₀²'s posterior CDF differs from R's by {cdfGap:E2}");
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void Readouts_MeetTheCorrectedGPrior_OnIshigami(bool declared)
    {
        if (!ROracle.IsAvailable) return;
        const int inputs = KhaosReadoutParityTests.Inputs, n = KhaosReadoutParityTests.Observations;
        (double[][] rows, InputDomain? domain) = Case(KhaosReadoutParityTests.Data(n, 20261107), inputs, declared);
        double[][] grid = KhaosReadoutParityTests.Grid();
        string gridFile = GbassOracle.WriteCsv(grid);
        JsonElement? reference;
        try { reference = GbassOracle.Run("khaos_readouts_oracle.R", rows, "gprior", gridFile, "6000", "3000", "11", "degree=8", "order=3", "corrected=1"); }
        finally { System.IO.File.Delete(gridFile); }
        if (reference is not { } r) return;

        BassData data = ToData(rows, inputs, domain);
        var sampling = new BassSampling { Chains = 2, MasterSeed = 11, BurnIn = 10000, Samples = 10000, ReportEvidence = false };
        var clock = System.Diagnostics.Stopwatch.StartNew();
        BassFit fit = Bass.Fit(GConfig(8, 3, n, null, sampling), data);
        double seconds = clock.Elapsed.TotalSeconds;
        BassDraw[] draws = [.. fit.Draws.SelectMany(s => s)];
        double engineM = draws.Average(d => d.Structure.TermIds.Length), engineS2 = draws.Average(d => d.Sigma2);
        double rM = r.GetProperty("M").GetProperty("mean").GetDouble(), rS2 = r.GetProperty("s2").GetProperty("mean").GetDouble();
        double[] lnG = [.. draws.Select(d => -Math.Log(d.Tau)).Order()];
        double engineG = Math.Exp(lnG[lnG.Length / 2]), rG = r.GetProperty("g0sq").GetProperty("q50").GetDouble();

        var cases = new PredictionCases([.. Enumerable.Range(0, inputs).Select(j => new ContinuousColumn($"x{j}", [.. grid.Select(g => g[j])]))]);
        double[,] mean = fit.Predict(data, cases, new BassPredictOptions { Predictive = false }).Mean();
        double[] rSurface = GbassOracle.Numbers(r.GetProperty("surface").GetProperty("mean"));
        double surfaceGap = Math.Sqrt(Enumerable.Range(0, grid.Length).Average(i => (mean[i, 0] - rSurface[i]) * (mean[i, 0] - rSurface[i])));
        double truthGap = Math.Sqrt(Enumerable.Range(0, grid.Length).Average(i => Math.Pow(mean[i, 0] - KhaosReadoutParityTests.Ishigami(grid[i]), 2)));
        double[] engineT = fit.Sensitivity(data).Decompositions[0].MeanTotalIndices();
        double[] rT = GbassOracle.Numbers(r.GetProperty("sobol").GetProperty("T_f_mean"));
        double sobolGap = engineT.Zip(rT, (x, y) => Math.Abs(x - y)).Max();

        output.WriteLine($"{(declared ? "raw inputs, [0, 1] declared" : "inputs rescaled to [0, 1]")}, modified g-prior: engine {seconds:F1} s, R {r.GetProperty("seconds").GetDouble():F1} s; " +
                         $"M {engineM:F2} against {rM:F2}; σ² {engineS2:G4} against {rS2:G4} (noise 0.25); g₀²'s median {engineG:G4} against {rG:G4}");
        output.WriteLine($"the surface on 16 points: engine against R's mean, root mean square {surfaceGap:F3}; engine against the function {truthGap:F3} (its sd 3.72)");
        output.WriteLine($"total indices over Var f: engine {string.Join(", ", engineT.Select((t, j) => $"x{j + 1} {t:F3}"))}; R {string.Join(", ", rT.Select((t, j) => $"x{j + 1} {t:F3}"))}; worst |Δ| {sobolGap:F3}");
        output.WriteLine($"moves {string.Join(", ", fit.Diagnostics.MoveKinds.Entries.Where(e => e.Attempts > 0).Select(e => $"{e.Kind} {e.Accepted}/{e.Attempts}"))}");
        Assert.True(Math.Abs(engineM - rM) < 2.0, $"M {engineM:F2} against {rM:F2}");
        Assert.True(Math.Abs(engineS2 / rS2 - 1.0) < 0.15, $"σ² {engineS2:G4} against {rS2:G4}");
        Assert.True(surfaceGap < 0.3, $"the surfaces differ by {surfaceGap:F3}");
        Assert.True(sobolGap < 0.04, $"the total indices differ by up to {sobolGap:F3}");
    }
}
