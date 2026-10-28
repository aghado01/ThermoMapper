using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// BASS's collapsed target (stage 4's S5) against the stage-4 oracle's exact posterior over an enumerable structure
/// space (<c>bass_enumeration_oracle.R</c>, <c>posterior</c>), whose closed forms the oracle checks against quadrature
/// over σ² of the Gaussian integral on the explicit uncentred design. Each of the oracle's states is loaded into a
/// <see cref="BassDesign"/> from its emitted terms, and the target's marginal is integrated over u = ln τ by this
/// test's own composite Gauss–Legendre rule:
/// <list type="bullet">
/// <item>every intercept convention, σ² prior and coefficient mode on fixture A: admission, ‖r‖², Q, the structure
/// prior (D01, D02) and the marginal (F-O1, F-O9);</item>
/// <item>fixture C's categorical predictor: the mixed order prior (D01, decided 2026-10-28) and its marginals;</item>
/// <item>fixture E's exact fits: the Jeffreys exclusion refuses exactly the oracle's;</item>
/// <item>fixture B's rungs: log Z_b and the slope's mean on the integrated and linear forms (D05);</item>
/// <item>R's own statistics under its zero centre at a fixed τ (<c>bass_statistics_oracle.R</c>; parity).</item>
/// </list>
/// The fixtures follow the stage-4 oracle brief's recipes with this side's noise: the oracle reads the data it is
/// given. Opt-in: it runs where the R toolchain and the BASS checkout are present.
/// </summary>
public sealed class BassTargetParityTests
{
    private static readonly string BassSource = Environment.GetEnvironmentVariable("BASS_R_SOURCE") ?? "D:/aghado01/packages/BASS";
    private readonly ITestOutputHelper _out;
    public BassTargetParityTests(ITestOutputHelper output) => _out = output;

    private static bool Available => ROracle.IsAvailable && Directory.Exists(BassSource);

    private static double Hinge(double x, double knot, int sign) => Math.Max(sign * (x - knot), 0.0);

    private static double Noise(Xoshiro256PlusPlus rng)
        => Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());

    // Rows of (y, x₁, x₂[, z]).
    private static double[][] FixtureA(bool exact = false)
    {
        var rng = new Xoshiro256PlusPlus(seed: 11);
        var rows = new List<double[]>();
        for (int rep = 0; rep < 3; rep++)
            for (int j = 0; j < 3; j++)
                for (int i = 0; i < 3; i++)
                {
                    double x1 = i * 0.5, x2 = j * 0.5;
                    double y = exact ? 1 + 2 * Hinge(x1, 0.5, 1) : 1 + 2 * x1 + 3 * Hinge(x1, 0.5, 1) * Hinge(x2, 0.5, -1) + 0.3 * Noise(rng);
                    rows.Add([y, x1, x2]);
                }
        return [.. rows];
    }

    private static double[][] FixtureB()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7);
        var rows = new List<double[]>();
        for (int rep = 0; rep < 2; rep++)
            for (int j = 0; j < 4; j++)
                for (int i = 0; i < 4; i++)
                {
                    double x1 = i / 3.0, x2 = j / 3.0;
                    rows.Add([1 + 2 * Hinge(x1, 1 / 3.0, 1) + 3 * Hinge(x1, 1 / 3.0, 1) * Hinge(x2, 2 / 3.0, -1) + 0.25 * Noise(rng), x1, x2]);
                }
        return [.. rows];
    }

    private static double[][] FixtureC()
    {
        var rng = new Xoshiro256PlusPlus(seed: 21);
        var rows = new List<double[]>();
        for (int z = 0; z < 3; z++)
            for (int j = 0; j < 3; j++)
                for (int i = 0; i < 3; i++)
                {
                    double x1 = i * 0.5, x2 = j * 0.5;
                    rows.Add([1 + 1.5 * x1 + 1.5 * Hinge(x2, 0.5, 1) * (z == 0 ? 1 : 0) + 0.4 * Noise(rng), x1, x2, z]);
                }
        return [.. rows];
    }

    private static BassData Data(double[][] rows, bool categorical = false)
    {
        var predictors = new List<BassPredictor>
        {
            new ContinuousPredictor("x1", [.. rows.Select(r => r[1])]),
            new ContinuousPredictor("x2", [.. rows.Select(r => r[2])]),
        };
        if (categorical) predictors.Add(new CategoricalPredictor("z", [.. rows.Select(r => (int)r[3])], ["0", "1", "2"]));
        return new BassData(predictors, new ScalarResponse([.. rows.Select(r => r[0])]));
    }

    private static JsonElement Posterior(double[][] rows, params string[] options)
    {
        string csv = Path.Combine(Path.GetTempPath(), $"bass_posterior_{Guid.NewGuid():N}.csv");
        File.WriteAllText(csv, Csv(rows));
        try
        {
            return ROracle.Run("oracles/bass_enumeration_oracle.R", csv, ["posterior", .. options]);
        }
        finally
        {
            File.Delete(csv);
        }
    }

    private static string Csv(IEnumerable<double[]> rows)
    {
        var text = new StringBuilder();
        foreach (double[] row in rows) text.AppendLine(string.Join(",", row.Select(v => v.ToString("R", CultureInfo.InvariantCulture))));
        return text.ToString();
    }

    // The oracle's terms, in the C# conventions it emits.
    private static BassTerm[] Terms(JsonElement r) => [.. r.GetProperty("terms").EnumerateArray().Select(t => new BassTerm(
        Items(t.GetProperty("cont")).Select(f => new HingeFactor(f[0].GetInt32(), f[1].GetInt32(), f[2].GetInt32())),
        Items(t.GetProperty("cat")).Select(f => new CategoricalFactor(f[0].GetInt32(), Items(f[1]).Select(l => l.GetInt32())))))];

    // R's split names the states, so jsonlite writes them as an object whose values are in state order.
    private static int[][] States(JsonElement r)
    {
        JsonElement states = r.GetProperty("states");
        IEnumerable<JsonElement> each = states.ValueKind == JsonValueKind.Object
            ? states.EnumerateObject().Select(p => p.Value)
            : states.EnumerateArray();
        return [.. each.Select(s => Items(s).Select(i => i.GetInt32()).ToArray())];
    }

    // jsonlite writes an empty vector as {} and a one-element one, unboxed, as a scalar.
    private static IEnumerable<JsonElement> Items(JsonElement e) => e.ValueKind switch
    {
        JsonValueKind.Array => e.EnumerateArray(),
        JsonValueKind.Object when !e.EnumerateObject().Any() => [],
        JsonValueKind.Number => [e],
        _ => throw new InvalidOperationException($"Expected an array, found {e.ValueKind}."),
    };

    private static double[] Doubles(JsonElement e) => [.. e.EnumerateArray().Select(v => v.GetDouble())];

    private sealed record Side(BassData Data, BassConfig Config, BassResolved Resolved, BassStructurePrior Prior, BassBasis Basis,
                               BassStructureFit[] Fits, double[] LogPrior);

    // Our side of an oracle run: each state loaded into one design, admitted, with its fit and its prior over F.
    private static Side Ours(BassData data, BassConfig config, BassTerm[] terms, int[][] states)
    {
        BassResolved resolved = BassResolved.Resolve(config, data);
        BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
        var prior = BassStructurePrior.From(config, resolved, inputs);
        var basis = new BassBasis(inputs, 1);
        var design = new BassDesign(basis, data, BassAdmission.From(config, resolved));
        var fits = new BassStructureFit[states.Length];
        var logPrior = new double[states.Length];
        for (int s = 0; s < states.Length; s++)
        {
            BassTerm[] structure = [.. states[s].Select(i => terms[i])];
            BassRefusal? refusal = design.Load(structure);
            Assert.True(refusal is null, $"state {s} ({string.Join(" + ", structure.Select(t => t.ToString()))}) refused: {refusal}");
            double termPrior = structure.Sum(prior.LogPrior);
            fits[s] = BassStructureFit.Of(design, termPrior);
            logPrior[s] = prior.LogCollectionWeight(structure.Length) + termPrior;
        }
        double z = LogSumExp(logPrior);
        for (int s = 0; s < logPrior.Length; s++) logPrior[s] -= z;
        return new(data, config, resolved, prior, basis, fits, logPrior);
    }

    private static BassConfig Config(GPriorCentre centre = GPriorCentre.Mean, bool jeffreys = false, bool fixedG = false, int maxInteraction = 2)
        => new()
        {
            Structure = new BassStructureSpec { MaxInteraction = maxInteraction, Support = new SupportRule(Minimum: 3) },
            Intercept = centre,
            Sigma2 = jeffreys ? Sigma2Prior.Jeffreys : Sigma2Prior.ScaledToResponse(1e-3),
            Coefficients = fixedG ? new FixedG() : new SampledG(),
        };

    // ln ∫ p(u) exp(ℓ_b(S, eᵘ)) du over the g mode's τ, or ℓ_b at the held τ; with the u-nodes' log weights when asked.
    private static double LogMarginal(BassTarget target, BassStructureFit fit, double b, List<(double U, double LogW)> nodes,
                                      double[]? weights = null, double[]? slopes = null)
    {
        if (!target.SamplesTau)
        {
            BassEvaluation held = target.Evaluate(fit);
            if (slopes is not null)
            {
                Array.Clear(weights!);
                Array.Clear(slopes);
                weights![0] = 1.0;
                slopes[0] = held.LogLikelihoodSlope(b);
            }
            return held.LogLikelihood(b);
        }
        var terms = new double[nodes.Count];
        for (int i = 0; i < nodes.Count; i++)
        {
            BassEvaluation ev = target.Evaluate(fit, nodes[i].U);
            terms[i] = nodes[i].LogW + target.LogTauDensity(nodes[i].U) + ev.LogLikelihood(b);
            if (slopes is not null) slopes[i] = ev.LogLikelihoodSlope(b);
        }
        double total = LogSumExp(terms);
        if (weights is not null)
            for (int i = 0; i < terms.Length; i++) weights[i] = Math.Exp(terms[i] - total);
        return total;
    }

    private static List<(double U, double LogW)> TauNodes(int n)
    {
        double centre = Math.Log(1.0 / n);
        return [.. Nodes(centre - 80, centre + 8).Select(p => (p.U, Math.Log(p.W)))];
    }

    // Composite 12-point Gauss–Legendre on panels of width ½ over [lo, hi], its nodes by Newton's method on P₁₂.
    private static IEnumerable<(double U, double W)> Nodes(double lo, double hi)
    {
        const int order = 12;
        var x = new double[order];
        var w = new double[order];
        for (int i = 0; i < order; i++)
        {
            double z = Math.Cos(Math.PI * (i + 0.75) / (order + 0.5));
            for (int iteration = 0; iteration < 100; iteration++)
            {
                (double value, double slope) = Legendre(order, z);
                double next = z - value / slope;
                bool done = Math.Abs(next - z) <= 1e-16;
                z = next;
                if (done) break;
            }
            double d = Legendre(order, z).Slope;
            x[i] = z;
            w[i] = 2 / ((1 - z * z) * d * d);
        }
        int panels = (int)Math.Ceiling((hi - lo) / 0.5);
        double width = (hi - lo) / panels;
        for (int p = 0; p < panels; p++)
            for (int i = 0; i < order; i++)
                yield return (lo + width * (p + 0.5 * (x[i] + 1)), 0.5 * width * w[i]);
    }

    private static (double Value, double Slope) Legendre(int n, double z)
    {
        double p0 = 1, p1 = z;
        for (int k = 2; k <= n; k++) (p0, p1) = (p1, ((2 * k - 1) * z * p1 - (k - 1) * p0) / k);
        return (p1, n * (z * p1 - p0) / (z * z - 1));
    }

    private static double LogSumExp(IReadOnlyList<double> values)
    {
        double max = values.Max();
        if (double.IsNegativeInfinity(max)) return max;
        double sum = 0;
        foreach (double v in values) sum += Math.Exp(v - max);
        return max + Math.Log(sum);
    }

    [Theory]
    [InlineData("mean")]
    [InlineData("flat")]
    [InlineData("zero")]
    [InlineData("jeffreys")]
    [InlineData("flat + jeffreys")]
    [InlineData("fixed g")]
    public void Marginals_MatchTheOraclesPosterior(string option)
    {
        if (!Available) return;
        (string[] args, BassConfig config) = option switch
        {
            "mean" => (Array.Empty<string>(), Config()),
            "flat" => (["centre=flat"], Config(GPriorCentre.Flat)),
            "zero" => (["centre=zero"], Config(GPriorCentre.Zero)),
            "jeffreys" => (["sigma2_prior=jeffreys"], Config(jeffreys: true)),
            "flat + jeffreys" => (["centre=flat", "sigma2_prior=jeffreys"], Config(GPriorCentre.Flat, jeffreys: true)),
            "fixed g" => (["coef=fixed"], Config(fixedG: true)),
            _ => throw new ArgumentException(option),
        };
        double[][] rows = FixtureA();
        JsonElement r = Posterior(rows, ["M_max=2", .. args]);
        JsonElement p = r.GetProperty("posterior");
        Side ours = Ours(Data(rows), config, Terms(r), States(r));
        var target = new BassTarget(config, ours.Resolved, ours.Prior);
        Compare(option, r, ours, target, config.Intercept);
        Assert.Equal(p.GetProperty("response").GetProperty("absolute").GetBoolean(), config.Sigma2.IsProper && config.Intercept != GPriorCentre.Flat);
    }

    // Per state: ‖r‖², Q in the convention's own design, the prior over F, and the marginal (absolute, or relative to
    // the first state's under an improper prior).
    private void Compare(string name, JsonElement r, Side ours, BassTarget target, GPriorCentre centre)
    {
        JsonElement p = r.GetProperty("posterior");
        double[] rss = Doubles(p.GetProperty("rss")), q = Doubles(p.GetProperty("Q"));
        double[] logPrior = Doubles(p.GetProperty("log_prior")), logMarginal = Doubles(p.GetProperty("log_marginal"));
        bool absolute = p.GetProperty("response").GetProperty("absolute").GetBoolean();
        int n = ours.Resolved.Observations;
        double level = centre == GPriorCentre.Zero ? n * ours.Resolved.ResponseMean * ours.Resolved.ResponseMean : 0.0;
        double scale = n * ours.Resolved.ResponseMeanSquare;
        List<(double U, double LogW)> nodes = TauNodes(n);
        var mine = new double[rss.Length];
        double worstRss = 0, worstQ = 0, worstPrior = 0, worstMarginal = 0;
        for (int s = 0; s < rss.Length; s++)
        {
            BassStructureFit fit = ours.Fits[s];
            worstRss = Math.Max(worstRss, Math.Abs(fit.ResidualSquares - rss[s]) / scale);
            worstQ = Math.Max(worstQ, Math.Abs(fit.FittedSquares + level - q[s]) / Math.Max(scale, q[s]));
            worstPrior = Math.Max(worstPrior, Math.Abs(ours.LogPrior[s] - logPrior[s]));
            mine[s] = LogMarginal(target, fit, 1.0, nodes);
        }
        for (int s = 0; s < rss.Length; s++)
        {
            double gap = absolute ? mine[s] - logMarginal[s] : (mine[s] - mine[0]) - (logMarginal[s] - logMarginal[0]);
            worstMarginal = Math.Max(worstMarginal, Math.Abs(gap));
        }
        _out.WriteLine($"{name}: {rss.Length} states; ‖r‖² {worstRss:E1}, Q {worstQ:E1} (of ‖y − level‖²), prior {worstPrior:E1}, marginal {worstMarginal:E1}{(absolute ? "" : " (relative)")}");
        Assert.True(worstRss <= 1e-12, $"{name}: ‖r‖² off by {worstRss:E2}");
        Assert.True(worstQ <= 1e-12, $"{name}: Q off by {worstQ:E2}");
        Assert.True(worstPrior <= 1e-12, $"{name}: the prior off by {worstPrior:E2}");
        Assert.True(worstMarginal <= 1e-9, $"{name}: the marginal off by {worstMarginal:E2}");
    }

    [Fact]
    public void Categorical_PriorAndMarginals_MatchTheOracle()
    {
        if (!Available) return;
        double[][] rows = FixtureC();
        JsonElement r = Posterior(rows, "M_max=2", "K_max=1", "categorical=1");
        BassConfig config = Config(maxInteraction: 1);
        Side ours = Ours(Data(rows, categorical: true), config, Terms(r), States(r));
        Assert.Equal(0, ours.Prior.MinimumOrder(BassPart.Categorical));
        Compare("categorical", r, ours, new BassTarget(config, ours.Resolved, ours.Prior), GPriorCentre.Mean);
    }

    [Fact]
    public void ExactFits_AreTheOraclesExclusion()
    {
        if (!Available) return;
        double[][] rows = FixtureA(exact: true);
        JsonElement all = Posterior(rows, "M_max=2");
        JsonElement kept = Posterior(rows, "M_max=2", "sigma2_prior=jeffreys");
        BassTerm[] terms = Terms(all);
        var keptStates = States(kept).Select(s => string.Join(",", s)).ToHashSet();

        BassConfig config = Config(jeffreys: true);
        Side ours = Ours(Data(rows), config, terms, States(all));
        var target = new BassTarget(config, ours.Resolved, ours.Prior);
        int[][] states = States(all);
        int refused = 0;
        for (int s = 0; s < states.Length; s++)
        {
            bool exact = target.Refusal(ours.Fits[s]) == BassRefusal.ExactFit;
            refused += exact ? 1 : 0;
            Assert.Equal(!keptStates.Contains(string.Join(",", states[s])), exact);
        }
        int expected = kept.GetProperty("posterior").GetProperty("n_exact_fit_refused").GetInt32();
        _out.WriteLine($"exact fits: {refused} refused of {states.Length}, the oracle's {expected}");
        Assert.Equal(expected, refused);
        Assert.True(refused > 0);
    }

    // Under flat + Jeffreys the integrated form is undefined at b = 0 and the evidence is null, but the rungs' laws, and
    // so the slope's mean, are defined wherever the form is.
    [Theory]
    [InlineData("defaults")]
    [InlineData("flat + jeffreys")]
    public void Rungs_MatchTheOraclesLadder(string option)
    {
        if (!Available) return;
        bool improper = option != "defaults";
        double[][] rows = FixtureB();
        double[] rungs = [0.0, 0.1, 0.5, 1.0];
        string[] args = ["M_max=2", "rungs=" + string.Join(",", rungs.Select(b => b.ToString(CultureInfo.InvariantCulture))),
                         .. improper ? new[] { "centre=flat", "sigma2_prior=jeffreys" } : []];
        JsonElement r = Posterior(rows, args);
        BassConfig config = improper ? Config(GPriorCentre.Flat, jeffreys: true) : Config();
        Side ours = Ours(Data(rows), config, Terms(r), States(r));
        List<(double U, double LogW)> nodes = TauNodes(ours.Resolved.Observations);
        JsonElement tempered = r.GetProperty("posterior").GetProperty("tempered");
        foreach ((string key, TemperedForm form) in new[] { ("integrated", TemperedForm.Integrated), ("linear", TemperedForm.Linear) })
        {
            var target = new BassTarget(config, ours.Resolved, ours.Prior, form);
            JsonElement[] reference = [.. tempered.GetProperty(key).EnumerateArray()];
            for (int k = 0; k < rungs.Length; k++)
            {
                double b = rungs[k];
                Assert.Equal(reference[k].GetProperty("defined").GetBoolean(), target.Defines(b));
                if (!target.Defines(b))
                {
                    Assert.Throws<ArgumentOutOfRangeException>(() => target.Evaluate(ours.Fits[0], nodes[0].U).LogLikelihood(b));
                    _out.WriteLine($"{option}, {key}, b = {b}: undefined on both sides");
                    continue;
                }
                var lm = new double[ours.Fits.Length];
                var meanSlope = new double[ours.Fits.Length];
                var weights = new double[nodes.Count];
                var slopes = new double[nodes.Count];
                for (int s = 0; s < ours.Fits.Length; s++)
                {
                    lm[s] = LogMarginal(target, ours.Fits[s], b, nodes, weights, slopes);
                    meanSlope[s] = weights.Zip(slopes, (w, d) => w * d).Sum();
                }
                double logZ = LogSumExp([.. ours.LogPrior.Zip(lm, (a, c) => a + c)]);
                double slope = 0;
                for (int s = 0; s < lm.Length; s++) slope += Math.Exp(ours.LogPrior[s] + lm[s] - logZ) * meanSlope[s];
                JsonElement expectedZ = reference[k].GetProperty("log_Z");
                JsonElement expectedSlope = reference[k].GetProperty("mean_slope");
                _out.WriteLine($"{option}, {key}, b = {b}: log Z {logZ:F9} against {(expectedZ.ValueKind == JsonValueKind.Null ? "null" : expectedZ.GetDouble().ToString("F9"))}; E[slope] {slope:F6} against {(expectedSlope.ValueKind == JsonValueKind.Null ? "null" : expectedSlope.GetDouble().ToString("F6"))}");
                if (improper) Assert.Equal(JsonValueKind.Null, expectedZ.ValueKind);
                else Assert.True(Math.Abs(logZ - expectedZ.GetDouble()) <= 1e-8, $"{key}, b = {b}: log Z {logZ:R} against {expectedZ.GetDouble():R}");
                if (expectedSlope.ValueKind != JsonValueKind.Null)
                    Assert.True(Math.Abs(slope - expectedSlope.GetDouble()) <= 1e-7 * Math.Max(1.0, Math.Abs(slope)),
                                $"{key}, b = {b}: E[slope] {slope:R} against {expectedSlope.GetDouble():R}");
            }
        }
    }

    [Fact]
    public void ZeroCentre_MatchesRsStatistics()
    {
        if (!Available) return;
        // The generating structure on fixture A, under R's convention: the uncentred design [1, B] with y as given.
        double[][] rows = FixtureA();
        BassData data = Data(rows);
        BassConfig config = Config(GPriorCentre.Zero);
        BassResolved resolved = BassResolved.Resolve(config, data);
        BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
        var prior = BassStructurePrior.From(config, resolved, inputs);
        var basis = new BassBasis(inputs, 1);
        BassTerm[] structure = [new BassTerm([new HingeFactor(0, 1, 1)]), new BassTerm([new HingeFactor(0, 1, 1), new HingeFactor(1, -1, 1)])];
        var design = new BassDesign(basis, data, BassAdmission.From(config, resolved));
        Assert.Null(design.Load(structure));
        BassStructureFit fit = BassStructureFit.Of(design, structure.Sum(prior.LogPrior));

        int n = rows.Length;
        var columns = structure.Select(t => { var c = new double[n]; basis.CaseColumn(t, inputs.Cases, c); return c; }).ToArray();
        string csv = Path.Combine(Path.GetTempPath(), $"bass_statistics_{Guid.NewGuid():N}.csv");
        File.WriteAllText(csv, Csv(Enumerable.Range(0, n).Select(i => new[] { rows[i][0], 1.0 }.Concat(columns.Select(c => c[i])).ToArray())));
        const double tau = 0.05;
        try
        {
            JsonElement r = ROracle.Run("oracles/bass_statistics_oracle.R", csv, "scalar", $"tau={tau.ToString(CultureInfo.InvariantCulture)}");
            double level = n * resolved.ResponseMean * resolved.ResponseMean;
            double qf = r.GetProperty("qf").GetDouble(), yty = r.GetProperty("yty").GetDouble();
            double rssDirect = r.GetProperty("check").GetProperty("rss_direct").GetDouble();
            double rateDirect = r.GetProperty("s2_rate_direct").GetDouble(), rate = r.GetProperty("s2_rate").GetDouble();
            double s = new BassTarget(config, resolved, prior).Marginal(fit, tau).S;
            _out.WriteLine($"qf {qf:R} against {fit.FittedSquares + level:R}; rss {rssDirect:R} against {fit.ResidualSquares:R}; S/2 {s / 2:R} against R's {rateDirect:R} (as written {rate:R})");
            Assert.True(Math.Abs(fit.FittedSquares + level - qf) <= 1e-13 * qf);
            Assert.True(Math.Abs(design.ResponseSquares + level - yty) <= 1e-13 * yty);
            Assert.True(Math.Abs(fit.ResidualSquares - rssDirect) <= 1e-11 * yty);
            Assert.True(Math.Abs(s / 2 - rateDirect) <= 1e-12 * rateDirect);
        }
        finally
        {
            File.Delete(csv);
        }
    }
}
