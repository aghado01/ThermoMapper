using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// S10's leave-one-out check of the pointwise contents (D04; owed since chip B's report): PSIS-LOO on the retained and
/// the collapsed content, scored unit by unit against the exact leave-one-out predictive by enumeration of F. It is
/// BARS's design (<see cref="LooBenchmark"/>) with an exact reference in place of a sampled one.
/// <list type="bullet">
///   <item><b>Data:</b> fixture A's design (x₁, x₂ on {0, ½, 1}, three replicates, N = 27, support minimum 3), on which F
///   enumerates. Three responses on it: fixture A's surface, 1 + 2x₁ + 3(x₁ − ½)₊(½ − x₂)₊ + 0.3ε; the same with an
///   outlier, +1.5 at observation 13; and with a corner spike, +0.9 at the three replicates of (1, 1), which only an
///   interaction at the corner fits.</item>
///   <item><b>Reference: the target PSIS estimates.</b> The prior is held at the full data's (the g-prior, the centre's
///   level, σ²'s scale and the support gate on all of x), as PSIS holds it. Both contents then estimate log p(yᵢ | y₋ᵢ),
///   and since p(y₋ᵢ) = Σ_S p(S)·p(y | S)/p(yᵢ | y₋ᵢ, S), it is exactly logΣ_S w_S − logΣ_S w_S/p(yᵢ | y₋ᵢ, S) with
///   w_S = p(S)·p(y | S). The leave-out term is the collapsed content's own, which <c>BassPointwiseTests</c> checks against
///   dense references. Under sampled g the sums run over u = ln τ as well, by quadrature.</item>
///   <item><b>Fits:</b> four chains of 1,000 burn-in and 2,000 draws.</item>
///   <item><b>Scores:</b> Δᵢ = PSIS − reference and z = Δ/mcse, the PSIS estimate's own error (the reference has none).
///   Flags are k̂ above 0.7, the threshold at 8,000 draws.</item>
/// </list>
/// <para>
/// Pre-registered verdicts, pooled over the three data sets, for each content: (a) its k̂ is <b>pessimistic</b> if at
/// least 90% of the units it flags have |z| ≤ 3; (b) its k̂ is a <b>false negative</b> if at least 10% of the units it
/// passes have |z| &gt; 3. The tests assert that (b) holds for neither content and report (a).
/// </para>
/// <para>
/// Both are qualification runs (<see cref="BassQualification.Enabled"/>), whose verdicts the fitting checklist records. Under
/// fixed g (about a minute) the collapsed ratios 1/p(yᵢ | y₋ᵢ, S) repeat exactly within a structure, the clumping that
/// made BARS's collapsed k̂ pessimistic; sampled g, the default, is heavier: F holds 464,094 structures here, each
/// with its quadrature in u, so the structures run in parallel. A total's Monte Carlo error is NaN whenever a unit is
/// flagged, as loo reports it.
/// </para>
/// </summary>
public sealed class BassQualificationLooTests
{
    private readonly ITestOutputHelper _out;
    public BassQualificationLooTests(ITestOutputHelper output) => _out = output;

    private const double KThreshold = 0.7;
    private static readonly string[] Kinds = ["surface", "outlier", "spike"];

    [Fact]
    public void Psis_UnderFixedG_IsNoFalseNegative()
    {
        if (!BassQualification.Enabled) return;   // heavy
        Verdicts(fixedG: true);
    }

    [Fact]
    public void Psis_UnderSampledG_IsNoFalseNegative()
    {
        if (!BassQualification.Enabled) return;   // heavy
        Verdicts(fixedG: false);
    }

    private void Verdicts(bool fixedG)
    {
        var scores = new Dictionary<string, List<(double Z, double K)>> { ["retained"] = [], ["collapsed"] = [] };
        foreach (string kind in Kinds)
        {
            BassData data = Data(kind);
            BassConfig config = Config(fixedG);
            BassFit fit = Bass.Fit(config, data);
            double[] reference = ExactLoo(BassModel.Create(config, data));
            double total = reference.Sum();
            foreach ((string content, PointwiseContent which) in new[] { ("retained", PointwiseContent.RetainedParameters), ("collapsed", PointwiseContent.CollapsedPredictive) })
            {
                PredictiveCheck loo = PredictiveChecks.PsisLoo(fit.PointwiseLogDensity(data, content: which));
                int flagged = 0;
                for (int u = 0; u < reference.Length; u++)
                {
                    double z = (loo.PerUnit[u] - reference[u]) / loo.MonteCarloErrors![u];
                    scores[content].Add((z, loo.ParetoK![u]));
                    if (loo.ParetoK[u] > KThreshold) flagged++;
                }
                double worst = Enumerable.Range(0, reference.Length).Max(u => Math.Abs(loo.PerUnit[u] - reference[u]));
                _out.WriteLine($"{(fixedG ? "fixed" : "sampled")} g, {kind}, {content}: elpd {loo.Elpd:F3} against {total:F3} " +
                               $"(± {loo.JointMonteCarloError:F3}), worst |Δ| {worst:F4}, {flagged} of {reference.Length} flagged, " +
                               $"largest k̂ {loo.ParetoK.Max():F2}");
            }
        }

        foreach ((string content, List<(double Z, double K)> s) in scores)
        {
            var passed = s.Where(p => !(p.K > KThreshold)).ToList();
            var flagged = s.Where(p => p.K > KThreshold).ToList();
            double missed = passed.Count == 0 ? 0.0 : passed.Count(p => Math.Abs(p.Z) > 3) / (double)passed.Count;
            string pessimism = flagged.Count == 0
                ? "none flagged"
                : $"{flagged.Count(p => Math.Abs(p.Z) <= 3)} of {flagged.Count} flagged within |z| ≤ 3 " +
                  (flagged.Count(p => Math.Abs(p.Z) <= 3) >= 0.9 * flagged.Count ? "(pessimistic)" : "(not pessimistic)");
            _out.WriteLine($"{content}: {passed.Count(p => Math.Abs(p.Z) > 3)} of {passed.Count} passed beyond |z| > 3; {pessimism}; " +
                           $"largest |z| among the passed {passed.Select(p => Math.Abs(p.Z)).DefaultIfEmpty(0).Max():F2}");
            Assert.True(missed < 0.10, $"The {content} content's k̂ passes {missed:P0} of units with |z| > 3: a false negative.");
        }
    }

    // log p(yᵢ | y₋ᵢ) for every observation, exact over F and, under sampled g, by quadrature in u: panels of width 2,
    // whose structure sums agreed with BassTargetTests' ½ to 1e-11 on 150 structures of each response, at a quarter of
    // the nodes. The structures run in parallel, each worker with its own term table; the sums over F run in order.
    private static double[] ExactLoo(BassModel model)
    {
        BassTarget target = model.Target;
        int n = model.Resolved.Observations;
        List<(double U, double W)> nodes = target.SamplesTau
            ? [.. BassTargetTests.Nodes(Math.Log(1.0 / n) - 80, Math.Log(1.0 / n) + 8, width: 2.0)]
            : [(double.NaN, 1.0)];
        List<(BassTerm[] Terms, double LogPrior)> f = BassSweepTests.EnumerateF(model);
        var numerators = new double[f.Count];
        var denominators = new double[f.Count, n];
        Parallel.For(0, f.Count, () => new BassTermTable(), (k, _, table) =>
        {
            (BassTerm[] terms, double logPrior) = f[k];
            var logw = new double[nodes.Count];
            var column = new double[nodes.Count];
            BassDesign design = model.NewDesign();
            Assert.Null(design.Load(terms));
            BassStructureFit fit = BassStructureFit.Of(design, 0.0);
            int[] ids = [.. terms.Select(table.Intern).Order()];
            var draws = new BassDraw[nodes.Count];
            for (int i = 0; i < nodes.Count; i++)
            {
                (double u, double w) = nodes[i];
                logw[i] = target.SamplesTau
                    ? logPrior + Math.Log(w) + target.LogTauDensity(u) + target.Evaluate(fit, u).LogMarginal
                    : logPrior + target.Evaluate(fit).LogMarginal;
                draws[i] = new BassDraw(new BassStructure(ids), target.SamplesTau ? Math.Exp(u) : target.HeldTau!.Value, 1.0,
                                        new double[ids.Length + 1], 1.0);
            }
            numerators[k] = LogSumExp(logw);
            double[,] leaveOut = BassPointwise.LogDensity(model, table, [draws], content: PointwiseContent.CollapsedPredictive).Values;
            for (int unit = 0; unit < n; unit++)
            {
                for (int i = 0; i < nodes.Count; i++) column[i] = logw[i] - leaveOut[i, unit];
                denominators[k, unit] = LogSumExp(column);
            }
            return table;
        }, _ => { });
        double numerator = LogSumExp(numerators);
        var reference = new double[n];
        var structures = new double[f.Count];
        for (int unit = 0; unit < n; unit++)
        {
            for (int k = 0; k < f.Count; k++) structures[k] = denominators[k, unit];
            reference[unit] = numerator - LogSumExp(structures);
        }
        return reference;
    }

    private static double LogSumExp(ReadOnlySpan<double> values)
    {
        double top = double.NegativeInfinity;
        foreach (double v in values) top = Math.Max(top, v);
        if (double.IsNegativeInfinity(top)) return top;
        double sum = 0.0;
        foreach (double v in values) sum += Math.Exp(v - top);
        return top + Math.Log(sum);
    }

    private static BassConfig Config(bool fixedG) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 2, Support = new SupportRule(Minimum: 3) },
        Coefficients = fixedG ? BassCoefficientPrior.G(TauPrior.Held()) : BassCoefficientPrior.G(),
        Sampling = new BassSampling { Chains = 4, MasterSeed = 23, BurnIn = 1000, Samples = 2000 },
    };

    // Fixture A's design and surface (its seed and row order: replicate, then x₂, then x₁), perturbed by `kind`.
    internal static BassData Data(string kind)
    {
        var rng = new Xoshiro256PlusPlus(seed: 11);
        var x1 = new List<double>();
        var x2 = new List<double>();
        var y = new List<double>();
        for (int rep = 0; rep < 3; rep++)
            for (int j = 0; j < 3; j++)
                for (int i = 0; i < 3; i++)
                {
                    double a = i * 0.5, b = j * 0.5;
                    double noise = Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());
                    double v = 1 + 2 * a + 3 * Math.Max(a - 0.5, 0) * Math.Max(0.5 - b, 0) + 0.3 * noise;
                    if (kind == "spike" && i == 2 && j == 2) v += 0.9;
                    x1.Add(a);
                    x2.Add(b);
                    y.Add(v);
                }
        if (kind == "outlier") y[13] += 1.5;
        return new BassData([new ContinuousPredictor("x1", [.. x1]), new ContinuousPredictor("x2", [.. x2])], new ScalarResponse([.. y]));
    }
}
