using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// <c>LOO_BENCH</c> (seeds per data set): brute-force leave-one-out for Normal BARS, against which the retained and the
/// collapsed PSIS-LOO are scored unit by unit. It asks whether the collapsed content's high k̂ is pessimistic, or the
/// retained content's low k̂ a false negative (a unit whose retained estimate is wrong while k̂ says it is fine).
/// <list type="bullet">
///   <item><b>Data:</b> n = 30 on (i + ½)/n. A sine with noise sd 0.15; the same with a narrow bump (width 0.03, about one
///   point spacing), which needs local knots and so high leverage; the sine with one outlier of 1.0 at x₁₇.</item>
///   <item><b>Full fit:</b> default BARS, four chains of 1,000 burn-in and 2,000 draws.</item>
///   <item><b>Reference: the same target as PSIS.</b> Both contents estimate leave-one-out with the prior held at the full
///   data's: the g-prior's G, τ = 1/n and median centre, and the support gate on all of x. That knot posterior is
///   p(M | y₋ᵢ) ∝ p(M)·p(y | M)/p(yᵢ | y₋ᵢ, M). BARS's own moves sample it here, on a target that subtracts the collapsed
///   leave-out term from BARS's. Four chains start at the full fit's last draws and run 1,000 burn-in and 2,000 draws.
///   The reference is the mean of p(yᵢ | y₋ᵢ, M) over those draws; its Monte Carlo error is the mean's, through the chains'
///   mean ESS. The leave-out term is exact (<see cref="CollapsedPredictiveTests"/> checks it against the dense
///   multivariate t).</item>
///   <item><b>Scores:</b> Δᵢ = PSIS − reference, and z = Δ/√(mcse²_PSIS + mcse²_ref). Flags are k̂ above 0.7, the threshold
///   at 8,000 draws.</item>
/// </list>
/// <para>
/// Pre-registered verdicts, pooled over data sets and seeds:
/// (a) <b>the collapsed k̂ is pessimistic</b> if at least 90% of the collapsed-flagged units have |z| ≤ 3;
/// (b) <b>the retained k̂ is a false negative</b> if at least 10% of the units it passes have |z| > 3, or if its mean signed Δ
/// on the units only the collapsed content flags exceeds three standard errors of that mean (systematic optimism).
/// </para>
/// <para>
/// The first design reran BARS on y₋ᵢ. Its pre-registered control, the prior's data dependence δᵢ, came out at up to 0.43
/// nats against a 0.01 limit, largest at the endpoints and the bump and outlier units: removing a point rebuilds the
/// g-prior, so the rerun answers a different question from PSIS. This design replaced it; the verdict rules are
/// unchanged.
/// </para>
/// </summary>
public sealed class LooBenchmark
{
    private readonly ITestOutputHelper _out;
    public LooBenchmark(ITestOutputHelper output) => _out = output;

    private const int N = 30, Burn = 1000, Samples = 2000;
    private static readonly SplineBasis Basis = new(3);
    private static readonly IComplexityPrior Prior = new PoissonPrior(4.0);
    private static readonly IKnotKernel Kernel = new LocalBetaKernel(50.0);

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static (double[] X, double[] Y) Data(string shape, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x = Enumerable.Range(0, N).Select(i => (i + 0.5) / N).ToArray();
        double[] y = x.Select(v => 1.5 + Math.Sin(2.0 * Math.PI * v) + 0.15 * Gaussian(rng)).ToArray();
        if (shape == "bump") for (int i = 0; i < N; i++) y[i] += 1.2 * Math.Exp(-Math.Pow((x[i] - 0.62) / 0.03, 2));
        if (shape == "outlier") y[17] += 1.0;
        return (x, y);
    }

    // log p(yᵢ | y₋ᵢ, M) with the g-prior held at the full data's: BARS's Normal block at τ = 1/n, the median centre.
    private static double LeaveOut(KnotConfig knots, double[] x, double[] y, int i)
    {
        double[,] z = Basis.Design(knots, x);
        GPriorConditional c = BandedGPrior.Fit(new BandedDesign(z), y, DesignAdmission.Admits, GPriorCentre.Median)!.At(1.0 / y.Length);
        return c.LeaveOutLogDensity(null, z, y, null, [i]);
    }

    // BARS's target without unit i: its knot posterior divided by p(yᵢ | y₋ᵢ, M), which is p(M | y₋ᵢ) under the full prior.
    private sealed class LeaveOneOutTarget(SplineTarget inner, double[] x, double[] y, int i) : ISamplingTarget<KnotConfig>
    {
        public ITargetEvaluation Evaluate(KnotConfig state)
        {
            ITargetEvaluation e = inner.Evaluate(state);
            if (!(Temperature.Posterior.LogDensity(e) > double.NegativeInfinity)) return e;
            return new Shifted(e, -LeaveOut(state, x, y, i));
        }

        private sealed class Shifted(ITargetEvaluation inner, double shift) : ITargetEvaluation
        {
            public double LogPrior => inner.LogPrior;
            public double LogLikelihood(double inverseTemperature) => inverseTemperature == 1.0
                ? inner.LogLikelihood(1.0) + shift
                : throw new NotSupportedException("The leave-one-out target is untempered.");
        }
    }

    // The reference: four chains on the leave-one-out target, started at the full fit's last draws.
    private static (double Estimate, double Error) Reference(BarsResult full, double[] x, double[] y, int i, int seed)
    {
        var target = new LeaveOneOutTarget(new SplineTarget(Basis, new WeightedNormalModel(), Prior, x, y), x, y, i);
        var logp = new double[full.Draws.Count][];
        for (int c = 0; c < full.Draws.Count; c++)
        {
            var density = new KnotBirthDensity(Kernel, 1.0, UniformKnotBase.Instance);
            IRjMove<KnotConfig>[] moves = [new KnotBirthMove(density), new KnotDeathMove(density), new KnotRelocateMove(Kernel)];
            var chain = new ReversibleJumpChain<KnotConfig>(moves, target, full.Draws[c][^1].Knots, new Xoshiro256PlusPlus(seed + 7919 * c));
            for (int s = 0; s < Burn; s++) chain.Step();
            logp[c] = new double[Samples];
            KnotConfig? last = null;
            double value = double.NaN;
            for (int s = 0; s < Samples; s++)
            {
                KnotConfig m = chain.Step();
                if (!ReferenceEquals(m, last)) { value = LeaveOut(m, x, y, i); last = m; }
                logp[c][s] = value;
            }
        }
        return LogMean(logp);
    }

    // log of the mean of exp(values) over chains, with its Monte Carlo error through the chains' mean ESS.
    private static (double Estimate, double Error) LogMean(double[][] chains)
    {
        double max = chains.SelectMany(c => c).Max();
        double[][] p = chains.Select(c => c.Select(v => Math.Exp(v - max)).ToArray()).ToArray();
        double[] all = p.SelectMany(c => c).ToArray();
        double mean = all.Average(), variance = all.Sum(v => (v - mean) * (v - mean)) / (all.Length - 1);
        EssDiagnostic ess = ChainDiagnostics.MeanEss(p);
        double n = ess.State == DiagnosticState.Estimated ? ess.Value : all.Length;
        return (max + Math.Log(mean), Math.Sqrt(Math.Log(1.0 + variance / (mean * mean) / n)));
    }

    private sealed record Unit(string Data, int Seed, int Index, double KRet, double KCol, double DRet, double DCol,
                               double ZRet, double ZCol, double RefError);

    [Fact]
    public void BruteForceLoo_ScoresTheRetainedAndCollapsedPsis()
    {
        if (!int.TryParse(Environment.GetEnvironmentVariable("LOO_BENCH"), out int seeds) || seeds <= 0) return;   // heavy
        var units = new List<Unit>();
        foreach (string shape in new[] { "sine", "bump", "outlier" })
            for (int seed = 1; seed <= seeds; seed++)
            {
                (double[] x, double[] y) = Data(shape, seed);
                var bars = new BarsEnsemble(Basis, new WeightedNormalModel(), Prior, Kernel);
                BarsResult full = bars.Run(x, y, grid: x, chains: 4, masterSeed: 100 * seed, burn: Burn, samples: Samples);
                PredictiveCheck ret = PredictiveChecks.PsisLoo(bars.PointwiseLogDensity(full, x, y));
                PredictiveCheck col = PredictiveChecks.PsisLoo(bars.PointwiseLogDensity(full, x, y, content: PointwiseContent.CollapsedPredictive));
                for (int i = 0; i < N; i++)
                {
                    (double reference, double refError) = Reference(full, x, y, i, 1000 * seed + i);
                    double dRet = ret.PerUnit[i] - reference, dCol = col.PerUnit[i] - reference;
                    units.Add(new Unit(shape, seed, i, ret.ParetoK![i], col.ParetoK![i], dRet, dCol,
                                       dRet / Math.Sqrt(ret.MonteCarloErrors![i] * ret.MonteCarloErrors[i] + refError * refError),
                                       dCol / Math.Sqrt(col.MonteCarloErrors![i] * col.MonteCarloErrors[i] + refError * refError), refError));
                }
                Unit[] these = units.Where(u => u.Data == shape && u.Seed == seed).ToArray();
                _out.WriteLine($"{shape} seed {seed}: flagged retained {these.Count(u => u.KRet > 0.7)}, collapsed {these.Count(u => u.KCol > 0.7)}; " +
                               $"rms z retained {Rms(these.Select(u => u.ZRet)):F2}, collapsed {Rms(these.Select(u => u.ZCol)):F2}; " +
                               $"rms Δ retained {Rms(these.Select(u => u.DRet)):F4}, collapsed {Rms(these.Select(u => u.DCol)):F4}; mean ref mcse {these.Average(u => u.RefError):F4}");
                foreach (Unit u in these.Where(u => u.KRet > 0.7 || u.KCol > 0.7 || Math.Abs(u.ZRet) > 3 || Math.Abs(u.ZCol) > 3))
                    _out.WriteLine($"   unit {u.Index,2}: k̂ ret {u.KRet,6:F2} col {u.KCol,6:F2} | Δ ret {u.DRet,8:F4} col {u.DCol,8:F4} | z ret {u.ZRet,6:F2} col {u.ZCol,6:F2} | ref mcse {u.RefError:F4}");
            }

        Unit[] colFlagged = units.Where(u => u.KCol > 0.7).ToArray();
        Unit[] retPassed = units.Where(u => u.KRet <= 0.7).ToArray();
        Unit[] retFlagged = units.Where(u => u.KRet > 0.7).ToArray();
        Unit[] onlyCol = units.Where(u => u.KCol > 0.7 && u.KRet <= 0.7).ToArray();
        double meanOnly = onlyCol.Length > 0 ? onlyCol.Average(u => u.DRet) : double.NaN;
        double seOnly = onlyCol.Length > 1 ? Math.Sqrt(onlyCol.Sum(u => (u.DRet - meanOnly) * (u.DRet - meanOnly)) / (onlyCol.Length - 1) / onlyCol.Length) : double.NaN;
        _out.WriteLine($"pooled: {units.Count} units");
        _out.WriteLine($"(a) collapsed-flagged {colFlagged.Length}: |z| ≤ 3 in {colFlagged.Count(u => Math.Abs(u.ZCol) <= 3)} ({Share(colFlagged, u => Math.Abs(u.ZCol) <= 3):P0}); rms z {Rms(colFlagged.Select(u => u.ZCol)):F2}, rms Δ {Rms(colFlagged.Select(u => u.DCol)):F4}");
        _out.WriteLine($"(b) retained-passed {retPassed.Length}: |z| > 3 in {retPassed.Count(u => Math.Abs(u.ZRet) > 3)} ({Share(retPassed, u => Math.Abs(u.ZRet) > 3):P0}); " +
                       $"flagged only by collapsed {onlyCol.Length}: mean Δ ret {meanOnly:F4} ± {seOnly:F4}");
        _out.WriteLine($"retained-flagged {retFlagged.Length}: rms z retained {Rms(retFlagged.Select(u => u.ZRet)):F2}, collapsed {Rms(retFlagged.Select(u => u.ZCol)):F2}; " +
                       $"rms Δ retained {Rms(retFlagged.Select(u => u.DRet)):F4}, collapsed {Rms(retFlagged.Select(u => u.DCol)):F4}");
        _out.WriteLine($"all units: rms z retained {Rms(units.Select(u => u.ZRet)):F2}, collapsed {Rms(units.Select(u => u.ZCol)):F2}; " +
                       $"rms Δ retained {Rms(units.Select(u => u.DRet)):F4}, collapsed {Rms(units.Select(u => u.DCol)):F4}; |z| > 3 retained {units.Count(u => Math.Abs(u.ZRet) > 3)}, collapsed {units.Count(u => Math.Abs(u.ZCol) > 3)}");
    }

    private static double Rms(IEnumerable<double> v) { double[] a = v.ToArray(); return a.Length == 0 ? double.NaN : Math.Sqrt(a.Average(x => x * x)); }

    private static double Share(Unit[] set, Func<Unit, bool> test) => set.Length == 0 ? double.NaN : set.Count(test) / (double)set.Length;
}
