using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Rng;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;
using static Maths.Samplers.Tests.FlatHistogramToys;

namespace Maths.Samplers.Tests;

/// <summary>
/// The warm start (N17 Phase 1, step 4). WHAM pools a tempering pilot's rung draws into ln g, at the fixed point of
/// the multiple-histogram equations, and SAMC's limit on that density of states is the start. On both toys the start
/// matches θ* from the exact density of states within the pilots' own spread, and learning from it reaches the flat
/// start's accuracy in fewer iterations. The pilot reads a traced ladder's slopes as ℓ₁, refusing a target whose ℓ_b
/// is not linear in b; the window, b_lo and K come from it.
/// </summary>
public sealed class WhamTests
{
    private readonly ITestOutputHelper _out;
    public WhamTests(ITestOutputHelper output) => _out = output;

    private const double Depth = 10.0;
    private static readonly double[] WellLadder = ParallelTempering<double>.GeometricLadder(8, 0.05);

    // A double-well pilot: one ladder from b = 1 down to b_lo = 0.05, traced for the given scans after a burn-in.
    internal static TemperingPilot WellPilot(int seed, int scans = 20_000, int ladders = 1)
    {
        int[] seeds = SeedTree.Derive(seed, ladders);
        var traced = new List<ParallelTempering<double>>();
        foreach (int s in seeds)
        {
            var pt = new ParallelTempering<double>(_ => [new GaussianStep(0.25)], DoubleWell(Depth), 1.0, WellLadder, TemperingPath.Likelihood, s);
            for (int i = 0; i < 2_000; i++) pt.Step();
            pt.TraceRungs();
            for (int i = 0; i < scans; i++) pt.Step();
            traced.Add(pt);
        }
        return TemperingPilot.FromLadders(traced);
    }

    // A Potts pilot: one ladder of eight even rungs from b = 1 to the prior, traced.
    internal static TemperingPilot PottsPilot(int seed, int scans = 20_000)
    {
        Potts lattice = PottsFlatHistogramTests.Lattice;
        double[] ladder = Enumerable.Range(0, 8).Select(i => 1.0 - i / 7.0).ToArray();
        var pt = new ParallelTempering<int[]>(_ => [new PottsFlatHistogramTests.SingleSiteFlip(lattice)], lattice.Target, new int[lattice.Sites],
                                              ladder, TemperingPath.Likelihood, seed);
        for (int i = 0; i < 2_000; i++) pt.Step();
        pt.TraceRungs();
        for (int i = 0; i < scans; i++) pt.Step();
        return TemperingPilot.FromLadders([pt]);
    }

    // The mean over nodes of z², z the pilots' mean departure from the exact over its standard error from their spread,
    // each pilot's ln w centred over the nodes first: about ν/(ν − 2) for ν = pilots − 1 when the spread covers the error.
    private static (double MeanZ2, double Rms) Agreement(IReadOnlyList<MulticanonicalWeights> starts, MulticanonicalWeights exact, bool[] use)
    {
        int nodes = exact.Nodes, p = starts.Count;
        double Centre(IReadOnlyList<double> v)
        {
            double s = 0.0;
            int n = 0;
            for (int k = 0; k < nodes; k++) if (use[k]) { s += v[k]; n++; }
            return s / n;
        }
        double exactCentre = Centre(exact.LogWeights);
        double[] centres = starts.Select(w => Centre(w.LogWeights)).ToArray();
        double z2 = 0.0, rms = 0.0;
        int counted = 0;
        for (int k = 0; k < nodes; k++)
        {
            if (!use[k]) continue;
            double[] d = Enumerable.Range(0, p).Select(i => starts[i].LogWeights[k] - centres[i] - (exact.LogWeights[k] - exactCentre)).ToArray();
            double mean = d.Average();
            double se = Math.Sqrt(d.Sum(v => (v - mean) * (v - mean)) / (p - 1) / p);
            z2 += mean * mean / (se * se);
            rms += d.Sum(v => v * v) / p;
            counted++;
        }
        return (z2 / counted, Math.Sqrt(rms / counted));
    }

    [Fact]
    public void DoubleWell_WarmStart_MatchesTheExactLimit_WithinThePilotsSpread()
    {
        MulticanonicalWeights shape = MulticanonicalWeights.Flat(-30.0, -0.02, 0.05, 31);
        MulticanonicalWeights exact = SamcLearningTests.WellFixedPoint(new FlatHistogramWindow(-30.0, -0.02), 0.05, 31, out bool[] empty);
        int[] seeds = SeedTree.Derive(51, 8);
        var starts = new List<MulticanonicalWeights>();
        var unweighted = new List<MulticanonicalWeights>();
        var inefficiencies = new List<string>();
        foreach (int s in seeds)
        {
            TemperingPilot pilot = WellPilot(s);
            WhamEstimate estimate = Wham.Estimate(pilot);
            Assert.True(estimate.Residual < 1e-8, $"WHAM's fixed point holds to {estimate.Residual:E1}");
            starts.Add(Wham.WarmStart(pilot, shape));
            inefficiencies.Add(string.Join("/", pilot.Inefficiencies.Select(g => g.ToString("F1"))));
            // The same draws with every g = 1, as for independent draws: a reading of what the weights buy.
            var flat = new TemperingPilot(pilot.InverseTemperatures, Enumerable.Range(0, pilot.Rungs).Select(i => pilot.Levels(i).ToArray()).ToArray());
            unweighted.Add(Wham.WarmStart(flat, shape));
        }
        bool[] full = Array.ConvertAll(empty, e => !e);
        (double z2, double rms) = Agreement(starts, exact, full);
        (double z2One, double rmsOne) = Agreement(unweighted, exact, full);
        _out.WriteLine($"[wham] double well, 8 pilots of 8 rungs × 20 000 scans: per-pilot centred RMS to θ* {rms:F3}; mean z² over nodes {z2:F2} (t₇: 1.4); " +
                       $"with every g = 1, {rmsOne:F3} and {z2One:F2}; g per rung, cold first, on the first pilots: {string.Join("; ", inefficiencies.Take(2))}");
        Assert.True(z2 < 3.0, $"the warm starts depart from θ* by more than their spread: mean z² {z2:F2}");
        Assert.True(rms < 0.25, $"a warm start departs from θ* by {rms:F3} RMS");
    }

    [Fact]
    public void Potts_WarmStart_MatchesTheExactLimit_WithinThePilotsSpread()
    {
        MulticanonicalWeights shape = MulticanonicalWeights.Flat(0.0, 27.0, 0.0, 19);
        MulticanonicalWeights exact = PottsFlatHistogramTests.PottsFixedPoint(shape, out bool[] empty);
        int[] seeds = SeedTree.Derive(53, 8);
        var starts = seeds.Select(s => Wham.WarmStart(PottsPilot(s), shape)).ToList();
        (double z2, double rms) = Agreement(starts, exact, Array.ConvertAll(empty, e => !e));
        _out.WriteLine($"[wham] Potts, 8 pilots of 8 rungs × 20 000 scans: per-pilot centred RMS to θ* {rms:F3}; mean z² over levels {z2:F2} (t₇: 1.4)");
        Assert.True(z2 < 3.0, $"the warm starts depart from θ* by more than their spread: mean z² {z2:F2}");
        Assert.True(rms < 0.25, $"a warm start departs from θ* by {rms:F3} RMS");
    }

    // The learned weights' centred RMS to θ* after a learning of the given length, from flat or from a pilot's WHAM start.
    internal static double WellError(long learning, TemperingPilot? pilot, int seed)
    {
        MulticanonicalWeights exact = SamcLearningTests.WellFixedPoint(new FlatHistogramWindow(-30.0, -0.02), 0.05, 31, out _);
        FlatHistogramSpec spec = SamcLearningTests.WellSpec(learning, 1) with { Start = pilot is null ? FlatHistogramStart.Flat : FlatHistogramStart.Wham };
        var result = MulticanonicalRun.Run(spec, DoubleWell(Depth), new WellModel(Depth, 0.25), seed, pilot, parallel: Serial);
        return CentredRms(result.Weights.LogWeights, exact.LogWeights);
    }

    // The mean over three learning seeds of the learned weights' RMS to θ*, from flat or from the pilot.
    internal static double MeanError(Func<long, TemperingPilot?, int, double> error, long learning, TemperingPilot? pilot, int seed)
        => SeedTree.Derive(seed, 3).Average(s => error(learning, pilot, s));

    // The comparison the warm-start tests assert, at any seeds: the flat start's error after its full length and the warm
    // start's after half of it, each averaged over three learning seeds, and the curves around them for the report.
    internal static (double Flat, double WarmHalf, double[] FlatCurve, double[] WarmCurve) WarmAgainstFlat(
        Func<long, TemperingPilot?, int, double> error, TemperingPilot pilot, long length, int seed)
    {
        double[] flat = [MeanError(error, length / 4, null, seed), MeanError(error, length / 2, null, seed), MeanError(error, length, null, seed)];
        double[] warm = [MeanError(error, 0, pilot, seed), MeanError(error, length / 4, pilot, seed), MeanError(error, length / 2, pilot, seed)];
        return (flat[2], warm[2], flat, warm);
    }

    [Fact]
    public void DoubleWell_LearningFromTheWarmStart_ReachesTheFlatStartsAccuracy_InHalfTheIterations()
    {
        // The gain starts at t_w = the pilot's draws over the chains, so the warm start skips the flat start's large-gain
        // discovery; on a toy that discovery is cheap, and the saving is about the pilot's own cost.
        TemperingPilot pilot = WellPilot(57);
        (double flat, double warmHalf, double[] f, double[] w) = WarmAgainstFlat(WellError, pilot, 32_000, 61);
        _out.WriteLine($"[wham] double well, mean centred RMS to θ* over 3 seeds: flat start 8k {f[0]:F3}, 16k {f[1]:F3}, 32k {f[2]:F3}; " +
                       $"warm start (pilot {pilot.Draws:N0} draws, the gain from t = {pilot.Draws / 4:N0}) 0 {w[0]:F3}, 8k {w[1]:F3}, 16k {w[2]:F3}");
        Assert.True(warmHalf < flat, $"16 000 iterations from the warm start ({warmHalf:F3}) should beat 32 000 from flat ({flat:F3})");
    }

    [Fact]
    public void Potts_LearningFromTheWarmStart_ReachesTheFlatStartsAccuracy_InHalfTheIterations()
    {
        TemperingPilot pilot = PottsPilot(59);
        (double flat, double warmHalf, double[] f, double[] w) = WarmAgainstFlat(PottsError, pilot, 100_000, 63);
        _out.WriteLine($"[wham] Potts, mean centred RMS to θ* over 3 seeds: flat start 25k {f[0]:F3}, 50k {f[1]:F3}, 100k {f[2]:F3}; " +
                       $"warm start (pilot {pilot.Draws:N0} draws, the gain from t = {pilot.Draws / 4:N0}) 0 {w[0]:F3}, 25k {w[1]:F3}, 50k {w[2]:F3}");
        Assert.True(warmHalf < flat, $"50 000 iterations from the warm start ({warmHalf:F3}) should beat 100 000 from flat ({flat:F3})");
    }

    internal static double PottsError(long learning, TemperingPilot? pilot, int seed)
    {
        FlatHistogramSpec spec = PottsFlatHistogramTests.PottsSpec(learning, 1, chains: 4) with { Start = pilot is null ? FlatHistogramStart.Flat : FlatHistogramStart.Wham };
        Potts lattice = PottsFlatHistogramTests.Lattice;
        var result = MulticanonicalRun.Run(spec, lattice.Target, new PottsFlatHistogramTests.PottsModel(lattice), seed, pilot, parallel: Serial);
        MulticanonicalWeights exact = PottsFlatHistogramTests.PottsFixedPoint(result.Weights, out bool[] empty);
        return CentredRms(result.Weights.LogWeights, exact.LogWeights, Array.ConvertAll(empty, e => !e));
    }

    [Fact]
    public void ARunFromThePilotAlone_ReadsItsWindowSlopeAndNodes_AndItsEvidenceIsExact()
    {
        TemperingPilot pilot = WellPilot(67, ladders: 2);
        Assert.Equal(2 * 8 * 20_000, pilot.Draws);
        FlatHistogramWindow window = Wham.Window(pilot, 0.01, 1.0);
        int nodes = Wham.Nodes(pilot, window);
        double[] hot = pilot.Levels(pilot.HotRung).ToArray(), cold = pilot.Levels(pilot.ColdRung).ToArray();
        double[] all = Enumerable.Range(0, pilot.Rungs).SelectMany(i => pilot.Levels(i).ToArray()).ToArray();
        Assert.True(window.Low >= all.Min() && window.Low <= hot.OrderBy(v => v).ElementAt(hot.Length / 100));
        Assert.True(window.High <= all.Max() && window.High >= cold.OrderBy(v => v).ElementAt(cold.Length * 99 / 100));
        double sd = Math.Sqrt(cold.Sum(v => (v - cold.Average()) * (v - cold.Average())) / (cold.Length - 1));
        Assert.Equal((int)Math.Ceiling((window.High - window.Low) / sd) + 1, nodes);

        var spec = new FlatHistogramSpec { Chains = 8, LearningIterations = 16_000, ProductionDraws = 40_000, Start = FlatHistogramStart.Wham };
        Assert.True(spec.NeedsPilot);
        Assert.Throws<ArgumentException>(() => MulticanonicalRun.Run(spec, DoubleWell(Depth), new WellModel(Depth, 0.25), 1));
        var result = MulticanonicalRun.Run(spec, DoubleWell(Depth), new WellModel(Depth, 0.25), 69, pilot);
        Assert.Equal((window.Low, window.High, 0.05, nodes), (result.Weights.Low, result.Weights.High, result.Weights.LowSlope, result.Weights.Nodes));
        Assert.Equal(pilot.Draws / 8, result.Learning.GainStart);
        Assert.Equal(10.0 * nodes, result.Learning.GainScale);
        LogEvidenceCurve curve = Reweighting.LogEvidence(result.Draws!, [0.05, 0.5, 1.0]);
        for (int g = 0; g < 3; g++)
        {
            double exact = WellLogEvidence(Depth, curve.InverseTemperatures[g]);
            _out.WriteLine($"[wham] from the pilot alone: window [{window.Low:F2}, {window.High:F3}], K = {nodes}; b = {curve.InverseTemperatures[g]}: " +
                           $"ln Z {curve.LogEvidence[g]:F4} ± {curve.MonteCarloError[g]:F4} against {exact:F4}");
            Assert.True(Math.Abs(curve.LogEvidence[g] - exact) < 4.5 * curve.MonteCarloError[g]);
        }
    }

    [Fact]
    public void Pilot_RefusesALadderWhoseLikelihoodIsNotLinearInB_AndBadRungs()
    {
        // ℓ_b = ln(1 + b·(e^{ℓ₁} − 1)), a collapsed-like likelihood whose slope in b is not ℓ₁.
        ISamplingTarget<double> curved = new CurvedTarget();
        var pt = new ParallelTempering<double>(_ => [new GaussianStep(0.3)], curved, 0.0, [1.0, 0.5, 0.2], TemperingPath.Likelihood, 3);
        pt.TraceRungs();
        for (int i = 0; i < 50; i++) pt.Step();
        Assert.Throws<ArgumentException>(() => TemperingPilot.FromLadders([pt]));

        var untraced = new ParallelTempering<double>(_ => [new GaussianStep(0.3)], DoubleWell(Depth), 1.0, [1.0, 0.5], TemperingPath.Likelihood, 3);
        Assert.Throws<ArgumentException>(() => TemperingPilot.FromLadders([untraced]));
        var powered = new ParallelTempering<double>(_ => [new GaussianStep(0.3)], DoubleWell(Depth), 1.0, [1.0, 0.5], TemperingPath.Powered, 3);
        powered.TraceRungs();
        powered.Step();
        Assert.Throws<ArgumentException>(() => TemperingPilot.FromLadders([powered]));

        Assert.Throws<ArgumentException>(() => new TemperingPilot([0.5, 0.2], [[1.0, 2.0], [1.0, 2.0]]));   // no rung at b = 1
        Assert.Throws<ArgumentException>(() => new TemperingPilot([1.0, 1.0], [[1.0, 2.0], [1.0, 2.0]]));
        Assert.Throws<ArgumentException>(() => new TemperingPilot([1.0, 0.5], [[1.0, 2.0], [double.NaN, 2.0]]));
        Assert.Throws<ArgumentException>(() => new TemperingPilot([1.0, 0.5], [[1.0, 2.0], [1.0]]));
    }

    [Fact]
    public void Wham_OnExactRungSamples_RecoversTheExactFreeEnergies()
    {
        // Rung draws placed at the Potts levels in exact proportion to each rung's law: WHAM's f_i is then ln Z(b_i) up to
        // a constant, to quadrature accuracy in the counts' rounding.
        long[] counts = PottsFlatHistogramTests.Lattice.Counts();
        double[] b = [1.0, 0.6, 0.3, 0.0];
        var levels = new List<double[]>();
        foreach (double bi in b)
        {
            double[] log = counts.Select((c, n) => c > 0 ? Math.Log(c) + bi * 1.5 * n : double.NegativeInfinity).ToArray();
            double z = LogSumExp(log);
            var draws = new List<double>();
            for (int n = 0; n < counts.Length; n++) for (int r = 0; r < (int)Math.Round(1e5 * Math.Exp(log[n] - z)); r++) draws.Add(1.5 * n);
            levels.Add([.. draws]);
        }
        WhamEstimate estimate = Wham.Estimate(new TemperingPilot(b, levels));
        double ExactF(double bi) => LogSumExp(counts.Select((c, n) => c > 0 ? Math.Log(c) + bi * 1.5 * n : double.NegativeInfinity).ToArray());
        for (int i = 0; i < b.Length; i++)
        {
            double expected = ExactF(b[i]) - ExactF(0.0);
            _out.WriteLine($"[wham] exact rungs: f({b[i]}) {estimate.LogNormalizers[i]:F4} against {expected:F4}");
            Assert.True(Math.Abs(estimate.LogNormalizers[i] - expected) < 2e-3);
        }
        Assert.True(estimate.Residual < 1e-8);
    }

    [Fact]
    public void Wham_WeighsEachRungByItsEffectiveDraws_AsFerrenbergAndSwendsen()
    {
        // Exact rung histograms give the exact f under any weights, so unequal inefficiencies leave them exact; one g on
        // every rung scales every weight alike and changes nothing; on sampled draws unequal g move the estimate.
        long[] counts = PottsFlatHistogramTests.Lattice.Counts();
        double[] b = [1.0, 0.6, 0.3, 0.0];
        var exactLevels = new List<double[]>();
        foreach (double bi in b)
        {
            double[] log = counts.Select((c, n) => c > 0 ? Math.Log(c) + bi * 1.5 * n : double.NegativeInfinity).ToArray();
            double z = LogSumExp(log);
            var draws = new List<double>();
            for (int n = 0; n < counts.Length; n++) for (int r = 0; r < (int)Math.Round(1e5 * Math.Exp(log[n] - z)); r++) draws.Add(1.5 * n);
            exactLevels.Add([.. draws]);
        }
        WhamEstimate plain = Wham.Estimate(new TemperingPilot(b, exactLevels));
        WhamEstimate weighted = Wham.Estimate(new TemperingPilot(b, exactLevels, [1.0, 3.0, 1.5, 7.0]));
        for (int i = 0; i < b.Length; i++) Assert.True(Math.Abs(weighted.LogNormalizers[i] - plain.LogNormalizers[i]) < 1e-3);
        Assert.Equal([exactLevels[0].Length / 1.0, exactLevels[1].Length / 3.0, exactLevels[2].Length / 1.5, exactLevels[3].Length / 7.0],
                     weighted.EffectiveDraws);

        TemperingPilot pilot = PottsPilot(71, scans: 4_000);
        double[][] levels = Enumerable.Range(0, pilot.Rungs).Select(i => pilot.Levels(i).ToArray()).ToArray();
        WhamEstimate one = Wham.Estimate(new TemperingPilot(pilot.InverseTemperatures, levels));
        WhamEstimate common = Wham.Estimate(new TemperingPilot(pilot.InverseTemperatures, levels, Enumerable.Repeat(4.0, pilot.Rungs).ToArray()));
        for (int i = 0; i < pilot.Rungs; i++) Assert.Equal(one.LogNormalizers[i], common.LogNormalizers[i], 9);
        double[] uneven = Enumerable.Range(0, pilot.Rungs).Select(i => 1.0 + i).ToArray();
        WhamEstimate moved = Wham.Estimate(new TemperingPilot(pilot.InverseTemperatures, levels, uneven));
        Assert.True(Enumerable.Range(0, pilot.Rungs).Max(i => Math.Abs(moved.LogNormalizers[i] - one.LogNormalizers[i])) > 1e-6);
        Assert.True(moved.Residual < 1e-8);
    }

    [Fact]
    public void FromLadders_ReadsEachRungsInefficiency_FromItsChains()
    {
        int[] seeds = SeedTree.Derive(73, 2);
        var ladders = new List<ParallelTempering<double>>();
        foreach (int s in seeds)
        {
            var pt = new ParallelTempering<double>(_ => [new GaussianStep(0.25)], DoubleWell(Depth), 1.0, WellLadder, TemperingPath.Likelihood, s);
            for (int i = 0; i < 1_000; i++) pt.Step();
            pt.TraceRungs();
            for (int i = 0; i < 5_000; i++) pt.Step();
            ladders.Add(pt);
        }
        TemperingPilot pilot = TemperingPilot.FromLadders(ladders);
        for (int i = 0; i < pilot.Rungs; i++)
        {
            double[][] chains = ladders.Select(l => l.RungTraces![i].Draws(Rjmcmc.RungFunctional.Slope)).ToArray();
            double ess = Ensemble.ChainDiagnostics.MeanEss(chains).Value;
            Assert.Equal(Math.Max(1.0, 10_000 / ess), pilot.Inefficiencies[i], 12);
            Assert.True(pilot.Inefficiencies[i] > 1.0, $"rung {i}'s random-walk draws should be correlated");
        }
        _out.WriteLine($"[wham] two ladders × 5 000 scans: g per rung, cold first, {string.Join(", ", pilot.Inefficiencies.Select(g => g.ToString("F1")))}");

        // Ladders traced for different numbers of scans are refused; so are ill-formed inefficiencies.
        ladders[1].Step();
        Assert.Throws<ArgumentException>(() => TemperingPilot.FromLadders(ladders));
        Assert.Throws<ArgumentException>(() => new TemperingPilot([1.0, 0.5], [[1.0, 2.0], [1.0, 2.0]], [1.0, 0.5]));
        Assert.Throws<ArgumentException>(() => new TemperingPilot([1.0, 0.5], [[1.0, 2.0], [1.0, 2.0]], [1.0, double.NaN]));
        Assert.Throws<ArgumentException>(() => new TemperingPilot([1.0, 0.5], [[1.0, 2.0], [1.0, 2.0]], [1.0]));
        Assert.Equal([1.0, 1.0], new TemperingPilot([1.0, 0.5], [[1.0, 2.0], [1.0, 2.0]]).Inefficiencies);
    }

    // ℓ_b = ln(1 + b(e^{ℓ₁} − 1)) with ℓ₁ = −x²/2 and a flat prior on [−3, 3]: not linear in b.
    private sealed class CurvedTarget : ISamplingTarget<double>
    {
        public ITargetEvaluation Evaluate(double x) => new Curved(Math.Abs(x) <= 3.0 ? 0.0 : double.NegativeInfinity, -0.5 * x * x);

        private sealed class Curved(double prior, double l1) : ITargetEvaluation
        {
            public double LogPrior => prior;
            public double LogLikelihood(double b) => Math.Log(1.0 + b * (Math.Exp(l1) - 1.0));
            public double LogLikelihoodSlope(double b) => (Math.Exp(l1) - 1.0) / (1.0 + b * (Math.Exp(l1) - 1.0));
        }
    }
}
