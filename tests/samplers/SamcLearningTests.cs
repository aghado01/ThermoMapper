using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.Rng;
using Maths.Samplers.FlatHistogram;
using Xunit;
using Xunit.Abstractions;
using static Maths.Samplers.Tests.FlatHistogramToys;

namespace Maths.Samplers.Tests;

/// <summary>
/// SAMC's learning and the multicanonical run (N17 Phase 1, step 2). The gain is Liang et al.'s t₀/max(t₀, t); a visit's
/// indicator is spread over its segment's two nodes by their hat functions, conserving Σ ln w; the limit is the θ* with
/// flat conditional hat occupancy, exact when ln g is piecewise linear on the nodes. On the double well, the learned
/// weights approach θ* from the prior's quadrature, closer the longer the learning; production's histogram is flat; the
/// shared weights give the same run under any parallelism bound; and the diagnostics flag a chain that never crosses.
/// </summary>
public sealed class SamcLearningTests
{
    private readonly ITestOutputHelper _out;
    public SamcLearningTests(ITestOutputHelper output) => _out = output;

    private const double Depth = 10.0;
    private static readonly FlatHistogramWindow WellWindow = new(-30.0, -0.02);
    private const double WellLowSlope = 0.05;
    private const int WellNodes = 31;

    [Fact]
    public void Gain_IsLiangsSequence()
    {
        Assert.Equal(1.0, SamcLearning.Gain(10.0, 1));
        Assert.Equal(1.0, SamcLearning.Gain(10.0, 10));
        Assert.Equal(10.0 / 11.0, SamcLearning.Gain(10.0, 11));
        Assert.Equal(1e-3, SamcLearning.Gain(10.0, 10_000));
        Assert.Throws<ArgumentOutOfRangeException>(() => SamcLearning.Gain(10.0, 0));
        var learning = new SamcLearning(MulticanonicalWeights.Flat(0.0, 1.0, 0.5, 3), gainScale: 4.0, gainStart: 7);
        Assert.Equal(0.5, learning.NextGain);   // t = 8
        learning.Update([0.5]);
        Assert.Equal(4.0 / 9.0, learning.NextGain);
        Assert.Throws<ArgumentOutOfRangeException>(() => new SamcLearning(MulticanonicalWeights.Flat(0.0, 1.0, 0.5, 3), 0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SamcLearning(MulticanonicalWeights.Flat(0.0, 1.0, 0.5, 3), 2.0, -1));
    }

    [Fact]
    public void Update_SpreadsEachVisitOverItsNodes_AveragesTheChains_AndConservesTheSum()
    {
        // Five nodes on [0, 4]: π = (1/8, 1/4, 1/4, 1/4, 1/8). One iteration at γ = 1 with three chains: a visit a quarter
        // of the way from node 1 to node 2, one at node 4 exactly, and one below the window, which does not update.
        var learning = new SamcLearning(MulticanonicalWeights.Flat(0.0, 4.0, 0.5, 5), gainScale: 10.0);
        MulticanonicalWeights w = learning.Update([1.25, 4.0, -3.0]);
        double[] visits = [0.0, 0.75, 0.25, 0.0, 1.0];
        double[] pi = [0.125, 0.25, 0.25, 0.25, 0.125];
        double sum = 0.0;
        for (int k = 0; k < 5; k++)
        {
            Assert.Equal(-(visits[k] - 2 * pi[k]) / 3.0, w.LogWeights[k], 14);
            sum += w.LogWeights[k];
        }
        Assert.Equal(0.0, sum, 14);
        Assert.Equal(1, learning.Iterations);

        // Only outside visits: no update, but the iteration counts.
        MulticanonicalWeights same = learning.Update([-1.0, 9.0]);
        Assert.Same(w, same);
        Assert.Equal(2, learning.Iterations);
    }

    [Fact]
    public void FixedPoint_RecoversAPiecewiseLinearLnG_UpToAConstant()
    {
        // A density of states whose ln g is piecewise linear on the nodes, given as a fine quadrature of ℓ, plus mass outside
        // the window, which plays no part.
        var shape = MulticanonicalWeights.Flat(-6.0, 2.0, 0.3, 9);
        double[] lnG = [3.0, 1.5, 2.5, 0.0, -1.0, -0.5, -4.0, -2.0, -6.0];
        const int cells = 80_000;
        var levels = new double[cells + 2];
        var masses = new double[cells + 2];
        double h = 8.0 / cells;
        var interpolant = shape.WithLogWeights(lnG);
        for (int j = 0; j < cells; j++)
        {
            levels[j] = -6.0 + (j + 0.5) * h;
            masses[j] = interpolant.LogWeight(levels[j]) + Math.Log(h);
        }
        (levels[cells], masses[cells]) = (-20.0, 5.0);
        (levels[cells + 1], masses[cells + 1]) = (7.0, 5.0);

        MulticanonicalWeights fixedPoint = SamcLearning.FixedPoint(shape, levels, masses, out bool[] empty);
        Assert.DoesNotContain(true, empty);
        double[] minusLnG = Array.ConvertAll(lnG, v => -v);
        double rms = CentredRms(fixedPoint.LogWeights, minusLnG);
        _out.WriteLine($"[fixed point] piecewise-linear ln g recovered to {rms:E2}");
        Assert.True(rms < 1e-6, $"θ* departs from ln g by {rms:E2}");
    }

    [Fact]
    public void FixedPoint_HandlesADensityOfStatesWiderThanADoublesRange()
    {
        // ln g = −3ℓ on [0, 1000], fine levels: its masses span 3 000 nats, past the 745 of a double's exponent, so a start
        // summed in linear space would underflow its far nodes. θ* is ln g itself, piecewise linear on any nodes.
        var shape = MulticanonicalWeights.Flat(0.0, 1000.0, 0.5, 14);
        const int cells = 100_000;
        var levels = new double[cells];
        var masses = new double[cells];
        for (int j = 0; j < cells; j++)
        {
            levels[j] = (j + 0.5) * 1000.0 / cells;
            masses[j] = -3.0 * levels[j];
        }
        MulticanonicalWeights fixedPoint = SamcLearning.FixedPoint(shape, levels, masses, out _);
        double[] expected = Enumerable.Range(0, 14).Select(k => 3.0 * shape.Node(k)).ToArray();
        double rms = CentredRms(fixedPoint.LogWeights, expected);
        _out.WriteLine($"[fixed point] ln g spanning 3 000 nats recovered to {rms:E2}");
        Assert.True(rms < 1e-4, $"θ* departs from ln g by {rms:E2}");
    }

    [Fact]
    public void FixedPoint_OnLevelsAtTheNodes_IsLnGLessLnPi_AndSkipsEmptyNodes()
    {
        // A lattice energy on the nodes: Liang et al.'s θ = C + ln g − ln π. Node 3 holds no mass, so the others take
        // π_k + d, d = π_3/(K − 1), and node 3 sits on the line through its neighbours.
        var shape = MulticanonicalWeights.Flat(0.0, 5.0, 0.0, 6);
        double[] lnG = [0.0, 2.0, 3.0, double.NegativeInfinity, 1.0, -2.0];
        double[] levels = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        MulticanonicalWeights fixedPoint = SamcLearning.FixedPoint(shape, levels, lnG, out bool[] empty);
        Assert.Equal([false, false, false, true, false, false], empty);
        double[] pi = shape.FlatOccupancy();
        double d = pi[3] / 5.0;
        var expected = new double[6];
        bool[] full = [true, true, true, false, true, true];
        for (int k = 0; k < 6; k++) expected[k] = full[k] ? -(lnG[k] - Math.Log(pi[k] + d)) : 0.0;
        Assert.True(CentredRms(fixedPoint.LogWeights, expected, full) < 1e-9);
        Assert.Equal(0.5 * (fixedPoint.LogWeights[2] + fixedPoint.LogWeights[4]), fixedPoint.LogWeights[3], 9);
    }

    [Fact]
    public void FixedPoint_RefusesAMeasureOutsideTheWindow()
    {
        var shape = MulticanonicalWeights.Flat(0.0, 1.0, 0.5, 3);
        Assert.Throws<ArgumentException>(() => SamcLearning.FixedPoint(shape, [2.0, 3.0], [0.0, 0.0], out _));
        Assert.Throws<ArgumentException>(() => SamcLearning.FixedPoint(shape, [0.5], [double.NaN], out _));
        Assert.Throws<ArgumentException>(() => SamcLearning.FixedPoint(shape, [0.5], [0.0, 1.0], out _));
    }

    // The double well's θ* over a window, from a 400 000-cell quadrature of the prior.
    internal static MulticanonicalWeights WellFixedPoint(FlatHistogramWindow window, double lowSlope, int nodes, out bool[] empty)
    {
        (_, double[] levels, double[] masses) = WellQuadrature(Depth, 400_000);
        return SamcLearning.FixedPoint(MulticanonicalWeights.Flat(window.Low, window.High, lowSlope, nodes), levels, masses, out empty);
    }

    internal static FlatHistogramSpec WellSpec(long learning, int production, int chains = 4, FlatHistogramWindow? window = null) => new()
    {
        Chains = chains,
        LearningIterations = learning,
        ProductionDraws = production,
        Window = window ?? WellWindow,
        LowSlope = WellLowSlope,
        Nodes = WellNodes,
    };

    // The learned weights' centred RMS departure from θ* after a learning of the given length.
    internal static double WellLearningError(long learning, int seed)
    {
        MulticanonicalWeights exact = WellFixedPoint(WellWindow, WellLowSlope, WellNodes, out _);
        var result = MulticanonicalRun.Run(WellSpec(learning, production: 1), DoubleWell(Depth), new WellModel(Depth, 0.25), seed, parallel: Serial);
        return CentredRms(result.Weights.LogWeights, exact.LogWeights);
    }

    [Fact]
    public void DoubleWell_LearnsSamcsLimit_CloserTheLongerItLearns()
    {
        double e1 = WellLearningError(4_000, seed: 21), e2 = WellLearningError(16_000, seed: 21), e3 = WellLearningError(64_000, seed: 21);
        _out.WriteLine($"[samc] double well, centred RMS of ln w − θ*: 4k {e1:F3}, 16k {e2:F3}, 64k {e3:F3}");
        Assert.True(e3 < e2 && e2 < e1, $"the error should fall with the learning's length: {e1:F3}, {e2:F3}, {e3:F3}");
        Assert.True(e3 < 0.15, $"after 64 000 iterations the learned ln w departs from θ* by {e3:F3}");
    }

    [Fact]
    public void DoubleWell_LearnedWeightsAreLnG_AwayFromItsSingularities()
    {
        // g has integrable singularities at ℓ = −c (the barrier's top, x = 0) and at ℓ → 0 (the wells' floors), where
        // the piecewise-linear θ* departs from ln g at the nodes. Away from them, within 2.5 spacings of −c and 1.5 of
        // the top, both θ* and the learned weights are ln g up to a constant.
        MulticanonicalWeights exact = WellFixedPoint(WellWindow, WellLowSlope, WellNodes, out _);
        var result = MulticanonicalRun.Run(WellSpec(64_000, 1), DoubleWell(Depth), new WellModel(Depth, 0.25), masterSeed: 5, parallel: Serial);
        var lnG = new double[WellNodes];
        var smooth = new bool[WellNodes];
        for (int k = 0; k < WellNodes; k++)
        {
            double node = exact.Node(k);
            lnG[k] = -WellLogDensityOfStates(node, Depth);   // as ln w
            smooth[k] = Math.Abs(node + Depth) > 2.5 * exact.Spacing && node < -1.5 * exact.Spacing;
        }
        double projection = CentredRms(exact.LogWeights, lnG, smooth);
        double learned = CentredRms(result.Weights.LogWeights, lnG, smooth);
        _out.WriteLine($"[samc] against the closed-form ln g at {smooth.Count(b => b)} of {WellNodes} nodes: θ* {projection:F3}, learned {learned:F3}");
        Assert.True(projection < 0.1, $"θ* departs from ln g by {projection:F3} away from the singularities");
        Assert.True(learned < 0.15, $"the learned ln w departs from ln g by {learned:F3} away from the singularities");
    }

    [Fact]
    public void DoubleWell_ProductionHistogramIsFlat_OnTheFrozenWeights_AndEveryChainTunnels()
    {
        // Two readings: the learned weights' own occupancy, exact by quadrature, against flat (the learning's accuracy),
        // and production's visits against that occupancy (production samples π_w on the frozen weights). Production's
        // per-node error runs at about one over the root of the round trips, some 600 here.
        var result = MulticanonicalRun.Run(WellSpec(64_000, 40_000), DoubleWell(Depth), new WellModel(Depth, 0.25), masterSeed: 5, parallel: Serial);
        FlatHistogramDiagnostics d = result.ProductionDiagnostics!;
        double[] predicted = WellOccupancy(result.Weights, Depth);
        double[] pi = result.Weights.FlatOccupancy();
        double learningFlatness = 0.0, samplingError = 0.0;
        for (int k = 0; k < WellNodes; k++)
        {
            learningFlatness = Math.Max(learningFlatness, Math.Abs(predicted[k] / pi[k] - 1.0));
            samplingError = Math.Max(samplingError, Math.Abs(d.Occupancy[k] / d.InWindow / predicted[k] - 1.0));
        }
        _out.WriteLine($"[samc] production: learned weights' exact flatness {learningFlatness:F3}, production against it {samplingError:F3}, " +
                       $"histogram's flatness {d.MaxFlatnessDeviation:F3}; in {d.InWindow} below {d.BelowWindow} above {d.AboveWindow}; " +
                       $"round trips {string.Join(",", d.RoundTrips)}; verdict {result.Production!.Verdict.Outcome}");
        Assert.True(learningFlatness < 0.2, $"the learned weights' occupancy departs from flat by {learningFlatness:F3}");
        Assert.True(samplingError < 0.2, $"production's occupancy departs from the frozen weights' by {samplingError:F3}");
        Assert.Empty(d.ChainsWithoutRoundTrip);
        Assert.True(result.Complete);
        Assert.Equal(40_000, result.Draws!.Draws);
        Assert.Same(result.Weights, result.Draws.Weights);
        foreach (IMulticanonicalChain<double> chain in result.Chains) Assert.Same(result.Weights, chain.CurrentEvaluation.Weights);
        Assert.Equal(result.Learning.Weights, result.Weights);
        Assert.Equal(64_000, result.Learning.Iterations);
    }

    [Fact]
    public void SharedWeights_GiveTheSameRun_UnderAnyParallelismBound()
    {
        FlatHistogramSpec spec = WellSpec(3_000, 500);
        MulticanonicalResult<double> Run(int bound) => MulticanonicalRun.Run(spec, DoubleWell(Depth), new WellModel(Depth, 0.25), masterSeed: 9,
            parallel: new ParallelOptions { MaxDegreeOfParallelism = bound });
        MulticanonicalResult<double> one = Run(1), two = Run(2), all = Run(-1);
        foreach (MulticanonicalResult<double> other in new[] { two, all })
        {
            Assert.Equal(one.Weights.LogWeights, other.Weights.LogWeights);
            for (int c = 0; c < spec.Chains; c++) Assert.Equal(one.Draws!.Levels(c).ToArray(), other.Draws!.Levels(c).ToArray());
            Assert.Equal(one.Production!.FunctionalMean, other.Production!.FunctionalMean);
        }
    }

    [Fact]
    public void OneChain_LearnsAlone()
    {
        MulticanonicalWeights exact = WellFixedPoint(WellWindow, WellLowSlope, WellNodes, out _);
        var result = MulticanonicalRun.Run(WellSpec(64_000, 1_000, chains: 1), DoubleWell(Depth), new WellModel(Depth, 0.25), masterSeed: 3, parallel: Serial);
        double error = CentredRms(result.Weights.LogWeights, exact.LogWeights);
        _out.WriteLine($"[samc] one chain: centred RMS {error:F3}; verdict {result.Production!.Verdict.Outcome}");
        Assert.True(error < 0.3, $"one chain's learned ln w departs from θ* by {error:F3}");
        Assert.Single(result.Chains);
    }

    [Fact]
    public void AChainThatNeverCrosses_IsFlagged()
    {
        int[] seeds = SeedTree.Derive(17, 2 * 4);
        var model = new WellModel(Depth, 0.25, stuckSeed: seeds[2]);
        var result = MulticanonicalRun.Run(WellSpec(4_000, 4_000), DoubleWell(Depth), model, masterSeed: 17, parallel: Serial);
        _out.WriteLine($"[samc] stuck chain: learning round trips {string.Join(",", result.LearningDiagnostics.RoundTrips)}, " +
                       $"production {string.Join(",", result.ProductionDiagnostics!.RoundTrips)}");
        Assert.Equal([2], result.ProductionDiagnostics.ChainsWithoutRoundTrip);
        Assert.Equal([2], result.LearningDiagnostics.ChainsWithoutRoundTrip);
    }

    [Fact]
    public void AWindowTooNarrow_IsLeftByTheExtensions_AndLearnsItsOwnLimit()
    {
        // The window covers only part of the range the chain reaches: ℓ₁ in [−12, −2]. Visits below and above take the
        // Boltzmann extensions, the learning ignores them, and the weights inside reach this window's θ*.
        var narrow = new FlatHistogramWindow(-12.0, -2.0);
        MulticanonicalWeights exact = WellFixedPoint(narrow, WellLowSlope, 11, out _);
        FlatHistogramSpec spec = WellSpec(64_000, 20_000, window: narrow) with { Nodes = 11 };
        var result = MulticanonicalRun.Run(spec, DoubleWell(Depth), new WellModel(Depth, 0.25), masterSeed: 13, parallel: Serial);
        double error = CentredRms(result.Weights.LogWeights, exact.LogWeights);
        FlatHistogramDiagnostics d = result.ProductionDiagnostics!;
        _out.WriteLine($"[samc] narrow window: centred RMS {error:F3}; production in {d.InWindow} below {d.BelowWindow} above {d.AboveWindow}, " +
                       $"flatness {d.MaxFlatnessDeviation:F3}");
        Assert.True(d.BelowWindow > 0 && d.AboveWindow > 0, "the chains should leave a narrow window on both sides");
        Assert.True(error < 0.15, $"the narrow window's learned ln w departs from θ* by {error:F3}");
        Assert.True(d.MaxFlatnessDeviation < 0.15);
    }

    [Fact]
    public void AWindowTooWide_LeavesItsEmptyNodesUnbounded_AndLearnsTheRest()
    {
        // ℓ₁ ≤ 0 everywhere, so a window reaching up to +6 has empty nodes at its top, those at 1 to 6, whose hats see no
        // mass (node 0's sees ℓ just below 0). Their ln w rises without bound while the nodes with mass reach θ* with
        // Liang et al.'s d.
        var wide = new FlatHistogramWindow(-30.0, 6.0);
        MulticanonicalWeights exact = WellFixedPoint(wide, WellLowSlope, 37, out bool[] empty);
        FlatHistogramSpec spec = WellSpec(64_000, 20_000, window: wide) with { Nodes = 37 };
        var result = MulticanonicalRun.Run(spec, DoubleWell(Depth), new WellModel(Depth, 0.25), masterSeed: 19, parallel: Serial);
        bool[] full = Array.ConvertAll(empty, e => !e);
        double error = CentredRms(result.Weights.LogWeights, exact.LogWeights, full);
        IReadOnlyList<double> learned = result.Weights.LogWeights;
        _out.WriteLine($"[samc] wide window: {Array.FindAll(empty, e => e).Length} empty nodes; centred RMS over the rest {error:F3}; " +
                       $"ln w at the empty nodes {string.Join(", ", learned.Skip(31).Select(v => v.ToString("F1")))}, at node 0 {learned[30]:F1}");
        Assert.Equal(Enumerable.Range(0, 37).Select(k => k > 30).ToArray(), empty);
        Assert.True(error < 0.15, $"the wide window's populated nodes depart from θ* by {error:F3}");
        for (int k = 31; k < 37; k++) Assert.True(learned[k] > learned[30] + 5.0, $"empty node {k}'s ln w should have risen");
        Assert.True(result.ProductionDiagnostics!.MaxFlatnessDeviation > 0.99, "the empty nodes leave the histogram unflat");
    }

    [Fact]
    public void RoundTrips_CountLowestBinToHighestAndBack_AndTheTraceHoldsTheLearning()
    {
        // A scripted chain walks the levels below; the window [0, 10] has 11 nodes, so the lowest bin is ℓ < 0.5 and the
        // highest ℓ > 9.5. It starts high (no trip yet), then goes low, high, low (one), high, middle, low (two),
        // middle, high.
        double[] script = [9.8, 5.0, 0.2, 4.0, 9.9, 3.0, -1.0, 12.0, 6.0, 0.4, 5.0, 9.7];
        var spec = new FlatHistogramSpec
        {
            Chains = 1, LearningIterations = script.Length, ProductionDraws = 1,
            Window = new FlatHistogramWindow(0.0, 10.0), LowSlope = 0.5, Nodes = 11, GainScale = 4.0,
        };
        var result = MulticanonicalRun.Run(spec, SamplingTarget.Separable<double>(_ => 0.0, x => x), new ScriptedModel(script), masterSeed: 1, parallel: Serial);
        FlatHistogramDiagnostics d = result.LearningDiagnostics;
        Assert.Equal([2L], d.RoundTrips);
        Assert.Equal(1, d.BelowWindow);
        Assert.Equal(1, d.AboveWindow);
        Assert.Equal(script.Length - 2, d.InWindow);
        double occupied = 0.0;
        foreach (double o in d.Occupancy) occupied += o;
        Assert.Equal(d.InWindow, occupied, 12);

        // The trace: the start, each iteration (a hundredth of 12 rounds to 1), the end once.
        Assert.Equal(script.Length + 1, result.WeightTrace.Count);
        Assert.Equal(0, result.WeightTrace[0].Iteration);
        Assert.Equal(new double[11], result.WeightTrace[0].Weights.LogWeights);
        Assert.Same(result.Weights, result.WeightTrace[^1].Weights);
        Assert.Equal(script.Length, result.WeightTrace[^1].Iteration);
    }

    [Fact]
    public void MulticanonicalChain_RefusesAChainOffTheWrappedTarget()
    {
        ISamplingTarget<double> inner = DoubleWell(Depth);
        var target = new MulticanonicalTarget<double>(inner, MulticanonicalWeights.Flat(-30.0, -0.02, 0.05, 31));
        var hot = new Rjmcmc.ReversibleJumpChain<double>([new GaussianStep(0.1)], target, 1.0, new Xoshiro256PlusPlus(1),
                                                          new Temperature(1.0, TemperingPath.Powered));
        Assert.Throws<ArgumentException>(() => new MulticanonicalChain<double>(hot, target, (in double _, Span<double> _) => { }));
        var plain = new Rjmcmc.ReversibleJumpChain<double>([new GaussianStep(0.1)], inner, 1.0, new Xoshiro256PlusPlus(1));
        Assert.Throws<ArgumentException>(() => new MulticanonicalChain<double>(plain, target, (in double _, Span<double> _) => { }));
    }

    // A chain that visits a fixed script of levels, ignoring the weights: for the diagnostics' bookkeeping.
    private sealed class ScriptedModel(double[] script) : IMulticanonicalModel<double>
    {
        public int FunctionalDim => 1;
        public int MonitoredDim => 0;
        public IMulticanonicalChain<double> StartChain(MulticanonicalTarget<double> target, int chainSeed, int readoutSeed) => new Scripted(script, target);
    }

    private sealed class Scripted(double[] script, MulticanonicalTarget<double> target) : IMulticanonicalChain<double>
    {
        private int _next;
        public MulticanonicalEvaluation CurrentEvaluation { get; private set; } = target.Evaluate(0.0);
        public long Accepted => _next;
        public long Attempts => _next;
        public double Step()
        {
            double x = script[_next++ % script.Length];
            CurrentEvaluation = target.Evaluate(x);
            return x;
        }
        public void Reweigh() => CurrentEvaluation = CurrentEvaluation.Reweighed(target.Weights);
        public void Accumulate(in double draw, Span<double> functionals) => functionals[0] = draw;
    }
}
