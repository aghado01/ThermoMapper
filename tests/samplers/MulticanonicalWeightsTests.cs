using System;
using Maths.Rng;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// The flat-histogram engine's weights and stage density (N17 Phase 1, step 1). ln w is piecewise linear between the
/// nodes and continues outside the window with slope b_lo below and 1 above, continuous at the nodes and at both edges.
/// The wrapped target's density is the inner prior plus ln w(ℓ₁), −∞ off the support with no likelihood work, and its
/// log likelihood is refused at b ≠ 1. A chain run on the wrapper samples π_w ∝ prior·w(ℓ₁) exactly.
/// </summary>
public sealed class MulticanonicalWeightsTests
{
    private static readonly double[] Values = [0.5, -1.25, 2.0, 0.0, -3.5];
    private static readonly MulticanonicalWeights W = new(low: -10.0, high: -2.0, lowSlope: 0.25, Values);

    [Fact]
    public void Nodes_HoldTheirValues_AndTheInterpolantIsLinearBetween()
    {
        Assert.Equal(5, W.Nodes);
        Assert.Equal(2.0, W.Spacing);
        for (int k = 0; k < W.Nodes; k++)
        {
            Assert.Equal(-10.0 + 2.0 * k, W.Node(k));
            Assert.Equal(Values[k], W.LogWeight(W.Node(k)), 12);
        }
        // A quarter of the way from node 1 to node 2.
        Assert.Equal(0.75 * -1.25 + 0.25 * 2.0, W.LogWeight(-7.5), 12);
        Assert.True(W.TryLocate(-7.5, out int segment, out double fraction));
        Assert.Equal(1, segment);
        Assert.Equal(0.25, fraction, 12);
        Assert.True(W.TryLocate(-2.0, out segment, out fraction));   // ℓ_hi: the last segment at u = 1
        Assert.Equal(3, segment);
        Assert.Equal(1.0, fraction);
        Assert.False(W.TryLocate(-10.0 - 1e-12, out _, out _));
        Assert.False(W.TryLocate(-2.0 + 1e-12, out _, out _));
        Assert.False(W.TryLocate(double.NaN, out _, out _));
    }

    [Fact]
    public void LogWeight_IsContinuous_AtEveryNode_AndAtBothEdges()
    {
        const double h = 1e-9;
        for (int k = 0; k < W.Nodes; k++)
        {
            double at = W.LogWeight(W.Node(k));
            Assert.Equal(at, W.LogWeight(W.Node(k) - h), 6);
            Assert.Equal(at, W.LogWeight(W.Node(k) + h), 6);
        }
        // The Boltzmann extensions: slope b_lo below the window, slope 1 above it.
        Assert.Equal(Values[0] + 0.25 * -5.0, W.LogWeight(-15.0), 12);
        Assert.Equal(Values[^1] + 3.0, W.LogWeight(1.0), 12);
        double belowSlope = (W.LogWeight(-10.0 - 1e-3) - W.LogWeight(-10.0 - 2e-3)) / 1e-3;
        double aboveSlope = (W.LogWeight(-2.0 + 2e-3) - W.LogWeight(-2.0 + 1e-3)) / 1e-3;
        Assert.Equal(0.25, belowSlope, 6);
        Assert.Equal(1.0, aboveSlope, 6);
    }

    [Fact]
    public void ZeroLikelihood_IsMinusInfinity_UnlessTheLowExtensionIsFlat()
    {
        Assert.Equal(double.NegativeInfinity, W.LogWeight(double.NegativeInfinity));
        var flatBelow = new MulticanonicalWeights(-10.0, -2.0, 0.0, Values);
        Assert.Equal(Values[0], flatBelow.LogWeight(double.NegativeInfinity));
        Assert.Equal(Values[0], flatBelow.LogWeight(-1e300));
        Assert.Equal(double.PositiveInfinity, W.LogWeight(double.PositiveInfinity));
        Assert.True(double.IsNaN(W.LogWeight(double.NaN)));
    }

    [Fact]
    public void FlatOccupancy_GivesTheEdgeNodesHalfAShare_AndSumsToOne()
    {
        double[] pi = W.FlatOccupancy();
        Assert.Equal([0.125, 0.25, 0.25, 0.25, 0.125], pi);
        double[] two = MulticanonicalWeights.Flat(0.0, 1.0, 0.5, 2).FlatOccupancy();
        Assert.Equal([0.5, 0.5], two);
    }

    [Fact]
    public void Construction_RefusesAnIllFormedWindow()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new MulticanonicalWeights(0.0, 0.0, 0.5, [0.0, 0.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new MulticanonicalWeights(1.0, 0.0, 0.5, [0.0, 0.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new MulticanonicalWeights(double.NegativeInfinity, 0.0, 0.5, [0.0, 0.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new MulticanonicalWeights(0.0, 1.0, -0.1, [0.0, 0.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new MulticanonicalWeights(0.0, 1.0, 1.5, [0.0, 0.0]));
        Assert.Throws<ArgumentException>(() => new MulticanonicalWeights(0.0, 1.0, 0.5, [0.0]));
        Assert.Throws<ArgumentException>(() => new MulticanonicalWeights(0.0, 1.0, 0.5, [0.0, double.NaN]));
        Assert.Throws<ArgumentOutOfRangeException>(() => MulticanonicalWeights.Flat(0.0, 1.0, 0.5, 1));
        Assert.Throws<ArgumentException>(() => W.WithLogWeights([0.0, 1.0]));
        Assert.Throws<ArgumentOutOfRangeException>(() => W.Node(5));
    }

    [Fact]
    public void WithLogWeights_KeepsTheWindow_AndCopiesTheValues()
    {
        double[] values = [1.0, 2.0, 3.0, 4.0, 5.0];
        MulticanonicalWeights other = W.WithLogWeights(values);
        values[0] = 99.0;
        Assert.Equal(1.0, other.LogWeights[0]);
        Assert.Equal((W.Low, W.High, W.LowSlope, W.Nodes), (other.Low, other.High, other.LowSlope, other.Nodes));
    }

    [Fact]
    public void Wrapper_IsTheInnerPriorPlusLnW_AndRefusesOtherTemperatures()
    {
        // A separable inner target: prior −x²/2 − 1 on x ≥ 0, ℓ₁ = −3x.
        int likelihoodCalls = 0;
        ISamplingTarget<double> inner = SamplingTarget.Separable<double>(
            x => x < 0.0 ? double.NegativeInfinity : -0.5 * x * x - 1.0, x => { likelihoodCalls++; return -3.0 * x; });
        var target = new MulticanonicalTarget<double>(inner, W);

        foreach (double x in new[] { 0.0, 0.4, 1.7, 2.5, 3.3, 5.0 })
        {
            MulticanonicalEvaluation e = target.Evaluate(x);
            Assert.Equal(-0.5 * x * x - 1.0, e.LogPrior);
            Assert.Equal(-3.0 * x, e.UntemperedLogLikelihood);
            Assert.Equal(W.LogWeight(-3.0 * x), e.LogWeight);
            Assert.Equal(e.LogPrior + W.LogWeight(-3.0 * x), Temperature.Posterior.LogDensity(e), 12);
            Assert.Same(W, e.Weights);
        }

        int before = likelihoodCalls;
        MulticanonicalEvaluation off = target.Evaluate(-1.0);
        Assert.Equal(before, likelihoodCalls);   // no likelihood work off the support
        Assert.Equal(double.NegativeInfinity, off.LogPrior);
        Assert.Equal(double.NegativeInfinity, off.UntemperedLogLikelihood);
        Assert.Equal(double.NegativeInfinity, off.LogWeight);
        Assert.Equal(double.NegativeInfinity, Temperature.Posterior.LogDensity(off));

        MulticanonicalEvaluation on = target.Evaluate(1.0);
        Assert.Throws<NotSupportedException>(() => on.LogLikelihood(0.5));
        Assert.Throws<NotSupportedException>(() => on.LogLikelihood(0.0));
        Assert.Throws<NotSupportedException>(() => new Temperature(0.5, TemperingPath.Likelihood).LogDensity(on));
        Assert.Throws<NotSupportedException>(() => ((ITargetEvaluation)on).LogLikelihoodSlope(1.0));
    }

    [Fact]
    public void Reweighed_ReadsTheCachedLevel_WithoutEvaluatingAgain()
    {
        int evaluations = 0;
        ISamplingTarget<double> inner = SamplingTarget.Separable<double>(_ => 0.0, x => { evaluations++; return -x; });
        var target = new MulticanonicalTarget<double>(inner, W);
        MulticanonicalEvaluation e = target.Evaluate(6.5);
        Assert.Equal(1, evaluations);
        MulticanonicalWeights other = W.WithLogWeights([0.0, 1.0, 2.0, 3.0, 4.0]);
        MulticanonicalEvaluation r = e.Reweighed(other);
        Assert.Equal(1, evaluations);
        Assert.Equal(-6.5, r.UntemperedLogLikelihood);
        Assert.Equal(other.LogWeight(-6.5), r.LogWeight);
        Assert.Same(e.Inner, r.Inner);
        Assert.Same(e, e.Reweighed(W));
    }

    [Fact]
    public void AChainOnTheWrapper_SamplesPriorTimesW()
    {
        // Ten states s = 0..9 under a prior ∝ (s + 1) with ℓ₁(s) = −1.5·s, a nearest-neighbour walk with wrap-around, and
        // weights whose window covers part of the range so that both extensions are visited.
        const int states = 10;
        ISamplingTarget<int> inner = SamplingTarget.Separable<int>(
            s => s is < 0 or >= states ? double.NegativeInfinity : Math.Log(s + 1.0), s => -1.5 * s);
        var weights = new MulticanonicalWeights(-10.0, -3.0, 0.4, [1.0, 4.0, 2.5, 6.0]);
        var target = new MulticanonicalTarget<int>(inner, weights);
        var chain = new ReversibleJumpChain<int>([new Ring(states)], target, 0, new Xoshiro256PlusPlus(11));

        var exact = new double[states];
        double total = 0.0;
        for (int s = 0; s < states; s++) total += exact[s] = (s + 1.0) * Math.Exp(weights.LogWeight(-1.5 * s));
        for (int s = 0; s < states; s++) exact[s] /= total;

        const int draws = 400_000;
        var counts = new int[states];
        for (int i = 0; i < 2_000; i++) chain.Step();
        for (int i = 0; i < draws; i++) counts[chain.Step()]++;
        for (int s = 0; s < states; s++)
        {
            double p = (double)counts[s] / draws;
            Assert.True(Math.Abs(p - exact[s]) < 6.0 * Math.Sqrt(exact[s] * (1 - exact[s]) / draws * 20.0) + 1e-3,
                        $"state {s}: {p:F4} against π_w's {exact[s]:F4}");
        }
    }

    // ±1 on a ring of states, a symmetric proposal.
    private sealed class Ring(int states) : IRjMove<int>
    {
        public string Key => "ring";
        public string ReverseKey => "ring";
        public double Weight(int state) => 1.0;
        public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng)
            => new Proposal<int>((current + (rng.NextDouble() < 0.5 ? 1 : states - 1)) % states, 0.0, 0.0);
    }
}
