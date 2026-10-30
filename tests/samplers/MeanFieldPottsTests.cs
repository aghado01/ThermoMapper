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
/// The mean-field Potts toy (N17 Phase 1's addition): the lumped occupation chain, its exact references and the engine on
/// it. The Hastings ratio on n reduces to the likelihood ratio; the partition sum of the density of states equals the
/// brute-force composition sum; the equal-weight b_t approaches the derived β_cJ = 2(q − 1)·ln(q − 1)/(q − 2) as N grows;
/// and at small N the learned weights reach SAMC's limit on the exact density of states, the reweighted law at b_t and
/// ln Z(b) match the exact ones.
/// </summary>
public sealed class MeanFieldPottsTests
{
    private readonly ITestOutputHelper _out;
    public MeanFieldPottsTests(ITestOutputHelper output) => _out = output;

    private const double Bound = 4.5;   // in the eight production chains' spread errors

    private static double Coupling(int q) => MeanFieldPotts.CriticalCoupling(q) / 0.6;

    [Fact]
    public void LumpedMove_HastingsRatioIsTheLikelihoodRatio()
    {
        // log π(n′) − log π(n) + log q(n′→n)/q(n→n′) = 0: the prior ratio n_k/(n_l + 1) cancels the proposal ratio.
        var model = new MeanFieldPotts(40, 4, 2.0);
        var move = new OccupationMove(model);
        var rng = new Xoshiro256PlusPlus(5);
        int[] n = [25, 9, 6, 0];
        for (int i = 0; i < 2_000; i++)
        {
            Proposal<int[]> p = move.Propose(n, rng)!.Value;
            Assert.Equal(40, p.Candidate.Sum());
            Assert.Equal(2, n.Zip(p.Candidate, (a, b) => Math.Abs(a - b)).Sum());
            Assert.Equal(0.0, model.LogPrior(p.Candidate) - model.LogPrior(n) + p.LogProposalRatio, 10);
            n = p.Candidate;
        }

        // The lumped chain at fixed b samples the exact law of ℓ₁.
        var small = new MeanFieldPotts(24, 3, Coupling(3));
        var exact = new MeanFieldExact(small);
        const double b = 0.55;
        var chain = new ReversibleJumpChain<int[]>([new OccupationMove(small)], small.Target, small.Disordered(), new Xoshiro256PlusPlus(7),
                                                   new Temperature(b, TemperingPath.Likelihood));
        double sum = 0.0, sq = 0.0;
        const int draws = 400_000;
        for (int i = 0; i < 20_000; i++) chain.Step();
        for (int i = 0; i < draws; i++)
        {
            double l = small.Level(chain.Step());
            sum += l;
            sq += l * l;
        }
        double mean = sum / draws, sd = Math.Sqrt(sq / draws - mean * mean);
        double exactMean = exact.Expectation(b, l => l);
        _out.WriteLine($"[meanfield] lumped chain, N = 24, q = 3, b = {b}: E[ℓ₁] {mean:F3} against {exactMean:F3} (sd {sd:F2})");
        Assert.True(Math.Abs(mean - exactMean) < 0.2 * sd, $"the lumped chain's mean ℓ₁ {mean:F3} departs from the exact {exactMean:F3}");
    }

    [Fact]
    public void DensityOfStates_ByPartitions_IsTheCompositionSum()
    {
        foreach ((int q, int n) in new[] { (3, 20), (4, 12), (5, 8), (3, 7) })
        {
            var model = new MeanFieldPotts(n, q, 1.0);
            (double[] levels, double[] logMasses) = model.DensityOfStates();
            var brute = new Dictionary<double, double>();
            var c = new int[q];
            void Visit(int k, int left)
            {
                if (k == q - 1)
                {
                    c[k] = left;
                    double l = model.Level(c);
                    brute[l] = (brute.TryGetValue(l, out double v) ? v : 0.0) + Math.Exp(model.LogPrior(c));
                    return;
                }
                for (int a = 0; a <= left; a++)
                {
                    c[k] = a;
                    Visit(k + 1, left - a);
                }
            }
            Visit(0, n);
            Assert.Equal(brute.Count, levels.Length);
            for (int j = 0; j < levels.Length; j++) Assert.Equal(brute[levels[j]], Math.Exp(logMasses[j]), 1e-13);
            Assert.Equal(1.0, logMasses.Sum(Math.Exp), 12);
        }
    }

    [Fact]
    public void TheEqualWeightPoint_ApproachesTheDerivedTransition()
    {
        foreach ((int q, int[] sizes, double tolerance) in new[] { (3, new[] { 128, 512, 2048 }, 2e-3), (4, new[] { 64, 128, 256 }, 1e-2) })
        {
            double closed = MeanFieldPotts.CriticalCoupling(q);
            var gaps = new List<double>();
            foreach (int n in sizes)
            {
                var exact = new MeanFieldExact(new MeanFieldPotts(n, q, Coupling(q)));
                double split = exact.Split();
                double bt = exact.EqualWeight(split);
                (int dis, int dip, int ord) = exact.Phases(bt, split);
                double[] smooth = exact.SmoothedLogLaw(bt);
                gaps.Add(Math.Abs(bt * Coupling(q) - closed));
                _out.WriteLine($"[meanfield] q = {q}, N = {n}: b_t = {bt:F4}, b_tJ = {bt * Coupling(q):F4} against {closed:F4}; " +
                               $"ln(peak/dip) {Math.Min(smooth[dis], smooth[ord]) - smooth[dip]:F2}");
            }
            for (int i = 1; i < gaps.Count; i++) Assert.True(gaps[i] < gaps[i - 1], $"q = {q}: b_tJ should approach the limit as N grows");
            Assert.True(gaps[^1] < tolerance * closed, $"q = {q}: b_tJ is {gaps[^1]:E2} from the limit at N = {sizes[^1]}");
        }
        Assert.Equal(4.0 * Math.Log(2.0), MeanFieldPotts.CriticalCoupling(3), 14);
        Assert.Equal(3.0 * Math.Log(3.0), MeanFieldPotts.CriticalCoupling(4), 14);
    }

    // A flat-start run at N with the window from the exact law's 1% quantile at b_lo = 0.3 to its 99% at b = 1. The
    // nodes are a b = 1 sd apart but at least twice the widest gap between levels in the window: near the ordered top a
    // small lattice's levels thin out (at N = 32, q = 4, the top two sit 5.3 nats apart), and a hat that sees mass only at
    // its edges makes θ* extreme there.
    private static (MulticanonicalResult<int[]> Run, MeanFieldExact Exact) Learn(int q, int n, long learning, int production, int seed)
    {
        var model = new MeanFieldPotts(n, q, Coupling(q));
        var exact = new MeanFieldExact(model);
        double low = exact.Quantile(0.3, 0.01), high = exact.Quantile(1.0, 0.99);
        double gap = 0.0;
        for (int j = 1; j < exact.Levels.Length; j++)
            if (exact.Levels[j - 1] >= low && exact.Levels[j] <= high) gap = Math.Max(gap, exact.Levels[j] - exact.Levels[j - 1]);
        int nodes = (int)Math.Ceiling((high - low) / Math.Max(exact.StandardDeviation(1.0), 2.0 * gap)) + 1;
        var spec = new FlatHistogramSpec
        {
            Chains = 8, LearningIterations = learning, ProductionDraws = production,
            Window = new FlatHistogramWindow(low, high), LowSlope = 0.3, Nodes = nodes,
        };
        return (MulticanonicalRun.Run(spec, model.Target, new MeanFieldModel(model), seed, parallel: Serial), exact);
    }

    [Theory]
    [InlineData(3, 32, 11)]
    [InlineData(3, 64, 13)]
    [InlineData(4, 32, 17)]
    public void Engine_OnTheMeanFieldToy_MatchesTheExactReferences(int q, int n, int seed)
    {
        // The learning's length grows as N², as the crossing does.
        (MulticanonicalResult<int[]> run, MeanFieldExact exact) = Learn(q, n, learning: 300L * n * n, production: 60_000, seed);

        // The learned weights against SAMC's limit on the exact density of states.
        MulticanonicalWeights limit = SamcLearning.FixedPoint(run.Weights, exact.Levels, exact.LogMasses, out bool[] empty);
        double rms = CentredRms(run.Weights.LogWeights, limit.LogWeights, Array.ConvertAll(empty, e => !e));

        // The reweighted law at b_t, in six bins of equal exact mass, and its ordered phase's mass and mean ℓ₁.
        double split = exact.Split(), bt = exact.EqualWeight(split);
        double[] edges = Enumerable.Range(1, 5).Select(i => exact.Quantile(bt, i / 6.0)).ToArray();
        int Bin(double l) => edges.Count(e => l > e);
        double[] law = exact.Law(bt);
        var exactBins = new double[6];
        for (int j = 0; j < law.Length; j++) exactBins[Bin(exact.Levels[j])] += law[j];
        double worst = 0.0;
        var checks = new List<(string Name, ReweightedEstimate Estimate, double Exact)>
        {
            ("P(ordered)", Reweighting.Expectation(run.Draws!, bt, (l, _) => l > split ? 1.0 : 0.0), exact.OrderedMass(bt, split)),
            ("E[ℓ₁]", Reweighting.Expectation(run.Draws!, bt, (l, _) => l), exact.Expectation(bt, l => l)),
        };
        for (int i = 0; i < 6; i++)
        {
            int bin = i;
            checks.Add(($"bin {i}", Reweighting.Expectation(run.Draws!, bt, (l, _) => Bin(l) == bin ? 1.0 : 0.0), exactBins[i]));
        }
        foreach ((string name, ReweightedEstimate e, double value) in checks)
        {
            double z = (e.Mean - value) / e.MonteCarloError;
            worst = Math.Max(worst, Math.Abs(z));
            Assert.True(Math.Abs(z) < Bound, $"q = {q}, N = {n}: {name} at b_t is {e.Mean:F4} ± {e.MonteCarloError:F4}, exactly {value:F4}");
        }

        // ln Z(b) on [b_lo, 1].
        double[] grid = [0.3, 0.45, bt, 0.75, 0.9, 1.0];
        LogEvidenceCurve curve = Reweighting.LogEvidence(run.Draws!, grid);
        double worstZ = 0.0;
        for (int g = 0; g < grid.Length; g++)
        {
            double z = (curve.LogEvidence[g] - exact.LogEvidence(grid[g])) / curve.MonteCarloError[g];
            worstZ = Math.Max(worstZ, Math.Abs(z));
            Assert.True(Math.Abs(z) < Bound, $"q = {q}, N = {n}: ln Z({grid[g]:F3}) is {curve.LogEvidence[g]:F4} ± {curve.MonteCarloError[g]:F4}, exactly {exact.LogEvidence(grid[g]):F4}");
        }
        _out.WriteLine($"[meanfield] q = {q}, N = {n}: K = {run.Weights.Nodes}, learned ln w against θ* {rms:F3} RMS; at b_t = {bt:F4}, worst |z| {worst:F2} " +
                       $"over the ordered mass, E[ℓ₁] and six bins; ln Z(b) worst |z| {worstZ:F2}; round trips {run.ProductionDiagnostics!.TotalRoundTrips}");
        Assert.True(rms < 0.15, $"q = {q}, N = {n}: the learned ln w departs from θ* by {rms:F3}");
    }
}
