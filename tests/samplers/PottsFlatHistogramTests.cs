using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;
using static Maths.Samplers.Tests.FlatHistogramToys;

namespace Maths.Samplers.Tests;

/// <summary>
/// The flat-histogram engine on toy 2 (N17 Phase 1, step 3): a 6-state Potts model on a periodic 3 × 3 lattice, qᴺ ≈
/// 1.0·10⁷ configurations under a uniform prior, ℓ₁ = J·(agreeing bonds) with J = 1.5, and single-site Metropolis
/// moves. Its exact density of states, by enumeration, gives an energy law with a disordered peak and an ordered one,
/// equally high at the transition, across a dip and the empty levels 13 and 15–17. The learned weights reach SAMC's limit on
/// the exact density of states, ln Z(b) matches the exact curve at every b on the grid, and the reweighted energy law
/// at the transition matches the exact one.
/// </summary>
public sealed class PottsFlatHistogramTests
{
    private readonly ITestOutputHelper _out;
    public PottsFlatHistogramTests(ITestOutputHelper output) => _out = output;

    private const double Coupling = 1.5;
    private const double Bound = 4.5;   // in the eight production chains' spread errors
    internal static readonly Potts Lattice = new(3, 3, 6, Coupling);
    private static readonly Lazy<long[]> ExactCounts = new(() => Lattice.Counts());

    // The window spans every level, a node on each, with a flat extension below.
    internal static FlatHistogramSpec PottsSpec(long learning, int production, int chains = 8) => new()
    {
        Chains = chains,
        LearningIterations = learning,
        ProductionDraws = production,
        Window = new FlatHistogramWindow(0.0, Coupling * Lattice.Bonds.Length),
        LowSlope = 0.0,
        Nodes = Lattice.Bonds.Length + 1,
    };

    // P_b(n) by enumeration.
    private static double[] ExactLaw(double b)
    {
        long[] counts = ExactCounts.Value;
        var log = new double[counts.Length];
        for (int n = 0; n < counts.Length; n++) log[n] = counts[n] > 0 ? Math.Log(counts[n]) + b * Coupling * n : double.NegativeInfinity;
        double z = LogSumExp(log);
        return Array.ConvertAll(log, l => Math.Exp(l - z));
    }

    private static double ExactLogEvidence(double b)
    {
        long[] counts = ExactCounts.Value;
        double logTotal = Lattice.Sites * Math.Log(Lattice.Q);
        var log = new double[counts.Length];
        for (int n = 0; n < counts.Length; n++) log[n] = counts[n] > 0 ? Math.Log(counts[n]) - logTotal + b * Coupling * n : double.NegativeInfinity;
        return LogSumExp(log);
    }

    // The transition's b, where the disordered peak (the most likely n ≤ 12) and the ordered one (n ≥ 14) are equally
    // high, the finite-size criterion for a first-order transition, by bisection.
    internal static double TransitionB()
    {
        double lo = 0.0, hi = 1.0;
        for (int i = 0; i < 60; i++)
        {
            double mid = 0.5 * (lo + hi);
            double[] law = ExactLaw(mid);
            if (law.Take(13).Max() > law.Skip(14).Max()) lo = mid;
            else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    [Fact]
    public void Enumeration_CountsEveryConfiguration_AndTheTransitionIsBimodal()
    {
        long[] counts = ExactCounts.Value;
        Assert.Equal((long)Math.Pow(6, 9), counts.Sum());
        Assert.Equal(6, counts[18]);              // the six uniform configurations
        Assert.Equal(9 * 6 * 5, counts[14]);      // one site of another colour breaks its four bonds
        foreach (int n in new[] { 13, 15, 16, 17 }) Assert.Equal(0, counts[n]);

        double bc = TransitionB();
        double[] law = ExactLaw(bc);
        int disordered = Enumerable.Range(0, 13).OrderByDescending(n => law[n]).First();
        int dipAt = Enumerable.Range(disordered + 1, 13 - disordered - 1).OrderBy(n => law[n]).First();
        double ordered = law.Skip(14).Max();
        double orderedMass = law.Skip(14).Sum();
        _out.WriteLine($"[potts] 3×3, q = 6, J = {Coupling}: b_c = {bc:F4} (βJ = {bc * Coupling:F4}); P(n) there: " +
                       string.Join(" ", Enumerable.Range(0, law.Length).Select(n => $"{n}:{law[n]:F4}")));
        _out.WriteLine($"[potts] disordered peak P({disordered}) = {law[disordered]:F4}, ordered peak {ordered:F4}, dip P({dipAt}) = {law[dipAt]:F4}: " +
                       $"peak/dip {Math.Min(law[disordered], ordered) / law[dipAt]:F2}, no mass at n = 13, 15–17; ordered mass {orderedMass:F3}");
        Assert.True(Math.Min(law[disordered], ordered) > 2.0 * law[dipAt], "the energy law at the transition should be bimodal");
    }

    [Fact]
    public void Potts_LearnsSamcsLimit_OnTheExactDensityOfStates()
    {
        MulticanonicalResult<int[]> run = MulticanonicalRun.Run(PottsSpec(200_000, 1, chains: 4), Lattice.Target, new PottsModel(Lattice), masterSeed: 7);
        MulticanonicalWeights exact = PottsFixedPoint(run.Weights, out bool[] empty);
        bool[] full = Array.ConvertAll(empty, e => !e);
        double error = CentredRms(run.Weights.LogWeights, exact.LogWeights, full);

        // On levels at the nodes, θ* is Liang et al.'s ln g − ln(π + d) exactly.
        long[] counts = ExactCounts.Value;
        double[] pi = run.Weights.FlatOccupancy();
        double d = Enumerable.Range(0, pi.Length).Where(k => empty[k]).Sum(k => pi[k]) / full.Count(f => f);
        double[] liang = Enumerable.Range(0, pi.Length).Select(k => empty[k] ? 0.0 : -(Math.Log(counts[k]) - Math.Log(pi[k] + d))).ToArray();
        _out.WriteLine($"[potts] learned ln w against θ*: centred RMS {error:F3} over {full.Count(f => f)} levels; θ* against ln g − ln(π + d): " +
                       $"{CentredRms(exact.LogWeights, liang, full):E1}; ln g spans {Math.Log(counts.Max()) - Math.Log(6):F1} nats");
        Assert.Equal([13, 15, 16, 17], Enumerable.Range(0, empty.Length).Where(k => empty[k]));
        Assert.True(CentredRms(exact.LogWeights, liang, full) < 1e-8);
        Assert.True(error < 0.1, $"the learned ln w departs from θ* by {error:F3}");
    }

    [Fact]
    public void Potts_LogEvidence_AndTheTransitionsEnergyLaw_MatchTheEnumeration()
    {
        // The ratio estimator needs the prior normalized: q^{−N} on each of the qᴺ configurations.
        Assert.Equal(-Lattice.Sites * Math.Log(Lattice.Q), Lattice.Target.Evaluate(new int[Lattice.Sites]).LogPrior);
        Assert.Equal((long)Math.Pow(Lattice.Q, Lattice.Sites), ExactCounts.Value.Sum());

        MulticanonicalResult<int[]> run = MulticanonicalRun.Run(PottsSpec(200_000, 100_000), Lattice.Target, new PottsModel(Lattice), masterSeed: 11);
        MulticanonicalDraws draws = run.Draws!;

        double[] grid = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0];
        LogEvidenceCurve curve = Reweighting.LogEvidence(draws, grid);
        for (int g = 0; g < grid.Length; g++)
        {
            double exact = ExactLogEvidence(grid[g]);
            double z = g == 0 ? 0.0 : (curve.LogEvidence[g] - exact) / curve.MonteCarloError[g];
            _out.WriteLine($"[potts] b = {grid[g]}: ln Z {curve.LogEvidence[g]:F4} ± {curve.MonteCarloError[g]:F4} against {exact:F4} (z = {z:F2})");
            Assert.True(Math.Abs(z) < Bound, $"b = {grid[g]}: ln Z is {curve.LogEvidence[g]:F4} ± {curve.MonteCarloError[g]:F4}, exactly {exact:F4}");
        }

        double bc = TransitionB();
        double[] law = ExactLaw(bc);
        double worst = 0.0;
        for (int n = 0; n < law.Length; n++)
        {
            double level = Coupling * n;
            ReweightedEstimate p = Reweighting.Expectation(draws, bc, (l, _) => l == level ? 1.0 : 0.0);
            if (law[n] == 0.0)
            {
                Assert.Equal(0.0, p.Mean);
                continue;
            }
            double z = (p.Mean - law[n]) / p.MonteCarloError;
            worst = Math.Max(worst, Math.Abs(z));
            Assert.True(Math.Abs(z) < Bound, $"P(n = {n}) at b_c is {p.Mean:F4} ± {p.MonteCarloError:F4}, exactly {law[n]:F4}");
        }
        WeightedDraws atOne = Reweighting.Weigh(draws, 1.0);
        _out.WriteLine($"[potts] the energy law at b_c = {bc:F4}: worst |z| {worst:F2} over the occupied levels; Kish ESS at b = 1 " +
                       $"{atOne.EffectiveSampleSize:F0} of {atOne.Weights.Length}; round trips {string.Join(",", run.ProductionDiagnostics!.RoundTrips)}");
    }

    [Fact]
    public void Exploration_TunnellingAgainstParallelTempering_AtEqualWork()
    {
        // Exploration, not a claim: the flat-histogram chains' round trips between the window's ends, against a tempering
        // ladder's round trips between b = 0 and b = 1, per million single-site proposals. The ladder is adapted in its
        // burn-in, and both spend the same proposals.
        const long learning = 100_000;
        const int production = 100_000, chains = 4;
        MulticanonicalResult<int[]> run = MulticanonicalRun.Run(PottsSpec(learning, production, chains), Lattice.Target, new PottsModel(Lattice), masterSeed: 13);
        long samcWork = (learning + production) * chains;
        double samcRate = 1e6 * run.ProductionDiagnostics!.TotalRoundTrips / (production * (double)chains);

        double[] ladder = Enumerable.Range(0, 8).Select(i => 1.0 - i / 7.0).ToArray();
        var pt = new ParallelTempering<int[]>(_ => [new SingleSiteFlip(Lattice)], Lattice.Target, new int[Lattice.Sites], ladder,
                                              TemperingPath.Likelihood, masterSeed: 13);
        long scans = samcWork / ladder.Length, burn = scans / 2;
        for (long s = 0; s < burn; s++)
        {
            pt.Step();
            if ((s + 1) % (burn / 5) == 0) pt.AdaptLadder();
        }
        long before = pt.RoundTrips;
        for (long s = burn; s < scans; s++) pt.Step();
        double ptRate = 1e6 * (pt.RoundTrips - before) / ((scans - burn) * (double)ladder.Length);
        _out.WriteLine($"[potts] exploration at {samcWork:N0} proposals each: flat histogram {samcRate:F1} round trips per 10⁶ production proposals " +
                       $"({run.ProductionDiagnostics.TotalRoundTrips} over {chains} chains); tempering {ptRate:F1} per 10⁶ " +
                       $"({pt.RoundTrips - before} over {ladder.Length} rungs, adapted to {string.Join(", ", pt.InverseTemperatures.Select(b => b.ToString("F2")))}); " +
                       $"swap acceptance {(double)pt.SwapAccepts / pt.SwapAttempts:F2}");
        Assert.True(run.ProductionDiagnostics.TotalRoundTrips > 0);
    }

    // SAMC's limit on the exact density of states, over the run's window.
    internal static MulticanonicalWeights PottsFixedPoint(MulticanonicalWeights shape, out bool[] empty)
    {
        long[] counts = ExactCounts.Value;
        double logTotal = Lattice.Sites * Math.Log(Lattice.Q);
        double[] levels = Enumerable.Range(0, counts.Length).Select(n => Coupling * n).ToArray();
        double[] masses = counts.Select(c => c > 0 ? Math.Log(c) - logTotal : double.NegativeInfinity).ToArray();
        return SamcLearning.FixedPoint(shape, levels, masses, out empty);
    }

    /// <summary>The Potts chains: single-site Metropolis on a colour array; the functional is ℓ₁.</summary>
    internal sealed class PottsModel(Potts lattice) : IMulticanonicalModel<int[]>
    {
        public int FunctionalDim => 1;
        public int MonitoredDim => 1;
        public IReadOnlyList<string> FunctionalNames => ["l1"];
        public string Fingerprint => new RunFingerprint().Add(lattice.Width).Add(lattice.Height).Add(lattice.Q).Add(lattice.Coupling).ToString();

        public IMulticanonicalChain<int[]> StartChain(MulticanonicalTarget<int[]> target, int chainSeed, int readoutSeed)
            => new MulticanonicalChain<int[]>(new ReversibleJumpChain<int[]>([new SingleSiteFlip(lattice)], target, new int[lattice.Sites],
                                                                             new Xoshiro256PlusPlus(chainSeed)), target, Reduce);

        public IMulticanonicalChain<int[]> RestoreChain(MulticanonicalTarget<int[]> target, IChainState state)
            => MulticanonicalChain<int[]>.Restore([new SingleSiteFlip(lattice)], target, state, Reduce);

        private void Reduce(in int[] colours, Span<double> functionals) => functionals[0] = lattice.Coupling * lattice.Agreeing(colours);
    }

    /// <summary>One site, uniformly, to one of the other q − 1 colours, uniformly: a symmetric proposal.</summary>
    internal sealed class SingleSiteFlip(Potts lattice) : IRjMove<int[]>
    {
        public string Key => "flip";
        public string ReverseKey => "flip";
        public double Weight(int[] state) => 1.0;

        public Proposal<int[]>? Propose(int[] current, Xoshiro256PlusPlus rng)
        {
            var next = (int[])current.Clone();
            int site = rng.NextInt(lattice.Sites);
            next[site] = (current[site] + 1 + rng.NextInt(lattice.Q - 1)) % lattice.Q;
            return new Proposal<int[]>(next, 0.0, 0.0);
        }
    }
}
