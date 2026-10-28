using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Samplers.Tests;

/// <summary>
/// Per-replica blocks in parallel tempering. The state is (k, θ): a discrete index k moved only by reversible-jump moves,
/// and a continuous θ moved only by a slice block at the replica's own temperature. Each rung's law, its index law and
/// θ's conditional means, matches its tempered target by enumeration over k and quadrature over θ; each pair's swap
/// rejection matches its expectation under the product of the rung laws, which a swap out of balance would move; a
/// resumed ladder continues bit for bit; a block that draws nothing leaves the ladder's trajectory exactly; and the same
/// block serves a lone chain through a <see cref="BlockScan{TState}"/>.
/// </summary>
public sealed class ParallelTemperingBlockTests
{
    private readonly ITestOutputHelper _out;
    public ParallelTemperingBlockTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(TemperingPath.Likelihood)]
    [InlineData(TemperingPath.Powered)]
    public void EveryRung_SamplesItsTemperedLaw_TheIndexByMovesAndTheCoordinateByTheBlock(TemperingPath path)
    {
        double[] betas = path == TemperingPath.Likelihood ? [1.0, 0.5, 0.2, 0.0] : [1.0, 0.5, 0.2, 0.05];
        var target = new IndexAndCoordinate();
        var pt = new ParallelTempering<State>(_ => Moves, target, new State(0, 0.0), betas, path, masterSeed: 17,
                                              blocksForReplica: _ => [new CoordinateSlice(target)]);
        Assert.True(pt.HasBlocks);
        for (int i = 0; i < 5_000; i++) pt.Step();

        const int n = 200_000;
        var counts = new double[betas.Length, IndexAndCoordinate.Indices];
        var sums = new double[betas.Length, IndexAndCoordinate.Indices];
        for (int i = 0; i < n; i++)
        {
            pt.Step();
            for (int r = 0; r < betas.Length; r++)
            {
                State s = pt.Replicas[r].Current;
                counts[r, s.K] += 1.0;
                sums[r, s.K] += s.Theta;
            }
        }

        for (int r = 0; r < betas.Length; r++)
        {
            (double[] law, double[] conditionalMean) = IndexAndCoordinate.Law(betas[r], path);
            for (int k = 0; k < IndexAndCoordinate.Indices; k++)
            {
                double observed = counts[r, k] / n;
                Assert.True(Math.Abs(observed - law[k]) < 0.01, $"{path} b = {betas[r]}, P(k = {k}): {observed:F4} vs {law[k]:F4}");
                double mean = sums[r, k] / counts[r, k];
                Assert.True(Math.Abs(mean - conditionalMean[k]) < 0.03, $"{path} b = {betas[r]}, E[θ | k = {k}]: {mean:F4} vs {conditionalMean[k]:F4}");
            }
            IReadOnlyList<MoveStat> stats = pt.BlockStats(r);
            Assert.Equal(5_000 + n, stats[0].Attempts);
        }
        _out.WriteLine($"[pt-block] {path}: swap acceptance {(double)pt.SwapAccepts / pt.SwapAttempts:F3}");
    }

    [Fact]
    public void Swaps_KeepTheProductOfTheRungLaws_TheirRejectionMatchesItsExpectation()
    {
        // In stationarity the replicas are independent, each at its rung's law, and every scan records each pair's
        // rejection 1 − min(1, e^{log α}) on them before any swap. Its long-run mean is then the expectation under the
        // product law, by quadrature over both rungs' (k, θ); stale evaluations or an unbalanced swap would move it.
        double[] betas = [1.0, 0.6, 0.3, 0.1];
        const TemperingPath path = TemperingPath.Likelihood;
        var target = new IndexAndCoordinate();
        var burnIn = new ParallelTempering<State>(_ => Moves, target, new State(1, 0.0), betas, path, masterSeed: 29,
                                                  blocksForReplica: _ => [new CoordinateSlice(target)]);
        for (int i = 0; i < 5_000; i++) burnIn.Step();
        // PairStats averages over every scan, so the measured run restarts from the burnt-in ladder with them cleared.
        var restarted = ParallelTempering<State>.Restore(_ => Moves, target, burnIn.Snapshot() with
        {
            Scans = 0, RejectionSum = new double[betas.Length - 1], RoundRejectionSum = new double[betas.Length - 1], RoundScans = 0,
        }, blocksForReplica: _ => [new CoordinateSlice(target)]);
        for (int i = 0; i < 200_000; i++) restarted.Step();

        IReadOnlyList<SwapPairStat> stats = restarted.PairStats;
        for (int p = 0; p < betas.Length - 1; p++)
        {
            double expected = IndexAndCoordinate.ExpectedRejection(betas[p], betas[p + 1], path);
            _out.WriteLine($"[pt-block] pair {betas[p]}–{betas[p + 1]}: rejection {stats[p].MeanRejection:F4} vs {expected:F4}");
            Assert.True(Math.Abs(stats[p].MeanRejection - expected) < 0.01, $"pair {p}: {stats[p].MeanRejection:F4} vs {expected:F4}");
        }
    }

    [Fact]
    public void Resume_IsBitIdentical_AndTheBlocksMustMatchTheSnapshot()
    {
        double[] betas = [1.0, 0.4, 0.1];
        var target = new IndexAndCoordinate();
        IReadOnlyList<IBlock<ReversibleJumpChain<State>>> Blocks(int replica) => [new CoordinateSlice(target)];
        var straight = new ParallelTempering<State>(_ => Moves, target, new State(0, 0.0), betas, TemperingPath.Likelihood, 5,
                                                    SwapSchedule.StochasticEvenOdd, blocksForReplica: Blocks);
        var first = new ParallelTempering<State>(_ => Moves, target, new State(0, 0.0), betas, TemperingPath.Likelihood, 5,
                                                 SwapSchedule.StochasticEvenOdd, blocksForReplica: Blocks);
        for (int i = 0; i < 400; i++)
        {
            straight.Step();
            first.Step();
        }
        LadderSnapshot<State> snapshot = first.Snapshot();
        Assert.NotNull(snapshot.BlockRngs);
        Assert.Throws<ArgumentException>(() => ParallelTempering<State>.Restore(_ => Moves, target, snapshot));
        ParallelTempering<State> resumed = ParallelTempering<State>.Restore(_ => Moves, target, snapshot, blocksForReplica: Blocks);
        for (int i = 0; i < 400; i++)
        {
            straight.Step();
            resumed.Step();
            for (int r = 0; r < betas.Length; r++)
            {
                Assert.Equal(straight.Replicas[r].Current.K, resumed.Replicas[r].Current.K);
                Assert.Equal(BitConverter.DoubleToInt64Bits(straight.Replicas[r].Current.Theta),
                             BitConverter.DoubleToInt64Bits(resumed.Replicas[r].Current.Theta));
            }
        }
        for (int r = 0; r < betas.Length; r++) Assert.Equal(straight.BlockStats(r), resumed.BlockStats(r));
        Assert.Equal(straight.SwapAccepts, resumed.SwapAccepts);

        var plain = new ParallelTempering<State>(_ => Moves, target, new State(0, 0.0), betas, TemperingPath.Likelihood, 5);
        Assert.Throws<ArgumentException>(() => ParallelTempering<State>.Restore(_ => Moves, target, plain.Snapshot(), blocksForReplica: Blocks));
        Assert.Empty(plain.BlockStats(0));
    }

    [Fact]
    public void ABlockThatDrawsNothing_LeavesTheLaddersTrajectoryBitForBit_AndKeepsItsCadence()
    {
        // The block streams come after the ladder's others, so the replica, swap and reference streams are those of a
        // ladder without blocks. The block runs on every second scan, and not at all on the exact reference rung.
        double[] betas = [1.0, 0.5, 0.2, 0.0];
        var target = new IndexAndCoordinate();
        var reference = new PriorDraw(target);
        var plain = new ParallelTempering<State>(_ => Moves, target, new State(0, 0.5), betas, TemperingPath.Likelihood, 3,
                                                 reference: reference);
        var idle = new ParallelTempering<State>(_ => Moves, target, new State(0, 0.5), betas, TemperingPath.Likelihood, 3,
                                                reference: reference, blocksForReplica: _ => [new Idle(every: 2)]);
        const int scans = 1_000;
        for (int i = 0; i < scans; i++)
        {
            plain.Step();
            idle.Step();
            for (int r = 0; r < betas.Length; r++) Assert.Equal(plain.Replicas[r].Current, idle.Replicas[r].Current);
        }
        Assert.Equal(plain.SwapAccepts, idle.SwapAccepts);
        for (int r = 0; r < betas.Length - 1; r++) Assert.Equal(scans / 2, idle.BlockStats(r)[0].Attempts);
        Assert.Equal(0, idle.BlockStats(betas.Length - 1)[0].Attempts);
    }

    [Fact]
    public void TheSameBlock_ServesALoneChain_ThroughABlockScanBesideItsMoves()
    {
        var target = new IndexAndCoordinate();
        var chain = new ReversibleJumpChain<State>(Moves, target, new State(2, 0.0), new Xoshiro256PlusPlus(41));
        var scan = new BlockScan<ReversibleJumpChain<State>>([new MoveStep<State>(), new CoordinateSlice(target)], chain,
                                                             new Xoshiro256PlusPlus(43), c => c);
        for (int i = 0; i < 5_000; i++) scan.Step();
        const int n = 200_000;
        var counts = new double[IndexAndCoordinate.Indices];
        var sums = new double[IndexAndCoordinate.Indices];
        for (int i = 0; i < n; i++)
        {
            State s = scan.Step().Current;
            counts[s.K] += 1.0;
            sums[s.K] += s.Theta;
        }
        (double[] law, double[] conditionalMean) = IndexAndCoordinate.Law(1.0, TemperingPath.Likelihood);
        for (int k = 0; k < IndexAndCoordinate.Indices; k++)
        {
            Assert.True(Math.Abs(counts[k] / n - law[k]) < 0.01, $"P(k = {k}): {counts[k] / n:F4} vs {law[k]:F4}");
            Assert.True(Math.Abs(sums[k] / counts[k] - conditionalMean[k]) < 0.03, $"E[θ | k = {k}]");
        }
        Assert.Equal(["moves", "theta"], scan.BlockStats().Select(s => s.Key));
    }

    [Fact]
    public void BlockSets_AreValidated()
    {
        var target = new IndexAndCoordinate();
        double[] betas = [1.0, 0.5];
        Assert.Throws<ArgumentException>(() => new ParallelTempering<State>(_ => Moves, target, new State(0, 0.0), betas,
            TemperingPath.Likelihood, 1, blocksForReplica: _ => null!));
        Assert.Throws<ArgumentException>(() => new ParallelTempering<State>(_ => Moves, target, new State(0, 0.0), betas,
            TemperingPath.Likelihood, 1, blocksForReplica: _ => [new Idle(every: 0)]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ParallelTempering<State>(_ => Moves, target, new State(0, 0.0), betas,
            TemperingPath.Likelihood, 1, blocksForReplica: _ => []).BlockStats(2));
    }

    private static readonly IRjMove<State>[] Moves = [new IndexStep("up", "down", +1), new IndexStep("down", "up", -1)];

    private readonly record struct State(int K, double Theta);

    // k ∈ {0, 1, 2, 3} with prior p_k, θ ~ N(0, 2²); ℓ_b = 3·ln(1 + b·a_k) − b·c_k − 2b(θ − μ_k)², not linear in b. Off
    // the index range the state is off the support.
    private sealed class IndexAndCoordinate : ISamplingTarget<State>
    {
        public const int Indices = 4;
        private static readonly double[] Prior = [0.4, 0.2, 0.3, 0.1];
        private static readonly double[] A = [0.0, 6.0, 2.0, 15.0];
        private static readonly double[] C = [0.0, 3.0, 0.5, 6.0];
        private static readonly double[] Mu = [-1.0, 0.0, 0.5, 1.5];

        private static double LogPrior(State s) => Math.Log(Prior[s.K]) - s.Theta * s.Theta / 8.0 - Math.Log(2.0 * Math.Sqrt(2.0 * Math.PI));
        private static double LogLikelihood(State s, double b)
            => 3.0 * Math.Log(1.0 + b * A[s.K]) - b * C[s.K] - 2.0 * b * (s.Theta - Mu[s.K]) * (s.Theta - Mu[s.K]);

        public ITargetEvaluation Evaluate(State state) => new Evaluation(state);

        private static double LogDensity(State s, double b, TemperingPath path)
            => path == TemperingPath.Likelihood ? LogPrior(s) + LogLikelihood(s, b) : b * (LogPrior(s) + LogLikelihood(s, 1.0));

        // θ by the trapezoid on [−12, 12], where every rung's conditional has all but a negligible tail.
        private const int Points = 4_801;
        private static double ThetaAt(int i) => -12.0 + 24.0 * i / (Points - 1);

        /// <summary>P_b(k) and E_b[θ | k].</summary>
        public static (double[] Law, double[] ConditionalMean) Law(double b, TemperingPath path)
        {
            var mass = new double[Indices];
            var first = new double[Indices];
            for (int k = 0; k < Indices; k++)
                for (int i = 0; i < Points; i++)
                {
                    double theta = ThetaAt(i), w = (i == 0 || i == Points - 1 ? 0.5 : 1.0) * Math.Exp(LogDensity(new State(k, theta), b, path));
                    mass[k] += w;
                    first[k] += w * theta;
                }
            double total = mass.Sum();
            return (mass.Select(m => m / total).ToArray(), first.Zip(mass, (f, m) => f / m).ToArray());
        }

        /// <summary>E[1 − min(1, e^{log α})] under π_a ⊗ π_c for a swap between rungs a and c.</summary>
        public static double ExpectedRejection(double a, double c, TemperingPath path)
        {
            const int grid = 481;
            double h = 24.0 / (grid - 1);
            var da = new double[Indices, grid];
            var dc = new double[Indices, grid];
            for (int k = 0; k < Indices; k++)
                for (int i = 0; i < grid; i++)
                {
                    var s = new State(k, -12.0 + i * h);
                    da[k, i] = LogDensity(s, a, path);
                    dc[k, i] = LogDensity(s, c, path);
                }
            double za = 0.0, zc = 0.0;
            for (int k = 0; k < Indices; k++) for (int i = 0; i < grid; i++) { za += Math.Exp(da[k, i]); zc += Math.Exp(dc[k, i]); }
            double expected = 0.0;
            for (int k = 0; k < Indices; k++)
                for (int i = 0; i < grid; i++)
                {
                    double px = Math.Exp(da[k, i]) / za;   // rung a holds x = (k, θ_i)
                    if (px < 1e-300) continue;
                    for (int l = 0; l < Indices; l++)
                        for (int j = 0; j < grid; j++)
                        {
                            double py = Math.Exp(dc[l, j]) / zc;   // rung c holds y = (l, θ_j)
                            double logA = da[l, j] + dc[k, i] - da[k, i] - dc[l, j];
                            expected += px * py * (logA >= 0.0 ? 0.0 : 1.0 - Math.Exp(logA));
                        }
                }
            return expected;
        }

        private sealed class Evaluation(State s) : ITargetEvaluation
        {
            public double LogPrior => s.K is >= 0 and < Indices ? IndexAndCoordinate.LogPrior(s) : double.NegativeInfinity;
            public double LogLikelihood(double inverseTemperature)
                => s.K is >= 0 and < Indices ? IndexAndCoordinate.LogLikelihood(s, inverseTemperature) : double.NegativeInfinity;
        }
    }

    // θ by one slice step at the replica's own temperature, the index held: π_b-invariant, handed back with its evaluation.
    private sealed class CoordinateSlice(ISamplingTarget<State> target) : IBlock<ReversibleJumpChain<State>>
    {
        public string Name => "theta";

        public bool Update(ReversibleJumpChain<State> chain, Xoshiro256PlusPlus rng)
        {
            State s = chain.Current;
            Temperature t = chain.Temperature;
            double theta = s.Theta, logTarget = chain.CurrentLogDensity;
            SliceSampler.Step(ref theta, ref logTarget, 1.0, 20, v => t.LogDensity(target.Evaluate(s with { Theta = v })), rng);
            State next = s with { Theta = theta };
            chain.SetState(next, target.Evaluate(next));
            return true;
        }
    }

    // Draws nothing and changes nothing.
    private sealed class Idle(int every) : IBlock<ReversibleJumpChain<State>>
    {
        public string Name => "idle";
        public int Every => every;
        public bool Update(ReversibleJumpChain<State> chain, Xoshiro256PlusPlus rng) => false;
    }

    // The likelihood path's b = 0 member exactly: k from its prior, θ from N(0, 2²).
    private sealed class PriorDraw(IndexAndCoordinate target) : IReferenceSampler<State>
    {
        private static readonly double[] Prior = [0.4, 0.2, 0.3, 0.1];
        public TemperingPath Path => TemperingPath.Likelihood;
        public long Attempts { get; private set; }

        public (State State, ITargetEvaluation Evaluation) Draw(Xoshiro256PlusPlus rng)
        {
            Attempts++;
            double u = rng.NextDouble(), cumulative = 0.0;
            int k = 0;
            for (; k < Prior.Length - 1; k++)
            {
                cumulative += Prior[k];
                if (u < cumulative) break;
            }
            var s = new State(k, 2.0 * Distributions.Normal.SampleStandard(rng));
            return (s, target.Evaluate(s));
        }
    }

    // k ± 1 with θ held; off the index range the candidate is off the support.
    private sealed class IndexStep(string key, string reverseKey, int delta) : IRjMove<State>
    {
        public string Key => key;
        public string ReverseKey => reverseKey;
        public double Weight(State state) => 1.0;
        public Proposal<State>? Propose(State current, Xoshiro256PlusPlus rng) => new(current with { K = current.K + delta }, 0.0, 0.0);
    }
}
