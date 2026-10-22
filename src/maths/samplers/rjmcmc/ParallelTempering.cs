using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace Maths.Samplers.Rjmcmc;

/// <summary>
/// A ladder's restart state: every replica's <see cref="ChainSnapshot{TState}"/> (its rung is its temperature), the
/// swap stream, the scan count, which state each rung holds, the round-trip bookkeeping and the swap statistics,
/// including those gathered since the last adaptation.
/// </summary>
public sealed record LadderSnapshot<TState>(
    IReadOnlyList<ChainSnapshot<TState>> Replicas, TemperingPath Path, SwapSchedule Schedule,
    ulong SwapRng0, ulong SwapRng1, ulong SwapRng2, ulong SwapRng3,
    long Scans, long SwapAttempts, long SwapAccepts, long RoundTrips, long AnnealedRestarts,
    int[] LabelAt, sbyte[] LastEnd, long[] PairAttempts, long[] PairAccepts, double[] RejectionSum,
    double[] RoundRejectionSum, long RoundScans);

/// <summary>
/// Parallel tempering (replica exchange) over a ladder of reversible-jump chains at inverse temperatures
/// 1 = β₀ &gt; β₁ &gt; … &gt; β_{L-1}, all on one declared <see cref="TemperingPath"/>: the powered posterior
/// <c>(prior · likelihood)^β</c> or the likelihood-tempered path <c>prior · likelihood^β</c>. After every scan, adjacent
/// replicas attempt Metropolis state swaps on the pairs the <see cref="SwapSchedule"/> chooses, each replica's density
/// evaluated at the other's state from the cached evaluations — so the hot replicas, which cross modal barriers freely,
/// ferry the cold replica (β = 1, the true posterior) between modes. The cold chain is the sampler; the rest are
/// scaffolding. For multimodal targets (e.g. a multi-peak T_c posterior) where a single chain would lock into one mode.
/// <para>The ladder reports what Syed et al. (2022) measure it by: round trips of states from the hot end to the cold
/// end and back, and each pair's rejection rate, from which <see cref="AdaptLadder"/> equalizes the rungs.</para>
/// As an <see cref="IChain{TState}"/> the ladder presents as a single chain whose draw is the cold state and whose
/// acceptance is the cold replica's, so the ensemble shell can drive a tempered replica like any other.
/// </summary>
public sealed class ParallelTempering<TState> : IChain<TState>
{
    private readonly ReversibleJumpChain<TState>[] _chains;
    private readonly Xoshiro256PlusPlus _swapRng;
    private readonly double[] _logA;              // this scan's log swap ratio per pair
    private readonly long[] _pairAttempts, _pairAccepts;
    private readonly double[] _rejectionSum;       // Σ (1 − α) per pair, every scan
    private readonly double[] _roundRejectionSum;  // the same since the last adaptation
    private long _roundScans;
    private readonly int[] _labelAt;               // which state (by its starting rung) each rung holds
    private readonly sbyte[] _lastEnd;             // per state: +1 last at the hot end, −1 at the cold end, 0 neither yet

    /// <summary>Adjacent-pair swap attempts so far.</summary>
    public long SwapAttempts { get; private set; }

    /// <summary>Adjacent-pair swaps accepted so far; <c>SwapAccepts / SwapAttempts</c> is the exchange rate.</summary>
    public long SwapAccepts { get; private set; }

    /// <param name="movesForReplica">Builds replica i's own move set (i = 0 is the cold replica). It is called once
    /// per replica, so any mutable move or kernel state is private to that replica.</param>
    /// <param name="target">The target; replica i samples it at <c>Temperature(β_i, path)</c>.</param>
    /// <param name="start">Common initial state for all replicas.</param>
    /// <param name="betas">The ladder: finite, strictly descending, cold value exactly 1 at index 0. The hot end is
    /// above 0 on the powered path, where β = 0 would be the bare reference measure; on the likelihood path β = 0 is
    /// the prior, which a target with an improper nuisance prior refuses. The array is not retained.</param>
    /// <param name="path">The tempering path every replica takes.</param>
    /// <param name="masterSeed">Seeds the per-replica streams and the swap stream via the <see cref="SeedTree"/>.</param>
    /// <param name="schedule">Which pairs each scan attempts; deterministic even-odd by default.</param>
    public ParallelTempering(Func<int, IReadOnlyList<IRjMove<TState>>> movesForReplica, ISamplingTarget<TState> target,
                             TState start, double[] betas, TemperingPath path, int masterSeed,
                             SwapSchedule schedule = SwapSchedule.DeterministicEvenOdd)
    {
        ArgumentNullException.ThrowIfNull(movesForReplica);
        ArgumentNullException.ThrowIfNull(target);
        ValidateLadder(betas, path);
        if (schedule is not (SwapSchedule.DeterministicEvenOdd or SwapSchedule.StochasticEvenOdd or SwapSchedule.Sequential))
            throw new ArgumentOutOfRangeException(nameof(schedule));

        Path = path;
        Schedule = schedule;
        int levels = betas.Length;
        int[] seeds = SeedTree.Derive(masterSeed, levels + 1);
        _chains = new ReversibleJumpChain<TState>[levels];
        for (int i = 0; i < levels; i++)
            _chains[i] = new ReversibleJumpChain<TState>(
                movesForReplica(i), target, start, new Xoshiro256PlusPlus(seeds[i]), new Temperature(betas[i], path));
        _swapRng = new Xoshiro256PlusPlus(seeds[levels]);

        _logA = new double[levels - 1];
        _pairAttempts = new long[levels - 1];
        _pairAccepts = new long[levels - 1];
        _rejectionSum = new double[levels - 1];
        _roundRejectionSum = new double[levels - 1];
        _labelAt = new int[levels];
        _lastEnd = new sbyte[levels];
        for (int i = 0; i < levels; i++) _labelAt[i] = i;
        _lastEnd[levels - 1] = +1;   // a state's first round trip starts at its first visit to the hot end
    }

    // A ladder rebuilt from a snapshot: the replicas and swap stream are restored, the bookkeeping copied after.
    private ParallelTempering(ReversibleJumpChain<TState>[] chains, Xoshiro256PlusPlus swapRng, TemperingPath path,
                              SwapSchedule schedule)
    {
        _chains = chains;
        _swapRng = swapRng;
        Path = path;
        Schedule = schedule;
        int levels = chains.Length;
        _logA = new double[levels - 1];
        _pairAttempts = new long[levels - 1];
        _pairAccepts = new long[levels - 1];
        _rejectionSum = new double[levels - 1];
        _roundRejectionSum = new double[levels - 1];
        _labelAt = new int[levels];
        _lastEnd = new sbyte[levels];
    }

    /// <summary>This ladder's restart state; the ladder itself is unchanged.</summary>
    public LadderSnapshot<TState> Snapshot()
    {
        var replicas = new ChainSnapshot<TState>[_chains.Length];
        for (int i = 0; i < replicas.Length; i++) replicas[i] = _chains[i].Snapshot();
        (ulong s0, ulong s1, ulong s2, ulong s3) = _swapRng.SaveState();
        return new LadderSnapshot<TState>(replicas, Path, Schedule, s0, s1, s2, s3, Scans, SwapAttempts, SwapAccepts,
            RoundTrips, AnnealedRestarts, (int[])_labelAt.Clone(), (sbyte[])_lastEnd.Clone(), (long[])_pairAttempts.Clone(),
            (long[])_pairAccepts.Clone(), (double[])_rejectionSum.Clone(), (double[])_roundRejectionSum.Clone(), _roundScans);
    }

    /// <summary>
    /// A ladder that continues exactly where <paramref name="snapshot"/>'s left off, over the move sets and target it
    /// ran with. Each replica's current state is evaluated once, against <paramref name="target"/>.
    /// </summary>
    public static ParallelTempering<TState> Restore(Func<int, IReadOnlyList<IRjMove<TState>>> movesForReplica,
                                                    ISamplingTarget<TState> target, LadderSnapshot<TState> snapshot)
    {
        ArgumentNullException.ThrowIfNull(movesForReplica);
        ArgumentNullException.ThrowIfNull(snapshot);
        int levels = snapshot.Replicas.Count;
        var chains = new ReversibleJumpChain<TState>[levels];
        for (int i = 0; i < levels; i++)
            chains[i] = ReversibleJumpChain<TState>.Restore(movesForReplica(i), target, snapshot.Replicas[i]);
        var ladder = new ParallelTempering<TState>(chains,
            Xoshiro256PlusPlus.Restore(snapshot.SwapRng0, snapshot.SwapRng1, snapshot.SwapRng2, snapshot.SwapRng3),
            snapshot.Path, snapshot.Schedule)
        {
            Scans = snapshot.Scans,
            SwapAttempts = snapshot.SwapAttempts,
            SwapAccepts = snapshot.SwapAccepts,
            RoundTrips = snapshot.RoundTrips,
            AnnealedRestarts = snapshot.AnnealedRestarts,
        };
        snapshot.LabelAt.CopyTo(ladder._labelAt, 0);
        snapshot.LastEnd.CopyTo(ladder._lastEnd, 0);
        snapshot.PairAttempts.CopyTo(ladder._pairAttempts, 0);
        snapshot.PairAccepts.CopyTo(ladder._pairAccepts, 0);
        snapshot.RejectionSum.CopyTo(ladder._rejectionSum, 0);
        snapshot.RoundRejectionSum.CopyTo(ladder._roundRejectionSum, 0);
        ladder._roundScans = snapshot.RoundScans;
        return ladder;
    }

    /// <summary>The path every replica takes.</summary>
    public TemperingPath Path { get; }

    /// <summary>Which pairs each scan attempts.</summary>
    public SwapSchedule Schedule { get; }

    /// <summary>The current cold-chain (β = 1) state — the posterior sample.</summary>
    public TState Cold => _chains[0].Current;

    /// <summary>The cold replica (β = 1) itself — the chain whose draws are the posterior, exposed so a consumer
    /// can read its acceptance / per-move diagnostics. Its state advances both by its own moves and by swaps.</summary>
    public ReversibleJumpChain<TState> ColdChain => _chains[0];

    /// <summary>Every replica, cold first, for per-rung diagnostics. Each one's <see cref="ReversibleJumpChain{TState}.Temperature"/>
    /// is its rung. Stepping a replica directly takes it out of the ladder's schedule.</summary>
    public IReadOnlyList<ReversibleJumpChain<TState>> Replicas => _chains;

    /// <summary>The current ladder, cold first.</summary>
    public double[] Betas
    {
        get
        {
            var b = new double[_chains.Length];
            for (int i = 0; i < b.Length; i++) b[i] = _chains[i].Temperature.InverseTemperature;
            return b;
        }
    }

    /// <summary>Scans so far: one <see cref="Step"/> each.</summary>
    public long Scans { get; private set; }

    /// <summary>
    /// Round trips so far, summed over states: a state that leaves the hot end, reaches the cold end and returns
    /// completes one (Syed et al.'s definition, with the hot end as the reference). <c>RoundTrips / Scans</c> is the
    /// round-trip rate, the measure of how well the ladder moves.
    /// </summary>
    public long RoundTrips { get; private set; }

    /// <summary>Arrivals at the cold end of states last seen at the hot end ("annealed restarts"): each carries what
    /// the hot rungs explored down to the posterior.</summary>
    public long AnnealedRestarts { get; private set; }

    /// <summary>Per adjacent pair, cold first: attempts, accepts and the mean rejection probability over every scan.</summary>
    public IReadOnlyList<SwapPairStat> PairStats
    {
        get
        {
            var stats = new SwapPairStat[_chains.Length - 1];
            for (int i = 0; i < stats.Length; i++)
                stats[i] = new SwapPairStat(_chains[i].Temperature.InverseTemperature, _chains[i + 1].Temperature.InverseTemperature,
                                            _pairAttempts[i], _pairAccepts[i], Scans > 0 ? _rejectionSum[i] / Scans : double.NaN);
            return stats;
        }
    }

    /// <summary>
    /// Advance every replica one step (local exploration), then attempt the scan's swaps. A swap of x_i and x_j is
    /// accepted with <c>log α = d_i(x_j) + d_j(x_i) − d_i(x_i) − d_j(x_j)</c>, d_i replica i's log density, which on
    /// the powered path reduces to <c>(β_i − β_j)(L_j − L_i)</c>. Every pair's rejection probability 1 − α is recorded
    /// every scan, attempted or not. Returns the cold state.
    /// </summary>
    public TState Step()
    {
        for (int i = 0; i < _chains.Length; i++) _chains[i].Step();
        Scans++;
        _roundScans++;
        int pairs = _chains.Length - 1;

        if (Schedule == SwapSchedule.Sequential)
        {
            // Every pair, cold to hot; each sees the swaps before it.
            for (int i = 0; i < pairs; i++)
            {
                double logA = Record(i);
                Attempt(i, logA);
            }
        }
        else
        {
            // The even and the odd pairs are each disjoint, so this scan's ratios are all taken before any swap.
            for (int i = 0; i < pairs; i++) _logA[i] = Record(i);
            int parity = Schedule == SwapSchedule.DeterministicEvenOdd
                ? (int)(Scans % 2 == 0 ? 0 : 1)
                : _swapRng.NextDouble() < 0.5 ? 0 : 1;
            for (int i = parity; i < pairs; i += 2) Attempt(i, _logA[i]);
        }

        TrackEnds();
        return _chains[0].Current;
    }

    /// <summary>
    /// Re-space the interior rungs so that adjacent rejection rates are equal (<see cref="TemperingLadder.Equalize"/>),
    /// from the rates observed since the last adaptation, and return those rates' barrier estimate Λ̂. The cold and hot
    /// ends stay where they are. Retempering is a query on each replica's cached evaluation: no state is re-evaluated.
    /// Adapt during burn-in only; a ladder that keeps changing does not leave the posterior invariant.
    /// </summary>
    public double AdaptLadder()
    {
        if (_roundScans == 0) throw new InvalidOperationException("No scans since the last adaptation.");
        int pairs = _chains.Length - 1;
        var rejection = new double[pairs];
        for (int i = 0; i < pairs; i++) rejection[i] = _roundRejectionSum[i] / _roundScans;
        double[] next = TemperingLadder.Equalize(Betas, rejection);
        ValidateLadder(next, Path);
        for (int i = 0; i < _chains.Length; i++) _chains[i].SetTemperature(new Temperature(next[i], Path));
        Array.Clear(_roundRejectionSum);
        _roundScans = 0;
        return TemperingLadder.Barrier(rejection);
    }

    // IChain view: the ladder advances by one swept Step; its acceptance is the cold replica's own move acceptance
    // (swaps move state via SetState and are tracked separately as SwapAccepts/SwapAttempts).
    long IChain<TState>.Accepted => _chains[0].Accepted;
    long IChain<TState>.Attempts => _chains[0].Attempts;

    /// <summary>Geometric inverse-temperature ladder from 1 (cold) down to <paramref name="betaMin"/>.</summary>
    public static double[] GeometricLadder(int levels, double betaMin)
    {
        if (levels < 2) throw new ArgumentOutOfRangeException(nameof(levels), "Need at least two levels.");
        if (!(betaMin > 0.0 && betaMin < 1.0)) throw new ArgumentOutOfRangeException(nameof(betaMin), "β_min ∈ (0,1).");
        var b = new double[levels];
        for (int i = 0; i < levels; i++)
            b[i] = Math.Pow(betaMin, i / (double)(levels - 1));   // 1 … β_min
        return b;
    }

    // Pair i's log swap ratio on the current states, with its rejection probability recorded.
    private double Record(int i)
    {
        ReversibleJumpChain<TState> a = _chains[i], b = _chains[i + 1];
        double logA = a.Temperature.LogDensity(b.CurrentEvaluation) + b.Temperature.LogDensity(a.CurrentEvaluation)
                    - a.CurrentLogDensity - b.CurrentLogDensity;
        double accept = logA >= 0.0 ? 1.0 : logA < 0.0 ? Math.Exp(logA) : 0.0;   // NaN: no swap
        _rejectionSum[i] += 1.0 - accept;
        _roundRejectionSum[i] += 1.0 - accept;
        return logA;
    }

    private void Attempt(int i, double logA)
    {
        SwapAttempts++;
        _pairAttempts[i]++;
        // Drawn on every attempt, so the swap stream does not advance by the decision: between states whose targets
        // agree to rounding, rounding decides a tie's sign, and it must not decide whether a uniform is consumed.
        double u = _swapRng.NextDouble();
        if (!(logA >= 0.0 || u < Math.Exp(logA))) return;

        SwapAccepts++;
        _pairAccepts[i]++;
        ReversibleJumpChain<TState> a = _chains[i], b = _chains[i + 1];
        TState xa = a.Current, xb = b.Current;
        ITargetEvaluation ea = a.CurrentEvaluation, eb = b.CurrentEvaluation;
        a.SetState(xb, eb);
        b.SetState(xa, ea);
        (_labelAt[i], _labelAt[i + 1]) = (_labelAt[i + 1], _labelAt[i]);
    }

    // A state at the hot end that was last at the cold end completes a round trip; one at the cold end that was last
    // at the hot end is an annealed restart. A state that has not yet been to the hot end is at neither.
    private void TrackEnds()
    {
        int hot = _chains.Length - 1;
        int atHot = _labelAt[hot], atCold = _labelAt[0];
        if (_lastEnd[atHot] == -1) RoundTrips++;
        _lastEnd[atHot] = +1;
        if (_lastEnd[atCold] == +1)
        {
            AnnealedRestarts++;
            _lastEnd[atCold] = -1;
        }
    }

    private static void ValidateLadder(double[] betas, TemperingPath path)
    {
        ArgumentNullException.ThrowIfNull(betas);
        if (betas.Length < 2) throw new ArgumentException("Tempering needs at least two levels.", nameof(betas));
        if (betas[0] != 1.0) throw new ArgumentException("The cold rung must be exactly β = 1.", nameof(betas));
        for (int i = 1; i < betas.Length; i++)
            if (!(betas[i] < betas[i - 1]))
                throw new ArgumentException("The ladder must be strictly descending.", nameof(betas));
        double hot = betas[^1];
        if (!(path == TemperingPath.Powered ? hot > 0.0 : hot >= 0.0))
            throw new ArgumentException(path == TemperingPath.Powered
                ? "On the powered path every β must be positive."
                : "Every β must be non-negative.", nameof(betas));
    }
}
