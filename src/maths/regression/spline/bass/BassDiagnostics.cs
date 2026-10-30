using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bars;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A BASS run's account (the contract shapes), over the sampling phase: everything after burn-in.
/// <list type="bullet">
/// <item><b>Convergence.</b> Rank-normalized R̂ with bulk and tail ESS (<see cref="ChainDiagnostics"/>) for M, τ (null
/// under fixed g, which holds it), each draw's log posterior, σ², and the fitted value f = a₀ + Σⱼ aⱼBⱼ of each draw at
/// <see cref="Fitted"/>'s points. <see cref="SamplesUsed"/> and <see cref="StoppedEarly"/> read the R̂ stop.</item>
/// <item><b>The verdict</b> (D10): the run's <see cref="ConvergenceVerdict"/> on three groups, <see cref="PredictiveGroup"/>
/// (σ² and the fitted values), <see cref="LevelGroup"/> (the log posterior) and <see cref="StructuralGroup"/> (M, and τ
/// when it is sampled), against the run's R̂ target or 1.01; and beside it the <see cref="ModeSeparation"/> readings,
/// which no stop reads.</item>
/// <item><b>Moves.</b> Attempts and acceptances per move type, and the moves' pooled acceptance rate; and per move kind
/// (<see cref="MoveKinds"/>), attempts, acceptances, null moves and refusals, F-M2's local steps and fixed-M family beside
/// D01's birth, death and redraw.</item>
/// <item><b>Support.</b> Refused candidates by reason (cap, support, numerical, variance inflation, exact fit), the
/// births whose draw came out the intercept, and each candidate term's support outcome by its orders, which keeps D01's
/// tilt visible.</item>
/// <item><b>The cap.</b> The share of retained draws at M_max (D02).</item>
/// <item><b>Tempering.</b> The ladders' readout when the run was tempered (<see cref="BassTemperingReadout"/>), null
/// otherwise. A tempered chain's moves, refusals and support account are its b = 1 replica's.</item>
/// <item><b>SMC</b> (stage 8). With <see cref="BassSampling.Smc"/> the population's readout is <see cref="Smc"/>
/// (<see cref="BassSmcReadout"/>), and the chain readings do not apply: a population has no order for R̂ or an
/// autocorrelation ESS to read, so every functional's diagnostic is unavailable, the verdict undetermined in every group,
/// and <see cref="Chains"/> empty. <see cref="ModeSeparation"/> reads the replicates' weighted means; the moves, support
/// and refusal tallies sum every particle's sweeps; <see cref="MassAtCap"/> is the weighted share at M_max; and
/// <see cref="SamplesUsed"/> is the particles per replicate.</item>
/// </list>
/// The totals are the chains' sums, and <see cref="Chains"/> keeps each chain's own, with its start.
/// </summary>
public sealed record BassDiagnostics(
    ConvergenceDiagnostic Count,
    ConvergenceDiagnostic? Tau,
    ConvergenceDiagnostic LogPosterior,
    ConvergenceDiagnostic Sigma2,
    IReadOnlyList<BassMonitoredValue> Fitted,
    IReadOnlyList<MoveStat> Moves,
    double AcceptanceRate,
    BassMoveTally MoveKinds,
    Tally<BassRefusal> Refusals,
    long EmptyBirths,
    BassSupportTally Support,
    double MassAtCap,
    int SamplesUsed,
    bool StoppedEarly,
    ConvergenceVerdict Verdict,
    BassModeSeparation ModeSeparation,
    IReadOnlyList<BassChainAccount> Chains,
    BassTemperingReadout? Tempering,
    BassSmcReadout? Smc = null)
{
    /// <summary>The predictive group: σ² and the fitted values. Do the chains predict the same surface and noise?</summary>
    public const string PredictiveGroup = "predictive";

    /// <summary>The level group: each draw's log posterior. Do the chains sit in regions of the same mass?</summary>
    public const string LevelGroup = "level";

    /// <summary>The structural group: M, and τ when it is sampled. Do the chains agree on the representation? In an
    /// overcomplete hinge basis M is not identified by the surface, so disagreement here beside agreement in the other
    /// groups reads as non-identifiability.</summary>
    public const string StructuralGroup = "structural";

    /// <summary>The groups an R̂ stop waits on under <paramref name="scope"/>.</summary>
    public static IReadOnlyList<string> StoppingGroups(BassStoppingScope scope) => scope switch
    {
        BassStoppingScope.All => [PredictiveGroup, LevelGroup, StructuralGroup],
        BassStoppingScope.Predictive => [PredictiveGroup],
        _ => throw new ArgumentOutOfRangeException(nameof(scope)),
    };

    /// <summary>
    /// The points whose fitted values are monitored when the config declares none
    /// (<see cref="BassSampling.MonitoredPoints"/>): up to three training cases spread through the data's order, at its
    /// quartile positions; for a direct functional response, each at the middle output location.
    /// </summary>
    public static IReadOnlyList<BassMonitoredPoint> DefaultPoints(BassResolved resolved, bool functional)
    {
        ArgumentNullException.ThrowIfNull(resolved);
        return DefaultPoints(resolved.Cases, functional ? resolved.Locations : null);
    }

    /// <summary>The defaults for <paramref name="cases"/> cases and, for a direct functional response,
    /// <paramref name="locations"/> output locations.</summary>
    public static IReadOnlyList<BassMonitoredPoint> DefaultPoints(int cases, int? locations)
    {
        if (cases < 1) throw new ArgumentOutOfRangeException(nameof(cases));
        int? location = locations is int m ? (m - 1) / 2 : null;
        var points = new List<BassMonitoredPoint>(3);
        foreach (double q in (ReadOnlySpan<double>)[0.25, 0.5, 0.75])
        {
            int i = (int)Math.Round((cases - 1) * q, MidpointRounding.AwayFromZero);
            if (points.Count == 0 || points[^1].Case != i) points.Add(new BassMonitoredPoint(i, location));
        }
        return points;
    }
}

/// <summary>
/// What a tempered run's ladders did over the sampling phase (D05, D10). Each chain froze its own rungs, at its own count
/// when the ladder sizes itself, so the pairs are read per chain, beside each chain's frozen ladder:
/// <list type="bullet">
/// <item>the round-trip and annealed-restart rates per scan, the chains' mean;</item>
/// <item><see cref="PairRejection"/>: each chain's adjacent pairs' mean swap rejection, cold first (Rao–Blackwellized:
/// 1 − α every scan, attempted or not), and <see cref="PairAcceptance"/>, their swap acceptance (accepted over attempted;
/// NaN for a pair DEO never reached, as in a one-scan sampling phase);</item>
/// <item><see cref="Barrier"/>, the chains' mean Λ̂, each the sum of its pairs' rejections, and
/// <see cref="MaxPairRejection"/>, the largest of any pair, which nears 1 on a ladder whose rungs barely communicate;</item>
/// <item><see cref="IdealRoundTripRate"/>, the chains' mean of the rate Syed et al.'s (2022) Theorem 1 predicts from each
/// chain's rejections under ideal local exploration, against which the observed rate measures the replicas' own mixing;</item>
/// <item><see cref="Ladders"/> and <see cref="Levels"/>, each chain's frozen ladder, cold first, and its count;
/// <see cref="AtMaxLevels"/>, whether a self-sizing chain's ladder reached <see cref="TemperingSpec.MaxLevels"/>;</item>
/// <item>with the reference rung, its sampler's attempts per exact draw, the rung's cost at equal work (NaN without it).</item>
/// </list>
/// </summary>
public sealed record BassTemperingReadout(
    TemperingPath Path, TemperedForm Form, SwapSchedule Schedule, bool ReferenceRung,
    double RoundTripsPerScan, double AnnealedRestartsPerScan, double[][] PairRejection, double[][] PairAcceptance, double Barrier,
    double MaxPairRejection, double IdealRoundTripRate, double[][] Ladders, int[] Levels, bool AtMaxLevels, double ReferenceAttemptsPerDraw)
{
    // The readout of `ladders`, one per chain, each past its burn-in.
    internal static BassTemperingReadout Of(IReadOnlyList<BassLadder> ladders)
    {
        TemperingSpec spec = ladders[0].Spec;
        int c = ladders.Count;
        double trips = 0.0, restarts = 0.0, referenceAttempts = 0.0, barrier = 0.0, ideal = 0.0, worst = 0.0;
        var rejection = new double[c][];
        var acceptance = new double[c][];
        var rungs = new double[c][];
        var levels = new int[c];
        bool atMax = false;
        for (int k = 0; k < c; k++)
        {
            ParallelTempering<BassState> ladder = ladders[k].Tempering;
            BassLadderMarks burn = ladders[k].BurnMarks ?? throw new InvalidOperationException("A ladder's burn-in has not ended.");
            double scans = ladder.Scans - burn.Scans;
            trips += (ladder.RoundTrips - burn.RoundTrips) / scans;
            restarts += (ladder.AnnealedRestarts - burn.AnnealedRestarts) / scans;
            IReadOnlyList<SwapPairStat> stats = ladder.PairStats;
            rejection[k] = new double[stats.Count];
            acceptance[k] = new double[stats.Count];
            for (int i = 0; i < stats.Count; i++)
            {
                rejection[k][i] = (stats[i].MeanRejection * ladder.PairScans - burn.Rejection[i]) / scans;
                long attempts = stats[i].Attempts - burn.PairAttempts[i];
                acceptance[k][i] = attempts > 0 ? (double)(stats[i].Accepts - burn.PairAccepts[i]) / attempts : double.NaN;
                worst = Math.Max(worst, rejection[k][i]);
            }
            barrier += TemperingLadder.Barrier(rejection[k]);
            if (spec.Schedule != SwapSchedule.Sequential) ideal += TemperingLadder.IdealRoundTripRate(rejection[k], spec.Schedule);
            if (ladders[k].Reference is { } reference) referenceAttempts += (reference.Attempts - burn.ReferenceAttempts) / scans;
            rungs[k] = ladder.InverseTemperatures;
            levels[k] = rungs[k].Length;
            atMax |= ladders[k].SizesItself && levels[k] == spec.MaxLevels;
        }
        return new BassTemperingReadout(spec.Path, ladders[0].Model.Target.Form, spec.Schedule, spec.ReferenceRung, trips / c, restarts / c,
                                        rejection, acceptance, barrier / c, worst,
                                        spec.Schedule == SwapSchedule.Sequential ? double.NaN : ideal / c, rungs, levels, atMax,
                                        spec.ReferenceRung ? referenceAttempts / c : double.NaN);
    }
}

/// <summary>
/// An SMC fit's population readout (stage 8), kept apart from chain R̂ and autocorrelation ESS.
/// <list type="bullet">
/// <item><see cref="Replicates"/>: each independent run's stage reports (ESS, conditional ESS, resampling, distinct
/// ancestors, families, acceptance and the normalizer's increments with its genealogy error), where it ended, its
/// normalizer, and its final population's Kish ESS and families.</item>
/// <item><see cref="Spread"/>: the between-run reading for each monitored quantity, the SD over replicates of their
/// weighted posterior means in units of the pooled posterior sd. Split R̂ means nothing for a population, which has no
/// order to split, and an unsplit R̂ would invite comparison with D10's 1.01, which assumes autocorrelated chains. No stop
/// reads it.</item>
/// <item>Particles share ancestry, so a Kish ESS, or a Monte Carlo error that treats them as independent (the pointwise
/// checks'), reads optimistic: each replicate's families say how many lineages remain, and the replicates' spread is the
/// honest error.</item>
/// </list>
/// </summary>
public sealed record BassSmcReadout(IReadOnlyList<BassSmcReplicate> Replicates, IReadOnlyList<BassSmcSpread> Spread)
{
    /// <summary>Whether every replicate reached b = 1 within its stage budget: only then are the draws the posterior.</summary>
    public bool ReachedPosterior
    {
        get
        {
            foreach (BassSmcReplicate r in Replicates)
                if (r.InverseTemperature != 1.0) return false;
            return true;
        }
    }

    /// <summary>The spread of the quantity named <paramref name="name"/> ("M", "tau", "log_posterior", "sigma2",
    /// "fitted[i]" or "fitted[i,r]").</summary>
    public BassSmcSpread Quantity(string name)
    {
        ArgumentNullException.ThrowIfNull(name);
        foreach (BassSmcSpread q in Spread)
            if (q.Name == name) return q;
        throw new KeyNotFoundException($"No monitored quantity is named \"{name}\".");
    }
}

/// <summary>
/// One SMC replicate's account: its stream's seed; its stage reports, the initial draw first; the b it ended at, 1 unless
/// its budget ran out; its log normalizer ln Ẑ with the genealogy error (NaN unless every resampling was multinomial);
/// its final population's Kish ESS and surviving families; the attempts its prior draws took; its kernel's account
/// (<see cref="BassSmcAccount"/>); and its weighted mass at M_max.
/// </summary>
public sealed record BassSmcReplicate(
    int Seed, IReadOnlyList<SmcStageReport> Stages, double InverseTemperature, double LogNormalizer, double LogNormalizerError,
    double EffectiveSampleSize, int Families, long PriorAttempts, BassSmcAccount Account, double MassAtCap)
{
    /// <summary>The stages taken after the initial draw.</summary>
    public int StageCount => Stages.Count - 1;
}

/// <summary>
/// One monitored quantity's between-run reading: each replicate's weighted posterior mean, the pooled posterior's mean and
/// sd over every replicate's draws (each replicate weighing 1/R), and <see cref="Spread"/>, the sample SD of the
/// replicates' means over that sd. NaN with one replicate; 0 for a quantity every draw holds at one value.
/// </summary>
public sealed record BassSmcSpread(string Name, double[] ReplicateMeans, double PosteriorMean, double PosteriorSd, double Spread);

/// <summary>
/// How far apart the chains sit (D10), reported beside the verdict: once chains sit in separate modes, R̂ and bulk ESS
/// saturate, and these measure the distance. No stop reads them; they have no calibrated threshold yet.
/// <list type="bullet">
/// <item><see cref="LogPosteriorGap"/>: the largest chain's mean log posterior minus the smallest's, in nats;</item>
/// <item><see cref="FittedDistance"/>: the largest Euclidean distance between two chains' mean fitted values at the
/// monitored points.</item>
/// </list>
/// Both are NaN with one chain.
/// </summary>
public sealed record BassModeSeparation(double LogPosteriorGap, double FittedDistance)
{
    // The readings from each chain's mean log posterior and mean fitted values, one row per chain.
    internal static BassModeSeparation Of(IReadOnlyList<double> logPosterior, IReadOnlyList<double[]> fitted)
    {
        int c = logPosterior.Count;
        if (c < 2) return new BassModeSeparation(double.NaN, double.NaN);
        double lo = double.PositiveInfinity, hi = double.NegativeInfinity, widest = 0.0;
        for (int i = 0; i < c; i++)
        {
            lo = Math.Min(lo, logPosterior[i]);
            hi = Math.Max(hi, logPosterior[i]);
            for (int j = i + 1; j < c; j++)
            {
                double squares = 0.0;
                for (int p = 0; p < fitted[i].Length; p++)
                {
                    double d = fitted[i][p] - fitted[j][p];
                    squares += d * d;
                }
                widest = Math.Max(widest, Math.Sqrt(squares));
            }
        }
        return new BassModeSeparation(hi - lo, widest);
    }
}

/// <summary>A monitored point: a training case and, for a direct functional response, an output location.</summary>
public readonly record struct BassMonitoredPoint(int Case, int? Location);

/// <summary>The convergence of the fitted value at <paramref name="Point"/>.</summary>
public sealed record BassMonitoredValue(BassMonitoredPoint Point, ConvergenceDiagnostic Diagnostic);

/// <summary>A chain's start: an exact draw of the prior on F, its terms and τ, and the attempts the draw took.</summary>
public sealed record BassStart(IReadOnlyList<BassTerm> Terms, double Tau, long Attempts);

/// <summary>One chain's account over the sampling phase, with its start: moves, by type and by kind, refusals, empty
/// births, support by order and its share of draws at the cap.</summary>
public sealed record BassChainAccount(
    BassStart Start, IReadOnlyList<MoveStat> Moves, BassMoveTally MoveKinds, Tally<BassRefusal> Refusals, long EmptyBirths,
    BassSupportTally Support, double MassAtCap);

/// <summary>One move kind's account: its proposals attempted, accepted, null and refused; the rest were rejected on their
/// ratio.</summary>
public readonly record struct BassKindCount(BassMoveKind Kind, long Attempts, long Accepted, long Null, long Refused)
{
    /// <summary>Proposals the ratio rejected: neither accepted, null nor refused.</summary>
    public long Rejected => Attempts - Accepted - Null - Refused;

    /// <summary>Accepted over attempted; NaN without an attempt.</summary>
    public double AcceptanceRate => Attempts > 0 ? (double)Accepted / Attempts : double.NaN;
}

/// <summary>
/// The structure moves by kind (F-M2's Phase 1): D01's birth, death and redraw, lever 2's local step and toggle, and the
/// fixed-M family's prune, grow and variable change, each with its attempts, acceptances, null moves and refusals.
/// <list type="bullet">
/// <item>An attempt is a proposal the chain weighs: a birth whose draw came out the intercept is no attempt
/// (<see cref="BassStructureMoves.EmptyBirths"/>), nor is a change with no term to pick.</item>
/// <item>A null move proposes nothing: a step off the knots, a toggle to an empty or full subset, a prune from an empty
/// part, a grow into a full one, a variable change with no unused variable, the family at M = 0.</item>
/// <item>A refusal is a candidate outside F, by the design's or the target's reason (tallied by reason in
/// <see cref="BassStructureMoves.Refusals"/>), or a prune of a term's only factor, which would leave the intercept.</item>
/// </list>
/// The kinds of one palette move add up to the chain's count for it (<see cref="MoveStat"/>): the redraw, the local step
/// and the toggle to the change's. A redraw that draws its factor back is the same state, accepted, as the chain counts it.
/// </summary>
public sealed class BassMoveTally : IMergeable<BassMoveTally>
{
    private static readonly int Kinds = Enum.GetValues<BassMoveKind>().Length;
    private readonly long[] _attempts = new long[Kinds], _accepted = new long[Kinds], _null = new long[Kinds], _refused = new long[Kinds];

    /// <summary>An empty tally.</summary>
    public BassMoveTally() { }

    private BassMoveTally(BassMoveTally source)
    {
        _attempts = (long[])source._attempts.Clone();
        _accepted = (long[])source._accepted.Clone();
        _null = (long[])source._null.Clone();
        _refused = (long[])source._refused.Clone();
    }

    /// <summary><paramref name="kind"/>'s account.</summary>
    public BassKindCount this[BassMoveKind kind]
        => new(kind, _attempts[(int)kind], _accepted[(int)kind], _null[(int)kind], _refused[(int)kind]);

    /// <summary>Every kind's account, in <see cref="BassMoveKind"/>'s order.</summary>
    public IReadOnlyList<BassKindCount> Entries => [.. Enum.GetValues<BassMoveKind>().Select(k => this[k])];

    internal void Attempted(BassMoveKind kind) => _attempts[(int)kind]++;
    internal void Accepted(BassMoveKind kind) => _accepted[(int)kind]++;
    internal void Nulled(BassMoveKind kind) => _null[(int)kind]++;
    internal void Refused(BassMoveKind kind) => _refused[(int)kind]++;

    /// <inheritdoc/>
    public void MergeFrom(BassMoveTally other)
    {
        ArgumentNullException.ThrowIfNull(other);
        for (int k = 0; k < Kinds; k++)
        {
            _attempts[k] += other._attempts[k];
            _accepted[k] += other._accepted[k];
            _null[k] += other._null[k];
            _refused[k] += other._refused[k];
        }
    }

    /// <summary>An independent copy.</summary>
    public BassMoveTally Copy() => new(this);

    /// <summary>Zeroes every count.</summary>
    public void Clear()
    {
        Array.Clear(_attempts);
        Array.Clear(_accepted);
        Array.Clear(_null);
        Array.Clear(_refused);
    }
}

/// <summary>One part and order's support outcomes: candidate terms of that order in that part, and those refused for
/// support.</summary>
public readonly record struct BassOrderSupport(BassPart Part, int Order, long Attempts, long Refused)
{
    /// <summary>The share refused; NaN without an attempt.</summary>
    public double RefusalRate => Attempts > 0 ? (double)Refused / Attempts : double.NaN;
}

/// <summary>
/// Candidate terms' support outcomes by their orders (D01's tilt): each term a birth or a change proposes counts once at
/// each part's order, and once more as refused there when the design refused it for support. Under the global
/// truncation, high orders fail the support count more often, so they get less mass than under per-order normalizers;
/// the refusal rate by order shows by how much.
/// </summary>
public sealed class BassSupportTally : IMergeable<BassSupportTally>
{
    private readonly long[][] _attempts = new long[3][];
    private readonly long[][] _refused = new long[3][];
    private readonly int[] _minimum = new int[3];   // each part's smallest order

    /// <summary>An empty tally over <paramref name="prior"/>'s parts and orders.</summary>
    public BassSupportTally(BassStructurePrior prior)
    {
        ArgumentNullException.ThrowIfNull(prior);
        for (int p = 0; p < 3; p++)
        {
            int orders = prior.Has((BassPart)p) ? prior.Cap((BassPart)p) + 1 : 0;
            _attempts[p] = new long[orders];
            _refused[p] = new long[orders];
            _minimum[p] = prior.MinimumOrder((BassPart)p);
        }
    }

    private BassSupportTally(BassSupportTally source)
    {
        for (int p = 0; p < 3; p++)
        {
            _attempts[p] = (long[])source._attempts[p].Clone();
            _refused[p] = (long[])source._refused[p].Clone();
            _minimum[p] = source._minimum[p];
        }
    }

    /// <summary>Counts <paramref name="term"/> at each part's order, as refused for support or not.</summary>
    public void Add(BassTerm term, bool refused)
    {
        ArgumentNullException.ThrowIfNull(term);
        for (int p = 0; p < 3; p++)
        {
            if (_attempts[p].Length == 0) continue;
            int order = BassStructurePrior.Order(term, (BassPart)p);
            _attempts[p][order]++;
            if (refused) _refused[p][order]++;
        }
    }

    /// <summary>Candidates whose <paramref name="part"/> has order <paramref name="order"/>.</summary>
    public long Attempts(BassPart part, int order) => _attempts[(int)part][order];

    /// <summary>Those of them refused for support.</summary>
    public long Refused(BassPart part, int order) => _refused[(int)part][order];

    /// <summary>Every part the fit has, by each order its terms can take.</summary>
    public IReadOnlyList<BassOrderSupport> Entries
    {
        get
        {
            var entries = new List<BassOrderSupport>();
            for (int p = 0; p < 3; p++)
                for (int k = _minimum[p]; k < _attempts[p].Length; k++) entries.Add(new BassOrderSupport((BassPart)p, k, _attempts[p][k], _refused[p][k]));
            return entries;
        }
    }

    /// <inheritdoc/>
    public void MergeFrom(BassSupportTally other)
    {
        ArgumentNullException.ThrowIfNull(other);
        for (int p = 0; p < 3; p++)
        {
            if (other._attempts[p].Length != _attempts[p].Length || other._minimum[p] != _minimum[p])
                throw new ArgumentException("The tallies count different orders.", nameof(other));
            for (int k = 0; k < _attempts[p].Length; k++)
            {
                _attempts[p][k] += other._attempts[p][k];
                _refused[p][k] += other._refused[p][k];
            }
        }
    }

    /// <summary>An independent copy.</summary>
    public BassSupportTally Copy() => new(this);

    /// <summary>Zeroes every count.</summary>
    public void Clear()
    {
        for (int p = 0; p < 3; p++)
        {
            Array.Clear(_attempts[p]);
            Array.Clear(_refused[p]);
        }
    }
}
