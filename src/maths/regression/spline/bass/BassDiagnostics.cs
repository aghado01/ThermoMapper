using System;
using System.Collections.Generic;
using Maths.Regression.Spline.Bars;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A BASS run's account (the contract shapes), over the sampling phase: everything after burn-in.
/// <list type="bullet">
/// <item><b>Convergence.</b> Rank-normalized R̂ with bulk and tail ESS (<see cref="ChainDiagnostics"/>) for M, τ (null
/// under fixed g, which holds it), σ², and the fitted value f = a₀ + Σⱼ aⱼBⱼ of each draw at <see cref="Fitted"/>'s
/// points. <see cref="SamplesUsed"/> and <see cref="StoppedEarly"/> read the R̂ stop.</item>
/// <item><b>Moves.</b> Attempts and acceptances per move type, and the moves' pooled acceptance rate.</item>
/// <item><b>Support.</b> Refused candidates by reason (cap, support, numerical, variance inflation, exact fit), the
/// births whose draw came out the intercept, and each candidate term's support outcome by its orders, which keeps D01's
/// tilt visible.</item>
/// <item><b>The cap.</b> The share of retained draws at M_max (D02).</item>
/// <item><b>Tempering.</b> The ladders' readout when the run was tempered (<see cref="BassTemperingReadout"/>), null
/// otherwise. A tempered chain's moves, refusals and support account are its b = 1 replica's.</item>
/// </list>
/// The totals are the chains' sums, and <see cref="Chains"/> keeps each chain's own, with its start.
/// </summary>
public sealed record BassDiagnostics(
    ConvergenceDiagnostic Count,
    ConvergenceDiagnostic? Tau,
    ConvergenceDiagnostic Sigma2,
    IReadOnlyList<BassMonitoredValue> Fitted,
    IReadOnlyList<MoveStat> Moves,
    double AcceptanceRate,
    Tally<BassRefusal> Refusals,
    long EmptyBirths,
    BassSupportTally Support,
    double MassAtCap,
    int SamplesUsed,
    bool StoppedEarly,
    IReadOnlyList<BassChainAccount> Chains,
    BassTemperingReadout? Tempering)
{
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
/// What a tempered run's ladders did over the sampling phase (D05, D10), averaged over chains as BARS's readout is: the
/// round-trip and annealed-restart rates per scan; each adjacent pair's mean swap rejection, cold first (Rao–Blackwellized:
/// 1 − α every scan, attempted or not), and its swap acceptance (accepted over attempted swaps, the chains' counts summed;
/// NaN for a pair DEO never reached, as in a one-scan sampling phase);
/// the barrier estimate Λ̂, the rejections' sum; the round-trip rate Syed et al.'s (2022) Theorem 1 predicts under ideal
/// local exploration, against which the observed rate measures the replicas' own mixing; each chain's frozen ladder, cold
/// first; and with the reference rung, its sampler's attempts per exact draw, the rung's cost at equal work (NaN without
/// it). The rungs are those burn-in left (<see cref="TemperingSpec.Adaptive"/>), so each chain's may differ.
/// </summary>
public sealed record BassTemperingReadout(
    TemperingPath Path, TemperedForm Form, SwapSchedule Schedule, bool ReferenceRung,
    double RoundTripsPerScan, double AnnealedRestartsPerScan, double[] PairRejection, double[] PairAcceptance, double Barrier,
    double IdealRoundTripRate, double[][] Ladders, double ReferenceAttemptsPerDraw)
{
    // The readout of `ladders`, one per chain, each past its burn-in.
    internal static BassTemperingReadout Of(IReadOnlyList<BassLadder> ladders)
    {
        TemperingSpec spec = ladders[0].Spec;
        int pairs = spec.Levels - 1, c = ladders.Count;
        double trips = 0.0, restarts = 0.0, referenceAttempts = 0.0;
        var rejection = new double[pairs];
        var attempts = new long[pairs];
        var accepts = new long[pairs];
        var rungs = new double[c][];
        for (int k = 0; k < c; k++)
        {
            ParallelTempering<BassState> ladder = ladders[k].Tempering;
            BassLadderMarks burn = ladders[k].BurnMarks ?? throw new InvalidOperationException("A ladder's burn-in has not ended.");
            double scans = ladder.Scans - burn.Scans;
            trips += (ladder.RoundTrips - burn.RoundTrips) / scans;
            restarts += (ladder.AnnealedRestarts - burn.AnnealedRestarts) / scans;
            IReadOnlyList<SwapPairStat> stats = ladder.PairStats;
            for (int i = 0; i < pairs; i++)
            {
                rejection[i] += (stats[i].MeanRejection * ladder.Scans - burn.Rejection[i]) / scans;
                attempts[i] += stats[i].Attempts - burn.PairAttempts[i];
                accepts[i] += stats[i].Accepts - burn.PairAccepts[i];
            }
            if (ladders[k].Reference is { } reference) referenceAttempts += (reference.Attempts - burn.ReferenceAttempts) / scans;
            rungs[k] = ladder.Betas;
        }
        var acceptance = new double[pairs];
        for (int i = 0; i < pairs; i++)
        {
            rejection[i] /= c;
            acceptance[i] = attempts[i] > 0 ? (double)accepts[i] / attempts[i] : double.NaN;
        }
        double ideal = spec.Schedule == SwapSchedule.Sequential ? double.NaN : TemperingLadder.IdealRoundTripRate(rejection, spec.Schedule);
        return new BassTemperingReadout(spec.Path, ladders[0].Model.Target.Form, spec.Schedule, spec.ReferenceRung, trips / c, restarts / c,
                                        rejection, acceptance, TemperingLadder.Barrier(rejection), ideal, rungs,
                                        spec.ReferenceRung ? referenceAttempts / c : double.NaN);
    }
}

/// <summary>A monitored point: a training case and, for a direct functional response, an output location.</summary>
public readonly record struct BassMonitoredPoint(int Case, int? Location);

/// <summary>The convergence of the fitted value at <paramref name="Point"/>.</summary>
public sealed record BassMonitoredValue(BassMonitoredPoint Point, ConvergenceDiagnostic Diagnostic);

/// <summary>A chain's start: an exact draw of the prior on F, its terms and τ, and the attempts the draw took.</summary>
public sealed record BassStart(IReadOnlyList<BassTerm> Terms, double Tau, long Attempts);

/// <summary>One chain's account over the sampling phase, with its start: moves, refusals, empty births, support by
/// order and its share of draws at the cap.</summary>
public sealed record BassChainAccount(
    BassStart Start, IReadOnlyList<MoveStat> Moves, Tally<BassRefusal> Refusals, long EmptyBirths, BassSupportTally Support,
    double MassAtCap);

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
