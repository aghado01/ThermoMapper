using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;

namespace Maths.Regression.Conformal;

/// <summary>A conditional diagnostic interval for the randomized PIT's probability, not a fitted response interval.</summary>
public readonly record struct ConformalPitBounds(double Lower, double Upper);

/// <summary>
/// One CDF-score ordering's errors and counts. Errors are on the probability scale, not density or log-density errors.
/// Fragile and PIT arrays are indexed by rank level k, including 0 and m; the two tails' fragile units are counted once.
/// Position arrays contain the stratum's training units in their unit order, followed by the candidate.
/// At a flagged candidate turn no band or PIT count is read. Flags are density-weight proxy diagnostics.
/// </summary>
public sealed record CdfMarginReading(double[] PositionErrors, double[]? CarriedErrors, bool CandidateFlagged, bool[]? FlaggedUnits,
    int Flagged, int[]? Fragile, int[]? DegenerateAccepted, int[]? PitLowerCounts, int[]? PitUpperCounts,
    int InflatedLowerRank, int InflatedUpperRank);

/// <summary>
/// Reduces CDF-score errors to two-tail union and PIT counts over a fixed/refined response grid. All numerical readouts are
/// diagnostics conditional on simultaneous error containment and grid coverage; Kept also requires proxy containment of
/// every population own-turn flag. Conjugate checks found proxy misses. These numbers are not guaranteed coverage bounds.
/// </summary>
public sealed class ConformalCdfMargin
{
    public int Units { get; }
    public MarginBounds? Bounds { get; }
    public InflationBounds Inflation { get; }
    public double Multiple { get; }
    public double ErrorLevel { get; }
    /// <summary>Training units' carried bounds before their held-ordering errors at the candidate are included.</summary>
    public double[] UnitCarriedErrors { get; }
    public IReadOnlyList<CdfMarginReading> Readings { get; }
    public int[] InflatedLowerRanks { get; }
    public int[] InflatedUpperRanks { get; }
    /// <summary>Largest union band plus D's inflated-accepted units, over kept grid values and probes, by level.</summary>
    public int[]? Counts { get; }
    public int? FlaggedCount { get; }
    /// <summary>Minimum lower and maximum upper PIT counts over kept grid values and probes, by level.</summary>
    public int[]? PitLowerCounts { get; }
    public int[]? PitUpperCounts { get; }
    public int RefinementProbes { get; }

    internal ConformalCdfMargin(int units, ConformalOptions options, double multiple, double[] unitCarried,
        CdfMarginReading[] readings, List<CdfMarginReading> probes)
    {
        Units = units; Bounds = options.MarginBounds; Inflation = options.InflationBounds; Multiple = multiple;
        ErrorLevel = options.ErrorLevel; UnitCarriedErrors = unitCarried; Readings = readings;
        InflatedLowerRanks = readings.Select(r => r.InflatedLowerRank).ToArray();
        InflatedUpperRanks = readings.Select(r => r.InflatedUpperRank).ToArray();
        RefinementProbes = probes.Count;
        if (Bounds is null) return;
        Counts = new int[units + 1]; PitLowerCounts = new int[units + 1]; PitUpperCounts = new int[units + 1];
        Array.Fill(PitLowerCounts, units);
        int flagged = 0, kept = 0;
        foreach (var r in readings.Concat(probes))
        {
            flagged = Math.Max(flagged, r.Flagged);
            if (r.CandidateFlagged) continue;
            kept++;
            for (int k = 0; k <= units; k++)
            {
                Counts[k] = Math.Max(Counts[k], r.Fragile![k] + (r.DegenerateAccepted?[k] ?? 0));
                PitLowerCounts[k] = Math.Min(PitLowerCounts[k], r.PitLowerCounts![k]);
                PitUpperCounts[k] = Math.Max(PitUpperCounts[k], r.PitUpperCounts![k]);
            }
        }
        if (kept == 0) Array.Clear(PitLowerCounts);
        FlaggedCount = Bounds == MarginBounds.Kept ? flagged : null;
    }

    internal double CoverageUnderError(int k, int units, bool inflated)
    {
        int d = FlaggedCount ?? 0;
        if (inflated) return 1.0 - (double)(2 * k + d) / units - ErrorLevel;
        return Counts is null ? double.NaN : 1.0 - (double)(2 * k + Counts[k] + d) / units
            - (Bounds == MarginBounds.Kept ? 2 : 1) * ErrorLevel;
    }

    /// <summary>
    /// Conditional diagnostic for P(PIT ≤ k/m), from doubled-shift counts. Kept charges the largest D again for unread
    /// orderings, in both ends of the grid reading, and 2ν; other modes pay ν. Exact endpoints are 0 and 1.
    /// </summary>
    public ConformalPitBounds Pit(int k)
    {
        if ((uint)k > (uint)Units) throw new ArgumentOutOfRangeException(nameof(k));
        if (PitLowerCounts is null) return new(double.NaN, double.NaN);
        if (k == 0 || k == Units) return new((double)k / Units, (double)k / Units);
        double charge = (double)(FlaggedCount ?? 0) / Units + (Bounds == MarginBounds.Kept ? 2 : 1) * ErrorLevel;
        return new(Math.Max(0, (double)PitLowerCounts[k] / Units - charge), Math.Min(1, (double)PitUpperCounts![k] / Units + charge));
    }

    /// <summary>The PIT diagnostic interpolated between rank levels, as pitCdf_interp; response-grid conditions still apply.</summary>
    public ConformalPitBounds PitAt(double u)
    {
        if (!(u >= 0 && u <= 1)) throw new ArgumentOutOfRangeException(nameof(u));
        int k = Math.Min((int)(u * Units), Units - 1);
        double f = u * Units - k;
        var a = Pit(k); var b = Pit(k + 1);
        return new((1 - f) * a.Lower + f * b.Lower, (1 - f) * a.Upper + f * b.Upper);
    }

    /// <summary>The largest diagnostic departure from uniformity over rank levels; interpolation adds no larger departure.</summary>
    public double PitDeparture
    {
        get
        {
            if (PitLowerCounts is null) return double.NaN;
            double most = 0;
            for (int k = 0; k <= Units; k++)
            {
                var b = Pit(k); double u = (double)k / Units;
                most = Math.Max(most, Math.Max(u - b.Lower, b.Upper - u));
            }
            return most;
        }
    }
}

internal static class CdfMarginCounts
{
    private static bool AtMost(CdfScore a, CdfScore b, double shift)
        => shift == 0 ? a.CompareTo(b) <= 0 : a.Minus(b) <= shift;
    private static bool Below(CdfScore a, CdfScore b, double shift, bool strict)
        => shift == 0 ? (strict ? a.CompareTo(b) < 0 : a.CompareTo(b) <= 0)
                      : (strict ? a.Minus(b) < -shift : a.Minus(b) <= -shift);

    internal static CdfMarginReading Read(CdfScore[] scores, double[] held, double[]? carried, bool ownFlag, bool[]? flags,
        MarginBounds? mode, InflationBounds inflation)
    {
        int m = scores.Length, candidate = m - 1;
        double maxHeld = held.Max();
        double[] inflationErrors = inflation == InflationBounds.Uniform ? Enumerable.Repeat(maxHeld, m).ToArray() : held;
        int low = 0, high = 0;
        for (int j = 0; j < m; j++)
        {
            double e = inflationErrors[j] + inflationErrors[candidate];
            if (AtMost(scores[j], scores[candidate], e)) low++;
            if (AtMost(scores[candidate], scores[j], e)) high++;
        }
        int d = flags?.Count(f => f) ?? 0;
        if (mode is null || ownFlag)
            return new(held, carried, ownFlag, flags, d, null, null, null, null, low, high);
        var band = new int[m + 1]; var accepted = mode == MarginBounds.Kept ? new int[m + 1] : null;
        var pitLow = new int[m + 1]; var pitHigh = new int[m + 1];
        double[] bounds = carried!;
        if (mode == MarginBounds.Uniform) bounds = Enumerable.Repeat(bounds.Max(), m).ToArray();
        for (int u = 0; u < m; u++)
        {
            int inflatedLow = 0, inflatedHigh = 0, shiftLow = 1, shiftHigh = 1, l = 0, upper = 0;
            for (int j = 0; j < m; j++)
            {
                double e = bounds[j] + bounds[u], doubled = 2 * e;
                if (AtMost(scores[j], scores[u], e)) inflatedLow++;
                if (AtMost(scores[u], scores[j], e)) inflatedHigh++;
                if (j != u && Below(scores[j], scores[u], doubled, false)) shiftLow++;
                if (j != u && Below(scores[u], scores[j], doubled, false)) shiftHigh++;
                if (Below(scores[j], scores[u], doubled, true)) l++;
                if (AtMost(scores[j], scores[u], doubled)) upper++;
            }
            bool flagged = flags?[u] ?? false;
            for (int k = 0; k <= m; k++)
            {
                if ((inflatedLow > k && shiftLow <= k) || (inflatedHigh > k && shiftHigh <= k)) band[k]++;
                if (accepted is not null && flagged && inflatedLow > k && inflatedHigh > k) accepted[k]++;
                if (!flagged && upper <= k) pitLow[k]++;
                if (flagged || l < k) pitHigh[k]++;
            }
        }
        return new(held, carried, ownFlag, flags, d, band, accepted, pitLow, pitHigh, low, high);
    }
}

// A CDF field is held by both tails. Six dense products accumulate its self-normalized influence variance under a batch
// of normalized weight columns. Near cancellation or underflow, a scaled norm of the centered influences is read directly.
internal sealed class CdfMomentField
{
    private readonly double[] _lower, _upper, _lowerSquared, _upperSquared;
    private readonly int _n, _s;
    internal CdfMomentField(double[,] logCdf, int[] rows, int[] columns)
    {
        _n = columns.Length; _s = rows.Length;
        _lower = new double[_n * _s]; _upper = new double[_n * _s];
        _lowerSquared = new double[_n * _s]; _upperSquared = new double[_n * _s];
        for (int i = 0; i < _n; i++)
            for (int d = 0; d < _s; d++)
            {
                var (f, q) = CdfScore.Of(logCdf[rows[d], columns[i]]); int at = i * _s + d;
                _lower[at] = f; _upper[at] = q; _lowerSquared[at] = f * f; _upperSquared[at] = q * q;
            }
    }

    internal (CdfScore[] Scores, double[] Errors) Measure(double[] weights, int count, double[] efficiency)
    {
        int size = _n * count;
        if (size == 0) return ([], []);
        var squared = new double[_s * count]; var q2 = new double[count];
        for (int d = 0; d < _s; d++) for (int b = 0; b < count; b++)
            q2[b] += squared[d * count + b] = weights[d * count + b] * weights[d * count + b];
        var w = DenseProduct.Pack(weights, _s, count); var w2 = DenseProduct.Pack(squared, _s, count);
        double[] Product(double[] left, PackedColumns right) { var p = new double[size]; DenseProduct.Multiply(left, _n, right, p); return p; }
        double[] lo = Product(_lower, w), hi = Product(_upper, w), lo1 = Product(_lower, w2), hi1 = Product(_upper, w2),
                 lo2 = Product(_lowerSquared, w2), hi2 = Product(_upperSquared, w2);
        var scores = new CdfScore[size]; var errors = new double[size];
        for (int i = 0; i < _n; i++) for (int b = 0; b < count; b++)
        {
            int at = i * count + b; scores[at] = new(lo[at], hi[at]);
            bool lower = lo[at] <= 0.5; double mean = lower ? lo[at] : hi[at];
            double first = lower ? lo1[at] : hi1[at], second = lower ? lo2[at] : hi2[at];
            double variance = second - 2 * mean * first + mean * mean * q2[b];
            if (!(efficiency[b] > 0)) errors[at] = double.PositiveInfinity;
            else if (variance > 1e-12 * (second + 2 * Math.Abs(mean * first) + mean * mean * q2[b]))
                errors[at] = Math.Sqrt(variance) / Math.Sqrt(efficiency[b]);
            else errors[at] = CenteredError(lower ? _lower : _upper, i * _s, weights, count, b, efficiency[b]);
        }
        return (scores, errors);
    }

    internal double OwnError(int field, double[] weights, int count, int b, CdfScore score, double efficiency)
        => CenteredError(score.Lower <= 0.5 ? _lower : _upper, field * _s, weights, count, b,
            efficiency);

    private double CenteredError(double[] values, int start, double[] weights, int count, int b, double efficiency)
    {
        if (!(efficiency > 0)) return double.PositiveInfinity;
        int pivot = 0;
        for (int d = 1; d < _s; d++) if (weights[d * count + b] > weights[pivot * count + b]) pivot = d;
        double anchor = values[start + pivot], centeredMean = 0;
        for (int d = 0; d < _s; d++) centeredMean = Math.FusedMultiplyAdd(weights[d * count + b], values[start + d] - anchor, centeredMean);
        double scale = 0, squares = 1; bool constant = true;
        for (int d = 0; d < _s; d++)
        {
            if (values[start + d] != values[start]) constant = false;
            double x = Math.Abs(weights[d * count + b] * ((values[start + d] - anchor) - centeredMean));
            if (x == 0) continue;
            if (x > scale) { double r = scale / x; squares = 1 + squares * r * r; scale = x; }
            else { double r = x / scale; squares += r * r; }
        }
        return constant || scale == 0 ? 0 : scale * (Math.Sqrt(squares) / Math.Sqrt(efficiency));
    }
}

internal sealed class CdfMarginScorer
{
    private readonly PointwiseLogDensity _training;
    private readonly int[] _rows, _members;
    private readonly ConformalOptions _options;
    private readonly CdfMomentField _units;
    private readonly CdfScorer _scorer;
    private readonly ConformalBayes.StratumScorer? _proxy;
    private readonly double[] _turnWeights, _turnEfficiency, _carried;
    private readonly double _multiple;
    private int N => _members.Length;
    private int S => _rows.Length;

    internal CdfMarginScorer(PointwiseLogDensity training, int[] rows, double[]? logWeights, int[] members, ConformalOptions options, CdfScorer scorer)
    {
        _training = training; _rows = rows; _members = members; _options = options;
        _scorer = scorer;
        _multiple = ConformalBayes.Multiple(options.ErrorLevel, N + 1);
        _units = new(training.LogCdf!, rows, members);
        if (options.MarginBounds == MarginBounds.Kept)
            _proxy = new(training.Values, rows, logWeights, members, training.ChainLengths, 0, _multiple,
                MarginBounds.Kept, degenerateEss: options.DegenerateEss, flagsOnly: true);
        _carried = options.MarginBounds is null ? [] : new double[N];
        _turnWeights = options.MarginBounds is null ? [] : new double[S * N];
        _turnEfficiency = options.MarginBounds is null ? [] : new double[N];
        if (options.MarginBounds is null) return;
        for (int i = 0; i < N; i++)
        {
            var p = new double[S]; double top = double.NegativeInfinity, weightTop = double.NegativeInfinity, total = 0;
            for (int d = 0; d < S; d++)
            {
                double value = training.Values[rows[d], members[i]];
                if (double.IsNaN(value) || double.IsPositiveInfinity(value))
                    throw new ArgumentException("A carried turn's log density is NaN or positive infinity.", nameof(training));
                top = Math.Max(top, value);
                weightTop = Math.Max(weightTop, value + (logWeights?[d] ?? 0));
            }
            if (double.IsNegativeInfinity(top)) throw new ArgumentException("A carried turn has zero density under every positive-weight draw.", nameof(training));
            for (int d = 0; d < S; d++)
            {
                p[d] = Math.Exp(training.Values[rows[d], members[i]] - top);
                total += _turnWeights[d * N + i] = Math.Exp(training.Values[rows[d], members[i]] + (logWeights?[d] ?? 0) - weightTop);
            }
            for (int d = 0; d < S; d++) _turnWeights[d * N + i] /= total;
            _turnEfficiency[i] = Efficiency(p);
        }
        Array.Copy(Carried(_units, N), _carried, N);
    }

    // Bound temporary moment matrices by fields × Batch, rather than fields × all training turns.
    private double[] Carried(CdfMomentField field, int fields)
    {
        var bounds = new double[fields];
        for (int first = 0; first < N; first += ConformalBayes.Batch)
        {
            int count = Math.Min(ConformalBayes.Batch, N - first);
            var weights = new double[S * count]; var efficiency = new double[count];
            for (int b = 0; b < count; b++) efficiency[b] = _turnEfficiency[first + b];
            for (int d = 0; d < S; d++) for (int b = 0; b < count; b++) weights[d * count + b] = _turnWeights[d * N + first + b];
            var turns = field.Measure(weights, count, efficiency);
            for (int j = 0; j < fields; j++) for (int b = 0; b < count; b++)
                if (!(_proxy?.FreeFlagged[first + b] ?? false)) bounds[j] = Math.Max(bounds[j], _multiple * turns.Errors[j * count + b]);
        }
        return bounds;
    }

    private double Efficiency(double[] p) => _training.ChainLengths is null ? 1 : PredictiveChecks.RelativeEfficiency(p, _training.ChainLengths);

    internal void Read(CandidateLogDensity candidates, int first, int count, double[] weights, CdfScore[] own, int at, CdfMarginReading[] into)
    {
        var efficiency = new double[count]; var ess = new double[count]; var supported = new bool[count];
        for (int b = 0; b < count; b++)
        {
            var p = new double[S]; double top = double.NegativeInfinity, squares = 0;
            for (int d = 0; d < S; d++) top = Math.Max(top, candidates.Values[_rows[d], first + b]);
            supported[b] = !double.IsNegativeInfinity(top);
            for (int d = 0; d < S; d++)
            {
                p[d] = supported[b] ? Math.Exp(candidates.Values[_rows[d], first + b] - top) : 0;
                squares += weights[d * count + b] * weights[d * count + b];
            }
            efficiency[b] = Efficiency(p); ess[b] = supported[b] ? 1 / squares : 0;
        }
        var held = _units.Measure(weights, count, efficiency);
        var field = new CdfMomentField(candidates.LogCdf!, _rows, Enumerable.Range(first, count).ToArray());
        double[] candidateCarried = _options.MarginBounds is null ? [] : Carried(field, count);
        for (int b = 0; b < count; b++)
        {
            var scores = new CdfScore[N + 1]; var errors = new double[N + 1];
            double[]? carried = _options.MarginBounds is null ? null : new double[N + 1];
            for (int i = 0; i < N; i++)
            {
                scores[i] = held.Scores[i * count + b]; errors[i] = supported[b] ? _multiple * held.Errors[i * count + b] : double.PositiveInfinity;
                if (carried is not null) carried[i] = Math.Max(errors[i], _carried[i]);
            }
            scores[N] = own[at + b];
            errors[N] = supported[b] ? _multiple * field.OwnError(b, weights, count, b, own[at + b], efficiency[b]) : double.PositiveInfinity;
            if (carried is not null)
            {
                carried[N] = Math.Max(errors[N], candidateCarried[b]);
            }
            bool ownFlag = false; bool[]? flags = null;
            if (_proxy is not null)
            {
                double[] w = Enumerable.Range(0, S).Select(d => weights[d * count + b]).ToArray();
                (bool flagged, bool[] units) = _proxy.Flags(w, ess[b], efficiency[b]); ownFlag = flagged;
                flags = [.. units, flagged];
            }
            into[at + b] = CdfMarginCounts.Read(scores, errors, carried, ownFlag, flags, _options.MarginBounds, _options.InflationBounds);
        }
    }

    internal ConformalCdfMargin Finish(PointwiseLogDensity training, CandidateLogDensity candidates, int t, CdfMarginReading[] readings)
    {
        var probes = new List<CdfMarginReading>();
        if (_options.MarginRefinement is { } refinement)
        {
            bool Same(CdfMarginReading a, CdfMarginReading b) => a.CandidateFlagged == b.CandidateFlagged && a.Flagged == b.Flagged
                && (a.Fragile ?? []).SequenceEqual(b.Fragile ?? []) && (a.DegenerateAccepted ?? []).SequenceEqual(b.DegenerateAccepted ?? [])
                && (a.PitLowerCounts ?? []).SequenceEqual(b.PitLowerCounts ?? []) && (a.PitUpperCounts ?? []).SequenceEqual(b.PitUpperCounts ?? []);
            var cells = new Stack<(double Lo, double Hi, CdfMarginReading Left, CdfMarginReading Right)>();
            double[] grid = candidates.Grids[t];
            for (int g = 1; g < grid.Length; g++) cells.Push((grid[g - 1], grid[g], readings[g - 1], readings[g]));
            while (cells.TryPop(out var cell))
            {
                if (Same(cell.Left, cell.Right) || cell.Hi - cell.Lo <= refinement.Tolerance) continue;
                double y = cell.Lo < 0 && cell.Hi > 0 ? (cell.Lo + cell.Hi) * .5 : cell.Lo + (cell.Hi - cell.Lo) * .5;
                if (!(y > cell.Lo && y < cell.Hi)) continue;
                CandidateLogDensity probe = refinement.At(t, y);
                ConformalBayes.Check(training, probe);
                if (probe.Cases != 1 || probe.Grids[0].Length != 1 || probe.Grids[0][0] != y || probe.StratumOf[0] != candidates.StratumOf[t] || probe.LogCdf is null)
                    throw new ArgumentException("A CDF margin probe must carry its log CDF at one response in the original stratum.", nameof(refinement));
                var own = new CdfScore[1]; var reading = new CdfMarginReading[1];
                _scorer.Score(probe.Values, probe.LogCdf, 0, 1, new int[1], new int[1], own, 0);
                Read(probe, 0, 1, _scorer.CopyWeights(1), own, 0, reading);
                probes.Add(reading[0]);
                cells.Push((cell.Lo, y, cell.Left, reading[0])); cells.Push((y, cell.Hi, reading[0], cell.Right));
            }
        }
        return new(N + 1, _options, _multiple, _carried, readings, probes);
    }
}
