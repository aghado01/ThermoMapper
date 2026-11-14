using System;
using System.Collections.Generic;

namespace Maths.Regression.Conformal;

/// <summary>The split distribution's ordered calibration PITs, with the draw weights and both tails of each score retained.</summary>
public sealed class SplitPitCalibration
{
    private readonly SplitDraws _draws;
    private readonly int[] _rows;
    private readonly double[] _weights;
    private readonly Dictionary<int, CdfScore[]> _sorted = [];
    private readonly Dictionary<int, int[]> _order = [];   // each stratum's calibration units, in their scores' sorted order
    private readonly int _units;

    internal SplitPitCalibration(CandidateLogDensity calibration)
    {
        SplitConformal.CheckCalibration(calibration);
        _draws = new SplitDraws(calibration);
        (_rows, double[]? logWeights) = PredictiveChecks.Rows(calibration.Weights, calibration.Draws);
        _weights = CdfScores.Weights(_rows, logWeights);
        var scores = new CdfScore[calibration.Cases];
        _units = scores.Length;
        Score(RequireCdf(calibration), 0, scores);
        var groups = new Dictionary<int, List<int>>();
        for (int i = 0; i < scores.Length; i++)
        {
            int h = calibration.StratumOf[i];
            if (!groups.TryGetValue(h, out var group)) groups.Add(h, group = []);
            group.Add(i);
        }
        foreach ((int h, List<int> group) in groups)
        {
            int[] order = [.. group];
            CdfScore[] sorted = Array.ConvertAll(order, i => scores[i]);
            Array.Sort(sorted, order);
            _sorted.Add(h, sorted);
            _order.Add(h, order);
        }
    }

    private static double[,] RequireCdf(CandidateLogDensity seam) =>
        seam.LogCdf ?? throw new ArgumentException("The predictive distribution requires a log CDF.", nameof(seam));

    private double[,] CheckCandidates(CandidateLogDensity candidates)
    {
        _draws.Check(candidates);
        return RequireCdf(candidates);
    }

    private void Score(double[,] logCdf, int first, Span<CdfScore> scores)
        => CdfScores.Fixed(logCdf, _rows, _weights, first, scores.Length, scores);

    private CdfScore[] Sorted(int stratum) => _sorted.TryGetValue(stratum, out var scores) ? scores : [];

    /// <summary>
    /// Read the split distribution on the supplied grids using the same calibration as the grid-free readouts. With
    /// <paramref name="weights"/> (<see cref="ShiftWeights"/>: one per calibration unit in the calibration's order, or per case, and one per
    /// candidate column) each case whose weights are not all one value reads its <see cref="ConformalDistribution.Weighted"/> ranks
    /// (slice 14's follow-up 1, <c>pit_weighted</c>).
    /// </summary>
    public ConformalDistributions OnGrid(CandidateLogDensity candidates, ShiftWeights? weights = null)
    {
        double[,] cdf = CheckCandidates(candidates);
        if (weights is not null)
        {
            if (weights.UnitCount != _units)
                throw new ArgumentException($"One weight per calibration unit is required: {weights.UnitCount} weights for {_units} units.", nameof(weights));
            if (weights.UnitsByCase is { } byCase && byCase.Count != candidates.Cases)
                throw new ArgumentException($"One set of unit weights per case is required: {byCase.Count} sets for {candidates.Cases} cases.", nameof(weights));
            if (weights.Candidates.Count != candidates.Values.GetLength(1))
                throw new ArgumentException($"One weight per candidate is required: {weights.Candidates.Count} weights for {candidates.Values.GetLength(1)} candidates.", nameof(weights));
        }
        var distributions = new ConformalDistribution[candidates.Cases];
        for (int t = 0; t < candidates.Cases; t++)
        {
            int h = candidates.StratumOf[t];
            CdfScore[] units = Sorted(h);
            double[] grid = candidates.Grids[t];
            var own = new CdfScore[grid.Length];
            Score(cdf, candidates.Offset(t), own);
            var below = new int[grid.Length];
            var tied = new int[grid.Length];
            for (int g = 0; g < grid.Length; g++) (below[g], tied[g]) = CdfScores.Rank(units, own[g]);
            distributions[t] = new ConformalDistribution(t, h, units.Length + 1, grid, below, tied, ConformalForm.Split,
                                                         CdfScores.Fall(below, tied), CdfScores.ScoreFall(own))
            {
                Weighted = weights is null ? null : Weigh(weights, t, h, candidates.Offset(t), below, tied),
            };
        }
        return new ConformalDistributions(distributions, ConformalForm.Split);
    }

    // Case t's weighted ranks: the stratum's unit weights summed along its sorted order, so the units below are the first B and those at or
    // below the first B + T; null where its weights are all one value, CB-1's distribution.
    private WeightedRanks? Weigh(ShiftWeights weights, int t, int h, int offset, int[] below, int[] tied)
    {
        IReadOnlyList<double> units = weights.UnitsFor(t);
        int[] order = _order.TryGetValue(h, out int[]? found) ? found : [];
        var running = new double[order.Length + 1];
        double common = order.Length > 0 ? units[order[0]] : weights.Candidates[offset];
        bool equal = true;
        for (int j = 0; j < order.Length; j++)
        {
            double w = units[order[j]];
            running[j + 1] = running[j] + w;
            if (w != common) equal = false;
        }
        for (int g = 0; g < below.Length; g++) if (weights.Candidates[offset + g] != common) equal = false;
        if (equal) return null;
        var weightBelow = new double[below.Length];
        var weightAtOrBelow = new double[below.Length];
        for (int g = 0; g < below.Length; g++)
        {
            weightBelow[g] = running[below[g]];
            weightAtOrBelow[g] = running[below[g] + tied[g]];
        }
        return CdfScores.Ranks(weightBelow, weightAtOrBelow, running[^1], weights.Candidates, offset);
    }

    /// <summary>
    /// The distribution at a fixed test input. At every y, including the initial probe at zero, <paramref name="candidateAt"/> supplies
    /// the same fit's seam with [y] as the selected case's grid. An optional inverse CDF of that case's fitted predictive supplies
    /// the initial root estimate. Crossing brackets use absolute response tolerance and the seam's own two-tail score throughout.
    /// With <paramref name="weights"/> (one per calibration unit, and one for the case: covariate shift's w(x), constant in y) the readouts
    /// are the weighted distribution's (slice 14's follow-up 1, <see cref="ConformalDistribution.Weighted"/>), equal weights giving the
    /// counts' bit for bit. Under response shift the case's weight moves with y; the grid form serves it.
    /// </summary>
    public SplitPredictiveDistribution AtCase(Func<double, CandidateLogDensity> candidateAt, int testCase = 0,
                                              Func<double, double>? predictiveQuantile = null, double tolerance = 1e-8,
                                              ShiftWeights? weights = null)
    {
        ArgumentNullException.ThrowIfNull(candidateAt);
        ArgumentOutOfRangeException.ThrowIfNegative(testCase);
        ConformalBisection.CheckTolerance(tolerance);
        if (weights is not null)
        {
            if (weights.UnitCount != _units)
                throw new ArgumentException($"One weight per calibration unit is required: {weights.UnitCount} weights for {_units} units.", nameof(weights));
            if (weights.Candidates.Count != 1 || weights.UnitsByCase is { Count: not 1 })
                throw new ArgumentException("One weight for the case is required, constant in y (covariate shift); the grid form serves weights that move with y.",
                                            nameof(weights));
        }
        int? stratum = null;
        CdfScore Read(double y)
        {
            CandidateLogDensity at = candidateAt(y);
            double[,] cdf = CheckCandidates(at);
            if (testCase >= at.Cases || at.Grids[testCase].Length != 1 || at.Grids[testCase][0] != y)
                throw new ArgumentException("The selected case must carry the requested y as its only candidate.", nameof(candidateAt));
            int h = at.StratumOf[testCase];
            stratum ??= h;
            if (h != stratum) throw new ArgumentException("The selected case changed stratum.", nameof(candidateAt));
            Span<CdfScore> value = stackalloc CdfScore[1];
            Score(cdf, at.Offset(testCase), value);
            return value[0];
        }
        CdfScore atZero = Read(0.0);
        int h = stratum!.Value;
        double[]? running = null;
        double own = 0.0;
        if (weights is not null)
        {
            // The stratum's weights summed along its sorted order, as the grid form sums them; none where all equal the case's.
            IReadOnlyList<double> units = weights.UnitsFor(0);
            int[] order = _order.TryGetValue(h, out int[]? found) ? found : [];
            own = weights.Candidates[0];
            running = new double[order.Length + 1];
            bool equal = true;
            for (int j = 0; j < order.Length; j++)
            {
                running[j + 1] = running[j] + units[order[j]];
                if (units[order[j]] != own) equal = false;
            }
            if (equal) running = null;
        }
        return new SplitPredictiveDistribution(h, Sorted(h), y => y == 0.0 ? atZero : Read(y), predictiveQuantile, tolerance, running, own);
    }
}

/// <summary>Grid-free readouts of the calibrated predictive CDF at one case, using its stratum's PIT order statistics.</summary>
public sealed class SplitPredictiveDistribution
{
    private readonly CdfScore[] _sorted;
    private readonly Func<double, CdfScore> _score;
    private readonly Func<double, double>? _quantile;
    private readonly double _tolerance;

    private readonly double[]? _running;   // the weight of the first j sorted PITs, where weighted
    private readonly double _own;

    internal SplitPredictiveDistribution(int stratum, CdfScore[] sorted, Func<double, CdfScore> score,
                                         Func<double, double>? quantile, double tolerance, double[]? running = null, double own = 0.0)
        => (Stratum, _sorted, _score, _quantile, _tolerance, _running, _own) = (stratum, sorted, score, quantile, tolerance, running, own);

    public int Stratum { get; }
    public int Units => _sorted.Length + 1;

    /// <summary>Whether the readouts are the weighted distribution's: the case's shift weights not all one value.</summary>
    public bool Weighted => _running is not null;

    // Q⁻ and Q⁺ with B PITs below and B + T at or below, as the grid form reads them.
    private (double Lower, double Upper) Shares(int below, int atOrBelow)
    {
        double total = _running![^1] + _own;
        return total > 0.0 ? (_running[below] / total, (_running[atOrBelow] + _own) / total) : (0.0, 1.0);
    }

    /// <summary>Q(y, τ) at any finite response. At the held-out response and uniform τ this is the randomized PIT.</summary>
    public double At(double y, double tau)
    {
        if (!double.IsFinite(y)) throw new ArgumentOutOfRangeException(nameof(y));
        if (!(tau >= 0.0 && tau <= 1.0)) throw new ArgumentOutOfRangeException(nameof(tau));
        (int below, int tied) = CdfScores.Rank(_sorted, _score(y));
        if (_running is null) return (below + tau * (tied + 1)) / Units;
        (double lower, double upper) = Shares(below, below + tied);
        return lower + tau * (upper - lower);
    }

    /// <summary>The infimum of responses where Q(y, ½) reaches p, including the two unbounded outer ranks.</summary>
    public double Quantile(double p)
    {
        if (!(p > 0 && p < 1)) throw new ArgumentOutOfRangeException(nameof(p));
        // The plateau just above the j-th PIT has Q̄ = (j + ½)/m, or the shares' mean where weighted. Compare as the grid does at
        // floating-point rank boundaries.
        int lo = 0, hi = Units;
        while (lo < hi)
        {
            int mid = lo + ((hi - lo) >> 1);
            if (Mean(mid, mid) < p) lo = mid + 1;
            else hi = mid;
        }
        if (lo == 0) return double.NegativeInfinity;
        if (lo == Units) return double.PositiveInfinity;
        CdfScore pit = _sorted[lo - 1];
        (int below, int tied) = CdfScores.Rank(_sorted, pit);
        bool strict = Mean(below, below + tied) < p;
        return Crossing(pit, strict).Upper;
    }

    // Q̄ with B PITs below and B + T at or below: (B + (T + 1)/2)/m, or (Q⁻ + Q⁺)/2 where weighted.
    private double Mean(int below, int atOrBelow)
    {
        if (_running is null) return (below + 0.5 * (atOrBelow - below + 1)) / Units;
        (double lower, double upper) = Shares(below, atOrBelow);
        return 0.5 * (lower + upper);
    }

    /// <summary>Closed equi-tailed ends. With k = ⌊αm/2⌋, the fitted CDF must reach PIT k and stay at or below PIT m − k
    /// (one-based). Round both ends outward by at most the response tolerance, including when the two PITs tie.</summary>
    public GridInterval Interval(double alpha)
    {
        SplitConformal.CheckAlpha(alpha);
        if (_running is not null) return WeightedInterval(alpha);
        int k = (int)Math.Floor(0.5 * alpha * Units);
        if (k == 0) return new(double.NegativeInfinity, double.PositiveInfinity);
        return new(Crossing(_sorted[k - 1], strict: false).Lower, Crossing(_sorted[Units - k - 1], strict: true).Upper);
    }

    // The weighted set's ends: y is held where Q⁺ > α/2, B + T at least the least j with (R_j + w_c)/W > α/2, and where 1 − Q⁻ > α/2, B below
    // the least j with 1 − R_j/W not above α/2; each test as the grid form reads it.
    private GridInterval WeightedInterval(double alpha)
    {
        double total = _running![^1] + _own;
        if (!(total > 0.0)) return new(double.NegativeInfinity, double.PositiveInfinity);
        int n = _sorted.Length, low = 0, high = 0;
        while (low < n && !((_running[low] + _own) / total > 0.5 * alpha)) low++;
        while (high <= n && 1.0 - _running[high] / total > 0.5 * alpha) high++;
        double lower = low == 0 ? double.NegativeInfinity : Crossing(_sorted[low - 1], strict: false).Lower;
        double upper = high > n ? double.PositiveInfinity : high == 0 ? double.NegativeInfinity : Crossing(_sorted[high - 1], strict: true).Upper;
        return new(lower, upper);
    }

    private GridInterval Crossing(CdfScore pit, bool strict)
    {
        // Find the first y with F >= pit (or F > pit). Strictness selects the far side of any tied plateau.
        if (!strict && pit.Lower == 0.0) return new(double.NegativeInfinity, double.NegativeInfinity);
        if (strict && pit.Upper == 0.0) return new(double.PositiveInfinity, double.PositiveInfinity);
        bool Above(double y) => strict ? _score(y).CompareTo(pit) > 0 : _score(y).CompareTo(pit) >= 0;

        double p = pit.Lower <= 0.5 ? pit.Lower : 1.0 - pit.Upper;
        bool useInverse = _quantile is not null && p > 0 && p < 1;
        double start = useInverse ? _quantile!(p) : 0.0;
        if (!double.IsFinite(start)) throw new InvalidOperationException("The fitted inverse returned a nonfinite interior quantile.");
        bool startAbove = Above(start);
        double radius = useInverse ? 0.5 * _tolerance : 1.0;
        double near = start;
        while (true)
        {
            double far = startAbove ? Math.Max(-double.MaxValue, start - radius) : Math.Min(double.MaxValue, start + radius);
            if (radius == double.MaxValue) far = startAbove ? -double.MaxValue : double.MaxValue;
            if (Above(far) != startAbove)
            {
                return ConformalBisection.Narrow(Math.Min(near, far), Math.Max(near, far), false, Above, _tolerance);
            }
            if (far == (startAbove ? -double.MaxValue : double.MaxValue))
            {
                double end = startAbove ? double.NegativeInfinity : double.PositiveInfinity;
                return new(end, end);
            }
            near = far;
            radius = Math.Min(double.MaxValue, Math.Max(1.0, 2.0 * radius));
        }
    }
}
