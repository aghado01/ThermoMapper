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

    internal SplitPitCalibration(CandidateLogDensity calibration)
    {
        SplitConformal.CheckCalibration(calibration);
        _draws = new SplitDraws(calibration);
        (_rows, double[]? logWeights) = PredictiveChecks.Rows(calibration.Weights, calibration.Draws);
        _weights = CdfScores.Weights(_rows, logWeights);
        var scores = new CdfScore[calibration.Cases];
        Score(RequireCdf(calibration), 0, scores);
        var groups = new Dictionary<int, List<CdfScore>>();
        for (int i = 0; i < scores.Length; i++)
        {
            int h = calibration.StratumOf[i];
            if (!groups.TryGetValue(h, out var group)) groups.Add(h, group = []);
            group.Add(scores[i]);
        }
        foreach ((int h, List<CdfScore> group) in groups)
        {
            CdfScore[] sorted = [.. group];
            Array.Sort(sorted);
            _sorted.Add(h, sorted);
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

    /// <summary>Read the split distribution on the supplied grids using the same calibration as the grid-free readouts.</summary>
    public ConformalDistributions OnGrid(CandidateLogDensity candidates)
    {
        double[,] cdf = CheckCandidates(candidates);
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
                                                         CdfScores.Fall(below, tied), CdfScores.ScoreFall(own));
        }
        return new ConformalDistributions(distributions, ConformalForm.Split);
    }

    /// <summary>
    /// The distribution at a fixed test input. At every y, including the initial probe at zero, <paramref name="candidateAt"/> supplies
    /// the same fit's seam with [y] as the selected case's grid. An optional inverse CDF of that case's fitted predictive supplies
    /// the initial root estimate. Crossing brackets use absolute response tolerance and the seam's own two-tail score throughout.
    /// </summary>
    public SplitPredictiveDistribution AtCase(Func<double, CandidateLogDensity> candidateAt, int testCase = 0,
                                              Func<double, double>? predictiveQuantile = null, double tolerance = 1e-8)
    {
        ArgumentNullException.ThrowIfNull(candidateAt);
        ArgumentOutOfRangeException.ThrowIfNegative(testCase);
        ConformalBisection.CheckTolerance(tolerance);
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
        return new SplitPredictiveDistribution(h, Sorted(h), y => y == 0.0 ? atZero : Read(y), predictiveQuantile, tolerance);
    }
}

/// <summary>Grid-free readouts of the calibrated predictive CDF at one case, using its stratum's PIT order statistics.</summary>
public sealed class SplitPredictiveDistribution
{
    private readonly CdfScore[] _sorted;
    private readonly Func<double, CdfScore> _score;
    private readonly Func<double, double>? _quantile;
    private readonly double _tolerance;

    internal SplitPredictiveDistribution(int stratum, CdfScore[] sorted, Func<double, CdfScore> score,
                                         Func<double, double>? quantile, double tolerance)
        => (Stratum, _sorted, _score, _quantile, _tolerance) = (stratum, sorted, score, quantile, tolerance);

    public int Stratum { get; }
    public int Units => _sorted.Length + 1;

    /// <summary>Q(y, τ) at any finite response. At the held-out response and uniform τ this is the randomized PIT.</summary>
    public double At(double y, double tau)
    {
        if (!double.IsFinite(y)) throw new ArgumentOutOfRangeException(nameof(y));
        if (!(tau >= 0.0 && tau <= 1.0)) throw new ArgumentOutOfRangeException(nameof(tau));
        (int below, int tied) = CdfScores.Rank(_sorted, _score(y));
        return (below + tau * (tied + 1)) / Units;
    }

    /// <summary>The infimum of responses where Q(y, ½) reaches p, including the two unbounded outer ranks.</summary>
    public double Quantile(double p)
    {
        if (!(p > 0 && p < 1)) throw new ArgumentOutOfRangeException(nameof(p));
        // The plateau just above the j-th PIT has Q̄ = (j + ½)/m. Compare as the grid does at floating-point rank boundaries.
        int lo = 0, hi = Units;
        while (lo < hi)
        {
            int mid = lo + ((hi - lo) >> 1);
            if ((mid + 0.5) / Units < p) lo = mid + 1;
            else hi = mid;
        }
        if (lo == 0) return double.NegativeInfinity;
        if (lo == Units) return double.PositiveInfinity;
        CdfScore pit = _sorted[lo - 1];
        (int below, int tied) = CdfScores.Rank(_sorted, pit);
        bool strict = (below + 0.5 * (tied + 1)) / Units < p;
        return Crossing(pit, strict).Upper;
    }

    /// <summary>Closed equi-tailed ends. With k = ⌊αm/2⌋, the fitted CDF must reach PIT k and stay at or below PIT m − k
    /// (one-based). Round both ends outward by at most the response tolerance, including when the two PITs tie.</summary>
    public GridInterval Interval(double alpha)
    {
        SplitConformal.CheckAlpha(alpha);
        int k = (int)Math.Floor(0.5 * alpha * Units);
        if (k == 0) return new(double.NegativeInfinity, double.PositiveInfinity);
        return new(Crossing(_sorted[k - 1], strict: false).Lower, Crossing(_sorted[Units - k - 1], strict: true).Upper);
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
