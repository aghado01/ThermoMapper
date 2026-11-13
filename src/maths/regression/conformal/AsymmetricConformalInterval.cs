using System.Collections.Generic;

namespace Maths.Regression.Conformal;

/// <summary>Closed asymmetric CQR ends and the two calibrated tail levels (<c>cqrAsym_coverage</c>).</summary>
public sealed record AsymmetricConformalInterval(int Case, int Stratum, int Units, int LowerLevel, int UpperLevel,
                                                double LowerThreshold, double UpperThreshold, double Lower, double Upper)
{
    public bool WholeSpace => WeightedAlpha is null ? LowerLevel == 0 && UpperLevel == 0
        : double.IsNegativeInfinity(LowerThreshold) && double.IsNegativeInfinity(UpperThreshold);
    public bool Empty => Lower > Upper;
    public double Width => Empty ? 0.0 : Upper - Lower;
    public bool Contains(double y) => y >= Lower && y <= Upper;
    public double LowerMiscoverageBound => WeightedAlpha is double a ? a / 2 : (double)LowerLevel / Units;
    public double UpperMiscoverageBound => WeightedAlpha is double a ? a / 2 : (double)UpperLevel / Units;
    public double Coverage => WeightedAlpha is double a ? 1.0 - a : 1.0 - (LowerLevel + UpperLevel) / (double)Units;

    /// <summary>The weighted rule's α where each tail is the weighted rule's at α/2 (slice 14, <see cref="ShiftWeights"/>), the thresholds
    /// the case's weighted ones and the levels CB-1's at the same m; null for CB-1's intervals.</summary>
    public double? WeightedAlpha { get; init; }
}

/// <summary>Asymmetric CQR intervals with the two calibration samples ranked at α/2.</summary>
public sealed record AsymmetricConformalIntervals(IReadOnlyList<AsymmetricConformalInterval> Intervals, double Alpha,
                                                 SplitCalibration LowerCalibration, SplitCalibration UpperCalibration);
