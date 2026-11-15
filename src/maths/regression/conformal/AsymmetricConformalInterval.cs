using System;
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

    /// <summary>The interval's own miscoverage bound (<see cref="TrainingConditional"/>): a miss falls below or above, two disjoint events,
    /// each tail's law Beta(k, m − k) at its own k, so each tail's bound at confidence (1 + β)/2 and their sum hold together with chance at
    /// least <paramref name="beta"/> over the calibration units; null under the weighted rule. The union holds whatever ties the two
    /// tails' scores, and so is loose: at m 100 and β 0.9 a Gaussian fixture passed it 0.85% of the time against 10%. The symmetric CQR
    /// interval, one score, reads the law exactly (<see cref="ConformalInterval.TrainingConditionalMiscoverage"/>).</summary>
    public double? TrainingConditionalMiscoverage(double beta)
    {
        if (WeightedAlpha is not null) return null;
        if (!(beta > 0.0 && beta < 1.0)) throw new ArgumentOutOfRangeException(nameof(beta), "β lies in (0, 1).");
        double each = (1.0 + beta) / 2.0;
        return Math.Min(1.0, TrainingConditional.Miscoverage(Units, LowerLevel, each) + TrainingConditional.Miscoverage(Units, UpperLevel, each));
    }

    /// <summary>The weighted rule's α where each tail is the weighted rule's at α/2 (slice 14, <see cref="ShiftWeights"/>), the thresholds
    /// the case's weighted ones and the levels CB-1's at the same m; null for CB-1's intervals.</summary>
    public double? WeightedAlpha { get; init; }
}

/// <summary>Asymmetric CQR intervals with the two calibration samples ranked at α/2.</summary>
public sealed record AsymmetricConformalIntervals(IReadOnlyList<AsymmetricConformalInterval> Intervals, double Alpha,
                                                 SplitCalibration LowerCalibration, SplitCalibration UpperCalibration);
