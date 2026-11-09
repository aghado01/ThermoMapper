using System.Collections.Generic;

namespace Maths.Regression.Conformal;

/// <summary>Closed asymmetric CQR ends and the two calibrated tail levels (<c>cqrAsym_coverage</c>).</summary>
public sealed record AsymmetricConformalInterval(int Case, int Stratum, int Units, int LowerLevel, int UpperLevel,
                                                double LowerThreshold, double UpperThreshold, double Lower, double Upper)
{
    public bool WholeSpace => LowerLevel == 0 && UpperLevel == 0;
    public bool Empty => Lower > Upper;
    public double Width => Empty ? 0.0 : Upper - Lower;
    public bool Contains(double y) => y >= Lower && y <= Upper;
    public double LowerMiscoverageBound => (double)LowerLevel / Units;
    public double UpperMiscoverageBound => (double)UpperLevel / Units;
    public double Coverage => 1.0 - (LowerLevel + UpperLevel) / (double)Units;
}

/// <summary>Asymmetric CQR intervals with the two calibration samples ranked at α/2.</summary>
public sealed record AsymmetricConformalIntervals(IReadOnlyList<AsymmetricConformalInterval> Intervals, double Alpha,
                                                 SplitCalibration LowerCalibration, SplitCalibration UpperCalibration);
