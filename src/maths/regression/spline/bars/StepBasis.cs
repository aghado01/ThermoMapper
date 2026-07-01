using System;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Piecewise-constant carrier realization (Green 1995, §4): the interior knots are changepoints splitting [0,1]
/// into k+1 segments, each with its own level. The design is the segment indicator matrix, so a fit is a step
/// function — the home of the discontinuity machinery, for integer/level observables over a swept parameter
/// (e.g. b₁(T)). Shares the carrier, moves, marginal, and ensemble with <see cref="SplineBasis"/>; only the
/// design differs. The continuous-peak readout does not apply here (the ensemble skips it for step bases).
/// </summary>
public sealed class StepBasis : IBasis
{
    /// <summary>k+1 segment levels.</summary>
    public int Dimension(int knotCount) => knotCount + 1;

    public double[,] Design(KnotConfig config, double[] x)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(x);
        double[] knots = config.InteriorKnots;
        int nu = knots.Length + 1;
        var z = new double[x.Length, nu];
        for (int i = 0; i < x.Length; i++)
            z[i, Segment(knots, x[i])] = 1.0;
        return z;
    }

    public double Evaluate(KnotConfig config, double[] coef, double x)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(coef);
        return coef[Segment(config.InteriorKnots, x)];
    }

    // Segment index = number of changepoints at or below x (knots are sorted ascending).
    private static int Segment(double[] knots, double x)
    {
        int seg = 0;
        while (seg < knots.Length && x >= knots[seg]) seg++;
        return seg;
    }
}
