namespace Maths.Regression.Spline;

/// <summary>
/// Realizes a <see cref="KnotConfig"/> as a linear model — the Transform from the (k,ξ) Field to the design the
/// observation likelihood sees. The carrier (knots/changepoints), the moves, the marginal, and the ensemble are
/// all basis-agnostic; only the design construction differs. <see cref="SplineBasis"/> gives continuous
/// piecewise polynomials (knots); <see cref="Bars.StepBasis"/> gives piecewise constants (changepoints) for
/// step-function inference (e.g. integer/level observables like b₁(T)).
/// </summary>
public interface IBasis
{
    /// <summary>Number of basis functions (design columns) for a config with <paramref name="knotCount"/> interior knots.</summary>
    int Dimension(int knotCount);

    /// <summary>Build the <c>m × ν</c> design matrix for <paramref name="config"/> at the points <paramref name="x"/> (each in [0,1]).</summary>
    double[,] Design(KnotConfig config, double[] x);

    /// <summary>Evaluate the realized function <c>Σ coef_j B_j(x)</c> at a single point.</summary>
    double Evaluate(KnotConfig config, double[] coef, double x);
}
