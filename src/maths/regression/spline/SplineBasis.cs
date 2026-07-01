using System;

namespace Maths.Regression.Spline;

/// <summary>
/// Realizes a <see cref="KnotConfig"/> as a B-spline design matrix — the Transform from the (k,ξ) Field to the
/// linear model the observation likelihood sees. Fixed degree (the model hyperparameter); each row is the
/// Cox–de Boor basis (<see cref="BSpline.EvaluateBasis"/>) evaluated at one design point over the clamped knot
/// vector. The basis dimension is ν = k + degree + 1.
/// </summary>
public sealed class SplineBasis : IBasis
{
    private readonly int _degree;

    public SplineBasis(int degree = 3)
    {
        if (degree < 1) throw new ArgumentOutOfRangeException(nameof(degree), "Degree must be at least 1.");
        _degree = degree;
    }

    /// <summary>Number of basis functions (design columns) for a config with <paramref name="knotCount"/> interior knots.</summary>
    public int Dimension(int knotCount) => knotCount + _degree + 1;

    /// <summary>The spline degree (3 = cubic).</summary>
    public int Degree => _degree;

    /// <summary>Evaluate the spline <c>Σ coef_j B_j(x)</c> at a single point <paramref name="x"/> ∈ [0,1].</summary>
    public double Evaluate(KnotConfig config, double[] coef, double x)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(coef);
        double[] knots = BSpline.MakeClampedKnots(config.InteriorKnots, _degree);
        double[] row = BSpline.EvaluateBasis(x, knots, _degree);
        double f = 0.0;
        for (int j = 0; j < row.Length; j++) f += row[j] * coef[j];
        return f;
    }

    /// <summary>
    /// Build the <c>m × ν</c> design matrix for <paramref name="config"/> at design points
    /// <paramref name="x"/> (each in [0,1]).
    /// </summary>
    public double[,] Design(KnotConfig config, double[] x)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(x);

        double[] knots = BSpline.MakeClampedKnots(config.InteriorKnots, _degree);
        int nu = knots.Length - _degree - 1;
        var z = new double[x.Length, nu];
        for (int i = 0; i < x.Length; i++)
        {
            double[] row = BSpline.EvaluateBasis(x[i], knots, _degree);
            for (int j = 0; j < nu; j++)
                z[i, j] = row[j];
        }
        return z;
    }
}
