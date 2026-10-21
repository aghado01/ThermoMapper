using System;
using Maths.LinAlg;
using Maths.Regression.GPrior;

namespace Maths.Regression.Spline;

/// <summary>
/// Builds the collapsed g-prior block from a banded spline design: the Gram matrices go straight into band storage
/// (<see cref="BandedDesign"/>), are factored by a root-free <see cref="BandCholesky"/> (LDLᵀ, so the block can
/// draw), and the explicit residual comes from <see cref="BandedDesign.Apply"/>. The block's closed forms are the
/// same for any design layer; this is the banded one BARS uses.
/// </summary>
public static class BandedGPrior
{
    /// <summary>The unweighted block: G = ZᵀZ factored, β_OLS, <c>‖y − Zβ_OLS‖²</c> and <c>Q = ‖Zβ_OLS‖²</c>.</summary>
    public static CollapsedGPrior Fit(BandedDesign design, double[] y)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        int m = design.Rows, nu = design.Dimension;
        if (y.Length != m) throw new ArgumentException("Response length must match the design row count.", nameof(y));

        var band = new double[design.Bandwidth + 1, nu];
        var b = new double[nu];
        design.Accumulate(null, band, y, b);
        var gram = new BandCholesky(nu, design.Bandwidth, BandFactorization.Ldlt);
        gram.DecomposeBanded(band);
        double[] beta = gram.Solve(b);   // β_OLS = G⁻¹Zᵀy; no inverse is formed

        var f = new double[m];
        design.Apply(beta, f);
        double rss = 0.0, fitted = 0.0;
        for (int i = 0; i < m; i++)
        {
            double r = y[i] - f[i];
            rss += r * r;
            fitted += f[i] * f[i];
        }
        return CollapsedGPrior.FromLeastSquares(gram, beta, rss, fitted, m);
    }

    /// <summary>
    /// The weighted block at prior precision <paramref name="tau"/>: G = ZᵀZ and <c>P_a = Zᵀ(W + τI)Z</c> factored,
    /// the posterior mean <c>β̃ = P_a⁻¹ZᵀWy</c>, and <c>S = ‖y − Zβ̃‖²_W + τ‖Zβ̃‖²</c>, a sum of non-negative terms.
    /// </summary>
    public static CollapsedGPrior Fit(BandedDesign design, double[] y, double[] weights, double tau)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(weights);
        int m = design.Rows, nu = design.Dimension;
        if (y.Length != m) throw new ArgumentException("Response length must match the design row count.", nameof(y));
        if (weights.Length != m) throw new ArgumentException("Weights length must match the design row count.", nameof(weights));

        var gBand = new double[design.Bandwidth + 1, nu];
        design.Accumulate(null, gBand, null, null);
        var gram = new BandCholesky(nu, design.Bandwidth, BandFactorization.Ldlt);
        gram.DecomposeBanded(gBand);

        var wPlus = new double[m];
        double sumLogW = 0.0;
        for (int i = 0; i < m; i++)
        {
            wPlus[i] = weights[i] + tau;
            sumLogW += Math.Log(weights[i]);
        }
        var pBand = new double[design.Bandwidth + 1, nu];
        design.Accumulate(wPlus, pBand, null, null);
        var h = new double[nu];
        design.AccumulateRhs(weights, y, h);
        var posterior = new BandCholesky(nu, design.Bandwidth, BandFactorization.Ldlt);
        posterior.DecomposeBanded(pBand);
        double[] beta = posterior.Solve(h);

        var f = new double[m];
        design.Apply(beta, f);
        double weightedRss = 0.0, fitted = 0.0;
        for (int i = 0; i < m; i++)
        {
            double r = y[i] - f[i];
            weightedRss += weights[i] * r * r;
            fitted += f[i] * f[i];
        }
        return CollapsedGPrior.FromWeighted(gram, posterior, beta, weightedRss + tau * fitted, sumLogW, m, tau);
    }
}
