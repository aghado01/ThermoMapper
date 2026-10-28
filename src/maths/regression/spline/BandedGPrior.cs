using System;
using Maths.LinAlg;
using Maths.Regression.GPrior;

namespace Maths.Regression.Spline;

/// <summary>
/// Builds the collapsed g-prior block from a banded spline design: the Gram matrices go straight into band storage
/// (<see cref="BandedDesign"/>), are factored by a root-free <see cref="BandCholesky"/> (LDLᵀ, so the block can
/// draw), and the explicit residual comes from <see cref="BandedDesign.Apply"/>. The block's closed forms are the
/// same for any design layer; this is the banded one BARS uses.
/// <para>The caller's admission rule judges each factor as soon as it exists, before anything is solved with it: a
/// factor that breaks down can solve to NaN, which no least-squares statistic may be formed from. A refused factor
/// gives no block (null).</para>
/// <para>With <see cref="GPriorCentre.Median"/> or <see cref="GPriorCentre.Mean"/> the prior is centred on the response's
/// level, <c>μ₀ = c·G⁻¹Zᵀ1</c> with c its median or mean, and the statistics are those of <c>ỹ = y − Zμ₀</c>. Both are
/// unweighted even when weights enter, so the prior does not depend on the robust weights a sampler redraws.</para>
/// </summary>
public static class BandedGPrior
{
    /// <summary>The unweighted block: G = ZᵀZ factored, β_OLS, <c>‖ỹ − Zβ_OLS‖²</c> and <c>Q = ‖Zβ_OLS‖²</c> for the
    /// response centred by <paramref name="centre"/>; null when <paramref name="admits"/> refuses G's factor. A caller
    /// that already holds the centre's <see cref="Level"/> passes it as <paramref name="level"/>, which saves the median's
    /// sort.</summary>
    public static CollapsedGPrior? Fit(BandedDesign design, double[] y, Func<ISpdFactor, bool> admits,
                                       GPriorCentre centre = GPriorCentre.Zero, double? level = null)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(admits);
        RefuseFlat(centre);
        int m = design.Rows, nu = design.Dimension;
        if (y.Length != m) throw new ArgumentException("Response length must match the design row count.", nameof(y));

        var band = new double[design.Bandwidth + 1, nu];
        var b = new double[nu];
        design.Accumulate(null, band, y, b);
        var gram = new BandCholesky(nu, design.Bandwidth, BandFactorization.Ldlt, strict: true);
        gram.DecomposeBanded(band);
        if (!admits(gram)) return null;
        (double[]? priorMean, double[] response) = Centre(design, gram, y, centre, level);
        if (priorMean is not null)
        {
            Array.Clear(b);
            design.AccumulateRhs(null, response, b);
        }
        double[] beta = gram.Solve(b);   // β_OLS = G⁻¹Zᵀỹ; no inverse is formed

        var f = new double[m];
        design.Apply(beta, f);
        double rss = 0.0, fitted = 0.0;
        for (int i = 0; i < m; i++)
        {
            double r = response[i] - f[i];
            rss += r * r;
            fitted += f[i] * f[i];
        }
        return CollapsedGPrior.FromLeastSquares(gram, beta, rss, fitted, m, priorMean);
    }

    /// <summary>
    /// The weighted block at prior precision <paramref name="tau"/>: G = ZᵀZ and <c>P_a = Zᵀ(W + τI)Z</c> factored,
    /// the posterior mean <c>β̃ = P_a⁻¹ZᵀWy</c>, and <c>S = ‖y − Zβ̃‖²_W + τ‖Zβ̃‖²</c>, a sum of non-negative terms;
    /// null when <paramref name="admits"/> refuses G.
    /// <para>Admission reads G alone, so it depends on the design and not on the weights, which a robust sampler redraws:
    /// a weight update can never move the current knots out of the support. P_a ≥ τG is positive definite whenever G is,
    /// and its conditioning is within a factor of about (w_max + τ)/(w_min + τ) of G's; a P_a that still fails to factor
    /// is a numerical defect, raised as an exception.</para>
    /// </summary>
    public static CollapsedGPrior? Fit(BandedDesign design, double[] y, double[] weights, double tau,
                                       Func<ISpdFactor, bool> admits, GPriorCentre centre = GPriorCentre.Zero,
                                       double? level = null)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(weights);
        ArgumentNullException.ThrowIfNull(admits);
        RefuseFlat(centre);
        int m = design.Rows, nu = design.Dimension;
        if (y.Length != m) throw new ArgumentException("Response length must match the design row count.", nameof(y));
        if (weights.Length != m) throw new ArgumentException("Weights length must match the design row count.", nameof(weights));

        var gBand = new double[design.Bandwidth + 1, nu];
        design.Accumulate(null, gBand, null, null);
        var gram = new BandCholesky(nu, design.Bandwidth, BandFactorization.Ldlt, strict: true);
        gram.DecomposeBanded(gBand);
        if (!admits(gram)) return null;
        (double[]? priorMean, double[] response) = Centre(design, gram, y, centre, level);

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
        design.AccumulateRhs(weights, response, h);
        var posterior = new BandCholesky(nu, design.Bandwidth, BandFactorization.Ldlt, strict: true);
        posterior.DecomposeBanded(pBand);
        if (!posterior.Status.PositiveDefinite)
            throw new InvalidOperationException(
                $"P_a = Zᵀ(W + τI)Z failed to factor at pivot {posterior.FailedAt} although G was admitted; P_a ≥ τG must be positive definite.");
        double[] beta = posterior.Solve(h);

        var f = new double[m];
        design.Apply(beta, f);
        double weightedRss = 0.0, fitted = 0.0;
        for (int i = 0; i < m; i++)
        {
            double r = response[i] - f[i];
            weightedRss += weights[i] * r * r;
            fitted += f[i] * f[i];
        }
        return CollapsedGPrior.FromWeighted(gram, posterior, beta, weightedRss + tau * fitted, sumLogW, m, tau, priorMean);
    }

    // The prior mean and the response it leaves: none and y at zero; for a level c, the median or the mean,
    // μ₀ = c·G⁻¹Zᵀ1 and y − Zμ₀. The constant is kept exactly where it can be. A constant response takes its value as c,
    // not a rounded average. A basis whose functions sum to one (B-splines, steps) holds the constant as μ₀ = c·1 with no
    // solve, so a constant response centres to exactly zero. Another basis projects; when its span holds the constant (the projection
    // reproduces 1 to 1e-10), a constant response's centred value is zero by definition and is set so, rather than left
    // as the projection's roundoff, which the improper σ² prior would read as a near-perfect fit.
    private static (double[]? PriorMean, double[] Response) Centre(BandedDesign design, ISpdFactor gram, double[] y,
                                                                  GPriorCentre centre, double? level)
    {
        if (centre == GPriorCentre.Zero) return (null, y);
        if (centre is not (GPriorCentre.Median or GPriorCentre.Mean)) throw new ArgumentOutOfRangeException(nameof(centre));
        int m = design.Rows, nu = design.Dimension;
        bool constant = Constant(y);
        double mean = level ?? Level(y, centre);

        var unit = new double[nu];
        Array.Fill(unit, 1.0);
        var rowSums = new double[m];
        design.Apply(unit, rowSums);
        bool sumsToOne = true;
        for (int i = 0; i < m && sumsToOne; i++) sumsToOne = Math.Abs(rowSums[i] - 1.0) <= 1e-12;

        var response = new double[m];
        if (sumsToOne)
        {
            var levelMean = new double[nu];
            Array.Fill(levelMean, mean);
            for (int i = 0; i < m; i++) response[i] = y[i] - mean;
            return (levelMean, response);
        }

        var ones = new double[m];
        Array.Fill(ones, 1.0);
        var columnSums = new double[nu];
        design.AccumulateRhs(null, ones, columnSums);
        double[] projection = gram.Solve(columnSums);   // G⁻¹Zᵀ1
        var projected = new double[m];
        design.Apply(projection, projected);            // P·1
        bool spansConstant = true;
        for (int i = 0; i < m && spansConstant; i++) spansConstant = Math.Abs(projected[i] - 1.0) <= 1e-10;
        var priorMean = new double[nu];
        for (int j = 0; j < nu; j++) priorMean[j] = mean * projection[j];
        if (!(constant && spansConstant))
            for (int i = 0; i < m; i++) response[i] = y[i] - mean * projected[i];
        return (priorMean, response);
    }

    /// <summary>
    /// The level <paramref name="centre"/> centres the response on: zero, its median or its mean, and a constant
    /// response's own value exactly. It reads the data alone, so it is the same for every design. About it, the
    /// response's mean square <c>Σ(yᵢ − c)²/m</c> is the scale <see cref="Sigma2Prior.ScaledToResponse"/> takes.
    /// </summary>
    public static double Level(double[] y, GPriorCentre centre)
    {
        ArgumentNullException.ThrowIfNull(y);
        if (y.Length == 0) throw new ArgumentException("The response is empty.", nameof(y));
        RefuseFlat(centre);
        return centre switch
        {
            GPriorCentre.Zero => 0.0,
            GPriorCentre.Median or GPriorCentre.Mean when Constant(y) => y[0],
            GPriorCentre.Median => Median(y),
            GPriorCentre.Mean => Average(y),
            _ => throw new ArgumentOutOfRangeException(nameof(centre)),
        };
    }

    // The flat intercept needs the constant split out of the basis, which B-splines do not do, and the block's
    // flat-intercept form, which it does not have yet; BASS selects it as an option once the block does.
    internal static void RefuseFlat(GPriorCentre centre)
    {
        if (centre == GPriorCentre.Flat)
            throw new ArgumentException("The flat intercept (GPriorCentre.Flat) needs the constant split out of the basis "
                                        + "and the block's flat-intercept form; BARS has neither.", nameof(centre));
    }

    private static bool Constant(double[] y)
    {
        for (int i = 1; i < y.Length; i++)
            if (y[i] != y[0]) return false;
        return true;
    }

    private static double Median(double[] y)
    {
        var sorted = (double[])y.Clone();
        Array.Sort(sorted);
        int half = sorted.Length / 2;
        return sorted.Length % 2 == 1 ? sorted[half] : 0.5 * (sorted[half - 1] + sorted[half]);
    }

    private static double Average(double[] y)
    {
        double s = 0.0;
        foreach (double v in y) s += v;
        return s / y.Length;
    }
}
