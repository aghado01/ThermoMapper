using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// The composition of a retained draw of several outputs (E34, D04's composition carried to K outputs), each piece exact given
/// (structure, τ): Σ | Y ∼ IW(ν₀ + N, Ψ₀ + S) from the K × K block (<see cref="OutputMarginal.SigmaPosterior"/>), then the centred
/// columns' coefficients A | Σ ∼ MN(Ā, (τ + 1)⁻¹G_c⁻¹, Σ) with Ā = β̂/(τ + 1), vec stacking A's columns so Cov(vec A) = Σ ⊗ (τ + 1)⁻¹G_c⁻¹,
/// and the intercepts apart from them, α | Σ ∼ N(m, Σ·c₀/N): m = ℓ + (ȳ − ℓ)/(τ + 1) and c₀ = 1/(τ + 1) inside the g-prior, m = ȳ and
/// c₀ = 1 under the flat intercept. The stream is read as the scalar draw reads it: Σ's gamma variates and normals, then each output's
/// M normals through the factor's innovation, mixed across the outputs by Σ's factor, then the intercepts' K normals; so at K = 1 the
/// draw is the scalar one, to rounding.
/// </summary>
internal static class BassOutputDraws
{
    /// <summary>Σ and the centred coefficients at <paramref name="design"/>'s structure and τ = <paramref name="tau"/>: row 0 the
    /// intercepts, then one row per term in the design's order, one column per output.</summary>
    internal static (double[,] Sigma, double[,] Centred) Sample(BassTarget target, BassDesign design, BassStructureFit fit, double tau, Xoshiro256PlusPlus rng,
                                                              InverseWishart? prior = null)
    {
        BassResolvedOutputs outputs = target.Outputs ?? throw new InvalidOperationException("The target holds one response.");
        if (outputs.Rates is not null && prior is null)
            throw new ArgumentException("A draw under Huang and Wand's prior needs the conditional inverse Wishart at its state's rates.", nameof(prior));
        int k = outputs.Count, m = design.Count, n = design.Observations;
        double precision = tau / target.CoefficientConstant;
        double[,] sigma = target.OutputsMarginal(fit, tau).SigmaPosterior(prior ?? outputs.Prior).Sample(rng);
        double[,] root = Root(sigma);
        if (fit.Outputs?.Prior is { } pencil)
            return (sigma, SamplePencil(pencil, outputs, root, precision, n, m, rng));
        double shrink = 1.0 / (precision + 1.0), scale = Math.Sqrt(shrink);
        double[,] beta = design.OutputCoefficients ?? throw new InvalidOperationException("The design holds one response.");
        ISpdFactor factor = design.Factor;

        var innovations = new double[k][];
        var z = new double[m];
        for (int l = 0; l < k; l++)
        {
            var v = new double[m];
            if (m > 0)
            {
                for (int j = 0; j < m; j++) z[j] = Normal.SampleStandard(rng);
                factor.SampleInnovation(z, v);
            }
            innovations[l] = v;
        }
        var centred = new double[m + 1, k];
        for (int j = 0; j < m; j++)
            for (int a = 0; a < k; a++)
            {
                double mixed = 0.0;
                for (int l = 0; l <= a; l++) mixed += root[a, l] * innovations[l][j];
                centred[j + 1, a] = shrink * beta[j, a] + scale * mixed;
            }

        bool flat = target.Intercept == GPriorCentre.Flat;
        double interceptScale = Math.Sqrt((flat ? 1.0 : shrink) / n);
        var u = new double[k];
        for (int l = 0; l < k; l++) u[l] = Normal.SampleStandard(rng);
        for (int a = 0; a < k; a++)
        {
            double mean = outputs.Means[a], level = outputs.Levels[a];
            double centre = flat ? mean : level + shrink * (mean - level);
            double mixed = 0.0;
            for (int l = 0; l <= a; l++) mixed += root[a, l] * u[l];
            centred[0, a] = centre + interceptScale * mixed;
        }
        return (sigma, centred);
    }

    private static double[,] SamplePencil(OutputPencil pencil, BassResolvedOutputs outputs, double[,] root, double tau,
                                          int rows, int terms, Xoshiro256PlusPlus rng)
    {
        var (factor, mean) = pencil.Conditional(tau, 1);
        int p = factor.Dimension, k = outputs.Count, offset = p == terms ? 0 : 1;
        var innovations = new double[k][];
        var z = new double[p];
        for (int l = 0; l < k; l++)
        {
            for (int j = 0; j < p; j++) z[j] = Normal.SampleStandard(rng);
            innovations[l] = new double[p];
            factor.SampleInnovation(z, innovations[l]);
        }
        var centred = new double[terms + 1, k];
        for (int j = 0; j < p; j++)
            for (int a = 0; a < k; a++)
            {
                double value = mean[j, a];
                for (int l = 0; l <= a; l++) value += root[a, l] * innovations[l][j];
                centred[j + 1 - offset, a] = value;
            }
        if (offset == 0)
        {
            var intercept = new double[k];
            for (int l = 0; l < k; l++) intercept[l] = Normal.SampleStandard(rng);
            for (int a = 0; a < k; a++)
            {
                double value = outputs.Means[a];
                for (int l = 0; l <= a; l++) value += root[a, l] * intercept[l] / Math.Sqrt(rows);
                centred[0, a] = value;
            }
        }
        return centred;
    }

    // The lower Cholesky factor of a draw of Σ, positive definite by construction.
    private static double[,] Root(double[,] sigma)
    {
        int k = sigma.GetLength(0);
        var l = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = sigma[i, j];
                for (int r = 0; r < j; r++) sum -= l[i, r] * l[j, r];
                if (i == j)
                {
                    if (!(sum > 0.0)) throw new InvalidOperationException("A draw of Σ is not positive definite in doubles.");
                    l[i, i] = Math.Sqrt(sum);
                }
                else l[i, j] = sum / l[j, j];
            }
        return l;
    }
}
