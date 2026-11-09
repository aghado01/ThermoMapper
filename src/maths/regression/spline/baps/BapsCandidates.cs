using System;
using System.Collections.Generic;
using Maths.Distributions;

namespace Maths.Regression.Spline.Baps;

public sealed partial class Baps
{
    /// <summary>
    /// BAPS's candidate density (E16's seam, the conformal layer's slice 5): log p(y | z_t, draw) = N(y; z_tβ, σ²) for every retained
    /// draw and every value y on each new case's grid, beside <see cref="PointwiseLogDensity"/>'s matrix over the same draws in the same
    /// rows, with its chain lengths and latent report, so the conformal layer (<see cref="Conformal.ConformalBayes"/>) reads both. A
    /// candidate is one observation; its design row is the fitted basis at the new case, evaluated as the training design was. On
    /// <paramref name="cdf"/> each draw's log CDF, ln Φ((y − z_tβ)/σ), rides beside the density.
    /// </summary>
    /// <param name="result">A run of this model, whose retained draws are read.</param>
    /// <param name="design">The new cases' design rows (cases × ν), of the basis the model was built on.</param>
    /// <param name="grids">Each case's candidate values, ascending.</param>
    /// <param name="strata">Each case's stratum in the training density's partition; null puts every case in stratum 0.</param>
    /// <param name="cdf">Adds each draw's log CDF at each candidate, for the conformal predictive distribution.</param>
    public CandidateLogDensity CandidateLogDensity(BapsResult result, double[,] design, IReadOnlyList<double[]> grids, int[]? strata = null,
                                                   bool cdf = false)
    {
        ArgumentNullException.ThrowIfNull(result);
        ArgumentNullException.ThrowIfNull(design);
        if (design.GetLength(1) != _nu) throw new ArgumentException($"The design has {design.GetLength(1)} columns; the model's basis has {_nu}.", nameof(design));
        int cases = design.GetLength(0);
        int[] offsets = CandidateSeam.Offsets(grids, cases);
        var values = new double[result.Draws, offsets[cases]];
        double[,]? logCdfs = cdf ? new double[result.Draws, offsets[cases]] : null;
        int row = 0;
        foreach (IReadOnlyList<BapsDraw> chain in result.PosteriorDraws)
            foreach (BapsDraw d in chain)
            {
                if (d.Coefficients.Length != _nu)
                    throw new ArgumentException("The result's draws do not fit this model's basis.", nameof(result));
                double logNorm = -0.5 * Math.Log(2.0 * Math.PI * d.Sigma2), sigma = Math.Sqrt(d.Sigma2);
                for (int t = 0; t < cases; t++)
                {
                    double fit = Fit(design, t, d.Coefficients);
                    double[] grid = grids[t];
                    for (int g = 0; g < grid.Length; g++)
                    {
                        double r = grid[g] - fit;
                        values[row, offsets[t] + g] = logNorm - 0.5 * r * r / d.Sigma2;
                        if (logCdfs is not null) logCdfs[row, offsets[t] + g] = Normal.StandardLogCdf(r / sigma);
                    }
                }
                row++;
            }
        return new CandidateLogDensity(values, grids, Latents(), strata, ChainLengths.Of(result.PosteriorDraws)) { LogCdf = logCdfs };
    }

    // z_iβ, summed in column order: the training terms' and the candidates' one arithmetic.
    private static double Fit(double[,] z, int i, double[] beta)
    {
        double fit = 0.0;
        for (int j = 0; j < beta.Length; j++) fit += z[i, j] * beta[j];
        return fit;
    }

    // A draw's coefficients and σ² are global, so every unit, a candidate too, conditions on them.
    private static List<LatentReport> Latents() =>
    [
        new("coefficients β", LatentScope.Global, LatentTreatment.Conditioned),
        new("noise variance σ²", LatentScope.Global, LatentTreatment.Conditioned),
    ];
}
