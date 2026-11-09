using System;
using System.Collections.Generic;
using Maths.Distributions;

namespace Maths.Regression.Spline.Baps;

public sealed partial class HierarchicalPSpline
{
    /// <summary>
    /// The hierarchical P-spline's candidate density (E16's seam, the conformal layer's slice 5): log p(y | z_t, draw) = N(y; z_tβ_r, σ²),
    /// a new point on a seen replicate r's curve, for every retained draw and every value y on each new case's grid, in
    /// <see cref="PointwiseLogDensity"/>'s rows with its chain lengths. It conditions on the drawn curve β_r, as a training unit that
    /// splits its replicate does, and its latent report says so. On <paramref name="cdf"/> each draw's log CDF, ln Φ((y − z_tβ_r)/σ),
    /// rides beside the density.
    /// <para><b>The split form takes it</b> (<see cref="Conformal.SplitConformal"/>, CB-2): calibration points and test points on seen curves,
    /// each at an input drawn as the others are, are exchangeable with one another whatever the fit. <b>The full step's guarantee does not
    /// reach it.</b> The model takes one design for every replicate, so each training point's input is held by every curve, while a new
    /// point's input is held by none: the candidate is not exchangeable with the training points (CB-1's premise), and no partition
    /// restores it. A point on a new curve would integrate that curve while the training units condition on theirs (CB-4), and a whole
    /// new curve on the shared design, which is exchangeable with the training curves, is a vector, E17's.</para>
    /// </summary>
    /// <param name="result">A run of this model with replicate draws retained (<c>retainReplicates: true</c>).</param>
    /// <param name="design">The new cases' design rows (cases × ν), of the basis the model was built on.</param>
    /// <param name="replicates">Each case's replicate, numbered as the responses were.</param>
    /// <param name="grids">Each case's candidate values, ascending.</param>
    /// <param name="strata">Each case's stratum, the calibration's numbering; null puts every case in stratum 0.</param>
    /// <param name="cdf">Adds each draw's log CDF at each candidate, for the conformal predictive distribution.</param>
    public CandidateLogDensity CandidateLogDensity(HierarchicalResult result, double[,] design, int[] replicates, IReadOnlyList<double[]> grids,
                                                   int[]? strata = null, bool cdf = false)
    {
        ArgumentNullException.ThrowIfNull(result);
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(replicates);
        if (design.GetLength(1) != _nu) throw new ArgumentException($"The design has {design.GetLength(1)} columns; the model's basis has {_nu}.", nameof(design));
        int cases = design.GetLength(0);
        if (replicates.Length != cases) throw new ArgumentException("One replicate per case is required.", nameof(replicates));
        for (int t = 0; t < cases; t++)
            if (replicates[t] < 0 || replicates[t] >= _responses.Length)
                throw new ArgumentOutOfRangeException(nameof(replicates),
                    $"Case {t}'s replicate {replicates[t]} is not one the model fitted. A point on a new curve integrates that curve under the "
                    + "population's law, while the training units that split a curve condition on it (CB-4); a whole new curve is E17's.");
        int[] offsets = CandidateSeam.Offsets(grids, cases);

        int draws = 0;
        foreach (IReadOnlyList<HierarchicalDraw> chain in result.PosteriorDraws) draws += chain.Count;
        var values = new double[draws, offsets[cases]];
        double[,]? logCdfs = cdf ? new double[draws, offsets[cases]] : null;
        int row = 0;
        foreach (IReadOnlyList<HierarchicalDraw> chain in result.PosteriorDraws)
            foreach (HierarchicalDraw d in chain)
            {
                if (d.Replicates is null)
                    throw new InvalidOperationException("A point on a seen curve conditions on that curve, which needs retained replicate draws (retainReplicates: true).");
                double s2 = d.NoiseVariance, sigma = Math.Sqrt(s2);
                for (int t = 0; t < cases; t++)
                {
                    double f = CurveAt(design, t, d.Replicates[replicates[t]]);
                    double[] grid = grids[t];
                    for (int g = 0; g < grid.Length; g++)
                    {
                        double e = grid[g] - f;
                        values[row, offsets[t] + g] = -0.5 * Math.Log(2.0 * Math.PI * s2) - 0.5 * e * e / s2;
                        if (logCdfs is not null) logCdfs[row, offsets[t] + g] = Normal.StandardLogCdf(e / sigma);
                    }
                }
                row++;
            }
        return new CandidateLogDensity(values, grids, Latents(LatentTreatment.Conditioned), strata, ChainLengths.Of(result.PosteriorDraws))
        {
            LogCdf = logCdfs,
        };
    }

    // z_iβ, summed in column order: the training terms' and the candidates' one arithmetic.
    private static double CurveAt(double[,] z, int i, double[] beta)
    {
        double f = 0.0;
        for (int j = 0; j < beta.Length; j++) f += z[i, j] * beta[j];
        return f;
    }

    // The replicate curves as the units treated them; the population curve and the two variances are global, always conditioned on.
    private static LatentReport[] Latents(LatentTreatment curves) =>
    [
        new LatentReport("replicate curves β_r", LatentScope.Group, curves),
        new LatentReport("population curve β_0", LatentScope.Global, LatentTreatment.Conditioned),
        new LatentReport("noise variance σ²", LatentScope.Global, LatentTreatment.Conditioned),
        new LatentReport("deviation variance τ_u²", LatentScope.Global, LatentTreatment.Conditioned),
    ];
}
