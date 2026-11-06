using System;
using System.Collections.Generic;
using Maths.Regression.Conformal;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's candidate density (E16's seam, the conformal layer's slice 1): log p(y | x_t, draw) for every retained draw and every
/// value y on each new case's grid, beside <see cref="BassPointwise.LogDensity"/>'s training matrix over the same draws, so the
/// conformal layer (<see cref="ConformalBayes"/>) reads both. The draws are pooled in chain order with the training density's
/// chain lengths, draw weights and latent report. The surface f at the cases is prediction's (<see cref="BassSurface"/>, the
/// training basis at new inputs). The density is the training terms' own: N(y | f, σ²), or under a noise law the law's with each
/// observation's latent integrated out (OQ17's scope rule, CB-4's default), at the draw's σ² and its latents' globals. A
/// candidate is one observation of a scalar response, so it pairs with the default units. A law whose density is an integral over
/// the latent (<see cref="IBassNoiseLaw.DensityByQuadrature"/>) costs draws × candidates integrals, so it is taken on
/// <c>quadrature</c> alone, as the pointwise path takes it. Thread-safe: it reads the model and the table and writes only its own
/// buffers.
/// </summary>
internal static class BassCandidates
{
    public static CandidateLogDensity LogDensity(BassModel model, BassTermTable terms, IReadOnlyList<IReadOnlyList<BassDraw>> draws,
                                                 IReadOnlyList<double>? weights, PredictionCases cases, IReadOnlyList<double[]> grids,
                                                 int[]? strata, UnseenLevel unseenLevels, bool quadrature = false)
    {
        ArgumentNullException.ThrowIfNull(cases);
        ArgumentNullException.ThrowIfNull(grids);
        if (model.Data.Response is FunctionalResponse)
            throw new NotSupportedException("A direct functional fit's unit is a curve; a candidate here is one observation of a scalar response. Whole curves are E17's.");
        IBassNoiseLaw? law = model.Observation.Law;
        if (law is { DensityByQuadrature: true } && !quadrature)
            throw new NotSupportedException("This noise's density is an integral over each observation's latent, draws × candidates integrals of hundreds "
                                            + "of nodes each; pass quadrature: true to take it.");
        if (law is null && !model.Working.IsUnweighted)
            throw new NotSupportedException("A weighted fit's terms read each observation's weight, and a new case carries none, so its candidates are not read as its units are.");
        (List<BassDraw> pooled, int[] chains) = BassDrawFit.Pool(draws, terms, false, nameof(draws));
        if (law is null && pooled.Exists(d => !ReferenceEquals(d.WorkingOn(model), model.Working)))
            throw new NotSupportedException("The draws' working responses moved, so the training terms condition on each observation's offset and weights; "
                                            + "a candidate has neither, and the conformal layer reads every unit by one density (CB-4).");
        (CaseInputs inputs, _) = cases.Prepare(model.Transform, model.Inputs, unseenLevels);
        var surface = new BassSurface(model.Basis, terms, inputs, null);
        int n = surface.Cases;
        if (grids.Count != n) throw new ArgumentException($"One grid per case is required: {grids.Count} grids for {n} cases.", nameof(grids));
        var offsets = new int[n + 1];
        for (int t = 0; t < n; t++)
            offsets[t + 1] = checked(offsets[t] + (grids[t] ?? throw new ArgumentException($"Case {t}'s grid is null.", nameof(grids))).Length);

        var values = new double[pooled.Count, offsets[n]];
        var f = new double[n];
        for (int s = 0; s < pooled.Count; s++)
        {
            BassDraw draw = pooled[s];
            surface.Evaluate(draw.Structure.TermIds, draw.Coefficients, f);
            double sigma2 = draw.Sigma2;
            ObservationLatents? globals = law is null ? null : draw.LatentsOn(model);
            double normalizer = -0.5 * Math.Log(2.0 * Math.PI * sigma2);
            for (int t = 0; t < n; t++)
            {
                double[] grid = grids[t];
                int at = offsets[t];
                for (int g = 0; g < grid.Length; g++)
                {
                    double e = grid[g] - f[t];
                    values[s, at + g] = law is null ? normalizer - 0.5 * e * e / sigma2 : law.LogDensity(e, sigma2, globals);
                }
            }
        }
        List<LatentReport> latents = BassPointwise.Latents(model, pooled, collapsed: false);
        return weights is null
            ? new CandidateLogDensity(values, grids, latents, strata, chains)
            : new CandidateLogDensity(values, grids, latents, strata) { Weights = weights };
    }
}
