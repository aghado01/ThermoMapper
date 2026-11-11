using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's candidate density (E16's seam, the conformal layer's slice 1): log p(y | x_t, draw) for every retained draw and every
/// value y on each new case's grid, beside <see cref="BassPointwise.LogDensity"/>'s training matrix over the same draws, so the
/// conformal layer (<see cref="ConformalBayes"/>) reads both. The draws are pooled in chain order with the training density's
/// chain lengths, draw weights and latent report. The surface f at the cases is prediction's (<see cref="BassSurface"/>, the
/// training basis at new inputs). The density is the training terms' own: N(y | f, σ²), or under a noise law the law's with each
/// observation's latent integrated out (OQ17's scope rule, CB-4's default), at the draw's σ² and its latents' globals. A fit with
/// case weights (<see cref="CaseWeights"/>) reads each new case's weight as its training rows read theirs, N(y | f, σ²/w) at a
/// precision and w copies of N(y | f, σ²) when the weights count rows, in the pointwise terms' arithmetic, so a training case's
/// candidate is its column bit for bit; a working response that moves per draw stays refused. A
/// candidate is one observation of a scalar response, so it pairs with the default units. Under random effects (E21's slice 5) a case
/// in a level some training case took reads the draw's effects, f + Σₖ z̃ₖ·u_{k,g}, as the training terms read theirs, so a training
/// case's candidate is its column bit for bit; a case in a new group (a level no training case took, placed outside every subset, or
/// its grouping's column left out) integrates that group's effect from its prior, N(y | f + the known effects, σ²(1 + Σ s_k²·z̃_k²))
/// under Gaussian noise, and under a noise the mixture over v of N(σβv, σ²(c·v + Σ s_k²·z̃_k²)), by quadrature in ln v on <c>quadrature</c>
/// (<see cref="BassWidenedMixture"/>). Such a call serves the split form only, and its strata lie inside cells (<see cref="BassGroupCells"/>). A law whose density is an integral over
/// the latent (<see cref="IBassNoiseLaw.DensityByQuadrature"/>) costs draws × candidates integrals, so it is taken on
/// <c>quadrature</c> alone, as the pointwise path takes it. On <c>cdf</c> each draw's log CDF rides beside the density, from the same law
/// at the draw (<see cref="IBassNoiseLaw.At"/>) or N(f, σ²), and <see cref="UnitLogCdf"/> gives the training units' (E16, slice 4).
/// Under the latent-posterior reading (the conformal layer's slice 12b) a noise law's candidates also carry
/// <see cref="CandidateLogDensity.ConditionedValues"/>: each candidate's latent drawn per draw from its exact conditional given its residual
/// net of f and the known groups' effects, N(e | σβv′, σ²c·v′), beside the training terms at their retained latents.
/// Thread-safe: it reads the model and the table and writes only its own buffers.
/// </summary>
internal static class BassCandidates
{
    public static CandidateLogDensity LogDensity(BassModel model, BassTermTable terms, IReadOnlyList<IReadOnlyList<BassDraw>> draws,
                                                 IReadOnlyList<double>? weights, PredictionCases cases, IReadOnlyList<double[]> grids,
                                                 int[]? strata, UnseenLevel unseenLevels, bool quadrature = false, bool cdf = false,
                                                 IReadOnlyList<double>? caseWeights = null,
                                                 LatentTreatment observationLatents = LatentTreatment.Integrated, int seed = 0)
    {
        ArgumentNullException.ThrowIfNull(cases);
        ArgumentNullException.ThrowIfNull(grids);
        bool conditioned = BassPointwise.ConditionsLatents(model, observationLatents);
        if (conditioned && cdf)
            throw new ArgumentException("The CDF is read with each observation's latent integrated, as the predictive distribution reads it.", nameof(cdf));
        if (model.Data.Response is FunctionalResponse)
            throw new NotSupportedException("A direct functional fit's unit is a curve; a candidate here is one observation of a scalar response. Whole curves are E17's.");
        IBassNoiseLaw? law = model.Observation.Law;
        if (law is { DensityByQuadrature: true } && !quadrature)
            throw new NotSupportedException("This noise's density is an integral over each observation's latent, draws × candidates integrals of hundreds "
                                            + "of nodes each; pass quadrature: true to take it.");
        bool weighted = law is null && !model.Working.IsUnweighted, counts = model.Working.Meaning == WeightMeaning.Replication;
        if (weighted && caseWeights is null)
            throw new ArgumentException($"A weighted fit's terms read each observation's weight ({model.Working.Meaning}), so each new case needs its own.", nameof(caseWeights));
        if (!weighted && caseWeights is not null)
            throw new ArgumentException("This fit's terms read no weight, so a new case's would state nothing: give none.", nameof(caseWeights));
        (List<BassDraw> pooled, int[] chains) = BassDrawFit.Pool(draws, terms, false, nameof(draws));
        RandomEffectsBlock? block = model.RandomEffects;
        // Under random effects the training terms read the response at f + u_g; the whitening is the block's integral, not an offset.
        if (law is null && block is null && BassPointwise.Moved(model, pooled))
            throw new NotSupportedException("The draws' working responses moved, so the training terms condition on each observation's offset and weights; "
                                            + "a candidate has neither, and the conformal layer reads every unit by one density (CB-4).");
        (CaseInputs inputs, _) = cases.Prepare(model.Transform, model.Inputs, unseenLevels, model.Data.Groupings);
        var surface = new BassSurface(model.Basis, terms, inputs, null);
        int n = surface.Cases;
        if (grids.Count != n) throw new ArgumentException($"One grid per case is required: {grids.Count} grids for {n} cases.", nameof(grids));
        int[][]? groups = null;
        bool[]? fresh = null;
        if (block is not null)
        {
            BassGroupCells cells = BassGroupCells.Of(model, cases, unseenLevels);
            int[] cellOf = cells.Candidates!;
            if (strata is null && Math.Max(cells.TrainingCells, cellOf.Length == 0 ? 0 : cellOf.Max() + 1) > 1)
                throw new ArgumentException("Under random effects cases are exchangeable only within a cell, one level of every grouping jointly, and these "
                                            + "cases and the training cases span more than one: give each case its cell's stratum (BassFit.GroupCells).", nameof(strata));
            if (strata is not null)
            {
                if (strata.Length != n) throw new ArgumentException("One stratum per case is required.", nameof(strata));
                BassGroupCells.CheckStrata(strata, t => [cellOf[t]], nameof(strata));
            }
            groups = [.. Enumerable.Range(0, block.Count).Select(k => cases.GroupsOf(block.Layout(k), block.Groupings[k], unseenLevels))];
            fresh = cells.NewGroup!;
            if (conditioned && fresh.Any(x => x))
                throw new NotSupportedException("A new group's effect is integrated from its prior, and beside a candidate's drawn latent its conditional "
                                                + "given y is not the latent's alone; the latent-posterior reading takes cases in known groups.");
            if (law is not null && fresh.Any(x => x) && !quadrature)
                throw new NotSupportedException("A new group under a noise is the mixture over v of N(σβv, σ²(c·v + Σ s²z̃²)), an integral over v for every "
                                                + "draw and candidate (BassWidenedMixture); pass quadrature: true to take it.");
        }
        double[]? w = caseWeights is null ? null : Checked(caseWeights, n);
        var offsets = new int[n + 1];
        for (int t = 0; t < n; t++)
            offsets[t + 1] = checked(offsets[t] + (grids[t] ?? throw new ArgumentException($"Case {t}'s grid is null.", nameof(grids))).Length);

        var values = new double[pooled.Count, offsets[n]];
        double[,]? logCdfs = cdf ? new double[pooled.Count, offsets[n]] : null;
        // Under the latent-posterior reading, each candidate's latent drawn from its exact conditional given its residual, in draw, case and grid
        // order on one stream, and its term at that latent.
        double[,]? conditionedValues = conditioned ? new double[pooled.Count, offsets[n]] : null;
        Xoshiro256PlusPlus? rng = conditioned ? new Xoshiro256PlusPlus(seed) : null;
        var f = new double[n];
        double[]? widened = fresh is not null && fresh.Any(x => x) ? new double[n] : null;
        for (int s = 0; s < pooled.Count; s++)
        {
            BassDraw draw = pooled[s];
            surface.Evaluate(draw.Structure.TermIds, draw.Coefficients, f);
            if (groups is not null)
            {
                // The known groups' effects as the training fit adds them (BassDrawFit.Evaluate); a new group's prior variance s_k²·z̃_k² in σ²'s units.
                BassGroupDraw[] effects = block!.Ordered(draw.Groups ?? throw new ArgumentException("A random-effect fit's draw carries its named groups.", nameof(draws)));
                if (widened is not null) Array.Clear(widened);
                for (int k = 0; k < groups.Length; k++)
                    for (int t = 0; t < n; t++)
                    {
                        double z = block.Layout(k).Multiplier(inputs, t);
                        if (groups[k][t] >= 0) f[t] += z*effects[k].Effects[groups[k][t]];
                        else widened![t] += effects[k].Scale * effects[k].Scale * z * z;
                    }
            }
            double sigma2 = draw.Sigma2;
            ObservationLatents? globals = law is null ? null : draw.LatentsOn(model);
            Func<double, double, double, double>? drawLatent = null;
            double unitVariance = 0.0, unitDrift = 0.0;
            if (conditioned)
            {
                var mixture = (VarianceMeanMixture)model.Observation;
                drawLatent = ExactDraw(mixture.MixingAt(globals), rng!);
                unitVariance = sigma2 * mixture.VarianceScale;
                unitDrift = Math.Sqrt(sigma2) * (globals?.Skew ?? 0.0);
            }
            IContinuousDistribution? noise = cdf && law is not null ? law.At(sigma2, globals) : null;
            double normalizer = -0.5 * Math.Log(2.0 * Math.PI * sigma2), sigma = Math.Sqrt(sigma2);
            for (int t = 0; t < n; t++)
            {
                double[] grid = grids[t];
                int at = offsets[t];
                if (fresh is not null && fresh[t] && law is not null)
                {
                    // A new group under a noise: the mixture over v of N(σβv, σ²(c·v + Σ s_k²·z̃_k²)), by quadrature in ln v.
                    var mixture = (VarianceMeanMixture)model.Observation;
                    IMixingLaw mixing = mixture.MixingAt(globals);
                    double unit = sigma2 * mixture.VarianceScale, drift = Math.Sqrt(sigma2) * (globals?.Skew ?? 0.0), extra = sigma2 * widened![t];
                    for (int g = 0; g < grid.Length; g++)
                    {
                        double e = grid[g] - f[t];
                        values[s, at + g] = BassWidenedMixture.LogDensity(mixing, e, unit, drift, extra);
                        if (logCdfs is not null) logCdfs[s, at + g] = BassWidenedMixture.LogCdf(mixing, e, unit, drift, extra);
                    }
                    continue;
                }
                if (fresh is not null && fresh[t])
                {
                    // A new group under Gaussian noise: N(y | f + the known effects, σ²(1 + Σ s_k²·z̃_k²)).
                    double variance = sigma2 * (1.0 + widened![t]), wide = -0.5 * Math.Log(2.0 * Math.PI * variance), sd = Math.Sqrt(variance);
                    for (int g = 0; g < grid.Length; g++)
                    {
                        double e = grid[g] - f[t];
                        values[s, at + g] = BassPointwise.NormalTerm(e, variance, wide);
                        if (logCdfs is not null) logCdfs[s, at + g] = Normal.StandardLogCdf(e / sd);
                    }
                    continue;
                }
                for (int g = 0; g < grid.Length; g++)
                {
                    double e = grid[g] - f[t];
                    values[s, at + g] = law is not null ? law.LogDensity(e, sigma2, globals)
                        : w is null ? BassPointwise.NormalTerm(e, sigma2, normalizer) : BassPointwise.WeightedNormalTerm(e, sigma2, normalizer, w[t], counts);
                    if (drawLatent is not null)
                        conditionedValues![s, at + g] = BassPointwise.ConditionedTerm(e, drawLatent(e, unitVariance, unitDrift), unitVariance, unitDrift);
                    if (logCdfs is not null)
                        logCdfs[s, at + g] = noise is not null ? LogCdf(noise, e)
                            : w is null ? Normal.StandardLogCdf(e / sigma) : BassPointwise.WeightedNormalLogCdf(e, sigma, w[t], counts);
                }
            }
        }
        List<LatentReport> latents = BassPointwise.Latents(model, pooled, collapsed: false,
                                                           conditioned ? LatentTreatment.Conditioned : LatentTreatment.Integrated);
        bool anyNew = fresh is not null && fresh.Any(x => x);
        if (anyNew)
        {
            // A new group's effect is integrated where the training units' are conditioned on (agreed with the conformal layer, E21's slice 5).
            LatentTreatment treatment = groups!.All(group => group.All(g => g < 0)) ? LatentTreatment.Integrated : LatentTreatment.Mixed;
            for (int i = 0; i < latents.Count; i++)
                if (latents[i].Name == BassPointwise.GroupEffects) latents[i] = latents[i] with { Treatment = treatment };
        }
        CandidateLogDensity density = weights is null
            ? new CandidateLogDensity(values, grids, latents, strata, chains) { LogCdf = logCdfs, ConditionedValues = conditionedValues }
            : new CandidateLogDensity(values, grids, latents, strata) { Weights = weights, LogCdf = logCdfs, ConditionedValues = conditionedValues };
        if (!anyNew) return density;
        return density with
        {
            NewGroup = fresh,
            Forms = ConformalForms.Split,
            FormsReason = "a new group's effect is integrated from its prior while the training units' terms condition on theirs, so the full step's "
                          + "in-sample scores would read the candidate by another density, which undercovers (CB-4).",
        };
    }

    /// <summary>A latent v from its exact conditional given the residual (e, s², m) at likelihood power 1: Student t's, the exponential's
    /// and the GIG's draws. A law whose conditional is a slice step from the current value has none, and is refused.</summary>
    private static Func<double, double, double, double> ExactDraw(IMixingLaw mixing, Xoshiro256PlusPlus rng) => mixing switch
    {
        StudentTMixing t => (e, s2, m) => t.SampleConditional(e, s2, m, 1.0, rng),
        ExponentialMixing x => (e, s2, m) => x.SampleConditional(e, s2, m, 1.0, rng),
        GigMixing g => (e, s2, m) => g.SampleConditional(e, s2, m, 1.0, rng),
        _ => throw new NotSupportedException("This noise's latent has no exact draw from its conditional, only a step from its current value, so the "
                                             + "latent-posterior reading cannot draw a candidate's latent."),
    };

    /// <summary>
    /// Each draw's log CDF at each training unit's response, ln F(y_u | x_u, θ_s), draws × units in the pooled draws' rows, for the
    /// training density's <see cref="PointwiseLogDensity.LogCdf"/>: the surface at the training cases (<see cref="BassDrawFit"/>), under
    /// the noise's law at the draw (<see cref="IBassNoiseLaw.At"/>) or N(f, σ²), each unit one observation. Under case weights it is
    /// N(f, σ²/w) at a precision and N(f, σ²) when the weights count rows, each copy one observation of that law. The candidates' log CDF
    /// reads the same law, so every unit is read by one distribution (CB-4). The fits the candidate density refuses are refused here too.
    /// </summary>
    public static double[,] UnitLogCdf(BassModel model, BassTermTable terms, IReadOnlyList<IReadOnlyList<BassDraw>> draws, UnitPartition units)
    {
        ArgumentNullException.ThrowIfNull(units);
        if (model.Data.Response is FunctionalResponse)
            throw new NotSupportedException("A direct functional fit's unit is a curve; a CDF orders one observation of a scalar response.");
        int n = model.Data.Response.Observations;
        if (units.Observations != n) throw new ArgumentException("The partition must cover the model's N observations.", nameof(units));
        for (int u = 0; u < units.Units; u++)
            if (units.Members[u].Length != 1) throw new ArgumentException($"Unit {u} holds {units.Members[u].Length} observations; a CDF orders one value.", nameof(units));
        IBassNoiseLaw? law = model.Observation.Law;
        (List<BassDraw> pooled, _) = BassDrawFit.Pool(draws, terms, false, nameof(draws));
        bool effects = model.RandomEffects is not null;
        if (law is null && !effects && BassPointwise.Moved(model, pooled))
            throw new NotSupportedException("The draws' working responses moved, so the training terms condition on each observation's offset and weights; "
                                            + "a candidate has neither, and the conformal layer reads every unit by one distribution (CB-4).");
        var fit = new BassDrawFit(model, terms);
        // Under random effects the fit reads the response itself at f + u_g, as the training terms do; the working response is whitened.
        ReadOnlySpan<double> y = law is null && !effects ? model.Working.Values : model.Data.Response.Values;
        double[]? weights = law is null && !effects ? model.Working.WeightArray : null;
        bool counts = model.Working.Meaning == WeightMeaning.Replication;
        var logCdfs = new double[pooled.Count, units.Units];
        var f = new double[n];
        for (int s = 0; s < pooled.Count; s++)
        {
            BassDraw draw = pooled[s];
            fit.Evaluate(draw, f);
            double sigma = Math.Sqrt(draw.Sigma2);
            IContinuousDistribution? noise = law?.At(draw.Sigma2, draw.LatentsOn(model));
            for (int u = 0; u < units.Units; u++)
            {
                int i = units.Members[u][0];
                double e = y[i] - f[i];
                logCdfs[s, u] = noise is not null ? LogCdf(noise, e)
                    : weights is null ? Normal.StandardLogCdf(e / sigma) : BassPointwise.WeightedNormalLogCdf(e, sigma, weights[i], counts);
            }
        }
        return logCdfs;
    }

    // The new cases' weights, one positive and finite value per case, copied.
    private static double[] Checked(IReadOnlyList<double> caseWeights, int cases)
    {
        if (caseWeights.Count != cases)
            throw new ArgumentException($"One weight per case is required: {caseWeights.Count} weights for {cases} cases.", nameof(caseWeights));
        var w = new double[cases];
        for (int t = 0; t < cases; t++)
        {
            w[t] = caseWeights[t];
            if (!(w[t] > 0.0 && double.IsFinite(w[t]))) throw new ArgumentException($"Case {t}'s weight is {w[t]:R}; a weight is positive and finite.", nameof(caseWeights));
        }
        return w;
    }

    /// <summary>ln F(e) under a noise at a draw, each half on its own tail (<see cref="NoiseLogCdf"/>, which every engine's seam shares).</summary>
    internal static double LogCdf(IContinuousDistribution noise, double e) => NoiseLogCdf.Of(noise, e);
}
