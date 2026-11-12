using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>Joint posterior readouts of a several-output fit (E34), on the fit's original data and retained draw weights.</summary>
public static class BassMultiOutput
{
    internal static BassModel Model(BassFit fit, BassData data)
    {
        ArgumentNullException.ThrowIfNull(fit);
        var model = fit.ModelFor(data);
        if (data.Response is not MultiOutputResponse) throw new ArgumentException("This readout needs several outputs.", nameof(data));
        return model;
    }

    internal static double[] Functional(int outputs, IReadOnlyList<double> weights)
    {
        ArgumentNullException.ThrowIfNull(weights);
        if (weights.Count != outputs) throw new ArgumentException("One functional weight per output is required.", nameof(weights));
        var copy = weights.ToArray();
        if (copy.Any(x => !double.IsFinite(x))) throw new ArgumentException("The output functional must be finite.", nameof(weights));
        return copy;
    }

    internal static double[] Project(BassDraw draw, IReadOnlyList<double> functional)
    {
        var matrix = draw.OutputCoefficients ?? throw new ArgumentException("The draw has no output coefficients.", nameof(draw));
        if (matrix.GetLength(1) != functional.Count || matrix.GetLength(0) != draw.Structure.TermIds.Length + 1)
            throw new ArgumentException("The output coefficients must match the structure and functional.", nameof(draw));
        var result = new double[matrix.GetLength(0)];
        for (int j = 0; j < result.Length; j++) for (int k = 0; k < functional.Count; k++) result[j] += matrix[j, k] * functional[k];
        return result;
    }

    public static BassOutputPrediction Predict(BassFit fit, BassData data, PredictionCases cases, BassPredictOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(cases);
        BassModel model = Model(fit, data); options ??= new BassPredictOptions();
        if (options.OutputLocations is not null || options.BasisParts != BasisPredictiveParts.NoiseAndTruncation)
            throw new ArgumentException("Several outputs have named columns, without functional locations or truncation noise.", nameof(options));
        var pooled = fit.Draws.SelectMany(c => c).ToArray();
        var (selected, weights) = Draws3.Select(options.Draws, pooled.Length, fit.DrawWeights, nameof(options));
        var (inputs, outside) = cases.Prepare(model.Transform, model.Inputs, options.UnseenLevels, model.Data.Groupings);
        // Random intercepts (E34's slice 4): a case of an observed level reads its draw's u_g in R^K; one of the population, u = 0 in f, and its
        // own effect ~ N(0, s^2 Sigma) in the predictive draw, as scalar prediction reads one response's.
        RandomEffectsBlock? block = model.RandomEffects;
        int[][]? groups = block is null ? null : [.. Enumerable.Range(0, block.Count).Select(g => cases.GroupsOf(block.Layout(g), block.Groupings[g], options.UnseenLevels))];
        var surface = new BassSurface(model.Basis, fit.Terms, inputs, null);
        int n = surface.Cases, k = data.Response.Outputs;
        var latent = new double[selected.Length, n, k];
        var predictive = options.Predictive ? new double[selected.Length, n, k] : null;
        var sigma = new double[selected.Length][,];
        var rng = new Xoshiro256PlusPlus(options.Seed);
        var column = new double[n]; var coefficients = new double[0]; var z = new double[k];
        // A noise by row (E34's decision 6): a new case's outputs share one v from the law at the draw's tail, yᵢ = fᵢ + vγ + √(c·v)·L_Σz.
        var mixture = model.Observation as VarianceMeanMixture;
        for (int d = 0; d < selected.Length; d++)
        {
            BassDraw draw = pooled[selected[d]];
            if (draw.OutputCoefficients is not { } a || draw.Sigma is not { } s || a.GetLength(1) != k || a.GetLength(0) != draw.Structure.TermIds.Length + 1 || s.GetLength(0) != k || s.GetLength(1) != k)
                throw new ArgumentException("A several-output draw must carry its coefficients and covariance.", nameof(fit));
            sigma[d] = (double[,])s.Clone();
            Array.Resize(ref coefficients, a.GetLength(0));
            BassGroupDraw[]? effects = groups is null ? null : block!.Ordered(draw.Groups ?? throw new ArgumentException("A draw under random intercepts carries its effects.", nameof(fit)));
            for (int l = 0; l < k; l++)
            {
                for (int j = 0; j < coefficients.Length; j++) coefficients[j] = a[j, l];
                surface.Evaluate(draw.Structure.TermIds, coefficients, column);
                if (effects is not null)
                    for (int g = 0; g < groups!.Length; g++)
                        for (int i = 0; i < n; i++)
                            if (groups[g][i] >= 0) column[i] += block!.Layout(g).Multiplier(inputs, i) * effects[g].OutputEffects![l][groups[g][i]];
                for (int i = 0; i < n; i++) latent[d, i, l] = column[i];
            }
            if (predictive is null) continue;
            var root = MatrixTPredictive.Root(s);
            ObservationLatents? rowNoise = mixture is null ? null
                : draw.LatentsOn(model) ?? throw new ArgumentException("A draw under a noise carries its latents.", nameof(fit));
            IMixingLaw? law = rowNoise is null ? null : mixture!.MixingAt(rowNoise);
            double[]? skew = rowNoise?.OutputSkew;
            for (int i = 0; i < n; i++)
            {
                double v = law?.SamplePrior(rng) ?? 1.0, spread = law is null ? 1.0 : Math.Sqrt(mixture!.VarianceScale * v);
                for (int l = 0; l < k; l++) z[l] = Normal.SampleStandard(rng);
                for (int l = 0; l < k; l++)
                {
                    double noise = skew is null ? 0 : v * skew[l];
                    for (int j = 0; j <= l; j++) noise += spread * root[l, j] * z[j];
                    predictive[d, i, l] = latent[d, i, l] + noise;
                }
            }
            if (effects is null) continue;
            for (int g = 0; g < groups!.Length; g++)   // a population case's own effect, after every case's noise
                for (int i = 0; i < n; i++)
                {
                    if (groups[g][i] >= 0) continue;
                    for (int l = 0; l < k; l++) z[l] = Normal.SampleStandard(rng);
                    double scale = effects[g].Scale * block!.Layout(g).Multiplier(inputs, i);
                    for (int l = 0; l < k; l++) { double own = 0; for (int j = 0; j <= l; j++) own += root[l, j] * z[j]; predictive[d, i, l] += scale * own; }
                }
        }
        double[,]? rb = options.MeanSurface ? MeanSurface(model, fit, pooled, selected, weights, surface, groups, inputs) : null;
        double nu = model.Resolved.Outputs!.Prior.DegreesOfFreedom + data.Response.Cases - (model.Target.Intercept == GPriorCentre.Flat ? 1 : 0);
        return new BassOutputPrediction(((MultiOutputResponse)data.Response).Names, selected, weights, latent, predictive, sigma, outside, rb, nu > k + 1);
    }

    // Only coefficient means are needed here; no case-by-case row covariance is constructed.
    private static double[,] MeanSurface(BassModel model, BassFit fit, BassDraw[] pooled, int[] selected, double[] weights, BassSurface surface,
                                         int[][]? groups = null, CaseInputs? inputs = null)
    {
        int k = model.Data.Response.Outputs, rows = model.Data.Response.Cases;
        bool flat = model.Target.Intercept == GPriorCentre.Flat;
        if (!(model.Resolved.Outputs!.Prior.DegreesOfFreedom + rows - (flat ? 1 : 0) > k))
            throw new InvalidOperationException("The posterior surface has no finite mean.");
        var designs = new Dictionary<(string, WorkingResponse), (BassDesign Design, BassStructureFit Fit, double[] Means)>();
        var result = new double[surface.Cases, k]; var buffer = new double[surface.Cases];
        for (int d = 0; d < selected.Length; d++)
        {
            var draw = pooled[selected[d]]; var ids = draw.Structure.TermIds;
            WorkingResponse reading = draw.WorkingOn(model);   // under random intercepts, the reading at the draw's s (E34's slice 4)
            var key = (string.Join(',', ids), reading);
            if (!designs.TryGetValue(key, out var cached))
            {
                var design = model.NewDesign();
                if (design.Load(ids.Select(id => fit.Terms[id]).ToArray(), reading) is { } refusal) throw new ArgumentException($"The retained structure was refused ({refusal}).", nameof(fit));
                var means = new double[ids.Length]; var training = new double[rows];
                for (int j = 0; j < ids.Length; j++) { model.Basis.CaseColumn(fit.Terms[ids[j]], model.Inputs.Cases, training); means[j] = reading.ColumnMean(training, null); }
                cached = (design, BassStructureFit.Of(design, 0), means); designs.Add(key, cached);
            }
            double precision = draw.Tau / model.Target.CoefficientConstant;
            double[,]? pencilMean = cached.Fit.Outputs!.Prior?.Conditional(precision, 1).Mean;
            var coefficients = new double[ids.Length + 1]; int offset = flat ? 0 : 1;
            bool moved = !reading.IsResponseItself;
            double[] outputMeans = moved ? reading.OutputMeans.ToArray() : model.Resolved.Outputs.Means;
            var all = groups is null ? null : new double[ids.Length + 1, k];
            var columns = new double[k][];
            for (int l = 0; l < k; l++)
            {
                double mean = outputMeans[l], level = moved && reading.Latents is null ? (model.Target.Intercept == GPriorCentre.Zero ? 0.0 : mean) : model.Resolved.Outputs.Levels[l];
                coefficients[0] = flat ? mean : pencilMean is not null ? pencilMean[0, l] : level + (mean - level) / (1 + precision);
                for (int j = 0; j < ids.Length; j++)
                {
                    double value = pencilMean is not null ? pencilMean[j + offset, l] : cached.Design.OutputCoefficients![j, l] / (1 + precision);
                    coefficients[j + 1] = value; coefficients[0] -= cached.Means[j] * value;
                }
                surface.Evaluate(ids, coefficients, buffer);
                columns[l] = (double[])buffer.Clone();
                if (all is not null) for (int j = 0; j < coefficients.Length; j++) all[j, l] = coefficients[j];
            }
            if (all is not null)   // E[u | S, tau, Y], linear in the coefficients: the effects' conditional mean at their conditional means
            {
                RandomEffectsBlock intercepts = model.BlockAt(draw)!;   // at the draw's scales
                double[][][] u = intercepts.MeanOutputEffects(model, [.. ids.Select(id => fit.Terms[id])], all, reading);   // under a noise, at the draw's weights
                for (int l = 0; l < k; l++)
                    for (int g = 0; g < groups!.Length; g++)
                        for (int i = 0; i < surface.Cases; i++)
                            if (groups[g][i] >= 0) columns[l][i] += intercepts.Layout(g).Multiplier(inputs!, i) * u[l][g][groups[g][i]];
            }
            for (int l = 0; l < k; l++)
                for (int i = 0; i < surface.Cases; i++) result[i, l] += weights[d] * columns[l][i];
        }
        return result;
    }

    /// <summary>The matrix-t components conditional on each retained structure, τ and covariance rates, with A and Σ integrated out.</summary>
    public static BassMatrixTPrediction Conditionals(BassFit fit, BassData data, PredictionCases cases, BassPredictOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(cases);
        BassModel model = Model(fit, data); options ??= new BassPredictOptions();
        if (options.OutputLocations is not null) throw new ArgumentException("Several outputs have no output locations.", nameof(options));
        if (model.RandomEffects is not null)
            throw new NotSupportedException("Under random intercepts a new case's row correlates with the training rows through u, which the matrix-t's row "
                                            + "scale does not carry; Predict reads each draw's effects (E34's slice 4).");
        if (model.Observation is VarianceMeanMixture)
            throw new NotSupportedException("Under a noise by row a new case's outputs are a variance-mean mixture, not a matrix-t's row; Predict draws "
                                            + "them (E34's decision 6).");
        var pooled = fit.Draws.SelectMany(c => c).ToArray();
        var (selected, weights) = Draws3.Select(options.Draws, pooled.Length, fit.DrawWeights, nameof(options));
        var (inputs, outside) = cases.Prepare(model.Transform, model.Inputs, options.UnseenLevels);
        var surface = new BassSurface(model.Basis, fit.Terms, inputs, null);
        var components = new MatrixTPredictive[selected.Length];
        var designs = new Dictionary<string, (BassDesign Design, BassStructureFit Fit, double[] Means)>();
        int n = surface.Cases, k = data.Response.Outputs;
        for (int d = 0; d < selected.Length; d++)
        {
            BassDraw draw = pooled[selected[d]]; var ids = draw.Structure.TermIds;
            string key = string.Join(',', ids);
            if (!designs.TryGetValue(key, out var cached))
            {
                var design = model.NewDesign();
                if (design.Load(ids.Select(id => fit.Terms[id]).ToArray()) is { } refusal) throw new ArgumentException($"The retained structure was refused ({refusal}).", nameof(fit));
                var means = new double[ids.Length]; var training = new double[data.Response.Cases];
                for (int j = 0; j < ids.Length; j++) { model.Basis.CaseColumn(fit.Terms[ids[j]], model.Inputs.Cases, training); means[j] = model.Working.ColumnMean(training, null); }
                cached = (design, BassStructureFit.Of(design, 0), means); designs.Add(key, cached);
            }
            bool flat = model.Target.Intercept == GPriorCentre.Flat;
            double precision = draw.Tau / model.Target.CoefficientConstant;
            int offset = flat ? 0 : 1, p = ids.Length + offset;
            var z = new double[n, p]; var mean = new double[n, k];
            if (!flat) for (int i = 0; i < n; i++) z[i, 0] = 1;
            for (int j = 0; j < ids.Length; j++) { var column = surface.Columns(ids[j]).Case; for (int i = 0; i < n; i++) z[i, j + offset] = column[i] - cached.Means[j]; }
            ISpdFactor factor; double factorScale; double[,] beta;
            if (cached.Fit.Outputs!.Prior is { } pencil)
            { (factor, beta) = pencil.Conditional(precision, 1); factorScale = 1; }
            else
            {
                factor = cached.Design.Factor; factorScale = 1 / (1 + precision); beta = new double[p, k];
                for (int l = 0; l < k; l++)
                {
                    if (!flat) beta[0, l] = model.Resolved.Outputs!.Levels[l] + (model.Resolved.Outputs.Means[l] - model.Resolved.Outputs.Levels[l]) / (1 + precision);
                    for (int j = 0; j < ids.Length; j++) beta[j + offset, l] = cached.Design.OutputCoefficients![j, l] / (1 + precision);
                }
            }
            for (int i = 0; i < n; i++) for (int l = 0; l < k; l++)
            { mean[i, l] = flat ? model.Resolved.Outputs!.Means[l] : 0; for (int j = 0; j < p; j++) mean[i, l] += z[i, j] * beta[j, l]; }
            bool g = cached.Fit.Outputs.Prior is null;
            var solved = new double[n][];
            for (int i = 0; i < n; i++)
            { var row = new double[factor.Dimension]; for (int j = 0; j < row.Length; j++) row[j] = z[i, j + (g ? offset : 0)]; solved[i] = factor.Solve(row); }
            var v = new double[n, n];
            for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++)
            {
                double entry = i == j ? 1 : 0;
                if (flat || g) entry += 1.0 / (data.Response.Cases * (flat ? 1 : 1 + precision));
                for (int a = 0; a < factor.Dimension; a++) entry += factorScale * z[i, a + (g ? offset : 0)] * solved[j][a];
                v[i, j] = v[j, i] = entry;
            }
            InverseWishart prior = model.Resolved.Outputs!.Rates is { } rates
                ? rates.PriorAt((draw.OutputRates ?? throw new ArgumentException("A Huang–Wand draw must carry its covariance rates.", nameof(fit))).Select(x => Math.Log(x)).ToArray())
                : model.Resolved.Outputs.Prior;
            var posterior = model.Target.OutputsMarginal(cached.Fit, draw.Tau).SigmaPosterior(prior);
            components[d] = new MatrixTPredictive(mean, v, posterior.Scale, posterior.DegreesOfFreedom);
        }
        return new BassMatrixTPrediction(components, weights, selected, outside);
    }

    public static SobolAnalysis Sensitivity(BassFit fit, BassData data, IReadOnlyList<double> output, SensitivitySpec? spec = null)
    { Model(fit, data); return BassSensitivity.Projected(fit, data, Functional(data.Response.Outputs, output), spec); }

    public static CoActiveAnalysis CoActive(BassFit fit, BassData data, CoActiveSpec spec)
    { Model(fit, data); return BassActiveSubspaces.MultiOutput(fit, data, spec); }

    public static BassGradient Gradient(BassFit fit, BassData data, PredictionCases cases, IReadOnlyList<int>? selection = null)
    {
        ArgumentNullException.ThrowIfNull(cases);
        var model = Model(fit, data); var pooled = fit.Draws.SelectMany(c => c).ToArray();
        var (selected, weights) = Draws3.Select(selection, pooled.Length, fit.DrawWeights, nameof(selection));
        var (inputs, outside) = cases.Prepare(model.Transform, model.Inputs, UnseenLevel.Refuse);
        var surface = new BassSurface(model.Basis, fit.Terms, inputs, null);
        var continuous = model.Transform.Predictors.OfType<RangeScale>().ToArray();
        var partials = new double[continuous.Length][,,]; var buffer = new double[surface.Cases];
        for (int v = 0; v < continuous.Length; v++)
        {
            var values = new double[selected.Length, surface.Cases, data.Response.Outputs];
            for (int d = 0; d < selected.Length; d++) for (int k = 0; k < data.Response.Outputs; k++)
            {
                var draw = pooled[selected[d]]; var coefficients = new double[draw.Structure.TermIds.Length + 1];
                for (int j = 0; j < coefficients.Length; j++) coefficients[j] = draw.OutputCoefficients![j, k];
                surface.EvaluateDerivative(draw.Structure.TermIds, coefficients, v, buffer);
                for (int i = 0; i < surface.Cases; i++) values[d, i, k] = buffer[i] / (continuous[v].Max - continuous[v].Min);
            }
            partials[v] = values;
        }
        return new BassGradient(continuous.Select(x => x.Name).ToArray(), partials, weights, outside);
    }

    /// <summary>Noise-covariance and correlation summaries; correlation of E[Σ] and E[correlation(Σ)] remain separate.</summary>
    public static BassOutputCovarianceSummary Covariance(BassFit fit, BassData data, IReadOnlyList<double>? probabilities = null)
    {
        var model = Model(fit, data); var pooled = fit.Draws.SelectMany(c => c).ToArray();
        var (_, weights) = Draws3.Select(null, pooled.Length, fit.DrawWeights, nameof(fit));
        int k = data.Response.Outputs;
        var entries = new PosteriorSummary[k, k]; var correlations = new PosteriorSummary[k, k];
        double nu = model.Resolved.Outputs!.Prior.DegreesOfFreedom + data.Response.Cases - (model.Target.Intercept == GPriorCentre.Flat ? 1 : 0);
        var mean = nu > k + 1 ? new double[k, k] : null; var correlationOfMean = mean is null ? null : new double[k, k];
        for (int a = 0; a < k; a++) for (int b = 0; b < k; b++)
        {
            double[] values = [.. pooled.Select(d => d.Sigma![a, b])];
            double[] rho = [.. pooled.Select(d => d.Sigma![a, b] / Math.Sqrt(d.Sigma[a, a] * d.Sigma[b, b]))];
            entries[a, b] = PosteriorSummary.Of(values, weights, probabilities);
            correlations[a, b] = PosteriorSummary.Of(rho, weights, probabilities);
            if (mean is not null) mean[a, b] = entries[a, b].Mean;
        }
        if (mean is not null) for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) correlationOfMean![a, b] = mean[a, b] / Math.Sqrt(mean[a, a] * mean[b, b]);
        return new BassOutputCovarianceSummary(((MultiOutputResponse)data.Response).Names, mean, correlationOfMean, entries, correlations);
    }
}

public sealed record BassOutputCovarianceSummary(IReadOnlyList<string> Outputs, double[,]? Mean, double[,]? CorrelationOfMean,
                                                 PosteriorSummary[,] CovarianceEntries, PosteriorSummary[,] CorrelationEntries);
