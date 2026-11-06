using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// What a basis fit's predictive draws add to the latent output, R's <c>nugget</c> and <c>trunc.error</c>. Neither is
/// <see cref="BassBasisPrediction.Latent"/>, and R's defaults are <see cref="Noise"/>.
/// </summary>
public enum BasisPredictiveParts
{
    /// <summary>Component noise and a whole truncation residual: D07's full draw.</summary>
    NoiseAndTruncation,

    /// <summary>Component noise alone.</summary>
    Noise,

    /// <summary>A whole truncation residual alone.</summary>
    Truncation,
}

/// <summary>
/// What a prediction call asks for (D07): which draws, in which order; whether predictive draws join the latent ones;
/// the call's own seed, so the same seed reproduces the same draws; what an unseen level does; new output locations for
/// a direct functional fit; and the Rao–Blackwellized mean surface, off by default.
/// </summary>
public sealed record BassPredictOptions
{
    /// <summary>The draws to predict with, in the order wanted, as indices into the fit's draws pooled in chain order (a
    /// basis fit's joint draws); repeats allowed. Null takes every draw in order.</summary>
    public IReadOnlyList<int>? Draws { get; init; }

    /// <summary>Whether predictive draws, with their noise, join the latent ones.</summary>
    public bool Predictive { get; init; } = true;

    /// <summary>What a basis fit's predictive draws and covariance add to the latent output: both sources by default (D07).
    /// A direct fit, which has noise alone, takes only the default.</summary>
    public BasisPredictiveParts BasisParts { get; init; } = BasisPredictiveParts.NoiseAndTruncation;

    /// <summary>The call's seed for the noise (and a basis fit's truncation draws).</summary>
    public int Seed { get; init; }

    /// <summary>What a categorical level no training case took does: refused by default.</summary>
    public UnseenLevel UnseenLevels { get; init; } = UnseenLevel.Refuse;

    /// <summary>A direct functional fit's new output locations, one column per output variable, matched by name; null
    /// keeps the training locations. A scalar fit takes none, and a basis exists only at its training locations.</summary>
    public IReadOnlyList<ContinuousColumn>? OutputLocations { get; init; }

    /// <summary>Whether the Rao–Blackwellized mean surface rides beside the draws: each draw's E[f | structure, τ, y],
    /// averaged under the weights. A point estimate, never a draw.</summary>
    public bool MeanSurface { get; init; }
}

/// <summary>
/// A direct fit's prediction (D07, P07): latent draws f(x), or f(x, t) for a direct functional fit, and predictive draws
/// that add N(0, σ²) noise drawn independently for each draw, case and location, as draw × case × location arrays with
/// every axis kept, singletons included, in the requested draw order. A scalar fit has one location. Under a noise law
/// (<see cref="IBassObservationModel.Law"/>, E1) the three readings are the law's: the surface f, the mean f plus
/// <see cref="MeanShifts"/>, and predictive draws whose noise draws each observation's latent from its prior, then ζ.
/// <list type="bullet">
/// <item><see cref="CaseOutsideRange"/> flags each case outside the training range of some continuous predictor, and
/// <see cref="LocationOutsideRange"/> each location outside the output variables' training ranges; the hinges extrapolate
/// linearly there.</item>
/// <item><see cref="Weights"/> are the requested draws' weights, normalized: equal for MCMC, SMC's when it joins.</item>
/// <item><see cref="Quantile"/> is the predictive mixture's quantile over the weighted draws, not a quantile of each
/// draw's own interval.</item>
/// </list>
/// </summary>
public sealed class BassPrediction
{
    private readonly IContinuousDistribution[]? _noises;   // each draw's noise under a law; null for the Gaussian

    internal BassPrediction(int[] draws, double[] weights, double[] noiseVariances, double[,,] latent, double[,,]? predictive,
                            bool[] caseOutside, bool[] locationOutside, double[,]? meanSurface, IContinuousDistribution[]? noises = null,
                            double[]? meanShifts = null)
    {
        _noises = noises;
        MeanShifts = Array.AsReadOnly(meanShifts ?? new double[draws.Length]);
        Draws = Array.AsReadOnly(draws);
        Weights = Array.AsReadOnly(weights);
        NoiseVariances = Array.AsReadOnly(noiseVariances);
        Latent = latent;
        Predictive = predictive;
        CaseOutsideRange = Array.AsReadOnly(caseOutside);
        LocationOutsideRange = Array.AsReadOnly(locationOutside);
        MeanSurface = meanSurface;
    }

    /// <summary>The draws used, in order, as indices into the fit's pooled draws.</summary>
    public IReadOnlyList<int> Draws { get; }

    /// <summary>Each draw's weight, normalized to sum to 1.</summary>
    public IReadOnlyList<double> Weights { get; }

    /// <summary>Each draw's noise variance σ², the core's: under a noise law its scale, the noise's variance being the law's
    /// σ²(β²Var[v] + c·E[v]).</summary>
    public IReadOnlyList<double> NoiseVariances { get; }

    /// <summary>Each draw's E[ε], σβE[v] under a skewed law: the mean of y at x is f plus this. Zero for a symmetric noise.</summary>
    public IReadOnlyList<double> MeanShifts { get; }

    /// <summary>f under each draw: draw × case × location.</summary>
    public double[,,] Latent { get; }

    /// <summary>f plus noise, independent for each draw, case and location; null when not asked for.</summary>
    public double[,,]? Predictive { get; }

    /// <summary>Whether each case lies outside the training range of some continuous predictor.</summary>
    public IReadOnlyList<bool> CaseOutsideRange { get; }

    /// <summary>Whether each location lies outside an output variable's training range; all false at the training
    /// locations and for a scalar fit.</summary>
    public IReadOnlyList<bool> LocationOutsideRange { get; }

    /// <summary>The Rao–Blackwellized mean surface, case × location; null when not asked for.</summary>
    public double[,]? MeanSurface { get; }

    /// <summary>The draws.</summary>
    public int DrawCount => Latent.GetLength(0);

    /// <summary>The cases.</summary>
    public int Cases => Latent.GetLength(1);

    /// <summary>The locations.</summary>
    public int Locations => Latent.GetLength(2);

    /// <summary>The posterior mean of f, the weighted mean of the latent draws, case × location.</summary>
    public double[,] Mean() => Draws3.WeightedMean(Latent, Weights);

    /// <summary>The posterior mean of y's mean at each case and location, f plus each draw's <see cref="MeanShifts"/>: f's under a
    /// symmetric noise.</summary>
    public double[,] ObservationMean()
    {
        double[,] mean = Mean();
        double shift = 0.0;
        for (int d = 0; d < DrawCount; d++) shift += Weights[d] * MeanShifts[d];
        for (int i = 0; i < Cases; i++)
            for (int l = 0; l < Locations; l++) mean[i, l] += shift;
        return mean;
    }

    /// <summary>
    /// The <paramref name="probability"/> quantile of the predictive mixture at each case and location, case × location:
    /// Σₛ wₛ·N(fₛ, σₛ²), solved exactly (<see cref="PredictiveMixture.Quantile(ReadOnlySpan{double}, ReadOnlySpan{double}, ReadOnlySpan{double}, double)"/>);
    /// under a noise law, Σₛ wₛ·Fₛ(· − fₛ) with each draw's noise, its latent integrated out.
    /// </summary>
    public double[,] Quantile(double probability)
    {
        int s = DrawCount, n = Cases, m = Locations;
        if (_noises is { } noises)
        {
            var w = new double[s];
            for (int d = 0; d < s; d++) w[d] = Weights[d];
            var at = new double[s];
            var lawQ = new double[n, m];
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++)
                {
                    for (int d = 0; d < s; d++) at[d] = Latent[d, i, l];
                    lawQ[i, l] = PredictiveMixture.Quantile(at, noises, w, probability);
                }
            return lawQ;
        }
        var sds = new double[s];
        var weights = new double[s];
        for (int d = 0; d < s; d++)
        {
            sds[d] = Math.Sqrt(NoiseVariances[d]);
            weights[d] = Weights[d];
        }
        var means = new double[s];
        var q = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++)
            {
                for (int d = 0; d < s; d++) means[d] = Latent[d, i, l];
                q[i, l] = PredictiveMixture.Quantile(means, sds, weights, probability);
            }
        return q;
    }
}

/// <summary>Weighted reductions over a draw × case × location array's first axis.</summary>
internal static class Draws3
{
    public static double[,] WeightedMean(double[,,] values, IReadOnlyList<double> weights)
    {
        int s = values.GetLength(0), n = values.GetLength(1), m = values.GetLength(2);
        var mean = new double[n, m];
        for (int d = 0; d < s; d++)
        {
            double w = weights[d];
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++) mean[i, l] += w * values[d, i, l];
        }
        return mean;
    }

    // The requested draws, checked against the pool, and their weights normalized over the request; equal when the fit
    // carries none. Refusals name the caller's parameter.
    public static (int[] Draws, double[] Weights) Select(IReadOnlyList<int>? requested, int pooled, IReadOnlyList<double>? weights, string name = "options")
    {
        int[] draws;
        if (requested is null)
        {
            draws = new int[pooled];
            for (int d = 0; d < pooled; d++) draws[d] = d;
        }
        else
        {
            if (requested.Count == 0) throw new ArgumentException("Ask for at least one draw.", name);
            draws = new int[requested.Count];
            for (int d = 0; d < draws.Length; d++)
            {
                if ((uint)requested[d] >= (uint)pooled)
                    throw new ArgumentException($"Draw {requested[d]} is not among the fit's {pooled}.", name);
                draws[d] = requested[d];
            }
        }
        if (weights is not null && weights.Count != pooled)
            throw new ArgumentException($"The fit's {weights.Count} draw weights do not match its {pooled} draws.", name);
        var selected = new double[draws.Length];
        double total = 0.0;
        for (int d = 0; d < draws.Length; d++)
        {
            double w = weights is null ? 1.0 : weights[draws[d]];
            if (!(w >= 0.0) || double.IsPositiveInfinity(w)) throw new ArgumentException("The draw weights must be finite and not negative.", name);
            selected[d] = w;
            total += w;
        }
        if (!(total > 0.0)) throw new ArgumentException("The requested draws carry no weight.", name);
        for (int d = 0; d < selected.Length; d++) selected[d] /= total;
        return (draws, selected);
    }
}

/// <summary>
/// Direct prediction (D07; <see cref="BassFit.Predict"/>): the requested draws evaluated at the cases (and locations)
/// through the training basis, the noise drawn from the call's seed in draw, case, location order, and the
/// Rao–Blackwellized mean surface on request.
/// </summary>
internal static class BassDirectPrediction
{
    public static BassPrediction Predict(BassModel model, BassTermTable terms, IReadOnlyList<IReadOnlyList<BassDraw>> draws,
                                         IReadOnlyList<double>? drawWeights, PredictionCases cases, BassPredictOptions options)
    {
        ArgumentNullException.ThrowIfNull(cases);
        ArgumentNullException.ThrowIfNull(options);
        if (options.BasisParts != BasisPredictiveParts.NoiseAndTruncation)
            throw new ArgumentException("A direct fit has no truncation residuals: its predictive draws add noise, and BasisParts is a basis fit's.", nameof(options));
        (List<BassDraw> pooled, _) = BassDrawFit.Pool(draws, terms, options.MeanSurface, nameof(draws));
        (int[] selected, double[] weights) = Draws3.Select(options.Draws, pooled.Count, drawWeights);
        (CaseInputs caseInputs, bool[] caseOutside) = cases.Prepare(model.Transform, model.Inputs, options.UnseenLevels);

        OutputInputs? locations = null;
        bool[] locationOutside;
        if (model.Data.Response is FunctionalResponse functional)
        {
            IReadOnlyList<ContinuousColumn> given = options.OutputLocations
                ?? Array.ConvertAll([.. functional.OutputVariables], v => new ContinuousColumn(v.Name, v.Values.ToArray()));
            (locations, locationOutside) = PredictionCases.PrepareLocations(model.Transform, given);
        }
        else
        {
            if (options.OutputLocations is not null) throw new ArgumentException("A scalar fit has no output locations.", nameof(options));
            locationOutside = [false];
        }

        var surface = new BassSurface(model.Basis, terms, caseInputs, locations);
        int s = selected.Length, n = surface.Cases, m = surface.Locations;
        var latent = new double[s, n, m];
        var noise = new double[s];
        var f = new double[n * m];
        IBassNoiseLaw? law = model.Observation.Law;
        IContinuousDistribution[]? noises = law is null ? null : new IContinuousDistribution[s];
        double[]? shifts = law is null ? null : new double[s];
        for (int d = 0; d < s; d++)
        {
            BassDraw draw = pooled[selected[d]];
            noise[d] = draw.Sigma2;
            if (law is not null)
            {
                ObservationLatents? globals = draw.LatentsOn(model);
                noises![d] = law.At(draw.Sigma2, globals);
                shifts![d] = law.MeanShift(draw.Sigma2, globals);
            }
            surface.Evaluate(draw.Structure.TermIds, draw.Coefficients, f);
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++) latent[d, i, l] = f[i * m + l];
        }

        double[,,]? predictive = null;
        if (options.Predictive)
        {
            var rng = new Xoshiro256PlusPlus(options.Seed);
            predictive = new double[s, n, m];
            for (int d = 0; d < s; d++)
            {
                if (law is not null)
                {
                    ObservationLatents? globals = pooled[selected[d]].LatentsOn(model);
                    for (int i = 0; i < n; i++)
                        for (int l = 0; l < m; l++) predictive[d, i, l] = latent[d, i, l] + law.Sample(noise[d], globals, rng);
                    continue;
                }
                double sd = Math.Sqrt(noise[d]);
                for (int i = 0; i < n; i++)
                    for (int l = 0; l < m; l++) predictive[d, i, l] = latent[d, i, l] + sd * Normal.SampleStandard(rng);
            }
        }

        double[,]? meanSurface = null;
        if (options.MeanSurface)
        {
            var means = new BassConditionalMeans(model, terms);
            meanSurface = new double[n, m];
            for (int d = 0; d < s; d++)
            {
                BassDraw draw = pooled[selected[d]];
                surface.Evaluate(draw.Structure.TermIds, means.Coefficients(draw), f);
                for (int i = 0; i < n; i++)
                    for (int l = 0; l < m; l++) meanSurface[i, l] += weights[d] * f[i * m + l];
            }
        }
        return new BassPrediction(selected, weights, noise, latent, predictive, caseOutside, locationOutside, meanSurface, noises, shifts);
    }
}
