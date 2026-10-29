using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A basis fit's prediction (D07, P08), over the joint draws in the requested order, as draw × case × location arrays in
/// physical units with every axis kept:
/// <list type="bullet">
/// <item><b>Latent</b>: y_m + y_s ∘ K w(x), each component's coefficient surface under its draw (<see cref="Coefficients"/>).</item>
/// <item><b>Predictive</b>: y_m + y_s ∘ (K(w(x) + ε) + e), with ε_k ~ N(0, σ_k²) drawn for each draw and case (R draws one
/// per draw and recycles it across cases), and e a whole training residual vector, which keeps its correlation across
/// locations, drawn uniformly over the training cases independently for each draw and case, from the declared residual
/// set in standardized units, added before y_s.</item>
/// <item><see cref="Covariance"/>: the predictive mixture's covariance over the locations at one case,
/// diag(y_s)·(K(Cov[w] + E[diag σ²])Kᵀ + Cov[e])·diag(y_s), or the latent K·Cov[w]·Kᵀ under the same transform.</item>
/// <item><see cref="Quantile"/>: the weighted quantile of the predictive draws, the predictive mixture's quantile up to
/// Monte Carlo error.</item>
/// </list>
/// </summary>
public sealed class BassBasisPrediction
{
    private readonly LinearReconstruction _map;
    private readonly TruncationResiduals _truncation;

    internal BassBasisPrediction(LinearReconstruction map, TruncationResiduals truncation, int[] draws, double[] weights, double[,] noise,
                                 double[,,] coefficients, double[,,] latent, double[,,]? predictive, bool[] caseOutside, double[,]? meanSurface)
    {
        _map = map;
        _truncation = truncation;
        Draws = Array.AsReadOnly(draws);
        Weights = Array.AsReadOnly(weights);
        NoiseVariances = noise;
        Coefficients = coefficients;
        Latent = latent;
        Predictive = predictive;
        CaseOutsideRange = Array.AsReadOnly(caseOutside);
        MeanSurface = meanSurface;
    }

    /// <summary>The joint draws used, in order.</summary>
    public IReadOnlyList<int> Draws { get; }

    /// <summary>Each draw's weight, normalized: the joint draws weigh equally.</summary>
    public IReadOnlyList<double> Weights { get; }

    /// <summary>Each draw's component noise variances σ_k²: draw × component.</summary>
    public double[,] NoiseVariances { get; }

    /// <summary>The coefficient surfaces w_k(x): draw × case × component.</summary>
    public double[,,] Coefficients { get; }

    /// <summary>y_m + y_s ∘ K w(x): draw × case × location.</summary>
    public double[,,] Latent { get; }

    /// <summary>y_m + y_s ∘ (K(w(x) + ε) + e); null when not asked for.</summary>
    public double[,,]? Predictive { get; }

    /// <summary>Whether each case lies outside the training range of some continuous predictor.</summary>
    public IReadOnlyList<bool> CaseOutsideRange { get; }

    /// <summary>y_m + y_s ∘ K w̄(x), with w̄ each component's Rao–Blackwellized mean surface; null when not asked for.</summary>
    public double[,]? MeanSurface { get; }

    /// <summary>The draws.</summary>
    public int DrawCount => Latent.GetLength(0);

    /// <summary>The cases.</summary>
    public int Cases => Latent.GetLength(1);

    /// <summary>The locations.</summary>
    public int Locations => Latent.GetLength(2);

    /// <summary>The posterior mean of the latent output, case × location.</summary>
    public double[,] Mean() => Draws3.WeightedMean(Latent, Weights);

    /// <summary>The <paramref name="probability"/> quantile of the weighted predictive draws at each case and location,
    /// case × location; refused when the predictive draws were not asked for.</summary>
    public double[,] Quantile(double probability)
    {
        double[,,] draws = Predictive ?? throw new InvalidOperationException("The predictive draws were not asked for.");
        int s = DrawCount, n = Cases, m = Locations;
        var weights = new double[s];
        for (int d = 0; d < s; d++) weights[d] = Weights[d];
        var values = new double[s];
        var q = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++)
            {
                for (int d = 0; d < s; d++) values[d] = draws[d, i, l];
                q[i, l] = PredictiveMixture.Quantile(values, weights, probability);
            }
        return q;
    }

    /// <summary>
    /// The covariance over the locations at case <paramref name="case"/>, in physical units (m × m): the predictive
    /// mixture's, diag(y_s)·(K(Cov[w] + E[diag σ²])Kᵀ + Cov[e])·diag(y_s), with Cov[w] and E[σ²] under the draws' weights
    /// and Cov[e] under the uniform draw over the training residuals; or, when <paramref name="predictive"/> is false, the
    /// latent output's, diag(y_s)·K·Cov[w]·Kᵀ·diag(y_s).
    /// </summary>
    public double[,] Covariance(int @case, bool predictive = true)
    {
        if ((uint)@case >= (uint)Cases) throw new ArgumentOutOfRangeException(nameof(@case));
        int s = DrawCount, r = Coefficients.GetLength(2), m = Locations;
        var mean = new double[r];
        for (int d = 0; d < s; d++)
            for (int k = 0; k < r; k++) mean[k] += Weights[d] * Coefficients[d, @case, k];
        var inner = new double[r, r];   // Cov[w] (+ E[diag σ²])
        for (int d = 0; d < s; d++)
            for (int k = 0; k < r; k++)
            {
                double a = Coefficients[d, @case, k] - mean[k];
                for (int q = 0; q < r; q++) inner[k, q] += Weights[d] * a * (Coefficients[d, @case, q] - mean[q]);
                if (predictive) inner[k, k] += Weights[d] * NoiseVariances[d, k];
            }
        double[,] basis = _map.BasisArray;
        var kInner = new double[m, r];
        for (int l = 0; l < m; l++)
            for (int q = 0; q < r; q++)
            {
                double v = 0.0;
                for (int k = 0; k < r; k++) v += basis[l, k] * inner[k, q];
                kInner[l, q] = v;
            }
        double[,] standardized = predictive ? _truncation.Covariance() : new double[m, m];
        for (int l = 0; l < m; l++)
            for (int j = 0; j < m; j++)
            {
                double v = 0.0;
                for (int q = 0; q < r; q++) v += kInner[l, q] * basis[j, q];
                standardized[l, j] += v;
            }
        return _map.Physical(standardized);
    }
}

/// <summary>
/// Basis prediction (D07; <see cref="BassBasisFit.Predict"/>): each component's surface under its draw in each requested
/// joint draw, reconstructed; the predictive draws' component noise and truncation residuals from the call's seed, in
/// draw then case order, the r noises before the residual's case.
/// </summary>
internal static class BassBasisPredictor
{
    public static BassBasisPrediction Predict(BassBasisFit fit, PredictionCases cases, BassPredictOptions options)
    {
        ArgumentNullException.ThrowIfNull(cases);
        ArgumentNullException.ThrowIfNull(options);
        if (options.OutputLocations is not null)
            throw new ArgumentException("A basis exists only at its training locations: it has no evaluation at new output coordinates.", nameof(options));
        BassDecomposition decomposition = fit.Decomposition;
        LinearReconstruction map = decomposition.Map;
        BassJointDraws joint = fit.JointDraws;
        int r = decomposition.Rank, m = decomposition.Locations;
        (int[] selected, double[] weights) = Draws3.Select(options.Draws, joint.Count, null);
        BassModel first = fit.Model(0);
        (CaseInputs caseInputs, bool[] caseOutside) = cases.Prepare(first.Transform, first.Inputs, options.UnseenLevels);
        int s = selected.Length, n = caseInputs.Cases;

        var coefficients = new double[s, n, r];
        var noise = new double[s, r];
        double[,]? meanCoefficients = options.MeanSurface ? new double[n, r] : null;
        var f = new double[n];
        for (int k = 0; k < r; k++)
        {
            BassFit component = fit.Components[k];
            BassModel model = fit.Model(k);
            (List<BassDraw> pooled, _) = BassDrawFit.Pool(component.Draws, component.Terms, options.MeanSurface, "fit");
            var surface = new BassSurface(model.Basis, component.Terms, caseInputs, null);
            BassConditionalMeans? means = options.MeanSurface ? new BassConditionalMeans(model, component.Terms) : null;
            for (int d = 0; d < s; d++)
            {
                BassDraw draw = pooled[joint.Draw(k, selected[d])];
                noise[d, k] = draw.Sigma2;
                surface.Evaluate(draw.Structure.TermIds, draw.Coefficients, f);
                for (int i = 0; i < n; i++) coefficients[d, i, k] = f[i];
                if (means is null) continue;
                surface.Evaluate(draw.Structure.TermIds, means.Coefficients(draw), f);
                for (int i = 0; i < n; i++) meanCoefficients![i, k] += weights[d] * f[i];
            }
        }

        var latent = new double[s, n, m];
        var w = new double[r];
        var y = new double[m];
        for (int d = 0; d < s; d++)
            for (int i = 0; i < n; i++)
            {
                for (int k = 0; k < r; k++) w[k] = coefficients[d, i, k];
                map.Reconstruct(w, [], y);
                for (int l = 0; l < m; l++) latent[d, i, l] = y[l];
            }

        double[,,]? predictive = null;
        if (options.Predictive)
        {
            double[,] residuals = decomposition.Truncation.Values;
            int trainingCases = residuals.GetLength(0);
            var e = new double[m];
            var rng = new Xoshiro256PlusPlus(options.Seed);
            predictive = new double[s, n, m];
            for (int d = 0; d < s; d++)
                for (int i = 0; i < n; i++)
                {
                    for (int k = 0; k < r; k++) w[k] = coefficients[d, i, k] + Math.Sqrt(noise[d, k]) * Normal.SampleStandard(rng);
                    int drawn = rng.NextInt(trainingCases);
                    for (int l = 0; l < m; l++) e[l] = residuals[drawn, l];
                    map.Reconstruct(w, e, y);
                    for (int l = 0; l < m; l++) predictive[d, i, l] = y[l];
                }
        }

        double[,]? meanSurface = null;
        if (meanCoefficients is not null)
        {
            meanSurface = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                for (int k = 0; k < r; k++) w[k] = meanCoefficients[i, k];
                map.Reconstruct(w, [], y);
                for (int l = 0; l < m; l++) meanSurface[i, l] = y[l];
            }
        }
        return new BassBasisPrediction(map, decomposition.Truncation, selected, weights, noise, coefficients, latent, predictive, caseOutside, meanSurface);
    }
}
