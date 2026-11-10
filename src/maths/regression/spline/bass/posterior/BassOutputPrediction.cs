using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bass;

/// <summary>A simultaneous band calibrated by the weighted maximum standardized deviation of the joint draws.</summary>
public sealed record BassOutputBands(double[,] Lower, double[,] Upper, double Probability, bool Predictive, bool AcrossCases);

/// <summary>Several outputs' joint latent and noisy prediction draws, preserving their posterior weights and within-case covariance.</summary>
public sealed class BassOutputPrediction
{
    internal BassOutputPrediction(IReadOnlyList<string> outputs, int[] draws, double[] weights, double[,,] latent, double[,,]? predictive,
                                  double[][,] sigma, bool[] outside, double[,]? meanSurface, bool finiteCovariance)
    { Outputs = outputs; Draws = draws; Weights = weights; Latent = latent; Predictive = predictive; Sigma = sigma; CaseOutsideRange = outside; MeanSurface = meanSurface; HasFiniteCovariance = finiteCovariance; }
    public IReadOnlyList<string> Outputs { get; }
    public IReadOnlyList<int> Draws { get; }
    public IReadOnlyList<double> Weights { get; }
    public double[,,] Latent { get; }
    public double[,,]? Predictive { get; }
    public IReadOnlyList<double[,]> Sigma { get; }
    public IReadOnlyList<bool> CaseOutsideRange { get; }
    public double[,]? MeanSurface { get; }
    public int Cases => Latent.GetLength(1);
    public bool HasFiniteCovariance { get; }

    public double[,] Mean()
    {
        var mean = new double[Cases, Outputs.Count];
        for (int d = 0; d < Draws.Count; d++) for (int i = 0; i < Cases; i++) for (int k = 0; k < Outputs.Count; k++) mean[i, k] += Weights[d] * Latent[d, i, k];
        return mean;
    }

    /// <summary>Cov(f) plus E[Σ] for a noisy prediction, or Cov(f) alone; correlations follow from this joint covariance.</summary>
    public double[,] Covariance(int @case, bool predictive = true)
    {
        if (!HasFiniteCovariance) throw new InvalidOperationException("The posterior has no finite covariance; use the predictive draws or quantiles.");
        return Covariance(@case, predictive, Mean());
    }

    private double[,] Covariance(int @case, bool predictive, double[,] mean)
    {
        int k = Outputs.Count; var covariance = new double[k, k];
        for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) for (int d = 0; d < Draws.Count; d++)
            covariance[a, b] += Weights[d] * ((Latent[d, @case, a] - mean[@case, a]) * (Latent[d, @case, b] - mean[@case, b]) + (predictive ? Sigma[d][a, b] : 0));
        return covariance;
    }

    public double[,] Correlation(int @case, bool predictive = true)
    {
        var covariance = Covariance(@case, predictive); int k = Outputs.Count; var correlation = new double[k, k];
        for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) correlation[a, b] = covariance[a, b] / Math.Sqrt(covariance[a, a] * covariance[b, b]);
        return correlation;
    }

    /// <summary>Marginal predictive quantiles of Σ_d w_d N(f_dk,Σ_dkk), with observation noise integrated analytically.</summary>
    public double[,] Quantile(double probability)
    {
        int s = Draws.Count; var result = new double[Cases, Outputs.Count];
        var at = new double[s]; var sd = new double[s]; var weights = Weights.ToArray();
        for (int k = 0; k < Outputs.Count; k++)
        {
            for (int d = 0; d < s; d++) sd[d] = Math.Sqrt(Sigma[d][k, k]);
            for (int i = 0; i < Cases; i++) { for (int d = 0; d < s; d++) at[d] = Latent[d, i, k]; result[i, k] = PredictiveMixture.Quantile(at, sd, weights, probability); }
        }
        return result;
    }

    public BassOutputBands JointBands(double probability = 0.95, bool predictive = true, bool acrossCases = false)
    {
        if (!(probability > 0 && probability < 1)) throw new ArgumentOutOfRangeException(nameof(probability));
        var values = predictive ? Predictive ?? throw new InvalidOperationException("Joint noisy bands need predictive draws.") : Latent;
        var centre = Mean(); var sd = new double[Cases, Outputs.Count];
        for (int i = 0; i < Cases; i++) { var covariance = Covariance(i, predictive, centre); for (int k = 0; k < Outputs.Count; k++) sd[i, k] = Math.Sqrt(covariance[k, k]); }
        var lower = new double[Cases, Outputs.Count]; var upper = new double[Cases, Outputs.Count];
        for (int group = 0; group < (acrossCases ? 1 : Cases); group++)
        {
            var maxima = new DrawSample();
            for (int d = 0; d < Draws.Count; d++)
            {
                double maximum = 0;
                for (int i = acrossCases ? 0 : group; i < (acrossCases ? Cases : group + 1); i++) for (int k = 0; k < Outputs.Count; k++)
                    maximum = Math.Max(maximum, sd[i, k] > 0 ? Math.Abs(values[d, i, k] - centre[i, k]) / sd[i, k] : 0);
                maxima.Add(maximum, Weights[d]);
            }
            double q = maxima.Quantile(probability);
            for (int i = acrossCases ? 0 : group; i < (acrossCases ? Cases : group + 1); i++) for (int k = 0; k < Outputs.Count; k++)
            { lower[i, k] = centre[i, k] - q * sd[i, k]; upper[i, k] = centre[i, k] + q * sd[i, k]; }
        }
        return new(lower, upper, probability, predictive, acrossCases);
    }

    /// <summary>The weighted event frequency P(all outputs in their bounds), per case, over the joint noisy draws.</summary>
    public double[] SpecificationProbability(double[,] lower, double[,] upper)
    {
        ArgumentNullException.ThrowIfNull(lower); ArgumentNullException.ThrowIfNull(upper);
        var values = Predictive ?? throw new InvalidOperationException("Specification probabilities need predictive draws.");
        if (lower.GetLength(0) != Cases || upper.GetLength(0) != Cases || lower.GetLength(1) != Outputs.Count || upper.GetLength(1) != Outputs.Count)
            throw new ArgumentException("Bounds must be case by output.");
        for (int i = 0; i < Cases; i++) for (int k = 0; k < Outputs.Count; k++)
            if (double.IsNaN(lower[i, k]) || double.IsNaN(upper[i, k]) || lower[i, k] > upper[i, k]) throw new ArgumentException("Bounds must be ordered and not NaN.");
        var probability = new double[Cases];
        for (int i = 0; i < Cases; i++) for (int d = 0; d < Draws.Count; d++)
        { bool inside = true; for (int k = 0; k < Outputs.Count; k++) inside &= values[d, i, k] >= lower[i, k] && values[d, i, k] <= upper[i, k]; if (inside) probability[i] += Weights[d]; }
        return probability;
    }
}

/// <summary>The weighted matrix-t predictive conditional on retained structures, precisions and covariance rates.</summary>
public sealed record BassMatrixTPrediction(IReadOnlyList<MatrixTPredictive> Components, IReadOnlyList<double> Weights,
                                           IReadOnlyList<int> Draws, IReadOnlyList<bool> CaseOutsideRange)
{
    public double LogDensity(double[,] values)
    {
        double largest = double.NegativeInfinity; var logs = new double[Components.Count];
        for (int d = 0; d < logs.Length; d++) { logs[d] = Weights[d] == 0 ? double.NegativeInfinity : Math.Log(Weights[d]) + Components[d].LogDensity(values); largest = Math.Max(largest, logs[d]); }
        if (double.IsNegativeInfinity(largest)) return largest;
        double sum = 0; for (int d = 0; d < logs.Length; d++) sum += Math.Exp(logs[d] - largest);
        return largest + Math.Log(sum);
    }

    public double[,] Mean()
    {
        if (Components.Any(c => !c.HasMean)) throw new InvalidOperationException("A predictive component has no finite mean.");
        var result = new double[Components[0].Cases, Components[0].Outputs];
        for (int d = 0; d < Components.Count; d++) { var mean = Components[d].Location; for (int i = 0; i < result.GetLength(0); i++) for (int k = 0; k < result.GetLength(1); k++) result[i, k] += Weights[d] * mean[i, k]; }
        return result;
    }
}
