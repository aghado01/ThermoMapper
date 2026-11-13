using System;
using System.Collections.Generic;
using Maths.LinAlg;

namespace Maths.Regression.Conformal;

/// <summary>
/// The likelihood ratio w(x) = dP̃_X/dP_X between a target law of the inputs and the source law, estimated by telling their samples apart
/// (the conformal layer's slice 14 follow-up): a logistic regression of "target" on the features, its odds corrected for the two samples'
/// sizes, w(x) = odds(x)·n_source/n_target (the classifier's probability ratio; Tibshirani, Barber, Candès and Ramdas 2019, §3). The features
/// are centred and scaled on the pooled samples, so the ridge on the slopes, ridge·N/2·‖β‖² over the N pooled points, reads alike whatever
/// their units; it keeps the fit finite where the samples separate and costs nothing otherwise. Fitted by Newton's method with step halving
/// on the penalized log likelihood. The logistic model is right exactly when log w is linear in the features: give log x for a ratio
/// proportional to a power of x, x for a Gaussian mean shift.
/// <para>The weighted rule takes the estimate in place of the true ratio and loses at most the expected total variation between the data
/// set's normalized true and estimated weights (<c>sum_weight_refused_estimated_le</c>): the estimate's quality is the coverage it costs,
/// which no estimate here certifies.</para>
/// </summary>
public sealed class DensityRatio
{
    private readonly double[] _beta, _mean, _scale;
    private readonly double _offset;

    private DensityRatio(double[] beta, double[] mean, double[] scale, double offset, int iterations)
    {
        _beta = beta;
        _mean = mean;
        _scale = scale;
        _offset = offset;
        Iterations = iterations;
    }

    /// <summary>Fits the ratio from <paramref name="source"/> (rows the source law's inputs' features) and <paramref name="target"/> (the
    /// target law's), each at least one row and the same width.</summary>
    /// <param name="ridge">The slopes' penalty per pooled point on the standardized features; positive.</param>
    /// <param name="maxIterations">Newton steps allowed before the fit is refused as not converged.</param>
    /// <param name="tolerance">The largest standardized gradient entry per pooled point at which the fit stops.</param>
    public static DensityRatio Logistic(double[,] source, double[,] target, double ridge = 1e-4, int maxIterations = 100, double tolerance = 1e-10)
    {
        ArgumentNullException.ThrowIfNull(source);
        ArgumentNullException.ThrowIfNull(target);
        int n0 = source.GetLength(0), n1 = target.GetLength(0), p = source.GetLength(1);
        if (n0 == 0 || n1 == 0) throw new ArgumentException("Each sample needs at least one row.", n0 == 0 ? nameof(source) : nameof(target));
        if (target.GetLength(1) != p) throw new ArgumentException("The two samples must have the same features.", nameof(target));
        if (!(ridge > 0.0 && double.IsFinite(ridge))) throw new ArgumentOutOfRangeException(nameof(ridge), "The ridge is positive and finite.");
        ArgumentOutOfRangeException.ThrowIfLessThan(maxIterations, 1);
        int n = n0 + n1, d = p + 1;
        var mean = new double[p];
        var scale = new double[p];
        for (int j = 0; j < p; j++)
        {
            double sum = 0.0;
            for (int i = 0; i < n; i++) sum += Feature(source, target, n0, i, j);
            mean[j] = sum / n;
            double squares = 0.0;
            for (int i = 0; i < n; i++) { double e = Feature(source, target, n0, i, j) - mean[j]; squares += e * e; }
            scale[j] = Math.Sqrt(squares / n);
            if (!(scale[j] > 0.0 && double.IsFinite(scale[j])))
                throw new ArgumentException($"Feature {j} is constant or not finite across the samples, so it tells nothing apart.", nameof(source));
        }
        // The standardized design, an intercept column first; the labels 0 for the source, 1 for the target.
        var x = new double[n, d];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = 1.0;
            for (int j = 0; j < p; j++) x[i, j + 1] = (Feature(source, target, n0, i, j) - mean[j]) / scale[j];
        }
        double penalty = ridge * n;
        var beta = new double[d];
        beta[0] = Math.Log((double)n1 / n0);
        double Objective(double[] b)
        {
            double l = 0.0;
            for (int i = 0; i < n; i++)
            {
                double eta = Eta(x, i, b);
                l += (i >= n0 ? eta : 0.0) - Softplus(eta);
            }
            for (int j = 1; j < d; j++) l -= 0.5 * penalty * b[j] * b[j];
            return l;
        }
        var factor = new DenseCholesky(d);
        var gradient = new double[d];
        var hessian = new double[d, d];
        double current = Objective(beta);
        for (int iteration = 1; iteration <= maxIterations; iteration++)
        {
            Array.Clear(gradient);
            Array.Clear(hessian);
            for (int i = 0; i < n; i++)
            {
                double eta = Eta(x, i, beta), prob = Logistic(eta), w = prob * (1.0 - prob), residual = (i >= n0 ? 1.0 : 0.0) - prob;
                for (int a = 0; a < d; a++)
                {
                    gradient[a] += x[i, a] * residual;
                    for (int c = 0; c <= a; c++) hessian[a, c] += w * x[i, a] * x[i, c];
                }
            }
            double largest = 0.0;
            for (int j = 0; j < d; j++)
            {
                if (j > 0)
                {
                    gradient[j] -= penalty * beta[j];
                    hessian[j, j] += penalty;
                }
                largest = Math.Max(largest, Math.Abs(gradient[j]));
            }
            if (largest <= tolerance * n) return new DensityRatio(beta, mean, scale, Math.Log((double)n0 / n1), iteration - 1);
            factor.Decompose(hessian);
            if (factor.FailedAt >= 0) throw new InvalidOperationException($"The logistic fit's Hessian is not positive definite at step {iteration}.");
            double[] step = factor.Solve(gradient);
            // A Newton step below rounding on the standardized scale is the optimum: its gain is below the objective's own rounding, which
            // the halving would read as a loss.
            double stepSize = 0.0;
            foreach (double v in step) stepSize = Math.Max(stepSize, Math.Abs(v));
            if (stepSize <= 1e-12) return new DensityRatio(beta, mean, scale, Math.Log((double)n0 / n1), iteration - 1);
            double size = 1.0, slack = 1e-12 * (1.0 + Math.Abs(current));
            for (int halving = 0; halving < 50; halving++)
            {
                var trial = new double[d];
                for (int j = 0; j < d; j++) trial[j] = beta[j] + size * step[j];
                double value = Objective(trial);
                if (value >= current - slack)
                {
                    beta = trial;
                    current = value;
                    break;
                }
                size *= 0.5;
            }
        }
        throw new InvalidOperationException($"The logistic fit did not converge in {maxIterations} steps.");
    }

    /// <summary>The Newton steps the fit took.</summary>
    public int Iterations { get; }

    /// <summary>ln w(x) at <paramref name="features"/>: the fitted log odds plus ln(n_source/n_target).</summary>
    public double LogRatio(ReadOnlySpan<double> features)
    {
        if (features.Length != _mean.Length) throw new ArgumentException($"{_mean.Length} features are required.", nameof(features));
        double eta = _beta[0];
        for (int j = 0; j < _mean.Length; j++) eta += _beta[j + 1] * (features[j] - _mean[j]) / _scale[j];
        return eta + _offset;
    }

    /// <summary>w(x) at <paramref name="features"/>.</summary>
    public double Ratio(ReadOnlySpan<double> features) => Math.Exp(LogRatio(features));

    /// <summary>w at each row of <paramref name="features"/>.</summary>
    public double[] Ratios(double[,] features)
    {
        ArgumentNullException.ThrowIfNull(features);
        var row = new double[features.GetLength(1)];
        var ratios = new double[features.GetLength(0)];
        for (int i = 0; i < ratios.Length; i++)
        {
            for (int j = 0; j < row.Length; j++) row[j] = features[i, j];
            ratios[i] = Ratio(row);
        }
        return ratios;
    }

    /// <summary>The weighted rule's weights (<see cref="ShiftWeights"/>) from the estimate: each unit's at its features, and each test case's,
    /// spread over its grid where <paramref name="grids"/> is given (the density forms) and one per case otherwise (the closed forms).</summary>
    public ShiftWeights Weights(double[,] units, double[,] cases, IReadOnlyList<double[]>? grids = null)
    {
        double[] unitWeights = Ratios(units), caseWeights = Ratios(cases);
        return grids is null ? new ShiftWeights(unitWeights, caseWeights) : ShiftWeights.ByCase(unitWeights, caseWeights, grids);
    }

    private static double Feature(double[,] source, double[,] target, int n0, int i, int j) => i < n0 ? source[i, j] : target[i - n0, j];

    private static double Eta(double[,] x, int i, double[] b)
    {
        double eta = 0.0;
        for (int j = 0; j < b.Length; j++) eta += x[i, j] * b[j];
        return eta;
    }

    private static double Logistic(double eta) => eta >= 0.0 ? 1.0 / (1.0 + Math.Exp(-eta)) : Math.Exp(eta) / (1.0 + Math.Exp(eta));

    // ln(1 + e^η) without overflow.
    private static double Softplus(double eta) => eta > 0.0 ? eta + Math.Log(1.0 + Math.Exp(-eta)) : Math.Log(1.0 + Math.Exp(eta));
}
