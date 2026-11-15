using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline;

namespace Maths.Regression.Conformal;

/// <summary>
/// The likelihood ratio w(x) = dP̃_X/dP_X between a target law of the inputs and the source law, estimated by telling their samples apart
/// (the conformal layer's slice 14 follow-up): a logistic regression of "target" on the features, its odds corrected for the two samples'
/// sizes, w(x) = odds(x)·n_source/n_target (the classifier's probability ratio; Tibshirani, Barber, Candès and Ramdas 2019, §3). The features
/// are centred and scaled on the pooled samples, so the ridge on the slopes, ridge·N/2·‖β‖² over the N pooled points, reads alike whatever
/// their units; it keeps the fit finite where the samples separate and costs nothing otherwise. Fitted by Newton's method with step halving
/// on the penalized log likelihood. The logistic model is right exactly when log w is linear in the features: give log x for a ratio
/// proportional to a power of x, x for a Gaussian mean shift. Where that form is not known, <see cref="Spline"/> expands each input into
/// cubic B-splines first, so log w may be any smooth sum of one function per input.
/// <para>The weighted rule takes the estimate in place of the true ratio and loses at most the expected total variation between the data
/// set's normalized true and estimated weights (<c>sum_weight_refused_estimated_le</c>): the estimate's quality is the coverage it costs,
/// which no estimate here certifies.</para>
/// </summary>
public sealed class DensityRatio
{
    private readonly double[] _beta, _mean, _scale;
    private readonly double _offset;
    private readonly SplineExpansion? _expansion;   // each input's spline columns, where the fit took them (Spline)

    private DensityRatio(double[] beta, double[] mean, double[] scale, double offset, int iterations, double ridge, SplineExpansion? expansion = null)
    {
        _beta = beta;
        _mean = mean;
        _scale = scale;
        _offset = offset;
        Iterations = iterations;
        Ridge = ridge;
        _expansion = expansion;
    }

    /// <summary>
    /// The ratio with each input expanded into cubic B-splines (the ratio estimator's follow-up: a nonparametric ratio): an additive
    /// logistic model, ln w(x) = Σ_j s_j(x_j) up to the size correction, where <see cref="Logistic"/> is linear in the features given.
    /// Each input's B-splines run over the pooled samples' range, <paramref name="knots"/> interior knots at the pooled quantiles, the
    /// first of each input's columns dropped since they sum to one, which the intercept spans; an input beyond the range holds its edge
    /// value. The ridge on the standardized columns is chosen from <paramref name="ridges"/> by <paramref name="folds"/>-fold
    /// cross-validated log loss of the classifier, each row in fold i mod folds of its sample; a ridge whose fit is refused in some fold is
    /// passed over; a single ridge is taken as given. The fit then reads every row at that ridge (<see cref="Ridge"/>). A column constant
    /// over the pooled samples is left out (a few-valued input leaves some), and one constant within a fold's training rows sits idle there,
    /// its slope held at zero by the ridge. Interactions between inputs are not modelled.
    /// </summary>
    /// <param name="source">The source law's sample, one row of inputs.</param>
    /// <param name="target">The target law's sample, the same inputs.</param>
    /// <param name="knots">Interior knots per input, at least 1.</param>
    /// <param name="ridges">The ridges to choose from, each positive and finite; null takes 10⁻⁴ to 10 by decades.</param>
    /// <param name="folds">The cross-validation's folds, at least 2 and at most either sample's rows.</param>
    /// <param name="maxIterations">Newton steps allowed per fit.</param>
    /// <param name="tolerance">The largest standardized gradient entry per pooled point at which a fit stops.</param>
    public static DensityRatio Spline(double[,] source, double[,] target, int knots = 6, IReadOnlyList<double>? ridges = null, int folds = 5,
                                      int maxIterations = 100, double tolerance = 1e-10)
    {
        ArgumentNullException.ThrowIfNull(source);
        ArgumentNullException.ThrowIfNull(target);
        int n0 = source.GetLength(0), n1 = target.GetLength(0), p = source.GetLength(1);
        if (target.GetLength(1) != p) throw new ArgumentException("The two samples must have the same features.", nameof(target));
        ArgumentOutOfRangeException.ThrowIfLessThan(knots, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(folds, 2);
        if (n0 < folds || n1 < folds) throw new ArgumentException($"Each sample needs at least {folds} rows, one per fold.", n0 < folds ? nameof(source) : nameof(target));
        double[] grid = ridges is null ? [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0] : [.. ridges];
        if (grid.Length == 0) throw new ArgumentException("At least one ridge is required.", nameof(ridges));
        foreach (double r in grid)
            if (!(r > 0.0 && double.IsFinite(r))) throw new ArgumentOutOfRangeException(nameof(ridges), "Each ridge is positive and finite.");
        var expansion = new SplineExpansion(source, target, knots);
        double[,] s = expansion.Expand(source), t = expansion.Expand(target);
        double chosen = grid[0];
        if (grid.Length > 1)
        {
            double best = double.PositiveInfinity;
            foreach (double ridge in grid)
            {
                double loss = 0.0;
                for (int fold = 0; fold < folds && double.IsFinite(loss); fold++)
                {
                    DensityRatio fit;
                    try { fit = Fit(Rows(s, fold, folds, false), Rows(t, fold, folds, false), ridge, maxIterations, tolerance, idleConstants: true); }
                    catch (Exception e) when (e is InvalidOperationException or ArgumentException) { loss = double.PositiveInfinity; break; }
                    // The classifier's held-out log loss: −ln(1 − σ(η)) at a source row, −ln σ(η) at a target row.
                    foreach ((double[,] held, bool isTarget) in new[] { (Rows(s, fold, folds, true), false), (Rows(t, fold, folds, true), true) })
                    {
                        var row = new double[held.GetLength(1)];
                        for (int i = 0; i < held.GetLength(0); i++)
                        {
                            for (int j = 0; j < row.Length; j++) row[j] = held[i, j];
                            double eta = fit.Linear(row);
                            loss += isTarget ? Softplus(-eta) : Softplus(eta);
                        }
                    }
                }
                if (loss < best) { best = loss; chosen = ridge; }
            }
            if (double.IsPositiveInfinity(best)) throw new InvalidOperationException("No ridge's fit was accepted in every fold.");
        }
        DensityRatio whole = Fit(s, t, chosen, maxIterations, tolerance, idleConstants: false);
        return new DensityRatio(whole._beta, whole._mean, whole._scale, whole._offset, whole.Iterations, chosen, expansion);
    }

    // The rows of fold `fold` (held) or the rest, row i in fold i mod folds.
    private static double[,] Rows(double[,] rows, int fold, int folds, bool held)
    {
        int n = rows.GetLength(0), w = rows.GetLength(1), count = 0;
        for (int i = 0; i < n; i++) if ((i % folds == fold) == held) count++;
        var picked = new double[count, w];
        for (int i = 0, at = 0; i < n; i++)
        {
            if ((i % folds == fold) != held) continue;
            for (int j = 0; j < w; j++) picked[at, j] = rows[i, j];
            at++;
        }
        return picked;
    }

    /// <summary>Fits the ratio from <paramref name="source"/> (rows the source law's inputs' features) and <paramref name="target"/> (the
    /// target law's), each at least one row and the same width.</summary>
    /// <param name="source">The source law's sample, one row of features per input.</param>
    /// <param name="target">The target law's sample, its rows the same features.</param>
    /// <param name="ridge">The slopes' penalty per pooled point on the standardized features; positive.</param>
    /// <param name="maxIterations">Newton steps allowed before the fit is refused as not converged.</param>
    /// <param name="tolerance">The largest standardized gradient entry per pooled point at which the fit stops.</param>
    public static DensityRatio Logistic(double[,] source, double[,] target, double ridge = 1e-4, int maxIterations = 100, double tolerance = 1e-10)
        => Fit(source, target, ridge, maxIterations, tolerance, idleConstants: false);

    // The logistic fit; under idleConstants a constant feature is centred to zero and its slope stays at zero, where the public fit refuses it.
    private static DensityRatio Fit(double[,] source, double[,] target, double ridge, int maxIterations, double tolerance, bool idleConstants)
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
            if (idleConstants && scale[j] == 0.0) { scale[j] = 1.0; continue; }
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
            if (largest <= tolerance * n) return new DensityRatio(beta, mean, scale, Math.Log((double)n0 / n1), iteration - 1, ridge);
            factor.Decompose(hessian);
            if (factor.FailedAt >= 0) throw new InvalidOperationException($"The logistic fit's Hessian is not positive definite at step {iteration}.");
            double[] step = factor.Solve(gradient);
            // A Newton step below rounding on the standardized scale is the optimum: its gain is below the objective's own rounding, which
            // the halving would read as a loss.
            double stepSize = 0.0;
            foreach (double v in step) stepSize = Math.Max(stepSize, Math.Abs(v));
            if (stepSize <= 1e-12) return new DensityRatio(beta, mean, scale, Math.Log((double)n0 / n1), iteration - 1, ridge);
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

    /// <summary>The ridge the fit took: <see cref="Logistic"/>'s as given, <see cref="Spline"/>'s as cross-validation chose it.</summary>
    public double Ridge { get; }

    /// <summary>The inputs <see cref="LogRatio"/> reads: the features as given, which a spline fit expands itself.</summary>
    public int Inputs => _expansion?.Inputs ?? _mean.Length;

    /// <summary>ln w(x) at <paramref name="features"/>: the fitted log odds plus ln(n_source/n_target), the features expanded first where the
    /// fit took a spline expansion (<see cref="Spline"/>).</summary>
    public double LogRatio(ReadOnlySpan<double> features)
    {
        if (features.Length != Inputs) throw new ArgumentException($"{Inputs} features are required.", nameof(features));
        if (_expansion is null) return Linear(features) + _offset;
        var expanded = new double[_expansion.Width];
        _expansion.Expand(features, expanded);
        return Linear(expanded) + _offset;
    }

    // The fitted log odds at the standardized columns `columns`.
    private double Linear(ReadOnlySpan<double> columns)
    {
        double eta = _beta[0];
        for (int j = 0; j < _mean.Length; j++) eta += _beta[j + 1] * (columns[j] - _mean[j]) / _scale[j];
        return eta;
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

/// <summary>
/// Each input's cubic B-spline columns for <see cref="DensityRatio.Spline"/>: over the pooled samples' range, interior knots at the pooled
/// quantiles (repeats dropped), the first column of each input left out since an input's B-splines sum to one, and any column constant over
/// the pooled samples left out too. An input beyond the range holds its edge value.
/// </summary>
internal sealed class SplineExpansion
{
    private static readonly SplineBasis Basis = new(3);
    private readonly double[] _low, _span;
    private readonly KnotConfig[] _knots;
    private readonly int[][] _kept;   // each input's columns kept, by index into its B-splines

    internal SplineExpansion(double[,] source, double[,] target, int knots)
    {
        int n0 = source.GetLength(0), n1 = target.GetLength(0), p = source.GetLength(1);
        _low = new double[p];
        _span = new double[p];
        _knots = new KnotConfig[p];
        _kept = new int[p][];
        var pooled = new double[n0 + n1];
        for (int j = 0; j < p; j++)
        {
            for (int i = 0; i < n0; i++) pooled[i] = source[i, j];
            for (int i = 0; i < n1; i++) pooled[n0 + i] = target[i, j];
            Array.Sort(pooled);
            double low = pooled[0], high = pooled[^1];
            if (!(high > low) || !double.IsFinite(low) || !double.IsFinite(high))
                throw new ArgumentException($"Input {j} is constant or not finite across the samples, so it tells nothing apart.", nameof(source));
            var interior = new List<double>();
            for (int q = 1; q <= knots; q++)
            {
                double u = (pooled[(int)Math.Floor(q * (pooled.Length - 1) / (double)(knots + 1))] - low) / (high - low);
                if (u > 0.0 && u < 1.0 && (interior.Count == 0 || u > interior[^1])) interior.Add(u);
            }
            _low[j] = low;
            _span[j] = high - low;
            _knots[j] = new KnotConfig([.. interior]);
            // The columns that vary over the pooled samples, the first passed over. Some always vary: the identity is a combination of an
            // input's B-splines, so they cannot all be constant over two distinct values.
            double[,] design = Basis.Design(_knots[j], [.. pooled.Select(v => (v - low) / (high - low))]);
            var kept = new List<int>();
            for (int c = 1; c < design.GetLength(1); c++)
            {
                for (int i = 1; i < pooled.Length; i++)
                    if (design[i, c] != design[0, c]) { kept.Add(c); break; }
            }
            _kept[j] = [.. kept];
            Width += kept.Count;
        }
    }

    internal int Inputs => _low.Length;

    internal int Width { get; }

    internal void Expand(ReadOnlySpan<double> inputs, Span<double> into)
    {
        int at = 0;
        for (int j = 0; j < _low.Length; j++)
        {
            double u = Math.Clamp((inputs[j] - _low[j]) / _span[j], 0.0, 1.0);
            if (double.IsNaN(u))
            {
                foreach (int _ in _kept[j]) into[at++] = double.NaN;
                continue;
            }
            double[,] design = Basis.Design(_knots[j], [u]);
            foreach (int c in _kept[j]) into[at++] = design[0, c];
        }
    }

    internal double[,] Expand(double[,] rows)
    {
        int n = rows.GetLength(0);
        var expanded = new double[n, Width];
        var row = new double[_low.Length];
        var into = new double[Width];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < row.Length; j++) row[j] = rows[i, j];
            Expand(row, into);
            for (int c = 0; c < into.Length; c++) expanded[i, c] = into[c];
        }
        return expanded;
    }
}
