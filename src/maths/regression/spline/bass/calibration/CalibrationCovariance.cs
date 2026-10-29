using System;
using System.Text.Json.Serialization;
using Maths.LinAlg;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// The likelihood's covariance formulation (D09), declared rather than read from the discrepancy's shape (R branches on
/// <c>ncol(discrep.mat) &lt; length(y)</c>, C-O1). Shapes and positive definiteness are checked when it is built
/// (<see cref="CalibrationCovariance.Create"/>), and a mismatch is refused, never resolved by switching formulations.
/// </summary>
[JsonPolymorphic(TypeDiscriminatorPropertyName = "Kind")]
[JsonDerivedType(typeof(FunctionalFormulation), "Functional")]
[JsonDerivedType(typeof(MultivariateFormulation), "Multivariate")]
public abstract record CovarianceFormulation
{
    private protected CovarianceFormulation() { }

    /// <summary>The truncation term (<see cref="TruncationTerm"/>).</summary>
    public TruncationTerm Truncation { get; init; } = TruncationTerm.ResidualSet;
}

/// <summary>How the emulator's truncation residuals enter the covariance (D07, D09).</summary>
public enum TruncationTerm
{
    /// <summary>
    /// D07's residual set whole: (1/n)·Σᵢ (y_s ∘ eᵢ)(y_s ∘ eᵢ)ᵀ, the covariance of D07's uniform draw over the n residual
    /// vectors taken about zero, as that draw is, so the residuals keep their correlation across outputs. A low-rank term
    /// of rank at most n in the functional formulation. The default.
    /// </summary>
    ResidualSet,

    /// <summary>
    /// R's type 1: each output's sample variance of the residuals (centred, n − 1), y_s² ∘ var(e), on the diagonal. An
    /// approximation that over-counts a dense output grid when the misfit is smooth.
    /// </summary>
    Diagonal,
}

/// <summary>
/// R's type 1, for dense outputs: Σ = diag(σ²_{g(j)}) + the emulator's noise + the truncation term, a diagonal plus a low-
/// rank term solved through <see cref="DiagonalPlusLowRank"/> without an m × m matrix, with the discrepancy δ = μ_δ + D·v
/// on a supplied basis and v drawn from its Gaussian conditional (<see cref="FunctionalDiscrepancy"/>).
/// </summary>
public sealed record FunctionalFormulation : CovarianceFormulation
{
    /// <summary>The discrepancy; null for none.</summary>
    public FunctionalDiscrepancy? Discrepancy { get; init; }
}

/// <summary>
/// The functional discrepancy δ = μ_δ + D·v, v ~ N(0, s_v²·I) (D09): <paramref name="Basis"/> D is m × q, row j for output j;
/// <paramref name="Scale"/> s_v has no default and is refused when missing (R fixes it at 1, which ties the posterior to
/// D's units whenever D is a dimensionless basis; declaring 1 reproduces R); <paramref name="Mean"/> μ_δ, m values, null
/// for zero.
/// </summary>
public sealed record FunctionalDiscrepancy(double[][] Basis, double? Scale, double[]? Mean = null);

/// <summary>
/// R's type 2, for short outputs: Σ = Σ_trunc + diag(σ_{g(j)})·R·diag(σ_{g(j)}) + Σ_δ + the emulator's noise, on the dense
/// Cholesky factor, with <paramref name="Correlation"/> R over the outputs (m × m, unit diagonal, positive definite), a
/// fixed discrepancy covariance Σ_δ (null for none) and a fixed mean shift μ_δ (null for zero).
/// </summary>
public sealed record MultivariateFormulation(double[][] Correlation, double[][]? DiscrepancyCovariance = null, double[]? DiscrepancyMean = null)
    : CovarianceFormulation;

/// <summary>
/// A formulation bound to its measurement model and emulator, validated: what stays fixed (the truncation term, the
/// discrepancy's declaration) is built once, and <see cref="Factor"/> assembles Σ at the group variances and an emulator
/// draw's noise, in physical units throughout (D07). Immutable; concurrent chains may share it.
/// </summary>
public abstract class CalibrationCovariance
{
    private protected CalibrationCovariance(MeasurementModel measurement, double[]? mean)
    {
        Measurement = measurement;
        Mean = mean;
    }

    /// <summary>The measurement model the variances belong to.</summary>
    public MeasurementModel Measurement { get; }

    /// <summary>The outputs m.</summary>
    public int Outputs => Measurement.Outputs;

    /// <summary>The fixed mean shift μ_δ, null for zero.</summary>
    public double[]? Mean { get; }

    /// <summary>The functional discrepancy's basis D (m × q), null when there is none to sample.</summary>
    public virtual double[,]? DiscrepancyBasis => null;

    /// <summary>The functional discrepancy's scale s_v; 0 when there is none to sample.</summary>
    public virtual double DiscrepancyScale => 0.0;

    /// <summary>Σ at <paramref name="groupVariances"/> and the draw's <paramref name="noise"/>, factored.</summary>
    public abstract CovarianceFactor Factor(ReadOnlySpan<double> groupVariances, EmulatorNoise noise);

    /// <summary>
    /// The formulation <paramref name="declaration"/> declares, over <paramref name="measurement"/>'s outputs and
    /// <paramref name="emulator"/>'s draws. Refused: a measurement model over another output count than the emulator's;
    /// a discrepancy, correlation or mean of the wrong shape, not finite, not symmetric or not positive definite; a
    /// missing or nonpositive s_v; the diagonal truncation term on an emulator without truncation residuals.
    /// </summary>
    public static CalibrationCovariance Create(CovarianceFormulation declaration, MeasurementModel measurement, ICalibrationEmulator emulator)
    {
        ArgumentNullException.ThrowIfNull(declaration);
        ArgumentNullException.ThrowIfNull(measurement);
        ArgumentNullException.ThrowIfNull(emulator);
        int m = measurement.Outputs;
        if (emulator.Outputs != m)
            throw new ArgumentException($"The emulator predicts {emulator.Outputs} outputs and the measurement model covers {m}.", nameof(emulator));
        if (declaration.Truncation is not (TruncationTerm.ResidualSet or TruncationTerm.Diagonal))
            throw new ArgumentException("The truncation term is of no known kind.", nameof(declaration));
        double[,]? residuals = emulator.TruncationResiduals();
        if (residuals is not null)
        {
            if (residuals.GetLength(1) != m) throw new ArgumentException($"The truncation residuals have {residuals.GetLength(1)} columns, not {m}.", nameof(emulator));
            foreach (double v in residuals)
                if (!double.IsFinite(v)) throw new ArgumentException("The truncation residuals are not finite.", nameof(emulator));
        }
        if (declaration.Truncation == TruncationTerm.Diagonal)
        {
            if (residuals is null) throw new ArgumentException("The diagonal truncation term needs the emulator's truncation residuals, and it has none.", nameof(declaration));
            if (residuals.GetLength(0) < 2) throw new ArgumentException("The diagonal truncation term's sample variance needs at least two residual vectors.", nameof(declaration));
        }
        return declaration switch
        {
            FunctionalFormulation f => new FunctionalCovariance(f, measurement, residuals),
            MultivariateFormulation v => new MultivariateCovariance(v, measurement, residuals),
            _ => throw new ArgumentException("The covariance formulation is of no known kind.", nameof(declaration)),
        };
    }

    // y_s² ∘ var(e) per output: R's diag(cov(t(trunc.error))), centred with n − 1, from the physical residuals.
    private protected static double[] SampleVariance(double[,] residuals)
    {
        int n = residuals.GetLength(0), m = residuals.GetLength(1);
        var variance = new double[m];
        var column = new double[n];
        for (int j = 0; j < m; j++)
        {
            for (int i = 0; i < n; i++) column[i] = residuals[i, j];
            double mean = PairwiseSum.Of(column) / n;
            for (int i = 0; i < n; i++) column[i] = (column[i] - mean) * (column[i] - mean);
            variance[j] = PairwiseSum.Of(column) / (n - 1);
        }
        return variance;
    }

    private protected static double[] Vector(double[]? values, int m, string name, string parameter)
    {
        if (values is null) throw new ArgumentException($"{name} is missing.", parameter);
        if (values.Length != m) throw new ArgumentException($"{name} has {values.Length} values; the outputs are {m}.", parameter);
        foreach (double v in values)
            if (!double.IsFinite(v)) throw new ArgumentException($"{name} is not finite.", parameter);
        return (double[])values.Clone();
    }

    private protected static double[,] Matrix(double[][]? rows, int m, int? columns, string name, string parameter)
    {
        if (rows is null) throw new ArgumentException($"{name} is missing.", parameter);
        if (rows.Length != m) throw new ArgumentException($"{name} has {rows.Length} rows; the outputs are {m}.", parameter);
        int q = columns ?? rows[0]?.Length ?? 0;
        if (q < 1) throw new ArgumentException($"{name} has no column.", parameter);
        var a = new double[m, q];
        for (int i = 0; i < m; i++)
        {
            if (rows[i] is null || rows[i].Length != q) throw new ArgumentException($"{name}'s rows must each hold {q} values.", parameter);
            for (int j = 0; j < q; j++)
            {
                if (!double.IsFinite(rows[i][j])) throw new ArgumentException($"{name} is not finite.", parameter);
                a[i, j] = rows[i][j];
            }
        }
        return a;
    }

    // A symmetric positive definite m × m matrix, checked by its factor.
    private protected static double[,] SymmetricPositiveDefinite(double[][]? rows, int m, string name, string parameter)
    {
        double[,] a = Matrix(rows, m, m, name, parameter);
        double scale = 0.0;
        for (int i = 0; i < m; i++) scale = Math.Max(scale, Math.Abs(a[i, i]));
        for (int i = 0; i < m; i++)
            for (int j = 0; j < i; j++)
                if (Math.Abs(a[i, j] - a[j, i]) > 1e-12 * scale)
                    throw new ArgumentException($"{name} is not symmetric: entries ({i}, {j}) and ({j}, {i}) differ.", parameter);
        var factor = new DenseCholesky(m);
        factor.Decompose(a);
        if (!factor.Status.PositiveDefinite) throw new ArgumentException($"{name} is not positive definite.", parameter);
        return a;
    }
}

/// <summary>The functional formulation built (<see cref="FunctionalFormulation"/>).</summary>
internal sealed class FunctionalCovariance : CalibrationCovariance
{
    private readonly double[]? _truncationDiagonal;     // the diagonal option's y_s² ∘ var(e)
    private readonly double[,]? _truncationColumns;     // the residual set's second moment as columns: its thin SVD's Q·diag(λ), m × rank
    private readonly double[,]? _discrepancy;
    private readonly double _scale;
    private Assembled? _assembled;                       // the last component factor, with the truncation columns beside it

    private sealed record Assembled(double[,] Factor, double[,] Columns);   // swapped whole, so concurrent readers see a pair

    public FunctionalCovariance(FunctionalFormulation declaration, MeasurementModel measurement, double[,]? residuals)
        : base(measurement, declaration.Discrepancy?.Mean is { } mean ? Vector(mean, measurement.Outputs, "The discrepancy's mean", nameof(declaration)) : null)
    {
        int m = measurement.Outputs;
        if (residuals is not null)
        {
            if (declaration.Truncation == TruncationTerm.Diagonal) _truncationDiagonal = SampleVariance(residuals);
            else _truncationColumns = Compressed(residuals);
        }
        if (declaration.Discrepancy is { } d)
        {
            _discrepancy = Matrix(d.Basis, m, null, "The discrepancy's basis", nameof(declaration));
            if (d.Scale is not { } s) throw new ArgumentException("The discrepancy's scale s_v has no default: declare it (1 reproduces R).", nameof(declaration));
            if (!(s > 0.0 && double.IsFinite(s))) throw new ArgumentException($"The discrepancy's scale s_v must be positive and finite, not {s}.", nameof(declaration));
            _scale = s;
        }
    }

    public override double[,]? DiscrepancyBasis => _discrepancy;
    public override double DiscrepancyScale => _scale;

    /// <summary>The residual set's term (1/n)·Σᵢ rᵢrᵢᵀ as columns C with CCᵀ equal to it: the thin SVD of [r₁ … r_n]/√n,
    /// Q·diag(λ) over its numerical rank, at most min(m, n) columns, taken once. A direction below the rank tolerance,
    /// max(m, n)·ε·λ₁, contributes under the term's own rounding. Null when the residuals are zero.</summary>
    internal static double[,]? Compressed(double[,] residuals)
    {
        int n = residuals.GetLength(0), m = residuals.GetLength(1);
        double root = Math.Sqrt(n);
        var a = new double[m, n];
        bool any = false;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
            {
                a[j, i] = residuals[i, j] / root;
                any |= a[j, i] != 0.0;
            }
        if (!any) return null;
        ThinSvdResult svd = ThinSvd.Decompose(a);
        int rank = svd.NumericalRank;
        if (rank == 0) return null;
        var columns = new double[m, rank];
        for (int k = 0; k < rank; k++)
            for (int j = 0; j < m; j++) columns[j, k] = svd.LeftVectors[k][j] * svd.SingularValues[k];
        return columns;
    }

    public override CovarianceFactor Factor(ReadOnlySpan<double> groupVariances, EmulatorNoise noise)
    {
        ArgumentNullException.ThrowIfNull(noise);
        int m = Outputs;
        var diagonal = new double[m];
        Measurement.OutputVariances(groupVariances, diagonal);
        if (_truncationDiagonal is not null)
            for (int j = 0; j < m; j++) diagonal[j] += _truncationDiagonal[j];

        double[,]? columns = null;
        double[] weights = [];
        switch (noise)
        {
            case IndependentNoise independent:
                CheckVariance(independent.Variance);
                for (int j = 0; j < m; j++) diagonal[j] += independent.Variance;
                columns = _truncationColumns;
                if (columns is not null)
                {
                    weights = new double[columns.GetLength(1)];
                    Array.Fill(weights, 1.0);
                }
                break;
            case ComponentNoise component:
                (columns, weights) = WithComponents(component);
                break;
            default:
                throw new ArgumentException("The emulator's noise is of no known kind.", nameof(noise));
        }
        return columns is null ? new DiagonalFactor(diagonal) : new LowRankFactor(new DiagonalPlusLowRank(diagonal, columns, weights));
    }

    // [y_s ∘ K | (y_s ∘ E)/√n] with weights [σ²_k | 1], the concatenation kept while the draws share their factor.
    private (double[,] Columns, double[] Weights) WithComponents(ComponentNoise noise)
    {
        double[,] factor = noise.Factor ?? throw new ArgumentException("The component noise has no factor.", nameof(noise));
        int m = Outputs, r = factor.GetLength(1), n = _truncationColumns?.GetLength(1) ?? 0;
        if (factor.GetLength(0) != m) throw new ArgumentException($"The component noise's factor has {factor.GetLength(0)} rows, not {m}.", nameof(noise));
        if (noise.Variances is null || noise.Variances.Length != r) throw new ArgumentException($"The component noise needs {r} variances.", nameof(noise));
        foreach (double v in noise.Variances) CheckVariance(v);
        double[,] columns;
        if (n == 0) columns = factor;
        else if (_assembled is { } cached && ReferenceEquals(cached.Factor, factor)) columns = cached.Columns;
        else
        {
            columns = new double[m, r + n];
            for (int j = 0; j < m; j++)
            {
                for (int k = 0; k < r; k++) columns[j, k] = factor[j, k];
                for (int i = 0; i < n; i++) columns[j, r + i] = _truncationColumns![j, i];
            }
            _assembled = new Assembled(factor, columns);
        }
        var weights = new double[r + n];
        noise.Variances.CopyTo(weights, 0);
        for (int i = r; i < weights.Length; i++) weights[i] = 1.0;
        return (columns, weights);
    }

    internal static void CheckVariance(double v)
    {
        if (!(v >= 0.0 && double.IsFinite(v))) throw new ArgumentException($"An emulator noise variance is {v}; each must be finite and not negative.", "noise");
    }
}

/// <summary>The multivariate formulation built (<see cref="MultivariateFormulation"/>).</summary>
internal sealed class MultivariateCovariance : CalibrationCovariance
{
    private readonly double[,] _fixed;          // Σ_trunc + Σ_δ
    private readonly double[,] _correlation;

    public MultivariateCovariance(MultivariateFormulation declaration, MeasurementModel measurement, double[,]? residuals)
        : base(measurement, declaration.DiscrepancyMean is { } mean ? Vector(mean, measurement.Outputs, "The discrepancy's mean", nameof(declaration)) : null)
    {
        int m = measurement.Outputs;
        _correlation = SymmetricPositiveDefinite(declaration.Correlation, m, "The correlation", nameof(declaration));
        for (int j = 0; j < m; j++)
            if (Math.Abs(_correlation[j, j] - 1.0) > 1e-12) throw new ArgumentException($"The correlation's diagonal entry {j} is {_correlation[j, j]}, not 1.", nameof(declaration));
        _fixed = declaration.DiscrepancyCovariance is null
            ? new double[m, m]
            : SymmetricPositiveDefinite(declaration.DiscrepancyCovariance, m, "The discrepancy's covariance", nameof(declaration));
        if (residuals is null) return;
        if (declaration.Truncation == TruncationTerm.Diagonal)
        {
            double[] variance = SampleVariance(residuals);
            for (int j = 0; j < m; j++) _fixed[j, j] += variance[j];
            return;
        }
        int n = residuals.GetLength(0);
        for (int j = 0; j < m; j++)
            for (int l = 0; l <= j; l++)
            {
                double s = 0.0;
                for (int i = 0; i < n; i++) s += residuals[i, j] * residuals[i, l];
                _fixed[j, l] += s / n;
                if (l != j) _fixed[l, j] += s / n;
            }
    }

    public override CovarianceFactor Factor(ReadOnlySpan<double> groupVariances, EmulatorNoise noise)
    {
        ArgumentNullException.ThrowIfNull(noise);
        int m = Outputs;
        var sd = new double[m];
        Measurement.OutputVariances(groupVariances, sd);
        for (int j = 0; j < m; j++) sd[j] = Math.Sqrt(sd[j]);
        var sigma = new double[m, m];
        for (int j = 0; j < m; j++)
            for (int l = 0; l <= j; l++) sigma[j, l] = _fixed[j, l] + sd[j] * _correlation[j, l] * sd[l];
        switch (noise)
        {
            case IndependentNoise independent:
                FunctionalCovariance.CheckVariance(independent.Variance);
                for (int j = 0; j < m; j++) sigma[j, j] += independent.Variance;
                break;
            case ComponentNoise component:
                double[,] u = component.Factor ?? throw new ArgumentException("The component noise has no factor.", nameof(noise));
                int r = u.GetLength(1);
                if (u.GetLength(0) != m || component.Variances is null || component.Variances.Length != r)
                    throw new ArgumentException($"The component noise needs an {m}-row factor and one variance per column.", nameof(noise));
                foreach (double v in component.Variances) FunctionalCovariance.CheckVariance(v);
                for (int j = 0; j < m; j++)
                    for (int l = 0; l <= j; l++)
                    {
                        double s = 0.0;
                        for (int k = 0; k < r; k++) s += u[j, k] * component.Variances[k] * u[l, k];
                        sigma[j, l] += s;
                    }
                break;
            default:
                throw new ArgumentException("The emulator's noise is of no known kind.", nameof(noise));
        }
        var factor = new DenseCholesky(m);
        factor.Decompose(sigma);   // the lower triangle
        if (!factor.Status.PositiveDefinite) throw new InvalidOperationException("The multivariate covariance did not factor as positive definite.");
        return new DenseFactor(factor);
    }
}

/// <summary>Σ factored: its log-determinant, solves and the quadratic form bᵀΣ⁻¹b, the likelihood's pieces.</summary>
public abstract class CovarianceFactor
{
    private protected CovarianceFactor() { }

    /// <summary>The order m.</summary>
    public abstract int Dimension { get; }

    /// <summary>ln|Σ|.</summary>
    public abstract double LogDet { get; }

    /// <summary>bᵀΣ⁻¹b.</summary>
    public abstract double InverseQuadraticForm(ReadOnlySpan<double> b);

    /// <summary>Σ⁻¹b, in place.</summary>
    public abstract void SolveInPlace(Span<double> x);

    /// <summary>Σ⁻¹b.</summary>
    public double[] Solve(ReadOnlySpan<double> b)
    {
        double[] x = b.ToArray();
        SolveInPlace(x);
        return x;
    }
}

internal sealed class DiagonalFactor : CovarianceFactor
{
    private readonly double[] _diagonal;

    public DiagonalFactor(double[] diagonal)
    {
        double logDet = 0.0;
        foreach (double d in diagonal)
        {
            if (!(d > 0.0 && double.IsFinite(d))) throw new ArgumentException($"A diagonal entry is {d}; each must be finite and positive.", nameof(diagonal));
            logDet += Math.Log(d);
        }
        _diagonal = diagonal;
        LogDet = logDet;
    }

    public override int Dimension => _diagonal.Length;
    public override double LogDet { get; }

    public override double InverseQuadraticForm(ReadOnlySpan<double> b)
    {
        double q = 0.0;
        for (int i = 0; i < _diagonal.Length; i++) q += b[i] * b[i] / _diagonal[i];
        return q;
    }

    public override void SolveInPlace(Span<double> x)
    {
        for (int i = 0; i < _diagonal.Length; i++) x[i] /= _diagonal[i];
    }
}

internal sealed class LowRankFactor(DiagonalPlusLowRank sigma) : CovarianceFactor
{
    public override int Dimension => sigma.Dimension;
    public override double LogDet => sigma.LogDet;
    public override double InverseQuadraticForm(ReadOnlySpan<double> b) => sigma.InverseQuadraticForm(b);
    public override void SolveInPlace(Span<double> x) => sigma.SolveInPlace(x);
}

internal sealed class DenseFactor(DenseCholesky factor) : CovarianceFactor
{
    public override int Dimension => factor.Dimension;
    public override double LogDet => factor.LogDet;
    public override double InverseQuadraticForm(ReadOnlySpan<double> b) => factor.InverseQuadraticForm(b);
    public override void SolveInPlace(Span<double> x) => factor.SolveInPlace(x);
}
