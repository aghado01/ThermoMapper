using System;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Rng;

namespace Maths.Regression.GPrior;

/// <summary>The matrix-t mixture Y | Σ ~ MN(M, V, Σ), Σ ~ IW(ν, Ψ), in the inverse Wishart's degrees of freedom.</summary>
public sealed class MatrixTPredictive
{
    private readonly double[,] _mean, _rows, _scale, _rowRoot;
    private readonly DenseCholesky _rowFactor;
    private readonly InverseWishart _sigma;
    private readonly double _constant;

    public MatrixTPredictive(double[,] location, double[,] rowScale, double[,] columnScale, double nu)
    {
        ArgumentNullException.ThrowIfNull(location);
        ArgumentNullException.ThrowIfNull(rowScale);
        ArgumentNullException.ThrowIfNull(columnScale);
        Cases = location.GetLength(0); Outputs = location.GetLength(1);
        if (Cases == 0 || Outputs == 0 || rowScale.GetLength(0) != Cases || rowScale.GetLength(1) != Cases || columnScale.GetLength(0) != Outputs)
            throw new ArgumentException("The location, row scale and output scale must have matching nonempty axes.");
        foreach (double x in location) if (!double.IsFinite(x)) throw new ArgumentException("The location must be finite.", nameof(location));
        for (int i = 0; i < Cases; i++) for (int j = 0; j < Cases; j++)
            if (!double.IsFinite(rowScale[i, j]) || rowScale[i, j] != rowScale[j, i]) throw new ArgumentException("The row scale must be finite and symmetric.", nameof(rowScale));
        _mean = (double[,])location.Clone(); _rows = (double[,])rowScale.Clone(); _scale = (double[,])columnScale.Clone();
        _sigma = InverseWishart.Of(nu, columnScale);
        _rowFactor = new DenseCholesky(Cases); _rowFactor.Decompose(rowScale);
        if (_rowFactor.FailedAt >= 0) throw new ArgumentException("The row scale must be positive definite.", nameof(rowScale));
        _rowRoot = Root(rowScale);
        Nu = nu;
        _constant = SpecialFunctions.LogMultivariateGamma(Outputs, (nu + Cases) / 2) - SpecialFunctions.LogMultivariateGamma(Outputs, nu / 2)
                    - 0.5 * Cases * Outputs * Math.Log(Math.PI) - 0.5 * Outputs * _rowFactor.LogDet + 0.5 * nu * _sigma.LogDeterminantScale;
    }

    public int Cases { get; }
    public int Outputs { get; }
    /// <summary>The mixing inverse Wishart's ν; one output-vector marginal has ν − K + 1 degrees of freedom.</summary>
    public double Nu { get; }
    public double[,] Location => (double[,])_mean.Clone();
    public double[,] RowScale => (double[,])_rows.Clone();
    public double[,] ColumnScale => (double[,])_scale.Clone();
    public double[,]? ExpectedSigma => _sigma.Mean;
    public bool HasMean => Nu > Outputs;
    public bool HasCovariance => Nu > Outputs + 1;

    public StudentT Marginal(int @case, int output)
        => StudentT.LocationScale(Nu - Outputs + 1, _mean[@case, output], Math.Sqrt(_rows[@case, @case] * _scale[output, output] / (Nu - Outputs + 1)));

    /// <summary>Cov(Y_ia,Y_jb), or null where its second moment does not exist.</summary>
    public double? Covariance(int i, int a, int j, int b)
        => HasCovariance ? _rows[i, j] * _scale[a, b] / (Nu - Outputs - 1) : null;

    public double LogDensity(double[,] values)
    {
        ArgumentNullException.ThrowIfNull(values);
        if (values.GetLength(0) != Cases || values.GetLength(1) != Outputs) throw new ArgumentException("The observation has the location's shape.", nameof(values));
        var solved = new double[Outputs][];
        for (int a = 0; a < Outputs; a++)
        {
            var e = new double[Cases];
            for (int i = 0; i < Cases; i++) { if (!double.IsFinite(values[i, a])) return double.NegativeInfinity; e[i] = values[i, a] - _mean[i, a]; }
            solved[a] = e; _rowFactor.ForwardSolveInPlace(solved[a]);
        }
        var total = (double[,])_scale.Clone();
        for (int a = 0; a < Outputs; a++) for (int b = 0; b <= a; b++)
        {
            double value = total[a, b];
            for (int i = 0; i < Cases; i++) value += solved[a][i] * solved[b][i];
            total[a, b] = total[b, a] = value;
        }
        var factor = new DenseCholesky(Outputs); factor.Decompose(total);
        if (factor.FailedAt >= 0) throw new InvalidOperationException("The predictive scale lost positive definiteness.");
        return _constant - 0.5 * (Nu + Cases) * factor.LogDet;
    }

    /// <summary>A joint draw, sharing Σ across all requested cases and the matrix-normal row covariance across their coefficients.</summary>
    public double[,] Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        double[,] column = Root(_sigma.Sample(rng));
        var z = new double[Cases, Outputs]; var row = new double[Cases, Outputs]; var sample = (double[,])_mean.Clone();
        for (int i = 0; i < Cases; i++) for (int a = 0; a < Outputs; a++) z[i, a] = Normal.SampleStandard(rng);
        for (int i = 0; i < Cases; i++) for (int a = 0; a < Outputs; a++) for (int j = 0; j <= i; j++) row[i, a] += _rowRoot[i, j] * z[j, a];
        for (int i = 0; i < Cases; i++) for (int a = 0; a < Outputs; a++) for (int b = 0; b <= a; b++) sample[i, a] += row[i, b] * column[a, b];
        return sample;
    }

    internal static double[,] Root(double[,] matrix)
    {
        int n = matrix.GetLength(0); var root = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++)
        {
            double value = matrix[i, j]; for (int k = 0; k < j; k++) value -= root[i, k] * root[j, k];
            if (i == j) { if (!(value > 0)) throw new ArgumentException("The matrix must be positive definite."); root[i, j] = Math.Sqrt(value); }
            else root[i, j] = value / root[j, j];
        }
        return root;
    }
}
