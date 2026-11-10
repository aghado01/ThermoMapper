using System;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Tests;

// Dense conditioning, independent of the engine's pencils and prediction readouts.
internal sealed record MultiOutputPredictiveReference(double[,] Mean, double[,] Precision, double[,] Scale,
                                                     double Nu, bool Flat, double[] ColumnMeans, double[] OutputMeans, int Rows)
{
    internal static MultiOutputPredictiveReference Fit(MultiOutputCase c, ReferenceSigmaPrior prior)
    {
        if (c.InverseTemperature != 1) throw new ArgumentException("The reference predicts at the posterior.");
        double[,] x = c.FlatIntercept ? ReferenceAlgebra.Centred(c.X) : c.X;
        double[,] y = c.FlatIntercept ? ReferenceAlgebra.Centred(c.Y) : c.Y;
        double[,] gram = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(x), x);
        double[,] p = ReferenceAlgebra.Sum(gram, c.PriorPrecision);
        double[,] rhs = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(x), y);
        rhs = ReferenceAlgebra.Sum(rhs, ReferenceAlgebra.Multiply(c.PriorPrecision, c.PriorMean ?? new double[c.Columns, c.Outputs]));
        var mean = ReferenceAlgebra.Solve(ReferenceAlgebra.Cholesky(p)!, rhs);
        var (nu, scale) = MatrixTReference.SigmaPosterior(MatrixTReference.Collapse(c), prior);
        var xm = new double[c.Columns];
        var ym = new double[c.Outputs];
        for (int i = 0; i < c.Rows; i++)
        {
            for (int j = 0; j < xm.Length; j++) xm[j] += c.X[i, j] / c.Rows;
            for (int k = 0; k < ym.Length; k++) ym[k] += c.Y[i, k] / c.Rows;
        }
        return new(mean, p, scale, nu, c.FlatIntercept, xm, ym, c.Rows);
    }

    internal (double[,] Mean, double[,] RowScale) At(double[,] cases, bool predictive = true)
    {
        int n = cases.GetLength(0), p = cases.GetLength(1), k = Mean.GetLength(1);
        var z = (double[,])cases.Clone();
        if (Flat) for (int i = 0; i < n; i++) for (int j = 0; j < p; j++) z[i, j] -= ColumnMeans[j];
        double[,] location = ReferenceAlgebra.Multiply(z, Mean);
        if (Flat) for (int i = 0; i < n; i++) for (int l = 0; l < k; l++) location[i, l] += OutputMeans[l];
        var v = ReferenceAlgebra.Multiply(ReferenceAlgebra.Multiply(z, ReferenceAlgebra.Inverse(Precision)), ReferenceAlgebra.Transpose(z));
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) v[i, j] += (Flat ? 1.0 / Rows : 0) + (predictive && i == j ? 1 : 0);
        return (location, v);
    }

    internal double LogDensity(double[,] cases, double[,] values)
    {
        var (mean, v) = At(cases);
        return MatrixTReference.MatrixTLogDensity(values, mean, v, Scale, Nu);
    }

    internal double[,] Sample(double[,] cases, Xoshiro256PlusPlus rng)
    {
        var (mean, v) = At(cases);
        var row = ReferenceAlgebra.Cholesky(v)!;
        var column = ReferenceAlgebra.Cholesky(ReferenceInverseWishart.Sample(Nu, Scale, rng))!;
        int n = mean.GetLength(0), k = mean.GetLength(1);
        var z = new double[n, k];
        for (int i = 0; i < n; i++) for (int l = 0; l < k; l++) z[i, l] = Normal.SampleStandard(rng);
        var noise = ReferenceAlgebra.Multiply(ReferenceAlgebra.Multiply(row, z), ReferenceAlgebra.Transpose(column));
        return ReferenceAlgebra.Sum(mean, noise);
    }
}

// Fixed mean-centred intercept structure. A conjugate Gibbs reference for the HW covariance mixture:
// r_k | Sigma ~ Gamma((nu+K)/2, rate (Sigma^-1)_kk + 1/(nu A_k^2)).
internal sealed class CovarianceGibbsReference(double[,] residual, int rows, double[] scales, bool huangWand)
{
    private readonly double[] _rates = Array.ConvertAll(scales, a => a * a);
    private readonly double[,]? _fixedRoot = huangWand ? null : ReferenceInverseWishart.SamplingRoot(Scaled(residual, rows));
    internal double Nu => rows + scales.Length - 1 + (huangWand ? 2 : 0.002);
    internal ReadOnlySpan<double> Rates => _rates;

    private static double[,] Scaled(double[,] r, int n)
    {
        var scale = (double[,])r.Clone();
        for (int i = 0; i < scale.GetLength(0); i++) for (int j = 0; j < scale.GetLength(0); j++) scale[i, j] *= 1 + 0.002 / n;
        return scale;
    }

    internal double[,] Step(Xoshiro256PlusPlus rng)
    {
        if (_fixedRoot is { } root) return ReferenceInverseWishart.SampleFromRoot(Nu, root, rng);
        var psi = (double[,])residual.Clone();
        for (int k = 0; k < scales.Length; k++) psi[k, k] += 2 * _rates[k];
        double[,] sigma = ReferenceInverseWishart.Sample(Nu, psi, rng);
        double[,] inverse = ReferenceAlgebra.Inverse(sigma);
        for (int k = 0; k < scales.Length; k++)
            _rates[k] = Gamma.ShapeRate((2.0 + scales.Length) / 2, inverse[k, k] + 1 / (2 * scales[k] * scales[k])).Sample(rng);
        return sigma;
    }
}
