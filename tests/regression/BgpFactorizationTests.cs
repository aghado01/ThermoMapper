using System;
using Maths.LinAlg;
using Maths.Regression.Bgp;
using Xunit;

namespace Maths.Regression.Tests;

public sealed class BgpFactorizationTests
{
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Admission_RetainsTheHistoricalPivotBoundary(bool unknownScale)
    {
        double[,] x = { { 0.0 }, { 1.0 } };
        double[] y = { 0.0, 1.0 };
        // A zero kernel makes the raw pivot exactly the nugget: test both adjacent representable values.
        foreach (double pivot in new[] { Math.BitDecrement(1e-12), 1e-12, Math.BitIncrement(1e-12) })
        {
            var gp = unknownScale
                ? GpRegression.WithUnknownScale(x, y, new ConstantKernel(0.0))
                : new GpRegression(x, y, pivot, new ConstantKernel(0.0));
            GpFit Fit() => unknownScale ? gp.Fit(0.1, pivot) : gp.Fit(0.1);
            if (pivot <= 1e-12) Assert.Throws<InvalidOperationException>(() => Fit());
            else Assert.True(double.IsFinite(Fit().LogMarginal));
        }
        // A nearly rank-one covariance exercises a small positive Schur pivot after cancellation.
        foreach (double nugget in new[] { 4e-13, 6e-13 })
        {
            var gp = unknownScale
                ? GpRegression.WithUnknownScale(x, y, new ConstantKernel(1.0))
                : new GpRegression(x, y, nugget, new ConstantKernel(1.0));
            GpFit Fit() => unknownScale ? gp.Fit(0.1, nugget) : gp.Fit(0.1);
            if (nugget == 4e-13) Assert.Throws<InvalidOperationException>(() => Fit());
            else Assert.True(double.IsFinite(Fit().LogMarginal));
        }
    }

    [Fact]
    public void PrecisionBlocks_MatchTheCovarianceInverse_WithUnorderedMembers()
    {
        double[,] covariance = { { 2.0, 0.4, -0.2, 0.1 }, { 0.4, 1.0, 0.3, -0.1 },
                                 { -0.2, 0.3, 1.5, 0.2 }, { 0.1, -0.1, 0.2, 0.8 } };
        var factor = new DenseCholesky(4);
        factor.Decompose(covariance);
        var precision = new GpPrecision(factor);
        int[] members = { 3, 0, 2 };
        double[,] block = precision.Block(members);
        for (int b = 0; b < members.Length; b++)
        {
            var rhs = new double[4];
            rhs[members[b]] = 1.0;
            double[] column = factor.Solve(rhs);
            for (int a = 0; a < members.Length; a++)
                Assert.InRange(Math.Abs(block[a, b] - column[members[a]]), 0.0, 1e-14);
        }
    }

    private sealed class ConstantKernel(double value) : IGpKernel
    {
        public double Evaluate(double squaredDistance, double bandwidth) => value;
    }
}
