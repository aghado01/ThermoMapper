using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class MultiOutputPredictiveReferenceTests(ITestOutputHelper output)
{
    [Fact]
    public void ConjugateCovarianceGibbsMeetsIndependentLogRateQuadrature()
    {
        var response = (Maths.Regression.Spline.Bass.MultiOutputResponse)MultiOutputTemperingReference.Data().Response;
        int n = response.Cases;
        var y = new double[n, 2];
        for (int i = 0; i < n; i++) for (int k = 0; k < 2; k++) y[i, k] = response[i, k];
        var centred = ReferenceAlgebra.Centred(y);
        var r = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(centred), centred);
        double[] scales = [1, 0.7];
        var nodes = BassTargetTests.Nodes(-30, 8, 1).ToArray();
        var logs = new double[nodes.Length * nodes.Length];
        var rates = new (double A, double B)[logs.Length];
        int at = 0;
        foreach (var a in nodes) foreach (var b in nodes)
        {
            double ra = Math.Exp(a.U), rb = Math.Exp(b.U);
            var total = new double[,] { { r[0, 0] + 2 * ra, r[0, 1] }, { r[1, 0], r[1, 1] + 2 * rb } };
            double logDet = ReferenceAlgebra.LogDeterminant(ReferenceAlgebra.Cholesky(total)!);
            logs[at] = Math.Log(a.W * b.W) + 2 * (a.U + b.U) - ra / 2 - rb / 0.98 - 0.5 * (n + 3) * logDet;
            rates[at++] = (ra, rb);
        }
        double logZ = BassTemperingTests.LogSumExp(logs), ea = 0, eb = 0;
        for (int i = 0; i < logs.Length; i++) { double w = Math.Exp(logs[i] - logZ); ea += w * rates[i].A; eb += w * rates[i].B; }
        var oracle = new CovarianceGibbsReference(r, n, scales, true);
        var rng = new Xoshiro256PlusPlus(20261214);
        for (int i = 0; i < 2000; i++) oracle.Step(rng);
        var series = new[] { new double[60_000], new double[60_000] };
        for (int i = 0; i < series[0].Length; i++) { oracle.Step(rng); series[0][i] = oracle.Rates[0]; series[1][i] = oracle.Rates[1]; }
        double[] truth = [ea, eb];
        for (int j = 0; j < 2; j++)
        {
            double z = (series[j].Average() - truth[j]) / KhaosChainTests.BatchError(series[j], 50);
            output.WriteLine($"rate {j}, mean {series[j].Average():G8}, quadrature {truth[j]:G8}, z={z:F3}");
            Assert.InRange(Math.Abs(z), 0, 4.5);
        }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void DenseConditioningAgreesWithTheAugmentedMarginalRatio(bool flat)
    {
        var rng = new Xoshiro256PlusPlus(20261211);
        double worst = 0;
        for (int trial = 0; trial < 100; trial++)
        {
            int n = 7, p = 2, k = 1 + trial % 3, q = 1 + trial % 2;
            var x = new double[n, p]; var y = new double[n, k];
            var z = new double[q, p]; var values = new double[q, k];
            for (int i = 0; i < n; i++) { for (int j = 0; j < p; j++) x[i, j] = Normal.SampleStandard(rng); for (int l = 0; l < k; l++) y[i, l] = Normal.SampleStandard(rng); }
            for (int i = 0; i < q; i++) { for (int j = 0; j < p; j++) z[i, j] = Normal.SampleStandard(rng); for (int l = 0; l < k; l++) values[i, l] = Normal.SampleStandard(rng); }
            var psi = new double[k, k]; for (int l = 0; l < k; l++) psi[l, l] = 0.3 + l;
            var prior = new ReferenceSigmaPrior(k + 3, psi);
            var c = new MultiOutputCase(y, x, new double[,] { { 0.7, 0.2 }, { 0.2, 1.4 } }, null, 1, flat);
            var allX = new double[n + q, p]; var allY = new double[n + q, k];
            for (int i = 0; i < n + q; i++) { for (int j = 0; j < p; j++) allX[i, j] = i < n ? x[i, j] : z[i - n, j]; for (int l = 0; l < k; l++) allY[i, l] = i < n ? y[i, l] : values[i - n, l]; }
            double ratio = MatrixTReference.LogMarginal(MatrixTReference.Collapse(c with { X = allX, Y = allY }), prior) - MatrixTReference.LogMarginal(MatrixTReference.Collapse(c), prior);
            double predictive = MultiOutputPredictiveReference.Fit(c, prior).LogDensity(z, values);
            worst = Math.Max(worst, Math.Abs(ratio - predictive));
        }
        output.WriteLine($"flat={flat}: worst augmented-ratio error {worst:G5}");
        Assert.InRange(worst, 0, 2e-12);
    }

    [Fact]
    public void SingleCasePredictiveIsAMultivariateStudentT_AndDrawsHaveItsMoments()
    {
        var model = Maths.Regression.Spline.Bass.BassModel.Create(MultiOutputCoefficientTests.Config("ridge"), MultiOutputCoefficientTests.Data(2));
        var d = model.NewDesign(); Assert.Null(d.Load(MultiOutputCoefficientTests.Terms("ridge")));
        var c = MultiOutputCoefficientTests.DenseCase(model, d, 0.4, 1);
        var prior = new ReferenceSigmaPrior(model.Resolved.Outputs!.Prior.DegreesOfFreedom, model.Resolved.Outputs.Prior.Scale);
        var reference = MultiOutputPredictiveReference.Fit(c, prior);
        var z = new double[1, c.Columns]; for (int j = 0; j < c.Columns; j++) z[0, j] = c.X[0, j];
        var (mean, row) = reference.At(z);
        double df = reference.Nu - 1;
        var value = new double[,] { { mean[0, 0] + 0.7, mean[0, 1] - 0.2 } };
        var e = new[] { 0.7, -0.2 };
        var root = ReferenceAlgebra.Cholesky(reference.Scale)!;
        double quad = ReferenceAlgebra.ForwardSolve(root, e).Sum(t => t * t) / row[0, 0];
        double expected = SpecialFunctions.LogGamma((df + 2) / 2) - SpecialFunctions.LogGamma(df / 2) - Math.Log(Math.PI * row[0, 0])
                          - 0.5 * ReferenceAlgebra.LogDeterminant(root) - (df + 2) / 2 * Math.Log(1 + quad);
        Assert.InRange(Math.Abs(expected - reference.LogDensity(z, value)), 0, 1e-12);
        var rng = new Xoshiro256PlusPlus(20261212);
        var moments = new RunningMoments[3];
        for (int sample = 0; sample < 30_000; sample++)
        {
            var a = reference.Sample(z, rng);
            moments[0].Add(a[0, 0]); moments[1].Add(a[0, 0] * a[0, 0]); moments[2].Add((a[0, 0] - mean[0, 0]) * (a[0, 1] - mean[0, 1]));
        }
        double[] truth = [mean[0, 0], mean[0, 0] * mean[0, 0] + row[0, 0] * reference.Scale[0, 0] / (reference.Nu - 3), row[0, 0] * reference.Scale[0, 1] / (reference.Nu - 3)];
        for (int j = 0; j < 3; j++) Assert.InRange(Math.Abs(moments[j].Mean - truth[j]) / Math.Sqrt(moments[j].Variance / moments[j].Count), 0, 4.5);
    }
}
