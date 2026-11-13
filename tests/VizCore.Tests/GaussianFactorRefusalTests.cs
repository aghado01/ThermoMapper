using System;
using System.IO;
using Clustering.Statistical.GMM;
using Viz;
using Viz.Renderers;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The callers of <see cref="Maths.LinAlg.CholeskyDecomposition"/> on a covariance with no finite factor: GMM's component
/// cache, the Bhattacharyya distance and the scene export's ellipsoids each refuse it, where they read wrong finite
/// numbers or exported a meaningless factor before (shared-code defects brief, slice 2).
/// </summary>
public sealed class GaussianFactorRefusalTests
{
    // A squared-exponential covariance on n points of [0, 1]: its spectrum falls below rounding, and at length 0.1 on 50
    // points the floored factor overflows.
    private static double[,] Kernel(int n, double length)
    {
        var k = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double r = (i - j) / (double)(n - 1) / length;
                k[i, j] = k[j, i] = Math.Exp(-0.5 * r * r);
            }
        return k;
    }

    private static double[,] Identity(int n)
    {
        var a = new double[n, n];
        for (int i = 0; i < n; i++) a[i, i] = 1.0;
        return a;
    }

    private static GaussianComponent Component(double[,] covariance, bool update = true)
    {
        int n = covariance.GetLength(0);
        var c = new GaussianComponent(n) { Weight = 0.5 };
        Array.Copy(covariance, c.Covariance, covariance.Length);
        if (update) c.UpdateCache();
        return c;
    }

    [Fact]
    public void AComponentWithNoFiniteFactor_RefusesItsCache_AndFactorsOnceRegularized()
    {
        const int n = 50;
        GaussianComponent c = Component(Kernel(n, 0.1), update: false);
        Assert.Throws<InvalidOperationException>(c.UpdateCache);
        Assert.True(double.IsNaN(c.LogNormalizationFactor));
        Assert.True(double.IsNaN(c.EvaluateLogPdf(new double[n])));
        Assert.Throws<InvalidOperationException>(() => c.Sample(new Random(1)));

        // The M-step's ridge lifts every pivot above the floor.
        for (int i = 0; i < n; i++) c.Covariance[i, i] += 1e-6;
        c.UpdateCache();
        GaussianComponent fresh = Component(c.Covariance);
        Assert.Equal(fresh.LogNormalizationFactor, c.LogNormalizationFactor);
        Assert.Equal(fresh.CovarianceInverse, c.CovarianceInverse);
        Assert.True(double.IsFinite(c.EvaluateLogPdf(new double[n])));
    }

    [Fact]
    public void TheDistance_RefusesACovarianceWithNoFiniteFactor()
    {
        const int n = 50;
        double[,] kernel = Kernel(n, 0.1), identity = Identity(n);
        var mean = new double[n];
        // The average factors; the second covariance does not, and its LogDet read a finite −1156 before.
        Assert.Equal("cov2", Assert.Throws<ArgumentException>(() => BhattacharyyaCoefficient.Distance(mean, identity, mean, kernel)).ParamName);
        Assert.Equal("cov1", Assert.Throws<ArgumentException>(() => BhattacharyyaCoefficient.Distance(mean, kernel, mean, identity)).ParamName);
        Assert.Throws<ArgumentException>(() => BhattacharyyaCoefficient.Distance(mean, kernel, mean, kernel));

        GaussianComponent good = Component(identity), bad = Component(Kernel(n, 0.1), update: false);
        Assert.Throws<InvalidOperationException>(bad.UpdateCache);
        Assert.Equal("b", Assert.Throws<ArgumentException>(() => BhattacharyyaCoefficient.Distance(good, bad)).ParamName);
        Assert.Equal(0.0, BhattacharyyaCoefficient.Distance(good, good), 12);
    }

    [Fact]
    public void TheExport_RefusesAnEllipsoidWithNoFiniteFactor()
    {
        // Rank one at 1e200 with one entry 2ε above it: L stays finite, but its last pivot overflows to −∞ under the floor.
        double big = 1e200, off = big * (1 + 2 * Math.ScaleB(1.0, -52));
        double[] covariance = [big, big, big, big, big, off, big, off, big];
        var layer = new GaussianLayer("overflow", new double[3], covariance, [1.0], 1, 3);
        var scene = new ScenePackage(new PointCloud(new double[3], 1, 3), [], [], [], [], [layer], [], null, 0, SceneRenderHints.Default);
        var refusal = Assert.Throws<InvalidOperationException>(() => new JsonExportRenderTarget().Render(scene, new MemoryStream()));
        Assert.Contains("component 0", refusal.Message);

        var fine = new GaussianLayer("fine", new double[3], [2.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 1.0], [1.0], 1, 3);
        var ok = new ScenePackage(new PointCloud(new double[3], 1, 3), [], [], [], [], [fine], [], null, 0, SceneRenderHints.Default);
        var output = new MemoryStream();
        new JsonExportRenderTarget().Render(ok, output);
        Assert.True(output.Length > 0);
    }
}
