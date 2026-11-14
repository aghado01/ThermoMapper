using System;
using System.IO;
using System.Text.Json;
using Clustering.Statistical.GMM;
using Viz;
using Viz.Renderers;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The callers of <see cref="Maths.LinAlg.CholeskyDecomposition"/> on a covariance with no finite factor: GMM's component
/// cache and the Bhattacharyya distance refuse it, where they read wrong finite numbers before; the scene export draws a
/// singular block thin, as it draws the dimensions a layer lacks, and refuses an indefinite one (shared-code defects brief, slice 2).
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

    // The exported factor of one 3-D component, from the scene's JSON.
    private static double[] ExportedFactor(double[] covariance)
    {
        var layer = new GaussianLayer("g", new double[3], covariance, [1.0], 1, 3);
        var scene = new ScenePackage(new PointCloud(new double[3], 1, 3), [], [], [], [], [layer], [], null, 0, SceneRenderHints.Default);
        var output = new MemoryStream();
        new JsonExportRenderTarget().Render(scene, output);
        using var json = JsonDocument.Parse(output.ToArray());
        JsonElement factor = json.RootElement.GetProperty("gaussian_layers")[0].GetProperty("cholesky_l");
        var l = new double[9];
        for (int i = 0; i < 9; i++) l[i] = factor[i].GetDouble();
        return l;
    }

    // L·Lᵀ against Σ + ridge·I, entrywise relative to Σ's largest variance.
    private static double Gap(double[] l, double[] covariance, double ridge)
    {
        double largest = Math.Max(covariance[0], Math.Max(covariance[4], covariance[8])), gap = 0.0;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
            {
                double s = 0.0;
                for (int k = 0; k < 3; k++) s += l[r * 3 + k] * l[c * 3 + k];
                gap = Math.Max(gap, Math.Abs(s - covariance[r * 3 + c] - (r == c ? ridge : 0.0)) / largest);
            }
        return gap;
    }

    [Fact]
    public void TheExport_DrawsASingularBlockThin_AndRefusesAnIndefiniteOne()
    {
        // A block of rank two takes the missing dimensions' ε, relative to its largest variance: a thin ellipsoid.
        double[] flat = [2.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0];
        double[] l = ExportedFactor(flat);
        Assert.All(l, v => Assert.True(double.IsFinite(v)));
        Assert.True(Gap(l, flat, 2e-8) <= 1e-15, $"{Gap(l, flat, 2e-8):E2}");
        Assert.Equal(Math.Sqrt(2e-8), l[8], 15);

        // Rank one at 1e200 with one entry 2ε above it, which overflowed the floored factor: now a thin needle.
        double big = 1e200, off = big * (1 + 2 * Math.ScaleB(1.0, -52));
        double[] needle = [big, big, big, big, big, off, big, off, big];
        double[] n = ExportedFactor(needle);
        Assert.All(n, v => Assert.True(double.IsFinite(v)));
        Assert.True(Gap(n, needle, 1e-8 * big) <= 1e-12, $"{Gap(n, needle, 1e-8 * big):E2}");

        // A positive definite block exports its own factor, as before.
        double[] fine = [2.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 1.0];
        Assert.True(Gap(ExportedFactor(fine), fine, 0.0) <= 1e-15);

        // An indefinite block is no ellipsoid, at any ε.
        double[] indefinite = [1.0, 2.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        var refusal = Assert.Throws<InvalidOperationException>(() => ExportedFactor(indefinite));
        Assert.Contains("component 0", refusal.Message);
    }

    [Fact]
    public void AWarmStartWithASingularCovariance_IsRefused()
    {
        var gmm = new GaussianMixtureModel(1, 2);
        Assert.Throws<InvalidOperationException>(() => gmm.InitializeWithParameters([[0.0, 0.0]], [new double[,] { { 1, 1 }, { 1, 1 } }], [1.0]));
        gmm.InitializeWithParameters([[0.0, 0.0]], [new double[,] { { 1, 0.5 }, { 0.5, 1 } }], [1.0]);
        Assert.True(double.IsFinite(gmm.Components[0].LogNormalizationFactor));
    }
}
