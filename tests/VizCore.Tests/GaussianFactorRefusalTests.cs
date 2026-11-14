using System;
using System.IO;
using System.Linq;
using System.Text.Json;
using Clustering.Statistical.GMM;
using Viz;
using Viz.Renderers;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The callers of a Gaussian's covariance factor (<see cref="Maths.LinAlg.DenseCholesky"/>) on a covariance with none: GMM's component
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

        // The same block at 1e-8 of the scale stays as thin: its ε is relative below unit scale too.
        double[] tiny = [.. flat.Select(v => 1e-8 * v)];
        double[] t = ExportedFactor(tiny);
        Assert.True(Gap(t, tiny, 2e-16) <= 1e-15, $"{Gap(t, tiny, 2e-16):E2}");
        Assert.True(Math.Abs(t[8] - Math.Sqrt(2e-16)) <= 1e-15 * Math.Sqrt(2e-16), $"{t[8]:R}");

        // An all-zero block takes the missing dimensions' ε as it stands: a sphere of radius 1e-4.
        double[] zero = ExportedFactor(new double[9]);
        Assert.Equal(1e-4, zero[0], 15);
        Assert.Equal(1e-4, zero[8], 15);

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

    // ── The components on DenseCholesky ──────────────────────────────────────────────────────────────────────────────

    // A positive definite covariance with off-diagonal entries: B·Bᵀ + I/2 for a seeded B.
    private static double[,] Spd(int n, int seed)
    {
        var rng = new Random(seed);
        var b = new double[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) b[i, j] = rng.NextDouble() - 0.5;
        var a = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                for (int k = 0; k < n; k++) a[i, j] += b[i, k] * b[j, k];
                if (i == j) a[i, j] += 0.5;
            }
        return a;
    }

    // The retired CholeskyDecomposition's factor, inverse and sample, verbatim.
    private static (double[,] L, double[,] Inverse, double LogDet) EarlierFactor(double[,] sigma)
    {
        int n = sigma.GetLength(0);
        var l = new double[n, n];
        var lInv = new double[n, n];
        double logDetL = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = sigma[i, j];
                for (int k = 0; k < j; k++) sum -= l[i, k] * l[j, k];
                if (i == j) { l[i, i] = Math.Sqrt(sum); logDetL += Math.Log(l[i, i]); }
                else l[i, j] = sum / l[j, j];
            }
        for (int i = 0; i < n; i++)
        {
            lInv[i, i] = 1.0 / l[i, i];
            for (int j = 0; j < i; j++)
            {
                double sum = 0.0;
                for (int k = j; k < i; k++) sum -= l[i, k] * lInv[k, j];
                lInv[i, j] = sum / l[i, i];
            }
        }
        var inverse = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = i; j < n; j++)
            {
                double sum = 0.0;
                for (int k = j; k < n; k++) sum += lInv[k, i] * lInv[k, j];
                inverse[i, j] = inverse[j, i] = sum;
            }
        return (l, inverse, 2.0 * logDetL);
    }

    private static double[] EarlierSample(double[,] l, Random rng, double[] mean)
    {
        int dim = mean.Length;
        double[] z = new double[dim];
        for (int i = 0; i < dim - 1; i += 2)
        {
            double u1 = 1.0 - rng.NextDouble(), u2 = 1.0 - rng.NextDouble();
            double mag = Math.Sqrt(-2.0 * Math.Log(u1));
            z[i] = mag * Math.Cos(2.0 * Math.PI * u2);
            z[i + 1] = mag * Math.Sin(2.0 * Math.PI * u2);
        }
        if (dim % 2 == 1)
        {
            double u1 = 1.0 - rng.NextDouble(), u2 = 1.0 - rng.NextDouble();
            z[dim - 1] = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
        double[] x = new double[dim];
        for (int i = 0; i < dim; i++)
        {
            x[i] = mean[i];
            for (int j = 0; j <= i; j++) x[i] += l[i, j] * z[j];
        }
        return x;
    }

    [Fact]
    [Trait("Tier", "Pinned")]
    public void AComponentsCacheAndSample_AreTheRetiredClassesBitForBit()
    {
        foreach (int n in new[] { 1, 2, 3, 7, 12 })
        {
            double[,] sigma = Spd(n, 40 + n);
            GaussianComponent c = Component(sigma);
            for (int i = 0; i < n; i++) c.Mean[i] = 0.25 * i - 1.0;
            (double[,] l, double[,] inverse, double logDet) = EarlierFactor(sigma);
            Assert.Equal(-0.5 * (n * Math.Log(2.0 * Math.PI) + logDet), c.LogNormalizationFactor);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) Assert.Equal(BitConverter.DoubleToInt64Bits(inverse[i, j]), BitConverter.DoubleToInt64Bits(c.CovarianceInverse[i, j]));
            for (int s = 0; s < 5; s++)
                Assert.Equal(Array.ConvertAll(EarlierSample(l, new Random(s), c.Mean), BitConverter.DoubleToInt64Bits),
                             Array.ConvertAll(c.Sample(new Random(s)), BitConverter.DoubleToInt64Bits));
        }
    }

    [Fact]
    public void AComponentWhoseInverseOverflows_IsRefused()
    {
        // [1e-320] factors, but its precision passes double.MaxValue: the component is left unusable, not half-published.
        GaussianComponent c = Component(new double[,] { { 1e-320 } }, update: false);
        Assert.Contains("inverse", Assert.Throws<InvalidOperationException>(c.UpdateCache).Message);
        Assert.True(double.IsNaN(c.LogNormalizationFactor));
        Assert.True(double.IsNaN(c.CovarianceInverse[0, 0]));
    }

    [Fact]
    public void TheDistance_ReadsItsMeanTermThroughTheFactor()
    {
        foreach (int n in new[] { 1, 3, 8 })
        {
            double[,] c1 = Spd(n, 70 + n), c2 = Spd(n, 90 + n);
            double[] m1 = new double[n], m2 = new double[n];
            for (int i = 0; i < n; i++) { m1[i] = 0.3 * i; m2[i] = 1.0 - 0.2 * i; }
            var bar = new double[n, n];
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) bar[i, j] = 0.5 * (c1[i, j] + c2[i, j]);
            (_, double[,] inverse, double logDetBar) = EarlierFactor(bar);
            double mean = 0.0;
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) mean += (m1[i] - m2[i]) * inverse[i, j] * (m1[j] - m2[j]);
            double expected = mean / 8.0 + 0.5 * (logDetBar - 0.5 * (EarlierFactor(c1).LogDet + EarlierFactor(c2).LogDet));
            double distance = BhattacharyyaCoefficient.Distance(m1, c1, m2, c2);
            Assert.True(Math.Abs(distance - expected) <= 1e-13 * Math.Max(1.0, Math.Abs(expected)), $"{distance:R} against {expected:R}");
            GaussianComponent a = Component(c1), b = Component(c2);
            Array.Copy(m1, a.Mean, n);
            Array.Copy(m2, b.Mean, n);
            Assert.True(Math.Abs(BhattacharyyaCoefficient.Distance(a, b) - expected) <= 1e-13 * Math.Max(1.0, Math.Abs(expected)));
        }
    }
}
