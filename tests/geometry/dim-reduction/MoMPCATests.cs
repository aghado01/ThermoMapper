using System;
using System.Linq;
using Maths.LinAlg;
using Maths.Geometry.DimReduction;
using Xunit;

namespace Maths.Geometry.Tests;

public sealed class MoMPCATests
{
    /// <summary>
    /// Three clean blocks share a mean ≈ (5, 0) and a high-variance axis ≈ (1, 0); two corrupted
    /// "nodes" sit at mean (5, 20) along the (1, 1) axis. The scale-calibrated product median on
    /// Rᵖ × Gr(r, p) resists the minority of corrupted blocks in BOTH factors — recovering the clean
    /// axis and a mean whose y-coordinate is pulled back toward 0, not toward the outliers' 20.
    /// </summary>
    [Fact]
    public void ComputeMoM_CorruptedBlocks_RecoverCleanMeanAndSubspace()
    {
        const int perBlock = 12;
        double[][] data = new double[5 * perBlock][];
        int idx = 0;
        idx = AddBlock(data, idx, perBlock, 5.0, 0.0, 0.0);    // clean
        idx = AddBlock(data, idx, perBlock, 5.2, 0.3, 0.04);   // clean
        idx = AddBlock(data, idx, perBlock, 4.8, -0.3, -0.04); // clean
        idx = AddBlock(data, idx, perBlock, 5.0, 20.0, 1.0);   // corrupted node
        idx = AddBlock(data, idx, perBlock, 5.0, 20.0, 1.0);   // corrupted node

        MoMPcaResult mom = MoMPCA.ComputeMoM(data, kBlocks: 5, numComponents: 1);

        double[] axis = mom.Basis[0];
        Assert.True(
            Math.Abs(axis[0]) > 0.9,
            $"Expected the robust axis near (1,0); |x|={Math.Abs(axis[0]):F4} (corrupted (1,1) would give ~0.707).");
        Assert.InRange(mom.Mean[1], -3.0, 3.0);   // robust mean resists the y=20 outlier nodes
    }

    /// <summary>
    /// Cut-locus regression guard mirroring DistributedSpredTests: the corrupted leading block
    /// spans xz while every clean block spans xy — principal angle exactly π/2 (YᵀZ singular),
    /// so a blockFrames[0] warm start would stall the preliminary Grassmann median on the
    /// corrupted subspace and poison the scale calibration and joint median downstream. The
    /// medoid warm start must recover the clean majority.
    /// </summary>
    [Fact]
    public void ComputeMoM_CorruptedFirstBlock_RecoversCleanSubspace()
    {
        const int perBlock = 24;
        double[][] data = new double[5 * perBlock][];
        int idx = 0;
        idx = AddCircleBlock(data, idx, perBlock, UnitX, UnitZ); // corrupted leading block: xz plane
        idx = AddCircleBlock(data, idx, perBlock, UnitX, UnitY); // clean majority: xy plane
        idx = AddCircleBlock(data, idx, perBlock, UnitX, UnitY);
        idx = AddCircleBlock(data, idx, perBlock, UnitX, UnitY);
        idx = AddCircleBlock(data, idx, perBlock, UnitY, UnitZ); // corrupted trailing block: yz plane

        MoMPcaResult mom = MoMPCA.ComputeMoM(data, kBlocks: 5, numComponents: 2);

        double toXy = SubspaceDistance(mom.Basis, [UnitX, UnitY]);

        Assert.InRange(toXy, 0.0, 0.1);
        Assert.True(toXy < SubspaceDistance(mom.Basis, [UnitX, UnitZ]));
        Assert.True(toXy < SubspaceDistance(mom.Basis, [UnitY, UnitZ]));
    }

    /// <summary>
    /// In-domain contrast: the corrupted leading plane tilted 0.3 rad off xz keeps its principal
    /// angles to the clean frames strictly below π/2, so Weiszfeld is not stalled — a
    /// blockFrames[0] warm start lands near tolerance (≈0.13) rather than hard-stalled at π/2 as
    /// in ComputeMoM_CorruptedFirstBlock_RecoversCleanSubspace. Full recovery within tolerance
    /// still requires the medoid warm start, whose preliminary median seeds the scale calibration
    /// and joint median cleanly.
    /// </summary>
    [Fact]
    public void ComputeMoM_TiltedCorruptedFirstBlock_RecoversCleanSubspace()
    {
        const int perBlock = 24;
        double[] tilted = [0.0, Math.Sin(0.3), Math.Cos(0.3)];
        double[][] data = new double[5 * perBlock][];
        int idx = 0;
        idx = AddCircleBlock(data, idx, perBlock, UnitX, tilted); // corrupted, 0.3 rad off xz
        idx = AddCircleBlock(data, idx, perBlock, UnitX, UnitY);  // clean majority: xy plane
        idx = AddCircleBlock(data, idx, perBlock, UnitX, UnitY);
        idx = AddCircleBlock(data, idx, perBlock, UnitX, UnitY);
        idx = AddCircleBlock(data, idx, perBlock, UnitY, UnitZ);  // corrupted trailing block: yz plane

        MoMPcaResult mom = MoMPCA.ComputeMoM(data, kBlocks: 5, numComponents: 2);

        double toXy = SubspaceDistance(mom.Basis, [UnitX, UnitY]);

        Assert.InRange(toXy, 0.0, 0.1);
        Assert.True(toXy < SubspaceDistance(mom.Basis, [UnitX, UnitZ]));
        Assert.True(toXy < SubspaceDistance(mom.Basis, [UnitY, UnitZ]));
    }

    /// <summary>
    /// The estimate is (μ, U) and the calibrated scale: an orthonormal basis of the subspace factor and α̂ strictly
    /// inside (0, 2), on the corrupted-blocks fixture.
    /// </summary>
    [Fact]
    public void ComputeMoM_ReturnsAnOrthonormalBasis_AndTheCalibratedScale()
    {
        const int perBlock = 24;
        double[] tilted = [0.0, Math.Sin(0.3), Math.Cos(0.3)];
        double[][] data = new double[5 * perBlock][];
        int idx = 0;
        idx = AddCircleBlock(data, idx, perBlock, UnitX, tilted);
        idx = AddCircleBlock(data, idx, perBlock, UnitX, UnitY);
        idx = AddCircleBlock(data, idx, perBlock, UnitX, UnitY);
        idx = AddCircleBlock(data, idx, perBlock, UnitX, UnitY);
        AddCircleBlock(data, idx, perBlock, UnitY, UnitZ);

        MoMPcaResult mom = MoMPCA.ComputeMoM(data, kBlocks: 5, numComponents: 2);

        Assert.Equal(3, mom.Mean.Length);
        Assert.Equal(2, mom.Basis.Length);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j <= i; j++)
            {
                double dot = 0.0;
                for (int t = 0; t < 3; t++) dot += mom.Basis[i][t] * mom.Basis[j][t];
                Assert.True(Math.Abs(dot - (i == j ? 1.0 : 0.0)) <= 1e-12, $"basis[{i}]·basis[{j}] = {dot:R}");
            }
        Assert.InRange(mom.Alpha, 1e-6, 2.0 - 1e-6);
        Assert.Equal(2, mom.Variances.Length);
        Assert.True(mom.Variances[0] >= mom.Variances[1] && mom.Variances[1] >= 0.0);
        foreach (double[] axis in mom.Basis)
        {
            int largest = 0;
            for (int t = 1; t < 3; t++) if (Math.Abs(axis[t]) > Math.Abs(axis[largest])) largest = t;
            Assert.True(axis[largest] > 0.0);
        }
    }

    /// <summary>
    /// Every block is an ellipse in the xy plane, so the robust subspace is xy. Three clean blocks have semi-axes near
    /// (2, 1), variances near (2, 0.5); two corrupted blocks stretch y to 20, variance near 200. Over all the data y
    /// carries the most variance, so plain PCA puts y first. The median of the blocks' projected covariances keeps the
    /// clean variances, so x comes first. The geometric median is pulled toward the outliers by about the clean blocks'
    /// own spread, not by the outliers' size: the clean x-variances spread over 1.9–2.1, which puts the y-variance near
    /// 0.56, where all the data give about 76.
    /// </summary>
    [Fact]
    public void ComputeMoM_InflatedVarianceBlocks_KeepTheCleanAxesAndVariances()
    {
        const int perBlock = 40;
        double[][] data = new double[5 * perBlock][];
        int idx = 0;
        idx = AddEllipseBlock(data, idx, perBlock, 2.00, 1.00);
        idx = AddEllipseBlock(data, idx, perBlock, 2.05, 1.02);
        idx = AddEllipseBlock(data, idx, perBlock, 1.95, 0.98);
        idx = AddEllipseBlock(data, idx, perBlock, 2.00, 20.0);
        AddEllipseBlock(data, idx, perBlock, 2.00, 19.0);

        PcaResult plain = Pca.Compute(data, numComponents: 2, center: true, whiten: false);
        Assert.True(Math.Abs(plain.Components[0][1]) > 0.99, "plain PCA should put y first here");

        MoMPcaResult mom = MoMPCA.ComputeMoM(data, kBlocks: 5, numComponents: 2);
        Assert.True(mom.Basis[0][0] > 0.999, $"first axis {string.Join(", ", mom.Basis[0])}");
        Assert.True(mom.Basis[1][1] > 0.999, $"second axis {string.Join(", ", mom.Basis[1])}");
        Assert.InRange(mom.Variances[0], 1.8, 2.2);
        Assert.InRange(mom.Variances[1], 0.4, 0.65);
        Assert.True(mom.Variances[1] < 0.01 * plain.Eigenvalues[0], $"robust {mom.Variances[1]:R} vs all-data {plain.Eigenvalues[0]:R}");
    }

    /// <summary>
    /// Rescaling the data rescales the robust mean and variances and keeps the subspace: the medians inside MoMPCA are
    /// scale-equivariant. Eight blocks in R⁴ with jittered means and tilted frames, two of them corrupted; α̂ ≈ 0.15 at
    /// scale 1, so the calibration stays inside its clamp [1e-6, 2 − 1e-6] at both scales tested, and the calibrated
    /// metrics differ only by a constant factor. The solver with absolute constants missed the 1e-3 case by 0.012 in
    /// the basis.
    /// </summary>
    [Theory]
    [InlineData(1e-3)]
    [InlineData(1e2)]
    public void ComputeMoM_RescaledData_GivesTheRescaledMeanAndTheSameSubspace(double s)
    {
        double[][] data = TiltedBlocks();
        double[][] scaled = data.Select(x => x.Select(v => s * v).ToArray()).ToArray();

        MoMPcaResult reference = MoMPCA.ComputeMoM(data, kBlocks: 8, numComponents: 2);
        MoMPcaResult rescaled = MoMPCA.ComputeMoM(scaled, kBlocks: 8, numComponents: 2);

        Assert.InRange(rescaled.Alpha, 1e-6 * (1 + 1e-9), (2.0 - 1e-6) * (1 - 1e-12));
        double balance = reference.Alpha / (2.0 - reference.Alpha);
        double rescaledBalance = rescaled.Alpha * s * s / (2.0 - rescaled.Alpha);
        Assert.True(Math.Abs(rescaledBalance / balance - 1.0) <= 1e-9, $"calibration balance {rescaledBalance:R} vs {balance:R}");
        for (int j = 0; j < 4; j++)
            Assert.True(Math.Abs(rescaled.Mean[j] / s - reference.Mean[j]) <= 1e-12 * Math.Abs(reference.Mean[j]) + 1e-12,
                $"mean[{j}]: {rescaled.Mean[j] / s:R} vs {reference.Mean[j]:R}");
        for (int c = 0; c < 2; c++)
        {
            for (int j = 0; j < 4; j++)
                Assert.True(Math.Abs(rescaled.Basis[c][j] - reference.Basis[c][j]) <= 1e-12, $"basis[{c}][{j}]");
            Assert.True(Math.Abs(rescaled.Variances[c] / (s * s) / reference.Variances[c] - 1.0) <= 1e-12, $"variance[{c}]");
        }
    }

    [Fact]
    public void ComputeMoM_SingleBlock_IsTheFullDataPca()
    {
        var rng = new Random(99);
        double[][] data = new double[24][];
        for (int i = 0; i < data.Length; i++)
            data[i] = new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble() };

        MoMPcaResult mom = MoMPCA.ComputeMoM(data, kBlocks: 1, numComponents: 2);
        PcaResult pca = Pca.Compute(data, numComponents: 2, center: true, whiten: false);

        Assert.Equal(pca.Mean, mom.Mean);
        Assert.Equal(pca.Components, mom.Basis);
        Assert.Equal(pca.Eigenvalues, mom.Variances);
        Assert.True(double.IsNaN(mom.Alpha));
    }

    // Eight blocks of 30 in R⁴: means jittered around (3, −2, 1, 5) by 0.3, frames tilted from span{e₁, e₂} by about
    // 0.12; blocks 7 and 8 corrupted (means by 4, frames by 0.8). Within a block, axis scales 2 and 1, noise 0.05.
    private static double[][] TiltedBlocks()
    {
        var rng = new Random(424242);
        const int p = 4, r = 2, perBlock = 30, blocks = 8;
        var grass = new GrassmannManifold(p, r);
        var u0 = new double[p * r];
        u0[0] = 1.0; u0[p + 1] = 1.0;
        var data = new double[blocks * perBlock][];
        for (int k = 0; k < blocks; k++)
        {
            bool corrupted = k >= 6;
            double[] mean = [3.0, -2.0, 1.0, 5.0];
            for (int j = 0; j < p; j++) mean[j] += (corrupted ? 4.0 : 0.3) * Gauss(rng);
            var tangent = new double[p * r];
            for (int t = 0; t < tangent.Length; t++) tangent[t] = (corrupted ? 0.8 : 0.12) * Gauss(rng);
            for (int col = 0; col < r; col++)          // horizontal at u0: zero the rows inside span{e₁, e₂}
                for (int row = 0; row < r; row++) tangent[col * p + row] = 0.0;
            var frame = new double[p * r];
            grass.ExpMap(u0, tangent, frame);
            for (int i = 0; i < perBlock; i++)
            {
                double a = 2.0 * Gauss(rng), b = Gauss(rng), e = 0.05 * Gauss(rng), f = 0.05 * Gauss(rng);
                var x = new double[p];
                for (int j = 0; j < p; j++) x[j] = mean[j] + a * frame[j] + b * frame[p + j];
                x[2] += e;
                x[3] += f;
                data[k * perBlock + i] = x;
            }
        }
        return data;
    }

    private static double Gauss(Random rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double[] UnitX => [1.0, 0.0, 0.0];
    private static double[] UnitY => [0.0, 1.0, 0.0];
    private static double[] UnitZ => [0.0, 0.0, 1.0];

    // A contiguous block of `count` points on the unit circle spanning {u, v} (orthonormal) in R³.
    private static int AddCircleBlock(double[][] data, int start, int count, double[] u, double[] v)
    {
        for (int i = 0; i < count; i++)
        {
            double t = 2.0 * Math.PI * i / count;
            var point = new double[3];
            for (int j = 0; j < 3; j++)
                point[j] = Math.Cos(t) * u[j] + Math.Sin(t) * v[j];
            data[start + i] = point;
        }
        return start + count;
    }

    // A contiguous block of `count` points on the ellipse with semi-axes a along x and b along y, centred at 0.
    private static int AddEllipseBlock(double[][] data, int start, int count, double a, double b)
    {
        for (int i = 0; i < count; i++)
        {
            double t = 2.0 * Math.PI * i / count;
            data[start + i] = [a * Math.Cos(t), b * Math.Sin(t), 0.0];
        }
        return start + count;
    }

    private static double SubspaceDistance(double[][] a, double[][] b)
    {
        var grass = new GrassmannManifold(ambientN: 3, subspaceR: 2);
        return grass.Distance(PackFrame(a), PackFrame(b));
    }

    // Rows (r vectors of length 3) → column-major 3×r Grassmann frame.
    private static double[] PackFrame(double[][] rows)
    {
        var frame = new double[3 * rows.Length];
        for (int c = 0; c < rows.Length; c++)
            Array.Copy(rows[c], 0, frame, c * 3, 3);
        return frame;
    }

    // A contiguous block of `count` points with mean (meanX, meanY), spread along (1, slope).
    private static int AddBlock(double[][] data, int start, int count, double meanX, double meanY, double slope)
    {
        double mid = (count - 1) / 2.0;
        for (int i = 0; i < count; i++)
        {
            double t = i - mid;
            data[start + i] = new[] { meanX + t, meanY + slope * t };
        }
        return start + count;
    }
}
