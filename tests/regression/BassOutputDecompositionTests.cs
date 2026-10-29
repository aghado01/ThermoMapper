using System;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's response decomposition (stage 5, O1; D07) against independent computations: the training reconstruction
/// y = y_m + y_s ∘ (K W + e) to rounding, the residual's orientation on an asymmetric fixture, R's convention (K = U_r D_r,
/// W = V_rᵀ with unit rows), the rank rule (≥ the share, capped at the numerical rank, full retention at ρ), constant
/// columns, the sign convention, leave-one-case-out residuals against explicit refits in output space under both
/// centres, the edge cases (n = 2, r = n − 1, full rank with m &lt; n, one location), a supplied nonorthogonal basis, the
/// truncation covariance against a dense computation, and the refusals.
/// </summary>
public sealed class BassOutputDecompositionTests
{
    private readonly ITestOutputHelper _out;
    public BassOutputDecompositionTests(ITestOutputHelper output) => _out = output;

    // n curves at m locations whose phase and amplitude move with two inputs, with positive spikes on every third case
    // so the residuals are skewed; z holds the inputs.
    internal static double[,] Curves(int n, int m, int seed, out double[,] z, bool spikes = true)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        z = new double[n, 2];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            z[i, 0] = rng.NextDouble();
            z[i, 1] = rng.NextDouble();
            for (int l = 0; l < m; l++)
            {
                double t = m == 1 ? 0.5 : l / (m - 1.0);
                y[i, l] = Math.Sin(2 * Math.PI * (t - 0.3 * z[i, 0])) * (1 + z[i, 1]) + 0.5 * z[i, 1] * t * t
                        + (spikes && i % 3 == 0 ? 2 * Math.Exp(-Math.Pow((t - 0.2 - 0.07 * i) / 0.04, 2)) : 0.0)
                        + 0.01 * Normal.SampleStandard(rng);
            }
        }
        return y;
    }

    private static double[,] Curves(int n, int m, int seed = 3, bool spikes = true) => Curves(n, m, seed, out _, spikes);

    private static double MaxAbs(double[,] a) => a.Cast<double>().Max(Math.Abs);

    // The training reconstruction with the in-sample residual.
    private static double ReconstructionError(BassDecomposition d, double[,] y)
    {
        int n = d.Cases, m = d.Locations;
        double worst = 0.0;
        var w = new double[d.Rank];
        var e = new double[m];
        var r = new double[m];
        for (int i = 0; i < n; i++)
        {
            for (int k = 0; k < d.Rank; k++) w[k] = d.Coefficient(k, i);
            for (int l = 0; l < m; l++) e[l] = d.Truncation.InSample![i, l];
            d.Map.Reconstruct(w, e, r);
            for (int l = 0; l < m; l++) worst = Math.Max(worst, Math.Abs(r[l] - y[i, l]));
        }
        return worst;
    }

    // Leave-one-case-out by explicit refits in output space: the others centred by their mean (or the fixed centre),
    // their SVD, the rank min(r, the count above the full data's tolerance), case i projected off it.
    private static double[,] ExplicitLeaveOneOut(BassDecomposition d, double[,] y, double tolerance, out int[] ranks, bool fixedCentre = false)
    {
        int n = d.Cases, m = d.Locations;
        var e = new double[n, m];
        ranks = new int[n];
        for (int i = 0; i < n; i++)
        {
            var centre = new double[m];
            for (int l = 0; l < m; l++)
                centre[l] = fixedCentre ? d.Map.Offsets[l] : Enumerable.Range(0, n).Where(j => j != i).Average(j => y[j, l]);
            var others = new double[m, n - 1];
            for (int j = 0, c = 0; j < n; j++)
            {
                if (j == i) continue;
                for (int l = 0; l < m; l++) others[l, c] = (y[j, l] - centre[l]) / d.Map.Scales[l];
                c++;
            }
            ThinSvdResult svd = ThinSvd.Decompose(others);
            int rank = Math.Min(d.Rank, svd.SingularValues.Count(s => s > tolerance));
            ranks[i] = rank;
            var x = new double[m];
            for (int l = 0; l < m; l++) x[l] = (y[i, l] - centre[l]) / d.Map.Scales[l];
            for (int q = 0; q < rank; q++)
            {
                double p = 0.0;
                for (int l = 0; l < m; l++) p += svd.LeftVectors[q][l] * x[l];
                for (int l = 0; l < m; l++) x[l] -= p * svd.LeftVectors[q][l];
            }
            for (int l = 0; l < m; l++) e[i, l] = x[l];
        }
        return e;
    }

    [Fact]
    public void Pca_ReconstructsTheTrainingResponses_InRsConvention()
    {
        double[,] y = Curves(8, 40);
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.9, ScaleByColumnSd = true });
        _out.WriteLine($"r = {d.Rank} of ρ = {d.Pca!.NumericalRank}, share {d.Pca.ShareReached:F4}, reconstruction {ReconstructionError(d, y):E1}");
        Assert.InRange(ReconstructionError(d, y), 0, 1e-13 * MaxAbs(y));
        // W = V_rᵀ: unit, orthogonal rows; K = U_r D_r: orthogonal columns of norm d_k.
        for (int k = 0; k < d.Rank; k++)
        {
            for (int q = 0; q <= k; q++)
            {
                double ww = Enumerable.Range(0, d.Cases).Sum(i => d.Coefficient(k, i) * d.Coefficient(q, i));
                double kk = Enumerable.Range(0, d.Locations).Sum(l => d.Map.Basis(l, k) * d.Map.Basis(l, q));
                Assert.Equal(k == q ? 1.0 : 0.0, ww, 12);
                Assert.Equal(k == q ? d.Pca.SingularValues[k] * d.Pca.SingularValues[k] : 0.0, kk, 10);
            }
            // Each row of W sums to zero: the centred responses have the constant in their null space.
            Assert.Equal(0.0, Enumerable.Range(0, d.Cases).Sum(i => d.Coefficient(k, i)), 12);
            // The sign convention: each component's largest entry of K is positive.
            double largest = Enumerable.Range(0, d.Locations).Select(l => d.Map.Basis(l, k)).MaxBy(Math.Abs);
            Assert.True(largest > 0);
        }
        // The scale is each column's sd (n − 1 denominator).
        for (int l = 0; l < d.Locations; l++)
        {
            double mean = Enumerable.Range(0, 8).Average(i => y[i, l]);
            double sd = Math.Sqrt(Enumerable.Range(0, 8).Sum(i => (y[i, l] - mean) * (y[i, l] - mean)) / 7);
            Assert.Equal(mean, d.Map.Offsets[l], 13);
            Assert.Equal(sd, d.Map.Scales[l], 13);
        }
    }

    [Fact]
    public void Residual_IsObservedMinusReconstructed_InStandardizedUnits()
    {
        // Positive spikes the leading components miss: observed minus reconstructed is skewed to the right.
        double[,] y = Curves(8, 40);
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.9, ScaleByColumnSd = true });
        double[,] e = d.Truncation.InSample!;
        int n = d.Cases, m = d.Locations;
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++)
            {
                double reconstructed = Enumerable.Range(0, d.Rank).Sum(k => d.Map.Basis(l, k) * d.Coefficient(k, i));
                Assert.Equal((y[i, l] - d.Map.Offsets[l]) / d.Map.Scales[l] - reconstructed, e[i, l], 13);
            }
        double[] all = e.Cast<double>().ToArray();
        double mean = all.Average();
        double skew = all.Average(v => Math.Pow(v - mean, 3)) / Math.Pow(all.Average(v => (v - mean) * (v - mean)), 1.5);
        _out.WriteLine($"skewness of e: {skew:F2}");
        Assert.True(skew > 1.0);
    }

    [Fact]
    public void RankRule_TakesTheSmallestRankReachingTheShare()
    {
        double[,] y = Curves(12, 30);
        foreach (double share in new[] { 0.5, 0.9, 0.99, 0.999 })
        {
            BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = share });
            double[] s = d.Pca!.SingularValues;
            double total = s.Sum(v => v * v);
            double Share(int r) => s.Take(r).Sum(v => v * v) / total;
            Assert.True(Share(d.Rank) >= share);
            Assert.True(d.Rank == 1 || Share(d.Rank - 1) < share);
            Assert.Equal(Share(d.Rank), d.Pca.ShareReached, 14);
        }
        // An explicit rank overrides the share, capped at ρ.
        BassDecomposition three = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 3 });
        Assert.Equal(3, three.Rank);
        Assert.Null(three.Pca!.VarianceShare);
        BassDecomposition many = BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 50 });
        Assert.Equal(11, many.Rank);
        Assert.True(many.Pca!.RankCapped);
    }

    [Fact]
    public void FullRetention_IsTheNumericalRank_WithZeroInSampleTruncation()
    {
        double[,] y = Curves(5, 40);
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 1.0 });
        Assert.Equal(4, d.Pca!.NumericalRank);   // the centred rank, n − 1
        Assert.Equal(4, d.Rank);
        Assert.InRange(MaxAbs(d.Truncation.InSample!), 0, d.Pca.RankTolerance);
        Assert.Null(d.Truncation.StabilityRatio);
        // The leave-one-out ranks are capped at the other cases' rank, n − 2, and their residuals are not zero.
        Assert.All(d.Truncation.LeaveOneOutRanks!, r => Assert.Equal(3, r));
        Assert.True(d.Truncation.LeaveOneOutRms > 0.01);
        Assert.Equal(TruncationResidual.LeaveOneCaseOut, d.Truncation.Kind);
        Assert.Same(d.Truncation.LeaveOneOut, d.Truncation.Values);
    }

    [Fact]
    public void FullRank_WithFewerLocationsThanCases()
    {
        // m < n: the centred responses have rank m, and full retention reconstructs every case.
        double[,] y = Curves(9, 3);
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 1.0 });
        Assert.Equal(3, d.Pca!.NumericalRank);
        Assert.Equal(3, d.Rank);
        Assert.InRange(MaxAbs(d.Truncation.InSample!), 0, d.Pca.RankTolerance);
        // Eight other cases in three dimensions keep rank 3, so the leave-one-out residuals vanish too.
        Assert.All(d.Truncation.LeaveOneOutRanks!, r => Assert.Equal(3, r));
        Assert.InRange(MaxAbs(d.Truncation.LeaveOneOut!), 0, 1e-12);
    }

    [Fact]
    public void ConstantColumn_IsCentred_AndScaledByOne()
    {
        double[,] y = Curves(8, 20);
        for (int i = 0; i < 8; i++) y[i, 4] = 7.0;
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.95, ScaleByColumnSd = true });
        Assert.Equal([4], d.Pca!.ConstantColumns);
        Assert.Equal(1.0, d.Map.Scales[4]);
        Assert.Equal(7.0, d.Map.Offsets[4]);
        for (int k = 0; k < d.Rank; k++) Assert.InRange(Math.Abs(d.Map.Basis(4, k)), 0, 1e-14);
        for (int i = 0; i < 8; i++)
        {
            Assert.InRange(Math.Abs(d.Truncation.InSample![i, 4]), 0, 1e-14);
            Assert.InRange(Math.Abs(d.Truncation.LeaveOneOut![i, 4]), 0, 1e-14);
        }
        Assert.InRange(ReconstructionError(d, y), 0, 1e-13 * MaxAbs(y));
    }

    [Theory]
    [InlineData(5, 40, 0.99, false)]
    [InlineData(10, 40, 0.99, false)]
    [InlineData(10, 40, 1.0, false)]
    [InlineData(12, 6, 0.9, false)]
    [InlineData(8, 30, 0.95, true)]
    public void LeaveOneOut_EqualsExplicitRefits(int n, int m, double share, bool scaled)
    {
        double[,] y = Curves(n, m, seed: n + m);
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = share, ScaleByColumnSd = scaled });
        double[,] reference = ExplicitLeaveOneOut(d, y, d.Pca!.RankTolerance, out int[] ranks);
        double gap = 0.0;
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++) gap = Math.Max(gap, Math.Abs(reference[i, l] - d.Truncation.LeaveOneOut![i, l]));
        _out.WriteLine($"n = {n}, m = {m}, r = {d.Rank}: gap {gap:E1}, ratio {d.Truncation.StabilityRatio:F3}, ranks {string.Join(" ", ranks)}");
        Assert.Equal(ranks, d.Truncation.LeaveOneOutRanks);
        Assert.InRange(gap, 0, 1e-11 * MaxAbs(d.Truncation.LeaveOneOut!));
    }

    [Fact]
    public void LeaveOneOut_UnderAFixedCentre_EqualsExplicitRefits()
    {
        double[,] y = Curves(7, 25);
        var centre = Enumerable.Range(0, 25).Select(l => 0.1 * l).ToArray();
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.95, Centre = centre });
        double[,] reference = ExplicitLeaveOneOut(d, y, d.Pca!.RankTolerance, out int[] ranks, fixedCentre: true);
        Assert.Equal(ranks, d.Truncation.LeaveOneOutRanks);
        for (int i = 0; i < 7; i++)
            for (int l = 0; l < 25; l++) Assert.Equal(reference[i, l], d.Truncation.LeaveOneOut![i, l], 11);
        Assert.InRange(ReconstructionError(d, y), 0, 1e-13 * MaxAbs(y));
    }

    [Fact]
    public void TwoCases_LeaveOneOutIsTheirDifference()
    {
        double[,] y = Curves(2, 30);
        BassDecomposition d = BassDecomposition.FromPca(y);
        Assert.Equal(1, d.Rank);
        Assert.Equal([0, 0], d.Truncation.LeaveOneOutRanks!);
        for (int l = 0; l < 30; l++)
        {
            Assert.Equal(y[0, l] - y[1, l], d.Truncation.LeaveOneOut![0, l], 13);
            Assert.Equal(y[1, l] - y[0, l], d.Truncation.LeaveOneOut![1, l], 13);
        }
    }

    [Fact]
    public void OneLocation_DecomposesToOneComponent()
    {
        double[,] y = Curves(6, 1);
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { Residual = TruncationResidual.InSample });
        Assert.Equal(1, d.Rank);
        Assert.Equal(TruncationResidual.InSample, d.Truncation.Kind);
        Assert.InRange(MaxAbs(d.Truncation.Values), 0, 1e-14);
        Assert.InRange(ReconstructionError(d, y), 0, 1e-14);
    }

    [Fact]
    public void StabilityRatio_IsTheRmsRatio()
    {
        double[,] y = Curves(10, 40);
        BassDecomposition d = BassDecomposition.FromPca(y);
        double Rms(double[,] e) => Math.Sqrt(e.Cast<double>().Average(v => v * v));
        Assert.Equal(Rms(d.Truncation.LeaveOneOut!) / Rms(d.Truncation.InSample!), d.Truncation.StabilityRatio!.Value, 12);
        Assert.True(d.Truncation.StabilityRatio > 1);
    }

    [Fact]
    public void SuppliedBasis_Nonorthogonal_ReconstructsWithItsResiduals()
    {
        int n = 6, m = 12, r = 2;
        var rng = new Xoshiro256PlusPlus(4);
        var k = new double[m, r];
        for (int l = 0; l < m; l++)
        {
            k[l, 0] = 1.0;
            k[l, 1] = l / (m - 1.0) + 0.3;   // not orthogonal to the first
        }
        var w = new double[r, n];
        for (int i = 0; i < n; i++)
        {
            w[0, i] = Normal.SampleStandard(rng);
            w[1, i] = Normal.SampleStandard(rng);
        }
        double[] offsets = [.. Enumerable.Range(0, m).Select(l => 0.5 * l)];
        double[] scales = [.. Enumerable.Range(0, m).Select(l => 1.0 + 0.1 * l)];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++) y[i, l] = offsets[l] + scales[l] * (k[l, 0] * w[0, i] + k[l, 1] * w[1, i] + 0.05 * Normal.SampleStandard(rng));
        BassDecomposition d = BassDecomposition.FromBasis(k, w, offsets, scales, responses: y);
        Assert.Null(d.Pca);
        Assert.Null(d.Truncation.LeaveOneOut);
        Assert.Null(d.Truncation.StabilityRatio);
        Assert.Equal(TruncationResidual.InSample, d.Truncation.Kind);
        Assert.InRange(ReconstructionError(d, y), 0, 1e-13 * MaxAbs(y));
        Assert.Equal(w[1, 3], d.CoefficientResponse(1)[3]);

        // Residuals as supplied.
        var given = new double[n, m];
        given[2, 5] = 0.25;
        BassDecomposition s = BassDecomposition.FromBasis(k, w, offsets, scales, residuals: given);
        Assert.Equal(TruncationResidual.Supplied, s.Truncation.Kind);
        Assert.Equal(0.25, s.Truncation.Values[2, 5]);
        Assert.Null(s.Truncation.InSample);
    }

    [Fact]
    public void TruncationCovariance_MatchesADenseComputation()
    {
        double[,] y = Curves(9, 15);
        BassDecomposition d = BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0.9, ScaleByColumnSd = true });
        double[,] covariance = d.TruncationCovariance();
        // Dense: the residuals in physical units, y_s ∘ eᵢ, and their covariance under a uniform draw over the cases.
        int n = 9, m = 15;
        var physical = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++) physical[i, l] = d.Map.Scales[l] * d.Truncation.Values[i, l];
        for (int l = 0; l < m; l++)
            for (int j = 0; j < m; j++)
            {
                double ml = Enumerable.Range(0, n).Average(i => physical[i, l]), mj = Enumerable.Range(0, n).Average(i => physical[i, j]);
                double dense = Enumerable.Range(0, n).Average(i => (physical[i, l] - ml) * (physical[i, j] - mj));
                Assert.Equal(dense, covariance[l, j], 13);
            }
    }

    [Fact]
    public void Refusals_SayWhy()
    {
        double[,] y = Curves(6, 10);
        Assert.Contains("two cases", Assert.Throws<ArgumentException>(() => BassDecomposition.FromPca(Curves(1, 10))).Message);
        Assert.Throws<ArgumentOutOfRangeException>(() => BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => BassDecomposition.FromPca(y, new BassPcaSpec { VarianceShare = 1.01 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 0 }));
        Assert.Contains("not both", Assert.Throws<ArgumentException>(() =>
            BassDecomposition.FromPca(y, new BassPcaSpec { ScaleByColumnSd = true, Scale = new double[10] })).Message);
        Assert.Contains("positive", Assert.Throws<ArgumentException>(() =>
            BassDecomposition.FromPca(y, new BassPcaSpec { Scale = new double[10] })).Message);
        Assert.Contains("one value per location", Assert.Throws<ArgumentException>(() =>
            BassDecomposition.FromPca(y, new BassPcaSpec { Centre = new double[3] })).Message);
        Assert.Contains("builds its own", Assert.Throws<ArgumentException>(() =>
            BassDecomposition.FromPca(y, new BassPcaSpec { Residual = TruncationResidual.Supplied })).Message);
        var flat = new double[4, 5];
        Assert.Contains("no component", Assert.Throws<ArgumentException>(() => BassDecomposition.FromPca(flat)).Message);
        y[2, 3] = double.NaN;
        Assert.Contains("non-finite", Assert.Throws<ArgumentException>(() => BassDecomposition.FromPca(y)).Message);

        var k = new double[10, 1];
        k[0, 0] = 1;
        var w = new double[1, 6];
        Assert.Contains("exactly one", Assert.Throws<ArgumentException>(() =>
            BassDecomposition.FromBasis(k, w, new double[10], Enumerable.Repeat(1.0, 10).ToArray())).Message);
        Assert.Contains("6 × 10", Assert.Throws<ArgumentException>(() =>
            BassDecomposition.FromBasis(k, w, new double[10], Enumerable.Repeat(1.0, 10).ToArray(), responses: new double[5, 10])).Message);
        Assert.Contains("one row per component", Assert.Throws<ArgumentException>(() =>
            BassDecomposition.FromBasis(k, new double[2, 6], new double[10], Enumerable.Repeat(1.0, 10).ToArray(), responses: new double[6, 10])).Message);
    }
}
