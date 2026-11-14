using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Regression.Subspaces;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's draw weights over their whole finite range. <c>[1, 1]</c>, <c>[1e308, 1e308]</c> and <c>[5e307, 5e307]</c> are one
/// posterior, so every reader of <see cref="BassFit.DrawWeights"/> (and of the subspace posteriors' weights) gives them one
/// reading: the selection, the summaries, the Sobol means, the resamplers, the residual check's posterior mean, and the
/// projection report's term count. A power-of-two multiple reads bit for bit, ordinary weights keep their earlier bits, and
/// zero weights drop out as they did (shared-code defects brief, slice 3).
/// </summary>
public sealed class DrawWeightRangeTests(ITestOutputHelper output)
{
    // The ordinary row, the row whose sum overflows, the row whose sum is finite but whose weighted sums overflow, and exact
    // power-of-two multiples at both ends. Subnormal weights cannot hold a pattern's thirds exactly, so they appear only where
    // their ratios are exact.
    private static readonly double[] Scales = [1.0, 1e308, 5e307, Math.ScaleB(1.0, 1000), Math.ScaleB(1.0, -1000)];

    private static bool PowerOfTwo(double scale) => (BitConverter.DoubleToInt64Bits(scale) & 0xFFFFFFFFFFFFFL) == 0 && double.IsNormal(scale);

    private static void AssertSame(double expected, double actual, double scale, string what)
    {
        if (PowerOfTwo(scale))
            Assert.True(BitConverter.DoubleToInt64Bits(expected) == BitConverter.DoubleToInt64Bits(actual), $"{what} at {scale:R}: {actual:R} against {expected:R}.");
        else
            Assert.True(Math.Abs(actual - expected) <= 1e-11 * Math.Max(1.0, Math.Abs(expected)), $"{what} at {scale:R}: {actual:R} against {expected:R}.");
    }

    // ── The selection and the summary, the readers the review measured ──────────────────────────────────────────────

    [Theory]
    [InlineData(1e308)]
    [InlineData(5e307)]
    [InlineData(1e-320)]
    public void TheSelection_ReadsEqualWeightsAsEqual(double scale)
    {
        (_, double[] weights) = Draws3.Select(null, 2, [scale, scale]);
        Assert.Equal([0.5, 0.5], weights);
        // Normalized over the requested draws, repeats included.
        (int[] draws, double[] repeated) = Draws3.Select([0, 0, 1], 2, [scale, scale / 2]);
        Assert.Equal([0, 0, 1], draws);
        Assert.Equal(0.4, repeated[0], 15);
        Assert.Equal(0.4, repeated[1], 15);
        Assert.Equal(0.2, repeated[2], 15);
        Assert.Throws<ArgumentException>(() => Draws3.Select([1], 2, [scale, 0.0]));
    }

    [Theory]
    [InlineData(1e308)]
    [InlineData(5e307)]
    [InlineData(1e-320)]
    public void TheSummary_ReadsTwoAndFourAsMeanThreeSdOne(double scale)
    {
        PosteriorSummary s = PosteriorSummary.Of([2.0, 4.0], [scale, scale], [0.5]);
        Assert.Equal(3.0, s.Mean, 14);
        Assert.Equal(1.0, s.Sd, 14);
        Assert.Equal(3.0, s.Quantiles[0], 14);
        // A zero weight still drops its draw out of the count and the range.
        PosteriorSummary one = PosteriorSummary.Of([2.0, 4.0], [scale, 0.0], [0.5]);
        Assert.Equal(1, one.Draws);
        Assert.Equal(2.0, one.Min);
        Assert.Equal(2.0, one.Max);
        Assert.Equal(2.0, one.Mean, 14);
    }

    [Fact]
    public void ThePowerOfTwoMultiple_ReadsBitForBit()
    {
        double big = Math.ScaleB(1.0, 1000);
        Assert.Equal(Draws3.Select(null, 3, [1.0, 2.0, 5.0]).Weights, Draws3.Select(null, 3, [big, 2 * big, 5 * big]).Weights);
        PosteriorSummary a = PosteriorSummary.Of([2.0, 4.0, 7.0], [1.0, 2.0, 5.0]), b = PosteriorSummary.Of([2.0, 4.0, 7.0], [big, 2 * big, 5 * big]);
        Assert.Equal(a.Mean, b.Mean);
        Assert.Equal(a.Sd, b.Sd);
        Assert.Equal(a.Quantiles, b.Quantiles);
        Assert.Equal(Whitening.Normalized([1.0, 2.0, 5.0], 3, "w"), Whitening.Normalized([big, 2 * big, 5 * big], 3, "w"));
    }

    // The selection and the summary as they stood before 2026-11-13, verbatim, on ordinary weights: bit for bit.
    [Fact]
    [Trait("Tier", "Pinned")]
    public void OrdinaryWeights_KeepTheEarlierFormsBits()
    {
        var rng = new Xoshiro256PlusPlus(1113);
        for (int trial = 0; trial < 2000; trial++)
        {
            int n = 1 + (int)(rng.NextUInt64() % 50);
            var w = new double[n];
            var x = new double[n];
            for (int d = 0; d < n; d++)
            {
                w[d] = rng.NextDouble() < 0.1 ? 0.0 : Math.Exp(-12.0 * rng.NextDouble());   // SMC-like: below 1, spanning e⁻¹²
                x[d] = 10.0 * (rng.NextDouble() - 0.5);
            }
            if (w.All(v => v == 0.0)) w[0] = 0.5;
            int[] requested = [.. Enumerable.Range(0, n).Select(_ => (int)(rng.NextUInt64() % (ulong)n))];
            if (requested.All(d => w[d] == 0.0)) requested[0] = Array.FindIndex(w, v => v > 0.0);

            double[] selected = Draws3.Select(requested, n, w).Weights;
            double[] earlier = [.. requested.Select(d => w[d])];
            double total = 0.0;
            foreach (double v in earlier) total += v;
            for (int d = 0; d < earlier.Length; d++) earlier[d] /= total;
            Assert.Equal(earlier.Select(BitConverter.DoubleToInt64Bits), selected.Select(BitConverter.DoubleToInt64Bits));

            PosteriorSummary s = PosteriorSummary.Of(x, w);
            double sum = 0.0, mass = 0.0;
            for (int d = 0; d < n; d++) { mass += w[d]; sum += w[d] * x[d]; }
            double mean = sum / mass, squares = 0.0;
            for (int d = 0; d < n; d++) squares += w[d] * (x[d] - mean) * (x[d] - mean);
            Assert.Equal(BitConverter.DoubleToInt64Bits(mean), BitConverter.DoubleToInt64Bits(s.Mean));
            Assert.Equal(BitConverter.DoubleToInt64Bits(Math.Sqrt(squares / mass)), BitConverter.DoubleToInt64Bits(s.Sd));
        }
    }

    // ── The subspace posteriors ──────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void TheSubspacePosteriors_ReadOnePosteriorAtEveryScale()
    {
        double[][,] c = [new double[,] { { 2, 1 }, { 1, 3 } }, new double[,] { { 1, 0 }, { 0, 4 } }, new double[,] { { 5, 2 }, { 2, 1 } }];
        var pairs = new[] { new CoActivePair(c[0], c[1], c[2]), new CoActivePair(c[1], c[2], c[0]) };
        SubspacePosterior reference = SubspacePosterior.Of(c, [1.0, 2.0, 3.0], 1);
        CoActivePosterior coReference = CoActivePosterior.Of(pairs, [1.0, 3.0]);
        double[] normalized = Whitening.Normalized([1.0, 2.0, 3.0], 3, "w");
        foreach (double scale in Scales.Where(s => s * 3 <= double.MaxValue).Append(1e-320))
        {
            double[] w = Whitening.Normalized([scale, 2 * scale, 3 * scale], 3, "w");
            for (int d = 0; d < 3; d++) AssertSame(normalized[d], w[d], scale, $"normalized weight {d}");
            SubspacePosterior p = SubspacePosterior.Of(c, [scale, 2 * scale, 3 * scale], 1);
            CoActivePosterior q = CoActivePosterior.Of(pairs, [scale, 3 * scale]);
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                {
                    AssertSame(reference.MeanC[i, j], p.MeanC[i, j], scale, "mean C");
                    AssertSame(coReference.MeanCf[i, j], q.MeanCf[i, j], scale, "co-active mean C_f");
                }
        }
    }

    // ── The resamplers ───────────────────────────────────────────────────────────────────────────────────────────────

    // The calibration's outer draws read only the count and the weights.
    private sealed class WeightsOnly(double[] weights) : ICalibrationEmulator
    {
        public IReadOnlyList<CalibrationParameter> Parameters => throw new NotSupportedException();
        public int DrawCount => weights.Length;
        public IReadOnlyList<double>? DrawWeights => weights;
        public int Outputs => throw new NotSupportedException();
        public EmulatorNoise Noise(int draw) => throw new NotSupportedException();
        public double[,]? TruncationResiduals() => throw new NotSupportedException();
        public void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output) => throw new NotSupportedException();
    }

    [Fact]
    public void TheCalibrationsOuterDraws_PickTheSameIndicesAtEveryScale()
    {
        int[] Draws(double scale) => NestedCalibrationSampler.OuterDraws(new WeightsOnly([scale, scale / 3 * 2, scale, 0.0, scale / 3]), 400, new Xoshiro256PlusPlus(17));
        int[] reference = Draws(1.0);
        Assert.DoesNotContain(3, reference);
        foreach (double scale in Scales) Assert.Equal(reference, Draws(scale));
        Assert.Throws<ArgumentException>(() => Draws(0.0));
    }

    [Fact]
    public void TheJointDraws_ResampleTheSameAncestorsAtEveryScale()
    {
        BassFit fit = SensitivityFixture.FitScalar(1);
        int total = fit.Draws.Sum(c => c.Count);
        BassJointDraws Pair(double scale) => BassJointDraws.Pair([fit with { DrawWeights = [.. Enumerable.Range(0, total).Select(d => scale / 3 * (1 + d % 3))] }], 5);
        BassJointDraws reference = Pair(1.0);
        Assert.True(reference.Resampled);
        foreach (double scale in Scales)
        {
            BassJointDraws j = Pair(scale);
            for (int d = 0; d < reference.Count; d++) Assert.Equal(reference.Draw(0, d), j.Draw(0, d));
        }
    }

    // ── The fit's readers ────────────────────────────────────────────────────────────────────────────────────────────

    private static BassFit Weighted(BassFit fit, double scale)
    {
        int total = fit.Draws.Sum(c => c.Count);
        return fit with { DrawWeights = [.. Enumerable.Range(0, total).Select(d => scale / 3 * (1 + d % 3))] };
    }

    [Fact]
    public void TheFitsSummaryPredictionAndSensitivity_ReadOnePosteriorAtEveryScale()
    {
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = SensitivityFixture.FitScalar(1);
        BassFit ordinary = Weighted(fit, 1.0);
        BassSummary summary = ordinary.Summary();
        double[,] mean = ordinary.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false }).Mean();
        var decomposition = ordinary.Sensitivity(data).Decompositions[0];
        double[] main = decomposition.MeanMainIndices();
        SobolSummary sobol = SobolSummary.Of(decomposition);
        foreach (double scale in Scales)
        {
            BassFit weighted = Weighted(fit, scale);
            BassSummary s = weighted.Summary();
            AssertSame(summary.Posterior.Count.Mean, s.Posterior.Count.Mean, scale, "the term count's mean");
            AssertSame(summary.Posterior.Sigma2.Sd, s.Posterior.Sigma2.Sd, scale, "σ²'s sd");
            double[,] m = weighted.Predict(data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false }).Mean();
            for (int i = 0; i < m.GetLength(0); i++) AssertSame(mean[i, 0], m[i, 0], scale, $"the predictive mean at case {i}");
            var d = weighted.Sensitivity(data).Decompositions[0];
            double[] dm = d.MeanMainIndices();
            for (int k = 0; k < main.Length; k++) AssertSame(main[k], dm[k], scale, $"main index {k}'s mean");
            SobolSummary ss = SobolSummary.Of(d);
            for (int k = 0; k < sobol.MainIndex.Count; k++) AssertSame(sobol.MainIndex[k].Mean, ss.MainIndex[k].Mean, scale, $"main index {k}'s summary");
        }
        output.WriteLine($"Term count mean {summary.Posterior.Count.Mean:R}; {main.Length} main indices.");
    }

    [Fact]
    public void TheResidualCheck_ReadsOnePosteriorMeanAtEveryScale()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        var index = new ResidualIndex.AlongVariable(0);
        ResidualDependence reference = Weighted(fit, 1.0).ResidualDependenceFor(data, index);
        foreach (double scale in Scales)
        {
            ResidualDependence r = Weighted(fit, scale).ResidualDependenceFor(data, index);
            for (int i = 0; i < reference.PerCase.Length; i++) AssertSame(reference.PerCase[i], r.PerCase[i], scale, $"case {i}'s statistic");
        }
    }

    [Fact]
    public void TheProjectionReport_ReadsOneMeanTermCountAtEveryScale()
    {
        BassProjectedFitTests.RidgeFits f = BassProjectedFitTests.Fixture.Value;
        double reference = f.Refit.Compare(Weighted(f.Native, 1.0), f.Data, f.Analysis, new ProjectionComparisonSpec()).Native.BasisFunctions;
        foreach (double scale in Scales)
            AssertSame(reference, f.Refit.Compare(Weighted(f.Native, scale), f.Data, f.Analysis, new ProjectionComparisonSpec()).Native.BasisFunctions, scale,
                       "the mean term count");
    }
}
