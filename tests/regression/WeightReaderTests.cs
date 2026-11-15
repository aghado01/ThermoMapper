using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The shared weight reader, <see cref="PredictiveChecks"/>' <c>Rows</c>, across the weights' whole finite range: its normalized log
/// weights against an exact reference (the weights summed as integers), the ordinary range bit for bit against the reader's earlier
/// form, positive rescaling, zero weights, and the predictive, CDF and carried-error readings built on them. The acceptance is conformal
/// slice 10's parked follow-up (issues/regression/shared-code-defects-brief.md, slice 1).
/// </summary>
public sealed class WeightReaderTests(ITestOutputHelper output)
{
    private const double SmallestNormal = 2.2250738585072014e-308;   // 2^-1022
    private const double Ln2Hi = 6.93147180369123816490e-01;          // ln 2 split so that k·Ln2Hi is exact for |k| < 2^21
    private const double Ln2Lo = 1.90821492927058770002e-10;

    // A weight in integer units of 2^-1074, the smallest subnormal: exact for every finite non-negative double.
    private static BigInteger Units(double w)
    {
        long bits = BitConverter.DoubleToInt64Bits(w);
        int exponent = (int)((bits >> 52) & 0x7FF);
        long mantissa = bits & 0xFFFFFFFFFFFFFL;
        return exponent == 0 ? new BigInteger(mantissa) : new BigInteger(mantissa | (1L << 52)) << (exponent - 1);
    }

    // The exact reference: ln(W_d / ΣW) for each positive weight, the sum exact in integers, the ratio taken to at least 62 bits by
    // integer division and read through a double where it is normal, else as ln Q − k ln 2 with ln 2 split.
    internal static double[] Reference(double[] weights)
    {
        BigInteger total = BigInteger.Zero;
        foreach (double w in weights) total += Units(w);
        var log = new double[weights.Count(w => w > 0.0)];
        int k = 0;
        foreach (double w in weights)
        {
            if (!(w > 0.0)) continue;
            BigInteger own = Units(w);
            long shift = 64 - (own.GetBitLength() - total.GetBitLength());
            double q = (double)((own << (int)shift) / total);
            double ratio = Math.ScaleB(q, -(int)shift);
            log[k++] = ratio >= SmallestNormal ? Math.Log(ratio) : Math.Log(q) - shift * Ln2Hi - shift * Ln2Lo;
        }
        return log;
    }

    // The reader's form before 2026-11-13, verbatim: ln(w / Σw), the sum taken in draw order.
    private static double[] EarlierForm(double[] weights)
    {
        double total = 0.0;
        foreach (double w in weights) total += w;
        return [.. weights.Where(w => w > 0.0).Select(w => Math.Log(w / total))];
    }

    // The rounding the reader's normal-ratio path may carry: the sum's (n − 1) roundings, the division's and the log's, each within
    // 2^-53, bounded by 2^-52 apiece; the log-difference path's are within that on its results of magnitude above 708.
    private static double Tolerance(double reference, int n) => (n + 4) * 2.220446049250313e-16 * Math.Max(1.0, Math.Abs(reference));

    private static double LogSumExp(params double[] values) => PredictiveChecks.LogSumExp(values);

    private static void AssertRead(double[] weights, double[] logWeights, int[] rows, ref double worst)
    {
        double[] reference = Reference(weights);
        int[] positive = [.. Enumerable.Range(0, weights.Length).Where(d => weights[d] > 0.0)];
        Assert.Equal(positive, rows);
        Assert.Equal(reference.Length, logWeights.Length);
        for (int k = 0; k < reference.Length; k++)
        {
            Assert.True(double.IsFinite(logWeights[k]), $"Weight {rows[k]} of [{string.Join(", ", weights.Select(w => w.ToString("R")))}] read {logWeights[k]}.");
            double error = Math.Abs(logWeights[k] - reference[k]), tolerance = Tolerance(reference[k], weights.Length);
            Assert.True(error <= tolerance, $"Weight {rows[k]} read {logWeights[k]:R} against {reference[k]:R}, off by {error:E3} past {tolerance:E3}.");
            worst = Math.Max(worst, error / tolerance);
        }
    }

    public static TheoryData<double[]> Extremes => new()
    {
        new[] { 1e300, 1e-300 },
        new[] { 1e308, 1e308 },
        new[] { double.MaxValue, double.MaxValue, double.MaxValue },
        new[] { double.MaxValue, double.Epsilon },
        new[] { double.Epsilon, double.Epsilon },
        new[] { 1.0, double.Epsilon },
        new[] { 3.0, double.Epsilon },
        new[] { SmallestNormal, double.Epsilon, 1e-320 },
        new[] { 0.0, 1e300, 0.0, 1e-300, 0.0 },
    };

    [Theory]
    [MemberData(nameof(Extremes))]
    public void ExtremeWeightsReadTheirExactNormalizedLogs(double[] weights)
    {
        (int[] rows, double[]? log) = PredictiveChecks.Rows(weights, weights.Length);
        double worst = 0.0;
        AssertRead(weights, log!, rows, ref worst);
    }

    [Fact]
    public void TheParkedInputsKeepTheirRowsAndValues()
    {
        (int[] rows, double[]? log) = PredictiveChecks.Rows([1e300, 1e-300], 2);
        Assert.Equal([0, 1], rows);
        Assert.Equal(0.0, log![0]);
        Assert.Equal(-1381.551055796427, log[1], 9);   // −600 ln 10

        (rows, log) = PredictiveChecks.Rows([1e308, 1e308], 2);
        Assert.Equal([0, 1], rows);
        Assert.Equal(Math.Log(0.5), log![0]);
        Assert.Equal(Math.Log(0.5), log[1]);
        Assert.Equal(0.5, Math.Exp(log[0]));
    }

    [Fact]
    [Trait("Tier", "Pinned")]
    public void OrdinaryWeightsKeepTheEarlierFormsBits()
    {
        var rng = new Xoshiro256PlusPlus(20261113);
        for (int trial = 0; trial < 4000; trial++)
        {
            int n = 1 + (int)(rng.NextUInt64() % 64);
            int band = 1 + (int)(rng.NextUInt64() % 120);
            var weights = new double[n];
            for (int d = 0; d < n; d++)
                weights[d] = rng.NextDouble() < 0.1 ? 0.0 : Math.ScaleB(1.0 + rng.NextDouble(), (int)(rng.NextUInt64() % (ulong)band) - band / 2);
            if (weights.All(w => w == 0.0)) weights[0] = 1.0;
            (_, double[]? log) = PredictiveChecks.Rows(weights, n);
            double[] earlier = EarlierForm(weights);
            for (int k = 0; k < earlier.Length; k++)
                Assert.Equal(BitConverter.DoubleToInt64Bits(earlier[k]), BitConverter.DoubleToInt64Bits(log![k]));
        }
    }

    // A weight anywhere in the finite range: zeros, the smallest subnormals, values near double.MaxValue, and otherwise a uniform
    // binary exponent within [low, high) with a uniform mantissa.
    private static double DrawWeight(Xoshiro256PlusPlus rng, int low, int high)
    {
        double u = rng.NextDouble();
        if (u < 0.08) return 0.0;
        if (u < 0.12) return double.Epsilon * (1 + (int)(rng.NextUInt64() % 4));
        if (u < 0.16) return double.MaxValue * (1.0 - 1e-3 * rng.NextDouble());
        return Math.ScaleB(1.0 + rng.NextDouble(), low + (int)(rng.NextUInt64() % (ulong)(high - low)));
    }

    private static double[] DrawWeights(Xoshiro256PlusPlus rng)
    {
        int n = 2 + (int)(rng.NextUInt64() % 40);
        int low = -1074 + (int)(rng.NextUInt64() % 2098);
        int high = low + 1 + (int)(rng.NextUInt64() % (ulong)(1024 - low));
        var weights = new double[n];
        for (int d = 0; d < n; d++) weights[d] = DrawWeight(rng, low, high);
        if (weights.All(w => w == 0.0)) weights[0] = 1.0;
        return weights;
    }

    [Fact]
    public void ASweepOverTheWholeRangeMatchesTheExactReference()
    {
        var rng = new Xoshiro256PlusPlus(611);
        double worst = 0.0;
        int fallback = 0;
        for (int trial = 0; trial < 4000; trial++)
        {
            double[] weights = DrawWeights(rng);
            (int[] rows, double[]? log) = PredictiveChecks.Rows(weights, weights.Length);
            AssertRead(weights, log!, rows, ref worst);
            if (EarlierForm(weights).Any(l => !(l >= Math.Log(SmallestNormal)))) fallback++;
        }
        output.WriteLine($"Largest error over its tolerance: {worst:F3}; {fallback} of 4000 vectors had a ratio the earlier form lost.");
        Assert.True(fallback > 1000, "The sweep must reach the range the earlier form misread.");
    }

    [Fact]
    public void APositiveRescalingLeavesTheLogWeights()
    {
        var rng = new Xoshiro256PlusPlus(7);
        int exact = 0;
        for (int trial = 0; trial < 4000; trial++)
        {
            double[] weights = DrawWeights(rng);
            double[] positive = [.. weights.Where(w => w > 0.0)];
            int top = Math.ILogB(positive.Max()), bottom = Math.ILogB(positive.Min());
            // A power of two that keeps every weight normal when their range allows, so the rescaled weights are exact.
            int lowest = -1022 - bottom, highest = 1023 - top;
            int k = lowest >= highest ? 0 : lowest + (int)(rng.NextUInt64() % (ulong)(highest - lowest));
            double[] scaled = [.. weights.Select(w => Math.ScaleB(w, k))];
            if (scaled.Zip(weights).Any(p => Math.ScaleB(p.First, -k) != p.Second)) continue;   // a subnormal lost bits
            (_, double[]? a) = PredictiveChecks.Rows(weights, weights.Length);
            (_, double[]? b) = PredictiveChecks.Rows(scaled, scaled.Length);
            bool normal = positive.All(w => w >= SmallestNormal) && scaled.Where(w => w > 0.0).All(w => w >= SmallestNormal)
                       && a!.All(l => l >= Math.Log(SmallestNormal) + 1.0);
            for (int i = 0; i < a!.Length; i++)
            {
                if (normal) Assert.Equal(BitConverter.DoubleToInt64Bits(a[i]), BitConverter.DoubleToInt64Bits(b![i]));
                else Assert.True(Math.Abs(a[i] - b![i]) <= 2 * Tolerance(a[i], weights.Length), $"Rescaled by 2^{k}: {a[i]:R} against {b[i]:R}.");
            }
            if (normal) exact++;

            // A scale that is not a power of two moves each weight by a rounding: the log weights agree within it.
            double c = 0.75 + 0.5 * rng.NextDouble();
            double[] rounded = [.. positive.Select(w => w * c)];
            if (rounded.Any(w => !(w >= SmallestNormal) || double.IsInfinity(w)) || positive.Any(w => w < SmallestNormal)) continue;
            (_, double[]? r) = PredictiveChecks.Rows(rounded, rounded.Length);
            (_, double[]? p) = PredictiveChecks.Rows(positive, positive.Length);
            for (int i = 0; i < p!.Length; i++)
                Assert.True(Math.Abs(p[i] - r![i]) <= 3 * Tolerance(p[i], positive.Length), $"Rescaled by {c:R}: {p[i]:R} against {r[i]:R}.");
        }
        output.WriteLine($"{exact} of 4000 rescalings compared bit for bit.");
        Assert.True(exact > 1000);
    }

    [Fact]
    public void ZeroWeightsDropOutAndLeaveTheOthersBits()
    {
        var rng = new Xoshiro256PlusPlus(19);
        for (int trial = 0; trial < 2000; trial++)
        {
            double[] weights = DrawWeights(rng);
            double[] positive = [.. weights.Where(w => w > 0.0)];
            (int[] rows, double[]? log) = PredictiveChecks.Rows(weights, weights.Length);
            (_, double[]? dense) = PredictiveChecks.Rows(positive, positive.Length);
            Assert.Equal(positive.Length, rows.Length);
            Assert.All(rows, d => Assert.True(weights[d] > 0.0));
            Assert.Equal(dense!.Select(BitConverter.DoubleToInt64Bits), log!.Select(BitConverter.DoubleToInt64Bits));
        }
    }

    [Fact]
    public void TheDensitiesAdmitEveryFiniteWeightAndRefuseTheRest()
    {
        var values = new double[2, 1];
        PointwiseLogDensity Density(double[] w) => new(values, UnitPartition.EachObservation(1), ConformalBayesTests.Theta) { Weights = w };
        foreach (double[] w in new[] { [1e308, 1e308], [double.MaxValue, double.MaxValue], [double.Epsilon, 0.0], new[] { 0.0, 1e-300 } })
            Assert.Same(w, Density(w).Weights);
        foreach (double[] w in new[] { [0.0, 0.0], [double.NaN, 1.0], [double.PositiveInfinity, 1.0], new[] { -1.0, 2.0 } })
            Assert.Throws<ArgumentException>(() => Density(w));
    }

    // A tilt that makes the first parked input's tiny row material: the dominant draw's terms sit e^-gap below the tiny one's.
    [Fact]
    public void TheTinyRowEntersTheWeightedPredictiveUnderATilt()
    {
        double[] weights = [1e300, 1e-300];
        double[] reference = Reference(weights);
        double gap = reference[0] - reference[1];
        const int columns = 7;
        var values = new double[2, columns];
        for (int c = 0; c < columns; c++) values[0, c] = -gap + (c - 3);
        (int[] rows, double[]? log) = PredictiveChecks.Rows(weights, 2);
        var predictive = new double[columns];
        SplitConformal.LogPredictive(values, rows, log, 0, columns, predictive);
        for (int c = 0; c < columns; c++)
        {
            double expected = LogSumExp(reference[0] + values[0, c], reference[1] + values[1, c]);
            Assert.Equal(expected, predictive[c], 9);
        }
        Assert.Equal(reference[1] + Math.Log(2.0), predictive[3], 9);   // both terms equal: twice the tiny row's

        // WAIC over three draws, two dominant: its lppd reads the tilted tiny row, its variance the dominant pair.
        double[] three = [1e300, 1e300, 1e-300];
        double[] w3 = Reference(three);
        double lift = w3[0] - w3[2];
        var ll = new double[,] { { -lift - 0.3, -lift + 0.4 }, { -lift - 0.8, -lift - 0.1 }, { 0.0, 0.0 } };
        PredictiveCheck waic = PredictiveChecks.Waic(new PointwiseLogDensity(ll, UnitPartition.EachObservation(2), ConformalBayesTests.Theta) { Weights = three });
        double[] linear = [.. w3.Select(Math.Exp)];
        double reliability = 1.0 - linear.Sum(v => v * v);
        for (int u = 0; u < 2; u++)
        {
            double lppd = LogSumExp(w3[0] + ll[0, u], w3[1] + ll[1, u], w3[2] + ll[2, u]);
            double mean = Enumerable.Range(0, 3).Sum(d => linear[d] * ll[d, u]);
            double variance = Enumerable.Range(0, 3).Sum(d => linear[d] * (ll[d, u] - mean) * (ll[d, u] - mean)) / reliability;
            Assert.Equal(lppd - variance, waic.PerUnit[u], 9);
            Assert.True(lppd - LogSumExp(w3[0] + ll[0, u], w3[1] + ll[1, u]) > 0.3);   // the tiny row is material
        }
    }

    // Slice 10's carried CDF error, with the first parked input in place of [1, ε]: unit 1's bound under the turn tilted by unit 0.
    [Fact]
    public void TheTinyRowEntersTheCarriedCdfError()
    {
        double[] weights = [1e300, 1e-300];
        double[] reference = Reference(weights);
        double tilt = reference[0] - reference[1] - 0.5;
        foreach (int k in new[] { 0, -20, 20 })
        {
            double[] scaled = [.. weights.Select(w => Math.ScaleB(w, k))];
            double[,] ll = { { -tilt, 0, 0 }, { 0, 0, 0 } };
            double[,] fc = { { Math.Log(.25), Math.Log(.25), Math.Log(.25) }, { Math.Log(.75), Math.Log(.75), Math.Log(.75) } };
            var training = new PointwiseLogDensity(ll, UnitPartition.EachObservation(3), ConformalBayesTests.Theta) { LogCdf = fc, Weights = scaled };
            var candidates = new CandidateLogDensity(new double[2, 1], [[0.0]], ConformalBayesTests.Theta)
                { LogCdf = new double[,] { { Math.Log(.5) }, { Math.Log(.5) } }, Weights = scaled };
            var d = ConformalBayes.Distribution(training, candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried }).Distributions[0];
            double a = Math.Exp(reference[0] - tilt - reference[1]), w0 = a / (1 + a), w1 = 1 / (1 + a);
            double mean = w0 * .25 + w1 * .75;
            double expected = d.Margin!.Multiple * Math.Sqrt(w0 * w0 * (.25 - mean) * (.25 - mean) + w1 * w1 * (.75 - mean) * (.75 - mean));
            Assert.True(double.IsFinite(d.Margin.Multiple));
            Assert.Equal(expected, d.Margin.UnitCarriedErrors[1], 12);
            Assert.InRange(w1, .37, .38);
        }
    }

    // The second parked input reads as two weights of ½, so every reading equals the one under [1, 1], bit for bit.
    [Fact]
    public void AnOverflowingSumReadsAsItsRatios()
    {
        double[] big = [1e308, 1e308], unit = [1.0, 1.0];
        var ll = new double[,] { { -1.0, -2.5 }, { -1.7, -0.4 } };
        double[,] fc = { { Math.Log(.2), Math.Log(.6) }, { Math.Log(.7), Math.Log(.1) } };
        PointwiseLogDensity Density(double[] w) => new(ll, UnitPartition.EachObservation(2), ConformalBayesTests.Theta) { LogCdf = fc, Weights = w };
        Assert.Equal(PredictiveChecks.Waic(Density(unit)).PerUnit, PredictiveChecks.Waic(Density(big)).PerUnit);

        (int[] rows, double[]? log) = PredictiveChecks.Rows(big, 2);
        double[] linear = CdfScores.Weights(rows, log);
        Assert.Equal([0.5, 0.5], linear);
        var cdf = new CdfScore[2];
        CdfScores.Fixed(fc, rows, linear, 0, 2, cdf);
        Assert.Equal(0.5 * .2 + 0.5 * .7, cdf[0].Lower, 15);
        Assert.Equal(0.5 * .6 + 0.5 * .1, cdf[1].Lower, 15);

        double[] Carried(double[] w)
        {
            var candidates = new CandidateLogDensity(new double[2, 1], [[0.0]], ConformalBayesTests.Theta)
                { LogCdf = new double[,] { { Math.Log(.5) }, { Math.Log(.5) } }, Weights = w };
            return ConformalBayes.Distribution(Density(w), candidates, new ConformalOptions { MarginBounds = MarginBounds.Carried }).Distributions[0].Margin!.UnitCarriedErrors;
        }
        Assert.Equal(Carried(unit), Carried(big));
    }

    // ---- the guard (the defects brief's slice 10) ----

    // A list whose enumerator yields `shift` rows more or fewer than its Count, from its `from`-th enumeration on; its indexer is
    // true. Rows' fault cannot be called up, so the guard is tested by a list that miscounts on purpose.
    private sealed class MiscountingList(double[] values, int shift, int from) : IReadOnlyList<double>
    {
        private int _enumerations;
        public int Count => values.Length;
        public double this[int index] => values[index];
        public IEnumerator<double> GetEnumerator()
        {
            int n = _enumerations++ >= from ? Math.Max(0, values.Length + shift) : values.Length;
            for (int i = 0; i < n; i++) yield return values[i % values.Length];
        }
        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
    }

    // [1, ε] takes the scaled branch, so Rows passes over it twice; the first enumeration is its first pass, the second its scaled pass.
    [Theory]
    [InlineData(-1, 0, "first", 1)]
    [InlineData(1, 0, "first", 3)]
    [InlineData(-2, 1, "scaled", 0)]   // the run's fault: the scaled pass read nothing
    [InlineData(-1, 1, "scaled", 1)]
    [InlineData(1, 1, "scaled", 3)]
    public void APassThatReadsAnotherCountThrows_NamingItself(int shift, int from, string pass, int read)
    {
        var weights = new MiscountingList([1.0, double.Epsilon], shift, from);
        var thrown = Assert.Throws<InvalidOperationException>(() => PredictiveChecks.Rows(weights, 2));
        output.WriteLine(thrown.Message);
        Assert.Contains("PredictiveChecks.Rows", thrown.Message);
        Assert.Contains($"its {pass} pass over the draw weights read {read} rows where the list's Count is 2", thrown.Message);
        Assert.Contains(typeof(MiscountingList).FullName!, thrown.Message);
        Assert.Contains("Keep the run's folder", thrown.Message);
        Assert.Contains("tell the second coordinating session", thrown.Message);
    }

    [Fact]
    public void ARightCountReadsAsBefore_AndAnOrdinaryRowHasNoScaledPass()
    {
        // A list that miscounts only from its second enumeration: ordinary weights never take a second pass, so it reads as the array.
        double[] ordinary = [1.0, 2.0, 3.0];
        (int[] rows, double[]? log) = PredictiveChecks.Rows(new MiscountingList(ordinary, -3, 1), 3);
        (int[] arrayRows, double[]? arrayLog) = PredictiveChecks.Rows(ordinary, 3);
        Assert.Equal(arrayRows, rows);
        Assert.Equal(arrayLog, log);
        // A list that counts right through both passes reads the scaled branch as the array does, bit for bit.
        double[] tiny = [1.0, double.Epsilon];
        (rows, log) = PredictiveChecks.Rows(new MiscountingList(tiny, 0, 0), 2);
        (arrayRows, arrayLog) = PredictiveChecks.Rows(tiny, 2);
        Assert.Equal(arrayRows, rows);
        Assert.Equal(arrayLog, log);
    }
}
