using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Sensitivity;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The sampled-sensitivity brief's slice 0 gates on the pick-freeze estimator (<see cref="PickFreeze"/>), against the closed forms of
/// <see cref="SobolReferences"/>: on chosen cases, replicated over seeds, each piece's mean within its spread of the truth and its stated
/// error within 15% of its spread; the error falling as N^(−1/2); a zero-variance function; shift invariance, with the uncentred form's
/// cost beside it; and a seeded sweep over far ranges of scale, offset, boxes and small indices.
/// </summary>
public sealed class PickFreezeTests(ITestOutputHelper output)
{
    private static SobolReference Named(string name) => SobolReferences.All().Single(r => r.Name == name);

    // A and B, N rows each from the uniform law on the reference's box, one stream.
    internal static (double[,] A, double[,] B) Bases(SobolReference f, int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[,] Draw()
        {
            var x = new double[n, f.Dimension];
            for (int r = 0; r < n; r++)
                for (int v = 0; v < f.Dimension; v++) x[r, v] = f.Law.Lower[v] + (f.Law.Upper[v] - f.Law.Lower[v]) * rng.NextDouble();
            return x;
        }
        return (Draw(), Draw());
    }

    // Every pair of inputs, so each two-way main variance is read.
    private static PickFreezeDesign Design(SobolReference f, MainVarianceForm form = MainVarianceForm.Saltelli) =>
        new(f.Dimension, [.. Enumerable.Range(0, f.Dimension).SelectMany(i => Enumerable.Range(i + 1, f.Dimension - i - 1).Select(j => new[] { i, j }))], form);

    public static TheoryData<string, MainVarianceForm> NamesAndForms()
    {
        var cases = new TheoryData<string, MainVarianceForm>();
        foreach (SobolReference r in SobolReferences.All())
            foreach (MainVarianceForm form in new[] { MainVarianceForm.Saltelli, MainVarianceForm.TwoDifference }) cases.Add(r.Name, form);
        return cases;
    }

    internal sealed record Piece(string Name, double Truth, double Estimate, double Error);

    internal static List<Piece> Pieces(SobolReference f, PickFreezeDesign design, PickFreezeEstimate e)
    {
        double v = f.Variance;
        var pieces = new List<Piece> { new("mean", f.Mean, e.Mean, e.MeanError), new("V", v, e.Variance, e.VarianceError) };
        for (int k = 0; k < design.Blocks.Count; k++)
        {
            int[] u = design.Blocks[k];
            string label = string.Join(",", u.Select(i => i + 1));
            double main = f.MainVariance(u);
            pieces.Add(new($"V_{label}", main, e.Main[k], e.MainError[k]));
            pieces.Add(new($"S_{label}", main / v, e.MainIndex[k], e.MainIndexError[k]));
        }
        for (int i = 0; i < f.Dimension; i++)
        {
            double total = f.TotalVariance(i);
            pieces.Add(new($"T_{i + 1}", total, e.Total[i], e.TotalError[i]));
            pieces.Add(new($"ST_{i + 1}", total / v, e.TotalIndex[i], e.TotalIndexError[i]));
        }
        return pieces;
    }

    [Theory]
    [MemberData(nameof(NamesAndForms))]
    public void Pieces_AreUnbiasedWithinTheirStatedError(string name, MainVarianceForm form)
    {
        const int Replicates = 300, N = 2048;
        SobolReference f = Named(name);
        PickFreezeDesign design = Design(f, form);
        int[] seeds = SeedTree.Derive(20261115, Replicates);
        var runs = new List<List<Piece>>();
        foreach (int seed in seeds)
        {
            (double[,] a, double[,] b) = Bases(f, N, seed);
            PickFreezeEstimate e = PickFreeze.Estimate(design, f.Value, a, b);
            Assert.Equal(PickFreezeStatus.Defined, e.Status & PickFreezeStatus.ZeroVariance);
            runs.Add(Pieces(f, design, e));
        }
        Summary s = Summarize(runs, f.Variance);
        output.WriteLine($"{name}, {form}: {s.Tested} pieces tested, {s.Exact} exact zeros; worst |bias| / its error {s.WorstBias:F2} ({s.WorstBiasPiece}); "
                         + $"stated error / spread in [{s.LowestCalibration:F3}, {s.HighestCalibration:F3}] ({s.LowestPiece}, {s.HighestPiece})");
        foreach (string line in s.Lines) output.WriteLine(line);
        Assert.True(s.WorstBias <= 4.0, $"{name}: {s.WorstBiasPiece} is biased by {s.WorstBias:F2} of its error");
        Assert.True(s.LowestCalibration >= 0.85 && s.HighestCalibration <= 1.15,
                    $"{name}: stated error / spread runs from {s.LowestCalibration:F3} ({s.LowestPiece}) to {s.HighestCalibration:F3} ({s.HighestPiece})");
    }

    internal sealed record Summary(int Tested, int Exact, double WorstBias, string WorstBiasPiece, double LowestCalibration, string LowestPiece,
                                   double HighestCalibration, string HighestPiece, List<string> Lines);

    // Per piece over the replicates: the bias in units of the mean's error, sd(estimate)/√R, and the mean stated error over the spread.
    // A piece whose truth is 0 within the reference's rounding and whose every estimate is within rounding of 0 (its per-row terms vanish)
    // is an exact zero.
    internal static Summary Summarize(List<List<Piece>> runs, double variance)
    {
        int count = runs[0].Count, r = runs.Count, tested = 0, exact = 0;
        double worstBias = 0.0, low = double.PositiveInfinity, high = 0.0;
        string worstPiece = "", lowPiece = "", highPiece = "";
        var lines = new List<string>();
        for (int k = 0; k < count; k++)
        {
            string piece = runs[0][k].Name;
            double truth = runs[0][k].Truth;
            double[] est = [.. runs.Select(run => run[k].Estimate)];
            double[] err = [.. runs.Select(run => run[k].Error)];
            double scale = piece.StartsWith('S') ? 1.0 : variance;
            if (Math.Abs(truth) <= 1e-12 * scale && est.All(x => Math.Abs(x) <= 1e-10 * scale))
            {
                exact++;
                continue;
            }
            tested++;
            double mean = est.Average(), sd = Math.Sqrt(est.Sum(x => (x - mean) * (x - mean)) / (r - 1));
            double bias = Math.Abs(mean - truth) / (sd / Math.Sqrt(r)), calibration = err.Average() / sd;
            lines.Add($"  {piece,-8} truth {truth,12:G6} mean {mean,12:G6} sd {sd,10:E2} bias/err {bias,5:F2} stated/sd {calibration:F3}");
            if (bias > worstBias) (worstBias, worstPiece) = (bias, piece);
            if (calibration < low) (low, lowPiece) = (calibration, piece);
            if (calibration > high) (high, highPiece) = (calibration, piece);
        }
        return new Summary(tested, exact, worstBias, worstPiece, low, lowPiece, high, highPiece, lines);
    }

    [Fact]
    public void Error_FallsAsTheInverseRootOfN()
    {
        // The spread of each piece over 100 replicates at N = 256 … 16384, and the mean stated error: both slopes on log N near −½.
        const int Replicates = 100;
        int[] sizes = [256, 1024, 4096, 16384];
        foreach (SobolReference f in new SobolReference[] { SobolReferences.Ishigami(), SobolReferences.G() })
        {
            PickFreezeDesign design = new(f.Dimension);
            var spread = new List<double[]>();
            var stated = new List<double[]>();
            string[] names = [];
            foreach (int n in sizes)
            {
                var runs = new List<List<Piece>>();
                foreach (int seed in SeedTree.Derive(n, Replicates))
                {
                    (double[,] a, double[,] b) = Bases(f, n, seed);
                    runs.Add(Pieces(f, design, PickFreeze.Estimate(design, f.Value, a, b)));
                }
                names = [.. runs[0].Select(p => p.Name)];
                spread.Add([.. Enumerable.Range(0, names.Length).Select(k => Deviation([.. runs.Select(run => run[k].Estimate)]))]);
                stated.Add([.. Enumerable.Range(0, names.Length).Select(k => runs.Average(run => run[k].Error))]);
            }
            double[] logN = [.. sizes.Select(n => Math.Log(n))];
            for (int k = 0; k < names.Length; k++)
            {
                if (spread[0][k] == 0.0) continue;   // an exact zero: Ishigami's x₃ has no main effect and no pair term moves it
                double slope = Slope(logN, [.. spread.Select(s => Math.Log(s[k]))]), statedSlope = Slope(logN, [.. stated.Select(s => Math.Log(s[k]))]);
                output.WriteLine($"{f.Name} {names[k],-6}: spread slope {slope:F3}, stated slope {statedSlope:F3}");
                Assert.True(Math.Abs(slope + 0.5) <= 0.08, $"{f.Name} {names[k]}: the spread falls as N^{slope:F3}");
                Assert.True(Math.Abs(statedSlope + 0.5) <= 0.05, $"{f.Name} {names[k]}: the stated error falls as N^{statedSlope:F3}");
            }
        }
    }

    [Fact]
    public void ZeroVariance_LeavesTheIndicesUndefined()
    {
        foreach (double c in new[] { 0.0, 1.0, -3.7e5, 1e200 })
        {
            var f = new ConstantReference(c, 3);
            PickFreezeDesign design = Design(f);
            (double[,] a, double[,] b) = Bases(f, 64, 7);
            PickFreezeEstimate e = PickFreeze.Estimate(design, f.Value, a, b);
            Assert.Equal(PickFreezeStatus.ZeroVariance, e.Status);
            Assert.Equal(c, e.Mean);
            Assert.Equal(0.0, e.Variance);
            Assert.All(e.Main.Concat(e.MainError).Concat(e.Total).Concat(e.TotalError), x => Assert.Equal(0.0, x));
            Assert.All(e.MainIndex.Concat(e.MainIndexError).Concat(e.TotalIndex).Concat(e.TotalIndexError), x => Assert.True(double.IsNaN(x)));
        }
        // A variation below the tolerance times the mean's size is zero within rounding; one above it is read, however small.
        PickFreezeDesign one = new(2);
        (double[,] a2, double[,] b2) = Bases(new ConstantReference(0.0, 2), 64, 8);
        PickFreezeEstimate near = PickFreeze.Estimate(one, x => 1e6 * (1.0 + 1e-12 * x[0]), a2, b2);
        Assert.Equal(PickFreezeStatus.ZeroVariance, near.Status);
        PickFreezeEstimate small = PickFreeze.Estimate(one, x => 1e6 * (1.0 + 1e-8 * x[0]), a2, b2);
        Assert.Equal(PickFreezeStatus.Defined, small.Status);
        Assert.True(Math.Abs(small.TotalIndex[0] - 1.0) <= 4.0 * small.TotalIndexError[0], $"{small.TotalIndex[0]:R} ± {small.TotalIndexError[0]:R}");
        Assert.Equal(0.0, small.Total[1]);
    }

    [Fact]
    public void Estimates_AreShiftInvariant_AndTheUncentredFormIsNot()
    {
        // f + c reads as f to within the values' own rounding, for c up to 10⁸ standard deviations; Saltelli's uncentred form,
        // mean(f(B)(f(A_B^i) − f(A))), has a spread that grows with c, about 2c²·T_i/N in variance.
        SobolReference f = SobolReferences.Ishigami();
        PickFreezeDesign design = new(f.Dimension);
        double sd = Math.Sqrt(f.Variance);
        const int N = 2048;
        foreach (double shift in new[] { 10.0, 1e4, 1e8 })
        {
            double worst = 0.0;
            foreach (int seed in SeedTree.Derive(3, 8))
            {
                (double[,] a, double[,] b) = Bases(f, N, seed);
                PickFreezeEstimate e0 = PickFreeze.Estimate(design, f.Value, a, b);
                PickFreezeEstimate e1 = PickFreeze.Estimate(design, x => f.Value(x) + shift * sd, a, b);
                for (int k = 0; k < design.Blocks.Count; k++)
                {
                    worst = Math.Max(worst, Math.Abs(e1.Main[k] - e0.Main[k]) / e0.MainError[k]);
                    worst = Math.Max(worst, Math.Abs(e1.MainError[k] - e0.MainError[k]) / e0.MainError[k]);
                }
                for (int i = 0; i < f.Dimension; i++)
                    if (e0.TotalError[i] > 0.0) worst = Math.Max(worst, Math.Abs(e1.Total[i] - e0.Total[i]) / e0.TotalError[i]);
                worst = Math.Max(worst, Math.Abs(e1.Variance - e0.Variance) / e0.VarianceError);
            }
            output.WriteLine($"shift {shift:E0} sd: worst change / error {worst:E2}");
            Assert.True(worst <= 1e-4, $"shift {shift:E0}: {worst:E2}");
        }

        // The evidence for centring: the spread of V₁ over replicates, centred and uncentred, as the mean moves. The comparison is
        // registered in ComparativeSweep with its seed.
        foreach (double shift in new[] { 0.0, 10.0, 100.0 })
        {
            (double centred, double uncentred) = CentringSpread(shift, 11);
            output.WriteLine($"mean {shift} sd: spread of V₁ centred {centred:E3}, uncentred {uncentred:E3}, ratio {uncentred / centred:F2}");
            if (shift >= 100.0) Assert.True(20.0 * centred < uncentred, $"the uncentred form's spread is only {uncentred / centred:F2} times the centred one's at {shift} sd");
        }
    }

    // Ishigami's V₁ over 200 replicates at N = 2048, the function shifted by `shift` standard deviations: the spread of the centred
    // estimate and of Saltelli's uncentred form, mean(f(B)(f(A_B^1) − f(A))), the replicates' seeds derived from `seed`.
    internal static (double Centred, double Uncentred) CentringSpread(double shift, int seed)
    {
        const int N = 2048, Replicates = 200;
        SobolReference f = SobolReferences.Ishigami();
        PickFreezeDesign design = new(f.Dimension);
        double sd = Math.Sqrt(f.Variance);
        var centred = new List<double>();
        var uncentred = new List<double>();
        foreach (int replicate in SeedTree.Derive(seed, Replicates))
        {
            (double[,] a, double[,] b) = Bases(f, N, replicate);
            double F(ReadOnlySpan<double> x) => f.Value(x) + shift * sd;
            centred.Add(PickFreeze.Estimate(design, F, a, b).Main[0]);
            double[,] h = design.Hybrid(a, b, 0);
            double s = 0.0;
            for (int j = 0; j < N; j++) s += F(Row(b, j)) * (F(Row(h, j)) - F(Row(a, j)));
            uncentred.Add(s / N);
        }
        return (Deviation([.. centred]), Deviation([.. uncentred]));
    }

    [Theory]
    [InlineData(MainVarianceForm.Saltelli)]
    [InlineData(MainVarianceForm.TwoDifference)]
    public void SeededSweep_FarRanges(MainVarianceForm form)
    {
        // 300 cases: g-functions with indices down to 10⁻⁸ of the variance, Ishigami's family, and exponential ridge sums on boxes from
        // 10⁻³ to 10³ wide and centred out to ±10³; each scaled by 10^±100 and shifted by up to 10⁹ standard deviations. One estimate
        // per case at N = 4096, with up to three pairs; each piece's z = (estimate − truth)/stated error.
        var rng = new Xoshiro256PlusPlus(20261115);
        var z = new List<(double Z, string Where)>();
        int exact = 0, cases = 300;
        for (int c = 0; c < cases; c++)
        {
            SobolReference f = SweepCase(rng);
            int p = f.Dimension;
            var pairs = new List<int[]>();
            for (int t = 0; t < 3 && p >= 2; t++)
            {
                int i = rng.NextInt(p), j = rng.NextInt(p - 1);
                if (j >= i) j++;
                pairs.Add([Math.Min(i, j), Math.Max(i, j)]);
            }
            var design = new PickFreezeDesign(p, pairs, form);
            (double[,] a, double[,] b) = Bases(f, 4096, (int)rng.NextUInt64());
            PickFreezeEstimate e = PickFreeze.Estimate(design, f.Value, a, b);
            Assert.True(e.Status == PickFreezeStatus.Defined || e.Status == PickFreezeStatus.BelowZeroBeyondError, $"case {c} ({f.Name}): {e.Status}");
            foreach (Piece piece in Pieces(f, design, e))
            {
                double scale = piece.Name.StartsWith('S') ? 1.0 : piece.Name == "mean" ? Math.Abs(f.Mean) + Math.Sqrt(f.Variance) : f.Variance;
                if (Math.Abs(piece.Truth) <= 1e-12 * scale && Math.Abs(piece.Estimate) <= 1e-10 * scale)
                {
                    exact++;
                    continue;
                }
                Assert.True(piece.Error > 0.0, $"case {c} ({f.Name}) {piece.Name}: estimate {piece.Estimate:R}, truth {piece.Truth:R}, no error");
                z.Add(((piece.Estimate - piece.Truth) / piece.Error, $"case {c} ({f.Name}) {piece.Name}"));
            }
        }
        double mean = z.Average(t => t.Z), sd = Math.Sqrt(z.Sum(t => (t.Z - mean) * (t.Z - mean)) / (z.Count - 1));
        double beyond3 = z.Count(t => Math.Abs(t.Z) > 3.0) / (double)z.Count;
        var worst = z.MaxBy(t => Math.Abs(t.Z));
        output.WriteLine($"{form}: {z.Count} z-scores over {cases} cases ({exact} exact zeros): mean {mean:F3}, sd {sd:F3}, beyond 3: {beyond3:P2}, worst {worst.Z:F2} at {worst.Where}");
        Assert.True(Math.Abs(mean) <= 0.1, $"mean z {mean:F3}");
        Assert.True(sd >= 0.9 && sd <= 1.12, $"sd of z {sd:F3}");
        Assert.True(beyond3 <= 0.01, $"{beyond3:P2} beyond 3");
        Assert.True(Math.Abs(worst.Z) <= 6.0, $"worst z {worst.Z:F2} at {worst.Where}");
    }

    [Fact]
    public void Design_BlocksAndHybrids()
    {
        var design = new PickFreezeDesign(4, [[2, 0, 3], [1, 3]]);
        string[] blocks = [.. design.Blocks.Select(u => string.Join(",", u))];
        Assert.Equal(["0", "1", "2", "3", "0,2", "0,3", "1,3", "2,3", "0,2,3"], blocks);
        Assert.Equal(4096L * 11, design.Values(4096));
        Assert.Equal(8, design.IndexOf([0, 2, 3]));
        Assert.Equal(-1, design.IndexOf([1, 2]));
        double[,] a = { { 1, 2, 3, 4 }, { 5, 6, 7, 8 } }, b = { { -1, -2, -3, -4 }, { -5, -6, -7, -8 } };
        double[,] h = design.Hybrid(a, b, design.IndexOf([0, 3]));
        Assert.Equal(new double[,] { { -1, 2, 3, -4 }, { -5, 6, 7, -8 } }, h);

        Assert.Throws<ArgumentOutOfRangeException>(() => new PickFreezeDesign(0));
        Assert.Throws<ArgumentException>(() => new PickFreezeDesign(3, [[0, 0]]));
        Assert.Throws<ArgumentException>(() => new PickFreezeDesign(3, [[3]]));
        Assert.Throws<ArgumentException>(() => new PickFreezeDesign(3, [[]]));
        Assert.Throws<ArgumentException>(() => new PickFreezeDesign(12, [[.. Enumerable.Range(0, 12)]], budget: 100));
        double[] two = [1.0, 2.0];
        Assert.Throws<ArgumentException>(() => PickFreeze.Estimate(new PickFreezeDesign(1), [1.0], [2.0], [[3.0]]));
        Assert.Throws<ArgumentException>(() => PickFreeze.Estimate(new PickFreezeDesign(1), two, two, [[3.0, double.NaN]]));
        Assert.Throws<ArgumentException>(() => PickFreeze.Estimate(new PickFreezeDesign(2), two, two, [two]));

        // The two-difference form reads each block's reverse hybrid, B's rows with u's columns from A.
        var reversed = new PickFreezeDesign(4, [[2, 0, 3], [1, 3]], MainVarianceForm.TwoDifference);
        Assert.Equal(MainVarianceForm.TwoDifference, reversed.Form);
        Assert.Equal(4096L * 20, reversed.Values(4096));
        Assert.Equal(new double[,] { { 1, -2, -3, 4 }, { 5, -6, -7, 8 } }, reversed.Reverse(a, b, reversed.IndexOf([0, 3])));
        var one = new PickFreezeDesign(1, form: MainVarianceForm.TwoDifference);
        Assert.Throws<ArgumentException>(() => PickFreeze.Estimate(one, two, two, [two]));                // no reverse hybrids
        Assert.Throws<ArgumentException>(() => PickFreeze.Estimate(one, two, two, [two], [[1.0]]));      // a short one
        Assert.Throws<ArgumentException>(() => PickFreeze.Estimate(new PickFreezeDesign(1), two, two, [two], [two]));   // Saltelli's reads none
        Assert.Throws<ArgumentOutOfRangeException>(() => new PickFreezeDesign(2, form: (MainVarianceForm)5));
    }

    [Fact]
    public void TwoDifference_ReadsSmallIndicesToTheirOwnScale()
    {
        // The slice 0 probe's finding, kept as a gate: on the g-function's S = 7.2·10⁻⁵ input, the two-difference form's spread is a few
        // percent of the piece (Saltelli's, printed beside it, is about three times the piece); both are unbiased (the theory above).
        SobolReference f = SobolReferences.G();
        const int N = 4096, Replicates = 100;
        var spread = new Dictionary<MainVarianceForm, double>();
        foreach (MainVarianceForm form in new[] { MainVarianceForm.Saltelli, MainVarianceForm.TwoDifference })
        {
            var design = new PickFreezeDesign(f.Dimension, null, form);
            double[] v5 = [.. SeedTree.Derive(99, Replicates).Select(seed =>
            {
                (double[,] a, double[,] b) = Bases(f, N, seed);
                return PickFreeze.Estimate(design, f.Value, a, b).Main[4];
            })];
            spread[form] = Deviation(v5) / f.MainVariance([4]);
        }
        output.WriteLine($"V_5 (S = 7.2e-5) at N = {N}: relative spread Saltelli {spread[MainVarianceForm.Saltelli]:F3}, two-difference {spread[MainVarianceForm.TwoDifference]:F4}");
        Assert.True(spread[MainVarianceForm.TwoDifference] <= 0.06, $"{spread[MainVarianceForm.TwoDifference]:F4}");
    }

    [Fact]
    public void Estimate_FromValuesMeetsTheFunctionsForm()
    {
        // The values route, given f's values on the hybrids, is the function route bit for bit.
        SobolReference f = SobolReferences.G();
        var design = new PickFreezeDesign(f.Dimension, [[0, 1], [4, 5]]);
        (double[,] a, double[,] b) = Bases(f, 512, 4);
        double[] Values(double[,] x) => [.. Enumerable.Range(0, x.GetLength(0)).Select(j => f.Value(Row(x, j)))];
        PickFreezeEstimate byValues = PickFreeze.Estimate(design, Values(a), Values(b),
            [.. Enumerable.Range(0, design.Blocks.Count).Select(k => Values(design.Hybrid(a, b, k)))]);
        PickFreezeEstimate byFunction = PickFreeze.Estimate(design, f.Value, a, b);
        Assert.Equal(byFunction.Main, byValues.Main);
        Assert.Equal(byFunction.MainError, byValues.MainError);
        Assert.Equal(byFunction.Total, byValues.Total);
        Assert.Equal(byFunction.TotalIndexError, byValues.TotalIndexError);
        Assert.Equal(byFunction.Variance, byValues.Variance);
    }

    [Fact]
    public void BelowZero_IsFlaggedOnlyBeyondItsError()
    {
        // An additive pair's interaction is 0, and its estimate rounds about 0: never flagged. A hand-made negative far outside its error is.
        SobolReference f = SobolReferences.Ishigami();
        var design = new PickFreezeDesign(3, [[0, 1]]);
        (double[,] a, double[,] b) = Bases(f, 1024, 9);
        Assert.Equal(PickFreezeStatus.Defined, PickFreeze.Estimate(design, f.Value, a, b).Status);
        double[] A = [.. Enumerable.Range(0, 1024).Select(j => (double)(j % 7))];
        double[] B = [.. A.Select(x => -x)];
        // B = −A and the hybrid's values are A − B = 2A, so D = A and the main variance is Cov(−A, A) = −Var A, far below zero for its error.
        PickFreezeEstimate e = PickFreeze.Estimate(new PickFreezeDesign(1), A, B, [[.. A.Select((x, j) => x - B[j])]]);
        Assert.True(e.Main[0] < -3.0 * e.MainError[0]);
        Assert.Equal(PickFreezeStatus.BelowZeroBeyondError, e.Status);
    }

    // A random reference with its truth, scaled and shifted.
    internal static SobolReference SweepCase(Xoshiro256PlusPlus rng)
    {
        SobolReference inner;
        int kind = rng.NextInt(3);
        if (kind == 0)
        {
            int p = 1 + rng.NextInt(8);
            inner = new SobolGReference([.. Enumerable.Range(0, p).Select(_ => rng.NextDouble() < 0.2 ? 0.0 : Math.Pow(10.0, -3.0 + 7.0 * rng.NextDouble()))]);
        }
        else if (kind == 1)
            inner = new IshigamiReference(20.0 * rng.NextDouble(), rng.NextDouble() < 0.2 ? 0.0 : 0.5 * rng.NextDouble());
        else
        {
            int p = 1 + rng.NextInt(5), count = 1 + rng.NextInt(3);
            double[] center = [.. Enumerable.Range(0, p).Select(_ => (rng.NextDouble() < 0.5 ? -1.0 : 1.0) * Math.Pow(10.0, -3.0 + 6.0 * rng.NextDouble()))];
            double[] half = [.. Enumerable.Range(0, p).Select(_ => 0.5 * Math.Pow(10.0, -3.0 + 6.0 * rng.NextDouble()))];
            var law = new ReferenceBox([.. center.Select((m, j) => m - half[j])], [.. center.Select((m, j) => m + half[j])]);
            var ridges = new RidgeTerm[count];
            for (int k = 0; k < count; k++)
            {
                double[] w = new double[p];
                int first = rng.NextInt(p);
                for (int j = 0; j < p; j++) w[j] = j == first || rng.NextDouble() < 0.6 ? Normal(rng) : 0.0;
                double norm = Math.Sqrt(w.Sum(x => x * x));
                for (int j = 0; j < p; j++) w[j] /= norm;
                double reach = Enumerable.Range(0, p).Sum(j => Math.Abs(w[j]) * half[j]);   // the projection's half-range
                double mid = Enumerable.Range(0, p).Sum(j => w[j] * center[j]);
                ExponentialProfile profile = rng.NextDouble() < 0.5
                    ? ExponentialProfile.Exp((-2.0 + 4.0 * rng.NextDouble()) / reach, 0.5 + rng.NextDouble())
                    : rng.NextDouble() < 0.5 ? ExponentialProfile.Sin((0.5 + 5.5 * rng.NextDouble()) / reach) : ExponentialProfile.Cos((0.5 + 5.5 * rng.NextDouble()) / reach);
                ridges[k] = new RidgeTerm(profile, w, mid + reach * (-1.0 + 2.0 * rng.NextDouble()));
            }
            inner = new RidgeSobolReference($"ridges ({count} on {p} inputs)", new RidgeSum("sweep ridges", ridges), law);
        }
        double scale = Math.Pow(10.0, -100.0 + 200.0 * rng.NextDouble());
        double shift = (rng.NextDouble() < 0.5 ? -1.0 : 1.0) * Math.Pow(10.0, -2.0 + 11.0 * rng.NextDouble()) * Math.Sqrt(inner.Variance);
        return new AffineReference(inner, scale, scale * shift);
    }

    private static double Normal(Xoshiro256PlusPlus rng)
    {
        double u = 1.0 - rng.NextDouble(), v = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u)) * Math.Cos(2.0 * Math.PI * v);
    }

    internal static double[] Row(double[,] x, int j) => [.. Enumerable.Range(0, x.GetLength(1)).Select(v => x[j, v])];

    internal static double Deviation(double[] x)
    {
        double m = x.Average();
        return Math.Sqrt(x.Sum(v => (v - m) * (v - m)) / (x.Length - 1));
    }

    private static double Slope(double[] x, double[] y)
    {
        double mx = x.Average(), my = y.Average();
        return x.Zip(y).Sum(t => (t.First - mx) * (t.Second - my)) / x.Sum(t => (t - mx) * (t - mx));
    }
}

/// <summary>s·f + c for a reference f: every variance piece times s², the indices as f's.</summary>
internal sealed class AffineReference(SobolReference inner, double scale, double shift) : SobolReference($"{inner.Name} × {scale:E1} + {shift:E1}", inner.Law)
{
    public override double Value(ReadOnlySpan<double> x) => scale * inner.Value(x) + shift;

    public override double Mean => scale * inner.Mean + shift;

    public override double ClosedVariance(int[] u) => scale * scale * inner.ClosedVariance(u);

    public override double MainVariance(int[] u) => scale * scale * inner.MainVariance(u);

    public override double TotalVariance(int i) => scale * scale * inner.TotalVariance(i);

    public override double ConditionalMean(int[] u, ReadOnlySpan<double> values) => scale * inner.ConditionalMean(u, values) + shift;

    public override bool Varies(int v) => inner.Varies(v);
}
