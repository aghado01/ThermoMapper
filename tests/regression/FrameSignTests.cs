using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 7e gates (the frame's column signs). Every eigenvector the subspace readouts return carries <see cref="Pca.Orient(double[])"/>'s
/// sign, its largest-magnitude entry positive, the first on ties, so a solver that returns a column negated moves nothing:
/// <list type="bullet">
/// <item>the rule on a seeded sweep of symmetric matrices over twelve decades: each column negated in turn, the oriented vectors the same bits;</item>
/// <item>the paired form and the written-out copies it replaces (<see cref="BassDecomposition"/>'s components, <see cref="OutputPosterior"/>'s
/// left vectors and leading patterns) the same bits, ties included;</item>
/// <item>the posterior's vectors from a solve with a column negated the same bits, and the rule on every vector the posteriors return;</item>
/// <item>slice 7b's ridge: the frame from a solve with its column negated the same bits and the refit's draws the same; the frame negated as
/// a sign-blind copy would read it moves the draws, and the comparison refuses the analysis for it.</item>
/// </list>
/// </summary>
public sealed class FrameSignTests(ITestOutputHelper output)
{
    /// <summary>300 symmetric matrices, p from 2 to 10, entries scaled over 10⁻⁶ to 10⁶: the solver's vectors oriented, and the same
    /// vectors with each column negated in turn, oriented, are the same bits; each oriented column's largest entry is positive.</summary>
    [Fact]
    public void OrientColumns_ColumnNegated_SameBits()
    {
        var rng = new Xoshiro256PlusPlus(20261113);
        int columns = 0, turned = 0;
        for (int trial = 0; trial < 300; trial++)
        {
            int p = 2 + (int)(rng.NextDouble() * 9);
            double[,] raw = Solve(RandomSymmetric(rng, p, Math.Pow(10.0, -6.0 + 12.0 * rng.NextDouble())));
            double[,] oriented = (double[,])raw.Clone();
            Pca.OrientColumns(oriented);
            for (int j = 0; j < p; j++)
            {
                columns++;
                if (Largest(raw, j) < 0.0) turned++;
                Assert.True(Largest(oriented, j) > 0.0, $"trial {trial}, column {j}");
                double[,] negated = Negated(raw, j);
                Pca.OrientColumns(negated);
                AssertSameBits(oriented, negated);
            }
        }
        output.WriteLine($"{columns} columns, {turned} of them turned by the rule; every negated solve the same bits");
    }

    /// <summary>Ties go to the first entry; a zero column keeps its bits; the column form, the axis form and the paired form agree.</summary>
    [Fact]
    public void Orient_Ties_FirstEntryDecides()
    {
        foreach ((double[] axis, double[] expected) in new[]
        {
            (new[] { 0.6, -0.6, 0.5 }, new[] { 0.6, -0.6, 0.5 }),
            (new[] { -0.6, 0.6, 0.5 }, new[] { 0.6, -0.6, -0.5 }),
            (new[] { 0.0, -0.0, 0.0 }, new[] { 0.0, -0.0, 0.0 }),
            (new[] { 0.1, -0.9, 0.9 }, new[] { -0.1, 0.9, -0.9 }),
        })
        {
            double[] single = [.. axis];
            Pca.Orient(single);
            Assert.Equal(expected.Select(BitConverter.DoubleToInt64Bits), single.Select(BitConverter.DoubleToInt64Bits));
            var column = new double[axis.Length, 1];
            for (int i = 0; i < axis.Length; i++) column[i, 0] = axis[i];
            Pca.OrientColumns(column);
            for (int i = 0; i < axis.Length; i++) Assert.Equal(BitConverter.DoubleToInt64Bits(expected[i]), BitConverter.DoubleToInt64Bits(column[i, 0]));
            double[] paired = [.. axis], partner = [1.0, -2.0];
            Pca.Orient(paired, partner);
            Assert.Equal(single.Select(BitConverter.DoubleToInt64Bits), paired.Select(BitConverter.DoubleToInt64Bits));
            bool turned = !single.Select(BitConverter.DoubleToInt64Bits).SequenceEqual(axis.Select(BitConverter.DoubleToInt64Bits));
            Assert.Equal(turned ? new[] { -1.0, 2.0 } : new[] { 1.0, -2.0 }, partner);
        }
    }

    /// <summary>2,000 vector pairs whose entries are drawn from a few values, so ties and zeros are common: the paired form against
    /// <see cref="BassDecomposition"/>'s written-out loop, and the axis form against <see cref="OutputPosterior"/>'s (its Left and Signed
    /// held the same loop), the same bits.</summary>
    [Fact]
    public void PairedForm_AndTheWrittenOutCopies_SameBits()
    {
        var rng = new Xoshiro256PlusPlus(71113);
        double[] values = [-2.0, -1.0, -0.5, -0.0, 0.0, 0.5, 1.0, 2.0];
        for (int trial = 0; trial < 2_000; trial++)
        {
            int m = 1 + (int)(rng.NextDouble() * 6), n = 1 + (int)(rng.NextDouble() * 5);
            double[] left = [.. Enumerable.Range(0, m).Select(_ => values[(int)(rng.NextDouble() * values.Length)])];
            double[] right = [.. Enumerable.Range(0, n).Select(_ => values[(int)(rng.NextDouble() * values.Length)] * (1.0 + rng.NextDouble()))];

            double[] oldLeft = [.. left], oldRight = [.. right];   // BassDecomposition's loop before slice 7e
            int largest = 0;
            for (int l = 1; l < m; l++)
                if (Math.Abs(oldLeft[l]) > Math.Abs(oldLeft[largest])) largest = l;
            if (!(oldLeft[largest] >= 0.0))
            {
                for (int l = 0; l < m; l++) oldLeft[l] = -oldLeft[l];
                for (int i = 0; i < n; i++) oldRight[i] = -oldRight[i];
            }
            double[] newLeft = [.. left], newRight = [.. right];
            Pca.Orient(newLeft, newRight);
            Assert.Equal(oldLeft.Select(BitConverter.DoubleToInt64Bits), newLeft.Select(BitConverter.DoubleToInt64Bits));
            Assert.Equal(oldRight.Select(BitConverter.DoubleToInt64Bits), newRight.Select(BitConverter.DoubleToInt64Bits));

            double[] oldU = [.. left];   // OutputPosterior.Left's and Signed's loop before slice 7e
            int big = 0;
            for (int l = 1; l < m; l++)
                if (Math.Abs(oldU[l]) > Math.Abs(oldU[big])) big = l;
            if (oldU[big] < 0.0)
                for (int l = 0; l < m; l++) oldU[l] = -oldU[l];
            double[] newU = [.. left];
            Pca.Orient(newU);
            Assert.Equal(oldU.Select(BitConverter.DoubleToInt64Bits), newU.Select(BitConverter.DoubleToInt64Bits));
        }
    }

    /// <summary>
    /// Four draws of p = 5: E[C]'s and P̄'s vectors from a solve with each column negated in turn, oriented, are the posterior's bits, and
    /// so are each draw's from C and from a factor. Every vector the posteriors return keeps the rule: E[C]'s and P̄'s, each draw's, the
    /// co-active directions and the co-active projection mean's vectors and bases, and the rank posterior's.
    /// </summary>
    [Fact]
    public void Posteriors_FromANegatedSolve_SameBits()
    {
        var rng = new Xoshiro256PlusPlus(5113);
        const int p = 5, n = 40;
        double[][,] grads = [.. Enumerable.Range(0, 5).Select(_ => RandomMatrix(rng, n, p))];
        double[][,] cs = [.. grads.Take(4).Select(g => Gram(g, g))];
        double[] weights = [1.0, 2.0, 0.5, 1.5];
        SubspacePosterior posterior = SubspacePosterior.Of(cs, weights, 2);

        foreach ((double[,] m, double[,] vectors) in new[] { (posterior.MeanC, posterior.MeanCVectors), (posterior.Projection, posterior.ProjectionVectors) })
        {
            double[,] raw = Solve(m);
            for (int j = 0; j < p; j++)
            {
                double[,] negated = Negated(raw, j);
                Pca.OrientColumns(negated);
                AssertSameBits(vectors, negated);
            }
        }
        for (int s = 0; s < cs.Length; s++)
        {
            double[,] negated = Negated(Solve(cs[s]), s % p);
            Pca.OrientColumns(negated);
            AssertSameBits(posterior.Draws[s].Vectors, negated);

            DrawSubspace fromFactor = DrawSubspace.FromFactor(grads[s], cs[s], 2);
            ThinSvdResult svd = ThinSvd.Decompose(grads[s], leftVectors: false);
            var right = new double[p, svd.RightVectors.Length];
            for (int j = 0; j < svd.RightVectors.Length; j++)
                for (int i = 0; i < p; i++) right[i, j] = svd.RightVectors[j][i];
            double[,] flipped = Negated(right, s % svd.RightVectors.Length);
            Pca.OrientColumns(flipped);
            AssertSameBits(fromFactor.Vectors, flipped);
            AssertRule(fromFactor.Basis);
        }

        var pairs = new[] { new CoActivePair(cs[0], cs[1], Gram(grads[0], grads[1])), new CoActivePair(cs[2], cs[3], Gram(grads[2], grads[3])) };
        CoActivePosterior co = CoActivePosterior.Of(pairs, [1.0, 3.0]);
        CoActiveProjection mean = co.ProjectionMean(2);
        RankPosterior rank = RankPosterior.Of(posterior.Draws, weights, RankRule.TraceShare(0.9));
        foreach (double[,] vectors in new[] { posterior.MeanCVectors, posterior.ProjectionVectors, co.Directions, mean.Vectors, rank.ProjectionVectors }
                                      .Concat(posterior.Draws.Select(d => d.Vectors)).Concat(posterior.Draws.Select(d => d.Basis)).Concat(mean.Bases))
            AssertRule(vectors);
    }

    /// <summary>
    /// Slice 7b's ridge at d = 1. The frame from E[C]'s solve with any column negated, oriented, is the fit's frame bit for bit, and a
    /// projection read from it refits to the same draws. The frame negated, as a copy blind to the sign would read it, refits to other
    /// draws, and the comparison refuses the analysis for that refit while it accepts its own.
    /// </summary>
    [Fact]
    public void Ridge_FrameFromANegatedSolve_SameBitsAndDraws()
    {
        BassProjectedFitTests.RidgeFits f = BassProjectedFitTests.Fixture.Value;
        SubspacePosterior posterior = f.Analysis.Posterior;
        double[,] raw = Solve(posterior.MeanC);
        int p = raw.GetLength(0);
        output.WriteLine($"the solver's first column {(Largest(raw, 0) < 0.0 ? "has its largest entry negative: the rule turned it" : "already keeps the rule")}");
        for (int j = 0; j < p; j++)
        {
            double[,] negated = Negated(raw, j);
            Pca.OrientColumns(negated);
            AssertSameBits(posterior.MeanCVectors, negated);
        }

        double[,] turned = Negated(raw, 0);
        Pca.OrientColumns(turned);
        BassProjection projection = BassProjection.FromVectors(f.Native.Transform, f.Analysis, turned, 1, null);
        AssertSameBits(f.Projection.Frame, projection.Frame);
        BassProjectedFit refit = BassProjectedFit.Of(projection, f.Data, BassProjectedFitTests.Config());
        BassRunTests.AssertSameFit(f.Refit.Fit, refit.Fit);
        refit.Compare(f.Native, f.Data, f.Analysis, new ProjectionComparisonSpec());

        BassProjection reflected = BassProjection.FromVectors(f.Native.Transform, f.Analysis, Negated(posterior.MeanCVectors, 0), 1, null);
        BassProjectedFit moved = BassProjectedFit.Of(reflected, f.Data, BassProjectedFitTests.Config());
        BassDraw[] a = [.. f.Refit.Fit.Draws.SelectMany(c => c)], b = [.. moved.Fit.Draws.SelectMany(c => c)];
        int differ = Enumerable.Range(0, a.Length).Count(s => !a[s].Coefficients.SequenceEqual(b[s].Coefficients));
        output.WriteLine($"the frame negated: {differ} of {a.Length} draws' coefficients differ");
        Assert.True(differ > 0);
        Assert.Throws<ArgumentException>(() => moved.Compare(f.Native, f.Data, f.Analysis, new ProjectionComparisonSpec()));
    }

    // ── Helpers ──────────────────────────────────────────────────────────────────────────────────────────────────────────

    // The solve SubspacePosterior reads, before the rule.
    private static double[,] Solve(double[,] m)
    {
        int p = m.GetLength(0);
        var values = new double[p];
        var vectors = new double[p, p];
        DenseEigen.DecomposeSymmetricInto(m, values, vectors, tol: 1e-14);
        return vectors;
    }

    private static double[,] Negated(double[,] vectors, int column)
    {
        var copy = (double[,])vectors.Clone();
        for (int i = 0; i < copy.GetLength(0); i++) copy[i, column] = -copy[i, column];
        return copy;
    }

    // The column's largest-magnitude entry, the first on ties.
    private static double Largest(double[,] vectors, int column)
    {
        int best = 0;
        for (int i = 1; i < vectors.GetLength(0); i++)
            if (Math.Abs(vectors[i, column]) > Math.Abs(vectors[best, column])) best = i;
        return vectors[best, column];
    }

    private static void AssertRule(double[,] vectors)
    {
        for (int j = 0; j < vectors.GetLength(1); j++) Assert.True(Largest(vectors, j) >= 0.0, $"column {j}: {Largest(vectors, j):R}");
    }

    private static void AssertSameBits(double[,] expected, double[,] actual)
    {
        Assert.Equal(expected.GetLength(0), actual.GetLength(0));
        Assert.Equal(expected.GetLength(1), actual.GetLength(1));
        for (int i = 0; i < expected.GetLength(0); i++)
            for (int j = 0; j < expected.GetLength(1); j++)
                Assert.True(BitConverter.DoubleToInt64Bits(expected[i, j]) == BitConverter.DoubleToInt64Bits(actual[i, j]),
                            $"[{i}, {j}]: {expected[i, j]:R} against {actual[i, j]:R}");
    }

    private static double[,] RandomSymmetric(Xoshiro256PlusPlus rng, int p, double scale)
    {
        var m = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j <= i; j++) m[i, j] = m[j, i] = scale * (2.0 * rng.NextDouble() - 1.0);
        return m;
    }

    private static double[,] RandomMatrix(Xoshiro256PlusPlus rng, int rows, int columns)
    {
        var m = new double[rows, columns];
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < columns; c++) m[r, c] = (2.0 * rng.NextDouble() - 1.0) * (1.0 + c);
        return m;
    }

    // AᵀB / N for two N × p gradient matrices.
    private static double[,] Gram(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), p = a.GetLength(1);
        var g = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
            {
                double s = 0.0;
                for (int r = 0; r < n; r++) s += a[r, i] * b[r, j];
                g[i, j] = s / n;
            }
        return g;
    }
}
