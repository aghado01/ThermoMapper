using System;
using System.Collections.Generic;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using static Maths.LinAlg.Tests.FactorReferenceTests;

namespace Maths.LinAlg.Tests;

/// <summary>
/// <see cref="BorderedCholesky.Copy"/>: a copy taken after births, deaths, changes and rebuilds continues through the
/// same updates bit for bit with the original, rebuilding on the same cadence, and shares nothing with it; a pending
/// proposal stays behind.
/// </summary>
public sealed class BorderedCholeskyCopyTests
{
    private const int Rows = 60, Capacity = 16;

    [Fact]
    public void Copy_ContinuesBitForBit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 811);
        var original = new BorderedCholesky(Capacity, rebuildInterval: 25);
        var active = new List<double[]>();
        for (int step = 0; step < 160; step++) Update(original, active, rng);
        Assert.True(original.UpdatesSinceRebuild > 0);

        BorderedCholesky copy = original.Copy();
        var copied = new List<double[]>(active);
        AssertSame(original, copy);
        // The same updates from two equal streams.
        (ulong, ulong, ulong, ulong) words = rng.SaveState();
        var copyRng = Xoshiro256PlusPlus.Restore(words.Item1, words.Item2, words.Item3, words.Item4);
        for (int step = 0; step < 200; step++)
        {
            Update(original, active, rng);
            Update(copy, copied, copyRng);
            AssertSame(original, copy);
        }
    }

    [Fact]
    public void Copy_SharesNothing_AndLeavesAPendingProposal()
    {
        var rng = new Xoshiro256PlusPlus(seed: 823);
        var original = new BorderedCholesky(Capacity);
        var active = new List<double[]>();
        for (int step = 0; step < 6; step++) Update(original, active, rng);
        double[] c = Vector(Rows, rng);
        original.ProposeBirth(Border(active, c), Dot(c, c));

        BorderedCholesky copy = original.Copy();
        Assert.Null(copy.Pending);
        Assert.Throws<InvalidOperationException>(copy.Accept);
        double before = original.Factor(0, 0);
        copy.ProposeDeath(0);
        copy.Accept();
        Assert.Equal(before, original.Factor(0, 0));
        Assert.Equal(active.Count, original.Dimension);
        Assert.Equal(active.Count - 1, copy.Dimension);
        original.Accept();   // its own proposal is still pending
        Assert.Equal(active.Count + 1, original.Dimension);
    }

    // A random birth, death or change, accepted.
    private static void Update(BorderedCholesky f, List<double[]> active, Xoshiro256PlusPlus rng)
    {
        int n = active.Count;
        int kind = n < 3 ? 0 : n == Capacity ? 1 : rng.NextInt(3);
        if (kind == 0)
        {
            double[] c = Vector(Rows, rng);
            f.ProposeBirth(Border(active, c), Dot(c, c));
            active.Add(c);
        }
        else
        {
            int m = rng.NextInt(n);
            active.RemoveAt(m);
            if (kind == 1) f.ProposeDeath(m);
            else
            {
                double[] c = Vector(Rows, rng);
                f.ProposeChange(m, Border(active, c), Dot(c, c));
                active.Add(c);
            }
        }
        f.Accept();
    }

    private static void AssertSame(BorderedCholesky a, BorderedCholesky b)
    {
        Assert.Equal(a.Capacity, b.Capacity);
        Assert.Equal(a.RebuildInterval, b.RebuildInterval);
        Assert.Equal(a.Dimension, b.Dimension);
        Assert.Equal(a.UpdatesSinceRebuild, b.UpdatesSinceRebuild);
        Assert.Equal(a.LogDet, b.LogDet);
        Assert.Equal(a.Status, b.Status);
        Assert.Equal(a.LastDrift, b.LastDrift);
        for (int i = 0; i < a.Dimension; i++)
        {
            Assert.Equal(a.InverseDiagonal(i), b.InverseDiagonal(i));
            for (int j = 0; j <= i; j++)
            {
                Assert.Equal(a.Factor(i, j), b.Factor(i, j));
                Assert.Equal(a.Entry(i, j), b.Entry(i, j));
            }
        }
    }

    private static double[] Border(List<double[]> columns, double[] c)
    {
        var u = new double[columns.Count];
        for (int j = 0; j < u.Length; j++) u[j] = Dot(columns[j], c);
        return u;
    }

    private static double Dot(double[] a, double[] b) => PairwiseSum.Dot(a, b);
}
