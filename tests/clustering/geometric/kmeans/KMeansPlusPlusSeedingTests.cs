using System;
using System.Linq;
using Maths.Rng;
using Xunit;

namespace Clustering.Geometric.KMeans.Tests;

public class KMeansPlusPlusSeedingTests
{
    // Four blobs of 50 points on the corners of a square of half-side 50, unit spread: blob g holds rows 50g..50g+49.
    private static double[][] FarBlobs(int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var data = new double[200][];
        for (int i = 0; i < 200; i++)
        {
            int g = i / 50;
            data[i] = new[]
            {
                ((g & 1) == 1 ? 50.0 : -50.0) + 2.0 * rng.NextDouble() - 1.0,
                ((g & 2) == 2 ? 50.0 : -50.0) + 2.0 * rng.NextDouble() - 1.0,
            };
        }
        return data;
    }

    [Fact]
    public void Seeds_AreDistinctDataRows()
    {
        double[][] data = FarBlobs(seed: 1);
        double[][] seeds = KMeansPlusPlusSeeding.Seed(data, 20, new Xoshiro256PlusPlus(2));

        Assert.Equal(20, seeds.Length);
        int[] rows = seeds.Select(seed => Array.FindIndex(data, row => row.SequenceEqual(seed))).ToArray();
        Assert.All(rows, row => Assert.True(row >= 0));
        Assert.Equal(20, rows.Distinct().Count());
        Assert.All(seeds, seed => Assert.DoesNotContain(data, row => ReferenceEquals(row, seed)));
    }

    [Fact]
    public void Seeds_AreAFunctionOfTheGeneratorState()
    {
        double[][] data = FarBlobs(seed: 3);
        double[][] first = KMeansPlusPlusSeeding.Seed(data, 6, new Xoshiro256PlusPlus(4));
        double[][] second = KMeansPlusPlusSeeding.Seed(data, 6, new Xoshiro256PlusPlus(4));

        for (int c = 0; c < 6; c++) Assert.Equal(first[c], second[c]);
    }

    [Fact]
    public void FarSeparatedBlobs_TakeOneSeedEach()
    {
        // A second seed in a covered blob has probability about 1e-4 per draw at this separation.
        double[][] data = FarBlobs(seed: 5);
        double[][] seeds = KMeansPlusPlusSeeding.Seed(data, 4, new Xoshiro256PlusPlus(6));

        int[] blobs = seeds.Select(seed => (seed[0] > 0 ? 1 : 0) + (seed[1] > 0 ? 2 : 0)).ToArray();
        Assert.Equal(4, blobs.Distinct().Count());
    }

    [Fact]
    public void CoincidentRows_StillYieldKSeeds()
    {
        double[][] data = Enumerable.Range(0, 10).Select(_ => new[] { 1.5, -2.5 }).ToArray();
        double[][] seeds = KMeansPlusPlusSeeding.Seed(data, 3, new Xoshiro256PlusPlus(7));

        Assert.Equal(3, seeds.Length);
        Assert.All(seeds, seed => Assert.Equal(new[] { 1.5, -2.5 }, seed));
    }

    [Fact]
    public void MoreSeedsThanRows_IsRefused()
        => Assert.Throws<ArgumentOutOfRangeException>(
            () => KMeansPlusPlusSeeding.Seed(FarBlobs(seed: 8), 201, new Xoshiro256PlusPlus(9)));
}
