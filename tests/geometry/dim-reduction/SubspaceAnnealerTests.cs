using System;
using System.Threading;
using Maths.Geometry.DimReduction;
using Xunit;

namespace Maths.Geometry.Tests;

public sealed class SubspaceAnnealerTests
{
    // A fixed, deterministic subspace objective (no RNG inside), so any non-determinism in the
    // result can only originate from the sampler's own RNG.
    private static double Objective(double[][] projection)
    {
        double c = 0.0;
        for (int r = 0; r < projection.Length; r++)
            for (int j = 0; j < projection[r].Length; j++)
                c += (r + 1) * j * projection[r][j] * projection[r][j];
        return c;
    }

    [Fact]
    public void Compute_SameSeed_IsBitIdentical()
    {
        double[][] data = BuildData(samples: 40, dim: 4, seed: 11);

        double[][] a = SubspaceAnnealer.Compute(data, targetDim: 2, Objective, maxIters: 300, seed: 7);
        double[][] b = SubspaceAnnealer.Compute(data, targetDim: 2, Objective, maxIters: 300, seed: 7);

        Assert.Equal(a.Length, b.Length);
        for (int i = 0; i < a.Length; i++)
        {
            Assert.Equal(a[i].Length, b[i].Length);
            for (int j = 0; j < a[i].Length; j++)
                Assert.Equal(a[i][j], b[i][j]); // same seed -> same xoshiro stream -> identical bits
        }
    }

    [Fact]
    public void Compute_Result_RowsAreOrthonormal()
    {
        double[][] data = BuildData(samples: 40, dim: 4, seed: 11);

        double[][] proj = SubspaceAnnealer.Compute(data, targetDim: 2, Objective, maxIters: 300, seed: 7);

        for (int i = 0; i < proj.Length; i++)
        {
            double selfDot = 0.0;
            for (int j = 0; j < proj[i].Length; j++) selfDot += proj[i][j] * proj[i][j];
            Assert.InRange(selfDot, 1.0 - 1e-9, 1.0 + 1e-9);

            for (int k = i + 1; k < proj.Length; k++)
            {
                double crossDot = 0.0;
                for (int j = 0; j < proj[i].Length; j++) crossDot += proj[i][j] * proj[k][j];
                Assert.InRange(crossDot, -1e-9, 1e-9);
            }
        }
    }

    [Fact]
    public void Compute_CancellationDuringAnneal_Throws()
    {
        double[][] data = BuildData(samples: 20, dim: 4, seed: 11);
        using var cancellation = new CancellationTokenSource();
        int evaluations = 0;

        double CancelAfterFirstProposal(double[][] projection)
        {
            if (++evaluations == 2) cancellation.Cancel();
            return Objective(projection);
        }

        Assert.Throws<OperationCanceledException>(() =>
            SubspaceAnnealer.Compute(
                data,
                targetDim: 2,
                CancelAfterFirstProposal,
                maxIters: 100,
                seed: 7,
                cancellation.Token));

        Assert.Equal(2, evaluations);
    }

    private static double[][] BuildData(int samples, int dim, int seed)
    {
        var rng = new Random(seed);
        double[][] data = new double[samples][];
        for (int i = 0; i < samples; i++)
        {
            data[i] = new double[dim];
            for (int j = 0; j < dim; j++)
                data[i][j] = rng.NextDouble() * 2.0 - 1.0;
        }
        return data;
    }
}
