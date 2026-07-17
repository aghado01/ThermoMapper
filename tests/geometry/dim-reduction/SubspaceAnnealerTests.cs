using System;
using System.Threading;
using Maths.Geometry.DimReduction;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Geometry.Tests;

public sealed class SubspaceAnnealerTests
{
    private readonly ITestOutputHelper _out;
    public SubspaceAnnealerTests(ITestOutputHelper output) => _out = output;

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

        SubspaceAnnealerResult a = SubspaceAnnealer.Compute(data, targetDim: 2, Objective, maxIters: 300, seed: 7);
        SubspaceAnnealerResult b = SubspaceAnnealer.Compute(data, targetDim: 2, Objective, maxIters: 300, seed: 7);

        Assert.Equal(a.Objective, b.Objective); // same seed -> same xoshiro stream -> identical bits
        AssertBitIdentical(a.Projection, b.Projection);
    }

    [Fact]
    public void Compute_SameSeedSameOptions_IsBitIdentical()
    {
        // Step adaptation is state-dependent (it reacts to the acceptance history), but that history
        // is a deterministic function of the one seeded stream — options must not loosen the contract.
        double[][] data = BuildData(samples: 40, dim: 4, seed: 11);
        var options = new SubspaceAnnealerOptions
        {
            IsotropicFraction = 0.5,
            TargetAcceptance = 0.3,
            StepCeiling = 1.0,
        };

        SubspaceAnnealerResult a = SubspaceAnnealer.Compute(data, targetDim: 2, Objective, maxIters: 300, seed: 7, options);
        SubspaceAnnealerResult b = SubspaceAnnealer.Compute(data, targetDim: 2, Objective, maxIters: 300, seed: 7, options);

        Assert.Equal(a.Objective, b.Objective);
        AssertBitIdentical(a.Projection, b.Projection);
    }

    [Fact]
    public void Compute_ReturnedObjective_MatchesFreshEvaluation()
    {
        double[][] data = BuildData(samples: 40, dim: 4, seed: 11);

        SubspaceAnnealerResult result = SubspaceAnnealer.Compute(data, targetDim: 2, Objective, maxIters: 300, seed: 7);

        // Callers reuse the tracked value in place of re-evaluating, which is only sound while a
        // fresh deterministic evaluation at the returned projection reproduces it bit-for-bit.
        Assert.Equal(Objective(result.Projection), result.Objective);
    }

    [Theory]
    [InlineData(0.0)]   // pure two-plane Givens
    [InlineData(0.1)]   // default mixture
    [InlineData(1.0)]   // isotropic-only
    public void Compute_Result_RowsAreOrthonormal(double isotropicFraction)
    {
        double[][] data = BuildData(samples: 40, dim: 4, seed: 11);
        var options = new SubspaceAnnealerOptions { IsotropicFraction = isotropicFraction };

        double[][] proj = SubspaceAnnealer.Compute(data, targetDim: 2, Objective, maxIters: 300, seed: 7, options).Projection;

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
    public void Compute_HighCodimension_TwoPlaneDescends_WhereIsotropicStalls()
    {
        // Quadratic-on-Grassmann mobility fact on Gr(5, 200) (intrinsic dimension 975): the
        // objective is the chordal distance² to the target subspace span{e_0..e_4} — Σ sin²θ_i,
        // minimized at 0. The data steers the PCA warm start onto span{e_0..e_3, e_100}: four
        // retained directions already aligned, the fifth deep in the complement (the ISOLET
        // situation — a good warm start offset in a sliver of the 975 intrinsic dimensions).
        // A full-rank isotropic tangent damages the four aligned columns at O(step²) while the
        // recovery sliver's directional derivative thins like 1/√dim, so isotropic proposals
        // stall: fixed-length ones freeze outright, and the adaptive controller alone only buys
        // acceptance by shrinking toward diffusion scale. A two-plane Givens rotation confines
        // the move to one column — proposals on aligned columns are rejected without damage and
        // proposals on the offset column keep an O(1) improving fraction — so it descends.
        const int d = 200, k = 5;
        double[][] data = OffsetColumnData(samples: 300, dim: d, seed: 5);

        static double DistanceToTarget(double[][] projection)
        {
            double overlap = 0.0;
            for (int r = 0; r < projection.Length; r++)
                for (int j = 0; j < 5; j++)
                    overlap += projection[r][j] * projection[r][j];
            return projection.Length - overlap;   // k − Σ cos²θ_i = Σ sin²θ_i
        }

        // All arms share a Metropolis temperature commensurate with the per-move objective
        // increments so the anneal refines the warm start instead of melting it — the fact pins
        // proposal mobility, with the temperature schedule controlled for.
        const double temp = 0.001;
        SubspaceAnnealerResult twoPlane = SubspaceAnnealer.Compute(
            data, k, DistanceToTarget, maxIters: 2000, seed: 7,
            new SubspaceAnnealerOptions { InitialTemperature = temp });   // otherwise shipped defaults: two-plane primary
        SubspaceAnnealerResult isotropicFixed = SubspaceAnnealer.Compute(
            data, k, DistanceToTarget, maxIters: 2000, seed: 7,
            new SubspaceAnnealerOptions
            {
                IsotropicFraction = 1.0,
                InitialStep = 0.1,
                StepFloor = 0.1,
                StepCeiling = 0.1,
                InitialTemperature = temp,
            });
        SubspaceAnnealerResult isotropicAdaptive = SubspaceAnnealer.Compute(
            data, k, DistanceToTarget, maxIters: 2000, seed: 7,
            new SubspaceAnnealerOptions { IsotropicFraction = 1.0, InitialTemperature = temp });

        _out.WriteLine(
            $"start ≈ 1: two-plane={twoPlane.Objective:F6}  isotropic-fixed={isotropicFixed.Objective:F6}  " +
            $"isotropic-adaptive={isotropicAdaptive.Objective:F6}");

        // Measured at seed 7: two-plane 0.642, both isotropic arms exactly 1.000000 (bit-frozen).
        Assert.True(twoPlane.Objective < 0.8,
            $"Two-plane proposals should descend decisively from the offset start (≈ 1); got {twoPlane.Objective:F6}.");
        Assert.True(isotropicFixed.Objective > 0.99,
            $"Fixed-length isotropic proposals should stall at the offset start (≈ 1); got {isotropicFixed.Objective:F6} — " +
            "if they now descend, the mobility-contrast premise needs revisiting.");
        Assert.True(isotropicAdaptive.Objective > 0.99,
            $"Adaptive-step isotropic proposals should still stall (the controller alone is not the fix); " +
            $"got {isotropicAdaptive.Objective:F6}.");
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
                cancellationToken: cancellation.Token));

        Assert.Equal(2, evaluations);
    }

    private static void AssertBitIdentical(double[][] a, double[][] b)
    {
        Assert.Equal(a.Length, b.Length);
        for (int i = 0; i < a.Length; i++)
        {
            Assert.Equal(a[i].Length, b[i].Length);
            for (int j = 0; j < a[i].Length; j++)
                Assert.Equal(a[i][j], b[i][j]);
        }
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

    // Signal variance on coords 0..3 and 100 (distinct amplitudes fix the PCA ordering), faint noise
    // elsewhere: the PCA warm start lands on span{e_0..e_3, e_100} — four columns aligned with the
    // mobility fact's target span{e_0..e_4}, one deep in the complement.
    private static double[][] OffsetColumnData(int samples, int dim, int seed)
    {
        var rng = new Random(seed);
        double[][] data = new double[samples][];
        for (int i = 0; i < samples; i++)
        {
            // Exactly rank-5 data: PCA's top-5 spans exactly {e_0..e_3, e_100} however the sample
            // eigenvectors rotate within it. The 1.7-vs-0.7 amplitude gap keeps the sample
            // eigenvalues decisively separated — with a near-degenerate pair PCA smears the offset
            // direction across two components, and no single-column move can trade a mixture apart.
            var row = new double[dim];
            row[0] = 2.0 * (rng.NextDouble() * 2.0 - 1.0);
            row[1] = 1.9 * (rng.NextDouble() * 2.0 - 1.0);
            row[2] = 1.8 * (rng.NextDouble() * 2.0 - 1.0);
            row[3] = 1.7 * (rng.NextDouble() * 2.0 - 1.0);
            row[100] = 0.4 * (rng.NextDouble() * 2.0 - 1.0);
            data[i] = row;
        }
        return data;
    }
}
