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
            PairedFraction = 0.2,
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
    [InlineData(0.0, 0.0)]   // pure single-column Givens
    [InlineData(0.1, 0.0)]   // default mixture
    [InlineData(1.0, 0.0)]   // isotropic-only
    [InlineData(0.0, 1.0)]   // paired-only
    [InlineData(0.3, 0.3)]   // all three kinds
    public void Compute_Result_RowsAreOrthonormal(double isotropicFraction, double pairedFraction)
    {
        double[][] data = BuildData(samples: 40, dim: 4, seed: 11);
        var options = new SubspaceAnnealerOptions
        {
            IsotropicFraction = isotropicFraction,
            PairedFraction = pairedFraction,
        };

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

    // One row per validation rule, NaN cases prominent: the original relational checks were the
    // hole NaN slipped through (a NaN TargetAcceptance turns every step scale into NaN; a NaN
    // InitialTemperature silently freezes Metropolis into greedy descent). The record's ToString
    // identifies a failing row.
    public static TheoryData<SubspaceAnnealerOptions> InvalidOptions => new()
    {
        new() { IsotropicFraction = double.NaN },
        new() { IsotropicFraction = -0.1 },
        new() { IsotropicFraction = 1.1 },
        new() { TargetAcceptance = double.NaN },
        new() { TargetAcceptance = 0.0 },
        new() { TargetAcceptance = 1.0 },
        new() { InitialStep = double.NaN },
        new() { InitialStep = double.PositiveInfinity },
        new() { InitialStep = 0.0 },
        new() { StepFloor = double.NaN },
        new() { StepFloor = 0.0 },
        new() { StepCeiling = double.NaN },
        new() { StepCeiling = double.PositiveInfinity },
        new() { StepFloor = 0.5, StepCeiling = 0.1 },
        new() { InitialTemperature = double.NaN },
        new() { InitialTemperature = double.PositiveInfinity },
        new() { InitialTemperature = 0.0 },
        new() { CoolingRate = double.NaN },
        new() { CoolingRate = 0.0 },
        new() { CoolingRate = 1.1 },
        new() { PairedFraction = double.NaN },
        new() { PairedFraction = -0.1 },
        new() { PairedFraction = 1.1 },
        new() { IsotropicFraction = 0.6, PairedFraction = 0.6 },   // mixture over-full
    };

    [Theory]
    [MemberData(nameof(InvalidOptions))]
    public void Compute_InvalidOptions_ThrowsBeforeAnyEvaluation(SubspaceAnnealerOptions options)
    {
        double[][] data = BuildData(samples: 10, dim: 4, seed: 11);
        int evaluations = 0;
        double Counting(double[][] projection)
        {
            evaluations++;
            return Objective(projection);
        }

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SubspaceAnnealer.Compute(data, targetDim: 2, Counting, maxIters: 5, seed: 7, options));

        Assert.Equal(0, evaluations);
    }

    [Fact]
    public void Compute_PairedFractionWithTargetDimOne_Throws()
    {
        double[][] data = BuildData(samples: 10, dim: 4, seed: 11);
        var options = new SubspaceAnnealerOptions { PairedFraction = 0.5 };

        Assert.Throws<ArgumentException>(() =>
            SubspaceAnnealer.Compute(data, targetDim: 1, Objective, maxIters: 1, seed: 7, options));
    }

    [Fact]
    public void Compute_DegeneratePair_PairedMovesDescend_WhereSingleColumnCrawls()
    {
        // Finding 1's geometry made synthetic: rank-5 data on Gr(5, 200) whose 4th/5th variance
        // directions are the 45°-mixtures of the good direction e_3 and the complement defect
        // e_100, so the warm start spans {e_0..e_3, e_100} (objective 1 against the target
        // span{e_0..e_4}) with the defect smeared evenly across the near-degenerate pair. That
        // warm start is a saddle where strict single-column improvement is essentially impossible:
        // a rotation of either mixed column pays an O(1/2) loss of its good share against an
        // O(1/(d−k)) complement gain, and at 195 complement dimensions a draw carrying enough e_4
        // mass to beat that trade has vanishing probability (at small d the tail is fat enough to
        // leak — measured: d=30 lets single-column moves recover fully). The paired move's in-span
        // mixing angle reaches m ≈ e_100 and excises the defect with no good-share loss — and its
        // first success CONCENTRATES the pair basis ({≈e_3, ≈junk}), after which the ordinary
        // single-column share descends exactly as on a pure-offset column. The temperature is
        // pinned at 1e-6 (≪ any real increment) so thermal saddle diffusion — which SA temperature
        // legitimately provides, and which descends either arm at 1e-3 — cannot contribute: the
        // contrast isolates proposal geometry alone. The shared raised step floor keeps the
        // acceptance controller from buying acceptance by shrinking θ into the diffusion regime.
        const int d = 200, k = 5;
        double[][] data = DegeneratePairData(samples: 300, dim: d, seed: 5);

        static double DistanceToTarget(double[][] projection)
        {
            double overlap = 0.0;
            for (int r = 0; r < projection.Length; r++)
                for (int j = 0; j < 5; j++)
                    overlap += projection[r][j] * projection[r][j];
            return projection.Length - overlap;   // k − Σ cos²θ_i = Σ sin²θ_i
        }

        var singleOnly = new SubspaceAnnealerOptions
        {
            IsotropicFraction = 0.0,
            PairedFraction = 0.0,
            InitialStep = 0.5,
            StepFloor = 0.4,
            InitialTemperature = 1e-6,
        };
        var paired = singleOnly with { PairedFraction = 0.5 };

        SubspaceAnnealerResult single = SubspaceAnnealer.Compute(
            data, k, DistanceToTarget, maxIters: 12000, seed: 7, singleOnly);
        SubspaceAnnealerResult mixed = SubspaceAnnealer.Compute(
            data, k, DistanceToTarget, maxIters: 12000, seed: 7, paired);

        _out.WriteLine($"start ≈ 1: single-column={single.Objective:F6}  paired-mixture={mixed.Objective:F6}");

        // Measured at seed 7: single-column exactly 1.000000 (bit-frozen), paired 0.649 — and
        // super-linearly (0.951 at half budget): early descent is window-rate-limited because
        // accepted Δ≈0 uniform-φ rotations re-smear the pair basis between excision hits, until
        // enough e_4 share accumulates to open first-order channels and compound. The engine fact
        // pins the qualitative escape; efficiency at scale is the S0 probe's question.
        Assert.True(single.Objective > 0.99,
            $"Single-column moves should stay pinned at the smeared-defect saddle (≈ 1); got {single.Objective:F6} — " +
            "if they now descend, the degenerate-pair premise needs revisiting.");
        Assert.True(mixed.Objective < 0.8,
            $"The paired mixture should escape the saddle and descend decisively; got {mixed.Objective:F6}.");
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

    // Exact rank-5 data whose 4th/5th variance directions are the 45°-mixtures q1 = (e_3 + e_100)/√2
    // and q2 = (e_3 − e_100)/√2 with a deliberate ~10% eigengap between them: the gap pins PCA's
    // component choice to q1, q2 themselves (stable against sample noise, unlike an exactly
    // degenerate pair whose axes land at a sample-determined angle), so the warm start's last two
    // columns each carry half good direction e_3, half complement defect e_100 — maximal smearing,
    // the flat-eigentail geometry of mobility finding 1, deterministically.
    private static double[][] DegeneratePairData(int samples, int dim, int seed)
    {
        var rng = new Random(seed);
        const double invSqrt2 = 0.7071067811865476;
        var data = new double[samples][];
        for (int i = 0; i < samples; i++)
        {
            var row = new double[dim];
            row[0] = 2.0 * (rng.NextDouble() * 2.0 - 1.0);
            row[1] = 1.9 * (rng.NextDouble() * 2.0 - 1.0);
            row[2] = 1.8 * (rng.NextDouble() * 2.0 - 1.0);
            double a = 1.05 * (rng.NextDouble() * 2.0 - 1.0);   // weight on q1 = (e_3 + e_100)/√2
            double b = 0.95 * (rng.NextDouble() * 2.0 - 1.0);   // weight on q2 = (e_3 − e_100)/√2
            row[3] = invSqrt2 * (a + b);
            row[100] = invSqrt2 * (a - b);
            data[i] = row;
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
