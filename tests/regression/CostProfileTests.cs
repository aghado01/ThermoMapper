using System;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The cost profile's facts (<c>issues/regression/performance-cost-regression-samplers-brief.md</c>). Each runs one block of
/// <see cref="CostProfile"/> behind <c>BASS_QUALIFICATION</c> = 1 and returns at once without it; <c>COST_PROFILE_BLOCK</c>
/// selects blocks by prefix (semicolon-separated) and <c>COST_PROFILE_OUT</c> names the folder the rows go to. A block whose
/// rows are present is skipped, so a rerun resumes.
/// </summary>
public sealed class CostProfileTests
{
    [Fact] public void Machine() => CostProfile.Run("machine", CostKernels.MachineBlock);

    [Fact] public void Kernels_Band() => CostProfile.Run("kernels:band", CostKernels.BandBlock);

    [Fact] public void Kernels_Engines() => CostProfile.Run("kernels:engines", CostKernels.EngineBandsBlock);

    [Fact] public void Kernels_Dense() => CostProfile.Run("kernels:dense", CostKernels.DenseBlock);

    [Fact] public void Kernels_Bgp() => CostProfile.Run("kernels:bgp", CostKernels.BgpCholeskyBlock);

    [Fact] public void Kernels_Wide() => CostProfile.Run("kernels:wide", CostKernels.WideBlock);

    [Fact] public void Kernels_Bordered() => CostProfile.Run("kernels:bordered", CostKernels.BorderedBlock);

    [Fact] public void Kernels_LowRank() => CostProfile.Run("kernels:lowrank", CostKernels.LowRankBlock);

    [Fact] public void Kernels_Khatri() => CostProfile.Run("kernels:khatri", CostKernels.KhatriBlock);

    [Fact] public void Engine_Bars() => CostProfile.Run("engine:bars", CostEngines.BarsBlock);

    [Fact] public void Engine_Baps() => CostProfile.Run("engine:baps", CostEngines.BapsBlock);

    [Fact] public void Engine_Bgp() => CostProfile.Run("engine:bgp", CostEngines.BgpBlock);

    [Fact] public void Engine_Changepoint() => CostProfile.Run("engine:changepoint", CostEngines.ChangepointBlock);

    [Fact] public void Engine_GPrior() => CostProfile.Run("engine:gprior", CostEngines.GPriorBlock);

    [Fact] public void Sampler_ReversibleJump() => CostProfile.Run("sampler:rj", CostSamplers.RjBlock);

    [Fact] public void Sampler_Blocks() => CostProfile.Run("sampler:blocks", CostSamplers.BlocksBlock);

    [Fact] public void Sampler_Tempering() => CostProfile.Run("sampler:pt", CostSamplers.TemperingBlock);

    [Fact] public void Sampler_Smc() => CostProfile.Run("sampler:smc", CostSamplers.SmcBlock);

    [Fact] public void Sampler_SmcBass() => CostProfile.Run("sampler:smc-bass", CostSamplers.BassSmcBlock);

    [Fact] public void Sampler_Samc() => CostProfile.Run("sampler:samc", CostSamplers.SamcBlock);

    [Fact] public void Sampler_Wham() => CostProfile.Run("sampler:wham", CostSamplers.WhamBlock);

    [Fact] public void Sampler_Ladder() => CostProfile.Run("sampler:ladder", CostSamplers.LadderBlock);

    [Fact] public void Sampler_Ensemble() => CostProfile.Run("sampler:ensemble", CostSamplers.EnsembleBlock);

    [Fact] public void Alloc_Parallel() => CostProfile.Run("alloc:parallel", CostParallel.AllocationBlock);

    [Fact] public void MemoryBound() => CostProfile.Run("membound", CostParallel.MemoryBoundBlock);

    /// <summary>One process per setting, by environment variable: no discarded pass, so the warm-up reads from the first call.</summary>
    [Fact] public void Settings() => CostProfile.Run("settings", CostParallel.SettingsBlock, discardPass: false);

    [Fact] public void Engine_Levers() => CostProfile.Run("engine:levers", CostLevers.LeversBlock);

    /// <summary>The brief's checks that need no quiet machine, always on: the allocation counter reads a known allocation and
    /// <see cref="CollapsedGPrior.MarginalAt"/> reads nothing, and the counted units equal their closed forms at three sizes.</summary>
    [Fact]
    public void Checks_AllocationCounterAndCountedUnits()
    {
        long before = GC.GetAllocatedBytesForCurrentThread();
        CostKernels.ObjectSink = new byte[1000];
        long known = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.InRange(known, 1000, 1100);

        double[,] g = CostKernels.DenseSpd(20, 3);
        var factor = new DenseCholesky(20);
        factor.Decompose(g);
        CollapsedGPrior block = CollapsedGPrior.FromLeastSquares(factor, CostKernels.Random(20, 4), 5.0, 40.0, 1000, null);
        block.MarginalAt(0.001);   // warm
        before = GC.GetAllocatedBytesForCurrentThread();
        double sink = 0;
        for (int i = 0; i < 1000; i++) sink += block.MarginalAt(0.001).S;
        Assert.Equal(0, GC.GetAllocatedBytesForCurrentThread() - before);
        Assert.True(sink > 0);

        foreach ((int n, int p) in new[] { (64, 3), (200, 8), (500, 32) })
            Assert.InRange(Math.Abs(CostKernels.BandFactorTerms(n, p) - CostKernels.BandFactorTermsClosed(n, p)), 0, (long)p * p * p + (long)p * p);
        foreach (int m in new[] { 5, 50, 500 })
            Assert.Equal((long)m * (m - 1) * (m + 1) / 6 + (long)m * (m - 1) / 2, CostKernels.DenseFactorTerms(m));
        foreach (int n in new[] { 100, 400 })
        {
            long closed = (long)(n - 1) * n * (n + 1) / 6;
            Assert.Equal(closed, CostKernels.CholeskyStep1Terms(n));
            Assert.Equal(closed, CostKernels.CholeskyStep2Terms(n));
        }
    }
}
