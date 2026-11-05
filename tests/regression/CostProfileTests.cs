using System;
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
}
