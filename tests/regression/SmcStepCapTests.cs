using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Samplers;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;
using static Maths.Regression.Tests.BassTemperingTests;

namespace Maths.Regression.Tests;

/// <summary>
/// SMC's opt-in step cap (<see cref="SmcSpec.StepCap"/>, built 2026-11-09 for KHAOS's trap): each step from b &gt; 0 takes at most
/// r·b, the first step from b = 0 the conditional-ESS rule's alone.
/// <list type="bullet">
/// <item><b>The config:</b> r above 1 and finite, left out of the serialized form when null, the fingerprint moved when set.</item>
/// <item><b>The bound,</b> stage by stage, the steps it bound counted as the readout counts them, and the last step never a jump.</item>
/// <item><b>The evidence</b> against the enumeration on KHAOS's exact-law fixture, the cap with and without the residual birth.</item>
/// <item><b>Resume</b> bit for bit through the cap.</item>
/// <item><b>The trap</b> (qualification): a production fit on KHAOS's raw-ridge trap fixture with the cap and the residual birth.</item>
/// </list>
/// Without the cap the schedule is untouched: SMC's pinned baseline (<see cref="BassBitIdentityTests"/>) reads it.
/// </summary>
public sealed class SmcStepCapTests(ITestOutputHelper output)
{
    private static readonly SmcSpec Small = new() { Particles = 150, Replicates = 2, MutationSteps = 2, TargetConditionalEss = 0.9 };

    private static BassConfig Capped(BassConfig config, SmcSpec smc, double? cap, int seed = 7)
        => config with { Sampling = config.Sampling with { MasterSeed = seed, Smc = smc with { StepCap = cap } } };

    [Fact]
    public void Config_RefusesACapThatCannotReachOne_AndAConfigWithoutItStandsAsItWas()
    {
        BassConfig good = Capped(BassSweepTests.Config(maxBasis: 4), Small, 1.25);
        good.Validate();
        foreach (double bad in new[] { 1.0, 0.9, 0.0, -2.0, double.NaN, double.PositiveInfinity })
            Assert.Throws<ArgumentException>(() => Capped(BassSweepTests.Config(maxBasis: 4), Small, bad).Validate());
        BassConfig plain = Capped(BassSweepTests.Config(maxBasis: 4), Small, null);
        Assert.DoesNotContain("StepCap", plain.ToJson());
        Assert.Contains("\"StepCap\":1.25", good.ToJson());
        Assert.Equal(good, BassConfig.FromJson(good.ToJson()));
        BassData data = BassSweepTests.Data(large: true);
        Assert.NotEqual(Bass.Fingerprint(plain, data), Bass.Fingerprint(good, data));
    }

    [Theory]
    [InlineData(1.25)]
    [InlineData(1.05)]
    public void EveryStepFromAboveZero_StaysWithinTheCap_AndTheReadoutCountsTheBoundSteps(double cap)
    {
        BassData data = BassSweepTests.Data(large: true);
        BassFit plain = Bass.Fit(Capped(BassSweepTests.Config(maxBasis: 4), Small, null), data);
        BassFit capped = Bass.Fit(Capped(BassSweepTests.Config(maxBasis: 4), Small, cap), data);
        Assert.True(capped.Diagnostics.Smc!.ReachedPosterior);
        for (int c = 0; c < capped.Diagnostics.Smc.Replicates.Count; c++)
        {
            BassSmcReplicate x = capped.Diagnostics.Smc.Replicates[c], y = plain.Diagnostics.Smc!.Replicates[c];
            double[] b = [.. x.Stages.Select(s => s.Temperature.InverseTemperature)];
            Assert.Equal(0.0, b[0]);
            Assert.Equal(y.Stages[1].Temperature.InverseTemperature, b[1]);   // the first step is the rule's alone
            int bound = 0;
            for (int k = 2; k < b.Length; k++)
            {
                Assert.True(b[k] > b[k - 1] && b[k] <= Math.Min(1.0, cap * b[k - 1]), $"stage {k}: b {b[k]:R} after {b[k - 1]:R}");
                if (b[k] < 1.0 && b[k] == cap * b[k - 1]) bound++;
            }
            Assert.Equal(bound, x.CappedSteps);
            Assert.Null(y.CappedSteps);
            Assert.True(x.LastStepFrom >= 1.0 / cap, $"the last step left from b {x.LastStepFrom:R}");
            Assert.False(x.JumpedOverACollapse);
            output.WriteLine($"r = {cap}: replicate {c}, {x.StageCount} stages against {y.StageCount} uncapped, {x.CappedSteps} bound; " +
                             $"first step to b {b[1]:G3}, last from {x.LastStepFrom:G3}");
        }
        Assert.Contains(capped.Diagnostics.Smc.Replicates, x => x.CappedSteps > 0);
    }

    [Theory]
    [InlineData("raw ridge", false)]
    [InlineData("raw ridge", true)]
    [InlineData("modified g", true)]
    public void Evidence_MeetsTheEnumeration(string prior, bool residual)
    {
        BassData data = KhaosChainTests.Data();
        BassConfig config = KhaosReachTests.Config(prior);
        if (residual) config = ResidualBirthTests.WithResidual(config, new ResidualBirth(0.25));
        config = config with
        {
            Sampling = config.Sampling with
            {
                MasterSeed = 3, ReportEvidence = true, Smc = new SmcSpec { Particles = 2000, Replicates = 4, MutationSteps = 4, StepCap = 1.25 },
            },
        };
        RungLaw exact = Exact(BassModel.Create(config, data), 1.0)[0];
        BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        BassEvidence evidence = fit.Evidence!;
        BassSmcReadout readout = fit.Diagnostics.Smc!;
        output.WriteLine($"{prior}{(residual ? ", the residual birth" : "")}, r = 1.25: ln p(y) {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} " +
                         $"against {exact.LogZ:F4}; stages {string.Join(", ", readout.Replicates.Select(x => $"{x.StageCount} ({x.CappedSteps} bound)"))}");
        Assert.True(Math.Abs(evidence.LogEvidence - exact.LogZ) <= 4 * evidence.MonteCarloError + 0.01,
                    $"{evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against {exact.LogZ:F4}");
        Assert.Empty(readout.Flags);
        Assert.All(readout.Replicates, x => Assert.True(x.CappedSteps > 0));
    }

    [Fact]
    public void Resume_IsBitForBit_ThroughTheCap()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig Config(int maxStages) => Capped(BassSweepTests.Config(maxBasis: 4), Small with { MaxStages = maxStages }, 1.1, seed: 41);
        BassFit whole = Bass.Fit(Config(5000), data);
        var run = new BassRunOptions { Checkpoint = true };
        BassFit resumed = Bass.Fit(Config(3), data, run);
        int resumes = 0;
        while (!resumed.Diagnostics.Smc!.ReachedPosterior)
        {
            resumed = Bass.Resume(resumed.Checkpoint!, Config(3), data, run);
            Assert.True(++resumes <= 500, "the resumes never reached b = 1");
        }
        Assert.True(resumes > 1);
        Assert.Equal(BassBitIdentityTests.Hash(whole), BassBitIdentityTests.Hash(resumed));
        Assert.Equal(whole.Diagnostics.Smc!.Replicates.Select(x => x.CappedSteps), resumed.Diagnostics.Smc.Replicates.Select(x => x.CappedSteps));
        Assert.Equal(whole.Evidence!.LogEvidence, resumed.Evidence!.LogEvidence);
        output.WriteLine($"{resumes} resumes; bound steps {string.Join(", ", whole.Diagnostics.Smc.Replicates.Select(x => x.CappedSteps))}");
    }

    // A production fit on KHAOS's raw-ridge trap fixture (BassSmcMutationQualificationTests.Khaos: Ishigami, n 120, p 4, d_max 8),
    // N 500 at 5 sweeps a stage, four replicates: with the cap and the residual birth each reaches the posterior, its weighted σ²
    // within twice the long chains' 0.2245, and none jumps.
    [Fact]
    public void Trap_TheCapAndTheResidualBirth_ReachThePosterior()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassSmcMutationQualificationTests.Fixture f = BassSmcMutationQualificationTests.Khaos();
        BassConfig config = ResidualBirthTests.WithResidual(f.Config, new ResidualBirth(0.25)) is var c
            ? c with { Sampling = c.Sampling with { MasterSeed = 20261109, ReportEvidence = true, Smc = new SmcSpec { Particles = 500, Replicates = 4, MutationSteps = 5, StepCap = 1.25 } } }
            : throw new InvalidOperationException();
        BassFit fit = Bass.Fit(config, f.Data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        BassSmcReadout readout = fit.Diagnostics.Smc!;
        double[] weights = fit.DrawWeights!;
        int i = 0;
        for (int r = 0; r < fit.Draws.Count; r++)
        {
            double sigma2 = 0.0, m = 0.0, total = 0.0;
            foreach (BassDraw d in fit.Draws[r])
            {
                sigma2 += weights[i] * d.Sigma2;
                m += weights[i] * d.Structure.TermIds.Length;
                total += weights[i++];
            }
            BassSmcReplicate x = readout.Replicates[r];
            output.WriteLine($"replicate {r}: σ² {sigma2 / total:F4}, M {m / total:F2}, ln Ẑ {x.LogNormalizer:F2} ± {x.LogNormalizerError:F2}, " +
                             $"{x.StageCount} stages ({x.CappedSteps} bound), last from b {x.LastStepFrom:G3}");
            Assert.True(sigma2 / total < 2 * 0.2245, $"replicate {r} is trapped: σ² {sigma2 / total:F3}");
            Assert.False(x.JumpedOverACollapse);
        }
        output.WriteLine($"flags: {string.Join(" | ", readout.Flags)}");
    }
}
