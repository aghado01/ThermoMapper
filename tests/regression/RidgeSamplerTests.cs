using System;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The sampler a ridge basis takes when none is declared (the user's decision of 2026-11-15, E36): self-sizing parallel tempering, with the
/// exact reference rung where the priors and the tempered form admit it and an auto-picked hot end with its reason otherwise; declared
/// samplers win, declared plain chains carry a caution, the resolution is logged, a checkpoint resumes with the sampler that produced it,
/// and every other fit, and every existing config's form and fingerprint, stand.
/// </summary>
public sealed class RidgeSamplerTests(ITestOutputHelper output)
{
    // A ridge basis under BASS's default priors, which admit the reference rung: the mean centre, σ² scaled to the response, unit-information τ.
    private static BassConfig Ridge(int columns = 4, int samples = 40) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 6, Directional = new DirectionalBasis { Columns = columns } },
        Sampling = new BassSampling { Chains = 2, MasterSeed = 3761, BurnIn = 300, Samples = samples },
    };

    private static BassAutoPick? SamplerRow(BassModel model) => model.Resolved.AutoPicks.FirstOrDefault(p => p.Setting == "Sampler");

    [Fact]
    public void AnOmittedSampler_OnARidgeBasis_IsSelfSizingTempering_WithTheReferenceRung_Logged_AndTheFitKeepsTheDeclaredConfig()
    {
        BassData data = DirectionalChainTests.Data();
        BassConfig config = Ridge();
        BassModel model = BassModel.Create(config, data);
        TemperingSpec tempering = model.Config.Sampling.Tempering!;
        Assert.True(tempering.ReferenceRung);
        Assert.Null(tempering.Levels);   // self-sizing
        Assert.Same(config, model.Declared);
        BassAutoPick row = SamplerRow(model)!;
        Assert.Equal("parallel tempering, self-sizing, exact reference rung", row.Value);
        Assert.Same(row, model.Resolved.AutoPicks[0]);
        output.WriteLine($"{row.Setting}: {row.Value}; {row.Reason}");

        BassFit fit = Bass.Fit(config, data);
        Assert.Equal(config, fit.Config);   // the declared config, whose fingerprint the identity reads
        Assert.Equal(Bass.Fingerprint(config, data), fit.Identity.Fingerprint);
        Assert.NotNull(fit.Diagnostics.Tempering);
        Assert.True(fit.Diagnostics.Tempering!.ReferenceRung);
        Assert.Null(fit.Diagnostics.SamplerCaution);
        Assert.Contains(fit.Resolved.AutoPicks, p => p.Setting == "Sampler");
    }

    [Fact]
    public void WhereTheReferenceRungIsRefused_AHotEndIsAutoPicked_AndTheRefusalIsTheReason()
    {
        // Jeffreys' σ² leaves the integrated form undefined at b = 0, so no reference rung; a positive hot end is admitted.
        BassData data = DirectionalChainTests.Data();
        BassConfig config = Ridge() with { Sigma2 = Sigma2Prior.Jeffreys };
        BassModel model = BassModel.Create(config, data);
        TemperingSpec tempering = model.Config.Sampling.Tempering!;
        Assert.False(tempering.ReferenceRung);
        BassAutoPick row = SamplerRow(model)!;
        Assert.Equal("parallel tempering, self-sizing, auto-picked hot end", row.Value);
        Assert.Contains("no reference rung, since", row.Reason);
        Assert.Contains("ReferenceRung", row.Reason);
        Assert.DoesNotContain("(Parameter", row.Reason);
        Assert.Contains(model.Resolved.AutoPicks, p => p.Setting == "Tempering.HotEnd");
        output.WriteLine($"{row.Setting}: {row.Value}; {row.Reason}");
        BassFit fit = Bass.Fit(config, data);
        Assert.False(fit.Diagnostics.Tempering!.ReferenceRung);

        // The ladder takes the form the model runs: a caller naming the linear form, which reaches b = 0 under Jeffreys' σ², gets the rung.
        TemperingSpec linear = BassModel.Create(config, data, TemperedForm.Linear).Config.Sampling.Tempering!;
        Assert.Equal(TemperedForm.Linear, linear.Form);
        Assert.True(linear.ReferenceRung);
    }

    [Fact]
    public void DeclaredPlainChains_RunWithTheCaution_AndDeclaredSamplersWin()
    {
        BassData data = DirectionalChainTests.Data();
        BassConfig plain = Ridge() with { Sampling = Ridge().Sampling with { Sampler = BassSamplerChoice.PlainChains } };
        BassModel model = BassModel.Create(plain, data);
        Assert.Null(model.Config.Sampling.Tempering);
        Assert.Equal("plain chains", SamplerRow(model)!.Value);
        BassFit fit = Bass.Fit(plain, data);
        Assert.Null(fit.Diagnostics.Tempering);
        Assert.Contains("can hold M where they started", fit.Diagnostics.SamplerCaution);
        Assert.Contains("too narrow", fit.Diagnostics.SamplerCaution);
        output.WriteLine($"caution: {fit.Diagnostics.SamplerCaution}");

        // A declared SMC run or ladder is run as declared, with nothing resolved or logged.
        BassConfig smc = Ridge() with { Sampling = Ridge().Sampling with { Chains = 1, Smc = new SmcSpec { Particles = 50, Replicates = 1, MutationSteps = 2 } } };
        BassModel smcModel = BassModel.Create(smc, data);
        Assert.Same(smc, smcModel.Config);
        Assert.Null(SamplerRow(smcModel));
        var ladder = new TemperingSpec { Levels = 3, HotEnd = 0.3 };
        BassConfig tempered = Ridge() with { Sampling = Ridge().Sampling with { Tempering = ladder } };
        BassModel temperedModel = BassModel.Create(tempered, data);
        Assert.Same(ladder, temperedModel.Config.Sampling.Tempering);
        Assert.Null(SamplerRow(temperedModel));

        // Plain chains beside a declared sampler are refused.
        Assert.Throws<ArgumentException>(() => (tempered with { Sampling = tempered.Sampling with { Sampler = BassSamplerChoice.PlainChains } }).Validate());
        Assert.Throws<ArgumentException>(() => (plain with { Sampling = plain.Sampling with { Sampler = (BassSamplerChoice)7 } }).Validate());
    }

    [Fact]
    public void EveryOtherBasis_PassesThroughUntouched_AndTheConfigsFormStands()
    {
        BassData data = DirectionalChainTests.Data();
        // The one-column directional hinge, slice 2's basis, and BASS's hinges: plain chains, the config itself, nothing logged, no caution.
        foreach (BassConfig config in new[] { Ridge(columns: 1), Ridge() with { Structure = Ridge().Structure with { Directional = null } } })
        {
            BassModel model = BassModel.Create(config, data);
            Assert.Same(config, model.Config);
            Assert.Null(SamplerRow(model));
            Assert.Null(model.SamplerCaution);
            BassFit fit = Bass.Fit(config, data);
            Assert.Null(fit.Diagnostics.Tempering);
            Assert.Null(fit.Diagnostics.SamplerCaution);
        }
        // The choice is left out of the serialized form at its default, so no existing config's form or fingerprint moves; declared, it
        // round-trips.
        Assert.DoesNotContain("Sampler", Ridge().ToJson());
        BassConfig plain = Ridge() with { Sampling = Ridge().Sampling with { Sampler = BassSamplerChoice.PlainChains } };
        Assert.Contains("Sampler", plain.ToJson());
        Assert.Equal(BassSamplerChoice.PlainChains, BassConfig.FromJson(plain.ToJson()).Sampling.Sampler);
        Assert.NotEqual(Bass.Fingerprint(Ridge(), data), Bass.Fingerprint(plain, data));
    }

    [Fact]
    public void ACheckpoint_ResumesWithTheSamplerThatProducedIt()
    {
        BassData data = DirectionalChainTests.Data();

        // The automatic ladder: a half run's checkpoint resumed is the whole run, tempered, bit for bit.
        BassFit whole = Bass.Fit(Ridge(samples: 40), data);
        BassFit half = Bass.Fit(Ridge(samples: 20), data, new BassRunOptions { Checkpoint = true });
        Assert.True(BassEnsemble.Tempered(half.Checkpoint!));
        BassFit resumed = Bass.Resume(half.Checkpoint!, Ridge(samples: 20), data);
        Assert.NotNull(resumed.Diagnostics.Tempering);
        BassRunTests.AssertSameFit(whole, resumed);

        // A checkpoint whose chains ran plain under the same declared config (as one taken before the default changed would hold) resumes
        // plain, with the caution, and goes on as the plain run's own resume does.
        BassConfig plain = Ridge(samples: 20) with { Sampling = Ridge(samples: 20).Sampling with { Sampler = BassSamplerChoice.PlainChains } };
        EnsembleCheckpoint ran = Bass.Fit(plain, data, new BassRunOptions { Checkpoint = true }).Checkpoint!;
        Assert.False(BassEnsemble.Tempered(ran));
        EnsembleCheckpoint legacy = ran with { Identity = ran.Identity with { Fingerprint = Bass.Fingerprint(Ridge(samples: 20), data) } };
        BassFit continued = Bass.Resume(legacy, Ridge(samples: 20), data);
        Assert.Null(continued.Diagnostics.Tempering);
        Assert.NotNull(continued.Diagnostics.SamplerCaution);
        Assert.Contains(continued.Resolved.AutoPicks, p => p.Setting == "Sampler" && p.Reason.StartsWith("as the resumed checkpoint ran"));
        BassFit reference = Bass.Resume(ran, plain, data);
        Assert.Equal(reference.Draws.Count, continued.Draws.Count);
        for (int c = 0; c < reference.Draws.Count; c++)
            for (int i = 0; i < reference.Draws[c].Count; i++)
            {
                Assert.Equal(reference.Draws[c][i].Coefficients, continued.Draws[c][i].Coefficients);
                Assert.Equal(reference.Draws[c][i].Structure.TermIds, continued.Draws[c][i].Structure.TermIds);
            }
        output.WriteLine($"tempered resume: {whole.Draws.Sum(c => c.Count)} draws bit for bit; plain resume: {continued.Draws.Sum(c => c.Count)} draws as the plain run's own");
    }
}
