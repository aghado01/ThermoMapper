using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's SMC particles and kernel (stage 8, C4's step 2): a fork is an exact, independent copy of its parent, the
/// parent's design and later moves untouched by the fork's; and a tempered run gives the same population, normalizer,
/// stage reports and account whether its particles move one at a time or in parallel. Each stage's law against the
/// enumeration oracle is <c>BassSmcParityTests</c>'.
/// </summary>
public sealed class BassSmcTests
{
    private readonly ITestOutputHelper _out;
    public BassSmcTests(ITestOutputHelper output) => _out = output;

    private static BassModel Model() => BassModel.Create(BassSweepTests.Config(maxBasis: 4), BassSweepTests.Data(large: true));

    private static SmcStage<BassParticle> Stage(BassModel model, double b) => new(new BassParticleTarget(model.Target), new Temperature(b, TemperingPath.Likelihood));

    // A particle's design as it stands: its terms, its statistics and its coefficients.
    private static (BassTerm[] Terms, double Residual, double Fitted, double[] Beta, double LogTau) Design(BassParticle p)
    {
        BassDesign d = p.Moves.DesignAt(p.State);
        return ([.. d.Terms], d.ResidualSquares, d.FittedSquares, d.Coefficients.ToArray(), p.State.LogTau);
    }

    private static List<(BassTerm[] Terms, double LogTau)> Walk(BassSmcKernel kernel, BassParticle p, SmcStage<BassParticle> stage, int seed, int sweeps)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var path = new List<(BassTerm[], double)>();
        ITargetEvaluation e = stage.Target.Evaluate(p);
        for (int s = 0; s < sweeps; s++)
        {
            e = kernel.Move(p, e, stage, rng).Evaluation;
            path.Add(([.. p.State.Terms], p.State.LogTau));
        }
        return path;
    }

    private static void SamePath(List<(BassTerm[] Terms, double LogTau)> a, List<(BassTerm[] Terms, double LogTau)> b)
    {
        Assert.Equal(a.Count, b.Count);
        for (int s = 0; s < a.Count; s++)
        {
            Assert.Equal(a[s].Terms, b[s].Terms);
            Assert.Equal(a[s].LogTau, b[s].LogTau);
        }
    }

    [Fact]
    public void Fork_IsAnExactCopy_AndIndependentOfItsParent()
    {
        BassModel model = Model();
        var particles = new BassParticles(model);
        var kernel = new BassSmcKernel(model);
        SmcStage<BassParticle> stage = Stage(model, 0.6);
        BassParticle Parent()
        {
            BassParticle p = particles.Draw(new Xoshiro256PlusPlus(41));
            Walk(kernel, p, stage, seed: 42, sweeps: 60);
            return p;
        }

        // A fork moves exactly as its parent would from the same stream: the design is copied bit for bit.
        BassParticle parent = Parent(), twin = Parent();
        BassParticle fork = particles.Fork(parent);
        Assert.NotSame(parent.State, fork.State);
        Assert.NotSame(parent.Moves, fork.Moves);
        Assert.False(ReferenceEquals(parent.State.Terms, fork.State.Terms));
        var before = Design(parent);
        List<(BassTerm[] Terms, double LogTau)> forked = Walk(kernel, fork, stage, seed: 43, sweeps: 200);
        Assert.True(forked.Select(s => s.Terms.Length).Distinct().Count() > 1, "the fork's walk never changed its structure");

        // The fork's walk left the parent as it was, and the parent then moves as its twin, which was never forked.
        var after = Design(parent);
        Assert.Equal(before.Terms, after.Terms);
        Assert.Equal(before.Residual, after.Residual);
        Assert.Equal(before.Fitted, after.Fitted);
        Assert.Equal(before.Beta, after.Beta);
        Assert.Equal(before.LogTau, after.LogTau);
        List<(BassTerm[] Terms, double LogTau)> parentPath = Walk(kernel, parent, stage, seed: 43, sweeps: 200);
        SamePath(parentPath, forked);
        SamePath(Walk(kernel, twin, stage, seed: 43, sweeps: 200), parentPath);
    }

    [Fact]
    public void Run_InParallel_MatchesTheSerialRun()
    {
        BassModel model = Model();
        (SmcResult<BassParticle> Result, BassSmcAccount Account, long Attempts) Run(int parallelism)
        {
            var particles = new BassParticles(model);
            var kernel = new BassSmcKernel(model);
            var options = new SmcOptions { Particles = 96, MutationSteps = 2, MaxDegreeOfParallelism = parallelism };
            SmcResult<BassParticle> result = SmcSampler.RunTempered(new BassParticleTarget(model.Target), TemperingPath.Likelihood, 0.9, particles, kernel,
                                                                     options, seed: 2026);
            return (result, kernel.Account, particles.PriorAttempts);
        }
        var serial = Run(1);
        var parallel = Run(6);
        _out.WriteLine($"{serial.Result.Stages.Count - 1} stages, ln Ẑ {serial.Result.LogNormalizerRatio:F4} ± {serial.Result.LogNormalizerError:F4}, " +
                       $"acceptance {serial.Account.AcceptanceRate:F3}, prior attempts {serial.Attempts}");
        Assert.Equal(1.0, serial.Result.Stages[^1].Temperature.InverseTemperature);
        Assert.Equal(serial.Result.LogNormalizerRatio, parallel.Result.LogNormalizerRatio);
        Assert.Equal(serial.Result.Weights, parallel.Result.Weights);
        Assert.Equal(serial.Result.Stages, parallel.Result.Stages);
        for (int i = 0; i < serial.Result.Particles.Count; i++)
        {
            Assert.Equal(serial.Result.Particles[i].State.Terms, parallel.Result.Particles[i].State.Terms);
            Assert.Equal(serial.Result.Particles[i].State.LogTau, parallel.Result.Particles[i].State.LogTau);
        }
        Assert.Equal(serial.Attempts, parallel.Attempts);
        Assert.Equal(serial.Account.Moves, parallel.Account.Moves);
        Assert.Equal(serial.Account.Kinds.Entries, parallel.Account.Kinds.Entries);
        Assert.Equal(serial.Account.Refusals.Entries, parallel.Account.Refusals.Entries);
        Assert.Equal(serial.Account.Support.Entries, parallel.Account.Support.Entries);
        Assert.Equal(serial.Account.EmptyBirths, parallel.Account.EmptyBirths);
        // Every sweep is counted: the stages' moves over the particles of positive weight.
        Assert.Equal((serial.Result.Stages.Count - 1) * 96L * 2, serial.Account.Sweeps);
        Assert.Equal(serial.Account.Moves.Sum(m => m.Attempts), serial.Account.Kinds.Entries.Sum(k => k.Attempts));
    }

    [Fact]
    public void Kernel_RefusesAnotherModelsStage()
    {
        BassModel model = Model(), other = Model();
        var kernel = new BassSmcKernel(model);
        BassParticle p = new BassParticles(model).Draw(new Xoshiro256PlusPlus(3));
        SmcStage<BassParticle> foreign = Stage(other, 0.5);
        Assert.Throws<ArgumentException>(() => kernel.Move(p, foreign.Target.Evaluate(p), foreign, new Xoshiro256PlusPlus(4)));
    }
}
