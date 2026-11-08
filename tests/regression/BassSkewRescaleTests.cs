using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The skew's rescaling move (<see cref="BassSkewRescaleBlock"/>; slice 5's and slice 7's follow-up), β → β/κ and v → κ·v along the ridge
/// β·vᵢ holds:
/// <list type="bullet">
/// <item><b>Its statement on a frozen state:</b> from one state every move stays on the orbit (β₀·e⁻ᵗ, eᵗ·v₀), so its stationary law there is
/// the joint along the orbit with the map's Jacobian, ∝ π(Tₜ(x₀))·e^((N−1)t), read on a grid by the engine's own target; the move's chain
/// meets its moments. Its acceptance from x to x′ and back cancel.</item>
/// <item><b>The switch:</b> it needs a drawn skew and stays out of the serialized form when off.</item>
/// <item><b>The run:</b> a tempered fit with it resumes bit for bit.</item>
/// </list>
/// The joint law with the whole sweep, the measurement and the cost are behind <c>BASS_QUALIFICATION</c>.
/// </summary>
public sealed class BassSkewRescaleTests(ITestOutputHelper output)
{
    private static readonly BassNoise SkewedHorseshoe = BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5));
    private static readonly BassNoise SkewedGh = BassNoise.GeneralizedHyperbolic(new GigParameters(0.7, 1.2, 0.8), new BassNoisePrior(0.3, 0.5));

    private static BassConfig Config(BassNoise noise) => BassStudentTTests.Config(noise) with
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 6 },
        Sigma2 = Sigma2Prior.InverseGamma(3.0, 0.5),
    };

    [Fact]
    public void TheSwitch_NeedsADrawnSkew_AndIsSerializedOnlyWhenOn()
    {
        BassConfig plain = Config(SkewedHorseshoe), moving = Config(SkewedHorseshoe with { SkewRescale = true });
        Assert.DoesNotContain("SkewRescale", plain.ToJson());
        Assert.Contains("\"SkewRescale\":true", moving.ToJson());
        Assert.Equal(moving, BassConfig.FromJson(moving.ToJson()));
        BassData data = BassQuantileTests.Data();
        Assert.NotEqual(Bass.Fingerprint(plain, data), Bass.Fingerprint(moving, data));
        Assert.Throws<ArgumentException>(() => BassNoise.StudentT(4) with { SkewRescale = true });
        Assert.Throws<ArgumentException>(() => BassNoise.AsymmetricLaplace(0.9) with { SkewRescale = true });
        Assert.Throws<ArgumentException>(() => BassNoise.Horseshoe(BassNoisePrior.Held(0.5)) with { SkewRescale = true });

        BassFit fit = Bass.Fit(moving, data);
        MoveStat move = fit.Diagnostics.Chains[0].Blocks!.Single(k => k.Key == "skew rescale");
        output.WriteLine($"skew rescale: {move.Accepted} of {move.Attempts} accepted");
        Assert.True(move.Attempts > 0 && move.Accepted > 0);
    }

    public static IEnumerable<object[]> Laws() => [["horseshoe", 1.0], ["gh", 1.0], ["horseshoe", 0.5], ["gh", 0.5]];

    [Theory]
    [MemberData(nameof(Laws))]
    public void TheMove_KeepsItsOrbitsLaw_OnAFrozenState(string law, double b)
    {
        BassNoise noise = (law == "gh" ? SkewedGh : SkewedHorseshoe) with { SkewRescale = true };
        BassData data = BassQuantileTests.Data();
        BassModel model = BassModel.Create(Config(noise), data);
        var mixture = (VarianceMeanMixture)model.Observation;
        var rng = new Xoshiro256PlusPlus(20261108);
        // A structure and latents to freeze at: a short chain's.
        var start = new BassChain(model, [], 0.0, new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        for (int s = 0; s < 300; s++) start.Step();
        var moves = new BassStructureMoves(model, start.Current.Terms, start.Current.LogTau);
        var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(3),
                                                       new Temperature(b, TemperingPath.Likelihood));
        Assert.Null(BassStructureMoves.Follow(chain, start.Current.Fit.Working!));
        moves.DrawSweepFit(chain.Current, rng, model.Target, b);
        var block = new BassSkewRescaleBlock(mixture);
        BassSkewRescaleBlock.Proposal x0 = block.Read(chain);
        double h = x0.Scale, beta0 = x0.Skew, sd = noise.SkewPrior!.Sd, m = noise.SkewPrior.Mean, current = chain.CurrentLogDensity;
        int n = x0.V.Length;

        // The orbit's law, ln p(t) = ln π(Tₜ(x₀)) + (N − 1)·t, π the target with v's and β's priors, on a grid.
        double LogOrbit(double t)
        {
            double[] v = [.. x0.V.Select(vi => vi * Math.Exp(t))];
            double skew = beta0 * Math.Exp(-t);
            BassState candidate = moves.ProposeFollow(chain.Current, mixture.At(x0.Working, new ObservationLatents(v, skew, x0.Sigma2, x0.Latents.Gamma)));
            double target = chain.Temperature.LogDensity(chain.Target.Evaluate(candidate)) - current;
            return target + v.Sum(x0.Mixing.LogPrior) - x0.V.Sum(x0.Mixing.LogPrior) - (skew - m) * (skew - m) / (2.0 * sd * sd) + (n - 1) * t;
        }
        const int points = 801;
        double width = 10.0 * h;
        var grid = new double[points];
        var logs = new double[points];
        for (int k = 0; k < points; k++) logs[k] = LogOrbit(grid[k] = -width + 2.0 * width * k / (points - 1));
        double top = logs.Max(), mass = 0.0, first = 0.0, second = 0.0;
        for (int k = 0; k < points; k++)
        {
            double w = Math.Exp(logs[k] - top);
            mass += w;
            first += w * grid[k];
            second += w * grid[k] * grid[k];
        }
        double mean = first / mass, meanSquare = second / mass;
        Assert.True(Math.Exp(logs[0] - top) < 1e-10 && Math.Exp(logs[^1] - top) < 1e-10, "the grid holds the orbit's mass");

        // The move's chain from x₀: t = ln(β₀/β) at each step.
        const int steps = 60_000, batches = 40;
        var ts = new double[steps];
        int accepted = 0;
        for (int k = 0; k < steps; k++)
        {
            if (block.Update(chain, rng)) accepted++;
            ts[k] = Math.Log(beta0 / chain.Current.Fit.Working!.Latents!.Skew);
        }
        double Z(Func<double, double> g, double expected)
        {
            int size = steps / batches;
            double[] means = [.. Enumerable.Range(0, batches).Select(b => ts.Skip(b * size).Take(size).Average(g))];
            double average = means.Average(), se = Math.Sqrt(means.Sum(x => (x - average) * (x - average)) / (batches - 1) / batches);
            return (average - expected) / se;
        }
        double zMean = Z(t => t, mean), zSquare = Z(t => t * t, meanSquare);
        output.WriteLine($"{law} at b = {b}: h {h:F4}, the orbit's sd {Math.Sqrt(meanSquare - mean * mean):F4}; accepted {accepted} of {steps}; " +
                         $"mean t {ts.Average():F5} against {mean:F5} (z {zMean:F2}), mean t² z {zSquare:F2}");
        Assert.True(Math.Abs(zMean) < 4.0 && Math.Abs(zSquare) < 4.0, $"z {zMean:F2}, {zSquare:F2}");
        Assert.InRange(accepted / (double)steps, 0.15, 0.85);

        // From x to x′ and back the acceptances cancel: the reverse's proposal, its Jacobian and the priors are each other's.
        BassSkewRescaleBlock.Proposal x = block.Read(chain);
        double step = 0.7 * x.Scale;
        BassSkewRescaleBlock.Move forward = block.Propose(chain, x, step)!.Value;
        chain.SetState(forward.Candidate, forward.Evaluation);
        BassSkewRescaleBlock.Move back = block.Propose(chain, block.Read(chain), -step)!.Value;
        output.WriteLine($"{law}: ln A forward {forward.LogAcceptance:R}, back {back.LogAcceptance:R}");
        Assert.Equal(0.0, forward.LogAcceptance + back.LogAcceptance, 8);
    }

    [Fact]
    public void ATemperedFitWithTheMove_ResumesBitForBit()
    {
        BassData data = BassQuantileTests.Data();
        BassConfig whole = Config(SkewedHorseshoe with { SkewRescale = true }) with
        {
            Sampling = new BassSampling { Chains = 2, MasterSeed = 31, BurnIn = 100, Samples = 120, Thin = 2,
                                          Tempering = new TemperingSpec { Levels = 4, HotEnd = 0.1, ReferenceRung = true }, ReportEvidence = true },
        };
        BassConfig half = whole with { Sampling = whole.Sampling with { Samples = 60 } };
        BassFit full = Bass.Fit(whole, data);
        BassFit first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        BassFit resumed = Bass.Resume(first.Checkpoint!, half, data);
        BassRunTests.AssertSameFit(full, resumed);
        for (int c = 0; c < full.Draws.Count; c++)
            for (int d = 0; d < full.Draws[c].Count; d++)
            {
                ObservationLatents a = full.Draws[c][d].Latents!, b = resumed.Draws[c][d].Latents!;
                Assert.True(a.Values.SequenceEqual(b.Values));
                Assert.Equal(a.Skew, b.Skew);
                Assert.Equal(a.Sigma2, b.Sigma2);
            }
    }
}
