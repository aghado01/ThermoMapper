using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The interweaving step's worth (E1's slice 7), behind <c>BASS_QUALIFICATION</c> = 1: R̂ and the bulk and tail effective sizes of ln σ², β,
/// mean ln v and the term count M, per second of the chains' work, with and without <see cref="BassNoise.Interweave"/>, on one data set
/// under t, the asymmetric Laplace, a skewed generalized hyperbolic and a skewed horseshoe. Its law is the reach harness's Geweke cases
/// (<see cref="BassNoiseReachQualificationTests"/>).
/// </summary>
public sealed class BassInterweavingQualificationTests(ITestOutputHelper output)
{
    private static readonly string[] Names = ["ln sigma2", "beta", "mean ln v", "M"];

    [Theory]
    [InlineData(11)]
    [InlineData(23)]
    public void EffectiveSizesPerUnitOfWork_WithAndWithoutInterweaving(int masterSeed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int n = 300;
        var rng = new Xoshiro256PlusPlus(20261108);
        double[] x0 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())], x1 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        // Skewed noise outside every family here: 0.3·(E − 1), E ~ Exp(1).
        double[] y = [.. Enumerable.Range(0, n).Select(i => Math.Sin(2.0 * Math.PI * x0[i]) + x1[i] + 0.3 * (-Math.Log(1.0 - rng.NextDouble()) - 1.0))];
        var data = new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
        var noises = new (string Name, BassNoise Noise)[]
        {
            ("t(4)", BassNoise.StudentT(4)),
            ("AL(0.9)", BassNoise.AsymmetricLaplace(0.9)),
            ("skewed GH", BassNoise.GeneralizedHyperbolic(new GigParameters(0.7, 1.2, 0.8), new BassNoisePrior(0.0, 1.0))),
            // Slice 7's reading, the sweep without the move along β·v, which the skewed horseshoe runs by default since.
            ("skewed horseshoe", BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5)) with { SkewRescale = false }),
        };
        output.WriteLine($"Seed {masterSeed}. Work: the fit's thread cycles on one thread (MaxDegreeOfParallelism 1), indicative under load; 4 chains × 8,000 draws after 2,000, at most 6 additive terms.");
        foreach ((string name, BassNoise noise) in noises)
        {
            var perSecond = new double[2][];
            for (int w = 0; w < 2; w++)
            {
                bool interweave = w == 1;
                BassConfig config = BassStudentTTests.Config(noise with { Interweave = interweave }) with
                {
                    Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 6 },
                    Sampling = new BassSampling { Chains = 4, MasterSeed = masterSeed, BurnIn = 2000, Samples = 8000 },
                };
                var wall = Stopwatch.StartNew();
                long t0 = BassCost.Tick();
                BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
                double seconds = (BassCost.Tick() - t0) / BassCost.TicksPerSecond;
                wall.Stop();
                var cells = new List<string>();
                perSecond[w] = new double[Names.Length];
                for (int g = 0; g < Names.Length; g++)
                {
                    double[][] chains = [.. fit.Draws.Select(chain => chain.Select(d => Functional(d, g)).ToArray())];
                    ConvergenceDiagnostic diagnostic = ChainDiagnostics.Diagnose(chains);
                    double bulk = diagnostic.BulkEss.Value;
                    perSecond[w][g] = bulk / seconds;
                    cells.Add($"{Names[g]} R̂ {diagnostic.RHat.Value:F3} ESS {bulk:F0}/{diagnostic.TailEss.Value:F0} ({bulk / seconds:F1}/s)");
                }
                MoveStat[] block = [.. fit.Diagnostics.Chains[0].Blocks!.Where(k => k.Key == "interweaving")];
                output.WriteLine($"{name}, {(interweave ? "interwoven" : "plain")}: {seconds:F1} s ({wall.Elapsed.TotalSeconds:F1} s wall)" +
                                 (block.Length == 0 ? "" : $", interweaving accepted {block[0].Accepted} of {block[0].Attempts}") + $"; {string.Join("; ", cells)}");
            }
            output.WriteLine($"{name}: bulk ESS per second, interwoven over plain: " +
                             string.Join(", ", Names.Select((label, g) => $"{label} {perSecond[1][g] / perSecond[0][g]:F2}")));
        }
    }

    private const int N = 12;
    private const double Alpha0 = 8.0, Beta0 = 6.0;
    private static readonly double[] X = [.. Enumerable.Range(0, N).Select(i => i / (N - 1.0))];
    private static readonly string[] AloneNames = ["M", "ln tau", "ln sigma2", "ln v0", "mean ln v", "atan f0", "atan f6", "z0^2"];

    private static BassData Data(double[] y) => new([new ContinuousPredictor("x", X)], new ScalarResponse(y));

    private static BassNoise Noise(string law) => law switch
    {
        "gh" => BassNoise.GeneralizedHyperbolic(new GigParameters(0.7, 1.2, 0.8), BassNoisePrior.Held(0.5), reference: BassNoiseReference.None),
        _ => BassNoise.Horseshoe(BassNoisePrior.Held(0.5), reference: BassNoiseReference.None),
    };

    private static BassConfig Config(string law, bool interweave) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 2, Support = new SupportRule(SupportKind.Count, 2) },
        Coefficients = BassCoefficientPrior.G(TauPrior.UnitInformationGamma(5)),
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.InverseGamma(Alpha0, Beta0),
        Noise = Noise(law) with { Interweave = interweave },
        Sampling = new BassSampling { Chains = 1, MasterSeed = 1, BurnIn = 0, Samples = 1 },
    };

    private static WorkingResponse At(BassModel model, double[] v, double sigma2, double beta)
        => ((VarianceMeanMixture)model.Observation).At(model.Working, new ObservationLatents(v, beta, sigma2));

    // The step alone (E1's slice 7): β held at 0.5, so σ moves only by the interweaving step, the σ² block left out of the sweep; Geweke's
    // test at b = 1 on the reach harness's design. The same design with σ moved only by its own block reads the skewed horseshoe's v
    // mixing slowly (slice 7's report).
    [Theory]
    [InlineData("horseshoe", 20261107)]
    [InlineData("gh", 20261107)]
    public void Geweke_TheStepAlone_LeavesTheJointInvariant(string law, int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int draws = 400_000, steps = 600_000, batches = 50;
        const bool interweave = true;
        BassConfig config = Config(law, interweave);
        BassNoise noise = config.Noise!;
        double c = noise.VarianceScale, beta = noise.Skew;
        var sigmaPrior = InverseGamma.ShapeScale(Alpha0, Beta0);
        IMixingLaw mixing = noise.Mixing;
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(config, Data([.. X.Select(x => Math.Sin(7 * x))]));
        var sampler = new BassPriorSampler(prior);

        double[] Functions(BassState s, double sigma2, double[] v, double[] f, double[] y)
            => [s.Count, s.LogTau, Math.Log(sigma2), Math.Log(v[0]), v.Average(x => Math.Log(x)), Math.Atan(f[0]), Math.Atan(f[N / 2]),
                Math.Pow(y[0] - f[0] - Math.Sqrt(sigma2) * beta * v[0], 2) / (sigma2 * c * v[0])];
        double[] Respond(double[] f, double sigma2, double[] v)
            => [.. Enumerable.Range(0, N).Select(i => f[i] + Math.Sqrt(sigma2) * (beta * v[i] + Math.Sqrt(c * v[i]) * Normal.SampleStandard(rng)))];

        (BassState S, double Sigma2, double[] V, double[] F) Draw()
        {
            while (true)
            {
                BassState state = sampler.DrawState(rng).State;
                double[] v = [.. Enumerable.Range(0, N).Select(_ => mixing.SamplePrior(rng))];
                double held = sigmaPrior.Sample(rng);
                var chain = new BassChain(prior, state.Terms, state.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
                if (chain.Follow(At(prior, v, held, beta)) is not null) continue;
                BassSweepFit fit = chain.Structure.DrawSweepFit(chain.Current, rng, prior.Target, power: 0.0);
                return (chain.Current, held, v, fit.Fitted);
            }
        }
        var marginal = new double[draws][];
        for (int k = 0; k < draws; k++)
        {
            var d = Draw();
            marginal[k] = Functions(d.S, d.Sigma2, d.V, d.F, Respond(d.F, d.Sigma2, d.V));
        }

        var theta = Draw();
        var successive = new double[steps][];
        long refused = 0, blockRefused = 0;
        for (int t = 0; t < steps; t++)
        {
            double[] y = Respond(theta.F, theta.Sigma2, theta.V);
            BassModel model = BassModel.Create(config, Data(y));
            var moves = new BassStructureMoves(model, theta.S.Terms, theta.S.LogTau);
            var blockRng = new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue));
            var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
            if (BassStructureMoves.Follow(chain, At(model, theta.V, theta.Sigma2, beta)) is not null) { refused++; continue; }
            chain.Step();
            foreach (IBlock<ReversibleJumpChain<BassState>> block in BassSweepSlots.Blocks(model, BassSweepSlot.Tau))
            {
                if (block.Name == "scale") continue;   // σ moves only by the interweaving step
                if (!block.Update(chain, blockRng)) blockRefused++;
            }
            BassSweepFit fit = moves.SweepFit!;
            ObservationLatents latents = chain.Current.Fit.Working!.Latents!;
            theta = (chain.Current, latents.Sigma2!.Value, latents.ValueArray, fit.Fitted);
            successive[t] = Functions(theta.S, theta.Sigma2, theta.V, theta.F, y);
        }

        double worst = 0.0;
        for (int g = 0; g < AloneNames.Length; g++)
        {
            double[] a = [.. marginal.Select(r => r[g])], s = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1), meanB = s.Average();
            int size = s.Length / batches;
            double[] means = [.. Enumerable.Range(0, batches).Select(k => s.Skip(k * size).Take(size).Average())];
            double seB2 = means.Sum(m => (m - means.Average()) * (m - means.Average())) / (batches - 1) / batches;
            double z = (meanA - meanB) / Math.Sqrt(varA / a.Length + seB2);
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{AloneNames[g],-10} marginal-conditional {meanA,10:F5}  successive {meanB,10:F5}  z {z,6:F2}");
        }
        output.WriteLine($"{law}, the step alone: refused {refused}, blocks refused {blockRefused}");
        Assert.True(worst < 4.0, $"worst |z| {worst:F2}");
    }

    // ln σ², β, mean ln v and M of a retained draw.
    private static double Functional(BassDraw draw, int g) => g switch
    {
        0 => Math.Log(draw.Sigma2),
        1 => draw.Latents!.Skew,
        2 => draw.Latents!.ValueArray.Average(Math.Log),
        _ => draw.Structure.TermIds.Length,
    };
}
