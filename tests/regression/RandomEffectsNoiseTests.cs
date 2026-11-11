using System;
using System.Collections.Generic;
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
/// E21's slice 4, random effects under GBASS's noises (<c>issues/regression/bass/ridge-blocks-brief.md</c>, decision 7 as amended on
/// 2026-11-10): u is a latent the reading carries as the offset Z̃u, with σ² on the record; at b = 1, when the fit opts in
/// (<see cref="BassSampling.CollapseEffectsInMoves"/>), the structure step integrates it under the weights and draws it afresh.
/// <list type="bullet">
/// <item><b>The weighted whitening against dense algebra:</b> L = V_w^(−½)W^½ gives xᵀC⁻¹y with C = W⁻¹ + Z̃S²Z̃ᵀ, the GLS intercept's mass
/// and mean through c = L·1 (the convention centres u with unweighted counts, so c is not W^½·1), ln|V_w| = ln|C| + Σ ln wᵢ, and u's
/// conditional mean, at the weights and at the rung's bW; at W = I the unweighted block's whitening.</item>
/// <item><b>The design on the integrated reading is the GLS fit:</b> its residual squares, mean, mass and Σ ln w against dense algebra.</item>
/// <item><b>The model under a noise:</b> σ² on the record, u = 0 at the start, the sweep's blocks in their order, a fit's draws carrying
/// u and s; a tempered fit with a drawn scale and the reference rung resumes bit for bit, with the collapse and without; fits without a
/// block keep their draws.</item>
/// <item><b>The collapse is an opt-in:</b> every sampler's sweep without it is form B's, block by block; a chain, a ladder and SMC each
/// take it when asked; it is refused where it reads nothing and left out of the serialized form when off.</item>
/// </list>
/// The joint law by Geweke at b = 1 and below, the mixing measurement and recovery are behind <c>BASS_QUALIFICATION</c>
/// (<c>RandomEffectsNoiseQualificationTests</c>).
/// </summary>
public sealed class RandomEffectsNoiseTests(ITestOutputHelper output)
{
    private const int Cases = 30;

    // Two crossed intercept groupings and, where asked, a slope on x0 within the first, on 30 cases.
    private static (RandomEffectsBlock Block, double[,] Z, double[] ColumnScales) Fixture(bool slope)
    {
        var rng = new Xoshiro256PlusPlus(5);
        ContinuousPredictor[] predictors =
        [
            new("x0", [.. Enumerable.Range(0, Cases).Select(_ => rng.NextDouble())]),
            new("x1", [.. Enumerable.Range(0, Cases).Select(_ => rng.NextDouble())]),
        ];
        CategoricalPredictor[] groupings = [RandomInterceptTests.Grouping(Cases, 5, 11, "g1"), RandomInterceptTests.Grouping(Cases, 3, 12, "g2")];
        List<IBassRandomEffect> specs = [new BassRandomIntercept("g1", 0.7), new BassRandomIntercept("g2", 1.4)];
        if (slope) specs.Add(new BassRandomSlope("g1", "x0", RandomSlopeScale.Unit, 0.5));
        RandomEffectsBlock block = RandomEffectsBlock.Of(groupings, specs, predictors);
        int q = Enumerable.Range(0, block.Count).Sum(k => block.Layout(k).Dimension);
        var z = new double[Cases, q];
        var scales = new double[q];
        int column = 0;
        for (int k = 0; k < block.Count; k++)
        {
            RandomEffectGrouping layout = block.Layout(k);
            for (int j = 0; j < layout.Dimension; j++, column++)
            {
                scales[column] = block.Scales.Span[k];
                for (int i = 0; i < Cases; i++) z[i, column] = layout.Multiplier(i) * layout.Contrast(layout.GroupOf(i), j);
            }
        }
        return (block, z, scales);
    }

    // C = W⁻¹ + Z̃S²Z̃ᵀ, written out.
    private static double[,] Covariance(double[,] z, double[] scales, double[] weights)
    {
        int n = weights.Length, q = scales.Length;
        var c = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            c[i, i] = 1.0 / weights[i];
            for (int j = 0; j < n; j++)
                for (int k = 0; k < q; k++) c[i, j] += z[i, k] * scales[k] * scales[k] * z[j, k];
        }
        return c;
    }

    // The dense Cholesky factor of an SPD matrix, for solves and the log determinant.
    private static double[,] Cholesky(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = a[i, j];
                for (int k = 0; k < j; k++) sum -= l[i, k] * l[j, k];
                l[i, j] = i == j ? Math.Sqrt(sum) : sum / l[j, j];
            }
        return l;
    }

    private static double[] Solve(double[,] l, double[] b)
    {
        int n = b.Length;
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            double sum = b[i];
            for (int k = 0; k < i; k++) sum -= l[i, k] * y[k];
            y[i] = sum / l[i, i];
        }
        var x = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = y[i];
            for (int k = i + 1; k < n; k++) sum -= l[k, i] * x[k];
            x[i] = sum / l[i, i];
        }
        return x;
    }

    private static double LogDeterminant(double[,] l)
    {
        double sum = 0.0;
        for (int i = 0; i < l.GetLength(0); i++) sum += 2.0 * Math.Log(l[i, i]);
        return sum;
    }

    private static double Dot(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        double sum = 0.0;
        for (int i = 0; i < a.Length; i++) sum += a[i] * b[i];
        return sum;
    }

    private static double Rel(double a, double b) => Math.Abs(a - b) / Math.Max(1e-300, Math.Max(Math.Abs(a), Math.Abs(b)));

    [Theory]
    [InlineData(false, 1.0)]
    [InlineData(true, 1.0)]
    [InlineData(true, 0.3)]
    public void TheWeightedWhitening_IsTheDenseGls(bool slope, double power)
    {
        (RandomEffectsBlock block, double[,] z, double[] scales) = Fixture(slope);
        var rng = new Xoshiro256PlusPlus(44);
        double[] w = [.. Enumerable.Range(0, Cases).Select(_ => power * (0.3 + 2.7 * rng.NextDouble()))];
        double[] x = [.. Enumerable.Range(0, Cases).Select(_ => Normal.SampleStandard(rng))];
        double[] y = [.. Enumerable.Range(0, Cases).Select(_ => 2.0 + Normal.SampleStandard(rng))];
        double[] ones = [.. Enumerable.Repeat(1.0, Cases)];
        double[,] l = Cholesky(Covariance(z, scales, w));
        var whitening = new WeightedWhitening(block, w);

        double inner = Rel(Dot(whitening.Apply(x), whitening.Apply(y)), Dot(x, Solve(l, y)));
        double mass = Dot(ones, Solve(l, ones));
        double massGap = Rel(whitening.Mass, mass), meanGap = Rel(whitening.Mean(y), Dot(ones, Solve(l, y)) / mass);
        double logGap = Math.Abs(whitening.LogDeterminant - (LogDeterminant(l) + w.Sum(Math.Log)));
        // The direction is not W^½·1: the convention centres u with unweighted counts.
        double offAxis = Enumerable.Range(0, Cases).Max(i => Math.Abs(whitening.Direction[i] - Math.Sqrt(w[i])));

        // u's conditional mean given e: in contrasts S²Z̃ᵀC⁻¹e, each grouping's effects through its orthonormal contrasts.
        double[] e = [.. Enumerable.Range(0, Cases).Select(_ => Normal.SampleStandard(rng))];
        double[] solved = Solve(l, e);
        var coordinates = new double[scales.Length];
        for (int k = 0; k < scales.Length; k++)
            for (int i = 0; i < Cases; i++) coordinates[k] += scales[k] * scales[k] * z[i, k] * solved[i];
        double[][] mean = whitening.MeanEffects(e);
        double effectGap = 0.0;
        for (int k = 0, start = 0; k < block.Count; start += block.Layout(k).Dimension, k++)
        {
            double[] expected = block.Layout(k).Effects(coordinates.AsSpan(start, block.Layout(k).Dimension));
            for (int g = 0; g < expected.Length; g++) effectGap = Math.Max(effectGap, Math.Abs(mean[k][g] - expected[g]) / (1.0 + Math.Abs(expected[g])));
        }
        output.WriteLine($"slope {slope}, power {power}: xᵀC⁻¹y {inner:E2}, mass {massGap:E2}, mean {meanGap:E2}, ln|V_w| {logGap:E2}, " +
                         $"effects {effectGap:E2}; c departs from W^½·1 by {offAxis:E2}");
        Assert.True(inner < 1e-10 && massGap < 1e-10 && meanGap < 1e-10 && logGap < 1e-9 && effectGap < 1e-10);
        Assert.True(offAxis > 1e-3, "the GLS intercept direction differs from W^½·1");

        // At W = I, the unweighted block's whitening, c = 1 and the plain mean.
        var unit = new WeightedWhitening(block, ones);
        double[] copy = (double[])x.Clone();
        block.Whiten(copy);
        double[] applied = unit.Apply(x);
        Assert.True(Enumerable.Range(0, Cases).All(i => Math.Abs(applied[i] - copy[i]) <= 1e-12 * (1.0 + Math.Abs(copy[i]))));
        Assert.Equal(y.Average(), unit.Mean(y), 12);
    }

    // A t-noise model with two crossed groupings, the first's scale drawn, on 60 cases.
    private static (BassModel Model, BassConfig Config) NoiseModel(BassNoise? noise = null, int cases = 60)
    {
        BassData plain = BassQuantileTests.Data();
        int n = Math.Min(cases, plain.Response.Cases);
        var data = new BassData(plain.Predictors, plain.Response,
                                [RandomInterceptTests.Grouping(plain.Response.Cases, 5, 21, "g1"), RandomInterceptTests.Grouping(plain.Response.Cases, 3, 22, "g2")]);
        BassConfig config = BassStudentTTests.Config(noise ?? BassNoise.StudentT(4)) with
        {
            RandomIntercepts = [new BassRandomIntercept("g1", 0.8, ScalePrior: new BassGroupScalePrior()), new BassRandomIntercept("g2", 0.5)],
        };
        return (BassModel.Create(config, data), config);
    }

    [Fact]
    public void OnTheIntegratedReading_TheDesignIsTheGlsFit()
    {
        (BassModel model, _) = NoiseModel();
        var mixture = (VarianceMeanMixture)model.Observation;
        int n = model.Data.Response.Observations;
        var rng = new Xoshiro256PlusPlus(9);
        double[] v = [.. Enumerable.Range(0, n).Select(_ => 0.3 + 2.0 * rng.NextDouble())];
        WorkingResponse reading = mixture.At(model.Working, new ObservationLatents(v, 0.0, 0.8));
        WorkingResponse integrated = mixture.Integrated(reading);
        BassDesign design = model.NewDesign(integrated);
        var terms = new List<BassTerm>();
        for (int k = 0; k < 200 && design.Count < 3; k++)
        {
            var counts = new BassTermCounts(model.Prior);
            foreach (BassTerm t in design.Terms) counts.Add(t);
            if (model.Proposal.DrawBirth(rng, counts) is { } term && design.ProposeBirth(term).Admitted) design.Accept();
        }
        Assert.True(design.Count >= 2);

        // Dense GLS: P = [1, the raw columns], C = W⁻¹ + Z̃S²Z̃ᵀ with W = 1/(c·v).
        RandomEffectsBlock block = model.RandomEffects!;
        double[] w = [.. v.Select(vi => 1.0 / (mixture.VarianceScale * vi))];
        int q = Enumerable.Range(0, block.Count).Sum(k => block.Layout(k).Dimension);
        var z = new double[n, q];
        var scales = new double[q];
        for (int k = 0, column = 0; k < block.Count; k++)
            for (int j = 0; j < block.Layout(k).Dimension; j++, column++)
            {
                scales[column] = block.Scales.Span[k];
                for (int i = 0; i < n; i++) z[i, column] = block.Layout(k).Multiplier(i) * block.Layout(k).Contrast(block.Layout(k).GroupOf(i), j);
            }
        double[,] l = Cholesky(Covariance(z, scales, w));
        int p = design.Count + 1;
        var columns = new double[p][];
        columns[0] = [.. Enumerable.Repeat(1.0, n)];
        for (int j = 0; j < design.Count; j++)
        {
            columns[j + 1] = new double[n];
            model.Basis.CaseColumn(design.Terms[j], model.Inputs.Cases, columns[j + 1]);
        }
        double[] y = model.Data.Response.Values.ToArray();
        double[][] solvedColumns = [.. columns.Select(c => Solve(l, c))];
        var gram = new double[p, p];
        var cross = new double[p];
        for (int a = 0; a < p; a++)
        {
            cross[a] = Dot(solvedColumns[a], y);
            for (int b = 0; b < p; b++) gram[a, b] = Dot(columns[a], solvedColumns[b]);
        }
        double[] beta = Solve(Cholesky(gram), cross);
        double rss = Dot(y, Solve(l, y)) - Dot(beta, cross);
        double mass = Dot(columns[0], solvedColumns[0]);
        double rssGap = Rel(design.ResidualSquares, rss), meanGap = Rel(design.ResponseMean, Dot(solvedColumns[0], y) / mass);
        double massGap = Rel(design.Weighting.Mass, mass);
        double logGap = Math.Abs(design.Weighting.SumLogWeights - (w.Sum(Math.Log) - (LogDeterminant(l) + w.Sum(Math.Log))));
        output.WriteLine($"{design.Count} terms: ‖r‖² {rssGap:E2}, mean {meanGap:E2}, mass {massGap:E2}, Σ ln w − ln|V_w| {logGap:E2}");
        Assert.True(rssGap < 1e-10 && meanGap < 1e-10 && massGap < 1e-10 && logGap < 1e-9);
    }

    [Fact]
    public void UnderANoise_TheEffectsAreALatent_AndTheSweepReadsThem()
    {
        (BassModel model, BassConfig config) = NoiseModel();
        var mixture = (VarianceMeanMixture)model.Observation;
        Assert.True(model.EffectsAreLatent && mixture.RecordsSigma2);
        Assert.NotNull(model.Working.Latents!.Sigma2);
        Assert.Null(model.Working.Whitening);
        Assert.All(model.Working.Effects!.Effects, effects => Assert.All(effects, u => Assert.Equal(0.0, u)));
        string[] chain = [.. BassSweepSlots.Blocks(model).Select(b => b.Name)];
        string[] replica = [.. BassSweepSlots.Blocks(model, BassSweepSlot.Tau).Select(b => b.Name)];
        output.WriteLine($"chain: {string.Join(", ", chain)}");
        output.WriteLine($"replica: {string.Join(", ", replica)}");
        Assert.Equal("moves", chain[0]);   // form B by default: the move step on the likelihood given u
        Assert.Equal("tau", replica[0]);
        Assert.True(Array.IndexOf(chain, "group effects") == Array.IndexOf(chain, "coefficients") + 1, "u is drawn right after a");
        Assert.Contains("group scale g1", chain);
        Assert.DoesNotContain("group scale g2", chain);

        BassFit fit = Bass.Fit(config with { Sampling = new BassSampling { Chains = 2, MasterSeed = 7, BurnIn = 300, Samples = 300 } }, model.Data);
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        Assert.All(draws, d => Assert.Equal(2, d.Groups!.Count));
        Assert.True(draws.Select(d => d.Groups![0].Effects[0]).Distinct().Count() > draws.Length / 2, "u moves");
        Assert.True(draws.Select(d => d.Groups![0].Scale).Distinct().Count() > 10, "s moves");
        Assert.DoesNotContain(fit.Diagnostics.Chains[0].Blocks!, k => k.Key == "collapsed moves");
        // Each draw's effects keep the convention, and its reading rebuilds with its u.
        RandomEffectsBlock block = model.RandomEffects!;
        foreach (BassDraw d in draws.Take(20))
        {
            BassGroupDraw[] groups = block.Ordered(d.Groups!);
            for (int k = 0; k < block.Count; k++)
            {
                double centred = Enumerable.Range(0, block.Layout(k).Levels).Sum(g => block.Layout(k).Counts[g] * groups[k].Effects[g]);
                Assert.True(Math.Abs(centred) < 1e-9, $"Σ n_g u_g {centred:E2}");
            }
            WorkingResponse reading = d.WorkingOn(model);
            Assert.Equal(groups[0].Effects, reading.Effects!.Effects[0]);
        }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void ATemperedFitUnderANoise_WithADrawnScale_ResumesBitForBit(bool collapse)
    {
        (BassModel model, BassConfig config) = NoiseModel();
        var spec = new TemperingSpec { Levels = 3, HotEnd = 0.2, ReferenceRung = true };
        BassConfig whole = config with
        {
            Sigma2 = Sigma2Prior.InverseGamma(3.0, 0.5),
            Sampling = new BassSampling
            {
                Chains = 2, MasterSeed = 31, BurnIn = 100, Samples = 120, Thin = 2, Tempering = spec, ReportEvidence = true,
                CollapseEffectsInMoves = collapse,
            },
        };
        BassConfig half = whole with { Sampling = whole.Sampling with { Samples = 60 } };
        BassFit full = Bass.Fit(whole, model.Data);
        BassFit first = Bass.Fit(half, model.Data, new BassRunOptions { Checkpoint = true });
        BassFit resumed = Bass.Resume(first.Checkpoint!, half, model.Data);
        BassRunTests.AssertSameFit(full, resumed);
        for (int c = 0; c < full.Draws.Count; c++)
            for (int d = 0; d < full.Draws[c].Count; d++)
                for (int k = 0; k < 2; k++)
                {
                    Assert.Equal(full.Draws[c][d].Groups![k].Scale, resumed.Draws[c][d].Groups![k].Scale);
                    Assert.Equal(full.Draws[c][d].Groups![k].Effects, resumed.Draws[c][d].Groups![k].Effects);
                }
        Assert.True(double.IsFinite(full.Evidence!.LogEvidence));
    }

    // The b = 1 collapse is an opt-in (decision 7, amended 2026-11-10). Chains, ladder replicas and SMC particles all take their sweep from
    // BassSweepSlots.Blocks, so without the option each is form B's, block by block; with it, each sampler takes the collapsed step.
    [Fact]
    public void TheCollapse_IsAnOptIn_AndEveryFitWithoutIt_TakesFormBsSweep()
    {
        (BassModel model, BassConfig config) = NoiseModel();
        BassConfig opted = config with { Sampling = config.Sampling with { CollapseEffectsInMoves = true } };
        BassModel collapsing = BassModel.Create(opted, model.Data);
        static string[] Read(IEnumerable<IBlock<ReversibleJumpChain<BassState>>> blocks) => [.. blocks.Select(b => $"{b.GetType().Name}:{b.Name}")];
        foreach (BassSweepSlot first in new[] { BassSweepSlot.Structure, BassSweepSlot.Tau })
        {
            string[] plain = Read(BassSweepSlots.Blocks(model, first));
            Assert.Equal(Read(BassSweepSlots.Blocks(model, first, collapse: false)), plain);
            Assert.DoesNotContain(plain, b => b.StartsWith(nameof(BassCollapsedStructureBlock), StringComparison.Ordinal));
            string[] collapsed = Read(BassSweepSlots.Blocks(collapsing, first));
            Assert.Equal(Read(BassSweepSlots.Blocks(collapsing, first, collapse: true)), collapsed);
            Assert.StartsWith(nameof(BassCollapsedStructureBlock), collapsed[0], StringComparison.Ordinal);
            Assert.Equal(plain.Skip(first == BassSweepSlot.Structure ? 1 : 0), collapsed.Skip(1));
        }

        // Each sampler: off, stated or not, is one fit; on differs, so the option reaches it.
        var chains = new BassSampling { Chains = 2, MasterSeed = 7, BurnIn = 150, Samples = 150 };
        BassConfig proper = config with { Sigma2 = Sigma2Prior.InverseGamma(3.0, 0.5) };
        foreach ((string name, BassSampling sampling) in new (string, BassSampling)[]
                 {
                     ("chain", chains),
                     ("ladder", chains with { Tempering = new TemperingSpec { Levels = 3, HotEnd = 0.2, ReferenceRung = true } }),
                     ("SMC", new BassSampling { Chains = 1, MasterSeed = 7, Smc = new SmcSpec { Particles = 60, Replicates = 2, MutationSteps = 2 } }),
                 })
        {
            BassFit byDefault = Bass.Fit(proper with { Sampling = sampling }, model.Data);
            BassFit stated = Bass.Fit(proper with { Sampling = sampling with { CollapseEffectsInMoves = false } }, model.Data);
            BassFit on = Bass.Fit(proper with { Sampling = sampling with { CollapseEffectsInMoves = true } }, model.Data);
            BassRunTests.AssertSameFit(byDefault, stated);
            BassDraw[] a = [.. byDefault.Draws.SelectMany(c => c)], b = [.. on.Draws.SelectMany(c => c)];
            Assert.True(a.Zip(b).Any(p => !p.First.Groups![0].Effects.SequenceEqual(p.Second.Groups![0].Effects)), $"{name}: the option reaches the sampler");
            if (name == "chain")
            {
                Assert.Equal(Read(BassSweepSlots.Blocks(BassModel.Create(proper with { Sampling = sampling }, model.Data))).Select(k => k.Split(':')[1]),
                             byDefault.Diagnostics.Chains[0].Blocks!.Select(k => k.Key));
                MoveStat moves = on.Diagnostics.Chains[0].Blocks!.Single(k => k.Key == "collapsed moves");
                output.WriteLine($"collapsed moves accepted {moves.Accepted} of {moves.Attempts}");
                Assert.True(moves.Accepted > 0);
            }
        }

        // Off is left out of the serialized form, so configs and fingerprints from before the option read as they did; on round-trips.
        Assert.DoesNotContain(nameof(BassSampling.CollapseEffectsInMoves), config.ToJson());
        Assert.Contains($"\"{nameof(BassSampling.CollapseEffectsInMoves)}\":true", opted.ToJson());
        Assert.Equal(opted, BassConfig.FromJson(opted.ToJson()));
        // Refused where it reads nothing: under Gaussian noise u is always integrated, and without random effects there is no u.
        Assert.Contains("CollapseEffectsInMoves", Assert.Throws<ArgumentException>(() => (opted with { Noise = null }).Validate()).Message);
        Assert.Contains("CollapseEffectsInMoves", Assert.Throws<ArgumentException>(() => (opted with { RandomIntercepts = null }).Validate()).Message);
    }

    // Noise fits without a block keep their draws: pinned by hash on HEAD before the change.
    [Theory]
    [InlineData("t", "5471c4bcad7d0b09")]
    [InlineData("skewed horseshoe", "afd1aa3f399ea35f")]
    [InlineData("tempered normal-wald", "98b9820e14563b71")]
    public void ANoiseFitWithoutABlock_KeepsItsDraws(string fit, string expected)
    {
        BassData data = BassQuantileTests.Data();
        var sampling = new BassSampling { Chains = 2, MasterSeed = 2110, BurnIn = 200, Samples = 150 };
        BassConfig config = fit switch
        {
            "t" => BassStudentTTests.Config(BassNoise.StudentT(4)) with { Sampling = sampling },
            "skewed horseshoe" => BassStudentTTests.Config(BassNoise.Horseshoe(new BassNoisePrior(0.0, 0.5))) with
            {
                Sigma2 = Sigma2Prior.InverseGamma(3.0, 0.5), Sampling = sampling,
            },
            _ => BassStudentTTests.Config(BassNoise.NormalWald(skew: new BassNoisePrior(0.0, 1.0))) with
            {
                Sigma2 = Sigma2Prior.InverseGamma(3.0, 0.5),
                Sampling = sampling with { Tempering = new TemperingSpec { Levels = 3, HotEnd = 0.2, ReferenceRung = true } },
            },
        };
        BassFit result = Bass.Fit(config, data);
        var hash = new RunFingerprint().Add(BassBitIdentityTests.Hash(result));
        foreach (DrawStore<BassDraw> chain in result.Draws)
            foreach (BassDraw d in chain)
                if (d.Latents is { } latents) hash.Add(latents.Values).Add(latents.Skew).Add(latents.Sigma2 ?? 0.0);
        output.WriteLine($"{fit}: {hash}");
        Assert.Equal(expected, hash.ToString());
    }
}
