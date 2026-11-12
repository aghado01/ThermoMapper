using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E34's slice 4, stage 4: random intercepts for several outputs under a noise by row (decisions 6 and 7), U a latent on every rung (form B),
/// u_g ∈ ℝᴷ ~ N(0, s²Σ) with one s for every output, Σ integrated in the structure step with U's prior inside its integral as pseudo-rows.
/// <list type="bullet">
/// <item><b>The target is equivariant in the outputs:</b> Y → YT with γ → Tᵀγ and U → UT shifts it by −(bN − f + 1 + q)·ln|det T| at fixed
/// v and s, q the grouping's contrasts.</item>
/// <item><b>U's block draws its conditional:</b> mean and covariance across levels and outputs against the dense Gaussian conditional, the
/// prior s²Π on the observation-centred levels and the likelihood at bW, by Monte Carlo at b = 1 and 0.5.</item>
/// <item><b>The readers:</b> prediction adds an observed level's effects from the draw's U, and the mean surface reads the draw's weights.</item>
/// <item><b>The prior sampler</b> draws γ and U from one Σ and v from its law, their second moments against the closed forms; SMC and a ladder
/// with its reference rung run end to end; a restored chain continues bit for bit.</item>
/// <item><b>Geweke from exact starts</b> (qualification) under a drawn s: the t and the skewed t at K = 2, the skewed t at K = 1, and the
/// skewed t at a tilted b.</item>
/// <item><b>A ladder and SMC</b> (qualification) agree on the evidence with each other, and on posterior means with chains.</item>
/// <item><b>K = 1 against scalar BASS by law</b> (qualification): σ² on the record there, Σ integrated here, the same model; posterior means
/// of ln s, mean ln v and the surface agree within their Monte Carlo errors.</item>
/// </list>
/// </summary>
public sealed class MultiOutputNoiseEffectsTests(ITestOutputHelper output)
{
    private const int Seed = 2026_1112;

    private static string[] Names(int k) => [.. Enumerable.Range(0, k).Select(l => $"y{l}")];

    private static BassData Several(BassPredictor[] p, CategoricalPredictor g, double[,] y) => new(p, new MultiOutputResponse(y, Names(y.GetLength(1))), [g]);

    private static BassConfig Config(BassNoise noise, double s, bool drawn, int samples = 60, int seed = 17, int chains = 2, int burn = 100)
        => BassRunTests.Config(chains: chains, burn: burn, samples: samples, seed: seed, maxBasis: 6) with
        {
            Noise = noise,
            RandomIntercepts = [new BassRandomIntercept("g", s, ScalePrior: drawn ? new BassGroupScalePrior(1.0) : null)],
        };

    // [grouping 0][output][level] from a levels × K array.
    private static double[][][] Effects(double[,] u)
        => [[.. Enumerable.Range(0, u.GetLength(1)).Select(l => Enumerable.Range(0, u.GetLength(0)).Select(h => u[h, l]).ToArray())]];

    // The reading at latents v, γ and the tail, then U at scale s.
    private static WorkingResponse Reading(BassModel model, double[] v, double[]? gamma, double? tail, double[,] u, double s)
    {
        var mixture = (VarianceMeanMixture)model.Observation;
        WorkingResponse latents = mixture.At(model.Working, new ObservationLatents(v, 0.0, null, tail) { OutputSkew = gamma });
        return mixture.At(latents, RandomEffectsLatent.OfOutputs(model.RandomEffects!.At(0, s), Effects(u)));
    }

    // ── The target's equivariance in the outputs ──────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(GPriorCentre.Mean, 1.0, false)]
    [InlineData(GPriorCentre.Zero, 0.4, true)]
    [InlineData(GPriorCentre.Flat, 1.0, true)]
    public void Target_WithEffectsUnderANoise_IsEquivariantInTheOutputs(GPriorCentre centre, double b, bool drawn)
    {
        const int n = 40, k = 3, levels = 6;
        const double s = 0.8;
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = MultiOutputRandomEffectsTests.Data(n, k, levels, Seed);
        var rng = new Xoshiro256PlusPlus(Seed + 1);
        var map = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int c = 0; c < k; c++) map[a, c] = (a == c ? 1.0 : 0.0) + 0.35 * (2.0 * rng.NextDouble() - 1.0);
        double det = map[0, 0] * (map[1, 1] * map[2, 2] - map[1, 2] * map[2, 1]) - map[0, 1] * (map[1, 0] * map[2, 2] - map[1, 2] * map[2, 0])
                     + map[0, 2] * (map[1, 0] * map[2, 1] - map[1, 1] * map[2, 0]);
        double logDet = Math.Log(Math.Abs(det));
        BassConfig config = Config(BassNoise.SkewedT(4.0, new BassNoisePrior(0.0, 0.8)), s, drawn) with { Intercept = centre };
        BassModel model = BassModel.Create(config, Several(p, g, y)), image = BassModel.Create(config, Several(p, g, ReferenceAlgebra.Multiply(y, map)));
        double[] v = [.. Enumerable.Range(0, n).Select(_ => 0.3 + 2.5 * rng.NextDouble())];
        double[] gamma = [0.12, -0.07, 0.05];
        double[] movedGamma = [.. Enumerable.Range(0, k).Select(c => Enumerable.Range(0, k).Sum(a => gamma[a] * map[a, c]))];
        double[,] root = ReferenceAlgebra.Cholesky(new double[k, k].Select2((a, c) => a == c ? 0.05 : 0.01))!;
        double[,] u = MultiOutputRandomEffectsTests.CentredEffects(g.Codes.ToArray(), levels, k, s, root, rng);
        double[,] movedU = ReferenceAlgebra.Multiply(u, map);
        int q = model.RandomEffects!.Layout(0).Dimension;
        var sampler = new BassPriorSampler(model);
        var temperature = new Temperature(b, TemperingPath.Likelihood);
        double Value(BassModel m, BassTerm[] terms, double tau, double[] skew, double[,] effects, double scale)
        {
            BassDesign d = m.NewDesign();
            Assert.Null(d.Load(terms, Reading(m, v, skew, null, effects, scale)));
            var state = new BassState([.. d.Terms], BassStructureFit.Of(d, 0.0), Math.Log(tau), double.NaN);
            return temperature.LogDensity(m.Target.Evaluate(state));
        }
        double worst = 0.0, scaleMoves = 0.0;
        double[] taus = [1e-3, 0.05, 1.0, 20.0];
        for (int c = 0; c < 10; c++)
        {
            BassTerm[] terms = [.. sampler.DrawState(rng).State.Terms];
            double tau = taus[c % taus.Length];
            double at = Value(model, terms, tau, gamma, u, s), there = Value(image, terms, tau, movedGamma, movedU, s);
            double expected = -(b * n - (centre == GPriorCentre.Flat ? 1.0 : 0.0) + 1.0 + q) * logDet;
            worst = Math.Max(worst, Math.Abs(there - at - expected) / (1.0 + Math.Abs(at)));
            scaleMoves = Math.Max(scaleMoves, Math.Abs(Value(model, terms, tau, gamma, u, 2.0 * s) - at));
        }
        output.WriteLine($"{centre} at b = {b}, s {(drawn ? "drawn" : "held")}: the shift against −(bN − f + 1 + q)·ln|det T| to {worst:G3}; doubling s moves the target by up to {scaleMoves:G3}");
        Assert.True(worst < 1e-9, $"{worst:G3}");
        Assert.True(scaleMoves > 1e-3, "U's prior reaches the target");
    }

    // ── U's block against the dense conditional ───────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(1.0)]
    [InlineData(0.5)]
    public void TheEffectsBlock_DrawsTheirConditional(double b)
    {
        const int n = 40, k = 2, levels = 5, draws = 20_000;
        const double s = 0.7;
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = MultiOutputRandomEffectsTests.Data(n, k, levels, Seed + 10);
        BassModel model = BassModel.Create(Config(BassNoise.SkewedT(4.0, new BassNoisePrior(0.0, 1.0)), s, false), Several(p, g, y));
        var rng = new Xoshiro256PlusPlus(Seed + 11);
        int[] codes = g.Codes.ToArray();
        double[] v = [.. Enumerable.Range(0, n).Select(_ => 0.3 + 2.5 * rng.NextDouble())];
        double[] gamma = [0.1, -0.05];
        double[,] u0 = MultiOutputRandomEffectsTests.CentredEffects(codes, levels, k, s, ReferenceAlgebra.Cholesky(new double[,] { { 0.05, 0.01 }, { 0.01, 0.04 } })!, rng);
        BassTerm[] terms = [.. new BassPriorSampler(model).DrawState(rng).State.Terms];
        var moves = new BassStructureMoves(model, terms, 0.0);
        var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(3), new Temperature(b, TemperingPath.Likelihood));
        Assert.Null(BassStructureMoves.Follow(chain, Reading(model, v, gamma, null, u0, s)));
        BassOutputSweepFit fit = moves.DrawOutputSweepFit(chain.Current, rng, model.Target, b);
        double[,] sigma = fit.Sigma;
        // E = Y − vγᵀ − F, F the coefficients' part of the sweep's fit (it carries Z̃U₀).
        var e = new double[n, k];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < k; l++) e[i, l] = y[i, l] - v[i] * gamma[l] - (fit.Fitted[i * k + l] - u0[codes[i], l]);
        // The dense conditional: prior s²Π on the levels (Π the projector onto nᵀu = 0), e_l = Zu_l + ε, ε ~ N(0, diag(c·vᵢ/b)).
        var counts = new double[levels];
        foreach (int c in codes) counts[c]++;
        double nn = counts.Sum(c => c * c);
        double[,] prior = new double[levels, levels].Select2((h, j) => s * s * ((h == j ? 1.0 : 0.0) - counts[h] * counts[j] / nn));
        double[,] z = new double[n, levels].Select2((i, h) => codes[i] == h ? 1.0 : 0.0);
        double[,] pz = ReferenceAlgebra.Multiply(prior, ReferenceAlgebra.Transpose(z));
        double[,] gram = ReferenceAlgebra.Sum(ReferenceAlgebra.Multiply(z, pz), new double[n, n].Select2((i, j) => i == j ? v[i] / b : 0.0));
        double[,] gain = ReferenceAlgebra.Multiply(pz, ReferenceAlgebra.Inverse(gram));
        double[,] mean = ReferenceAlgebra.Multiply(gain, e);
        double[,] rows = ReferenceAlgebra.Sum(prior, ReferenceAlgebra.Multiply(ReferenceAlgebra.Multiply(gain, z), prior).Select2Negate());

        var block = new BassOutputGroupEffectsBlock(model);
        var sums = new double[levels, k];
        var products = new double[levels * k, levels * k];
        for (int t = 0; t < draws; t++)
        {
            Assert.True(block.Update(chain, rng));
            double[][] drawn = moves.DesignAt(chain.Current).Working.Effects!.OutputEffects![0];
            for (int h = 0; h < levels; h++)
                for (int l = 0; l < k; l++)
                {
                    double x = drawn[l][h] - mean[h, l];
                    sums[h, l] += x;
                    for (int h2 = 0; h2 < levels; h2++)
                        for (int l2 = 0; l2 < k; l2++) products[h * k + l, h2 * k + l2] += x * (drawn[l2][h2] - mean[h2, l2]);
                }
        }
        double zMean = 0.0, zCov = 0.0;
        for (int h = 0; h < levels; h++)
            for (int l = 0; l < k; l++)
            {
                double variance = rows[h, h] * sigma[l, l];
                zMean = Math.Max(zMean, Math.Abs(sums[h, l] / draws) / Math.Sqrt(variance / draws));
                for (int h2 = 0; h2 < levels; h2++)
                    for (int l2 = 0; l2 < k; l2++)
                    {
                        double expected = rows[h, h2] * sigma[l, l2], other = rows[h2, h2] * sigma[l2, l2];
                        double se = Math.Sqrt((variance * other + expected * expected) / draws);
                        zCov = Math.Max(zCov, Math.Abs(products[h * k + l, h2 * k + l2] / draws - expected) / se);
                    }
            }
        output.WriteLine($"b = {b}: {draws} draws of U, worst |z| on the means {zMean:F2}, on the covariances {zCov:F2}");
        Assert.True(zMean < 4.5 && zCov < 4.5, $"{zMean:F2}, {zCov:F2}");
    }

    // ── The readers ───────────────────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Prediction_AddsAnObservedLevelsEffects_AndTheMeanSurfaceRuns()
    {
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = MultiOutputRandomEffectsTests.Data(60, 2, 6, Seed + 20);
        BassData data = Several(p, g, y);
        BassFit fit = Bass.Fit(Config(BassNoise.SkewedT(4.0, new BassNoisePrior(0.0, 1.0)), 0.8, true, samples: 30), data);
        var at = new PredictionCases([new ContinuousColumn("x0", [0.3, 0.3]), new ContinuousColumn("x1", [0.6, 0.6]), new CategoricalColumn("g", ["L2", "L9"])]);
        BassOutputPrediction pr = BassMultiOutput.Predict(fit, data, at, new BassPredictOptions { UnseenLevels = UnseenLevel.OutsideEverySubset, Predictive = true, MeanSurface = true });
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        for (int d = 0; d < pr.Draws.Count; d++)
        {
            double[][] effects = pooled[pr.Draws[d]].Groups![0].OutputEffects!;
            for (int l = 0; l < 2; l++) Assert.Equal(effects[l][2], pr.Latent[d, 0, l] - pr.Latent[d, 1, l], 10);
        }
        Assert.True(pr.MeanSurface!.Cast<double>().All(double.IsFinite));
        Assert.All(pooled, draw => Assert.NotNull(draw.Latents?.OutputSkew));
    }

    // ── The prior sampler's latents, and the ensembles that start from them ───────────────────────────────────────────────────────

    [Fact]
    public void ThePriorSampler_DrawsTheLatentsFromTheirPrior()
    {
        const int n = 24, k = 2, levels = 6, draws = 80_000;
        const double nu0 = 20.0, s = 0.8;
        double[,] psi = { { 0.6, 0.15 }, { 0.15, 0.4 } };
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = MultiOutputRandomEffectsTests.Data(n, k, levels, Seed + 50);
        BassConfig config = Config(BassNoise.SkewedT(6.0, new BassNoisePrior(0.0, 1.0)), s, false) with
        {
            Intercept = GPriorCentre.Zero,
            OutputCovariance = OutputCovariancePrior.InverseWishart(nu0, psi),
        };
        BassModel model = BassModel.Create(config, Several(p, g, y));
        var sampler = new BassPriorSampler(model);
        var rng = new Xoshiro256PlusPlus(Seed + 51);
        int q = model.RandomEffects!.Layout(0).Dimension;
        // E[γγᵀ] = s_γ²·E[Σ] and E[Σ_h u_h u_hᵀ] = s²·q·E[Σ], E[Σ] = Ψ/(ν − K − 1); 1/v ~ Gamma(3, 3) under t₆'s mixing, mean 1.
        double[,] meanSigma = psi.Select2((a, d) => psi[a, d] / (nu0 - k - 1));
        (string Name, double Expected)[] rows =
        [
            ("g1^2", meanSigma[0, 0]), ("g2^2", meanSigma[1, 1]), ("g1 g2", meanSigma[1, 0]),
            ("U0'U0", s * s * q * meanSigma[0, 0]), ("U1'U1", s * s * q * meanSigma[1, 1]), ("U0'U1", s * s * q * meanSigma[1, 0]), ("mean 1/v", 1.0),
        ];
        var values = new double[rows.Length][];
        for (int f = 0; f < rows.Length; f++) values[f] = new double[draws];
        for (int t = 0; t < draws; t++)
        {
            sampler.DrawPrior(rng);
            WorkingResponse reading = sampler.LastDesign().Response!.Working;
            double[] gamma = reading.Latents!.OutputSkew!;
            double[][] u = reading.Effects!.OutputEffects![0];
            double[] read = [gamma[0] * gamma[0], gamma[1] * gamma[1], gamma[0] * gamma[1], u[0].Sum(x => x * x), u[1].Sum(x => x * x),
                             u[0].Zip(u[1]).Sum(x => x.First * x.Second), reading.Latents.ValueArray.Average(x => 1.0 / x)];
            for (int f = 0; f < rows.Length; f++) values[f][t] = read[f];
        }
        double worst = 0.0;
        for (int f = 0; f < rows.Length; f++)
        {
            double mean = values[f].Average(), sd = Math.Sqrt(values[f].Sum(x => (x - mean) * (x - mean)) / (draws - 1));
            double z = (mean - rows[f].Expected) / (sd / Math.Sqrt(draws));
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{rows[f].Name,-9} {mean:F5} against {rows[f].Expected:F5}, z {z:F2}");
        }
        Assert.True(worst < 4.0, $"{worst:F2}");
    }

    [Fact]
    public void SmcAndALadder_RunUnderANoiseWithEffects()
    {
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = MultiOutputRandomEffectsTests.Data(60, 2, 6, Seed + 60);
        BassData data = Several(p, g, y);
        // An informative Σ prior, as E21's ensembles declare an informative σ² prior: γ and U are drawn from Σ's prior at b = 0.
        BassConfig config = Config(BassNoise.SkewedT(6.0, new BassNoisePrior(0.0, 1.0)), 0.8, true) with
        {
            OutputCovariance = OutputCovariancePrior.InverseWishart(6.0, new double[,] { { 0.2, 0.05 }, { 0.05, 0.2 } }),
        };
        BassFit smc = Bass.Fit(config with { Sampling = config.Sampling with { Chains = 1, MasterSeed = 61, ReportEvidence = true, Smc = new SmcSpec { Particles = 200, Replicates = 2, MutationSteps = 2 } } }, data);
        BassFit ladder = Bass.Fit(config with { Sampling = config.Sampling with { Chains = 1, MasterSeed = 62, BurnIn = 50, Samples = 50, ReportEvidence = true,
                                                                                 Tempering = new TemperingSpec { Levels = 3, HotEnd = 0.1, ReferenceRung = true } } }, data);
        foreach (BassFit fit in new[] { smc, ladder })
        {
            output.WriteLine($"ln p(y) {fit.Evidence!.LogEvidence:F3} ± {fit.Evidence.MonteCarloError:F3}");
            Assert.True(double.IsFinite(fit.Evidence.LogEvidence));
            Assert.All(fit.Draws.SelectMany(c => c), d =>
            {
                Assert.NotNull(d.Latents?.OutputSkew);
                Assert.NotNull(d.Groups![0].OutputEffects);
                Assert.NotNull(d.Sigma);
            });
        }
    }

    [Fact]
    public void ARestoredChain_ContinuesBitForBit()
    {
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = MultiOutputRandomEffectsTests.Data(50, 2, 6, Seed + 70);
        BassModel model = BassModel.Create(Config(MultiOutputNoiseTests.Law("skewed-horseshoe"), 0.8, true), Several(p, g, y));
        var chain = new BassChain(model, [], 0.0, new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        for (int t = 0; t < 40; t++) chain.Step();
        BassChain restored = BassChain.Restore(model, chain.Snapshot());
        for (int t = 0; t < 40; t++)
        {
            BassState a = chain.Step(), b = restored.Step();
            Assert.Equal(a.Terms, b.Terms);
            Assert.Equal(a.LogTau, b.LogTau);
            WorkingResponse ra = chain.Structure.DesignAt(a).Working, rb = restored.Structure.DesignAt(b).Working;
            Assert.Equal(ra.Latents!.ValueArray, rb.Latents!.ValueArray);
            Assert.Equal(ra.Latents.OutputSkew, rb.Latents.OutputSkew);
            Assert.Equal(ra.Effects!.Block.Scales.ToArray(), rb.Effects!.Block.Scales.ToArray());
            for (int l = 0; l < 2; l++) Assert.Equal(ra.Effects.OutputEffects![0][l], rb.Effects.OutputEffects![0][l]);
        }
    }

    // ── Geweke from exact starts (qualification) ──────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData("t", 2, 1.0)]
    [InlineData("skewed-t", 2, 1.0)]
    [InlineData("skewed-t", 1, 1.0)]
    [InlineData("skewed-t", 2, 0.5)]
    public void Geweke_FromExactStarts_WithDrawnScaleUnderANoise(string law, int k, double b)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int starts = 40_000, sweeps = 10, n = 24, levels = 6;
        const double nu0 = 20.0;
        double[,] psi = k == 2 ? new double[,] { { 0.6, 0.15 }, { 0.15, 0.4 } } : new double[,] { { 0.6 } };
        BassNoise noise = MultiOutputNoiseTests.Law(law);
        double c = noise.VarianceScale;
        var rng = new Xoshiro256PlusPlus(Seed + 30);
        (BassPredictor[] p, CategoricalPredictor g, _) = MultiOutputRandomEffectsTests.Data(n, k, levels, Seed + 31);
        int[] codes = g.Codes.ToArray();
        var scalePrior = new BassGroupScalePrior(1.0);
        BassConfig config = BassRunTests.Config(chains: 1, maxBasis: 3) with
        {
            Structure = new BassStructureSpec { MaxBasis = 3, MaxInteraction = 2 },
            Intercept = GPriorCentre.Zero,
            OutputCovariance = OutputCovariancePrior.InverseWishart(nu0, psi),
            Noise = noise,
            RandomIntercepts = [new BassRandomIntercept("g", 1.0, ScalePrior: scalePrior)],
        };
        BassModel prior = BassModel.Create(config, Several(p, g, new double[n, k].Select2((i, l) => Math.Sin(i + l))));
        var sampler = new BassPriorSampler(prior);
        // The tempered joint's exact draw: Σ tilted by |Σ|^(N(1−b)/2), each v by v^(K(1−b)/2); s, γ and U from their priors.
        var sigmaTilted = InverseWishart.Of(nu0 - n * (1.0 - b), psi);
        Func<Xoshiro256PlusPlus, double> drawV = b == 1.0
            ? r => noise.MixingAt(null).SamplePrior(r)
            : law == "skewed-t" ? r => InverseGamma.ShapeScale(3.0 - 0.5 * k * (1.0 - b), 3.0).Sample(r)
            : throw new ArgumentException("A tilted law is written for the skewed t alone.", nameof(law));
        string[] names = ["M", "ln tau", "ln s", "ln S11", "ln S22", "atanh rho", "mean ln v", "g1", "atan u(0)", "ln sum u0^2", "atan f(0)", "energy(0)"];
        var deltas = new double[names.Length][];
        for (int f = 0; f < names.Length; f++) deltas[f] = new double[starts];
        int redrawn = 0, refused = 0;
        for (int r = 0; r < starts; r++)
        {
            BassState state = sampler.DrawState(rng).State;
            double s = scalePrior.Sample(rng);
            double[,] sigma = sigmaTilted.Sample(rng);
            double[,] root = MultiOutputRandomEffectsTests.MatrixTPredictiveRoot(sigma);
            double[] v = [.. Enumerable.Range(0, n).Select(_ => drawV(rng))];
            double[]? gamma = noise.SkewDrawn ? Mixed(root, noise.SkewPrior!.Sd, k, rng) : null;
            double[,] u = MultiOutputRandomEffectsTests.CentredEffects(codes, levels, k, s, root, rng);
            BassDesign d = prior.NewDesign();
            Assert.Null(d.Load(state.Terms));   // at v = 1, U = 0 and c = 1 the prior's Gram is the design's
            double[,] f0 = MultiOutputRandomEffectsTests.GewekeSurface(prior, d, root, Math.Exp(state.LogTau), n, k, rng);
            var y = new double[n, k];
            for (int i = 0; i < n; i++)
            {
                double[] e = Mixed(root, Math.Sqrt(c * v[i] / b), k, rng);
                for (int l = 0; l < k; l++) y[i, l] = f0[i, l] + u[codes[i], l] + v[i] * (gamma?[l] ?? 0.0) + e[l];
            }
            if (y.Cast<double>().Any(x => Math.Abs(x) > 1e6))   // a condition on y alone: the test stays exact
            {
                redrawn++;
                r--;
                continue;
            }
            double[] before = Read(state, s, sigma, v, gamma, u, f0, y);
            BassModel model = BassModel.Create(config, Several(p, g, y));
            var moves = new BassStructureMoves(model, state.Terms, state.LogTau);
            var blockRng = new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue));
            var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)),
                                                           new Temperature(b, TemperingPath.Likelihood));
            if (BassStructureMoves.Follow(chain, Reading(model, v, gamma, null, u, s)) is not null)   // θ outside F: the target is the prior on F
            {
                refused++;
                r--;
                continue;
            }
            List<IBlock<ReversibleJumpChain<BassState>>> blocks = BassSweepSlots.Blocks(model);
            for (int w = 0; w < sweeps; w++)
                foreach (IBlock<ReversibleJumpChain<BassState>> block in blocks) block.Update(chain, blockRng);
            WorkingResponse reading = moves.DesignAt(chain.Current).Working;
            RandomEffectsLatent latent = reading.Effects!;
            BassOutputSweepFit fit = moves.OutputSweepFit!;
            var effects = new double[n * k];
            latent.AddTo(effects);
            double[,] fitted = new double[n, k].Select2((i, l) => fit.Fitted[i * k + l] - effects[i * k + l]);
            double[,] after = new double[levels, k].Select2((h, l) => latent.OutputEffects![0][l][h]);
            double[] read = Read(chain.Current, latent.Block.Scales.Span[0], fit.Sigma, reading.Latents!.ValueArray, reading.Latents.OutputSkew, after, fitted, y);
            for (int f = 0; f < names.Length; f++) deltas[f][r] = read[f] - before[f];
        }
        double worst = 0.0;
        for (int f = 0; f < names.Length; f++)
        {
            double mean = deltas[f].Average(), sd = Math.Sqrt(deltas[f].Sum(x => (x - mean) * (x - mean)) / (starts - 1));
            double z = sd == 0.0 ? (mean == 0.0 ? 0.0 : double.PositiveInfinity) : mean / (sd / Math.Sqrt(starts));
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{names[f],-12} change {mean,10:F5} ± {sd / Math.Sqrt(starts):F5}  z {z,6:F2}");
        }
        output.WriteLine($"{law}, K = {k}, b = {b}, exact starts ({redrawn} drawn again for |y| > 10⁶, {refused} outside F): worst |z| {worst:F2}");
        Assert.True(worst < 4.0, $"{worst:F2}");

        // The structure, τ, s, Σ, v, γ, U (bounded, and its energy in logs: E21's funnel reading), the fit at a case, the first row's energy.
        double[] Read(BassState st, double s, double[,] sigma, double[] v, double[]? gamma, double[,] u, double[,] f, double[,] y)
        {
            double[,] root = MultiOutputRandomEffectsTests.MatrixTPredictiveRoot(sigma);
            var r0 = new double[k];
            for (int l = 0; l < k; l++) r0[l] = y[0, l] - f[0, l] - u[codes[0], l] - v[0] * (gamma?[l] ?? 0.0);
            double squares = Enumerable.Range(0, levels).Sum(h => u[h, 0] * u[h, 0]);
            return [st.Count, st.LogTau, Math.Log(s), Math.Log(sigma[0, 0]), k > 1 ? Math.Log(sigma[1, 1]) : 0.0,
                    k > 1 ? Math.Atanh(sigma[1, 0] / Math.Sqrt(sigma[0, 0] * sigma[1, 1])) : 0.0, v.Average(Math.Log), gamma?[0] ?? 0.0,
                    Math.Atan(u[codes[0], 0]), Math.Log(squares), Math.Atan(f[0, 0]), b * BassOutputLatentScaleBlock.Energy(root, r0) / (c * v[0])];
        }
    }

    private static double[] Mixed(double[,] root, double scale, int k, Xoshiro256PlusPlus rng)
    {
        var z = new double[k];
        for (int l = 0; l < k; l++) z[l] = Normal.SampleStandard(rng);
        var e = new double[k];
        for (int l = 0; l < k; l++)
        {
            for (int q = 0; q <= l; q++) e[l] += root[l, q] * z[q];
            e[l] *= scale;
        }
        return e;
    }

    // ── A ladder and SMC against each other and against chains (qualification) ─────────────────────────────────────────────────────

    [Fact]
    public void ALadderAndSmc_AgreeWithEachOtherAndWithChains()
    {
        if (!BassQualification.Enabled) return;   // heavy
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = MultiOutputRandomEffectsTests.Data(80, 2, 8, Seed + 80);
        BassData data = Several(p, g, y);
        BassConfig config = Config(BassNoise.SkewedT(6.0, new BassNoisePrior(0.0, 1.0)), 0.8, true) with
        {
            Structure = new BassStructureSpec { MaxInteraction = 1, MaxBasis = 6 },
            OutputCovariance = OutputCovariancePrior.InverseWishart(6.0, new double[,] { { 0.2, 0.05 }, { 0.05, 0.2 } }),
        };
        var run = new BassRunOptions { MaxDegreeOfParallelism = 4 };
        BassFit chains = Bass.Fit(config with { Sampling = config.Sampling with { Chains = 4, MasterSeed = 3, BurnIn = 3000, Samples = 6000 } }, data, run);
        BassFit ladder = Bass.Fit(config with
        {
            Sampling = config.Sampling with
            {
                Chains = 4, MasterSeed = 5, BurnIn = 2000, Samples = 4000, ReportEvidence = true,
                Tempering = new TemperingSpec { Levels = 8, HotEnd = 0.02, ReferenceRung = true },
            },
        }, data, run);
        BassFit smc = Bass.Fit(config with
        {
            Sampling = config.Sampling with { Chains = 1, MasterSeed = 7, ReportEvidence = true, Smc = new SmcSpec { Particles = 2000, Replicates = 4, MutationSteps = 4 } },
        }, data, run);
        (double Mean, double Error) Summary(BassFit fit, Func<BassDraw, double> read)
        {
            double[] means = [.. fit.Draws.Select(chain => chain.Average(read))];
            double[]? w = fit.DrawWeights;
            double spread = Math.Sqrt(means.Select(m => (m - means.Average()) * (m - means.Average())).Sum() / Math.Max(1, means.Length - 1) / means.Length);
            if (w is null) return (means.Average(), spread);
            BassDraw[] all = [.. fit.Draws.SelectMany(c => c)];
            return (all.Select((d, i) => w[i] * read(d)).Sum() / w.Sum(), spread);
        }
        BassEvidence a = ladder.Evidence!, s = smc.Evidence!;
        output.WriteLine($"ln p(y): ladder {a.LogEvidence:F3} ± {a.MonteCarloError:F3}, SMC {s.LogEvidence:F3} ± {s.MonteCarloError:F3}");
        Assert.True(Math.Abs(a.LogEvidence - s.LogEvidence) <= 4.0 * Math.Sqrt(a.MonteCarloError * a.MonteCarloError + s.MonteCarloError * s.MonteCarloError) + 0.05,
                    $"{a.LogEvidence:F3} against {s.LogEvidence:F3}");
        if (smc.Diagnostics.Smc is { } readout) Assert.Empty(readout.Flags);
        foreach ((string name, Func<BassDraw, double> read) in new (string, Func<BassDraw, double>)[]
                 {
                     ("ln s", d => Math.Log(d.Groups![0].Scale)), ("M", d => d.Structure.TermIds.Length),
                     ("mean ln v", d => d.Latents!.Values.ToArray().Average(Math.Log)), ("g1", d => d.Latents!.OutputSkew![0]),
                 })
        {
            (double cm, double ce) = Summary(chains, read);
            foreach ((string sampler, BassFit fit) in new[] { ("ladder", ladder), ("SMC", smc) })
            {
                (double m, double e) = Summary(fit, read);
                output.WriteLine($"{name}: {sampler} {m:F4} ± {e:F4} against chains' {cm:F4} ± {ce:F4}");
                Assert.True(Math.Abs(m - cm) <= 5.0 * Math.Sqrt(e * e + ce * ce) + 0.03, $"{name}, {sampler}: {m:F4} against {cm:F4}");
            }
        }
    }

    // ── K = 1 against scalar BASS by law (qualification) ──────────────────────────────────────────────────────────────────────────

    [Fact]
    public void KOne_WithDrawnScaleUnderANoise_IsScalarBass_ByLaw()
    {
        if (!BassQualification.Enabled) return;   // heavy
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = MultiOutputRandomEffectsTests.Data(120, 1, 8, Seed + 40);
        BassConfig config = Config(BassNoise.StudentT(4.0), 1.0, true, samples: 4000, seed: 41, chains: 4, burn: 2000);
        BassData one = new(p, new ScalarResponse([.. Enumerable.Range(0, 120).Select(i => y[i, 0])]), [g]), several = Several(p, g, y);
        var run = new BassRunOptions { MaxDegreeOfParallelism = 4 };
        BassFit a = Bass.Fit(config, one, run), b = Bass.Fit(config, several, run);
        var cases = new PredictionCases([new ContinuousColumn("x0", [0.2, 0.5, 0.8]), new ContinuousColumn("x1", [0.7, 0.3, 0.5]), new CategoricalColumn("g", ["L1", "L4", "L6"])]);
        double[,] ma = a.Predict(one, cases, new BassPredictOptions { MeanSurface = true }).MeanSurface!;
        double[,] mb = BassMultiOutput.Predict(b, several, cases, new BassPredictOptions { MeanSurface = true }).MeanSurface!;
        (double Mean, double Error) Summary(BassFit fit, Func<BassDraw, double> read)
        {
            double[] means = [.. fit.Draws.Select(chain => chain.Average(read))];
            double mean = means.Average();
            return (mean, Math.Sqrt(means.Sum(x => (x - mean) * (x - mean)) / (means.Length - 1) / means.Length));
        }
        var rows = new (string Name, Func<BassDraw, double> Read)[]
        {
            ("ln s", d => Math.Log(d.Groups![0].Scale)),
            ("mean ln v", d => d.Latents!.Values.ToArray().Average(Math.Log)),
            ("ln sigma2", d => Math.Log(d.Sigma is { } s ? s[0, 0] : d.Sigma2)),
            ("M", d => d.Structure.TermIds.Length),
        };
        double worst = 0.0;
        foreach ((string name, Func<BassDraw, double> read) in rows)
        {
            (double ua, double ea) = Summary(a, read);
            (double ub, double eb) = Summary(b, read);
            double z = (ua - ub) / Math.Sqrt(ea * ea + eb * eb);
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{name,-10} scalar {ua:F4} ± {ea:F4}, K = 1 {ub:F4} ± {eb:F4}, z {z:F2}");
        }
        double surface = Enumerable.Range(0, 3).Max(i => Math.Abs(ma[i, 0] - mb[i, 0]));
        output.WriteLine($"mean surfaces within {surface:F4}; worst |z| {worst:F2}");
        Assert.True(worst < 4.0, $"{worst:F2}");
        Assert.True(surface < 0.05, $"{surface:F4}");
    }
}

internal static class ArrayNegate
{
    // −a, a new array.
    internal static double[,] Select2Negate(this double[,] a) => a.Select2((i, j) => -a[i, j]);
}
