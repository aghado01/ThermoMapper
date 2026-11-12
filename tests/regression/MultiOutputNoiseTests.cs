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
/// E34's slice 4, stage 3: several outputs under a noise by row (decision 6), one vᵢ per case shared by its K outputs and a skew γ ∈ ℝᴷ with
/// γ | Σ ~ N(0, s²Σ), Σ integrated in the structure step, the prior apart from the weights and the intercept's level frozen apart from v.
/// <list type="bullet">
/// <item><b>K = 1 is scalar BASS</b> under the t, the Normal–Wald with its tail drawn and the horseshoe: a fit draw for draw, its latents
/// included, and the Rao–Blackwellized mean surface, under level centring and the flat intercept (scalar BASS's moments reference refuses
/// the zero centre; several outputs read it as the levels at 0).</item>
/// <item><b>The target is equivariant in the outputs:</b> Y → YT with γ → Tᵀγ shifts it by −(bN − f + 1)·ln|det T| at fixed v, f = 1 under
/// the flat intercept, which integrates one level, and 0 otherwise; the response-scaled Σ prior, the frozen levels and γ's pseudo-row move with
/// the map.</item>
/// <item><b>The γ·v move on a frozen state</b> keeps its orbit's law, read on a grid by the engine's own target, and its forward and reverse
/// acceptances cancel.</item>
/// <item><b>Geweke from exact starts</b> (qualification) under the t, the skewed t at b = 1 and at a tilted b, the Normal–Wald with its tail
/// and skew drawn, the skewed horseshoe with its rescaling, and the skewed t under Huang and Wand's Σ prior with its rates drawn.</item>
/// <item><b>Recovery from outside the model</b> (qualification): skewed heavy-tailed noise shared across outputs, surfaces outside the basis.</item>
/// <item><b>The readers and refusals:</b> prediction draws the noise by row; the matrix-t conditionals refuse a noise; the asymmetric
/// Laplace, another reference, interweaving, a skew prior off zero and a mixture built for another K are refused.</item>
/// </list>
/// </summary>
public sealed class MultiOutputNoiseTests(ITestOutputHelper output)
{
    private const int Seed = 2026_1111;

    // n cases on two inputs and K outputs from outside the model: smooth surfaces, a case's v from t_ν's mixing IG(ν/2, ν/2) shared by its
    // outputs, the skew `gamma` per unit v, and N(0, vΣ) noise, Σ's correlations ±rho alternating; the surfaces alone as the truth, and v.
    internal static (BassPredictor[] Predictors, double[,] Y, double[,] Truth, double[] V) Data(int n, int k, int seed, double nu = 4.0, double[]? gamma = null,
                                                                                                double rho = 0.4, double scale = 0.04)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x0 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())], x1 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        var sigma = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++) sigma[a, b] = scale * Math.Sqrt((1 + 0.5 * a) * (1 + 0.5 * b)) * (a == b ? 1.0 : rho * (((a + b) % 2) == 0 ? 1 : -1));
        double[,] root = ReferenceAlgebra.Cholesky(sigma)!;
        var mixing = InverseGamma.ShapeScale(0.5 * nu, 0.5 * nu);
        var y = new double[n, k];
        var truth = new double[n, k];
        var scales = new double[n];
        var z = new double[k];
        for (int i = 0; i < n; i++)
        {
            double[] f = [Math.Sin(3 * x0[i]) + 0.5 * x1[i], 2 * Math.Exp(-x1[i]) * x0[i], Math.Cos(2 * x1[i]) + x0[i] * x1[i]];
            double v = scales[i] = mixing.Sample(rng);
            for (int l = 0; l < k; l++) z[l] = Normal.SampleStandard(rng);
            for (int l = 0; l < k; l++)
            {
                double e = 0.0;
                for (int q = 0; q <= l; q++) e += root[l, q] * z[q];
                truth[i, l] = f[l % 3];
                y[i, l] = f[l % 3] + v * (gamma?[l] ?? 0.0) + Math.Sqrt(v) * e;
            }
        }
        return ([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], y, truth, scales);
    }

    private static string[] Names(int k) => [.. Enumerable.Range(0, k).Select(l => $"y{l}")];

    private static BassData Several(BassPredictor[] p, double[,] y) => new(p, new MultiOutputResponse(y, Names(y.GetLength(1))));

    // The laws by name: plain, with the Normal–Wald's tail drawn, and skewed with γ drawn.
    internal static BassNoise Law(string name) => name switch
    {
        "t" => BassNoise.StudentT(4.0),
        "t8" => BassNoise.StudentT(8.0),
        "nw" => BassNoise.NormalWald(gamma: new BassNoisePrior(2.0, 1.0)),
        "horseshoe" => BassNoise.Horseshoe(),
        "skewed-t" => BassNoise.SkewedT(6.0, new BassNoisePrior(0.0, 1.0)),
        "skewed-nw" => BassNoise.NormalWald(gamma: new BassNoisePrior(2.0, 1.0), skew: new BassNoisePrior(0.0, 0.7)),
        "skewed-horseshoe" => BassNoise.Horseshoe(skew: new BassNoisePrior(0.0, 1.0)) with { SkewRescale = true },
        _ => throw new ArgumentOutOfRangeException(nameof(name)),
    };

    // ── K = 1 is scalar BASS ──────────────────────────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData("t", GPriorCentre.Mean)]
    [InlineData("t", GPriorCentre.Flat)]
    [InlineData("nw", GPriorCentre.Mean)]
    [InlineData("horseshoe", GPriorCentre.Mean)]
    public void KOne_UnderANoise_IsScalarBass_DrawForDraw(string law, GPriorCentre centre)
    {
        (BassPredictor[] p, double[,] y, _, _) = Data(60, 1, Seed);
        BassConfig config = BassRunTests.Config(chains: 2, burn: 100, samples: 60, seed: 17, maxBasis: 6) with { Intercept = centre, Noise = Law(law) };
        BassData one = new(p, new ScalarResponse([.. Enumerable.Range(0, 60).Select(i => y[i, 0])])), several = Several(p, y);
        BassFit a = Bass.Fit(config, one), b = Bass.Fit(config, several);
        double worst = 0.0;
        int draws = 0;
        for (int c = 0; c < a.Draws.Count; c++)
            for (int t = 0; t < a.Draws[c].Count; t++, draws++)
            {
                BassDraw u = a.Draws[c][t], w = b.Draws[c][t];
                Assert.Equal(u.Structure.TermIds.Select(id => a.Terms[id]), w.Structure.TermIds.Select(id => b.Terms[id]));
                worst = Math.Max(worst, Relative(u.Tau, w.Tau));
                worst = Math.Max(worst, Relative(u.Sigma2, w.Sigma![0, 0]));
                for (int j = 0; j < u.Coefficients.Length; j++) worst = Math.Max(worst, Math.Abs(u.Coefficients[j] - w.OutputCoefficients![j, 0]) / (1 + Math.Abs(u.Coefficients[j])));
                ReadOnlySpan<double> vu = u.Latents!.Values, vw = w.Latents!.Values;
                for (int i = 0; i < vu.Length; i++) worst = Math.Max(worst, Relative(vu[i], vw[i]));
                if (u.Latents.Gamma is double tail) worst = Math.Max(worst, Relative(tail, w.Latents.Gamma!.Value));
            }
        var cases = new PredictionCases([new ContinuousColumn("x0", [0.2, 0.7]), new ContinuousColumn("x1", [0.4, 0.9])]);
        double[,] ma = a.Predict(one, cases, new BassPredictOptions { MeanSurface = true }).MeanSurface!;
        double[,] mb = BassMultiOutput.Predict(b, several, cases, new BassPredictOptions { MeanSurface = true }).MeanSurface!;
        double surface = Math.Max(Math.Abs(ma[0, 0] - mb[0, 0]), Math.Abs(ma[1, 0] - mb[1, 0]));
        output.WriteLine($"K = 1 under {law}, {centre}: {draws} draws to {worst:G3}, the mean surface to {surface:G3}");
        Assert.True(worst < 1e-9, $"the draws differ by {worst:G3}");
        Assert.True(surface < 1e-9, $"the mean surfaces differ by {surface:G3}");
    }

    // ── The target's equivariance in the outputs ──────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(GPriorCentre.Mean, 1.0)]
    [InlineData(GPriorCentre.Zero, 0.4)]
    [InlineData(GPriorCentre.Flat, 1.0)]
    public void Target_UnderANoise_IsEquivariantInTheOutputs(GPriorCentre centre, double b)
    {
        const int n = 40, k = 3;
        (BassPredictor[] p, double[,] y, _, _) = Data(n, k, Seed + 10, gamma: [0.1, -0.05, 0.08]);
        var rng = new Xoshiro256PlusPlus(Seed + 11);
        var map = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int c = 0; c < k; c++) map[a, c] = (a == c ? 1.0 : 0.0) + 0.35 * (2.0 * rng.NextDouble() - 1.0);
        double logDet = Math.Log(Math.Abs(Determinant3(map)));
        double[,] mapped = ReferenceAlgebra.Multiply(y, map);
        BassConfig config = BassRunTests.Config(chains: 1, maxBasis: 6) with { Intercept = centre, Noise = BassNoise.SkewedT(4.0, new BassNoisePrior(0.0, 0.8)) };
        BassModel model = BassModel.Create(config, Several(p, y)), image = BassModel.Create(config, Several(p, mapped));
        double[] v = [.. Enumerable.Range(0, n).Select(_ => 0.3 + 2.5 * rng.NextDouble())];
        double[] gamma = [0.12, -0.07, 0.05];
        double[] moved = [.. Enumerable.Range(0, k).Select(c => Enumerable.Range(0, k).Sum(a => gamma[a] * map[a, c]))];
        var sampler = new BassPriorSampler(model);
        var temperature = new Temperature(b, TemperingPath.Likelihood);
        double Value(BassModel m, BassTerm[] terms, double tau, double[] skew)
        {
            var mixture = (VarianceMeanMixture)m.Observation;
            WorkingResponse reading = mixture.At(m.Working, new ObservationLatents(v, 0.0, null, null) { OutputSkew = skew });
            BassDesign d = m.NewDesign();
            Assert.Null(d.Load(terms, reading));
            var state = new BassState([.. d.Terms], BassStructureFit.Of(d, 0.0), Math.Log(tau), double.NaN);
            return temperature.LogDensity(m.Target.Evaluate(state));
        }
        double worst = 0.0, skewMoves = 0.0;
        double[] taus = [1e-3, 0.05, 1.0, 20.0];
        for (int c = 0; c < 10; c++)
        {
            BassTerm[] terms = [.. sampler.DrawState(rng).State.Terms];
            double tau = taus[c % taus.Length];
            double at = Value(model, terms, tau, gamma), there = Value(image, terms, tau, moved);
            // The flat intercept integrates one level, whose Lebesgue measure moves by |det T|: the likelihood's shift is −(bN − 1)·ln|det T| there.
            worst = Math.Max(worst, Math.Abs(there - at + (b * n - (centre == GPriorCentre.Flat ? 1.0 : 0.0) + 1.0) * logDet) / (1.0 + Math.Abs(at)));
            skewMoves = Math.Max(skewMoves, Math.Abs(Value(model, terms, tau, [.. gamma.Select(g => 3.0 * g)]) - at));
        }
        output.WriteLine($"{centre} at b = {b}: the target's shift against −(bN − f + 1)·ln|det T| to {worst:G3}; tripling γ moves it by up to {skewMoves:G3}");
        Assert.True(worst < 1e-9, $"{worst:G3}");
        Assert.True(skewMoves > 1e-3, "γ reaches the target");
    }

    // ── The γ·v move on a frozen state ────────────────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(1.0)]
    [InlineData(0.5)]
    public void TheSkewRescale_KeepsItsOrbitsLaw_OnAFrozenState(double b)
    {
        const int k = 2;
        (BassPredictor[] p, double[,] y, _, _) = Data(60, k, Seed + 20, gamma: [0.15, -0.1]);
        BassModel model = BassModel.Create(BassRunTests.Config(chains: 1, maxBasis: 6) with { Noise = Law("skewed-horseshoe") }, Several(p, y));
        var mixture = (VarianceMeanMixture)model.Observation;
        var rng = new Xoshiro256PlusPlus(Seed + 21);
        // A structure and latents to freeze at: a short chain's, the whole sweep under the noise.
        var start = new BassChain(model, [], 0.0, new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        for (int s = 0; s < 300; s++) start.Step();
        var moves = new BassStructureMoves(model, start.Current.Terms, start.Current.LogTau);
        var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(3), new Temperature(b, TemperingPath.Likelihood));
        Assert.Null(BassStructureMoves.Follow(chain, start.Current.Fit.Working!));
        moves.DrawOutputSweepFit(chain.Current, rng, model.Target, b);
        var block = new BassOutputSkewRescaleBlock(mixture);
        WorkingResponse working = moves.DesignAt(chain.Current).Working;
        ObservationLatents latents = working.Latents!;
        double[] v0 = latents.ValueArray, gamma0 = latents.OutputSkew!;
        IMixingLaw mixing = mixture.MixingAt(latents);
        double h = BassOutputSkewRescaleBlock.Scale(v0, k, b, mixing), current = chain.CurrentLogDensity;
        int n = v0.Length;

        // The orbit's law, ln p(t) = ln π(Tₜ(x₀)) + (N − K)·t, π the target with Σ and the coefficients integrated, γ's prior inside it.
        double LogOrbit(double t)
        {
            double[] v = [.. v0.Select(vi => vi * Math.Exp(t))];
            double[] skew = [.. gamma0.Select(g => g * Math.Exp(-t))];
            BassState candidate = moves.ProposeFollow(chain.Current, mixture.At(working, new ObservationLatents(v, 0.0, null, latents.Gamma) { OutputSkew = skew }));
            double target = chain.Temperature.LogDensity(chain.Target.Evaluate(candidate)) - current;
            return target + v.Sum(mixing.LogPrior) - v0.Sum(mixing.LogPrior) + (n - k) * t;
        }
        const int points = 801;
        double width = 10.0 * h;
        var grid = new double[points];
        var logs = new double[points];
        for (int g = 0; g < points; g++) logs[g] = LogOrbit(grid[g] = -width + 2.0 * width * g / (points - 1));
        double top = logs.Max(), mass = 0.0, first = 0.0, second = 0.0;
        for (int g = 0; g < points; g++)
        {
            double w = Math.Exp(logs[g] - top);
            mass += w;
            first += w * grid[g];
            second += w * grid[g] * grid[g];
        }
        double mean = first / mass, meanSquare = second / mass;
        Assert.True(Math.Exp(logs[0] - top) < 1e-10 && Math.Exp(logs[^1] - top) < 1e-10, "the grid holds the orbit's mass");

        // The move's chain from x₀: t = ln(γ₀/γ) on the first output at each step.
        const int steps = 60_000, batches = 40;
        var ts = new double[steps];
        int accepted = 0;
        for (int s = 0; s < steps; s++)
        {
            if (block.Update(chain, rng))
            {
                accepted++;
                Assert.Same(chain.Current, moves.OutputSweepFit!.State);   // Σ and the coefficients drawn afresh at the new state
            }
            ts[s] = Math.Log(gamma0[0] / moves.DesignAt(chain.Current).Working.Latents!.OutputSkew![0]);
        }
        double Z(Func<double, double> f, double expected)
        {
            int size = steps / batches;
            double[] means = [.. Enumerable.Range(0, batches).Select(q => ts.Skip(q * size).Take(size).Average(f))];
            double average = means.Average(), se = Math.Sqrt(means.Sum(x => (x - average) * (x - average)) / (batches - 1) / batches);
            return (average - expected) / se;
        }
        double zMean = Z(t => t, mean), zSquare = Z(t => t * t, meanSquare);
        output.WriteLine($"b = {b}: h {h:F4}, the orbit's sd {Math.Sqrt(meanSquare - mean * mean):F4}; accepted {accepted} of {steps}; " +
                         $"mean t {ts.Average():F5} against {mean:F5} (z {zMean:F2}), mean t² z {zSquare:F2}");
        Assert.True(Math.Abs(zMean) < 4.0 && Math.Abs(zSquare) < 4.0, $"z {zMean:F2}, {zSquare:F2}");
        Assert.InRange(accepted / (double)steps, 0.15, 0.85);

        // From x to x′ and back the acceptances cancel: the reverse's proposal, the Jacobian and the priors are each other's.
        WorkingResponse here = moves.DesignAt(chain.Current).Working;
        double[] v1 = here.Latents!.ValueArray, g1 = here.Latents.OutputSkew!;
        double h1 = BassOutputSkewRescaleBlock.Scale(v1, k, b, mixing), step = 0.7 * h1;
        BassOutputSkewRescaleBlock.Move forward = block.Propose(chain, here, here.Latents, v1, g1, mixing, b, h1, step)!.Value;
        chain.SetState(forward.Candidate, forward.Evaluation);
        WorkingResponse there = moves.DesignAt(chain.Current).Working;
        double[] v2 = there.Latents!.ValueArray;
        BassOutputSkewRescaleBlock.Move back = block.Propose(chain, there, there.Latents, v2, there.Latents.OutputSkew!, mixing, b,
                                                             BassOutputSkewRescaleBlock.Scale(v2, k, b, mixing), -step)!.Value;
        output.WriteLine($"ln A forward {forward.LogAcceptance:R}, back {back.LogAcceptance:R}");
        Assert.Equal(0.0, forward.LogAcceptance + back.LogAcceptance, 8);
    }

    // ── Geweke from exact starts (qualification) ──────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData("t", 1.0)]
    [InlineData("skewed-t", 1.0)]
    [InlineData("skewed-t", 0.5)]
    [InlineData("skewed-nw", 1.0)]
    [InlineData("skewed-horseshoe", 1.0)]
    [InlineData("skewed-t+hw", 1.0)]
    public void Geweke_FromExactStarts_UnderANoise(string law, double b)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int starts = 40_000, sweeps = 10, n = 24, k = 2;
        const double nu0 = 20.0;
        double[,] psi = { { 0.6, 0.15 }, { 0.15, 0.4 } };
        bool hw = law.EndsWith("+hw", StringComparison.Ordinal);   // Huang and Wand's Σ prior at fixed scales, its rates drawn
        BassNoise noise = Law(hw ? law[..^3] : law);
        double c = noise.VarianceScale;
        var rng = new Xoshiro256PlusPlus(Seed + 30);
        (BassPredictor[] p, _, _, _) = Data(n, k, Seed + 31);
        BassConfig config = BassRunTests.Config(chains: 1, maxBasis: 3) with
        {
            Structure = new BassStructureSpec { MaxBasis = 3, MaxInteraction = 2 },
            Intercept = GPriorCentre.Zero,
            OutputCovariance = hw ? OutputCovariancePrior.HuangWand(2, [0.6, 0.5]) : OutputCovariancePrior.InverseWishart(nu0, psi),
            Noise = noise,
        };
        BassModel prior = BassModel.Create(config, Several(p, new double[n, k].Select2((i, l) => Math.Sin(i + l))));
        var sampler = new BassPriorSampler(prior);
        // The tempered joint's exact draw: Σ tilted by |Σ|^(N(1−b)/2), each v by v^(K(1−b)/2); at b = 1 the priors themselves.
        var sigmaTilted = InverseWishart.Of(nu0 - n * (1.0 - b), psi);
        Func<double?, Xoshiro256PlusPlus, double> drawV = b == 1.0
            ? (tail, r) => noise.MixingAt(tail).SamplePrior(r)
            : law == "skewed-t" ? (_, r) => InverseGamma.ShapeScale(3.0 - 0.5 * k * (1.0 - b), 3.0).Sample(r)
            : throw new ArgumentException("A tilted law is written for the skewed t alone.", nameof(law));
        string[] names = ["M", "ln tau", "ln S11", "ln S22", "atanh rho", "mean ln v", "g1", "g2", "tail", "atan f(0)", "atan f(1)", "energy(0)", "ln a1"];
        var deltas = new double[names.Length][];
        for (int f = 0; f < names.Length; f++) deltas[f] = new double[starts];
        int redrawn = 0, refused = 0;
        for (int r = 0; r < starts; r++)
        {
            BassState state = sampler.DrawState(rng).State;   // under Huang and Wand the rates from their prior too
            double[,] sigma = hw ? prior.Target.Outputs!.Rates!.PriorAt(state.OutputRateArray!).Sample(rng) : sigmaTilted.Sample(rng);
            double[,] root = MultiOutputRandomEffectsTests.MatrixTPredictiveRoot(sigma);
            double? tail = noise.GammaPrior switch
            {
                { IsHeld: false } g => Normal.MeanSd(g.Mean, g.Sd).SampleOn(rng, 0.0, double.PositiveInfinity),
                { } g => g.Mean,
                null => null,
            };
            double[] v = [.. Enumerable.Range(0, n).Select(_ => drawV(tail, rng))];
            double[]? gamma = noise.SkewDrawn ? Mixed(root, noise.SkewPrior!.Sd, k, rng) : null;
            BassDesign d = prior.NewDesign();
            Assert.Null(d.Load(state.Terms));   // at v = 1 and c = 1 the prior's Gram is the design's
            double tau = Math.Exp(state.LogTau);
            double[,] f0 = MultiOutputRandomEffectsTests.GewekeSurface(prior, d, root, tau, n, k, rng);
            var y = new double[n, k];
            for (int i = 0; i < n; i++)
            {
                double[] e = Mixed(root, Math.Sqrt(c * v[i] / b), k, rng);
                for (int l = 0; l < k; l++) y[i, l] = f0[i, l] + v[i] * (gamma?[l] ?? 0.0) + e[l];
            }
            // The chain starts at v = 1, where an output of the horseshoe's v ~ 10¹³ leaves Ψ₀ + S positive definite only to rounding. A start
            // whose y exceeds 10⁶ is drawn again: the condition reads y alone, so the sweeps still target p(θ | y) and the test stays exact on
            // the joint conditioned on it.
            if (y.Cast<double>().Any(x => Math.Abs(x) > 1e6))
            {
                redrawn++;
                r--;
                continue;
            }
            double[] before = Read(state, sigma, v, gamma, tail, f0, y);
            BassModel model = BassModel.Create(config, Several(p, y));
            var mixture = (VarianceMeanMixture)model.Observation;
            var moves = new BassStructureMoves(model, state.Terms, state.LogTau, double.NaN, state.OutputRateArray);
            var blockRng = new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue));
            var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target, moves.Start, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)),
                                                           new Temperature(b, TemperingPath.Likelihood));
            // A θ outside F (the horseshoe's v near 10⁻¹² puts one row's weight past the design's numerical floor) is drawn again: the target
            // is the prior on F, so the starts are its draws.
            if (BassStructureMoves.Follow(chain, mixture.At(model.Working, new ObservationLatents(v, 0.0, null, tail) { OutputSkew = gamma })) is not null)
            {
                refused++;
                r--;
                continue;
            }
            List<IBlock<ReversibleJumpChain<BassState>>> blocks = BassSweepSlots.Blocks(model);
            for (int s = 0; s < sweeps; s++)
                foreach (IBlock<ReversibleJumpChain<BassState>> block in blocks) block.Update(chain, blockRng);
            ObservationLatents latents = moves.DesignAt(chain.Current).Working.Latents!;
            BassOutputSweepFit fit = moves.OutputSweepFit!;
            var fitted = new double[n, k].Select2((i, l) => fit.Fitted[i * k + l]);
            double[] after = Read(chain.Current, fit.Sigma, latents.ValueArray, latents.OutputSkew, latents.Gamma, fitted, y);
            for (int f = 0; f < names.Length; f++) deltas[f][r] = after[f] - before[f];
        }
        double worst = 0.0;
        for (int f = 0; f < names.Length; f++)
        {
            double mean = deltas[f].Average(), sd = Math.Sqrt(deltas[f].Sum(x => (x - mean) * (x - mean)) / (starts - 1));
            double z = sd == 0.0 ? (mean == 0.0 ? 0.0 : double.PositiveInfinity) : mean / (sd / Math.Sqrt(starts));
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{names[f],-11} change {mean,10:F5} ± {sd / Math.Sqrt(starts):F5}  z {z,6:F2}");
        }
        output.WriteLine($"{law} at b = {b}, exact starts ({redrawn} drawn again for |y| > 10⁶, {refused} outside F): worst |z| {worst:F2}");
        Assert.True(worst < 4.0, $"{worst:F2}");

        // The functions read: the structure, τ, Σ, v, γ, the tail, the fit at two cases, the first row's tempered energy, and Huang and Wand's
        // first rate.
        double[] Read(BassState s, double[,] sigma, double[] v, double[]? gamma, double? tail, double[,] f, double[,] y)
        {
            double[,] root = MultiOutputRandomEffectsTests.MatrixTPredictiveRoot(sigma);
            var r0 = new double[k];
            for (int l = 0; l < k; l++) r0[l] = y[0, l] - f[0, l] - v[0] * (gamma?[l] ?? 0.0);
            return [s.Count, s.LogTau, Math.Log(sigma[0, 0]), Math.Log(sigma[1, 1]), Math.Atanh(sigma[1, 0] / Math.Sqrt(sigma[0, 0] * sigma[1, 1])),
                    v.Average(Math.Log), gamma?[0] ?? 0.0, gamma?[1] ?? 0.0, tail ?? 0.0, Math.Atan(f[0, 0]), Math.Atan(f[n / 2, 1]),
                    b * BassOutputLatentScaleBlock.Energy(root, r0) / (c * v[0]), s.OutputRateArray?[0] ?? 0.0];
        }
    }

    // L·z·scale for z standard normal: a draw of N(0, scale²·Σ) from Σ's factor L.
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

    // ── Recovery from outside the model (qualification) ───────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Recovery_UnderSkewedHeavyTails_SharedAcrossTheOutputs()
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int n = 200, k = 2;
        double[] gamma = [0.12, -0.08];
        (BassPredictor[] p, double[,] y, double[,] truth, double[] v) = Data(n, k, Seed + 40, nu: 6.0, gamma: gamma, rho: 0.5);
        BassData data = Several(p, y);
        BassConfig gaussian = BassRunTests.Config(chains: 4, burn: 2000, samples: 2000, seed: 41, maxBasis: 12);
        BassFit skewed = Bass.Fit(gaussian with { Noise = Law("skewed-t") }, data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        BassFit plain = Bass.Fit(gaussian, data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        var cases = new PredictionCases([new ContinuousColumn("x0", [.. ((ContinuousPredictor)p[0]).Values.ToArray()]),
                                         new ContinuousColumn("x1", [.. ((ContinuousPredictor)p[1]).Values.ToArray()])]);
        double Rmse(BassFit fit)
        {
            double[,] m = BassMultiOutput.Predict(fit, data, cases, new BassPredictOptions { MeanSurface = true }).MeanSurface!;
            double sum = 0.0;
            for (int i = 0; i < n; i++)
                for (int l = 0; l < k; l++) sum += (m[i, l] - truth[i, l]) * (m[i, l] - truth[i, l]);
            return Math.Sqrt(sum / (n * k));
        }
        double skewedRmse = Rmse(skewed), plainRmse = Rmse(plain);
        BassDraw[] draws = [.. skewed.Draws.SelectMany(c => c)];
        double[][] g = [.. Enumerable.Range(0, k).Select(l => draws.Select(d => d.Latents!.OutputSkew![l]).ToArray())];
        double[] rho = [.. draws.Select(d => d.Sigma![1, 0] / Math.Sqrt(d.Sigma[0, 0] * d.Sigma[1, 1]))];
        output.WriteLine($"surface RMSE: skewed t {skewedRmse:F4}, Gaussian {plainRmse:F4}");
        for (int l = 0; l < k; l++)
        {
            double[] sorted = [.. g[l].Order()];
            output.WriteLine($"γ{l}: mean {g[l].Average():F4}, 95% [{sorted[(int)(0.025 * sorted.Length)]:F4}, {sorted[(int)(0.975 * sorted.Length)]:F4}] against {gamma[l]}");
            Assert.InRange(gamma[l], sorted[(int)(0.005 * sorted.Length)], sorted[(int)(0.995 * sorted.Length)]);
        }
        // Σ's correlation is −0.5 (the generator alternates signs); the posterior centres on this data set's realized one, the standardized
        // Gaussian parts' (yᵢ − fᵢ − vᵢγ)/√vᵢ.
        double[][] e = [.. Enumerable.Range(0, k).Select(l => Enumerable.Range(0, n).Select(i => (y[i, l] - truth[i, l] - v[i] * gamma[l]) / Math.Sqrt(v[i])).ToArray())];
        double m0 = e[0].Average(), m1 = e[1].Average();
        double realized = e[0].Zip(e[1]).Sum(t => (t.First - m0) * (t.Second - m1)) / Math.Sqrt(e[0].Sum(x => (x - m0) * (x - m0)) * e[1].Sum(x => (x - m1) * (x - m1)));
        output.WriteLine($"ρ: mean {rho.Average():F3} against the realized {realized:F3} (Σ's −0.5)");
        Assert.InRange(rho.Average() - realized, -0.1, 0.1);
        Assert.True(skewedRmse < plainRmse, $"{skewedRmse:F4} against {plainRmse:F4}");
    }

    // ── The readers and refusals ──────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Prediction_DrawsTheNoiseByRow_AndTheConditionalsRefuseANoise()
    {
        const int k = 2;
        (BassPredictor[] p, double[,] y, _, _) = Data(80, k, Seed + 50, nu: 8.0);
        BassData data = Several(p, y);
        BassFit fit = Bass.Fit(BassRunTests.Config(chains: 1, burn: 100, samples: 40, seed: 51, maxBasis: 6) with { Noise = Law("t8") }, data);
        var rng = new Xoshiro256PlusPlus(Seed + 52);
        const int m = 500;
        var cases = new PredictionCases([new ContinuousColumn("x0", [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble())]),
                                         new ContinuousColumn("x1", [.. Enumerable.Range(0, m).Select(_ => rng.NextDouble())])]);
        BassOutputPrediction prediction = BassMultiOutput.Predict(fit, data, cases, new BassPredictOptions { Predictive = true, Seed = 53 });
        int draws = prediction.Sigma.Count;
        // E[(y − f)(y − f)ᵀ] = c·E[v]·Σ under t₈'s mixing, E[v] = 8/6, per draw.
        var expected = new double[k, k];
        var seen = new double[k, k];
        for (int d = 0; d < draws; d++)
            for (int a = 0; a < k; a++)
                for (int q = 0; q < k; q++)
                {
                    expected[a, q] += 8.0 / 6.0 * prediction.Sigma[d][a, q] / draws;
                    for (int i = 0; i < m; i++)
                        seen[a, q] += (prediction.Predictive![d, i, a] - prediction.Latent[d, i, a]) * (prediction.Predictive[d, i, q] - prediction.Latent[d, i, q]) / (draws * m);
                }
        double rhoSeen = seen[1, 0] / Math.Sqrt(seen[0, 0] * seen[1, 1]), rhoExpected = expected[1, 0] / Math.Sqrt(expected[0, 0] * expected[1, 1]);
        output.WriteLine($"noise variances {seen[0, 0]:G4}, {seen[1, 1]:G4} against {expected[0, 0]:G4}, {expected[1, 1]:G4}; ρ {rhoSeen:F3} against {rhoExpected:F3}");
        Assert.InRange(seen[0, 0] / expected[0, 0], 0.92, 1.08);
        Assert.InRange(seen[1, 1] / expected[1, 1], 0.92, 1.08);
        Assert.InRange(rhoSeen - rhoExpected, -0.04, 0.04);
        Assert.Throws<NotSupportedException>(() => BassMultiOutput.Conditionals(fit, data, cases));
    }

    [Fact]
    public void Refuses_WhatHasNoKOutputReading()
    {
        (BassPredictor[] p, double[,] y, _, _) = Data(40, 2, Seed + 60);
        BassData data = Several(p, y);
        BassConfig config = BassRunTests.Config(chains: 1, maxBasis: 4);
        Assert.Throws<ArgumentException>(() => BassModel.Create(config with { Noise = BassNoise.AsymmetricLaplace(0.5) }, data));
        Assert.Throws<ArgumentException>(() => BassModel.Create(config with { Noise = BassNoise.StudentT(4.0, reference: BassNoiseReference.Robust) }, data));
        Assert.Throws<ArgumentException>(() => BassModel.Create(config with { Noise = BassNoise.StudentT(4.0) with { Interweave = true } }, data));
        Assert.Throws<ArgumentException>(() => BassModel.Create(config with { Noise = BassNoise.SkewedT(4.0, new BassNoisePrior(0.3, 1.0)) }, data));
        Assert.Throws<ArgumentException>(() => BassModel.Create(config, data, new VarianceMeanMixture(BassNoise.StudentT(4.0))));
        BassModel model = BassModel.Create(config with { Noise = BassNoise.StudentT(4.0) }, data);
        Assert.Equal(2, ((VarianceMeanMixture)model.Observation).Outputs);
    }

    // ── Fixtures ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────

    private static double Determinant3(double[,] a)
        => a[0, 0] * (a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1]) - a[0, 1] * (a[1, 0] * a[2, 2] - a[1, 2] * a[2, 0]) + a[0, 2] * (a[1, 0] * a[2, 1] - a[1, 1] * a[2, 0]);

    private static double Relative(double a, double b) => Math.Abs(a - b) / Math.Max(1e-300, Math.Max(Math.Abs(a), Math.Abs(b)));
}
