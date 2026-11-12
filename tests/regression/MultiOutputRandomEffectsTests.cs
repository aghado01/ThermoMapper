using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E34's slice 4, its first stages: the K-output core on a reading with row weights or a random-effect block's whitening, and random
/// intercepts with one s for every output under Gaussian noise (decision 7): u_g ∈ ℝᴷ ~ N(0, s²Σ), the rows' covariance V ⊗ Σ.
/// <list type="bullet">
/// <item><b>K = 1 is scalar BASS</b> under a held and a drawn s, under each intercept convention: a fit draw for draw, its group effects
/// included, and the Rao–Blackwellized mean surface.</item>
/// <item><b>The target against a dense matrix-t</b> under V ⊗ Σ, V = I + s²ZΠZᵀ written out from its definition, and under W⁻¹ ⊗ Σ for row
/// weights: the rows transformed by a square root of V⁻¹ (or W^½), the g-prior on X_cᵀV⁻¹X_c, the Jacobian −(Kb/2)ln|V| (or (Kb/2)Σ ln wᵢ).</item>
/// <item><b>The follow</b> to a reading at another s is a fresh load there, and a restore after it continues bit for bit.</item>
/// <item><b>The group effects' conditional:</b> their mean and their covariance across the outputs, v_g·Σ, by Monte Carlo.</item>
/// <item><b>The readers:</b> prediction adds an observed level's effects and draws a population case's; the matrix-t conditionals refuse a
/// block.</item>
/// <item><b>Refusals:</b> random slopes, the normalized ridge with a block, and the b = 1 collapse under a noise.</item>
/// <item><b>Geweke from exact starts</b> under a drawn s (qualification), as E21's slice 4 gates its block.</item>
/// </list>
/// </summary>
public sealed class MultiOutputRandomEffectsTests(ITestOutputHelper output)
{
    private const int Seed = 2026_1111;

    // n cases on two inputs, a grouping of `levels` levels and K outputs: surfaces outside the model, u_g ~ N(0, s²Σ) and noise N(0, Σ).
    internal static (BassPredictor[] Predictors, CategoricalPredictor Grouping, double[,] Y) Data(int n, int k, int levels, int seed, double s = 0.7)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x0 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())], x1 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        CategoricalPredictor g = RandomInterceptTests.Grouping(n, levels, seed + 1);
        var sigma = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++) sigma[a, b] = 0.04 * Math.Sqrt((1 + 0.5 * a) * (1 + 0.5 * b)) * (a == b ? 1.0 : 0.4 * (((a + b) % 2) == 0 ? 1 : -1));
        double[,] root = ReferenceAlgebra.Cholesky(sigma)!;
        var effects = new double[levels, k];
        for (int h = 0; h < levels; h++) Mix(root, rng, k, s, (l, v) => effects[h, l] = v);
        var y = new double[n, k];
        int[] codes = g.Codes.ToArray();
        for (int i = 0; i < n; i++)
        {
            double[] f = [Math.Sin(3 * x0[i]) + 0.5 * x1[i], 2 * Math.Exp(-x1[i]) * x0[i], Math.Cos(2 * x1[i]) + x0[i] * x1[i]];
            int at = i;
            Mix(root, rng, k, 1.0, (l, v) => y[at, l] = f[l % 3] + effects[codes[at], l] + v);
        }
        return ([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], g, y);

        static void Mix(double[,] root, Xoshiro256PlusPlus rng, int k, double scale, Action<int, double> put)
        {
            var z = new double[k];
            for (int l = 0; l < k; l++) z[l] = Normal.SampleStandard(rng);
            for (int l = 0; l < k; l++)
            {
                double v = 0.0;
                for (int q = 0; q <= l; q++) v += root[l, q] * z[q];
                put(l, scale * v);
            }
        }
    }

    private static string[] Names(int k) => [.. Enumerable.Range(0, k).Select(l => $"y{l}")];

    private static BassData Several(BassPredictor[] p, CategoricalPredictor g, double[,] y) => new(p, new MultiOutputResponse(y, Names(y.GetLength(1))), [g]);

    private static BassData One(BassPredictor[] p, CategoricalPredictor g, double[,] y, int l)
        => new(p, new ScalarResponse([.. Enumerable.Range(0, y.GetLength(0)).Select(i => y[i, l])]), [g]);

    private static BassConfig Config(double s, bool drawn, GPriorCentre centre = GPriorCentre.Mean, int samples = 60, int seed = 17)
        => BassRunTests.Config(chains: 2, burn: 100, samples: samples, seed: seed, maxBasis: 6) with
        {
            Intercept = centre,
            RandomIntercepts = [new BassRandomIntercept("g", s, ScalePrior: drawn ? new BassGroupScalePrior(1.0) : null)],
        };

    // ── K = 1 is scalar BASS ──────────────────────────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(false, GPriorCentre.Mean)]
    [InlineData(false, GPriorCentre.Zero)]
    [InlineData(false, GPriorCentre.Flat)]
    [InlineData(true, GPriorCentre.Mean)]
    public void KOne_UnderABlock_IsScalarBass_DrawForDraw(bool drawn, GPriorCentre centre)
    {
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = Data(60, 1, 8, Seed);
        BassConfig config = Config(0.8, drawn, centre);
        BassData one = One(p, g, y, 0), several = Several(p, g, y);
        BassFit a = Bass.Fit(config, one), b = Bass.Fit(config, several);
        double worst = 0.0;
        int draws = 0;
        for (int c = 0; c < a.Draws.Count; c++)
            for (int t = 0; t < a.Draws[c].Count; t++, draws++)
            {
                BassDraw u = a.Draws[c][t], v = b.Draws[c][t];
                Assert.Equal(u.Structure.TermIds.Select(id => a.Terms[id]), v.Structure.TermIds.Select(id => b.Terms[id]));
                worst = Math.Max(worst, Relative(u.Tau, v.Tau));
                worst = Math.Max(worst, Relative(u.Sigma2, v.Sigma![0, 0]));
                for (int j = 0; j < u.Coefficients.Length; j++) worst = Math.Max(worst, Math.Abs(u.Coefficients[j] - v.OutputCoefficients![j, 0]) / (1 + Math.Abs(u.Coefficients[j])));
                BassGroupDraw gu = u.Groups![0], gv = v.Groups![0];
                worst = Math.Max(worst, Relative(gu.Scale, gv.Scale));
                for (int h = 0; h < gu.Effects.Length; h++) worst = Math.Max(worst, Math.Abs(gu.Effects[h] - gv.OutputEffects![0][h]) / (1 + Math.Abs(gu.Effects[h])));
            }
        var cases = new PredictionCases([new ContinuousColumn("x0", [0.2, 0.7]), new ContinuousColumn("x1", [0.4, 0.9]), new CategoricalColumn("g", ["L1", "L3"])]);
        double[,] ma = a.Predict(one, cases, new BassPredictOptions { MeanSurface = true }).MeanSurface!;
        double[,] mb = BassMultiOutput.Predict(b, several, cases, new BassPredictOptions { MeanSurface = true }).MeanSurface!;
        double surface = Math.Max(Math.Abs(ma[0, 0] - mb[0, 0]), Math.Abs(ma[1, 0] - mb[1, 0]));
        output.WriteLine($"K = 1 under a {(drawn ? "drawn" : "held")} s, {centre}: {draws} draws to {worst:G3}, the mean surface to {surface:G3}");
        Assert.True(worst < 1e-9, $"the draws differ by {worst:G3}");
        Assert.True(surface < 1e-9, $"the mean surfaces differ by {surface:G3}");
    }

    // ── The target against dense references ───────────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(2, GPriorCentre.Mean)]
    [InlineData(3, GPriorCentre.Mean)]
    [InlineData(3, GPriorCentre.Zero)]
    public void Target_UnderABlock_IsTheDenseMatrixT(int k, GPriorCentre centre)
    {
        const double s = 0.8;
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = Data(40, k, 6, Seed + k);
        BassModel model = BassModel.Create(Config(s, false, centre), Several(p, g, y));
        int n = y.GetLength(0);
        double[,] v = DenseV(g, s);
        double[,] lower = ReferenceAlgebra.Cholesky(v)!;
        double[,] root = ReferenceAlgebra.LowerInverse(lower);   // L⁻¹: L⁻¹V L⁻ᵀ = I
        double logDetV = ReferenceAlgebra.LogDeterminant(lower);
        double worst = 0.0;
        foreach ((List<BassTerm> terms, double tau, double b) in Cases(model, Seed + 10 * k))
        {
            BassDesign design = model.NewDesign();
            if (design.Load(terms) is not null) continue;
            double engine = model.Target.OutputsMarginal(BassStructureFit.Of(design, 0.0), tau, b).LogMarginal(model.Resolved.Outputs!.Prior);
            double[,] xc = CentredColumns(model, terms, null);
            double[] means = [.. Enumerable.Range(0, k).Select(l => Enumerable.Range(0, n).Average(i => y[i, l]))];
            MultiOutputCase reference = Transformed(y, xc, root, n, means, centre, tau, b);
            double dense = MatrixTReference.LogMarginal(MatrixTReference.Collapse(reference), Prior(model)) - b * 0.5 * k * logDetV;
            worst = Math.Max(worst, Math.Abs(engine - dense) / (1.0 + Math.Abs(dense)));
        }
        output.WriteLine($"K = {k}, {centre}, under V ⊗ Σ: the target against the dense matrix-t to {worst:G3}");
        Assert.True(worst < 1e-10, $"{worst:G3}");
    }

    [Theory]
    [InlineData(2, GPriorCentre.Mean)]
    [InlineData(3, GPriorCentre.Zero)]
    public void Target_UnderRowWeights_IsTheDenseMatrixT(int k, GPriorCentre centre)
    {
        (BassPredictor[] p, _, double[,] y) = Data(40, k, 6, Seed + 20 + k);
        BassData data = new(p, new MultiOutputResponse(y, Names(k)));
        BassModel model = BassModel.Create(BassRunTests.Config(chains: 1, maxBasis: 6) with { Intercept = centre }, data);
        int n = y.GetLength(0);
        var rng = new Xoshiro256PlusPlus(Seed + 30);
        double[] w = [.. Enumerable.Range(0, n).Select(_ => 0.3 + 2.0 * rng.NextDouble())];
        var reading = new WorkingResponse(data.Response, weights: w, priorGram: WeightedPriorGram.Weighted);
        var root = new double[n, n];
        for (int i = 0; i < n; i++) root[i, i] = Math.Sqrt(w[i]);
        double mass = w.Sum(), logWeights = w.Sum(Math.Log);
        double[] means = [.. Enumerable.Range(0, k).Select(l => Enumerable.Range(0, n).Sum(i => w[i] * y[i, l]) / mass)];
        double worst = 0.0;
        foreach ((List<BassTerm> terms, double tau, double b) in Cases(model, Seed + 40 + k))
        {
            BassDesign design = model.NewDesign(reading);
            if (design.Load(terms) is not null) continue;
            double engine = model.Target.OutputsMarginal(BassStructureFit.Of(design, 0.0), tau, b).LogMarginal(model.Resolved.Outputs!.Prior);
            double[,] xc = CentredColumns(model, terms, w);
            MultiOutputCase reference = Transformed(y, xc, root, mass, means, centre, tau, b);
            double dense = MatrixTReference.LogMarginal(MatrixTReference.Collapse(reference), Prior(model)) + b * 0.5 * k * logWeights;
            worst = Math.Max(worst, Math.Abs(engine - dense) / (1.0 + Math.Abs(dense)));
        }
        output.WriteLine($"K = {k}, {centre}, under W⁻¹ ⊗ Σ: the target against the dense matrix-t to {worst:G3}");
        Assert.True(worst < 1e-10, $"{worst:G3}");
    }

    // ── The follow, and a restore after it ────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Follow_ToAnotherScale_IsAFreshLoad_AndARestoreContinuesBitForBit()
    {
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = Data(50, 3, 7, Seed + 50);
        BassModel model = BassModel.Create(Config(0.8, true), Several(p, g, y));
        List<BassTerm> terms = Cases(model, Seed + 51).First(c => model.NewDesign().Load(c.Terms) is null && c.Terms.Count > 0).Terms;
        BassDesign design = model.NewDesign();
        Assert.Null(design.Load(terms));
        WorkingResponse there = model.RandomEffects!.At(0, 1.9).Rewhiten(model.BaseReading);
        DesignProposal followed = design.ProposeFollow(there);
        BassDesign fresh = model.NewDesign();
        Assert.Null(fresh.Load(terms, there));
        double gap = Math.Max(MaxGap(followed.Outputs!.Residual, fresh.OutputForms!.Residual), MaxGap(followed.Outputs.Fitted, fresh.OutputForms.Fitted));
        output.WriteLine($"the follow to s = 1.9 against a fresh load there: {gap:G3}");
        Assert.True(gap < 1e-9, $"{gap:G3}");
        design.Accept();
        Assert.Equal(0.0, MaxGap(design.OutputForms!.Residual, followed.Outputs.Residual));

        BassDesign restored = model.NewDesign();
        restored.Restore(design.Snapshot());
        Assert.Equal(0.0, MaxGap(design.OutputForms!.Residual, restored.OutputForms!.Residual));
        BassTerm candidate = Cases(model, Seed + 52).SelectMany(c => c.Terms).First(t => !terms.Contains(t));
        DesignProposal pa = design.ProposeBirth(candidate), pb = restored.ProposeBirth(candidate);
        Assert.Equal(pa.Refusal, pb.Refusal);
        if (pa.Refusal is null) Assert.Equal(0.0, MaxGap(pa.Outputs!.Residual, pb.Outputs!.Residual));
    }

    // ── The group effects' conditional ────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void GroupEffects_FollowTheirConditional_AcrossTheOutputs()
    {
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = Data(40, 2, 5, Seed + 60);
        BassModel model = BassModel.Create(Config(0.8, false), Several(p, g, y));
        List<BassTerm> terms = Cases(model, Seed + 61).First(c => model.NewDesign().Load(c.Terms) is null && c.Terms.Count > 0).Terms;
        RandomEffectsBlock block = model.RandomEffects!;
        var coefficients = new double[terms.Count + 1, 2];
        for (int j = 0; j <= terms.Count; j++) { coefficients[j, 0] = 0.1 * j; coefficients[j, 1] = -0.05 * j + 0.2; }
        double[,] sigma = { { 0.05, 0.02 }, { 0.02, 0.03 } };
        double[,] root = MatrixTPredictiveRoot(sigma);
        double[][][] mean = block.MeanOutputEffects(model, terms, coefficients);
        const int draws = 40_000;
        var rng = new Xoshiro256PlusPlus(Seed + 62);
        int levels = mean[0][0].Length;
        var sum = new double[levels, 2];
        var cross = new double[levels, 3];   // u0², u0u1, u1² per level about the conditional mean
        for (int t = 0; t < draws; t++)
        {
            BassGroupDraw drawn = block.DrawOutputGroups(model, terms, coefficients, root, rng)[0];
            for (int h = 0; h < levels; h++)
            {
                double a = drawn.OutputEffects![0][h] - mean[0][0][h], b = drawn.OutputEffects[1][h] - mean[1][0][h];
                sum[h, 0] += a; sum[h, 1] += b;
                cross[h, 0] += a * a; cross[h, 1] += a * b; cross[h, 2] += b * b;
            }
        }
        // The scalar block's spread at σ² = 1, v_g, from its own draws on a zero residual.
        var spread = new double[levels];
        var single = new Xoshiro256PlusPlus(Seed + 63);
        for (int t = 0; t < draws; t++)
        {
            BassGroupDraw drawn = block.DrawOutputGroups(model, terms, new double[terms.Count + 1, 2], new double[,] { { 1, 0 }, { 0, 1 } }, single)[0];
            double[][][] zero = block.MeanOutputEffects(model, terms, new double[terms.Count + 1, 2]);
            for (int h = 0; h < levels; h++) spread[h] += (drawn.OutputEffects![0][h] - zero[0][0][h]) * (drawn.OutputEffects[0][h] - zero[0][0][h]) / draws;
        }
        double worst = 0.0;
        for (int h = 0; h < levels; h++)
        {
            double v = spread[h];
            worst = Math.Max(worst, Math.Abs(sum[h, 0] / draws) / Math.Sqrt(v * sigma[0, 0] / draws));
            worst = Math.Max(worst, Math.Abs(sum[h, 1] / draws) / Math.Sqrt(v * sigma[1, 1] / draws));
            double[] expected = [v * sigma[0, 0], v * sigma[0, 1], v * sigma[1, 1]];
            for (int e = 0; e < 3; e++)
            {
                double se = Math.Sqrt(2.0 / draws) * v * Math.Max(sigma[0, 0], sigma[1, 1]) * 1.5;
                worst = Math.Max(worst, Math.Abs(cross[h, e] / draws - expected[e]) / se);
            }
        }
        output.WriteLine($"group effects: worst standardized gap {worst:F2} over {levels} levels");
        Assert.True(worst < 4.5, $"{worst:F2}");
    }

    // ── The readers ───────────────────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Prediction_AddsAnObservedLevelsEffects_AndTheConditionalsRefuseABlock()
    {
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = Data(60, 2, 6, Seed + 70);
        BassData data = Several(p, g, y);
        BassFit fit = Bass.Fit(Config(0.8, false, samples: 30), data);
        var at = new PredictionCases([new ContinuousColumn("x0", [0.3, 0.3]), new ContinuousColumn("x1", [0.6, 0.6]), new CategoricalColumn("g", ["L2", "L9"])]);
        BassOutputPrediction pr = BassMultiOutput.Predict(fit, data, at, new BassPredictOptions { UnseenLevels = UnseenLevel.OutsideEverySubset, Predictive = true });
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        for (int d = 0; d < pr.Draws.Count; d++)
        {
            double[][] effects = pooled[pr.Draws[d]].Groups![0].OutputEffects!;
            for (int l = 0; l < 2; l++) Assert.Equal(effects[l][2], pr.Latent[d, 0, l] - pr.Latent[d, 1, l], 10);
        }
        Assert.Throws<NotSupportedException>(() => BassMultiOutput.Conditionals(fit, data, at, new BassPredictOptions { UnseenLevels = UnseenLevel.OutsideEverySubset }));
    }

    // ── Refusals ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Refuses_Slopes_TheRidgeUnderABlock_AndANoise()
    {
        (BassPredictor[] p, CategoricalPredictor g, double[,] y) = Data(40, 2, 5, Seed + 80);
        BassData data = Several(p, g, y);
        BassConfig block = Config(0.8, false);
        Assert.Throws<ArgumentException>(() => BassModel.Create(block with { RandomIntercepts = null, RandomSlopes = [new BassRandomSlope("g", "x0", RandomSlopeScale.Unit, 0.5)] }, data));
        Assert.Throws<ArgumentException>(() => BassModel.Create(block with { Coefficients = new BassCoefficientPrior(new RidgeForm(), TauPrior.UnitInformationGamma()) }, data));
        // Under a noise the block is a latent on every rung (stage 4); the b = 1 collapse would need K outputs' weighted whitening.
        Assert.Throws<ArgumentException>(() => BassModel.Create(block with { Noise = BassNoise.StudentT(4.0), Sampling = block.Sampling with { CollapseEffectsInMoves = true } }, data));
    }

    // ── Geweke from exact starts (qualification) ──────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Geweke_FromExactStarts_UnderADrawnScale()
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int starts = 40_000, sweeps = 10, n = 24, k = 2, levels = 6;
        double[,] psi = { { 0.06, 0.015 }, { 0.015, 0.04 } };
        var rng = new Xoshiro256PlusPlus(Seed + 90);
        (BassPredictor[] p, CategoricalPredictor g, _) = Data(n, k, levels, Seed + 91);
        BassConfig Configured(double s) => BassRunTests.Config(chains: 1, maxBasis: 3) with
        {
            Structure = new BassStructureSpec { MaxBasis = 3, MaxInteraction = 2 },
            Intercept = GPriorCentre.Zero,
            OutputCovariance = OutputCovariancePrior.InverseWishart(5.0, psi),
            RandomIntercepts = [new BassRandomIntercept("g", s, ScalePrior: new BassGroupScalePrior(1.0))],
        };
        BassModel prior = BassModel.Create(Configured(1.0), Several(p, g, new double[n, k].Select2((i, l) => Math.Sin(i + l))));
        var sampler = new BassPriorSampler(prior);
        InverseWishart sigmaPrior = prior.Resolved.Outputs!.Prior;
        string[] names = ["M", "ln tau", "ln s", "ln S11", "ln S22", "atanh rho", "atan f0(0)", "atan u0(0)"];
        var deltas = new double[names.Length][];
        for (int f = 0; f < names.Length; f++) deltas[f] = new double[starts];
        int[] codes = g.Codes.ToArray();
        for (int t = 0; t < starts; t++)
        {
            (BassState state, _) = sampler.DrawPrior(rng);
            WorkingResponse reading = sampler.LastDesign().Response!.Working;   // the reading whitened at the drawn s
            double s = reading.Whitening!.Scales.Span[0];
            double[,] sigma = sigmaPrior.Sample(rng);
            double[,] root = MatrixTPredictiveRoot(sigma);
            BassDesign d = prior.NewDesign();
            Assert.Null(d.Load(state.Terms, reading));   // A's prior reads the Gram whitened at s
            int m = d.Count;
            double tau = Math.Exp(state.LogTau);
            double[,] f0 = GewekeSurface(prior, d, root, tau, n, k, rng);
            // u_g ~ N(0, s²Σ) on the observation-centred levels: Π z, Π the projector onto nᵀu = 0.
            double[,] u = CentredEffects(codes, levels, k, s, root, rng);
            var y = new double[n, k];
            for (int i = 0; i < n; i++)
            {
                var z = new double[k];
                for (int l = 0; l < k; l++) z[l] = Normal.SampleStandard(rng);
                for (int l = 0; l < k; l++)
                {
                    double e = 0.0;
                    for (int q = 0; q <= l; q++) e += root[l, q] * z[q];
                    y[i, l] = f0[i, l] + u[codes[i], l] + e;
                }
            }
            double[] before = [m, state.LogTau, Math.Log(s), Math.Log(sigma[0, 0]), Math.Log(sigma[1, 1]), Math.Atanh(sigma[1, 0] / Math.Sqrt(sigma[0, 0] * sigma[1, 1])),
                               Math.Atan(f0[0, 0]), Math.Atan(u[codes[0], 0])];
            BassModel model = BassModel.Create(Configured(s), Several(p, g, y));
            var chain = new BassChain(model, state.Terms, state.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
            BassState now = state;
            for (int w = 0; w < sweeps; w++) now = chain.Step();
            BassDesign design = chain.Structure.DesignAt(now);
            RandomEffectsBlock at = design.Working.Whitening!;
            (double[,] drawnSigma, double[,] centred) = BassOutputDraws.Sample(model.Target, design, now.Fit, Math.Exp(now.LogTau), rng);
            double[,] fitted = Fitted(model, design, centred);
            BassGroupDraw effects = at.DrawOutputGroups(model, design.Terms, Uncentred(model, design, centred), MatrixTPredictiveRoot(drawnSigma), rng)[0];
            double[] after = [now.Count, now.LogTau, Math.Log(at.Scales.Span[0]), Math.Log(drawnSigma[0, 0]), Math.Log(drawnSigma[1, 1]),
                              Math.Atanh(drawnSigma[1, 0] / Math.Sqrt(drawnSigma[0, 0] * drawnSigma[1, 1])), Math.Atan(fitted[0, 0]), Math.Atan(effects.OutputEffects![0][codes[0]])];
            for (int f = 0; f < names.Length; f++) deltas[f][t] = after[f] - before[f];
        }
        double worst = 0.0;
        for (int f = 0; f < names.Length; f++)
        {
            double mean = deltas[f].Average(), sd = Math.Sqrt(deltas[f].Sum(x => (x - mean) * (x - mean)) / (starts - 1));
            double z = mean / (sd / Math.Sqrt(starts));
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{names[f],-11} change {mean,10:F5} ± {sd / Math.Sqrt(starts):F5}  z {z,6:F2}");
        }
        output.WriteLine($"exact starts: worst |z| {worst:F2}");
        Assert.True(worst < 4.0, $"{worst:F2}");
    }

    // ── Fixtures ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────

    // Structures from the prior, each with a τ and a b.
    private static IEnumerable<(List<BassTerm> Terms, double Tau, double B)> Cases(BassModel model, int seed)
    {
        var sampler = new BassPriorSampler(model);
        var rng = new Xoshiro256PlusPlus(seed);
        double[] taus = [1e-3, 0.05, 1.0, 20.0], bs = [1.0, 0.4];
        for (int c = 0; c < 12; c++)
            yield return ([.. sampler.DrawState(rng).State.Terms], taus[c % taus.Length], bs[c % bs.Length]);
    }

    // V = I + s²ZΠZᵀ over the cases, Π the projector onto nᵀu = 0 (observation centring), written out.
    internal static double[,] DenseV(CategoricalPredictor grouping, double s)
    {
        int[] codes = grouping.Codes.ToArray();
        int n = codes.Length, l = codes.Max() + 1;
        var counts = new double[l];
        foreach (int c in codes) counts[c]++;
        double nn = counts.Sum(c => c * c);
        var v = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                v[i, j] = (i == j ? 1.0 : 0.0) + s * s * ((codes[i] == codes[j] ? 1.0 : 0.0) - counts[codes[i]] * counts[codes[j]] / nn);
        return v;
    }

    // The terms' raw case columns centred over the cases, under W's weighted means when given.
    private static double[,] CentredColumns(BassModel model, IReadOnlyList<BassTerm> terms, double[]? w)
    {
        int n = model.Data.Response.Cases, m = terms.Count;
        var x = new double[n, m];
        var column = new double[n];
        for (int j = 0; j < m; j++)
        {
            model.Basis.CaseColumn(terms[j], model.Inputs.Cases, column);
            double mean = w is null ? column.Average() : column.Select((v, i) => w[i] * v).Sum() / w.Sum();
            for (int i = 0; i < n; i++) x[i, j] = column[i] - mean;
        }
        return x;
    }

    // The g-form on rows transformed by `root` (a square root of the rows' precision): Y' = R·Y, X' = R·[1, X_c] (R·X_c alone under the flat
    // intercept), Λ = τ·blockdiag(mass, (R X_c)ᵀ(R X_c)), the intercept centred on the means (level centring) or at 0.
    private static MultiOutputCase Transformed(double[,] y, double[,] xc, double[,] root, double mass, double[] means, GPriorCentre centre, double tau, double b)
    {
        int n = y.GetLength(0), k = y.GetLength(1), m = xc.GetLength(1);
        double[,] ry = ReferenceAlgebra.Multiply(root, y), rx = ReferenceAlgebra.Multiply(root, xc);
        double[,] gram = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(rx), rx);
        var ones = new double[n, 1];
        for (int i = 0; i < n; i++) ones[i, 0] = 1.0;
        double[,] r1 = ReferenceAlgebra.Multiply(root, ones);
        var x = new double[n, m + 1];
        var lambda = new double[m + 1, m + 1];
        lambda[0, 0] = tau * mass;
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = r1[i, 0];
            for (int j = 0; j < m; j++) x[i, j + 1] = rx[i, j];
        }
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++) lambda[i + 1, j + 1] = tau * gram[i, j];
        if (centre == GPriorCentre.Zero) return new MultiOutputCase(ry, x, lambda, null, b);
        var a0 = new double[m + 1, k];
        for (int l = 0; l < k; l++) a0[0, l] = means[l];
        return new MultiOutputCase(ry, x, lambda, a0, b);
    }

    private static ReferenceSigmaPrior Prior(BassModel model)
    {
        InverseWishart iw = model.Resolved.Outputs!.Prior;
        return new ReferenceSigmaPrior(iw.DegreesOfFreedom, iw.Scale);
    }

    internal static double[,] MatrixTPredictiveRoot(double[,] sigma) => ReferenceAlgebra.Cholesky(sigma)!;

    // A prior draw's surface at the training cases: A | Σ, τ ~ MN(0, (τG_c)⁻¹, Σ) and the intercept from N(0, Σ/(τ·mass)) under level centring's
    // prior about the levels at 0 for the zero-mean design; the fitted values f = 1αᵀ + Z_cA.
    internal static double[,] GewekeSurface(BassModel model, BassDesign d, double[,] root, double tau, int n, int k, Xoshiro256PlusPlus rng)
    {
        int m = d.Count;
        var centred = new double[m + 1, k];
        var v = new double[k][];
        for (int l = 0; l < k; l++)
        {
            v[l] = new double[m];
            if (m == 0) continue;
            double[] z = [.. Enumerable.Range(0, m).Select(_ => Normal.SampleStandard(rng))];
            d.Factor.SampleInnovation(z, v[l]);
        }
        for (int j = 0; j < m; j++)
            for (int a = 0; a < k; a++)
            {
                double mixed = 0.0;
                for (int l = 0; l <= a; l++) mixed += root[a, l] * v[l][j];
                centred[j + 1, a] = mixed / Math.Sqrt(tau);
            }
        var u0 = new double[k];
        for (int l = 0; l < k; l++) u0[l] = Normal.SampleStandard(rng);
        for (int a = 0; a < k; a++)
        {
            double mixed = 0.0;
            for (int l = 0; l <= a; l++) mixed += root[a, l] * u0[l];
            centred[0, a] = mixed / Math.Sqrt(tau * n);
        }
        return Fitted(model, d, centred);
    }

    // u (levels × K) ~ N(0, s²Σ) with nᵀu = 0 per output: Π z mixed by Σ's factor.
    internal static double[,] CentredEffects(int[] codes, int levels, int k, double s, double[,] root, Xoshiro256PlusPlus rng)
    {
        var counts = new double[levels];
        foreach (int c in codes) counts[c]++;
        double nn = counts.Sum(c => c * c);
        var z = new double[levels, k];
        for (int h = 0; h < levels; h++)
            for (int l = 0; l < k; l++) z[h, l] = Normal.SampleStandard(rng);
        var u = new double[levels, k];
        for (int l = 0; l < k; l++)
        {
            double dot = 0.0;
            for (int h = 0; h < levels; h++) dot += counts[h] * z[h, l];
            for (int h = 0; h < levels; h++) z[h, l] -= counts[h] * dot / nn;
        }
        for (int h = 0; h < levels; h++)
            for (int a = 0; a < k; a++)
            {
                double mixed = 0.0;
                for (int l = 0; l <= a; l++) mixed += root[a, l] * z[h, l];
                u[h, a] = s * mixed;
            }
        return u;
    }

    // f = 1αᵀ + Z_cA at the design's structure, from centred coefficients (α first, then the terms in the design's order).
    internal static double[,] Fitted(BassModel model, BassDesign design, double[,] centred)
    {
        int n = model.Data.Response.Cases, k = centred.GetLength(1);
        double[,] xc = MultiOutputCoreTests.Columns(model, design);
        var f = new double[n, k];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < k; l++)
            {
                double v = centred[0, l];
                for (int j = 0; j < design.Count; j++) v += xc[i, j] * centred[j + 1, l];
                f[i, l] = v;
            }
        return f;
    }

    // The centred coefficients in the uncentred basis: α − Σⱼ βⱼ·z̄ⱼ, the raw columns' means.
    private static double[,] Uncentred(BassModel model, BassDesign design, double[,] centred)
    {
        int n = model.Data.Response.Cases, k = centred.GetLength(1), m = design.Count;
        var a = (double[,])centred.Clone();
        var column = new double[n];
        for (int j = 0; j < m; j++)
        {
            model.Basis.CaseColumn(design.Terms[j], model.Inputs.Cases, column);
            double mean = column.Average();
            for (int l = 0; l < k; l++) a[0, l] -= centred[j + 1, l] * mean;
        }
        return a;
    }

    private static double MaxGap(double[,] a, double[,] b)
    {
        double worst = 0.0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) worst = Math.Max(worst, Math.Abs(a[i, j] - b[i, j]) / (1.0 + Math.Abs(a[i, j])));
        return worst;
    }

    private static double Relative(double a, double b) => Math.Abs(a - b) / Math.Max(1e-300, Math.Max(Math.Abs(a), Math.Abs(b)));
}

internal static class ArrayFill
{
    // An n × k array filled by f(i, l).
    internal static double[,] Select2(this double[,] shape, Func<int, int, double> f)
    {
        var a = new double[shape.GetLength(0), shape.GetLength(1)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int l = 0; l < a.GetLength(1); l++) a[i, l] = f(i, l);
        return a;
    }
}
