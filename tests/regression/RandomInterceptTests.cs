using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Rng;
using Maths.Samplers;
using Xunit;
using Xunit.Abstractions;
using static Maths.Regression.Tests.BassTemperingTests;

namespace Maths.Regression.Tests;

/// <summary>
/// Random-effect intercepts as a ridge block (E21; <see cref="RandomInterceptBlock"/>), a general capability: any grouping, any coefficient
/// form, hinge and Legendre fits alike.
/// <list type="bullet">
/// <item><b>The block's algebra against dense matrices:</b> V^(−½) whitens V = I + s²ZΠZᵀ, keeps the constant column, and its ln|V|; the
/// group effects' conditional mean and covariance, and Σ_g n_g u_g = 0 in every draw.</item>
/// <item><b>The target against a dense multivariate t,</b> absolutely, the block's −½ln|V| included: y | S, τ ~ t with scale V + ΨS₀Ψᵀ under
/// σ² ~ IG(a, b), for the raw ridge and the modified g-prior, and without the block for the reference.</item>
/// <item><b>The engine's reach:</b> each rung's law from a ladder's replicas against the enumeration; SMC's and the ladder's ln p(y) against
/// it, SMC's flags clear; the flat histogram against it; the retained pointwise content against the Normal at f + u_g; on a hinge fit
/// (the g-form, τ drawn) and a Legendre fit (the modified g-prior, g₀² drawn). The collapsed content and the candidate density are
/// refused, with their reasons.</item>
/// <item><b>The convention and prediction:</b> the intercept carries ȳ; a case of an observed level reads f + u_g, one of the population
/// f, with its own effect in the predictive draw; an unseen level is refused or placed at the population.</item>
/// <item><b>The config and data:</b> refusals, and the serialized form, which a config without the block keeps.</item>
/// </list>
/// </summary>
public sealed class RandomInterceptTests(ITestOutputHelper output)
{
    // An unbalanced grouping of n cases on `levels` levels, each level taking at least two cases.
    internal static CategoricalPredictor Grouping(int n, int levels, int seed, string name = "g")
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var codes = new int[n];
        for (int i = 0; i < n; i++) codes[i] = i < 2 * levels ? i / 2 : (int)Math.Min(levels - 1, Math.Floor(levels * Math.Pow(rng.NextDouble(), 1.6)));
        var order = Enumerable.Range(0, n).OrderBy(_ => rng.NextDouble()).ToArray();
        return new CategoricalPredictor(name, [.. order.Select(i => codes[i])], [.. Enumerable.Range(0, levels).Select(l => $"L{l}")]);
    }

    internal static BassData Grouped(BassData data, int levels, int seed = 3) => new(data.Predictors, data.Response, [Grouping(data.Response.Cases, levels, seed)]);

    // Z (cases × observed levels), the counts and V = I + s²ZΠZᵀ, Π the projector onto nᵀu = 0, written out.
    private static (double[,] Z, double[] Counts, double[,] Pi, double[,] V) Dense(CategoricalPredictor grouping, double s)
    {
        int[] codes = grouping.Codes.ToArray();
        int n = codes.Length, l = codes.Max() + 1;
        var z = new double[n, l];
        var counts = new double[l];
        for (int i = 0; i < n; i++) { z[i, codes[i]] = 1.0; counts[codes[i]]++; }
        double nn = counts.Sum(c => c * c);
        var pi = new double[l, l];
        for (int a = 0; a < l; a++)
            for (int b = 0; b < l; b++) pi[a, b] = (a == b ? 1.0 : 0.0) - counts[a] * counts[b] / nn;
        var v = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) v[i, j] = (i == j ? 1.0 : 0.0) + s * s * pi[codes[i], codes[j]];
        return (z, counts, pi, v);
    }

    private static double[,] Multiply(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), k = a.GetLength(1), m = b.GetLength(1);
        var c = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
            {
                double sum = 0.0;
                for (int q = 0; q < k; q++) sum += a[i, q] * b[q, j];
                c[i, j] = sum;
            }
        return c;
    }

    private static double[,] Inverse(double[,] a)
    {
        int n = a.GetLength(0);
        var f = new DenseCholesky(n);
        f.Decompose(a);
        Assert.True(f.Status.PositiveDefinite);
        var inv = new double[n, n];
        for (int j = 0; j < n; j++)
        {
            var e = new double[n];
            e[j] = 1.0;
            double[] col = f.Solve(e);
            for (int i = 0; i < n; i++) inv[i, j] = col[i];
        }
        return inv;
    }

    private static double LogDet(double[,] a)
    {
        var f = new DenseCholesky(a.GetLength(0));
        f.Decompose(a);
        Assert.True(f.Status.PositiveDefinite);
        return f.LogDet;
    }

    // ---- the block's algebra -------------------------------------------------------------------------------------------------------

    [Fact]
    public void Whitening_IsVToTheMinusHalf_KeepsTheConstant_AndReadsLnV()
    {
        const double s = 1.7;
        CategoricalPredictor grouping = Grouping(30, 6, 11);
        RandomInterceptBlock block = RandomInterceptBlock.Of(grouping, new BassRandomIntercept("g", s));
        (_, _, _, double[,] v) = Dense(grouping, s);
        int n = 30;
        var w = new double[n, n];
        for (int j = 0; j < n; j++)
        {
            var e = new double[n];
            e[j] = 1.0;
            block.Whiten(e);
            for (int i = 0; i < n; i++) w[i, j] = e[i];
        }
        double[,] wvw = Multiply(Multiply(w, v), w);
        double worst = 0.0, asymmetry = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                worst = Math.Max(worst, Math.Abs(wvw[i, j] - (i == j ? 1.0 : 0.0)));
                asymmetry = Math.Max(asymmetry, Math.Abs(w[i, j] - w[j, i]));
            }
        double[] ones = [.. Enumerable.Repeat(1.0, n)];
        block.Whiten(ones);
        output.WriteLine($"W·V·W against I {worst:E2}, W against Wᵀ {asymmetry:E2}, W·1 against 1 {ones.Max(x => Math.Abs(x - 1.0)):E2}; " +
                         $"ln|V| {block.LogDeterminant:F12} against the dense {LogDet(v):F12}");
        Assert.True(worst < 1e-12 && asymmetry < 1e-13, $"{worst:E2}, {asymmetry:E2}");
        Assert.All(ones, x => Assert.Equal(1.0, x, 12));
        Assert.Equal(LogDet(v), block.LogDeterminant, 10);
    }

    [Fact]
    public void GroupEffects_FollowTheirConditional_AndCentreOnTheObservations()
    {
        const double s = 0.8, sigma2 = 0.7;
        CategoricalPredictor grouping = Grouping(30, 6, 12);
        RandomInterceptBlock block = RandomInterceptBlock.Of(grouping, new BassRandomIntercept("g", s));
        (double[,] z, double[] counts, double[,] pi, double[,] v) = Dense(grouping, s);
        int n = 30, l = counts.Length;
        var rng = new Xoshiro256PlusPlus(5);
        double[] r = [.. Enumerable.Range(0, n).Select(_ => Normal.SampleStandard(rng))];
        // Dense: E[u | r] = s²ΠZᵀV⁻¹r, Cov = σ²(s²Π − s⁴ΠZᵀV⁻¹ZΠ).
        double[,] vInv = Inverse(v), zt = new double[l, n];
        for (int i = 0; i < n; i++)
            for (int g = 0; g < l; g++) zt[g, i] = z[i, g];
        double[,] gain = Multiply(Multiply(pi, zt), vInv);   // ΠZᵀV⁻¹
        double[,] cov = Multiply(Multiply(gain, z), pi);
        var mean = new double[l];
        for (int g = 0; g < l; g++)
        {
            for (int i = 0; i < n; i++) mean[g] += s * s * gain[g, i] * r[i];
            for (int h = 0; h < l; h++) cov[g, h] = sigma2 * (s * s * pi[g, h] - s * s * s * s * cov[g, h]);
        }
        double[] sums = block.GroupSums(r), at = block.EffectMean(sums);
        double meanGap = Enumerable.Range(0, l).Max(g => Math.Abs(at[g] - mean[g]));
        const int draws = 200_000;
        var m1 = new double[l];
        var m2 = new double[l, l];
        double worstCentre = 0.0;
        for (int k = 0; k < draws; k++)
        {
            double[] u = block.DrawEffects(sums, sigma2, rng);
            worstCentre = Math.Max(worstCentre, Math.Abs(Enumerable.Range(0, l).Sum(g => counts[g] * u[g])));
            for (int g = 0; g < l; g++)
            {
                m1[g] += u[g] / draws;
                for (int h = 0; h < l; h++) m2[g, h] += (u[g] - mean[g]) * (u[h] - mean[h]) / draws;
            }
        }
        double scale = Enumerable.Range(0, l).Max(g => cov[g, g]), covGap = 0.0, drawGap = 0.0;
        for (int g = 0; g < l; g++)
        {
            drawGap = Math.Max(drawGap, Math.Abs(m1[g] - mean[g]) / Math.Sqrt(cov[g, g] / draws));
            for (int h = 0; h < l; h++) covGap = Math.Max(covGap, Math.Abs(m2[g, h] - cov[g, h]) / scale);
        }
        output.WriteLine($"the conditional mean against the dense {meanGap:E2}; over {draws:N0} draws the mean's worst z {drawGap:F2} and the covariance's " +
                         $"worst gap {covGap:F4} of the largest variance; Σ n_g u_g at most {worstCentre:E2}");
        Assert.True(meanGap < 1e-12, $"{meanGap:E2}");
        Assert.True(drawGap < 4.5, $"{drawGap:F2}");
        Assert.True(covGap < 0.02, $"{covGap:F4}");
        Assert.True(worstCentre < 1e-11, $"{worstCentre:E2}");
    }

    // ---- the target against a dense multivariate t ---------------------------------------------------------------------------------

    // ln p(y) under y | σ² ~ N(0, σ²C), σ² ~ IG(a, b): Γ(a + n/2)/Γ(a)·bᵃ(2π)^(−n/2)|C|^(−½)(b + yᵀC⁻¹y/2)^(−(a + n/2)).
    private static double DenseT(double[,] c, double[] y, double a, double b)
    {
        int n = y.Length;
        double[,] inv = Inverse(c);
        double q = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) q += y[i] * inv[i, j] * y[j];
        return SpecialFunctions.LogGamma(a + 0.5 * n) - SpecialFunctions.LogGamma(a) + a * Math.Log(b) - 0.5 * n * Math.Log(2 * Math.PI)
               - 0.5 * LogDet(c) - (a + 0.5 * n) * Math.Log(b + 0.5 * q);
    }

    [Theory]
    [InlineData("raw ridge", true)]
    [InlineData("raw ridge", false)]
    [InlineData("modified g", true)]
    [InlineData("modified g", false)]
    public void Target_IsTheDenseT_Absolutely(string prior, bool block)
    {
        const double s = 1.3, a = 2.0, b = 0.5;
        BassData data = Grouped(KhaosChainTests.Data(), 6);
        if (!block) data = new BassData(data.Predictors, data.Response);
        BassConfig config = KhaosReachTests.Config(prior) with { RandomIntercepts = block ? [new BassRandomIntercept("g", s)] : null };
        BassModel model = BassModel.Create(config, data);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        double[,] v = block ? Dense(data.Groupings[0], s).V : Identity(n);
        var form = config.Coefficients.Form as ModifiedGForm;
        List<BassTerm> terms = [.. model.Prior.Legendre!.Terms()];
        BassTerm[][] structures = [[], [terms[0]], [terms[3], terms[9]], [terms[5], terms[20]]];
        double worst = 0.0;
        foreach (BassTerm[] st in structures)
        {
            int p = st.Length + 1;
            var psi = new double[n, p];
            for (int i = 0; i < n; i++) psi[i, 0] = 1.0;
            var g = new double[p];
            g[0] = 1.0;
            for (int j = 1; j < p; j++)
            {
                var col = new double[n];
                model.Basis.CaseColumn(st[j - 1], model.Inputs.Cases, col);
                for (int i = 0; i < n; i++) psi[i, j] = col[i];
                g[j] = form?.Factor(st[j - 1]) ?? 1.0;
            }
            // The modified g-prior reads the whitened design, ΨᵀV⁻¹Ψ; the raw ridge, the identity.
            double[,] info = Multiply(Multiply(Transpose(psi), Inverse(v)), psi);
            foreach (double logTau in prior == "raw ridge" ? new[] { Math.Log(0.1) } : new[] { -4.0, -1.5, 1.0 })
            {
                BassState state = BassMixing.Loaded(model, st, logTau, out _) ?? throw new InvalidOperationException("outside F");
                double tau = Math.Exp(logTau);
                double[,] s0 = form is null ? Scaled(Identity(p), 1.0 / tau) : Scaled(DgInverseDg(Inverse(info), g), 1.0 / tau);
                double[,] c = Add(v, Multiply(Multiply(psi, s0), Transpose(psi)));
                double dense = DenseT(c, y, a, b), engine = model.Target.Evaluate(state).LogLikelihood(1.0);
                worst = Math.Max(worst, Math.Abs(engine - dense));
            }
        }
        output.WriteLine($"{prior}, {(block ? $"a block of s = {s}" : "no block")}: the target's ln p(y | S, τ) against the dense t, worst absolute gap {worst:E2}" +
                         (block ? $"; the block's −½ln|V| = {-0.5 * model.RandomEffects!.LogDeterminant:F4}" : ""));
        Assert.True(worst < 1e-8, $"{worst:E2}");
    }

    private static double[,] Identity(int n)
    {
        var m = new double[n, n];
        for (int i = 0; i < n; i++) m[i, i] = 1.0;
        return m;
    }

    private static double[,] Transpose(double[,] a)
    {
        var t = new double[a.GetLength(1), a.GetLength(0)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) t[j, i] = a[i, j];
        return t;
    }

    private static double[,] Scaled(double[,] a, double c)
    {
        var m = (double[,])a.Clone();
        for (int i = 0; i < m.GetLength(0); i++)
            for (int j = 0; j < m.GetLength(1); j++) m[i, j] *= c;
        return m;
    }

    private static double[,] Add(double[,] a, double[,] b)
    {
        var m = (double[,])a.Clone();
        for (int i = 0; i < m.GetLength(0); i++)
            for (int j = 0; j < m.GetLength(1); j++) m[i, j] += b[i, j];
        return m;
    }

    private static double[,] DgInverseDg(double[,] inv, double[] g)
    {
        var m = (double[,])inv.Clone();
        for (int i = 0; i < m.GetLength(0); i++)
            for (int j = 0; j < m.GetLength(1); j++) m[i, j] *= g[i] * g[j];
        return m;
    }

    // ---- the engine's reach ------------------------------------------------------------------------------------------------------------

    // A hinge fit (BassSweepTests' fixture, the g-form with τ drawn, the defaults' proper priors) and a Legendre fit (slice 0's exact-law
    // fixture, the modified g-prior with g₀² drawn), each with a block of s = 1.3 on six unbalanced levels.
    internal static (BassModel Model, BassConfig Config) Fixture(string basis, TemperingSpec? spec = null)
    {
        if (basis == "hinge")
        {
            BassConfig hinge = BassTemperingTests.Config(spec) with { RandomIntercepts = [new BassRandomIntercept("g", 1.3)] };
            return (BassModel.Create(hinge, Grouped(BassSweepTests.Data(), 4)), hinge);
        }
        BassConfig legendre = KhaosReachTests.Config("modified g") with { RandomIntercepts = [new BassRandomIntercept("g", 1.3)] };
        if (spec is not null) legendre = legendre with { Sampling = legendre.Sampling with { Tempering = spec } };
        return (BassModel.Create(legendre, Grouped(KhaosChainTests.Data(), 6)), legendre);
    }

    [Theory]
    [InlineData("hinge", true)]
    [InlineData("legendre", true)]
    [InlineData("legendre", false)]
    public void Ladder_EachRungSamplesItsTemperedLaw(string basis, bool reference)
    {
        double[] rungs = reference ? [1.0, 0.4, 0.1, 0.0] : [1.0, 0.4, 0.0];
        var spec = new TemperingSpec { Levels = rungs.Length, HotEnd = rungs[^1] > 0.0 ? rungs[^1] : rungs[^2], ReferenceRung = reference };
        (BassModel model, _) = Fixture(basis, spec);
        const int ladders = 4, batches = 20;
        (double[][,] count, double[][] logTau) = RungBatches(model, spec, rungs, ladders, batches);
        RungLaw[] exact = Exact(model, rungs);
        for (int r = 0; r < rungs.Length; r++)
        {
            for (int m = 0; m <= model.MaxBasis; m++)
            {
                (double mean, double error) = Summary([.. Enumerable.Range(0, ladders * batches).Select(k => count[r][m, k])]);
                output.WriteLine($"{basis}, b = {rungs[r]:G3}, M = {m}: {mean:F4} against {exact[r].Count[m]:F4} (± {error:F4})");
                Assert.True(Math.Abs(mean - exact[r].Count[m]) <= 5 * error + 0.004, $"b = {rungs[r]}, M = {m}: {mean:F4} against {exact[r].Count[m]:F4}");
            }
            (double u, double uError) = Summary(logTau[r]);
            Assert.True(Math.Abs(u - exact[r].LogTauMean) <= 5 * uError + 0.01, $"b = {rungs[r]}, E[ln τ] {u:F4} against {exact[r].LogTauMean:F4}");
        }
    }

    [Theory]
    [InlineData("hinge", "smc")]
    [InlineData("legendre", "smc")]
    [InlineData("hinge", "ladder")]
    [InlineData("legendre", "ladder")]
    public void Fit_MeetsTheEnumeratedEvidence(string basis, string sampler)
    {
        (BassModel model, BassConfig config) = Fixture(basis);
        BassSampling sampling = sampler == "smc"
            ? config.Sampling with { Chains = 1, MasterSeed = 3, ReportEvidence = true, Smc = new SmcSpec { Particles = 2000, Replicates = 4, MutationSteps = 4 }, Tempering = null }
            : config.Sampling with { Chains = 4, MasterSeed = 5, BurnIn = 1_000, Samples = 4_000, ReportEvidence = true, Tempering = new TemperingSpec { ReferenceRung = true } };
        config = config with { Sampling = sampling };
        BassFit fit = Bass.Fit(config, model.Data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        double exact = Exact(model, 1.0)[0].LogZ;
        BassEvidence evidence = fit.Evidence!;
        output.WriteLine($"{basis}, {sampler}: ln p(y) {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against the enumeration's {exact:F4}");
        Assert.True(Math.Abs(evidence.LogEvidence - exact) <= 4 * evidence.MonteCarloError + 0.01, $"{evidence.LogEvidence:F4} against {exact:F4}");
        if (fit.Diagnostics.Smc is { } readout) Assert.Empty(readout.Flags);
        // Every retained draw carries its group effects, centred on the observations.
        RandomEffectGrouping block = model.RandomEffects!.Layout(0);
        foreach (BassDraw d in fit.Draws.SelectMany(c => c))
        {
            Assert.Equal(block.Levels, d.Groups![0].Effects.Length);
            double centre = 0.0;
            for (int g = 0; g < block.Levels; g++) centre += block.Counts[g] * d.Groups![0].Effects[g];
            Assert.True(Math.Abs(centre) < 1e-9, $"Σ n_g u_g = {centre:E2}");
        }
    }

    [Fact]
    public void FlatHistogram_MatchesTheEnumeration()
    {
        (BassModel model, _) = Fixture("legendre");
        BassFlatHistogramTests.MatchesTheEnumeration(model, 0.0, "a Legendre fit with a block", output);
    }

    // ---- the pointwise density, prediction and the convention ----------------------------------------------------------------------

    [Fact]
    public void Pointwise_ReadsEachDrawsGroupEffects_AndRefusesWhatIsNotBuilt()
    {
        (BassModel model, BassConfig config) = Fixture("hinge");
        BassFit fit = Bass.Fit(config with { Sampling = config.Sampling with { Chains = 2, BurnIn = 300, Samples = 30, Thin = 3 } }, model.Data);
        PointwiseLogDensity retained = fit.PointwiseLogDensity(model.Data);
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        double[] y = model.Data.Response.Values.ToArray();
        RandomEffectGrouping block = model.RandomEffects!.Layout(0);
        double worst = 0.0;
        for (int s = 0; s < pooled.Length; s++)
        {
            for (int i = 0; i < y.Length; i++)
            {
                double mu = DenseFit(model, fit.Terms, pooled[s], i) + pooled[s].Groups![0].Effects[block.GroupOf(i)];
                double expected = Normal.MeanVariance(mu, pooled[s].Sigma2).LogPdf(y[i]);
                worst = Math.Max(worst, Math.Abs(expected - retained.Values[s, i]) / Math.Max(1.0, Math.Abs(expected)));
            }
        }
        output.WriteLine($"{pooled.Length} draws: the retained content against N(yᵢ | fᵢ + u_g, σ²), largest gap {worst:E2}");
        Assert.True(worst < 1e-10, $"{worst:E2}");
        // The group effects are reported as conditioned on, a group's; the whitening is the block's integral, not an offset.
        Assert.Contains(retained.Latents, r => r is { Scope: LatentScope.Group, Treatment: LatentTreatment.Conditioned } && r.Name.Contains("group effects"));
        Assert.DoesNotContain(retained.Latents, r => r.Name.Contains("offset"));
        Assert.Contains("random intercepts", Assert.Throws<NotSupportedException>(() => fit.PointwiseLogDensity(model.Data, content: PointwiseContent.CollapsedPredictive)).Message);
        Assert.Contains("random-intercept", Assert.Throws<NotSupportedException>(() =>
            fit.CandidateLogDensity(model.Data, PredictionCases.Of(model.Data), [.. y.Select(v => new[] { v })])).Message);
    }

    // f at case i from the draw's coefficients and the terms' columns, as evaluated.
    private static double DenseFit(BassModel model, BassTermTable terms, BassDraw draw, int i)
    {
        double f = draw.Coefficients[0];
        var col = new double[model.Data.Response.Cases];
        for (int q = 0; q < draw.Structure.TermIds.Length; q++)
        {
            model.Basis.CaseColumn(terms[draw.Structure.TermIds[q]], model.Inputs.Cases, col);
            f += draw.Coefficients[q + 1] * col[i];
        }
        return f;
    }

    [Fact]
    public void Prediction_ReadsTheCasesLevels_AndTheConventionKeepsTheLevel()
    {
        (BassModel model, BassConfig config) = Fixture("hinge");
        BassData data = model.Data;
        BassFit fit = Bass.Fit(config with { Sampling = config.Sampling with { Chains = 2, BurnIn = 300, Samples = 400, Thin = 3 } }, data);
        CategoricalPredictor grouping = data.Groupings[0];
        RandomEffectGrouping block = model.RandomEffects!.Layout(0);
        var options = new BassPredictOptions { Predictive = true, Seed = 9 };
        BassPrediction grouped = fit.Predict(data, PredictionCases.Of(data), options);
        BassPrediction population = fit.Predict(data, PredictionCases.Of(data.Predictors), options);
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        double worst = 0.0;
        for (int d = 0; d < pooled.Length; d++)
            for (int i = 0; i < data.Response.Cases; i++)
                worst = Math.Max(worst, Math.Abs(grouped.Latent[d, i, 0] - population.Latent[d, i, 0] - pooled[d].Groups![0].Effects[block.GroupOf(i)]));
        Assert.True(worst < 1e-12, $"{worst:E2}");
        // A population case's predictive draw carries its own effect: its spread about f is σ√(1 + s²) where a grouped case's is σ.
        double grouping2 = 0.0, population2 = 0.0, sigma2 = 0.0;
        int count = 0;
        for (int d = 0; d < pooled.Length; d++)
            for (int i = 0; i < data.Response.Cases; i++)
            {
                grouping2 += Math.Pow(grouped.Predictive![d, i, 0] - grouped.Latent[d, i, 0], 2);
                population2 += Math.Pow(population.Predictive![d, i, 0] - population.Latent[d, i, 0], 2);
                sigma2 += pooled[d].Sigma2;
                count++;
            }
        output.WriteLine($"predictive spread over σ²: grouped cases {grouping2 / sigma2:F3}, population cases {population2 / sigma2:F3} against 1 + s² = {1 + 1.3 * 1.3:F3}");
        Assert.InRange(grouping2 / sigma2, 0.94, 1.06);
        Assert.InRange(population2 / sigma2 / (1 + 1.3 * 1.3), 0.94, 1.06);
        // The training fit reads the prediction, f + u_g, and puts back no offset: the whitening is not one.
        BassTrainingFit training = fit.TrainingFit(data);
        double[,] mean = grouped.Mean();
        Assert.Null(training.Offset);
        for (int i = 0; i < data.Response.Cases; i++)
        {
            Assert.Equal(mean[i, 0], training.Mean[i, 0], 12);
            Assert.Equal(training.Observed[i, 0] - mean[i, 0], training.Residuals[i, 0], 12);
        }
        // An unseen level: refused by default, the population when asked.
        var unseen = new PredictionCases([.. PredictionCases.Of(data.Predictors).Columns, new CategoricalColumn(grouping.Name, [.. Enumerable.Repeat("new", data.Response.Cases)])]);
        Assert.Throws<ArgumentException>(() => fit.Predict(data, unseen));
        BassPrediction placed = fit.Predict(data, unseen, new BassPredictOptions { UnseenLevels = UnseenLevel.OutsideEverySubset, Predictive = false });
        Assert.Equal(population.Latent, placed.Latent);
        // The convention: the block has no level of its own, so each draw's mean fit over the cases is the intercept's and the terms'.
        double[] y = data.Response.Values.ToArray();
        double meanFit = 0.0;
        for (int d = 0; d < pooled.Length; d++)
        {
            double level = 0.0, terms = 0.0;
            for (int i = 0; i < y.Length; i++) { level += grouped.Latent[d, i, 0] / y.Length; terms += population.Latent[d, i, 0] / y.Length; }
            Assert.Equal(terms, level, 12);
            meanFit += level / pooled.Length;
        }
        output.WriteLine($"the fits' mean level {meanFit:F4} against ȳ {y.Average():F4}");
        Assert.True(Math.Abs(meanFit - y.Average()) < 0.05 * y.Select(v => Math.Abs(v - y.Average())).Average(), "the level strays from ȳ");
    }

    [Fact]
    public void MeanSurface_IsTheDensePosteriorMean()
    {
        // The raw ridge holds τ, so each draw's E[f + u_g | S, y] is the dense joint posterior's: C = V + ΨΨᵀ/τ, E[a] = ΨᵀC⁻¹y/τ and
        // E[u] = s²ΠZᵀC⁻¹y, neither reading σ². The engine's goes through the whitened design and the block's mean instead.
        const double s = 1.3, tau = 0.1;
        BassData data = Grouped(KhaosChainTests.Data(), 6);
        BassConfig config = KhaosReachTests.Config("raw ridge") with { RandomIntercepts = [new BassRandomIntercept("g", s)] };
        config = config with { Sampling = config.Sampling with { Chains = 2, BurnIn = 200, Samples = 15, Thin = 4 } };
        BassFit fit = Bass.Fit(config, data);
        BassModel model = BassModel.Create(config, data);
        BassPrediction prediction = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { MeanSurface = true });
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        (double[,] z, _, double[,] pi, double[,] v) = Dense(data.Groupings[0], s);
        int l = pi.GetLength(0);
        var dense = new double[n];
        for (int d = 0; d < prediction.DrawCount; d++)
        {
            BassDraw draw = pooled[prediction.Draws[d]];
            int p = draw.Structure.TermIds.Length + 1;
            var psi = new double[n, p];
            var col = new double[n];
            for (int i = 0; i < n; i++) psi[i, 0] = 1.0;
            for (int j = 1; j < p; j++)
            {
                model.Basis.CaseColumn(fit.Terms[draw.Structure.TermIds[j - 1]], model.Inputs.Cases, col);
                for (int i = 0; i < n; i++) psi[i, j] = col[i];
            }
            double[,] c = Add(v, Scaled(Multiply(psi, Transpose(psi)), 1.0 / tau));
            double[,] inverse = Inverse(c);
            var w = new double[n];   // C⁻¹y
            for (int i = 0; i < n; i++)
                for (int k = 0; k < n; k++) w[i] += inverse[i, k] * y[k];
            var a = new double[p];
            for (int j = 0; j < p; j++)
                for (int i = 0; i < n; i++) a[j] += psi[i, j] * w[i] / tau;
            var zw = new double[l];
            for (int i = 0; i < n; i++)
                for (int g = 0; g < l; g++) zw[g] += z[i, g] * w[i];
            var u = new double[l];
            for (int g = 0; g < l; g++)
                for (int h = 0; h < l; h++) u[g] += s * s * pi[g, h] * zw[h];
            for (int i = 0; i < n; i++)
            {
                double f = 0.0;
                for (int j = 0; j < p; j++) f += psi[i, j] * a[j];
                for (int g = 0; g < l; g++) f += z[i, g] * u[g];
                dense[i] += prediction.Weights[d] * f;
            }
        }
        double worst = Enumerable.Range(0, n).Max(i => Math.Abs(prediction.MeanSurface![i, 0] - dense[i]));
        output.WriteLine($"{prediction.DrawCount} draws: the mean surface at the training cases against the dense posterior mean, worst gap {worst:E2}");
        Assert.True(worst < 1e-9, $"{worst:E2}");
    }

    [Fact]
    public void Resume_KeepsTheGroupEffects()
    {
        (BassModel model, BassConfig config) = Fixture("hinge");
        config = config with { Sampling = config.Sampling with { Chains = 2, MasterSeed = 3, BurnIn = 200, Samples = 400 } };
        BassFit whole = Bass.Fit(config, model.Data);
        BassConfig half = config with { Sampling = config.Sampling with { Samples = 200 } };
        BassFit first = Bass.Fit(half, model.Data, new BassRunOptions { Checkpoint = true });
        BassRunTests.AssertSameFit(whole, Bass.Resume(first.Checkpoint!, half, model.Data));
        Assert.All(whole.Draws.SelectMany(c => c), d => Assert.Equal(model.RandomEffects!.Layout(0).Levels, d.Groups![0].Effects.Length));
    }

    // ---- the config and data ---------------------------------------------------------------------------------------------------------

    [Fact]
    public void Config_AndData_RefuseWhatIsNotBuilt_AndTheFormStands()
    {
        BassData data = BassSweepTests.Data();
        BassData grouped = Grouped(data, 4);
        BassConfig plain = BassSweepTests.Config(), block = plain with { RandomIntercepts = [new BassRandomIntercept("g", 1.0)] };
        Assert.DoesNotContain("RandomIntercept", plain.ToJson());
        Assert.Contains("\"RandomIntercepts\":[{\"Grouping\":\"g\",\"Scale\":1", block.ToJson());
        Assert.Equal(block.ToJson(), BassConfig.FromJson(block.ToJson()).ToJson());
        Assert.Equal(block, BassConfig.FromJson(block.ToJson()));
        Assert.NotEqual(Bass.Fingerprint(plain, data), Bass.Fingerprint(block, grouped));
        Assert.Contains("only random effects read", Assert.Throws<ArgumentException>(() => BassModel.Create(plain, grouped)).Message);
        Assert.Contains("has no data grouping", Assert.Throws<ArgumentException>(() => BassModel.Create(block, data)).Message);
        Assert.Contains("RandomIntercepts[0].Scale", Assert.Throws<ArgumentException>(() => (block with { RandomIntercepts = [new BassRandomIntercept("g", 0.0)] }).Validate()).Message);
        Assert.Contains("take no Noise", Assert.Throws<ArgumentException>(() => (block with { Noise = BassNoise.StudentT(5) }).Validate()).Message);
        Assert.Contains("used twice", Assert.Throws<ArgumentException>(() =>
            new BassData(data.Predictors, data.Response, [new CategoricalPredictor(data.Predictors[0].Name,
                                                                                [.. Enumerable.Range(0, data.Response.Cases).Select(i => i % 2)], ["a", "b"])])).Message);
        BassData functional = BassPointwiseTests.Data(BassPointwiseTests.Kind.Functional);
        var functionalGrouped = new BassData(functional.Predictors, functional.Response, [Grouping(functional.Response.Cases, 3, 1)]);
        Assert.Throws<ArgumentException>(() => BassModel.Create(block, functionalGrouped));
        var curves = new double[12, 5];
        for (int i = 0; i < 12; i++)
            for (int l = 0; l < 5; l++) curves[i, l] = Math.Sin(i + 0.7 * l) + 0.1 * i * l;
        BassDecomposition decomposition = BassDecomposition.FromPca(curves, new BassPcaSpec { VarianceShare = 0.9 });
        Assert.Contains("random intercepts", Assert.Throws<ArgumentException>(() => BassComponents.Fit(decomposition,
            [new ContinuousPredictor("x", [.. Enumerable.Range(0, 12).Select(i => i / 11.0)])], [.. Enumerable.Repeat(block, decomposition.Rank)], 1)).Message);
        // Without the block, a fit's draws carry no group effects and its config's form stands.
        Assert.Null(Bass.Fit(plain with { Sampling = plain.Sampling with { Chains = 1, BurnIn = 10, Samples = 3 } }, data).Draws[0][0].Groups);
    }
}
