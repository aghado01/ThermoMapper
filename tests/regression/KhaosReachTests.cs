using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Samplers;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;
using static Maths.Regression.Tests.BassTemperingTests;

namespace Maths.Regression.Tests;

/// <summary>
/// KHAOS on the engine's reach (slice 5): tempering, SMC, the flat histogram, the pointwise density and the conformal layer's candidate
/// density on Legendre terms. The fixture is slice 0's exact-law fixture (<see cref="KhaosChainTests"/>: p 4, degree 3, order 2, M_max 2,
/// n 40, 466 structures) with σ² ~ IG(2, 0.5), so b = 0 and the evidence are defined, under the raw ridge at khaos's τ² = 10 held and the
/// modified g-prior (ζ = 1) with g₀² drawn.
/// <list type="bullet">
/// <item><b>The law at each rung,</b> from a ladder's own replicas, with and without the reference rung, against the enumeration's tempered
/// law (<see cref="BassTemperingTests.Exact"/>): the count law at every rung, and E[ln τ] where g₀² is drawn.</item>
/// <item><b>SMC and the ladder:</b> each fit's ln p(y) against the enumeration's, and its count law against the posterior's; SMC under the
/// fixed and the adaptive count of sweeps, and a ladder that steps its rungs at once, draw for draw as in order.</item>
/// <item><b>The flat histogram</b> against the enumeration (<see cref="BassFlatHistogramTests.MatchesTheEnumeration"/>).</item>
/// <item><b>The pointwise density:</b> the retained content against the Normal at the explicit uncentred design, and the collapsed content
/// against explicit refits under KHAOS's prior written out: τI on the uncentred coefficients, or τ·D_g⁻¹ΨᵀΨD_g⁻¹ with Ψ = [1, B].</item>
/// <item><b>The candidate density</b> at the training cases is the pointwise column, under the Gaussian and Student t noise.</item>
/// </list>
/// </summary>
public sealed class KhaosReachTests(ITestOutputHelper output)
{
    private const double TauShape = 2.0, TauRate = 80.0;   // g₀² about n: Gamma(2, rate 80) on τ = 1/g₀²

    internal static BassConfig Config(string prior, int maxBasis = 2)
        => KhaosChainTests.Config(raw: true, maxBasis: maxBasis, sigma2: Sigma2Prior.InverseGamma(2.0, 0.5)) with
        {
            Coefficients = prior == "modified g" ? BassCoefficientPrior.ModifiedG(TauPrior.Gamma(TauShape, TauRate)) : BassCoefficientPrior.RawRidge(TauPrior.Held(0.1)),
        };

    // ---- the law at each rung ---------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData("raw ridge", false, new[] { 1.0, 0.4, 0.0 })]
    [InlineData("raw ridge", true, new[] { 1.0, 0.4, 0.1, 0.0 })]
    [InlineData("modified g", false, new[] { 1.0, 0.4, 0.0 })]
    [InlineData("modified g", true, new[] { 1.0, 0.3, 0.05, 0.0 })]
    public void Ladder_EachRungSamplesItsTemperedLaw(string prior, bool reference, double[] rungs)
    {
        var spec = new TemperingSpec { Levels = rungs.Length, HotEnd = rungs[^1] > 0.0 ? rungs[^1] : rungs[^2], ReferenceRung = reference };
        BassConfig config = Config(prior);
        BassModel model = BassModel.Create(config with { Sampling = config.Sampling with { Tempering = spec } }, KhaosChainTests.Data());
        const int ladders = 4, batches = 20;
        (double[][,] count, double[][] logTau) = RungBatches(model, spec, rungs, ladders, batches);
        RungLaw[] exact = Exact(model, rungs);
        for (int r = 0; r < rungs.Length; r++)
        {
            for (int m = 0; m <= model.MaxBasis; m++)
            {
                (double mean, double error) = Summary([.. Enumerable.Range(0, ladders * batches).Select(k => count[r][m, k])]);
                output.WriteLine($"{prior}{(reference ? ", reference" : "")}, b = {rungs[r]:G3}, M = {m}: {mean:F4} against {exact[r].Count[m]:F4} (± {error:F4})");
                Assert.True(Math.Abs(mean - exact[r].Count[m]) <= 5 * error + 0.004, $"b = {rungs[r]}, M = {m}: {mean:F4} against {exact[r].Count[m]:F4} (± {error:F4})");
            }
            if (!model.Target.SamplesTau) continue;
            if (rungs[r] == 0.0) Assert.Equal(SpecialFunctions.Digamma(TauShape) - Math.Log(TauRate), exact[r].LogTauMean, 1e-9);   // τ's prior at b = 0
            (double u, double uError) = Summary(logTau[r]);
            output.WriteLine($"{prior}{(reference ? ", reference" : "")}, b = {rungs[r]:G3}, E[ln τ]: {u:F4} against {exact[r].LogTauMean:F4} (± {uError:F4})");
            Assert.True(Math.Abs(u - exact[r].LogTauMean) <= 5 * uError + 0.01, $"b = {rungs[r]}, E[ln τ] {u:F4} against {exact[r].LogTauMean:F4} (± {uError:F4})");
        }
    }

    // ---- SMC and the ladder: the posterior and the evidence --------------------------------------------------------------------------

    [Theory]
    [InlineData("raw ridge", "smc")]
    [InlineData("modified g", "smc")]
    [InlineData("raw ridge", "ladder")]
    [InlineData("modified g", "ladder")]
    [InlineData("modified g", "smc, adaptive")]
    public void Fit_MeetsTheEnumeratedEvidence(string prior, string sampler)
    {
        BassData data = KhaosChainTests.Data();
        BassConfig config = Config(prior);
        BassSampling sampling = sampler.StartsWith("smc")
            ? config.Sampling with
            {
                MasterSeed = 3, ReportEvidence = true,
                Smc = new SmcSpec { Particles = 2000, Replicates = 4, MutationSteps = 4, AdaptiveMutation = sampler.EndsWith("adaptive") ? new AdaptiveMutation { StayProbability = 0.1 } : null },
            }
            : config.Sampling with { Chains = 4, MasterSeed = 5, BurnIn = 1_000, Samples = 4_000, ReportEvidence = true, Tempering = new TemperingSpec { ReferenceRung = true } };
        config = config with { Sampling = sampling };
        BassModel model = BassModel.Create(config, data);
        RungLaw exact = Exact(model, 1.0)[0];
        BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        BassEvidence evidence = fit.Evidence!;
        output.WriteLine($"{prior}, {sampler}: ln p(y) {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against the enumeration's {exact.LogZ:F4}");
        Assert.True(Math.Abs(evidence.LogEvidence - exact.LogZ) <= 4 * evidence.MonteCarloError + 0.01,
                    $"{prior}, {sampler}: {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against {exact.LogZ:F4}");

        // The count law: each chain's (each replicate's) weighted frequencies, their spread for the error.
        double[] weights = fit.DrawWeights ?? [.. Enumerable.Repeat(1.0, fit.Draws.Sum(c => c.Count))];
        var perChain = new double[fit.Draws.Count][];
        int i = 0;
        for (int c = 0; c < fit.Draws.Count; c++)
        {
            perChain[c] = new double[model.MaxBasis + 1];
            double total = 0.0;
            foreach (BassDraw d in fit.Draws[c])
            {
                perChain[c][d.Structure.TermIds.Length] += weights[i];
                total += weights[i++];
            }
            for (int m = 0; m <= model.MaxBasis; m++) perChain[c][m] /= total;
        }
        for (int m = 0; m <= model.MaxBasis; m++)
        {
            (double mean, double error) = Summary([.. perChain.Select(p => p[m])]);
            output.WriteLine($"{prior}, {sampler}, P(M = {m}): {mean:F4} against {exact.Count[m]:F4} (± {error:F4})");
            Assert.True(Math.Abs(mean - exact.Count[m]) <= 5 * error + 0.01, $"P(M = {m}): {mean:F4} against {exact.Count[m]:F4} (± {error:F4})");
        }
    }

    // A tempered KHAOS fit steps its rungs at once as it steps them in order, draw for draw (BassRunOptions.RungWorkers).
    [Fact]
    public void ATemperedFit_IsTheSameWithItsRungsAtOnce()
    {
        BassData data = KhaosChainTests.Data();
        BassConfig config = Config("modified g") with
        {
            Sampling = new BassSampling { Chains = 2, MasterSeed = 9, BurnIn = 300, Samples = 200, Tempering = new TemperingSpec { ReferenceRung = true } },
        };
        string inOrder = BassBitIdentityTests.Hash(Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 8 }));
        string atOnce = BassBitIdentityTests.Hash(Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 8, RungWorkers = -1 }));
        output.WriteLine($"KHAOS tempered: {inOrder} in order, {atOnce} at once");
        Assert.Equal(inOrder, atOnce);
    }

    // ---- the flat histogram ------------------------------------------------------------------------------------------------------------

    // At a held τ the density of states is F's 466 point masses, which leave SAMC's limit undefined on 12 nodes; 6 nodes hold it.
    [Theory]
    [InlineData("raw ridge", 6)]
    [InlineData("modified g", 12)]
    public void FlatHistogram_MatchesTheEnumeration(string prior, int nodes)
        => BassFlatHistogramTests.MatchesTheEnumeration(BassModel.Create(Config(prior), KhaosChainTests.Data()), 0.0, $"KHAOS, {prior}", output, nodes);

    // ---- the pointwise density ---------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData("raw ridge")]
    [InlineData("modified g")]
    public void Pointwise_BothContents_MatchDenseReferencesAndRefits(string prior)
    {
        BassData data = KhaosChainTests.Data();
        (BassModel model, BassTermTable table, List<List<BassDraw>> draws) = BassPointwiseTests.Run(Config(prior, maxBasis: 4), data, seed: 11, kept: 12);
        PointwiseLogDensity retained = BassPointwise.LogDensity(model, table, draws, null);
        PointwiseLogDensity collapsed = BassPointwise.LogDensity(model, table, draws, null, PointwiseContent.CollapsedPredictive);
        BassDraw[] pooled = [.. draws.SelectMany(c => c)];
        Assert.Contains(pooled, d => d.Structure.TermIds.Length > 0);
        if (prior == "modified g")   // a structure met at two τ, so the collapsed content takes each τ's own forms
            Assert.Contains(pooled.GroupBy(d => string.Join(",", d.Structure.TermIds)), g => g.Select(d => d.Tau).Distinct().Count() > 1);

        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        var form = model.Config.Coefficients.Form as ModifiedGForm;
        double worstRetained = 0.0, worstCollapsed = 0.0;
        for (int s = 0; s < pooled.Length; s++)
        {
            BassDraw d = pooled[s];
            double[] log = BassPointwiseTests.DenseRetained(model, table, d);
            // Ψ = [1, B] as evaluated, uncentred; KHAOS's prior precision on its coefficients, the zero centre.
            int p = d.Structure.TermIds.Length + 1;
            var psi = new double[n, p];
            for (int k = 0; k < n; k++) psi[k, 0] = 1.0;
            var g = new double[p];
            g[0] = 1.0;
            for (int j = 1; j < p; j++)
            {
                BassTerm term = table[d.Structure.TermIds[j - 1]];
                double[] column = BassPointwiseTests.Column(model, term);
                for (int k = 0; k < n; k++) psi[k, j] = column[k];
                g[j] = form?.Factor(term) ?? 1.0;
            }
            var lambda0 = new double[p, p];
            for (int a = 0; a < p; a++)
                for (int b = 0; b < p; b++)
                {
                    double gram = 0.0;
                    for (int k = 0; k < n; k++) gram += psi[k, a] * psi[k, b];
                    lambda0[a, b] = form is null ? (a == b ? d.Tau : 0.0) : d.Tau * gram / (g[a] * g[b]);
                }
            for (int k = 0; k < n; k++)
            {
                worstRetained = Math.Max(worstRetained, Math.Abs(log[k] - retained.Values[s, k]) / Math.Max(1.0, Math.Abs(log[k])));
                double refit = BassPointwiseTests.Refit(psi, y, [k], lambda0, new double[p], p, model.Target.Sigma2);
                worstCollapsed = Math.Max(worstCollapsed, Math.Abs(refit - collapsed.Values[s, k]) / Math.Max(1.0, Math.Abs(refit)));
            }
        }
        output.WriteLine($"{prior}: {pooled.Length} draws over {pooled.Select(d => string.Join(",", d.Structure.TermIds)).Distinct().Count()} structures, " +
                         $"{n} units; retained against the Normal {worstRetained:E2}, collapsed against the refits {worstCollapsed:E2}");
        Assert.True(worstRetained <= 1e-10, $"retained {worstRetained:E2}");
        Assert.True(worstCollapsed <= 1e-8, $"collapsed {worstCollapsed:E2}");
    }

    // ---- the candidate density ---------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData("gaussian")]
    [InlineData("t")]
    public void Candidates_AtTheTrainingCases_AreThePointwiseColumn(string law)
    {
        BassData data = KhaosChainTests.Data();
        BassConfig config = Config("modified g", maxBasis: 4) with
        {
            // The zero centre fixes the level, so the noise reads no reference (BassNoiseReference.None).
            Noise = law == "t" ? BassNoise.StudentT(3, 0.7, BassNoiseReference.None) : null,
            Sampling = new BassSampling { Chains = 2, MasterSeed = 3, BurnIn = 300, Samples = 60 },
        };
        BassFit fit = Bass.Fit(config, data);
        PointwiseLogDensity pointwise = fit.PointwiseLogDensity(data);
        double[] y = data.Response.Values.ToArray();
        CandidateLogDensity candidates = fit.CandidateLogDensity(data, PredictionCases.Of(data), [.. y.Select(v => new[] { v })]);
        Assert.Equal(pointwise.Draws, candidates.Draws);
        Assert.Equal(pointwise.ChainLengths, candidates.ChainLengths);
        double worst = 0.0;
        for (int s = 0; s < pointwise.Draws; s++)
            for (int i = 0; i < y.Length; i++)
            {
                double a = pointwise.Values[s, i], b = candidates.Values[s, i];
                worst = Math.Max(worst, Math.Abs(a - b) / Math.Max(1.0, Math.Abs(a)));
            }
        output.WriteLine($"{law}: {pointwise.Draws} draws × {y.Length} cases, largest relative gap {worst:E2}");
        Assert.True(worst <= 1e-12, $"{law}: {worst:E2}");
        Assert.Equal(y.Length, ConformalBayes.FullStep(pointwise, candidates).Sets.Count);
    }
}
