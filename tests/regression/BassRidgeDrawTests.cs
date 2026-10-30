using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// What reads the normalized ridge's marginal (D03's ridge, decided 2026-10-29): the composed draws (τ, σ², the
/// coefficients) on a fixed structure against a long joint-Gibbs reference that samples the coefficients, σ² and τ in
/// turn; and fits under the ridge, functional and mixed, sampled and held, hyper-g (improper at a = 1.5 too), tempered
/// with the reference rung and by SMC, and under the improper options, with the g-form under hyper-g/n beside them, each
/// with every view over its draws. The structure posterior and the evidence against the enumeration oracle's exact ones
/// are the parity tests'.
/// </summary>
public sealed class BassRidgeDrawTests
{
    private readonly ITestOutputHelper _out;
    public BassRidgeDrawTests(ITestOutputHelper output) => _out = output;

    private static BassModel Model() => BassRidgeTests.Model();

    private static double[][] CentredColumns(BassModel model, IReadOnlyList<BassTerm> terms) => BassRidgeTests.CentredColumns(model, terms);

    [Fact]
    public void FixedStructure_ComposedDraws_MatchAJointGibbsReference()
    {
        BassModel model = Model();
        BassTerm[] structure = BassCoefficientPriorTests.TwoTerms;
        var moves = new BassStructureMoves(model, structure, 0.0);
        var chain = new ReversibleJumpChain<BassState>([new BassCoefficientPriorTests.Stay()], model.Target, moves.Start, new Xoshiro256PlusPlus(7));
        var scan = new BlockScan<ReversibleJumpChain<BassState>>([new BassTauBlock(model.Target)], chain, new Xoshiro256PlusPlus(8), c => c);
        var composer = new Xoshiro256PlusPlus(9);
        var table = new BassTermTable();
        const int sweeps = 100_000;
        var tau = new double[sweeps];
        var sigma = new double[sweeps];
        var a = new double[3][];
        for (int j = 0; j < 3; j++) a[j] = new double[sweeps];
        for (int i = 0; i < 1_000; i++) scan.Step();
        for (int i = 0; i < sweeps; i++)
        {
            BassDraw draw = moves.Compose(scan.Step().Current, composer, table);
            tau[i] = draw.Tau;
            sigma[i] = draw.Sigma2;
            for (int j = 0; j < 3; j++) a[j][i] = draw.Coefficients[j];
        }
        (double[] gTau, double[] gSigma, double[][] gA) = JointGibbs(model, structure, table, 400_000, seed: 31);

        void Check(string name, double[] collapsed, double[] gibbs)
        {
            double c = collapsed.Average(), g = gibbs.Average();
            double error = Math.Sqrt(Math.Pow(BassCoefficientPriorTests.BatchError(collapsed), 2) + Math.Pow(BassCoefficientPriorTests.BatchError(gibbs), 2));
            _out.WriteLine($"{name}: collapsed {c:G6}, joint Gibbs {g:G6} (± {error:G2})");
            Assert.True(Math.Abs(c - g) <= 5 * error, $"{name}: collapsed {c:G6} against joint Gibbs {g:G6} (± {error:G2})");
        }
        Check("τ", tau, gTau);
        Check("ln τ", [.. tau.Select(v => Math.Log(v))], [.. gTau.Select(v => Math.Log(v))]);
        Check("σ²", sigma, gSigma);
        Check("σ⁴", [.. sigma.Select(v => v * v)], [.. gSigma.Select(v => v * v)]);
        for (int j = 0; j < 3; j++) Check($"a{j}", a[j], gA[j]);
        Check("a₁²", [.. a[1].Select(v => v * v)], [.. gA[1].Select(v => v * v)]);
        Check("a₁a₂", [.. a[1].Zip(a[2], (u, v) => u * v)], [.. gA[1].Zip(gA[2], (u, v) => u * v)]);
    }

    // Gibbs on (a, σ², τ) for a fixed structure under the ridge, level centring, σ² scaled and τ ~ Gamma(a_τ, rate a_τ), on
    // the explicit centred design X = [1, Z_c] with prior precision τD₁, D₁ = diag(1, d), d_j = ‖z_j‖²/N, about μ₀ = ȳe₀:
    //   a | σ², τ ~ N(A⁻¹(Xᵀy + τD₁μ₀), σ²A⁻¹), A = XᵀX + τD₁,
    //   σ² | a, τ ~ IG(c + (N + k)/2, cs² + ‖y − Xa‖²/2 + τ(a − μ₀)ᵀD₁(a − μ₀)/2),
    //   τ | a, σ² ~ Gamma(a_τ + k/2, rate a_τ + (a − μ₀)ᵀD₁(a − μ₀)/(2σ²)).
    // The coefficients are returned uncentred, in the composed draws' order.
    private static (double[] Tau, double[] Sigma2, double[][] A) JointGibbs(BassModel model, BassTerm[] structure, BassTermTable table,
                                                                            int sweeps, int seed)
    {
        int n = model.Resolved.Observations, k = structure.Length + 1;
        double[] y = [.. model.Data.Response.Values.ToArray()];
        double[][] z = CentredColumns(model, structure);
        var means = new double[structure.Length];
        for (int j = 0; j < structure.Length; j++)
        {
            var raw = new double[n];
            model.Basis.CaseColumn(structure[j], model.Inputs.Cases, raw);
            means[j] = raw.Average();
        }
        var x = new double[n, k];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = 1;
            for (int j = 0; j < structure.Length; j++) x[i, j + 1] = z[j][i];
        }
        var g = new double[k, k];
        var xty = new double[k];
        for (int r = 0; r < k; r++)
        {
            for (int i = 0; i < n; i++) xty[r] += x[i, r] * y[i];
            for (int c = 0; c < k; c++)
                for (int i = 0; i < n; i++) g[r, c] += x[i, r] * x[i, c];
        }
        var d = new double[k];
        for (int j = 0; j < k; j++) d[j] = g[j, j] / n;
        var mu = new double[k];
        mu[0] = model.Resolved.ResponseMean;
        InverseGamma sigmaPrior = model.Resolved.Sigma2!.Value;
        Gamma tauPrior = model.Resolved.Scale.Gamma!.Value;
        Assert.Equal(tauPrior.Shape, tauPrior.Rate);   // τ₀ = 1
        var rng = new Xoshiro256PlusPlus(seed);
        double tau = 1.0, sigma2 = 0.1;
        var a = new double[k];
        var w = new double[k];
        var factor = new DenseCholesky(k);
        var outTau = new double[sweeps];
        var outSigma = new double[sweeps];
        var outA = new double[k][];
        for (int j = 0; j < k; j++) outA[j] = new double[sweeps];
        int[] ids = [.. structure.Select(t => table.Intern(t))];
        int[] order = [.. Enumerable.Range(0, structure.Length).OrderBy(j => ids[j])];
        for (int it = -2_000; it < sweeps; it++)
        {
            var precision = (double[,])g.Clone();
            var rhs = new double[k];
            for (int j = 0; j < k; j++)
            {
                precision[j, j] += tau * d[j];
                rhs[j] = xty[j] + tau * d[j] * mu[j];
            }
            factor.Decompose(precision);
            double[] mean = factor.Solve(rhs);
            var normals = new double[k];
            for (int j = 0; j < k; j++) normals[j] = Normal.SampleStandard(rng);
            factor.SampleInnovation(normals, w);
            for (int j = 0; j < k; j++) a[j] = mean[j] + Math.Sqrt(sigma2) * w[j];
            double rss = 0, quad = 0;
            for (int i = 0; i < n; i++)
            {
                double fitted = 0;
                for (int j = 0; j < k; j++) fitted += x[i, j] * a[j];
                rss += (y[i] - fitted) * (y[i] - fitted);
            }
            for (int j = 0; j < k; j++) quad += d[j] * (a[j] - mu[j]) * (a[j] - mu[j]);
            sigma2 = InverseGamma.ShapeScale(sigmaPrior.Shape + 0.5 * (n + k), sigmaPrior.Scale + 0.5 * rss + 0.5 * tau * quad).Sample(rng);
            tau = Gamma.ShapeRate(tauPrior.Shape + 0.5 * k, tauPrior.Rate + 0.5 * quad / sigma2).Sample(rng);
            if (it < 0) continue;
            outTau[it] = tau;
            outSigma[it] = sigma2;
            double intercept = a[0];
            for (int q = 0; q < structure.Length; q++)
            {
                int j = order[q];
                outA[q + 1][it] = a[j + 1];
                intercept -= a[j + 1] * means[j];
            }
            outA[0][it] = intercept;
        }
        return (outTau, outSigma, outA);
    }

    public static TheoryData<string> FitCases => new()
    {
        "functional", "mixed", "held", "hyper-g", "improper hyper-g", "g, hyper-g/n", "tempered", "smc", "flat, jeffreys",
    };

    [Theory]
    [MemberData(nameof(FitCases))]
    public void Fits_UnderTheRidge_RunEveryView(string option)
    {
        BassData data = option switch
        {
            "functional" => BassRunTests.FunctionalData(),
            "mixed" => BassPriorSamplerTests.MixedData(),
            _ => BassSweepTests.Data(large: true),
        };
        BassCoefficientPrior coefficients = option switch
        {
            "held" => BassCoefficientPrior.Ridge(TauPrior.Held(2.0)),
            "hyper-g" => BassCoefficientPrior.Ridge(TauPrior.HyperG(3)),
            "improper hyper-g" => BassCoefficientPrior.Ridge(TauPrior.HyperG(1.5)),
            "g, hyper-g/n" => BassCoefficientPrior.G(TauPrior.HyperGOverN(3)),
            _ => BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()),
        };
        BassConfig config = new()
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 8 },
            Coefficients = coefficients,
            Sampling = new BassSampling { Chains = 2, MasterSeed = 5, BurnIn = 400, Samples = 120 },
        };
        config = option switch
        {
            "tempered" => config with { Sampling = config.Sampling with { Tempering = new TemperingSpec { ReferenceRung = true } } },
            "smc" => config with { Sampling = config.Sampling with { Smc = new SmcSpec { Particles = 200, Replicates = 2, MutationSteps = 2 } } },
            "flat, jeffreys" => config with { Intercept = GPriorCentre.Flat, Sigma2 = Sigma2Prior.Jeffreys },
            _ => config,
        };
        BassFit fit = Bass.Fit(config, data);
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        Assert.All(draws, d =>
        {
            Assert.True(d.Tau > 0 && double.IsFinite(d.Tau) && d.Sigma2 > 0 && double.IsFinite(d.Sigma2));
            Assert.All(d.Coefficients, v => Assert.True(double.IsFinite(v)));
        });
        if (option == "held") Assert.All(draws, d => Assert.Equal(2.0, d.Tau));
        else Assert.True(draws.Select(d => d.Tau).Distinct().Count() > 10, "τ is sampled");
        Assert.Contains(draws, d => d.Structure.TermIds.Length > 0);

        // Evidence where the run can give one: the ladder with its reference rung and SMC, under the proper priors; the
        // record carries the form and τ's prior on its unit information.
        bool evidence = option is "tempered" or "smc";
        Assert.Equal(evidence, fit.Evidence is not null);
        if (fit.Evidence is { } e)
        {
            Assert.True(double.IsFinite(e.LogEvidence));
            Assert.Equal(coefficients, e.Priors.Coefficients);
            Assert.Equal((coefficients.UnitObservations(fit.Resolved.Observations), false), (e.Priors.Scale.UnitObservations, e.Priors.Scale.IsHeld));
        }

        // Every view runs on the ridge's draws.
        BassSummary summary = fit.Summary();
        Assert.Equal(option != "held", summary.Posterior.Tau is not null);
        if (option == "smc") Assert.Throws<InvalidOperationException>(() => fit.Traces());   // a weighted population has no draw order
        else Assert.Equal(option != "held", fit.Traces().TauSampled);
        _ = fit.Metadata();
        PointwiseLogDensity retained = fit.PointwiseLogDensity(data);
        PointwiseLogDensity collapsed = fit.PointwiseLogDensity(data, content: PointwiseContent.CollapsedPredictive);
        foreach (PointwiseLogDensity density in new[] { retained, collapsed })
        {
            Assert.All(density.Values.Cast<double>(), v => Assert.True(double.IsFinite(v)));
            Assert.True(double.IsFinite(PredictiveChecks.PsisLoo(density).Elpd));
        }
        Assert.True(double.IsFinite(PredictiveChecks.Waic(retained).Elpd));   // WAIC reads the retained content only
        BassPrediction prediction = fit.Predict(data, PredictionCases.Of(data), new BassPredictOptions { MeanSurface = true, Seed = 3 });
        Assert.NotNull(prediction);
        _ = fit.TrainingFit(data);
        _ = fit.Sensitivity(data);
        _ = fit.Effects(data);
        Assert.Equal(option == "functional", fit.ResidualDependence is not null);
        _out.WriteLine($"{option}: {draws.Length} draws, mean M {draws.Average(d => d.Structure.TermIds.Length):F2}, E[τ] {draws.Average(d => d.Tau):G4}, " +
                       $"E[σ²] {draws.Average(d => d.Sigma2):G4}{(fit.Evidence is { } ev ? $", ln Z {ev.LogEvidence:F3} ± {ev.MonteCarloError:F3}" : "")}");
    }
}
