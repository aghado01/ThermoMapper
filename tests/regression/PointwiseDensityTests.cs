using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Changepoint;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Each engine's pointwise density (OQ17) against an independent computation of the rule it implements: a latent is
/// integrated out of a unit's term when everything that depends on it lies inside the unit, and conditioned on
/// otherwise. BARS's Normal terms and its robust weights' marginals (against numerical integration over the weight's
/// prior), HierarchicalPSpline's curve units (against an explicit dense multivariate normal) and observation units, and the
/// changepoint draws' segment levels.
/// </summary>
public sealed class PointwiseDensityTests
{
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double LogNormal(double y, double mean, double variance)
        => -0.5 * Math.Log(2.0 * Math.PI * variance) - 0.5 * (y - mean) * (y - mean) / variance;

    private static (double[] X, double[] Y) Bump(int n, int seed, bool outliers)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Exp(-40.0 * (x[i] - 0.45) * (x[i] - 0.45)) + 0.15 * Gaussian(rng) + (outliers && i % 11 == 0 ? 1.5 : 0.0);
        }
        return (x, y);
    }

    private static BarsEnsemble Bars(IObservationModel model, IWeightResampler? resampler = null) =>
        new(new SplineBasis(3), model, new PoissonPrior(5.0), new LocalBetaKernel(50.0), resampler: resampler);

    [Fact]
    public void Bars_NormalTerms_AreTheGaussianAtEachDraw_AndUnitsSumThem()
    {
        (double[] x, double[] y) = Bump(60, 3, outliers: false);
        BarsEnsemble ensemble = Bars(new WeightedNormalModel());
        BarsResult r = ensemble.Run(x, y, x, chains: 2, masterSeed: 4, burn: 200, samples: 100);
        PointwiseLogDensity each = ensemble.PointwiseLogDensity(r, x, y);
        Assert.Equal(200, each.Draws);
        Assert.Equal(60, each.Units.Units);

        var basis = new SplineBasis(3);
        BarsDraw[] draws = r.Draws.SelectMany(c => c).ToArray();
        foreach (int d in new[] { 0, 57, 199 })
        {
            double[,] z = basis.Design(draws[d].Knots, x);
            foreach (int i in new[] { 0, 30, 59 })
            {
                double eta = 0.0;
                for (int j = 0; j < draws[d].Coefficients.Length; j++) eta += z[i, j] * draws[d].Coefficients[j];
                Assert.Equal(LogNormal(y[i], eta, draws[d].Sigma2), each.Values[d, i], 12);
            }
        }

        PointwiseLogDensity blocks = ensemble.PointwiseLogDensity(r, x, y, UnitPartition.ConsecutiveBlocks(6, 10));
        for (int d = 0; d < 200; d += 37)
            for (int u = 0; u < 6; u++)
            {
                double sum = 0.0;
                for (int i = 10 * u; i < 10 * u + 10; i++) sum += each.Values[d, i];
                Assert.Equal(sum, blocks.Values[d, u], 10);
            }
        Assert.All(each.Latents, l => Assert.Equal(LatentTreatment.Conditioned, l.Treatment));
        Assert.Contains(each.Latents, l => l.Name == "noise variance σ²");
    }

    [Theory]
    [InlineData(0.1)]
    [InlineData(1.0)]
    [InlineData(5.0)]
    public void StudentTMarginal_IntegratesTheWeightUnderItsGammaPrior(double residual)
    {
        // ∫ N(r; 0, σ²/w) Gamma(w; ν/2, rate ν/2) dw by the trapezoid rule over ln w, against the closed form.
        const double nu = 4.0, sigma2 = 0.3;
        double sum = 0.0;
        const int steps = 200_000;
        double lo = -25.0, h = 50.0 / steps;
        for (int k = 0; k <= steps; k++)
        {
            double u = lo + k * h, w = Math.Exp(u);
            double logPrior = 0.5 * nu * Math.Log(0.5 * nu) - SpecialFunctions.LogGamma(0.5 * nu) + (0.5 * nu - 1.0) * u - 0.5 * nu * w;
            double f = Math.Exp(LogNormal(residual, 0.0, sigma2 / w) + logPrior + u);   // + u: dw = w du
            sum += (k == 0 || k == steps ? 0.5 : 1.0) * f;
        }
        Assert.Equal(Math.Log(sum * h), new StudentTWeights(nu).MarginalLogDensity(residual, sigma2), 9);
    }

    [Fact]
    public void ContaminatedMarginal_IsTheNormalMixture_AndIntegratesToOne()
    {
        var weights = new ContaminatedNormalWeights(outlierFraction: 0.1, inflation: 5.0);
        const double sigma2 = 0.3;
        foreach (double r in new[] { 0.0, 0.7, 4.0 })
        {
            double mixture = 0.9 * Math.Exp(LogNormal(r, 0.0, sigma2)) + 0.1 * Math.Exp(LogNormal(r, 0.0, 25.0 * sigma2));
            Assert.Equal(Math.Log(mixture), weights.MarginalLogDensity(r, sigma2), 12);
        }
        double total = 0.0, h = 1e-3;
        for (double r = -40.0; r <= 40.0; r += h) total += Math.Exp(weights.MarginalLogDensity(r, sigma2)) * h;
        Assert.Equal(1.0, total, 6);
    }

    [Fact]
    public void Bars_RobustWeights_AreIntegrated_AndHelpOnOutliers()
    {
        (double[] x, double[] y) = Bump(60, 5, outliers: true);
        BarsEnsemble plain = Bars(new WeightedNormalModel());
        BarsEnsemble robust = Bars(new WeightedNormalModel(), new StudentTWeights(4.0));
        BarsResult rp = plain.Run(x, y, x, chains: 2, masterSeed: 2, burn: 400, samples: 400);
        BarsResult rr = robust.Run(x, y, x, chains: 2, masterSeed: 2, burn: 400, samples: 400);
        PointwiseLogDensity dp = plain.PointwiseLogDensity(rp, x, y), dr = robust.PointwiseLogDensity(rr, x, y);

        BarsDraw first = rr.Draws[0][0];
        double[,] z = new SplineBasis(3).Design(first.Knots, x);
        double eta = 0.0;
        for (int j = 0; j < first.Coefficients.Length; j++) eta += z[11, j] * first.Coefficients[j];
        Assert.Equal(new StudentTWeights(4.0).MarginalLogDensity(y[11] - eta, first.Sigma2), dr.Values[0, 11], 12);
        Assert.Contains(dr.Latents, l => l.Name == "robust weights wᵢ" && l.Scope == LatentScope.Observation && l.Treatment == LatentTreatment.Integrated);

        // The Student-t error model predicts the held-out outliers far better than the Normal one.
        Assert.True(PredictiveChecks.PsisLoo(dr).Elpd > PredictiveChecks.PsisLoo(dp).Elpd + 5.0);
    }

    [Fact]
    public void Bars_PoissonTerms_AreThePoissonLogMass_AndGibbsLossDeclines()
    {
        var rng = new Xoshiro256PlusPlus(seed: 7);
        int n = 50;
        var x = new double[n];
        var counts = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            double lambda = 3.0 + 5.0 * Math.Exp(-40.0 * (x[i] - 0.5) * (x[i] - 0.5)), p = Math.Exp(-lambda), s = p, u = rng.NextDouble();
            int k = 0;
            while (u > s && k < 200) { k++; p *= lambda / k; s += p; }
            counts[i] = k;
        }
        BarsEnsemble poisson = Bars(new PoissonModel());
        BarsResult r = poisson.Run(x, counts, x, chains: 2, masterSeed: 3, burn: 200, samples: 50);
        PointwiseLogDensity density = poisson.PointwiseLogDensity(r, x, counts);
        BarsDraw d = r.Draws[1][7];
        double[,] z = new SplineBasis(3).Design(d.Knots, x);
        double eta = 0.0;
        for (int j = 0; j < d.Coefficients.Length; j++) eta += z[20, j] * d.Coefficients[j];
        Assert.Equal(counts[20] * eta - Math.Exp(eta) - SpecialFunctions.LogGamma(counts[20] + 1.0), density.Values[57, 20], 10);

        BarsEnsemble gibbs = Bars(new GibbsLossModel());
        BarsResult rg = gibbs.Run(x, counts, x, chains: 1, masterSeed: 3, burn: 50, samples: 10);
        Assert.Throws<NotSupportedException>(() => gibbs.PointwiseLogDensity(rg, x, counts));
    }

    // Two replicates of eight points on five B-spline coefficients.
    private static (double[,] Z, double[][] Y) Replicates()
    {
        const int n = 8;
        var x = new double[n];
        for (int i = 0; i < n; i++) x[i] = (i + 0.5) / n;
        var rng = new Xoshiro256PlusPlus(seed: 303);
        var y = new double[2][];
        for (int r = 0; r < 2; r++)
        {
            y[r] = new double[n];
            for (int i = 0; i < n; i++) y[r][i] = 1.5 + Math.Sin(2.0 * Math.PI * x[i]) + 0.3 * r * Math.Cos(Math.PI * x[i]) + 0.3 * Gaussian(rng);
        }
        return (new SplineBasis(3).Design(new KnotConfig(new[] { 0.5 }), x), y);
    }

    // log N(y; Zβ_0, σ²I + τ²ZZᵀ) through an explicit n × n Cholesky factor.
    private static double DenseMarginal(double[,] z, double[] y, double[] beta0, double sigma2, double tau2)
    {
        int n = y.Length, nu = z.GetLength(1);
        var c = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                double s = i == j ? sigma2 : 0.0;
                for (int k = 0; k < nu; k++) s += tau2 * z[i, k] * z[j, k];
                c[i, j] = s;
            }
        var l = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = c[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                l[i, j] = i == j ? Math.Sqrt(s) : s / l[j, j];
            }
        var w = new double[n];
        double logDet = 0.0, quad = 0.0;
        for (int i = 0; i < n; i++)
        {
            double mean = 0.0;
            for (int k = 0; k < nu; k++) mean += z[i, k] * beta0[k];
            double s = y[i] - mean;
            for (int k = 0; k < i; k++) s -= l[i, k] * w[k];
            w[i] = s / l[i, i];
            logDet += 2.0 * Math.Log(l[i, i]);
            quad += w[i] * w[i];
        }
        return -0.5 * (n * Math.Log(2.0 * Math.PI) + logDet + quad);
    }

    [Fact]
    public void Hierarchical_CurveUnits_IntegrateTheReplicate_ObservationUnitsConditionOnIt()
    {
        (double[,] z, double[][] y) = Replicates();
        var model = new HierarchicalPSpline(z, y, new DifferencePenalty(2), 2.0, 0.2);
        HierarchicalResult h = model.Run(chains: 2, masterSeed: 5, burn: 200, samples: 100);
        HierarchicalDraw[] draws = h.PosteriorDraws.SelectMany(c => c).ToArray();

        PointwiseLogDensity curves = model.PointwiseLogDensity(h, UnitPartition.ConsecutiveBlocks(2, 8));
        PointwiseLogDensity points = model.PointwiseLogDensity(h);
        foreach (int d in new[] { 0, 99, 150 })
            for (int r = 0; r < 2; r++)
            {
                HierarchicalDraw draw = draws[d];
                Assert.Equal(DenseMarginal(z, y[r], draw.Population, draw.NoiseVariance, draw.DeviationVariance), curves.Values[d, r], 9);
                for (int i = 0; i < 8; i++)
                {
                    double f = 0.0;
                    for (int k = 0; k < z.GetLength(1); k++) f += z[i, k] * draw.Replicates![r][k];
                    Assert.Equal(LogNormal(y[r][i], f, draw.NoiseVariance), points.Values[d, 8 * r + i], 12);
                }
            }
        Assert.Equal(LatentTreatment.Integrated, curves.Latents.Single(l => l.Name == "replicate curves β_r").Treatment);
        Assert.Equal(LatentTreatment.Conditioned, points.Latents.Single(l => l.Name == "replicate curves β_r").Treatment);

        // Replicate 0 whole, replicate 1 point by point: mixed.
        var mixed = new int[16];
        for (int i = 8; i < 16; i++) mixed[i] = i - 7;
        Assert.Equal(LatentTreatment.Mixed, model.PointwiseLogDensity(h, new UnitPartition(mixed)).Latents[0].Treatment);

        // Without replicate draws, curve units still work, since they integrate the curves out; split units cannot.
        HierarchicalResult lean = model.Run(chains: 1, masterSeed: 5, burn: 100, samples: 50, retainReplicates: false);
        Assert.Equal(50, model.PointwiseLogDensity(lean, UnitPartition.ConsecutiveBlocks(2, 8)).Draws);
        Assert.Throws<InvalidOperationException>(() => model.PointwiseLogDensity(lean));
    }

    [Fact]
    public void Changepoint_Levels_IntegratedForWholeSegments_ConditionedOtherwise()
    {
        var rng = new Xoshiro256PlusPlus(seed: 9);
        int n = 30;
        var y = new double[n];
        for (int i = 0; i < n; i++) y[i] = (i < 12 ? 0.0 : 2.0) + 0.4 * Gaussian(rng);
        (double s2, double nu, double rho2) = ExactChangepoint.DefaultHyperparameters(y);
        ChangepointResult result = ExactChangepoint.Infer(y, s2, nu, rho2, maxSegments: 6);
        IReadOnlyList<ChangepointDraw> draws = result.Posterior.Sample(200, seed: 4);

        // One unit holding everything integrates every level: the term is the sum of the draw's segment evidences.
        PointwiseLogDensity all = result.Posterior.PointwiseLogDensity(draws, new UnitPartition(new int[n]));
        PointwiseLogDensity each = result.Posterior.PointwiseLogDensity(draws);
        for (int d = 0; d < draws.Count; d += 17)
        {
            ChangepointDraw draw = draws[d];
            double expected = 0.0;
            int start = 0;
            foreach (int end in draw.SegmentEnds)
            {
                expected += result.Posterior.SegmentLogEvidence(start, end);
                start = end + 1;
            }
            Assert.Equal(expected, all.Values[d, 0], 10);

            // Point units: a level conditioned on, unless its segment is that one point.
            start = 0;
            for (int s = 0; s < draw.SegmentEnds.Length; s++)
            {
                for (int i = start; i <= draw.SegmentEnds[s]; i++)
                {
                    double term = start == draw.SegmentEnds[s] ? result.Posterior.SegmentLogEvidence(i, i) : LogNormal(y[i], draw.Levels[s], s2);
                    Assert.Equal(term, each.Values[d, i], 10);
                }
                start = draw.SegmentEnds[s] + 1;
            }
        }
        Assert.Equal(LatentTreatment.Integrated, all.Latents.Single(l => l.Name == "segment levels").Treatment);
        Assert.Contains(each.Latents.Single(l => l.Name == "segment levels").Treatment, new[] { LatentTreatment.Conditioned, LatentTreatment.Mixed });
    }
}
