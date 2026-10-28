using System;
using System.Linq;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The model-free consumers of a pointwise density (OQ17). The generalized Pareto fit recovers known shapes; PSIS-LOO
/// matches exact leave-one-out on a conjugate Gaussian linear model, where each held-out predictive is closed-form; WAIC
/// follows its definition and agrees with exact leave-one-out to O(1/n); unit partitions validate.
/// </summary>
public sealed class PredictiveChecksTests
{
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    [Theory]
    [InlineData(0.3)]
    [InlineData(0.0)]
    [InlineData(0.7)]
    public void GeneralizedParetoFit_RecoversTheShape(double shape)
    {
        var rng = new Xoshiro256PlusPlus(seed: 5);
        const int n = 4000;
        const double scale = 2.0;
        var x = new double[n];
        for (int i = 0; i < n; i++)
        {
            double u = rng.NextDouble();
            x[i] = shape == 0.0 ? -scale * Math.Log(1.0 - u) : scale * (Math.Pow(1.0 - u, -shape) - 1.0) / shape;
        }
        Array.Sort(x);
        (double k, double sigma) = PredictiveChecks.FitGeneralizedPareto(x);
        Assert.InRange(k, shape - 0.07, shape + 0.07);
        Assert.InRange(sigma, scale * 0.9, scale * 1.1);
    }

    [Fact]
    public void ParetoTails_ReportHowEachUnitWasTreated()
    {
        // posterior's gpdfit guard: the grid's scale is the tail's first quartile, so a tail whose lowest quarter sits at
        // its minimum has no fit.
        (double k, double sigma) = PredictiveChecks.FitGeneralizedPareto([0.0, 0.0, 0.0, 0.0, 0.1, 0.4, 0.9, 1.6]);
        Assert.True(double.IsNaN(k) && double.IsNaN(sigma));

        // Unit 0: 355 of 400 log ratios share one value and 45 exceed it, so the tail of 60 starts with 15 at the cutoff:
        // no fit, k̂ = +∞, as in loo. Unit 1 varies smoothly and is fitted. Twenty draws give a tail of 4: too short.
        var rng = new Xoshiro256PlusPlus(seed: 3);
        var values = new double[400, 2];
        for (int d = 0; d < 400; d++)
        {
            values[d, 0] = d < 355 ? -1.0 : -1.0 - 0.02 * (d - 354);
            values[d, 1] = -0.5 * Math.Pow(0.3 - 0.1 * Gaussian(rng), 2);
        }
        PredictiveCheck check = PredictiveChecks.PsisLoo(new PointwiseLogDensity(values, UnitPartition.EachObservation(2), []));
        Assert.Equal([ParetoTail.NoFit, ParetoTail.Fitted], check.ParetoTails);
        Assert.Equal(double.PositiveInfinity, check.ParetoK![0]);
        Assert.True(double.IsNaN(check.MonteCarloError), "an unfitted tail is above the threshold, so the total's error is withheld");

        var few = new double[20, 1];
        for (int d = 0; d < 20; d++) few[d, 0] = -0.5 * Math.Pow(Gaussian(rng), 2);
        PredictiveCheck shortTail = PredictiveChecks.PsisLoo(new PointwiseLogDensity(few, UnitPartition.EachObservation(1), []));
        Assert.Equal([ParetoTail.TooShort], shortTail.ParetoTails);
        Assert.Equal(double.PositiveInfinity, shortTail.ParetoK![0]);
    }

    // y = Xβ + ε with ε ~ N(0, σ²) known and β ~ N(0, τ²I): the posterior is N(m, V), V = (XᵀX/σ² + I/τ²)⁻¹, m = VXᵀy/σ²,
    // and leaving observation i out gives the predictive N(xᵢᵀm₋ᵢ, σ² + xᵢᵀV₋ᵢxᵢ) from the precision without xᵢxᵢᵀ/σ².
    private const double Sigma2 = 0.25, Tau2 = 4.0;

    private static (double[,] X, double[] Y) Linear(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[n, 3];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            double t = (i + 0.5) / n;
            x[i, 0] = 1.0; x[i, 1] = t; x[i, 2] = t * t;
            y[i] = 0.5 - 1.2 * t + 2.0 * t * t + Math.Sqrt(Sigma2) * Gaussian(rng) + (i == 7 ? 1.5 : 0.0);   // one outlier
        }
        return (x, y);
    }

    private static double[,] Inverse3(double[,] a)
    {
        double c00 = a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1], c01 = a[1, 2] * a[2, 0] - a[1, 0] * a[2, 2], c02 = a[1, 0] * a[2, 1] - a[1, 1] * a[2, 0];
        double det = a[0, 0] * c00 + a[0, 1] * c01 + a[0, 2] * c02;
        return new[,]
        {
            { c00 / det, (a[0, 2] * a[2, 1] - a[0, 1] * a[2, 2]) / det, (a[0, 1] * a[1, 2] - a[0, 2] * a[1, 1]) / det },
            { c01 / det, (a[0, 0] * a[2, 2] - a[0, 2] * a[2, 0]) / det, (a[0, 2] * a[1, 0] - a[0, 0] * a[1, 2]) / det },
            { c02 / det, (a[0, 1] * a[2, 0] - a[0, 0] * a[2, 1]) / det, (a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0]) / det },
        };
    }

    private static (double[,] Precision, double[] Linear) Posterior(double[,] x, double[] y, int skip)
    {
        var precision = new double[3, 3];
        var linear = new double[3];
        for (int j = 0; j < 3; j++) precision[j, j] = 1.0 / Tau2;
        for (int i = 0; i < y.Length; i++)
        {
            if (i == skip) continue;
            for (int a = 0; a < 3; a++)
            {
                linear[a] += x[i, a] * y[i] / Sigma2;
                for (int b = 0; b < 3; b++) precision[a, b] += x[i, a] * x[i, b] / Sigma2;
            }
        }
        return (precision, linear);
    }

    private static (PointwiseLogDensity Density, double[] ExactLoo) LinearCase(int n, int draws, int seed)
    {
        (double[,] x, double[] y) = Linear(n, seed);
        (double[,] precision, double[] linear) = Posterior(x, y, -1);
        double[,] v = Inverse3(precision);
        var m = new double[3];
        for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) m[a] += v[a, b] * linear[b];
        // Cholesky of V for exact posterior draws.
        var l = new double[3, 3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = v[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                l[i, j] = i == j ? Math.Sqrt(s) : s / l[j, j];
            }
        var rng = new Xoshiro256PlusPlus(seed + 100);
        var ll = new double[draws, n];
        for (int d = 0; d < draws; d++)
        {
            double z0 = Gaussian(rng), z1 = Gaussian(rng), z2 = Gaussian(rng);
            var beta = new double[3];
            for (int a = 0; a < 3; a++) beta[a] = m[a] + l[a, 0] * z0 + (a >= 1 ? l[a, 1] * z1 : 0.0) + (a == 2 ? l[a, 2] * z2 : 0.0);
            for (int i = 0; i < n; i++)
            {
                double mu = x[i, 0] * beta[0] + x[i, 1] * beta[1] + x[i, 2] * beta[2], r = y[i] - mu;
                ll[d, i] = -0.5 * Math.Log(2.0 * Math.PI * Sigma2) - 0.5 * r * r / Sigma2;
            }
        }
        var exact = new double[n];
        for (int i = 0; i < n; i++)
        {
            (double[,] pi, double[] li) = Posterior(x, y, i);
            double[,] vi = Inverse3(pi);
            double mean = 0.0, var = Sigma2;
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++)
                {
                    mean += x[i, a] * vi[a, b] * li[b];
                    var += x[i, a] * vi[a, b] * x[i, b];
                }
            exact[i] = -0.5 * Math.Log(2.0 * Math.PI * var) - 0.5 * (y[i] - mean) * (y[i] - mean) / var;
        }
        var report = new[] { new LatentReport("β", LatentScope.Global, LatentTreatment.Conditioned) };
        return (new PointwiseLogDensity(ll, UnitPartition.EachObservation(n), report), exact);
    }

    [Fact]
    public void PsisLoo_MatchesExactLeaveOneOut_OnAConjugateLinearModel()
    {
        (PointwiseLogDensity density, double[] exact) = LinearCase(n: 30, draws: 8000, seed: 3);
        PredictiveCheck loo = PredictiveChecks.PsisLoo(density);
        for (int i = 0; i < exact.Length; i++)
        {
            Assert.True(loo.ParetoK![i] < 0.7, $"observation {i}: k̂ {loo.ParetoK[i]}");
            Assert.True(Math.Abs(loo.PerUnit[i] - exact[i]) < 0.03, $"observation {i}: PSIS {loo.PerUnit[i]} against exact {exact[i]}");
        }
        Assert.InRange(loo.Elpd - exact.Sum(), -0.1, 0.1);
        Assert.True(loo.ParetoK![7] > loo.ParetoK.Where((_, i) => i != 7).Average(), "the outlier should have the heaviest tail");
        Assert.InRange(loo.EffectiveParameters, 2.0, 5.0);   // three coefficients, a little more for the outlier
    }

    [Fact]
    public void Waic_FollowsItsDefinition_AndAgreesWithExactLooToOrderOneOverN()
    {
        // The definition on a hand-checkable matrix: two draws, two units.
        var tiny = new PointwiseLogDensity(new[,] { { -1.0, -2.0 }, { -3.0, -2.0 } }, UnitPartition.EachObservation(2),
                                           Array.Empty<LatentReport>());
        PredictiveCheck w = PredictiveChecks.Waic(tiny);
        double lppd0 = Math.Log(0.5 * (Math.Exp(-1.0) + Math.Exp(-3.0))), p0 = 2.0;   // sample variance of {−1, −3}
        Assert.Equal(lppd0 - p0, w.PerUnit[0], 12);
        Assert.Equal(-2.0, w.PerUnit[1], 12);
        Assert.Equal(p0, w.EffectiveParameters, 12);

        (PointwiseLogDensity density, double[] exact) = LinearCase(n: 30, draws: 8000, seed: 3);
        Assert.InRange(PredictiveChecks.Waic(density).Elpd - exact.Sum(), -0.4, 0.4);
    }

    [Fact]
    public void UnitPartitions_Validate()
    {
        var p = new UnitPartition(new[] { 1, 0, 1, 2 }, new[] { 0, 0, 1 });
        Assert.Equal(3, p.Units);
        Assert.Equal(new[] { 0, 2 }, p.Members[1]);
        Assert.Equal(new[] { 0, 0, 1 }, p.StratumOf);
        Assert.Equal(new[] { 0, 0, 0, 1, 1, 1 }, UnitPartition.ConsecutiveBlocks(2, 3).UnitOf);
        Assert.Throws<ArgumentException>(() => new UnitPartition(new[] { 0, 2 }));            // unit 1 is empty
        Assert.Throws<ArgumentException>(() => new UnitPartition(new[] { 0, 1 }, new[] { 0 }));   // one stratum per unit
        Assert.Throws<ArgumentOutOfRangeException>(() => new UnitPartition(new[] { -1 }));
    }
}
