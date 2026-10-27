using System;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Standalone validation of the weighted-Normal marginal (DMGK eq. 6) and the end-to-end reversible-jump
/// stack on it. The marginal is the engine's spine; these tests pin the quadratic-form computation, the
/// weight path, the over-knotting penalty, and that carrier + moves + target + engine wire together.
/// </summary>
public sealed class WeightedNormalMarginalTests
{
    // Box–Muller standard normal from the project die.
    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = 1.0 - rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static double[] InSpan(double[,] z, double[] beta)
    {
        int m = z.GetLength(0), nu = z.GetLength(1);
        var y = new double[m];
        for (int i = 0; i < m; i++)
        {
            double s = 0.0;
            for (int j = 0; j < nu; j++) s += z[i, j] * beta[j];
            y[i] = s;
        }
        return y;
    }

    [Theory]
    [InlineData(GPriorCentre.Median)]
    [InlineData(GPriorCentre.Mean)]
    [InlineData(GPriorCentre.Zero)]
    public void Marginal_PerfectFit_ReducesToShrinkageResidual(GPriorCentre centre)
    {
        // When y lies exactly in span(Z), so does ỹ = y − Zμ₀, whose regression SS is then ỹᵀỹ, so a = ỹᵀỹ/(m+1) exactly:
        // ỹ = y − c with c the median or mean (the B-splines sum to one), and y under zero.
        var basis = new SplineBasis(degree: 1);
        var config = new KnotConfig(new[] { 0.5 });
        double[] x = { 0.1, 0.3, 0.5, 0.7, 0.9 };
        double[,] z = basis.Design(config, x);
        int m = x.Length, nu = z.GetLength(1);

        double[] beta = { 1.0, 2.0, -1.0 };   // length must equal ν = k + degree + 1 = 1 + 1 + 1
        Assert.Equal(nu, beta.Length);
        double[] y = InSpan(z, beta);

        double level = NormalCentreTests.Level(y, centre);
        double yy = 0.0;
        foreach (double v in y) yy += (v - level) * (v - level);
        double expectedA = yy / (m + 1);
        double expected = -0.5 * nu * Math.Log(m + 1) - 0.5 * m * Math.Log(expectedA);

        double actual = new WeightedNormalModel(centre, sigma2Prior: Sigma2Prior.Jeffreys).LogMarginalLikelihood(z, y, null);
        Assert.Equal(expected, actual, 6);

        // The default IG(g₁, g₁·ỹᵀỹ/m) integrates σ² in closed form: the exact fit's complete marginal.
        const double g1 = 1e-3;
        double g2 = g1 * yy / m;
        double complete = -0.5 * m * Math.Log(2.0 * Math.PI) - 0.5 * nu * Math.Log(m + 1) + g1 * Math.Log(g2)
                          + SpecialFunctions.LogGamma(g1 + 0.5 * m) - SpecialFunctions.LogGamma(g1) - (g1 + 0.5 * m) * Math.Log(g2 + 0.5 * expectedA);
        Assert.Equal(complete, new WeightedNormalModel(centre).LogMarginalLikelihood(z, y, null), 6);
    }

    [Fact]
    public void Marginal_ZeroResponse_IsRejected_WhileAnInterpolatingFitIsNot()
    {
        // a ≥ Q/(m+1), so the a ≤ 0 guard fires only for an identically zero response. A fit that interpolates
        // with m > ν keeps a = Q/(m+1) > 0 and a finite marginal; over-knotting (m ≤ ν) has its own rule.
        var basis = new SplineBasis(degree: 1);
        var config = new KnotConfig(new[] { 0.5 });
        double[] x = { 0.1, 0.3, 0.5, 0.7, 0.9 };
        double[,] z = basis.Design(config, x);
        var model = new WeightedNormalModel();

        Assert.Equal(double.NegativeInfinity, model.LogMarginalLikelihood(z, new double[x.Length], null));
        Assert.True(double.IsFinite(model.LogMarginalLikelihood(z, InSpan(z, new[] { 1.0, 2.0, -1.0 }), null)));
    }

    [Fact]
    public void Marginal_UnitWeights_MatchNullWeights()
    {
        var basis = new SplineBasis(degree: 3);
        var config = new KnotConfig(new[] { 0.4, 0.7 });
        var rng = new Xoshiro256PlusPlus(seed: 1);
        int m = 40;
        var x = new double[m];
        var y = new double[m];
        for (int i = 0; i < m; i++) { x[i] = (i + 0.5) / m; y[i] = Math.Sin(6 * x[i]) + 0.1 * Gaussian(rng); }

        double[,] z = basis.Design(config, x);
        var ones = new double[m];
        Array.Fill(ones, 1.0);

        var model = new WeightedNormalModel();
        Assert.Equal(model.LogMarginalLikelihood(z, y, null), model.LogMarginalLikelihood(z, y, ones), 10);
    }

    [Fact]
    public void Posterior_PrefersTrueModel_OverSpuriousKnots()
    {
        // Data from a 2-knot cubic with tiny noise; the (m+1)^(-ν/2) penalty must beat a 9-knot bloat.
        (double trueScore, double bloatScore) = TrueAndBloated(noiseSeed: 12345);

        Assert.True(trueScore > bloatScore, $"true {trueScore} should beat bloated {bloatScore}");
    }

    // The test's computation at any noise seed, which the seed sweep reruns (ComparativeSweep): the log posterior of
    // the true two knots and of nine spurious ones.
    internal static (double TrueScore, double BloatScore) TrueAndBloated(int noiseSeed)
    {
        var basis = new SplineBasis(degree: 3);
        var rng = new Xoshiro256PlusPlus(seed: noiseSeed);
        int m = 80;
        var x = new double[m];
        for (int i = 0; i < m; i++) x[i] = (i + 0.5) / m;

        var trueConfig = new KnotConfig(new[] { 0.3, 0.6 });
        double[,] zTrue = basis.Design(trueConfig, x);
        int nuTrue = zTrue.GetLength(1);
        var beta = new double[nuTrue];
        for (int j = 0; j < nuTrue; j++) beta[j] = 1.0 + 0.5 * j * (j % 2 == 0 ? 1 : -1);

        double[] y = InSpan(zTrue, beta);
        for (int i = 0; i < m; i++) y[i] += 0.01 * Gaussian(rng);

        var target = new SplineTarget(basis, new WeightedNormalModel(), new UniformPrior(50), x, y);
        double trueScore = target.LogPosterior(trueConfig);
        var bloat = new KnotConfig(new[] { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 });
        return (trueScore, target.LogPosterior(bloat));
    }

    [Fact]
    public void Chain_RunsEndToEnd_WithValidAcceptanceRate()
    {
        var basis = new SplineBasis(degree: 3);
        var noise = new Xoshiro256PlusPlus(seed: 7);
        int m = 60;
        var x = new double[m];
        for (int i = 0; i < m; i++) x[i] = (i + 0.5) / m;

        var trueConfig = new KnotConfig(new[] { 0.4 });
        double[,] zt = basis.Design(trueConfig, x);
        int nut = zt.GetLength(1);
        var beta = new double[nut];
        for (int j = 0; j < nut; j++) beta[j] = 1.0 + 0.3 * j;
        double[] y = InSpan(zt, beta);
        for (int i = 0; i < m; i++) y[i] += 0.05 * Gaussian(noise);

        var target = new SplineTarget(basis, new WeightedNormalModel(), new PoissonPrior(3.0), x, y);
        var kernel = new UniformKernel();
        var moves = new IRjMove<KnotConfig>[]
        {
            new KnotBirthMove(kernel), new KnotDeathMove(kernel), new KnotRelocateMove(kernel),
        };
        var chain = new ReversibleJumpChain<KnotConfig>(
            moves, target, new KnotConfig(Array.Empty<double>()), new Xoshiro256PlusPlus(seed: 99));

        for (int s = 0; s < 500; s++) chain.Step();

        Assert.True(chain.Attempts > 0);
        double accept = (double)chain.Accepted / chain.Attempts;
        Assert.InRange(accept, 0.0, 1.0);
        Assert.True(chain.Current.Count >= 0);
    }

    [Fact]
    public void Marginal_RejectsSingularDesign()
    {
        // Rank-deficient design (column 1 ≡ column 0) ⇒ singular Gram. The marginal must reject it (−∞) rather
        // than reward the degenerate "perfect" fit (the old Epsilon-guard turned this into a large positive).
        var rng = new Xoshiro256PlusPlus(seed: 3);
        const int m = 8, nu = 3;
        var z = new double[m, nu];
        var y = new double[m];
        for (int i = 0; i < m; i++)
        {
            double t = (i + 0.5) / m;
            z[i, 0] = t;
            z[i, 1] = t;          // identical to column 0 ⇒ rank-deficient
            z[i, 2] = t * t;
            y[i] = Math.Sin(t) + 0.05 * Gaussian(rng);
        }

        double lml = new WeightedNormalModel().LogMarginalLikelihood(z, y, null);
        Assert.True(double.IsNegativeInfinity(lml), $"singular design should be rejected (−∞), got {lml}");
    }
}
