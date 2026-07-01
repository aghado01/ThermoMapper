using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
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

    [Fact]
    public void Marginal_PerfectFit_ReducesToShrinkageResidual()
    {
        // When y lies exactly in span(Z), the regression SS equals yᵀy, so a = yᵀy/(m+1) exactly.
        var basis = new SplineBasis(degree: 1);
        var config = new KnotConfig(new[] { 0.5 });
        double[] x = { 0.1, 0.3, 0.5, 0.7, 0.9 };
        double[,] z = basis.Design(config, x);
        int m = x.Length, nu = z.GetLength(1);

        double[] beta = { 1.0, 2.0, -1.0 };   // length must equal ν = k + degree + 1 = 1 + 1 + 1
        Assert.Equal(nu, beta.Length);
        double[] y = InSpan(z, beta);

        double yy = 0.0;
        foreach (double v in y) yy += v * v;
        double expectedA = yy / (m + 1);
        double expected = -0.5 * nu * Math.Log(m + 1) - 0.5 * m * Math.Log(expectedA);

        double actual = new WeightedNormalModel().LogMarginalLikelihood(z, y, null);
        Assert.Equal(expected, actual, 6);
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
        var basis = new SplineBasis(degree: 3);
        var rng = new Xoshiro256PlusPlus(seed: 12345);
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
        double bloatScore = target.LogPosterior(bloat);

        Assert.True(trueScore > bloatScore, $"true {trueScore} should beat bloated {bloatScore}");
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
