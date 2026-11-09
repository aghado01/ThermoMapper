using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E34's references (the multi-output brief, slice 0): the collapsed marginal of K outputs on one basis, with A integrated against
/// its matrix-normal prior and Σ against its inverse Wishart (<see cref="MatrixTReference"/>), checked by routes that share none of
/// its algebra.
/// <list type="bullet">
/// <item>A's integral against the dense Gaussian on vec(Y) (<see cref="DenseVecReference"/>), with the Kronecker products built entry
/// by entry: an asymmetric Σ at K = 2 separates Σ ⊗ Λ₀⁻¹ from Λ₀⁻¹ ⊗ Σ, and a seeded sweep takes far ranges and near-singular Σ.</item>
/// <item>Σ's integral against quadrature at K = 2 and Monte Carlo over the prior's draws at K = 3; Σ's posterior against the
/// quadrature's moments.</item>
/// <item>The matrix-t density in the data's n × n row scale, untempered.</item>
/// <item>K = 1 against the scalar collapsed g-prior under decision 2's reduction, IW(2c, 2c·s²) = IG(c, c·s²), for the g-form, the
/// normalized ridge and a pencil, under level centring, the zero centre and the flat intercept, refusals included.</item>
/// <item>Equivariance: Y → YL with decision 2's prior shifts every structure's marginal by the same −(nb − f) ln|det L|.</item>
/// <item>The sources: Brown, Vannucci and Fearn's (1998) eqs. (20)–(22) over every subset of five columns, and Huang and Wand's
/// (2013) Property 3, the correlation law of an inverse Wishart with a diagonal scale.</item>
/// </list>
/// </summary>
public sealed class MatrixTReferenceTests(ITestOutputHelper output)
{
    private const int Seed = 20261108;

    // ── A's integral: the dense vec form ────────────────────────────────────────────────────────────────────────────────────

    public static TheoryData<double, bool> ChosenSettings => new()
    {
        { 1.0, false },
        { 0.4, false },
        { 1.0, true },
        { 0.4, true },
    };

    [Theory]
    [MemberData(nameof(ChosenSettings))]
    public void Collapse_IsTheDenseVecForm_UnderColumnStacking_AndNotTheSwappedKronecker(double b, bool flat)
    {
        // K = p = 2, so Λ₀⁻¹ ⊗ Σ has Σ ⊗ Λ₀⁻¹'s shape and only the order tells them apart. Σ has unequal variances and a
        // correlation, Λ₀ is no multiple of it, and the outputs' scales differ.
        (double[,] y, double[,] x) = TwoOutputData(flat);
        double[,] lambda = { { 1.7, 0.45 }, { 0.45, 0.6 } };
        double[,] a0 = { { 0.3, -1.1 }, { 0.2, 0.8 } };
        double[,] sigma = { { 0.7, 0.9 }, { 0.9, 2.6 } };
        var c = new MultiOutputCase(y, x, lambda, a0, b, flat);

        double collapse = MatrixTReference.LogLikelihoodGiven(MatrixTReference.Collapse(c), sigma);
        double dense = DenseVecReference.LogLikelihoodGiven(c, sigma);
        double swapped = DenseVecReference.LogLikelihoodGiven(c, sigma, KroneckerOrder.Swapped);
        output.WriteLine($"b {b}, flat {flat}: collapse {collapse:R}, dense {dense:R}, swapped {swapped:R}");
        Close(dense, collapse, 1e-12);
        Assert.True(Math.Abs(swapped - collapse) > 0.05, $"the swapped Kronecker reads {swapped:R} against {collapse:R}");
    }

    [Fact]
    public void Collapse_IsTheDenseVecForm_OverASeededSweep_FarRangesAndNearSingularSigma()
    {
        // Two regimes. Moderate: the response's units from 10⁻⁶ to 10⁶, b down to 10⁻⁶, cond(Σ) ≤ 10², a prior and columns within two
        // decades, judged at a fixed tolerance. Far: cond(Σ) to 10⁸, Λ₀ over eight decades and the columns over six, judged against
        // κ̂ = cond(Σ)·(1 + b·tr XΛ₀⁻¹Xᵀ), the condition the dense Cholesky of Cov(vec Y) works at, up to κ̂ = 10¹³, past which the
        // dense form reads nothing (a case beyond it is drawn again). The collapse reads p × p and K × K matrices and stays exact there.
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 1)[0]);
        double worstModerate = 0.0, worstFar = 0.0, largestKappa = 0.0;
        int far = 0, redrawn = 0;
        const int cases = 500;
        for (int index = 0; index < cases; index++)
        {
            bool farCase = index % 2 == 1;
            int k = 1 + rng.NextInt(4), p = rng.NextInt(5), n = 2 + rng.NextInt(8);
            bool flat = rng.NextDouble() < 0.4;
            double b = rng.NextDouble() switch { < 0.2 => 1.0, < 0.4 => Math.Pow(10.0, -2.0 - 4.0 * rng.NextDouble()), _ => 0.05 + 0.95 * rng.NextDouble() };
            double scale = Math.Pow(10.0, -6.0 + 12.0 * rng.NextDouble());
            double logCondition = farCase ? 8.0 * rng.NextDouble() : 2.0 * rng.NextDouble();
            double[,] sigma = RandomSpd(rng, k, logCondition, scale * scale);
            double[,] lambda = farCase
                ? RandomSpd(rng, p, 4.0 * rng.NextDouble(), Math.Pow(10.0, -4.0 + 8.0 * rng.NextDouble()))
                : RandomSpd(rng, p, 2.0 * rng.NextDouble(), Math.Pow(10.0, -1.0 + 2.0 * rng.NextDouble()));
            double[,] x = RandomMatrix(rng, n, p, 1.0);
            for (int j = 0; j < p; j++)
            {
                double s = farCase ? Math.Pow(10.0, -3.0 + 6.0 * rng.NextDouble()) : Math.Pow(10.0, -1.0 + 2.0 * rng.NextDouble());
                for (int i = 0; i < n; i++) x[i, j] *= s;
            }
            double[,]? a0 = rng.NextDouble() < 0.5 ? null : RandomMatrix(rng, p, k, scale);
            double[,] y = RandomMatrix(rng, n, k, scale);
            var c = new MultiOutputCase(y, x, lambda, a0, b, flat);
            double[,] spread = ReferenceAlgebra.Multiply(ReferenceAlgebra.Multiply(x, ReferenceAlgebra.Inverse(lambda)), ReferenceAlgebra.Transpose(x));
            double trace = 0.0;
            for (int i = 0; i < n; i++) trace += spread[i, i];
            double kappa = Math.Pow(10.0, logCondition) * (1.0 + b * trace);
            if (kappa > 1e13) { redrawn++; index--; continue; }

            double collapse = MatrixTReference.LogLikelihoodGiven(MatrixTReference.Collapse(c), sigma);
            double dense = DenseVecReference.LogLikelihoodGiven(c, sigma);
            double error = Math.Abs(collapse - dense) / Math.Max(1.0, Math.Abs(dense));
            if (farCase)
            {
                far++;
                worstFar = Math.Max(worstFar, error / kappa);
                largestKappa = Math.Max(largestKappa, kappa);
            }
            else worstModerate = Math.Max(worstModerate, error);
        }
        output.WriteLine($"{cases - far} moderate cases, worst relative gap {worstModerate:E2}; {far} far cases to κ̂ {largestKappa:E1}, "
                         + $"worst relative gap / κ̂ {worstFar:E2}; {redrawn} drawn again past κ̂ = 1e13");
        Assert.True(worstModerate < 1e-11, $"worst relative gap {worstModerate:E2}");
        Assert.True(worstFar < 1e-13, $"worst relative gap over κ̂ {worstFar:E2}");
    }

    [Fact]
    public void Marginal_IsTheMatrixTDensity_Untempered()
    {
        // The intercept inside the prior at b = 1: Y ~ MT(XA₀, I + XΛ₀⁻¹Xᵀ, Ψ₀, ν₀), read in the n × n row scale.
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 2)[1]);
        double worst = 0.0;
        for (int index = 0; index < 200; index++)
        {
            int k = 1 + rng.NextInt(4), p = rng.NextInt(5), n = 1 + rng.NextInt(9);
            double scale = Math.Pow(10.0, -4.0 + 8.0 * rng.NextDouble());
            double[,] lambda = RandomSpd(rng, p, 3.0 * rng.NextDouble(), Math.Pow(10.0, -3.0 + 6.0 * rng.NextDouble()));
            double[,] x = RandomMatrix(rng, n, p, 1.0);
            double[,] a0 = RandomMatrix(rng, p, k, scale);
            double[,] y = RandomMatrix(rng, n, k, scale);
            var prior = new ReferenceSigmaPrior(k - 1 + 0.01 + 20.0 * rng.NextDouble(), RandomSpd(rng, k, 3.0 * rng.NextDouble(), scale * scale));
            double collapse = MatrixTReference.LogMarginal(MatrixTReference.Collapse(new MultiOutputCase(y, x, lambda, a0, 1.0)), prior);

            double[,] location = ReferenceAlgebra.Multiply(x, a0);
            double[,] rowScale = ReferenceAlgebra.Multiply(ReferenceAlgebra.Multiply(x, ReferenceAlgebra.Inverse(lambda)), ReferenceAlgebra.Transpose(x));
            for (int i = 0; i < n; i++) rowScale[i, i] += 1.0;
            ReferenceAlgebra.Symmetrize(rowScale);
            double density = MatrixTReference.MatrixTLogDensity(y, location, rowScale, prior.Scale!, prior.Nu);
            worst = Math.Max(worst, Math.Abs(collapse - density) / Math.Max(1.0, Math.Abs(density)));
        }
        output.WriteLine($"200 cases: worst relative gap {worst:E2}");
        Assert.True(worst < 1e-11, $"worst relative gap {worst:E2}");
    }

    // ── Σ's integral: quadrature at K = 2, Monte Carlo at K = 3 ─────────────────────────────────────────────────────────────

    public static TheoryData<string> QuadratureCases => new() { "informative", "tempered", "decision 2, level centring", "flat intercept" };

    [Theory]
    [MemberData(nameof(QuadratureCases))]
    public void Marginal_MatchesQuadratureOverSigma_AtTwoOutputs(string name)
    {
        // The integrand is A's integral at Σ (the collapse's, which the dense vec form checks above) times the inverse Wishart's
        // density, whose own normalization the quadrature checks below; so this checks Σ's integral alone: Γ₂, the powers of |Σ|
        // and ln|Ψ₀ + R|. It is evaluated through the 2 × 2 kernel, checked here against those two at points over the box. The first
        // pass's reference, (Ψ₀ + b·n·Cov(Y))⁻¹, and every box centre are read from the data, not from the closed form.
        (MultiOutputCase c, ReferenceSigmaPrior prior) = QuadratureCase(name);
        MatrixTCollapse collapse = MatrixTReference.Collapse(c);
        double[,] psi = prior.Scale!, merged = ReferenceAlgebra.Sum(collapse.R, psi);
        double constant = -collapse.Count * Math.Log(2.0 * Math.PI) + collapse.LogDeterminantTerm
                          + 0.5 * prior.Nu * Math.Log(psi[0, 0] * psi[1, 1] - psi[1, 0] * psi[1, 0]) - prior.Nu * Math.Log(2.0)
                          - MatrixTReference.LogMultivariateGamma(2, 0.5 * prior.Nu);
        double LogF(double s11, double s21, double s22) => constant + SigmaQuadrature.Kernel(s11, s21, s22, merged, collapse.Count + prior.Nu + 3.0);
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 8)[7]);
        for (int i = 0; i < 50; i++)
        {
            double[,] sigma = RandomSpd(rng, 2, 3.0 * rng.NextDouble(), Math.Pow(10.0, -2.0 + 3.0 * rng.NextDouble()));
            Close(MatrixTReference.LogLikelihoodGiven(collapse, sigma) + ReferenceInverseWishart.LogDensity(sigma, prior.Nu, psi),
                  LogF(sigma[0, 0], sigma[1, 0], sigma[1, 1]), 1e-12);
        }

        double[,] spread = Covariance(c.Y);
        for (int a = 0; a < 2; a++)
            for (int d = 0; d < 2; d++) spread[a, d] = psi[a, d] + c.InverseTemperature * c.Rows * spread[a, d];
        (SigmaQuadrature.Result fine, SigmaQuadrature.Result check, SigmaQuadrature.Result first) = SigmaQuadrature.Adaptive(LogF, ReferenceAlgebra.Inverse(spread));
        double closed = MatrixTReference.LogMarginal(collapse, prior);
        (double nuPost, double[,] scalePost) = MatrixTReference.SigmaPosterior(collapse, prior);
        output.WriteLine($"{name}: closed {closed:R}, quadrature {fine.LogIntegral:R} (coarser {check.LogIntegral:R}; first pass {first.LogIntegral:R}), "
                         + $"faces {fine.BoundaryGap:F1} below the peak");
        output.WriteLine($"  E[Σ]: posterior IW ({scalePost[0, 0] / (nuPost - 3.0):R}, {scalePost[1, 0] / (nuPost - 3.0):R}, {scalePost[1, 1] / (nuPost - 3.0):R}), "
                         + $"quadrature ({fine.Mean[0, 0]:R}, {fine.Mean[1, 0]:R}, {fine.Mean[1, 1]:R})");

        Assert.True(fine.BoundaryGap > 40.0, $"the box's faces sit {fine.BoundaryGap:F1} below the peak");
        Assert.True(Math.Abs(fine.LogIntegral - check.LogIntegral) < 1e-12, $"the rule has not converged: {fine.LogIntegral:R} against {check.LogIntegral:R}");
        Assert.True(Math.Abs(fine.LogIntegral - closed) < 1e-11, $"closed {closed:R}, quadrature {fine.LogIntegral:R}");
        for (int a = 0; a < 2; a++)
            for (int d = 0; d < 2; d++)
                Close(scalePost[a, d] / (nuPost - 3.0), fine.Mean[a, d], 1e-10);
    }

    [Fact]
    public void InverseWishart_IntegratesToOne_AtTwoOutputs_AndIsTheInverseGammaAtOne()
    {
        // Checks the quadrature's Jacobian with the density's normalization, ν = 4.5 putting a heavy left tail in ln A₂₂, χ²(3.5).
        double[,] psi = { { 1.3, -0.7 }, { -0.7, 2.2 } };
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 9)[8]);
        foreach (double nu in new[] { 4.5, 11.0 })
        {
            double constant = 0.5 * nu * Math.Log(psi[0, 0] * psi[1, 1] - psi[1, 0] * psi[1, 0]) - nu * Math.Log(2.0) - MatrixTReference.LogMultivariateGamma(2, 0.5 * nu);
            double LogF(double s11, double s21, double s22) => constant + SigmaQuadrature.Kernel(s11, s21, s22, psi, nu + 3.0);
            for (int i = 0; i < 20; i++)
            {
                double[,] sigma = RandomSpd(rng, 2, 3.0 * rng.NextDouble(), Math.Pow(10.0, -1.0 + 2.0 * rng.NextDouble()));
                Close(ReferenceInverseWishart.LogDensity(sigma, nu, psi), LogF(sigma[0, 0], sigma[1, 0], sigma[1, 1]), 1e-13);
            }
            (SigmaQuadrature.Result fine, SigmaQuadrature.Result check, _) = SigmaQuadrature.Adaptive(LogF, ReferenceAlgebra.Inverse(new[,] { { 2.0, 0.9 }, { 0.9, 3.5 } }));
            output.WriteLine($"ν {nu}: ln ∫ IW = {fine.LogIntegral:E2} (coarser {check.LogIntegral:E2}), faces {fine.BoundaryGap:F1} below the peak; "
                             + $"E[Σ] ({fine.Mean[0, 0]:R}, {fine.Mean[1, 0]:R}, {fine.Mean[1, 1]:R}) against Ψ/(ν − 3)");
            Assert.True(fine.BoundaryGap > 40.0, $"faces {fine.BoundaryGap:F1} below the peak");
            Assert.True(Math.Abs(fine.LogIntegral) < 1e-12, $"ln ∫ IW = {fine.LogIntegral:R}");
        }
        foreach ((double nu, double psi1, double s2) in new[] { (0.002, 0.004, 0.3), (3.0, 2.5, 1.7), (40.0, 1e-6, 2e-8) })
            Close(InverseGamma.ShapeScale(0.5 * nu, 0.5 * psi1).LogPdf(s2), ReferenceInverseWishart.LogDensity(new[,] { { s2 } }, nu, new[,] { { psi1 } }), 1e-13);
    }

    public static TheoryData<double> MonteCarloTemperatures => new() { 1.0, 0.5 };

    [Theory]
    [MemberData(nameof(MonteCarloTemperatures))]
    public void Marginal_MatchesMonteCarloOverTheSigmaPrior_AtThreeOutputs(double b)
    {
        // Σ drawn from IW(ν₀, Ψ₀) by Bartlett's decomposition, so no Γ₃ enters this side; the mean of p_b(Y | Σ), A integrated at
        // each draw by the collapse (checked against the dense vec form on the first draws here, and over K ≤ 4 above), is the
        // marginal. Y comes from a surface and noise outside the model; Ψ₀ sits near its residual scales, so the average settles.
        const int n = 6, k = 3, draws = 300_000;
        double[,] x = new double[n, 2];
        double[,] y = new double[n, k];
        var noise = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 3)[2]);
        for (int i = 0; i < n; i++)
        {
            double t = -1.0 + 2.0 * i / (n - 1);
            x[i, 0] = 1.0;
            x[i, 1] = t;
            double shared = 0.4 * (noise.NextDouble() - 0.5);
            y[i, 0] = Math.Sin(2.0 * t) + shared + 0.2 * (noise.NextDouble() - 0.5);
            y[i, 1] = 3.0 * t * t - 1.0 - 2.0 * shared + 0.5 * (noise.NextDouble() - 0.5);
            y[i, 2] = 0.05 * (Math.Exp(t) - 1.25) + 0.03 * (noise.NextDouble() - 0.5);
        }
        double[,] lambda = { { 0.05, 0.01 }, { 0.01, 0.08 } };
        double[,] psi = { { 0.56, -0.48, 0.004 }, { -0.48, 11.2, 0.008 }, { 0.004, 0.008, 0.00096 } };
        var prior = new ReferenceSigmaPrior(12.0, psi);
        var c = new MultiOutputCase(y, x, lambda, null, b);
        MatrixTCollapse collapse = MatrixTReference.Collapse(c);
        double closed = MatrixTReference.LogMarginal(collapse, prior);

        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 11)[b == 1.0 ? 9 : 10]);
        double[,] root = ReferenceInverseWishart.SamplingRoot(psi);
        var logs = new double[draws];
        double max = double.NegativeInfinity;
        for (int s = 0; s < draws; s++)
        {
            double[,] sigma = ReferenceInverseWishart.SampleFromRoot(prior.Nu, root, rng);
            logs[s] = MatrixTReference.LogLikelihoodGiven(collapse, sigma);
            if (s < 100) Close(DenseVecReference.LogLikelihoodGiven(c, sigma), logs[s], 1e-11);
            max = Math.Max(max, logs[s]);
        }
        double sum = 0.0, sumSquares = 0.0;
        foreach (double l in logs)
        {
            double w = Math.Exp(l - max);
            sum += w;
            sumSquares += w * w;
        }
        double mean = sum / draws, variance = (sumSquares / draws - mean * mean) * draws / (draws - 1.0);
        double estimate = max + Math.Log(mean), se = Math.Sqrt(variance / draws) / mean;
        output.WriteLine($"b {b}: closed {closed:R}, Monte Carlo {estimate:R} ± {se:E2} ({draws} draws), z {(estimate - closed) / se:F2}");
        Assert.True(se < 0.01, $"the average has not settled: se {se:E2}");
        Assert.True(Math.Abs(estimate - closed) < 4.0 * se, $"closed {closed:R}, Monte Carlo {estimate:R} ± {se:E2}");
    }

    [Fact]
    public void InverseWishartSampler_HasItsMean()
    {
        const int draws = 200_000;
        double[,] psi = { { 0.6, -0.4, 0.01 }, { -0.4, 2.4, -0.02 }, { 0.01, -0.02, 0.004 } };
        const double nu = 9.0;
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 5)[4]);
        var sum = new double[3, 3];
        var sumSquares = new double[3, 3];
        for (int s = 0; s < draws; s++)
        {
            double[,] sigma = ReferenceInverseWishart.Sample(nu, psi, rng);
            for (int a = 0; a < 3; a++)
                for (int d = 0; d < 3; d++) { sum[a, d] += sigma[a, d]; sumSquares[a, d] += sigma[a, d] * sigma[a, d]; }
        }
        for (int a = 0; a < 3; a++)
            for (int d = 0; d < 3; d++)
            {
                double mean = sum[a, d] / draws, se = Math.Sqrt((sumSquares[a, d] / draws - mean * mean) / draws), expected = psi[a, d] / (nu - 4.0);
                Assert.True(Math.Abs(mean - expected) < 4.5 * se, $"E[Σ]({a}, {d}) = {mean:R} ± {se:E2} against {expected:R}");
            }
    }

    // ── K = 1: the scalar collapsed g-prior under decision 2's reduction ───────────────────────────────────────────────────

    private enum ScalarForm { Gram, Ridge, Pencil }

    [Fact]
    public void Marginal_AtOneOutput_IsTheScalarMarginal_UnderDecisionTwosReduction()
    {
        // IW(2c, 2c·s²) at K = 1 is IG(c, c·s²), the scalar default (D03), with s² the mean square about the centre level: ȳ under
        // level centring and the flat intercept, 0 under the zero centre. Where one side refuses (σ²'s posterior shape c + ν/2 not
        // positive), the other must too.
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 6)[5]);
        GPriorCentre[] centres = [GPriorCentre.Mean, GPriorCentre.Zero, GPriorCentre.Flat];
        double worst = 0.0;
        int compared = 0, refusedBoth = 0;
        for (int index = 0; index < 600; index++)
        {
            int n = 2 + rng.NextInt(11), m = rng.NextInt(Math.Min(4, n - 1) + 1);
            GPriorCentre centre = centres[rng.NextInt(3)];
            var form = (ScalarForm)rng.NextInt(3);
            double scale = Math.Pow(10.0, -6.0 + 12.0 * rng.NextDouble()), tau = Math.Pow(10.0, -6.0 + 10.0 * rng.NextDouble());
            double b = rng.NextDouble() switch { < 0.15 => 1.0, < 0.3 => centre == GPriorCentre.Flat ? 1.0 : 0.0, < 0.45 => Math.Pow(10.0, -4.0 * rng.NextDouble()) / n, _ => rng.NextDouble() };
            bool jeffreys = b > 0.0 && rng.NextDouble() < 0.2;
            double c = rng.NextDouble() < 0.5 ? 1e-3 : 0.5;

            double[,] xc = ReferenceAlgebra.Centred(RandomMatrix(rng, n, m, Math.Pow(10.0, -2.0 + 4.0 * rng.NextDouble())));
            var y = new double[n];
            double offset = scale * (-3.0 + 6.0 * rng.NextDouble());
            for (int i = 0; i < n; i++) y[i] = offset + scale * Normal.SampleStandard(rng);
            double[,]? pencilPrior = form == ScalarForm.Pencil ? RandomSpd(rng, m + (centre == GPriorCentre.Flat ? 0 : 1), 2.0 * rng.NextDouble(), 1.0) : null;
            if (pencilPrior is not null && centre != GPriorCentre.Flat) { pencilPrior[0, 0] *= n; }   // the intercept's direction near its Gram entry

            (double scalar, bool scalarRefused) = Attempt(() => ScalarLogMarginal(y, xc, centre, form, pencilPrior, tau, b, c, jeffreys));
            (double reference, bool referenceRefused) = Attempt(() => ReferenceLogMarginal(y, xc, centre, form, pencilPrior, tau, b, c, jeffreys));
            Assert.True(scalarRefused == referenceRefused,
                        $"case {index}: scalar {(scalarRefused ? "refused" : scalar.ToString("R"))}, reference {(referenceRefused ? "refused" : reference.ToString("R"))} ({centre}, {form}, n {n}, b {b:R}, c {c})");
            if (scalarRefused) { refusedBoth++; continue; }
            compared++;
            double error = Math.Abs(scalar - reference) / Math.Max(1.0, Math.Abs(scalar));
            Assert.True(error < 1e-11, $"case {index}: scalar {scalar:R}, reference {reference:R} ({centre}, {form}, n {n}, M {m}, τ {tau:E2}, b {b:R}, c {c}, Jeffreys {jeffreys})");
            worst = Math.Max(worst, error);
        }
        output.WriteLine($"{compared} compared, worst relative gap {worst:E2}; {refusedBoth} refused by both");
        Assert.True(refusedBoth > 0, "the sweep reached no refusal");
    }

    // ── Equivariance, b = 0 and refusals ───────────────────────────────────────────────────────────────────────────────────

    public static TheoryData<GPriorCentre, double> EquivarianceSettings => new()
    {
        { GPriorCentre.Mean, 1.0 }, { GPriorCentre.Mean, 0.3 }, { GPriorCentre.Zero, 1.0 }, { GPriorCentre.Flat, 1.0 }, { GPriorCentre.Flat, 0.3 },
    };

    [Theory]
    [MemberData(nameof(EquivarianceSettings))]
    public void Marginal_UnderDecisionTwo_ShiftsEveryStructureAlike_WhenTheOutputsAreMixed(GPriorCentre centre, double b)
    {
        // Y → YL takes S_Y to LᵀS_YL, so ln p_b(YL) = ln p_b(Y) − (nb − f) ln|det L| for every structure, and the structure posterior
        // is unchanged. L rescales by 10⁻³ to 10³, mixes, permutes and flips a sign.
        const int n = 9, k = 3;
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 7)[6]);
        double[,] y = RandomMatrix(rng, n, k, 1.0);
        for (int i = 0; i < n; i++) { y[i, 0] += 4.0; y[i, 2] = 0.5 * y[i, 0] + 0.1 * y[i, 2]; }
        double[,] l = { { 0.0, 1e3, 0.2 }, { -0.7, 0.4, 1.1 }, { 1e-3, 0.0, 0.9 } };
        double logDetL = Math.Log(Math.Abs(Determinant3(l)));
        double[,] yl = ReferenceAlgebra.Multiply(y, l);
        double[,] all = RandomMatrix(rng, n, 3, 1.0);
        double expected = -(n * b - (centre == GPriorCentre.Flat ? 1.0 : 0.0)) * logDetL;
        foreach (int m in new[] { 1, 3 })
        {
            var xc = new double[n, m];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++) xc[i, j] = all[i, j];
            xc = ReferenceAlgebra.Centred(xc);
            double shift = GFormLogMarginal(yl, xc, centre, 0.4, b) - GFormLogMarginal(y, xc, centre, 0.4, b);
            output.WriteLine($"M {m}: shift {shift:R}, −(nb − f) ln|det L| = {expected:R}");
            Assert.True(Math.Abs(shift - expected) < 1e-9 * Math.Max(1.0, Math.Abs(expected)), $"M {m}: shift {shift:R} against {expected:R}");
        }
    }

    [Fact]
    public void Marginal_IsZeroAtBZero_AndRefusesAnImproperSigmaPosterior()
    {
        (double[,] y, double[,] x) = TwoOutputData(false);
        double[,] lambda = { { 1.7, 0.45 }, { 0.45, 0.6 } };
        double[,] a0 = { { 0.3, -1.1 }, { 0.2, 0.8 } };
        var proper = new ReferenceSigmaPrior(1.002, new[,] { { 0.01, 0.002 }, { 0.002, 0.03 } });
        Assert.Equal(0.0, MatrixTReference.LogMarginal(MatrixTReference.Collapse(new MultiOutputCase(y, x, lambda, a0, 0.0)), proper));
        Assert.Equal(0.0, MatrixTReference.LogLikelihoodGiven(MatrixTReference.Collapse(new MultiOutputCase(y, x, lambda, a0, 0.0)), proper.Scale!));

        // Decision 3: proper at every b under a proper prior with the intercept inside; under an improper ν₀ a hot rung fails.
        Assert.True(double.IsFinite(MatrixTReference.LogMarginal(MatrixTReference.Collapse(new MultiOutputCase(y, x, lambda, a0, 1e-9)), proper)));
        int n = y.GetLength(0);
        var improper = new ReferenceSigmaPrior(0.0, null);
        Assert.Throws<InvalidOperationException>(() => MatrixTReference.LogMarginal(MatrixTReference.Collapse(new MultiOutputCase(y, x, lambda, a0, 0.9 / n)), improper));
        Assert.True(double.IsFinite(MatrixTReference.LogMarginal(MatrixTReference.Collapse(new MultiOutputCase(y, x, lambda, a0, 1.1 / n)), improper)));

        // The flat intercept spends one observation: ν₀ + nb − 1 > K − 1 fails at small b even under decision 2's proper ν₀ = K − 1 + 2c,
        // and its integral diverges at b = 0.
        (double[,] yf, double[,] xf) = TwoOutputData(true);
        Assert.Throws<InvalidOperationException>(() => MatrixTReference.LogMarginal(MatrixTReference.Collapse(new MultiOutputCase(yf, xf, lambda, null, 0.99 / n, true)), proper));
        Assert.True(double.IsFinite(MatrixTReference.LogMarginal(MatrixTReference.Collapse(new MultiOutputCase(yf, xf, lambda, null, 1.0 / n, true)), proper)));
        Assert.Throws<ArgumentOutOfRangeException>(() => MatrixTReference.Collapse(new MultiOutputCase(yf, xf, lambda, null, 0.0, true)));
        Assert.Throws<ArgumentException>(() => MatrixTReference.LogMarginal(MatrixTReference.Collapse(new MultiOutputCase(y, x, lambda, a0, 1.0)),
                                                                            new ReferenceSigmaPrior(0.9, proper.Scale)));
    }

    // ── The sources: Brown, Vannucci and Fearn (1998); Huang and Wand (2013) ─────────────────────────────────────────────

    [Fact]
    public void StructureDifferences_AreBrownVannucciFearnsEquation20()
    {
        // Their model, (4)–(7): Y − 1αᵀ − XB ∼ N(I_n, Σ) on centred X and Y, α ∼ N(α₀, hΣ) with h → ∞, B ∼ N(H_γ, Σ) by rows, and
        // Σ ∼ IW(δ; Q) in Dawid's shape, ν₀ = δ + q − 1. Their (20): g(γ) = (|H_γ||K_γ|)^(−q/2)·|Q_γ|^(−(n + δ + q − 1)/2), K_γ = X_γᵀX_γ + H_γ⁻¹,
        // and (21) Q_γ = Q + YᵀY − YᵀX_γK_γ⁻¹X_γᵀY, taken here as written. On centred columns level centring keeps α apart from B, so
        // ln p − ln g is one constant over γ, and Q_γ is Ψ₀ + R. Their g-prior (22), H_γ = c(X_γᵀX_γ)⁻¹, gives |H_γ||K_γ| = (c + 1)^(p_γ) and
        // Q_γ = Q + c/(c + 1)·(YᵀY − YᵀX_γ(X_γᵀX_γ)⁻¹X_γᵀY) + YᵀY/(c + 1).
        const int n = 12, q = 3, p = 5;
        const double delta = 3.0, c = 30.0, v1 = 2.5;
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 13)[12]);
        double[,] x = ReferenceAlgebra.Centred(RandomMatrix(rng, n, p, 1.0));
        var raw = new double[n, q];
        for (int i = 0; i < n; i++)
        {
            double shared = Normal.SampleStandard(rng);
            raw[i, 0] = Math.Sin(1.5 * x[i, 0]) + 0.4 * x[i, 2] * x[i, 2] + 0.3 * shared;
            raw[i, 1] = -0.8 * x[i, 0] + Math.Abs(x[i, 3]) + 0.5 * shared + 0.2 * Normal.SampleStandard(rng);
            raw[i, 2] = 0.05 * Math.Exp(x[i, 1]) + 0.02 * Normal.SampleStandard(rng);
        }
        double[,] y = ReferenceAlgebra.Centred(raw);
        double[,] qPrior = { { 0.5, 0.0, 0.0 }, { 0.0, 0.5, 0.0 }, { 0.0, 0.0, 0.5 } };
        var prior = new ReferenceSigmaPrior(delta + q - 1.0, qPrior);
        double[,] yy = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(y), y);

        double worstConstant = 0.0, worstForm = 0.0, worstShortcut = 0.0;
        foreach (bool gPrior in new[] { false, true })
        {
            double? constant = null;
            for (int mask = 0; mask < 1 << p; mask++)
            {
                int[] columns = [.. Enumerable.Range(0, p).Where(j => (mask >> j & 1) == 1)];
                int pg = columns.Length;
                var xg = new double[n, pg];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < pg; j++) xg[i, j] = x[i, columns[j]];
                double[,] gram = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(xg), xg);
                var hInverse = new double[pg, pg];
                for (int i = 0; i < pg; i++)
                    for (int j = 0; j < pg; j++) hInverse[i, j] = gPrior ? gram[i, j] / c : i == j ? 1.0 / v1 : 0.0;
                double[,] k = ReferenceAlgebra.Sum(gram, hInverse);
                double logHK = ReferenceAlgebra.LogDeterminant(ReferenceAlgebra.Cholesky(k)!) - ReferenceAlgebra.LogDeterminant(ReferenceAlgebra.Cholesky(hInverse)!);
                double[,] xy = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(xg), y);
                double[,] fitted = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(xy), ReferenceAlgebra.Solve(ReferenceAlgebra.Cholesky(k)!, xy));
                var qg = new double[q, q];
                for (int a = 0; a < q; a++)
                    for (int d = 0; d < q; d++) qg[a, d] = qPrior[a, d] + yy[a, d] - fitted[a, d];
                ReferenceAlgebra.Symmetrize(qg);
                double logG = -0.5 * q * logHK - 0.5 * (n + delta + q - 1.0) * ReferenceAlgebra.LogDeterminant(ReferenceAlgebra.Cholesky(qg)!);

                if (gPrior)
                {
                    worstShortcut = Math.Max(worstShortcut, Math.Abs(logHK - pg * Math.Log(c + 1.0)));
                    double[,] ols = pg == 0 ? new double[q, q] : ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(xy), ReferenceAlgebra.Solve(ReferenceAlgebra.Cholesky(gram)!, xy));
                    for (int a = 0; a < q; a++)
                        for (int d = 0; d < q; d++)
                        {
                            double shortcut = qPrior[a, d] + c / (c + 1.0) * (yy[a, d] - ols[a, d]) + yy[a, d] / (c + 1.0);
                            worstShortcut = Math.Max(worstShortcut, Math.Abs(shortcut - qg[a, d]) / Math.Max(1.0, Math.Abs(qg[a, d])));
                        }
                }

                // The reference: the intercept inside the prior at precision 1·Σ⁻¹ (any precision: on centred Y its factor is one constant).
                var design = new double[n, pg + 1];
                var lambda = new double[pg + 1, pg + 1];
                lambda[0, 0] = 1.0;
                for (int i = 0; i < n; i++)
                {
                    design[i, 0] = 1.0;
                    for (int j = 0; j < pg; j++) design[i, j + 1] = xg[i, j];
                }
                for (int i = 0; i < pg; i++)
                    for (int j = 0; j < pg; j++) lambda[i + 1, j + 1] = hInverse[i, j];
                MatrixTCollapse collapse = MatrixTReference.Collapse(new MultiOutputCase(y, design, lambda, null, 1.0));
                double gap = MatrixTReference.LogMarginal(collapse, prior) - logG;
                constant ??= gap;
                worstConstant = Math.Max(worstConstant, Math.Abs(gap - constant.Value));
                (_, double[,] posteriorScale) = MatrixTReference.SigmaPosterior(collapse, prior);
                for (int a = 0; a < q; a++)
                    for (int d = 0; d < q; d++)
                        worstForm = Math.Max(worstForm, Math.Abs(posteriorScale[a, d] - qg[a, d]) / Math.Max(1.0, Math.Abs(qg[a, d])));
            }
        }
        output.WriteLine($"2 × 32 structures: ln p − ln g varies by {worstConstant:E2}; Ψ₀ + R against Q_γ {worstForm:E2}; (22)'s shortcuts {worstShortcut:E2}");
        Assert.True(worstConstant < 1e-10, $"ln p − ln g varies by {worstConstant:E2}");
        Assert.True(worstForm < 1e-12, $"Ψ₀ + R against Q_γ {worstForm:E2}");
        Assert.True(worstShortcut < 1e-12, $"(22)'s shortcuts {worstShortcut:E2}");
    }

    public static TheoryData<int, double, int> HuangWandShapes => new() { { 2, 0.5, 13 }, { 2, 2.0, 14 }, { 3, 2.0, 15 }, { 3, 6.0, 16 } };

    [Theory]
    [MemberData(nameof(HuangWandShapes))]
    public void InverseWishartDraws_GiveHuangAndWandsCorrelationLaw(int k, double nu, int stream)
    {
        // Huang and Wand (2013), Property 3: under Σ | a ∼ IW(ν + K − 1, 2ν diag(1/a)) each correlation has density ∝ (1 − ρ²)^(ν/2 − 1),
        // uniform at ν = 2. A correlation is unchanged by a diagonal rescaling, so the law holds at each a: ρ² ∼ Beta(½, ν/2), with
        // E[ρ²] = 1/(1 + ν) and E[ρ⁴] = 3/((1 + ν)(3 + ν)). Decision 2's ν₀ = K − 1 + 2c is this form at ν = 2c.
        const int draws = 200_000;
        var scale = new double[k, k];
        double[] a = [0.3, 4.0, 1.7];
        for (int i = 0; i < k; i++) scale[i, i] = 2.0 * nu / a[i];
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 17)[stream]);
        double[,] root = ReferenceInverseWishart.SamplingRoot(scale);
        double sum2 = 0.0, sumSq2 = 0.0, sum4 = 0.0, sumSq4 = 0.0;
        for (int s = 0; s < draws; s++)
        {
            double[,] sigma = ReferenceInverseWishart.SampleFromRoot(nu + k - 1.0, root, rng);
            int i = s % (k - 1), j = k - 1;   // ρ₁ₖ, and ρ₂ₖ at K = 3: each correlation's law is the same
            double rho2 = sigma[i, j] * sigma[i, j] / (sigma[i, i] * sigma[j, j]), rho4 = rho2 * rho2;
            sum2 += rho2; sumSq2 += rho2 * rho2; sum4 += rho4; sumSq4 += rho4 * rho4;
        }
        double m2 = sum2 / draws, m4 = sum4 / draws;
        double se2 = Math.Sqrt((sumSq2 / draws - m2 * m2) / draws), se4 = Math.Sqrt((sumSq4 / draws - m4 * m4) / draws);
        double e2 = 1.0 / (1.0 + nu), e4 = 3.0 / ((1.0 + nu) * (3.0 + nu));
        output.WriteLine($"K {k}, ν {nu}: E[ρ²] {m2:F5} ± {se2:E1} against {e2:F5}; E[ρ⁴] {m4:F5} ± {se4:E1} against {e4:F5}");
        Assert.True(Math.Abs(m2 - e2) < 4.5 * se2, $"E[ρ²] {m2:R} ± {se2:E2} against {e2:R}");
        Assert.True(Math.Abs(m4 - e4) < 4.5 * se4, $"E[ρ⁴] {m4:R} ± {se4:E2} against {e4:R}");
    }

    // ── Cases ──────────────────────────────────────────────────────────────────────────────────────────────────────────────

    // Two outputs on unequal scales over n = 7; X holds the intercept and a covariate, or, under the flat intercept, two covariates.
    private static (double[,] Y, double[,] X) TwoOutputData(bool flat)
    {
        const int n = 7;
        var y = new double[n, 2];
        var x = new double[n, 2];
        for (int i = 0; i < n; i++)
        {
            double t = -1.0 + 2.0 * i / (n - 1);
            x[i, 0] = flat ? t * t : 1.0;
            x[i, 1] = t;
            y[i, 0] = 0.4 + Math.Cos(1.3 * t) + 0.15 * Math.Sin(7.0 * i);
            y[i, 1] = 3.0 * t - 0.8 * y[i, 0] + 0.6 * Math.Cos(5.0 * i);
        }
        return (y, x);
    }

    private static (MultiOutputCase Case, ReferenceSigmaPrior Prior) QuadratureCase(string name)
    {
        const int n = 8;
        var y = new double[n, 2];
        var t = new double[n, 1];
        for (int i = 0; i < n; i++)
        {
            double u = -1.0 + 2.0 * i / (n - 1);
            t[i, 0] = u;
            y[i, 0] = 1.5 + Math.Tanh(2.0 * u) + 0.3 * Math.Sin(11.0 * i);
            y[i, 1] = -0.5 * u * u + 0.4 * y[i, 0] + 0.2 * Math.Cos(13.0 * i);
        }
        double[,] xc = ReferenceAlgebra.Centred(t);
        double[,] informative = { { 1.2, -0.5 }, { -0.5, 3.0 } };
        return name switch
        {
            "informative" => (GFormCase(y, xc, GPriorCentre.Mean, 0.5, 1.0), new ReferenceSigmaPrior(5.0, informative)),
            "tempered" => (GFormCase(y, xc, GPriorCentre.Mean, 0.5, 0.5), new ReferenceSigmaPrior(5.0, informative)),
            "decision 2, level centring" => (GFormCase(y, xc, GPriorCentre.Mean, 1.0 / n, 1.0), ReferenceSigmaPrior.ScaledToResponse(y, ColumnMeans(y))),
            "flat intercept" => (GFormCase(y, xc, GPriorCentre.Flat, 0.5, 1.0), new ReferenceSigmaPrior(4.0, informative)),
            _ => throw new ArgumentOutOfRangeException(nameof(name)),
        };
    }

    // BASS's g-form on centred columns X_c (D03): the intercept inside the prior with Gram entry n, centred on ȳ (level centring) or 0
    // (the zero centre), or flat; Λ₀ = τ·blockdiag(n, X_cᵀX_c).
    private static MultiOutputCase GFormCase(double[,] y, double[,] xc, GPriorCentre centre, double tau, double b)
        => Case(y, xc, centre, ScalarForm.Gram, null, tau, b);

    private static double GFormLogMarginal(double[,] y, double[,] xc, GPriorCentre centre, double tau, double b)
    {
        double[] levels = centre == GPriorCentre.Zero ? new double[y.GetLength(1)] : ColumnMeans(y);
        return MatrixTReference.LogMarginal(MatrixTReference.Collapse(GFormCase(y, xc, centre, tau, b)), ReferenceSigmaPrior.ScaledToResponse(y, levels));
    }

    // The case under one of the scalar engine's three prior forms: the g-form, the normalized ridge (D_α = 1, D_β = diag(G_c)/n) or a
    // pencil's prior D, at precision τ.
    private static MultiOutputCase Case(double[,] y, double[,] xc, GPriorCentre centre, ScalarForm form, double[,]? pencilPrior, double tau, double b)
    {
        int n = y.GetLength(0), k = y.GetLength(1), m = xc.GetLength(1);
        double[,] gram = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(xc), xc);
        bool flat = centre == GPriorCentre.Flat;
        int offset = flat ? 0 : 1, p = m + offset;
        var lambda = new double[p, p];
        if (form == ScalarForm.Pencil)
        {
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) lambda[i, j] = tau * pencilPrior![i, j];
        }
        else
        {
            if (!flat) lambda[0, 0] = tau * (form == ScalarForm.Gram ? n : 1.0);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    lambda[i + offset, j + offset] = form == ScalarForm.Gram ? tau * gram[i, j] : i == j ? tau * gram[i, i] / n : 0.0;
        }
        if (flat) return new MultiOutputCase(y, xc, lambda, null, b, true);
        var x = new double[n, p];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = 1.0;
            for (int j = 0; j < m; j++) x[i, j + 1] = xc[i, j];
        }
        double[,]? a0 = null;
        if (centre == GPriorCentre.Mean)
        {
            a0 = new double[p, k];
            double[] means = ColumnMeans(y);
            for (int l = 0; l < k; l++) a0[0, l] = means[l];
        }
        return new MultiOutputCase(y, x, lambda, a0, b);
    }

    private static double ReferenceLogMarginal(double[] y, double[,] xc, GPriorCentre centre, ScalarForm form, double[,]? pencilPrior,
                                               double tau, double b, double c, bool jeffreys)
    {
        int n = y.Length;
        var column = new double[n, 1];
        for (int i = 0; i < n; i++) column[i, 0] = y[i];
        double level = centre == GPriorCentre.Zero ? 0.0 : ColumnMeans(column)[0];
        ReferenceSigmaPrior prior = jeffreys ? new ReferenceSigmaPrior(0.0, null) : ReferenceSigmaPrior.ScaledToResponse(column, [level], c);
        return MatrixTReference.LogMarginal(MatrixTReference.Collapse(Case(column, xc, centre, form, pencilPrior, tau, b)), prior);
    }

    // The engine's scalar marginal from the centred least-squares fit's statistics, by its own forms.
    private static double ScalarLogMarginal(double[] y, double[,] xc, GPriorCentre centre, ScalarForm form, double[,]? pencilPrior,
                                            double tau, double b, double c, bool jeffreys)
    {
        int n = y.Length, m = xc.GetLength(1);
        double mean = 0.0;
        foreach (double v in y) mean += v;
        mean /= n;
        double level = centre == GPriorCentre.Zero ? 0.0 : mean, s2 = 0.0;
        foreach (double v in y) s2 += (v - level) * (v - level);
        s2 /= n;
        InverseGamma? prior = jeffreys ? null : InverseGamma.ShapeScale(c, c * s2);

        double[,] gram = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(xc), xc);
        var cross = new double[m];
        var crossMatrix = new double[m, 1];
        for (int j = 0; j < m; j++)
        {
            for (int i = 0; i < n; i++) cross[j] += xc[i, j] * (y[i] - mean);
            crossMatrix[j, 0] = cross[j];
        }
        double rss = 0.0, fitted = 0.0;
        double[,] beta = m == 0 ? new double[0, 1] : ReferenceAlgebra.Solve(ReferenceAlgebra.Cholesky(gram)!, crossMatrix);
        for (int i = 0; i < n; i++)
        {
            double fit = 0.0;
            for (int j = 0; j < m; j++) fit += xc[i, j] * beta[j, 0];
            rss += (y[i] - mean - fit) * (y[i] - mean - fit);
            fitted += fit * fit;
        }

        GPriorMarginal marginal;
        switch (form)
        {
            case ScalarForm.Gram:
                marginal = CollapsedGPrior.CentredMarginalAt(m, rss, fitted, n, mean, centre, tau, b);
                break;
            case ScalarForm.Ridge:
                marginal = CollapsedGPrior.CentredRidgeMarginalAt(RidgeSpectrum.Of(gram, cross, n), rss, n, mean, centre, tau, b);
                break;
            default:
                bool inside = centre != GPriorCentre.Flat;
                int p = m + (inside ? 1 : 0), offset = inside ? 1 : 0;
                var likelihood = new double[p, p];
                var pencilCross = new double[p];
                if (inside)
                {
                    likelihood[0, 0] = n;
                    pencilCross[0] = n * (mean - level);
                }
                for (int i = 0; i < m; i++)
                {
                    pencilCross[i + offset] = cross[i];
                    for (int j = 0; j < m; j++) likelihood[i + offset, j + offset] = gram[i, j];
                }
                PriorPencil pencil = PriorPencil.Of(likelihood, pencilPrior!, pencilCross, m, inside) ?? throw new InvalidOperationException("The pencil is not positive definite.");
                marginal = CollapsedGPrior.CentredPencilMarginalAt(pencil, rss, n, GPriorWeights.Unit(n), centre, tau, b);
                break;
        }
        return marginal.LogMarginal(prior);
    }

    // ── Helpers ────────────────────────────────────────────────────────────────────────────────────────────────────────────

    private static (double Value, bool Refused) Attempt(Func<double> f)
    {
        try { return (f(), false); }
        catch (InvalidOperationException) { return (double.NaN, true); }
    }

    private static double[] ColumnMeans(double[,] a)
    {
        int n = a.GetLength(0), k = a.GetLength(1);
        var means = new double[k];
        for (int l = 0; l < k; l++)
        {
            for (int i = 0; i < n; i++) means[l] += a[i, l];
            means[l] /= n;
        }
        return means;
    }

    private static double[,] Covariance(double[,] y)
    {
        double[,] c = ReferenceAlgebra.Centred(y);
        double[,] s = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(c), c);
        int n = y.GetLength(0);
        for (int a = 0; a < s.GetLength(0); a++)
            for (int d = 0; d < s.GetLength(1); d++) s[a, d] /= n;
        return s;
    }

    private static double Determinant3(double[,] a)
        => a[0, 0] * (a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1]) - a[0, 1] * (a[1, 0] * a[2, 2] - a[1, 2] * a[2, 0])
           + a[0, 2] * (a[1, 0] * a[2, 1] - a[1, 1] * a[2, 0]);

    private static double[,] RandomMatrix(Xoshiro256PlusPlus rng, int rows, int columns, double scale)
    {
        var a = new double[rows, columns];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++) a[i, j] = scale * Normal.SampleStandard(rng);
        return a;
    }

    // Q·diag(d)·Qᵀ with Q a random rotation and d spread from scale down to scale·10^(−logCondition).
    private static double[,] RandomSpd(Xoshiro256PlusPlus rng, int k, double logCondition, double scale)
    {
        double[,] q = RandomMatrix(rng, k, k, 1.0);
        for (int j = 0; j < k; j++)
        {
            for (int r = 0; r < j; r++)
            {
                double dot = 0.0;
                for (int i = 0; i < k; i++) dot += q[i, j] * q[i, r];
                for (int i = 0; i < k; i++) q[i, j] -= dot * q[i, r];
            }
            double norm = 0.0;
            for (int i = 0; i < k; i++) norm += q[i, j] * q[i, j];
            norm = Math.Sqrt(norm);
            for (int i = 0; i < k; i++) q[i, j] /= norm;
        }
        var s = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int d = 0; d <= a; d++)
            {
                double sum = 0.0;
                for (int j = 0; j < k; j++) sum += q[a, j] * q[d, j] * scale * Math.Pow(10.0, -logCondition * (k == 1 ? 0.0 : (double)j / (k - 1)));
                s[a, d] = s[d, a] = sum;
            }
        return s;
    }

    private static void Close(double expected, double actual, double rel)
        => Assert.True(Math.Abs(actual - expected) <= rel * Math.Max(1.0, Math.Abs(expected)), $"expected {expected:R}, got {actual:R}");
}
