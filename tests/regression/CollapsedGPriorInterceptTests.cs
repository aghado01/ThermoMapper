using System;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The block's centred forms, the intercept apart from centred columns (BASS's D03): level centring, the zero centre and
/// Liang et al.'s flat intercept. Checked against references computed another way:
/// <list type="bullet">
/// <item>dense algebra on the explicit design [1, Z_c], and on the uncentred [1, H] for the two g-prior conventions,
/// whose marginal depends only on the column space;</item>
/// <item>the flat form's marginal and σ²'s posterior mean against quadrature over the intercept and σ², at several b;</item>
/// <item>Liang's and level centring's Bayes factors, τ integrated under Zellner–Siow, against F-O9's oracle on F-O8's grid
/// and against the g-forms on an explicit hinge design, with shift invariance;</item>
/// <item>draws, leave-out densities and readout variances against their dense counterparts.</item>
/// </list>
/// </summary>
public sealed class CollapsedGPriorInterceptTests
{
    private static readonly GPriorCentre[] Conventions = [GPriorCentre.Mean, GPriorCentre.Zero, GPriorCentre.Flat];

    public static TheoryData<double, double> Settings => new()
    {
        { 1.0 / 40.0, 1.0 },
        { 0.01, 1.0 },
        { 0.3, 0.25 },
        { 2.0, 0.7 },
    };

    [Theory]
    [MemberData(nameof(Settings))]
    public void CentredForms_MatchDenseAlgebra_AtAnyTauAndTemperature(double tau, double b)
    {
        Design d = Design.Hinges(offset: 2.0);
        InverseGamma proper = InverseGamma.ShapeScale(2.0, 0.05);
        foreach (GPriorCentre convention in Conventions)
        {
            GPriorConditional c = d.Block(convention).At(tau, b);
            Dense e = Dense.Of(d.Explicit(), d.Y, convention, d.Mean, tau, b);
            Close(e.S, c.S, 1e-11);
            Close(e.LogDeterminantTerm, c.LogDeterminantTerm, 1e-12);
            Close(e.Nu, c.Sigma2Posterior(null).Shape * 2.0, 1e-14);
            for (int j = 0; j <= d.M; j++) Close(e.Mean[j], c.Mean[j], 1e-11);
            Close(e.LogMarginal(proper), c.LogMarginal(proper), 1e-11);
            Close(e.LogMarginal(null), c.LogMarginal(null), 1e-11);

            if (convention == GPriorCentre.Flat) continue;
            // The g-prior conventions read the column space alone, so R's uncentred Gram over [1, H] gives the same.
            Dense u = Dense.Of(d.Uncentred(), d.Y, convention, d.Mean, tau, b);
            Close(u.LogMarginal(proper), c.LogMarginal(proper), 1e-10);
            Close(u.S, c.S, 1e-10);
        }
    }

    [Fact]
    public void LevelCentring_IsTheExplicitInterceptColumn_ThroughFromLeastSquares()
    {
        // FromLeastSquares on the factor of [1, Z_c]'s Gram, the statistics of ỹ and μ₀ = ȳ·e₀: the same block.
        Design d = Design.Hinges(offset: -7.0);
        double[,] x = d.Explicit();
        int p = d.M + 1;
        double[] yt = d.Y.Select(v => v - d.Mean).ToArray();
        var factor = new DenseCholesky(p);
        factor.Decompose(Gram(x));
        double[] beta = factor.Solve(CrossProduct(x, yt));
        (double rss, double q) = Fit(x, yt, beta);
        var priorMean = new double[p];
        priorMean[0] = d.Mean;
        CollapsedGPrior explicitBlock = CollapsedGPrior.FromLeastSquares(factor, beta, rss, q, d.N, priorMean);
        CollapsedGPrior centred = d.Block(GPriorCentre.Mean);
        Assert.Equal(GPriorCentre.Mean, centred.Intercept);
        Assert.Null(explicitBlock.Intercept);
        Assert.Equal(explicitBlock.Columns, centred.Columns);

        foreach ((double tau, double b) in new[] { (0.025, 1.0), (0.3, 0.25), (5.0, 0.0) })
        {
            GPriorConditional a = explicitBlock.At(tau, b), c = centred.At(tau, b);
            Close(a.S, c.S, 1e-12);
            Close(a.LogDeterminantTerm, c.LogDeterminantTerm, 1e-14);
            Close(a.LogMarginal(InverseGamma.ShapeScale(1e-3, 1e-3)), c.LogMarginal(InverseGamma.ShapeScale(1e-3, 1e-3)), 1e-12);
            for (int j = 0; j < p; j++) Close(a.Mean[j], c.Mean[j], 1e-12);
        }
    }

    public static TheoryData<double, double> FlatQuadratureSettings => new()
    {
        { 1.0, 0.0 },   // g₁ = 0: the improper 1/σ²
        { 1.0, 2.0 },
        { 0.5, 0.0 },
        { 0.5, 2.0 },
        { 0.2, 0.0 },
        { 0.2, 2.0 },
    };

    [Theory]
    [MemberData(nameof(FlatQuadratureSettings))]
    public void FlatMarginal_MatchesQuadratureOverTheInterceptAndSigma2(double b, double g1)
    {
        // ∫∫ p_b(y | α, σ²) π(σ²) dα dσ² with α under its flat prior (density 1), β's integral from dense algebra at each
        // α, and α and u = ln σ² by quadrature: the (2πσ²/(Nb))^½ and the shape (Nb − 1)/2 are what this checks.
        Design d = Design.Hinges(offset: 2.0);
        const double tau = 0.1;
        InverseGamma? prior = g1 > 0.0 ? InverseGamma.ShapeScale(g1, 0.05) : null;
        double[,] zc = d.Centred();
        double[,] gInverse = Inverse(Gram(zc));
        int n = d.N, m = d.M;

        double LogPrior(double s2) => prior is InverseGamma ig ? ig.LogPdf(s2) : -Math.Log(s2);
        double Energy(double alpha)
        {
            var r = d.Y.Select(v => v - alpha).ToArray();
            double[] h = CrossProduct(zc, r);
            double rr = r.Sum(v => v * v), hgh = 0.0;
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) hgh += h[i] * gInverse[i, j] * h[j];
            return b * rr - b * b / (tau + b) * hgh;
        }

        const int nu = 3000, na = 160;
        const double lo = -14.0, hi = 9.0;
        double hu = (hi - lo) / nu;
        var logs = new double[nu + 1];
        var weightedSigma2 = new double[nu + 1];
        for (int iu = 0; iu <= nu; iu++)
        {
            double u = lo + iu * hu, s2 = Math.Exp(u), half = 12.0 * Math.Sqrt(s2 / (n * b)), ha = 2.0 * half / na;
            var inner = new double[na + 1];
            for (int ia = 0; ia <= na; ia++)
            {
                double alpha = d.Mean - half + ia * ha;
                inner[ia] = -0.5 * n * b * Math.Log(2.0 * Math.PI * s2) + 0.5 * m * Math.Log(tau / (tau + b)) - Energy(alpha) / (2.0 * s2);
            }
            logs[iu] = LogTrapezoid(inner, ha) + LogPrior(s2) + u;   // + u: dσ² = σ² du
            weightedSigma2[iu] = logs[iu] + u;
        }
        double quadrature = LogTrapezoid(logs, hu);

        GPriorConditional c = d.Block(GPriorCentre.Flat).At(tau, b);
        Assert.True(Math.Abs(c.LogMarginal(prior) - quadrature) < 1e-8, $"log marginal {c.LogMarginal(prior):R} vs quadrature {quadrature:R}");
        double meanSigma2 = Math.Exp(LogTrapezoid(weightedSigma2, hu) - quadrature);
        Close(meanSigma2, c.Sigma2Posterior(prior).Mean, 1e-8);
    }

    // F-O9's grid, F-O8's (N, k, R²): ln Bayes factors against the intercept-only structure, from bass_intercept_oracle.R
    // (grid mode), Liang's and level centring's.
    public static TheoryData<int, int, double, double, double> OracleGrid => new()
    {
        { 40, 1, 0.05, -1.138514956796143, -0.70487073407550671 },
        { 40, 1, 0.2, 2.0594135302426606, 2.5053225027387063 },
        { 40, 1, 0.5, 10.788879721289103, 11.257362710480709 },
        { 40, 2, 0.25, 1.4802077357669088, 2.1559032456022651 },
        { 100, 1, 0.05, -0.047637447780182418, 0.39641170610310095 },
        { 100, 1, 0.2, 8.2950257050944494, 8.7446410222152888 },
        { 100, 1, 0.5, 31.104015384047102, 31.563157513048946 },
        { 100, 2, 0.25, 9.2006438442700258, 9.8870927135596638 },
        { 1000, 1, 0.05, 21.889333880980445, 22.3401398897775 },
        { 1000, 1, 0.2, 107.5575086193892, 108.00891836669309 },
        { 1000, 1, 0.5, 341.85569813266397, 342.30808756914689 },
        { 1000, 2, 0.25, 136.35767629988567, 137.05017948532907 },
    };

    [Theory]
    [MemberData(nameof(OracleGrid))]
    public void BayesFactors_UnderZellnerSiow_MatchFO9sOracle(int n, int k, double r2, double liang, double level)
    {
        // A design with exactly that R²: k cosine columns, orthogonal to the constant and to each other, and a response
        // at level 3 whose centred part is √R² along them and √(1 − R²) along a further cosine.
        var columns = new double[k][];
        for (int j = 0; j < k; j++) columns[j] = Cosine(n, j + 1);
        double[] signal = new double[n], noise = Unit(Cosine(n, k + 1));
        foreach (double[] col in columns.Select(Unit)) for (int i = 0; i < n; i++) signal[i] += col[i] / Math.Sqrt(k);
        var y = new double[n];
        for (int i = 0; i < n; i++) y[i] = 3.0 + Math.Sqrt(r2) * signal[i] + Math.Sqrt(1.0 - r2) * noise[i];
        var structure = new Design(y, columns);
        var intercept = new Design(y, []);

        Assert.True(Math.Abs(BayesFactor(structure, intercept, GPriorCentre.Flat) - liang) < 1e-8, "Liang");
        Assert.True(Math.Abs(BayesFactor(structure, intercept, GPriorCentre.Mean) - level) < 1e-8, "level centring");
    }

    [Fact]
    public void ShiftingY_MovesNeitherLevelCentringNorLiang_ButMovesTheZeroCentre_AndTheGFormsHold()
    {
        // F-O9's design mode on a hinge design of our own: the intercept alone, one hinge, two hinges.
        double[] Factors(double offset, GPriorCentre convention)
        {
            Design two = Design.Hinges(offset), one = new(two.Y, [two.Columns[0]]), none = new(two.Y, []);
            return [BayesFactor(one, none, convention), BayesFactor(two, none, convention)];
        }

        double[] liang0 = Factors(0.0, GPriorCentre.Flat), level0 = Factors(0.0, GPriorCentre.Mean), zero0 = Factors(0.0, GPriorCentre.Zero);
        foreach (double offset in new[] { 10.0, -250.0 })
            for (int s = 0; s < 2; s++)
            {
                Assert.True(Math.Abs(Factors(offset, GPriorCentre.Flat)[s] - liang0[s]) < 1e-8, $"Liang moved at {offset}");
                Assert.True(Math.Abs(Factors(offset, GPriorCentre.Mean)[s] - level0[s]) < 1e-8, $"level centring moved at {offset}");
            }
        Assert.True(Math.Abs(Factors(-250.0, GPriorCentre.Zero)[1] - zero0[1]) > 1.0, "the zero centre should move with y's offset");

        // The g-forms with each structure's centred R²: Liang's, and level centring's with −N/2 and the null's constant.
        Design hinges = Design.Hinges(0.0);
        for (int s = 0; s < 2; s++)
        {
            var design = new Design(hinges.Y, hinges.Columns.Take(s + 1).ToArray());
            double r2 = design.CentredFit().Q / design.CentredResponseSquares();
            int n = design.N, k = s + 1;
            double liang = GForm(n, u => (n - 1 - k) / 2.0 * Log1pExp(u) - (n - 1) / 2.0 * Math.Log(1.0 + Math.Exp(u) * (1.0 - r2)));
            double level = GForm(n, u => (n - 1 - k) / 2.0 * Log1pExp(u) - n / 2.0 * Math.Log(1.0 + Math.Exp(u) * (1.0 - r2)))
                         - GForm(n, u => -0.5 * Log1pExp(u));
            Assert.True(Math.Abs(liang0[s] - liang) < 1e-8, $"Liang's g-form, {k} hinges: {liang0[s]:R} vs {liang:R}");
            Assert.True(Math.Abs(level0[s] - level) < 1e-8, $"level centring's g-form, {k} hinges: {level0[s]:R} vs {level:R}");
        }
    }

    [Fact]
    public void FlatDraws_HaveTheirConditionalMoments_TheInterceptApartFromTheColumns()
    {
        // σ² | y ~ IG(g₁ + (Nb − 1)/2, g₂ + S/2); α | σ² ~ N(ȳ, σ²/(Nb)); β_c | σ² ~ N(b/(τ+b)·β̂_c, σ²G_c⁻¹/(τ+b)).
        Design d = Design.Hinges(offset: 2.0);
        const double tau = 0.05, b = 0.5;
        var prior = InverseGamma.ShapeScale(3.0, 0.4);
        GPriorConditional c = d.Block(GPriorCentre.Flat).At(tau, b);
        Dense e = Dense.Of(d.Explicit(), d.Y, GPriorCentre.Flat, d.Mean, tau, b);
        InverseGamma post = InverseGamma.ShapeScale(3.0 + 0.5 * (d.N * b - 1.0), 0.4 + 0.5 * e.S);
        Close(post.Mean, c.Sigma2Posterior(prior).Mean, 1e-12);
        Assert.Equal(1.0 / b, c.InterceptCovarianceScale);
        Close(d.Mean, c.Mean[0], 1e-15);

        var rng = new Xoshiro256PlusPlus(seed: 31);
        const int draws = 40_000;
        int p = d.M + 1;
        var theta = new double[p];
        var mean = new double[p];
        var second = new double[p, p];
        double s2Sum = 0.0;
        for (int s = 0; s < draws; s++)
        {
            s2Sum += c.Sample(prior, rng, theta);
            for (int i = 0; i < p; i++) { mean[i] += theta[i] / draws; for (int j = 0; j < p; j++) second[i, j] += theta[i] * theta[j] / draws; }
        }
        Assert.True(Math.Abs(s2Sum / draws - post.Mean) < 5.0 * Math.Sqrt(post.Variance / draws), $"E[σ²] {s2Sum / draws:R} vs {post.Mean:R}");
        for (int i = 0; i < p; i++)
        {
            double sd = Math.Sqrt(post.Mean * e.PInverse[i, i]);
            Assert.True(Math.Abs(mean[i] - e.Mean[i]) < 5.0 * sd / Math.Sqrt(draws), $"mean[{i}] {mean[i]:R} vs {e.Mean[i]:R}");
            for (int j = 0; j < p; j++)
            {
                double cov = second[i, j] - mean[i] * mean[j], expected = post.Mean * e.PInverse[i, j];
                double tol = 0.05 * post.Mean * Math.Sqrt(e.PInverse[i, i] * e.PInverse[j, j]);
                Assert.True(Math.Abs(cov - expected) < tol, $"cov[{i},{j}] {cov:R} vs {expected:R}");
            }
        }
    }

    [Fact]
    public void LeaveOutDensities_AreRatiosOfDenseMarginals_AndReadoutVariancesAreDense()
    {
        // ln p(y_U | y_rest) = ln p(y) − ln p(y_rest), each a dense marginal under the full data's prior (its Gram, τ
        // and centre), for one and two held-out observations; the readout variance E[σ² | y]·xᵀP⁻¹x.
        Design d = Design.Hinges(offset: 2.0);
        double[,] x = d.Explicit();
        const double tau = 0.08;
        InverseGamma proper = InverseGamma.ShapeScale(2.0, 0.05);
        foreach (GPriorCentre convention in Conventions)
        {
            GPriorConditional c = d.Block(convention).At(tau);
            Dense full = Dense.Of(x, d.Y, convention, d.Mean, tau, 1.0);
            foreach (int[] unit in new[] { new[] { 3 }, new[] { 5, 17 } })
            {
                int[] rest = Enumerable.Range(0, d.N).Except(unit).ToArray();
                double[,] xRest = Rows(x, rest);
                double[] yRest = rest.Select(i => d.Y[i]).ToArray();
                Dense partial = Dense.Of(xRest, yRest, convention, d.Mean, tau, 1.0, full.Lambda);
                foreach (InverseGamma? prior in new[] { proper, (InverseGamma?)null })
                    Close(full.LogMarginal(prior) - partial.LogMarginal(prior), c.LeaveOutLogDensity(prior, x, d.Y, null, unit), 1e-10);
            }

            double[] v = c.ConditionalVariance(proper, x);
            double eSigma2 = c.Sigma2Posterior(proper).Mean;
            for (int t = 0; t < d.N; t++)
            {
                double q = 0.0;
                for (int i = 0; i <= d.M; i++) for (int j = 0; j <= d.M; j++) q += x[t, i] * full.PInverse[i, j] * x[t, j];
                Close(eSigma2 * q, v[t], 1e-10);
            }
        }
    }

    [Fact]
    public void FlatMarginalAt_IsTheConditionalsBitForBit_WithoutAllocating()
    {
        Design d = Design.Hinges(offset: 2.0);
        CollapsedGPrior block = d.Block(GPriorCentre.Flat);
        InverseGamma proper = InverseGamma.ShapeScale(1e-3, 1e-3);
        foreach (double tau in new[] { 1e-6, 0.025, 1.0, 1e3 })
            foreach (double b in new[] { 0.1, 0.5, 1.0 })
                foreach (InverseGamma? prior in new[] { proper, (InverseGamma?)null })
                    Assert.Equal(BitConverter.DoubleToInt64Bits(block.At(tau, b).LogMarginal(prior)),
                                 BitConverter.DoubleToInt64Bits(block.MarginalAt(tau, b).LogMarginal(prior)));

        double sink = 0.0;
        double[] taus = [1e-4, 0.1, 10.0];
        for (int k = 0; k < 2; k++) sink += block.MarginalAt(0.3, 0.7).LogMarginal(proper);
        long before = GC.GetAllocatedBytesForCurrentThread();
        foreach (double tau in taus) sink += block.MarginalAt(tau, 0.7).LogMarginal(proper);
        Assert.Equal(before, GC.GetAllocatedBytesForCurrentThread());
        Assert.True(double.IsFinite(sink));
    }

    [Fact]
    public void FlatForm_RefusesBZero_AndAShapeItCannotIntegrate_AndTheMedianIsRefused()
    {
        Design d = Design.Hinges(offset: 2.0);
        CollapsedGPrior flat = d.Block(GPriorCentre.Flat);
        Assert.Equal(GPriorCentre.Flat, flat.Intercept);
        Assert.Equal(d.M + 1, flat.Columns);
        Assert.Throws<ArgumentOutOfRangeException>(() => flat.MarginalAt(0.1, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => flat.At(0.1, 0.0));

        // Nb = ½ leaves the improper prior a negative shape; a proper prior's g₁ + ν/2 stays positive there.
        GPriorMarginal edge = flat.MarginalAt(0.1, 0.5 / d.N);
        Assert.Throws<InvalidOperationException>(() => edge.LogMarginal(null));
        Assert.True(double.IsFinite(edge.LogMarginal(InverseGamma.ShapeScale(0.5, 0.05))));

        // Level centring keeps b = 0: a proper prior's marginal is exactly 0 there.
        Assert.Equal(0.0, d.Block(GPriorCentre.Mean).MarginalAt(0.1, 0.0).LogMarginal(InverseGamma.ShapeScale(2.0, 0.5)));

        (ISpdFactor factor, double[] beta, double rss, double q) = d.CentredFit();
        Assert.Throws<ArgumentException>(() => CollapsedGPrior.FromCentredLeastSquares(factor, beta, rss, q, d.N, d.Mean, GPriorCentre.Median));
        Assert.Throws<ArgumentException>(() => CollapsedGPrior.FromCentredLeastSquares(factor, [1.0], rss, q, d.N, d.Mean, GPriorCentre.Flat));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollapsedGPrior.FromCentredLeastSquares(factor, beta, rss, q, d.N, double.NaN, GPriorCentre.Flat));
    }

    // ln BF of a structure against the intercept-only one, τ ~ Gamma(½, rate N/2) (Zellner–Siow) integrated on v = ln τ,
    // under the improper σ² prior, whose constant every structure on the same response shares.
    private static double BayesFactor(Design structure, Design intercept, GPriorCentre convention)
    {
        CollapsedGPrior a = structure.Block(convention), z = intercept.Block(convention);
        int n = structure.N;
        double LogPrior(double v) => 0.5 * (v + Math.Log(n / 2.0)) - n / 2.0 * Math.Exp(v) - SpecialFunctions.LogGamma(0.5);
        double Integral(CollapsedGPrior block)
        {
            const int points = 60_000;
            const double lo = -50.0, hi = 10.0;
            double h = (hi - lo) / points;
            var logs = new double[points + 1];
            for (int i = 0; i <= points; i++)
            {
                double v = lo + i * h;
                logs[i] = LogPrior(v) + block.MarginalAt(Math.Exp(v)).LogMarginal(null);
            }
            return LogTrapezoid(logs, h);
        }
        return Integral(a) - Integral(z);
    }

    // ln ∫ p(g) exp(f(u)) dg over u = ln g, g ~ IG(½, N/2), as the oracle writes its g-forms.
    private static double GForm(int n, Func<double, double> f)
    {
        const int points = 120_000;
        const double lo = -60.0, hi = 60.0;
        double h = (hi - lo) / points;
        var logs = new double[points + 1];
        for (int i = 0; i <= points; i++)
        {
            double u = lo + i * h;
            logs[i] = 0.5 * Math.Log(n / 2.0) - SpecialFunctions.LogGamma(0.5) - 0.5 * u - n / 2.0 * Math.Exp(-u) + f(u);
        }
        return LogTrapezoid(logs, h);
    }

    private static double Log1pExp(double u) => u > 30.0 ? u + Math.Exp(-u) : Math.Log(1.0 + Math.Exp(u));

    private static double LogTrapezoid(double[] logs, double h)
    {
        double max = logs.Max(), sum = 0.0;
        for (int i = 0; i < logs.Length; i++) sum += (i == 0 || i == logs.Length - 1 ? 0.5 : 1.0) * Math.Exp(logs[i] - max);
        return max + Math.Log(sum * h);
    }

    private static double[] Cosine(int n, int frequency)
        => Enumerable.Range(0, n).Select(i => Math.Cos(Math.PI * frequency * (i + 0.5) / n)).ToArray();

    private static double[] Unit(double[] v)
    {
        double norm = Math.Sqrt(v.Sum(t => t * t));
        return v.Select(t => t / norm).ToArray();
    }

    private static void Close(double expected, double actual, double rel)
        => Assert.True(Math.Abs(actual - expected) <= rel * Math.Max(1.0, Math.Abs(expected)), $"expected {expected:R}, got {actual:R}");

    // A response and its uncentred columns, from which each convention's centred block is built as BASS's design would.
    private sealed class Design(double[] y, double[][] columns)
    {
        public double[] Y { get; } = y;
        public double[][] Columns { get; } = columns;
        public int N => Y.Length;
        public int M => Columns.Length;
        public double Mean => Y.Average();

        // 40 points, a sine with noise at the given level; hinges (x − 0.3)₊ and (0.6 − x)₊.
        public static Design Hinges(double offset)
        {
            var rng = new Xoshiro256PlusPlus(seed: 9);
            double[] x = Enumerable.Range(0, 40).Select(_ => rng.NextDouble()).OrderBy(v => v).ToArray();
            double[] y = x.Select(v => offset + 2.0 + 0.8 * Math.Sin(2.0 * Math.PI * v) + 0.3 * Normal.SampleStandard(rng)).ToArray();
            return new Design(y, [x.Select(v => Math.Max(v - 0.3, 0.0)).ToArray(), x.Select(v => Math.Max(0.6 - v, 0.0)).ToArray()]);
        }

        public double[,] Centred()
        {
            var z = new double[N, M];
            for (int j = 0; j < M; j++)
            {
                double mean = Columns[j].Average();
                for (int i = 0; i < N; i++) z[i, j] = Columns[j][i] - mean;
            }
            return z;
        }

        public double[,] Explicit() => WithConstant(Centred());

        public double[,] Uncentred()
        {
            var h = new double[N, M];
            for (int i = 0; i < N; i++) for (int j = 0; j < M; j++) h[i, j] = Columns[j][i];
            return WithConstant(h);
        }

        public double CentredResponseSquares() => Y.Sum(v => (v - Mean) * (v - Mean));

        public (ISpdFactor Factor, double[] Beta, double Rss, double Q) CentredFit()
        {
            double[,] z = Centred();
            double[] yt = Y.Select(v => v - Mean).ToArray();
            var factor = new DenseCholesky(M);
            factor.Decompose(Gram(z));
            double[] beta = factor.Solve(CrossProduct(z, yt));
            (double rss, double q) = Fit(z, yt, beta);
            return (factor, beta, rss, q);
        }

        public CollapsedGPrior Block(GPriorCentre convention)
        {
            (ISpdFactor factor, double[] beta, double rss, double q) = CentredFit();
            return CollapsedGPrior.FromCentredLeastSquares(factor, beta, rss, q, N, Mean, convention);
        }

        private double[,] WithConstant(double[,] z)
        {
            var x = new double[N, M + 1];
            for (int i = 0; i < N; i++)
            {
                x[i, 0] = 1.0;
                for (int j = 0; j < M; j++) x[i, j + 1] = z[i, j];
            }
            return x;
        }
    }

    // The dense reference on an explicit design X = [1, …]: p_b(y | θ, σ²) = (2πσ²)^(−nb/2) e^(−b‖y − Xθ‖²/2σ²), the
    // g-prior N(m₀, σ²Λ⁻¹) on its coordinates (all of them, or all but a flat intercept's, whose density is 1), with
    // Λ = τ·XᵀX there unless a full-data Λ is held. Then P = bXᵀX + Λ, S = b·yᵀy + m₀ᵀΛm₀ − rᵀP⁻¹r with r = bXᵀy + Λm₀,
    // the determinant term ½ln|Λ_g| − ½ln|P|, and ν = nb less one per flat coordinate.
    private sealed record Dense(double Nu, double LogDeterminantTerm, double S, double[] Mean, double[,] PInverse, double[,] Lambda)
    {
        public static Dense Of(double[,] x, double[] y, GPriorCentre convention, double level, double tau, double b, double[,]? lambda = null)
        {
            int n = y.Length, d = x.GetLength(1), first = convention == GPriorCentre.Flat ? 1 : 0;
            double[,] xtx = Gram(x);
            double[] xty = CrossProduct(x, y);
            if (lambda is null)
            {
                lambda = new double[d, d];
                for (int i = first; i < d; i++) for (int j = first; j < d; j++) lambda[i, j] = tau * xtx[i, j];
            }
            var m0 = new double[d];
            if (convention == GPriorCentre.Mean) m0[0] = level;

            var p = new double[d, d];
            var r = new double[d];
            double m0Lm0 = 0.0;
            for (int i = 0; i < d; i++)
            {
                for (int j = 0; j < d; j++)
                {
                    p[i, j] = b * xtx[i, j] + lambda[i, j];
                    r[i] += lambda[i, j] * m0[j];
                    m0Lm0 += m0[i] * lambda[i, j] * m0[j];
                }
                r[i] += b * xty[i];
            }
            double[,] pInv = Inverse(p);
            var mean = new double[d];
            double rPr = 0.0;
            for (int i = 0; i < d; i++) { for (int j = 0; j < d; j++) mean[i] += pInv[i, j] * r[j]; rPr += r[i] * mean[i]; }
            var lambdaG = new double[d - first, d - first];
            for (int i = first; i < d; i++) for (int j = first; j < d; j++) lambdaG[i - first, j - first] = lambda[i, j];
            double s = b * y.Sum(v => v * v) + m0Lm0 - rPr;
            double logDet = (d > first ? 0.5 * LogDet(lambdaG) : 0.0) - 0.5 * LogDet(p);
            return new Dense(n * b - first, logDet, s, mean, pInv, lambda);
        }

        public double LogMarginal(InverseGamma? prior)
        {
            double half = 0.5 * Nu, common = -half * Math.Log(2.0 * Math.PI) + LogDeterminantTerm;
            return prior is InverseGamma ig
                ? common + ig.Shape * Math.Log(ig.Scale) - SpecialFunctions.LogGamma(ig.Shape)
                  + SpecialFunctions.LogGamma(ig.Shape + half) - (ig.Shape + half) * Math.Log(ig.Scale + 0.5 * S)
                : common + SpecialFunctions.LogGamma(half) - half * Math.Log(0.5 * S);
        }
    }

    private static double[,] Gram(double[,] x)
    {
        int n = x.GetLength(0), d = x.GetLength(1);
        var g = new double[d, d];
        for (int i = 0; i < n; i++) for (int a = 0; a < d; a++) for (int c = 0; c < d; c++) g[a, c] += x[i, a] * x[i, c];
        return g;
    }

    private static double[] CrossProduct(double[,] x, double[] y)
    {
        int n = x.GetLength(0), d = x.GetLength(1);
        var h = new double[d];
        for (int i = 0; i < n; i++) for (int a = 0; a < d; a++) h[a] += x[i, a] * y[i];
        return h;
    }

    private static (double Rss, double Q) Fit(double[,] x, double[] y, double[] beta)
    {
        double rss = 0.0, q = 0.0;
        for (int i = 0; i < y.Length; i++)
        {
            double f = 0.0;
            for (int j = 0; j < beta.Length; j++) f += x[i, j] * beta[j];
            rss += (y[i] - f) * (y[i] - f);
            q += f * f;
        }
        return (rss, q);
    }

    private static double[,] Rows(double[,] x, int[] rows)
    {
        var r = new double[rows.Length, x.GetLength(1)];
        for (int i = 0; i < rows.Length; i++) for (int j = 0; j < x.GetLength(1); j++) r[i, j] = x[rows[i], j];
        return r;
    }

    // Gauss–Jordan with partial pivoting, a factorization of its own beside the block's Cholesky.
    private static double[,] Inverse(double[,] a0)
    {
        int n = a0.GetLength(0);
        var a = (double[,])a0.Clone();
        var inv = new double[n, n];
        for (int i = 0; i < n; i++) inv[i, i] = 1.0;
        for (int c = 0; c < n; c++)
        {
            int piv = c;
            for (int r = c + 1; r < n; r++) if (Math.Abs(a[r, c]) > Math.Abs(a[piv, c])) piv = r;
            for (int k = 0; k < n; k++) { (a[c, k], a[piv, k]) = (a[piv, k], a[c, k]); (inv[c, k], inv[piv, k]) = (inv[piv, k], inv[c, k]); }
            double d = a[c, c];
            for (int k = 0; k < n; k++) { a[c, k] /= d; inv[c, k] /= d; }
            for (int r = 0; r < n; r++)
            {
                if (r == c) continue;
                double f = a[r, c];
                for (int k = 0; k < n; k++) { a[r, k] -= f * a[c, k]; inv[r, k] -= f * inv[c, k]; }
            }
        }
        return inv;
    }

    private static double LogDet(double[,] a0)
    {
        int n = a0.GetLength(0);
        var a = (double[,])a0.Clone();
        double logDet = 0.0;
        for (int c = 0; c < n; c++)
        {
            int piv = c;
            for (int r = c + 1; r < n; r++) if (Math.Abs(a[r, c]) > Math.Abs(a[piv, c])) piv = r;
            for (int k = 0; k < n; k++) (a[c, k], a[piv, k]) = (a[piv, k], a[c, k]);
            logDet += Math.Log(Math.Abs(a[c, c]));
            for (int r = c + 1; r < n; r++)
            {
                double f = a[r, c] / a[c, c];
                for (int k = c; k < n; k++) a[r, k] -= f * a[c, k];
            }
        }
        return logDet;
    }
}
