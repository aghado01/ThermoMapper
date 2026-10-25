using System;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BAPS and <see cref="HierarchicalPSpline"/> on any <see cref="IBandPenalty"/>, checked with the tensor penalty of a
/// surface, whose bandwidth (order_x·ν_y) and nullity (order_x·order_y) differ from a difference penalty's single order.
/// Its roughness is βᵀPβ for P built from an explicit Kronecker product; BAPS's first coefficient draw follows
/// N(A⁻¹Zᵀy, A⁻¹) with A = ZᵀZ + P built independently; the smoothing variances' conditionals use the nullity; and the
/// locally adaptive scheme, which needs local terms of full row rank, refuses the tensor penalty.
/// </summary>
public sealed class PenaltySeamTests
{
    private const int Nu1 = 4;   // a cubic basis with no interior knots in each direction: 16 coefficients

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static (double[,] Z, double[] Xs, double[] Ys) Surface(int seed, int n)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var xs = new double[n];
        var ys = new double[n];
        for (int i = 0; i < n; i++) { xs[i] = rng.NextDouble(); ys[i] = rng.NextDouble(); }
        var none = new KnotConfig(Array.Empty<double>());
        TensorDesign t = TensorProductDesign.Build(new SplineBasis(3), none, xs, new SplineBasis(3), none, ys);
        return (t.Design, xs, ys);
    }

    private static double SurfaceAt(double x, double y) => Math.Sin(2.0 * x) + 0.5 * y * y;

    // P = D_xᵀD_x ⊗ I + I ⊗ D_yᵀD_y with second differences, coefficients flattened as j·ν_y + k.
    private static double[,] KroneckerPenalty()
    {
        var d = new double[Nu1 - 2, Nu1];
        for (int i = 0; i < Nu1 - 2; i++) { d[i, i] = 1.0; d[i, i + 1] = -2.0; d[i, i + 2] = 1.0; }
        var p1 = new double[Nu1, Nu1];
        for (int a = 0; a < Nu1; a++)
            for (int b = 0; b < Nu1; b++)
                for (int i = 0; i < Nu1 - 2; i++) p1[a, b] += d[i, a] * d[i, b];
        int nu = Nu1 * Nu1;
        var p = new double[nu, nu];
        for (int j = 0; j < Nu1; j++)
            for (int k = 0; k < Nu1; k++)
                for (int j2 = 0; j2 < Nu1; j2++)
                    for (int k2 = 0; k2 < Nu1; k2++)
                        p[j * Nu1 + k, j2 * Nu1 + k2] = p1[j, j2] * (k == k2 ? 1.0 : 0.0) + (j == j2 ? 1.0 : 0.0) * p1[k, k2];
        return p;
    }

    private static double[,] Cholesky(double[,] a)
    {
        int n = a.GetLength(0);
        var l = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = a[i, j];
                for (int k = 0; k < j; k++) s -= l[i, k] * l[j, k];
                l[i, j] = i == j ? Math.Sqrt(s) : s / l[j, j];
            }
        return l;
    }

    private static double[] Solve(double[,] l, double[] b)
    {
        int n = b.Length;
        var z = new double[n];
        for (int i = 0; i < n; i++) { double s = b[i]; for (int k = 0; k < i; k++) s -= l[i, k] * z[k]; z[i] = s / l[i, i]; }
        var x = new double[n];
        for (int i = n - 1; i >= 0; i--) { double s = z[i]; for (int k = i + 1; k < n; k++) s -= l[k, i] * x[k]; x[i] = s / l[i, i]; }
        return x;
    }

    [Fact]
    public void TensorRoughness_IsTheKroneckerQuadratic()
    {
        var penalty = new TensorPenalty(Nu1, Nu1);
        double[,] p = KroneckerPenalty();
        var rng = new Xoshiro256PlusPlus(seed: 2);
        for (int trial = 0; trial < 5; trial++)
        {
            double[] beta = Enumerable.Range(0, Nu1 * Nu1).Select(_ => Gaussian(rng)).ToArray();
            double q = 0.0;
            for (int i = 0; i < beta.Length; i++)
                for (int j = 0; j < beta.Length; j++) q += beta[i] * p[i, j] * beta[j];
            Assert.Equal(q, penalty.Roughness(beta), 10);
        }
        Assert.Equal(4, penalty.Nullity);
        Assert.Equal(2 * Nu1, penalty.Bandwidth);
    }

    [Fact]
    public void Baps_OnATensorPenalty_DrawsTheFirstSweepFromItsConditionals()
    {
        (double[,] z, double[] xs, double[] ys) = Surface(seed: 5, n: 120);
        var rng = new Xoshiro256PlusPlus(seed: 6);
        double[] y = xs.Select((x, i) => SurfaceAt(x, ys[i]) + 0.2 * Gaussian(rng)).ToArray();
        int n = y.Length, nu = Nu1 * Nu1;

        // From the fixed start σ² = τ² = 1: β ~ N(A⁻¹Zᵀy, A⁻¹) with A = ZᵀZ + P.
        double[,] p = KroneckerPenalty();
        var a = new double[nu, nu];
        var zty = new double[nu];
        for (int i = 0; i < nu; i++)
        {
            for (int t = 0; t < n; t++) zty[i] += z[t, i] * y[t];
            for (int j = 0; j < nu; j++)
            {
                double s = p[i, j];
                for (int t = 0; t < n; t++) s += z[t, i] * z[t, j];
                a[i, j] = s;
            }
        }
        double[,] l = Cholesky(a);
        double[] mean = Solve(l, zty);

        var baps = new Baps(z, y, new TensorPenalty(Nu1, Nu1));
        const int draws = 1500;
        var white = new double[draws][];
        var scaled = new double[draws];
        for (int s = 0; s < draws; s++)
        {
            BapsDraw d = baps.Run(BapsScheme.Gibbs, chains: 1, masterSeed: 1000 + s, burn: 0, samples: 2).PosteriorDraws[0][0];
            var dev = new double[nu];
            for (int j = 0; j < nu; j++) dev[j] = d.Coefficients[j] - mean[j];
            var w = new double[nu];
            for (int i = 0; i < nu; i++) for (int k = i; k < nu; k++) w[i] += l[k, i] * dev[k];   // Lᵀ(β − m) ~ N(0, I)
            white[s] = w;
            // τ² | β ~ IG(a + (ν − r)/2, b + βᵀPβ/2) with r the nullity, so (b + βᵀPβ/2)/τ² ~ Gamma(a + (ν − r)/2, 1).
            scaled[s] = (1e-3 + 0.5 * new TensorPenalty(Nu1, Nu1).Roughness(d.Coefficients)) / d.Tau2;
        }
        for (int i = 0; i < nu; i++)
        {
            double m = white.Average(w => w[i]), v = white.Average(w => (w[i] - m) * (w[i] - m));
            Assert.True(Math.Abs(m) * Math.Sqrt(draws) < 5.0, $"coordinate {i}: mean {m}");
            Assert.True(Math.Abs(v - 1.0) < 5.0 * Math.Sqrt(2.0 / draws), $"coordinate {i}: variance {v}");
        }
        double shape = 1e-3 + 0.5 * (nu - 4), meanScaled = scaled.Average();
        Assert.True(Math.Abs(meanScaled - shape) < 5.0 * Math.Sqrt(shape / draws), $"Gamma mean {meanScaled} against shape {shape}");
    }

    [Fact]
    public void Baps_LocallyAdaptiveScheme_RefusesAPenaltyWithoutLocalTerms()
    {
        (double[,] z, double[] xs, double[] ys) = Surface(seed: 5, n: 60);
        double[] y = xs.Select((x, i) => SurfaceAt(x, ys[i])).ToArray();
        var baps = new Baps(z, y, new TensorPenalty(Nu1, Nu1));
        Assert.Throws<ArgumentException>(() => baps.Run(BapsScheme.LocallyAdaptiveGibbs, chains: 1, burn: 1, samples: 2));
        Assert.Equal(2 * 60, baps.Run(BapsScheme.Collapsed, chains: 2, masterSeed: 3, burn: 50, samples: 60).Draws);
    }

    [Fact]
    public void Hierarchical_OnATensorPenalty_UsesItsNullityForTheSmoothingVariance()
    {
        (double[,] z, double[] xs, double[] ys) = Surface(seed: 8, n: 80);
        var rng = new Xoshiro256PlusPlus(seed: 9);
        var reps = new double[3][];
        for (int r = 0; r < 3; r++) reps[r] = xs.Select((x, i) => SurfaceAt(x, ys[i]) + 0.1 * r + 0.15 * Gaussian(rng)).ToArray();
        var model = new HierarchicalPSpline(z, reps, new TensorPenalty(Nu1, Nu1));
        var penalty = new TensorPenalty(Nu1, Nu1);

        // τ_0² | β_0 ~ IG(a + (ν − r)/2, b + β_0ᵀPβ_0/2): the first sweep's draws, scaled, are Gamma(a + (ν − r)/2, 1).
        const int draws = 1500;
        int nu = Nu1 * Nu1;
        var scaled = new double[draws];
        for (int s = 0; s < draws; s++)
        {
            HierarchicalDraw d = model.Run(chains: 1, masterSeed: 2000 + s, burn: 0, samples: 1, start: HierarchicalStart.Fixed)
                                      .PosteriorDraws[0][0];
            scaled[s] = (1e-3 + 0.5 * penalty.Roughness(d.Population)) / d.PopulationSmoothingVariance;
        }
        double shape = 1e-3 + 0.5 * (nu - 4), meanScaled = scaled.Average();
        Assert.True(Math.Abs(meanScaled - shape) < 5.0 * Math.Sqrt(shape / draws), $"Gamma mean {meanScaled} against shape {shape}");

        // And the fit recovers the population surface.
        HierarchicalResult h = model.Run(chains: 2, masterSeed: 4, burn: 400, samples: 400);
        double mse = 0.0;
        for (int i = 0; i < xs.Length; i++)
        {
            double f = 0.0;
            for (int j = 0; j < nu; j++) f += z[i, j] * h.PopulationCoefficients[j];
            mse += (f - SurfaceAt(xs[i], ys[i]) - 0.1) * (f - SurfaceAt(xs[i], ys[i]) - 0.1);
        }
        Assert.True(mse / xs.Length < 0.01, $"population surface MSE {mse / xs.Length}");
    }
}
