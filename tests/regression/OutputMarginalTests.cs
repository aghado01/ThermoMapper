using System;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E34's engine pieces for Σ (the multi-output brief, slice 1) against slice 0's references, which share none of their code:
/// <list type="bullet">
/// <item><see cref="InverseWishart"/>'s density against <see cref="ReferenceInverseWishart"/>'s, itself normalized by quadrature, and
/// the inverse gamma's at K = 1; its draws' mean, Huang and Wand's correlation law at a small ν, and the inverse gamma's draws at
/// K = 1, to rounding on the same stream.</item>
/// <item><see cref="OutputMarginal"/> against <see cref="MatrixTReference"/> over a seeded sweep: K to 4, the three intercept
/// conventions with a frozen level, τ over seven decades, b to 0, units 10^±6, decision 2's prior, a declared one and the improper
/// |Σ|^(−(K + 1)/2), refusals agreeing; Σ's posterior against the reference's; and at K = 1 the scalar block's marginal.</item>
/// </list>
/// </summary>
public sealed class OutputMarginalTests(ITestOutputHelper output)
{
    private const int Seed = 20261109;

    // ── The inverse Wishart ───────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void InverseWishart_Density_IsTheReferences_AndTheInverseGammaAtOneOutput()
    {
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 1)[0]);
        double worst = 0.0;
        for (int index = 0; index < 300; index++)
        {
            int k = 1 + rng.NextInt(4);
            double nu = k - 1 + 0.01 + 20.0 * rng.NextDouble(), scale = Math.Pow(10.0, -6.0 + 12.0 * rng.NextDouble());
            double[,] psi = Spd(rng, k, 4.0 * rng.NextDouble(), scale);
            double[,] sigma = Spd(rng, k, 4.0 * rng.NextDouble(), scale * Math.Pow(10.0, -2.0 + 4.0 * rng.NextDouble()));
            double engine = InverseWishart.Of(nu, psi).LogPdf(sigma), reference = ReferenceInverseWishart.LogDensity(sigma, nu, psi);
            worst = Math.Max(worst, Math.Abs(engine - reference) / Math.Max(1.0, Math.Abs(reference)));
            if (k == 1)
                Close(InverseGamma.ShapeScale(0.5 * nu, 0.5 * psi[0, 0]).LogPdf(sigma[0, 0]), engine, 1e-13);
        }
        output.WriteLine($"300 cases: worst relative gap {worst:E2}");
        Assert.True(worst < 1e-12, $"worst relative gap {worst:E2}");
        Assert.Equal(double.NegativeInfinity, InverseWishart.Of(3.0, new[,] { { 1.0, 0.0 }, { 0.0, 1.0 } }).LogPdf(new[,] { { 1.0, 2.0 }, { 2.0, 1.0 } }));
    }

    [Fact]
    public void InverseWishart_Refuses_AnImproperShape_AndAScaleNotPositiveDefinite()
    {
        double[,] unit = { { 1.0, 0.0 }, { 0.0, 1.0 } };
        Assert.Throws<ArgumentOutOfRangeException>(() => InverseWishart.Of(1.0, unit));
        Assert.Throws<ArgumentOutOfRangeException>(() => InverseWishart.Of(double.PositiveInfinity, unit));
        Assert.Throws<ArgumentException>(() => InverseWishart.Of(3.0, new[,] { { 1.0, 2.0 }, { 2.0, 1.0 } }));
        Assert.Throws<ArgumentException>(() => InverseWishart.Of(3.0, new[,] { { 1.0, 0.1 }, { 0.2, 1.0 } }));
        Assert.Throws<ArgumentException>(() => InverseWishart.Of(3.0, new double[2, 3]));
        Assert.Null(InverseWishart.Of(3.0, unit).Mean);
    }

    [Fact]
    public void InverseWishart_Draws_HaveTheirMean_AndAreTheInverseGammasAtOneOutput()
    {
        // The mean at K = 3; at K = 1 the draw reads the inverse gamma's gamma variate on the same stream, so the two agree to rounding.
        double[,] psi = { { 0.6, -0.4, 0.01 }, { -0.4, 2.4, -0.02 }, { 0.01, -0.02, 0.004 } };
        InverseWishart iw = InverseWishart.Of(9.0, psi);
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 2)[1]);
        const int draws = 200_000;
        var sum = new double[3, 3];
        var squares = new double[3, 3];
        for (int s = 0; s < draws; s++)
        {
            double[,] sigma = iw.Sample(rng);
            for (int a = 0; a < 3; a++)
                for (int d = 0; d < 3; d++) { sum[a, d] += sigma[a, d]; squares[a, d] += sigma[a, d] * sigma[a, d]; }
        }
        double[,] expected = iw.Mean!;
        for (int a = 0; a < 3; a++)
            for (int d = 0; d < 3; d++)
            {
                double mean = sum[a, d] / draws, se = Math.Sqrt((squares[a, d] / draws - mean * mean) / draws);
                Assert.True(Math.Abs(mean - expected[a, d]) < 4.5 * se, $"E[Σ]({a}, {d}) = {mean:R} ± {se:E2} against {expected[a, d]:R}");
            }

        foreach ((double nu, double scale) in new[] { (0.5, 3.0), (3.0, 1e-6), (40.0, 2e5) })
        {
            InverseWishart one = InverseWishart.Of(nu, new[,] { { scale } });
            InverseGamma ig = InverseGamma.ShapeScale(0.5 * nu, 0.5 * scale);
            var left = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 3)[2]);
            var right = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 3)[2]);
            for (int s = 0; s < 1000; s++) Close(ig.Sample(right), one.Sample(left)[0, 0], 1e-14);
        }
    }

    [Fact]
    public void InverseWishart_Draws_KeepHuangAndWandsCorrelationLaw_AtASmallShape()
    {
        // ν = 1.5 at K = 2 is Huang and Wand's ν = 0.5: ρ² ∼ Beta(½, ¼), E[ρ²] = 2/3. A tiny χ² pivot leaves the Wishart singular in doubles
        // there; the draw through Ψ's factor stays exact.
        InverseWishart iw = InverseWishart.Of(1.5, new[,] { { 2.0 / 0.3, 0.0 }, { 0.0, 2.0 / 4.0 } });
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 4)[3]);
        const int draws = 200_000;
        double sum = 0.0, squares = 0.0;
        for (int s = 0; s < draws; s++)
        {
            double[,] sigma = iw.Sample(rng);
            double rho2 = sigma[1, 0] * sigma[1, 0] / (sigma[0, 0] * sigma[1, 1]);
            Assert.True(double.IsFinite(rho2), "a draw is not finite");
            sum += rho2;
            squares += rho2 * rho2;
        }
        double mean = sum / draws, se = Math.Sqrt((squares / draws - mean * mean) / draws);
        output.WriteLine($"E[ρ²] {mean:F5} ± {se:E1} against {2.0 / 3.0:F5}");
        Assert.True(Math.Abs(mean - 2.0 / 3.0) < 4.5 * se, $"E[ρ²] {mean:R} ± {se:E2}");
    }

    // ── The K-output marginal ─────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Marginal_IsTheMatrixTReference_OverASeededSweep()
    {
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 5)[4]);
        GPriorCentre[] centres = [GPriorCentre.Mean, GPriorCentre.Zero, GPriorCentre.Flat];
        double worst = 0.0, worstPosterior = 0.0;
        int compared = 0, refused = 0;
        for (int index = 0; index < 600; index++)
        {
            int k = 1 + rng.NextInt(4), n = 3 + rng.NextInt(10), m = rng.NextInt(Math.Min(4, n - 2) + 1);
            GPriorCentre centre = centres[rng.NextInt(3)];
            double[]? levels = null;
            double scale = Math.Pow(10.0, -6.0 + 12.0 * rng.NextDouble()), tau = Math.Pow(10.0, -4.0 + 7.0 * rng.NextDouble());
            double b = rng.NextDouble() switch { < 0.15 => 1.0, < 0.3 => centre == GPriorCentre.Flat ? 1.0 : 0.0, < 0.45 => Math.Pow(10.0, -3.0 * rng.NextDouble()) / n, _ => rng.NextDouble() };
            double[,] xc = ReferenceAlgebra.Centred(Matrix(rng, n, m, Math.Pow(10.0, -2.0 + 4.0 * rng.NextDouble())));
            var y = new double[n, k];
            for (int i = 0; i < n; i++)
                for (int l = 0; l < k; l++) y[i, l] = scale * ((l + 1) * 1.5 + Normal.SampleStandard(rng) + (m > 0 ? 0.7 * xc[i, 0] : 0.0));
            if (centre == GPriorCentre.Mean && rng.NextDouble() < 0.25)
            {
                levels = new double[k];
                for (int l = 0; l < k; l++) levels[l] = scale * Normal.SampleStandard(rng);
            }
            double[] means = ColumnMeans(y), centres0 = centre == GPriorCentre.Zero ? new double[k] : levels ?? means;
            int kind = rng.NextInt(3);   // decision 2's prior, a declared one, or the improper |Σ|^(−(K + 1)/2)
            // Decision 2's S_Y and Jeffreys' S are singular when the cases cannot span the outputs; decision 2 refuses that S_Y.
            if (kind != 1 && n - 1 < k + 1) kind = 1;
            if (kind == 2 && b == 0.0) kind = 0;
            ReferenceSigmaPrior prior = kind switch
            {
                0 => ReferenceSigmaPrior.ScaledToResponse(y, centres0, rng.NextDouble() < 0.5 ? 1e-3 : 0.5),
                1 => new ReferenceSigmaPrior(k - 1 + 0.05 + 10.0 * rng.NextDouble(), Spd(rng, k, 3.0 * rng.NextDouble(), scale * scale)),
                _ => new ReferenceSigmaPrior(0.0, null),
            };
            InverseWishart? enginePrior = prior.Scale is { } psi ? InverseWishart.Of(prior.Nu, psi) : null;

            (double[,] r, double[,] f) = LeastSquaresForms(y, xc);
            (double reference, bool referenceRefused) = Attempt(() => MatrixTReference.LogMarginal(MatrixTReference.Collapse(MultiOutputCase.GForm(y, xc, centre, tau, b, levels)), prior));
            (double engine, bool engineRefused) = Attempt(() => OutputMarginal.CentredAt(m, r, f, n, means, centre, tau, b, levels).LogMarginal(enginePrior));
            Assert.True(referenceRefused == engineRefused,
                        $"case {index}: reference {(referenceRefused ? "refused" : reference.ToString("R"))}, engine {(engineRefused ? "refused" : engine.ToString("R"))} ({centre}, K {k}, n {n}, b {b:R})");
            if (engineRefused) { refused++; continue; }
            compared++;
            double gap = Math.Abs(engine - reference) / Math.Max(1.0, Math.Abs(reference));
            Assert.True(gap < 1e-10, $"case {index}: engine {engine:R}, reference {reference:R} ({centre}, K {k}, n {n}, M {m}, τ {tau:E2}, b {b:R}, prior {kind})");
            worst = Math.Max(worst, gap);
            if (b > 0.0)
            {
                InverseWishart posterior = OutputMarginal.CentredAt(m, r, f, n, means, centre, tau, b, levels).SigmaPosterior(enginePrior);
                (double nu, double[,] scalePosterior) = MatrixTReference.SigmaPosterior(MatrixTReference.Collapse(MultiOutputCase.GForm(y, xc, centre, tau, b, levels)), prior);
                Close(nu, posterior.DegreesOfFreedom, 1e-14);
                double[,] engineScale = posterior.Scale;
                for (int a = 0; a < k; a++)
                    for (int d = 0; d < k; d++)
                        worstPosterior = Math.Max(worstPosterior, Math.Abs(engineScale[a, d] - scalePosterior[a, d]) / Math.Sqrt(Math.Abs(scalePosterior[a, a] * scalePosterior[d, d])));
            }
        }
        output.WriteLine($"{compared} compared, worst relative gap {worst:E2}; Σ's posterior scale {worstPosterior:E2}; {refused} refused by both");
        Assert.True(worstPosterior < 1e-11, $"Σ's posterior scale {worstPosterior:E2}");
        Assert.True(refused > 0, "the sweep reached no refusal");
    }

    [Fact]
    public void Marginal_AtOneOutput_IsTheScalarBlocks()
    {
        // IW(ν₀, Ψ₀) at K = 1 is IG(ν₀/2, Ψ₀/2); the improper |Σ|^(−1) is Jeffreys'.
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 6)[5]);
        GPriorCentre[] centres = [GPriorCentre.Mean, GPriorCentre.Zero, GPriorCentre.Flat];
        double worst = 0.0;
        for (int index = 0; index < 400; index++)
        {
            int n = 3 + rng.NextInt(10), m = rng.NextInt(Math.Min(4, n - 2) + 1);
            GPriorCentre centre = centres[rng.NextInt(3)];
            double tau = Math.Pow(10.0, -4.0 + 7.0 * rng.NextDouble()), b = centre == GPriorCentre.Flat ? 0.5 + 0.5 * rng.NextDouble() : 0.01 + 0.99 * rng.NextDouble();
            double scale = Math.Pow(10.0, -6.0 + 12.0 * rng.NextDouble());
            double[,] xc = ReferenceAlgebra.Centred(Matrix(rng, n, m, 1.0));
            var y = new double[n, 1];
            for (int i = 0; i < n; i++) y[i, 0] = scale * (2.0 + Normal.SampleStandard(rng));
            double? level = centre == GPriorCentre.Mean && rng.NextDouble() < 0.25 ? scale * Normal.SampleStandard(rng) : null;
            (double[,] r, double[,] f) = LeastSquaresForms(y, xc);
            double mean = ColumnMeans(y)[0], nu0 = 0.002 + 5.0 * rng.NextDouble(), psi0 = scale * scale * (0.01 + rng.NextDouble());
            bool jeffreys = rng.NextDouble() < 0.25;
            GPriorMarginal scalar = CollapsedGPrior.CentredMarginalAt(m, r[0, 0], f[0, 0], n, mean, centre, tau, b, level: level);
            OutputMarginal outputs = OutputMarginal.CentredAt(m, r, f, n, [mean], centre, tau, b, level is double l ? [l] : null);
            double expected = scalar.LogMarginal(jeffreys ? null : InverseGamma.ShapeScale(0.5 * nu0, 0.5 * psi0));
            double actual = outputs.LogMarginal(jeffreys ? null : InverseWishart.Of(nu0, new[,] { { psi0 } }));
            double gap = Math.Abs(actual - expected) / Math.Max(1.0, Math.Abs(expected));
            Assert.True(gap < 1e-12, $"case {index}: scalar {expected:R}, outputs {actual:R} ({centre}, b {b:R}, Jeffreys {jeffreys})");
            worst = Math.Max(worst, gap);
            Close(scalar.S, outputs.S[0, 0], 1e-15);
            Close(scalar.DegreesOfFreedom, outputs.DegreesOfFreedom, 1e-15);
        }
        output.WriteLine($"400 cases: worst relative gap {worst:E2}");
    }

    // ── Helpers ────────────────────────────────────────────────────────────────────────────────────────────────────────────

    // The centred least-squares fit's K × K residual form R = EᵀE and fitted form F = B_cᵀG_c⁻¹B_c, Y centred on its means.
    internal static (double[,] Residual, double[,] Fitted) LeastSquaresForms(double[,] y, double[,] xc)
    {
        int n = y.GetLength(0), k = y.GetLength(1), m = xc.GetLength(1);
        double[,] yc = ReferenceAlgebra.Centred(y);
        if (m == 0) return (ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(yc), yc), new double[k, k]);
        double[,] gram = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(xc), xc);
        double[,] cross = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(xc), yc);
        double[,] beta = ReferenceAlgebra.Solve(ReferenceAlgebra.Cholesky(gram)!, cross);
        double[,] fit = ReferenceAlgebra.Multiply(xc, beta);
        var e = new double[n, k];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < k; l++) e[i, l] = yc[i, l] - fit[i, l];
        double[,] r = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(e), e);
        double[,] f = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(cross), beta);
        ReferenceAlgebra.Symmetrize(r);
        ReferenceAlgebra.Symmetrize(f);
        return (r, f);
    }

    private static (double Value, bool Refused) Attempt(Func<double> f)
    {
        try { return (f(), false); }
        catch (InvalidOperationException) { return (double.NaN, true); }
        catch (ArgumentOutOfRangeException) { return (double.NaN, true); }
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

    private static double[,] Matrix(Xoshiro256PlusPlus rng, int rows, int columns, double scale)
    {
        var a = new double[rows, columns];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++) a[i, j] = scale * Normal.SampleStandard(rng);
        return a;
    }

    // Q·diag(d)·Qᵀ, Q a random rotation, d from scale down to scale·10^(−logCondition).
    private static double[,] Spd(Xoshiro256PlusPlus rng, int k, double logCondition, double scale)
    {
        double[,] q = Matrix(rng, k, k, 1.0);
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
