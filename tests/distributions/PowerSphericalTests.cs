using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Distributions.Tests;

/// <summary>
/// E36's slice 0 gates on <see cref="PowerSpherical"/> (the BPPR brief). Its draws against its density, cell by cell on the circle and
/// the sphere, each cell's probability by quadrature of the density itself over cells set off-centre from μ; the cosine's law against
/// its closed form, (1 − μᵀx)/2 ~ Beta(β, α), and the tangent's against the uniform's, by Kolmogorov–Smirnov, on chosen κ and a seeded
/// sweep to κ = 10⁵ and d = 16 with means at the reflection's branch and next to ±e₁; the normalizer against the beta's through the
/// tangent-normal decomposition; and the proposal's symmetry, in the density bit for bit and in the chain it drives: started uniform,
/// one step leaves the uniform law in place and the pair exchangeable. The seeds are fixed, so a failure repeats.
/// </summary>
public sealed class PowerSphericalTests(ITestOutputHelper output)
{
    private const double KolmogorovCut = 2.6;   // P(√n·D > 2.6) ≈ 3 in a million per statistic
    private const double SignCut = 5.0;         // |z| of a sign count, about 6 in ten million
    private const double CellFloor = 1e-6;      // the smallest χ² p-value a cell test accepts
    private const double Eps = 2.220446049250313e-16;

    // ── The cosine and the tangent ──────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(2, 0.0)]
    [InlineData(2, 0.5)]
    [InlineData(2, 1000.0)]
    [InlineData(3, 0.0)]
    [InlineData(3, 4.0)]
    [InlineData(3, 1000.0)]
    [InlineData(3, 1e5)]
    [InlineData(5, 0.5)]
    [InlineData(5, 1000.0)]
    [InlineData(8, 30.0)]
    [InlineData(16, 1e5)]
    public void CosineAndTangent_FollowTheirLaws(int d, double kappa)
    {
        var rng = new Xoshiro256PlusPlus(1000 * d + (int)Math.Round(Math.Log10(kappa + 1.0) * 10));
        var law = PowerSpherical.Of(Unit(rng, d), kappa);
        Laws laws = Read(rng, law, 5000);
        output.WriteLine($"d = {d}, κ = {kappa}: {laws}");
        Assert.True(laws.Worst <= KolmogorovCut && laws.Sign <= SignCut, laws.ToString());
    }

    [Fact]
    public void CosineAndTangent_SeededSweep()
    {
        var rng = new Xoshiro256PlusPlus(20261114);
        var failures = new List<string>();
        double worst = 0.0, worstSign = 0.0, worstNorm = 0.0;
        string[] means = ["random", "a coordinate axis", "at the reflection's branch", "next to ±e₁"];
        for (int c = 0; c < 300; c++)
        {
            int d = 2 + rng.NextInt(15), kind = c % 4;
            double u = rng.NextDouble();
            double kappa = u < 0.1 ? 0.0 : u < 0.2 ? 1000.0 : Math.Pow(10.0, -3.0 + 8.0 * rng.NextDouble());
            double[] mu = kind switch
            {
                0 => Unit(rng, d),
                1 => Axis(d, rng.NextInt(d), rng.NextDouble() < 0.5 ? -1.0 : 1.0),
                2 => Tilted(rng, d, (rng.NextDouble() < 0.5 ? -1.0 : 1.0) * Math.Pow(10.0, -17.0 + 4.0 * rng.NextDouble())),
                _ => Tilted(rng, d, (rng.NextDouble() < 0.5 ? -1.0 : 1.0) * (1.0 - Math.Pow(10.0, -14.0 + 8.0 * rng.NextDouble()))),
            };
            var law = PowerSpherical.Of(mu, kappa);
            Laws laws = Read(rng, law, 2000);
            worst = Math.Max(worst, laws.Worst);
            worstSign = Math.Max(worstSign, laws.Sign);
            worstNorm = Math.Max(worstNorm, laws.NormError);
            if (laws.Worst > KolmogorovCut || laws.Sign > SignCut || !(laws.NormError <= 8 * d * Eps))
                failures.Add($"case {c}: d = {d}, κ = {kappa:R}, mean {means[kind]}: {laws}");
        }
        output.WriteLine($"300 cases: worst √n·D {worst:F2}, worst sign |z| {worstSign:F2}, worst |‖x‖² − 1| {worstNorm:E2}");
        Assert.True(failures.Count == 0, string.Join("\n", failures.Take(20)));
    }

    // ── The draws against the density, on cells ─────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(2, 0.0)]
    [InlineData(2, 0.5)]
    [InlineData(2, 3.0)]
    [InlineData(2, 30.0)]
    [InlineData(2, 1000.0)]
    [InlineData(2, 1e5)]
    [InlineData(3, 0.0)]
    [InlineData(3, 0.5)]
    [InlineData(3, 3.0)]
    [InlineData(3, 30.0)]
    [InlineData(3, 1000.0)]
    [InlineData(3, 1e5)]
    public void Draws_MatchTheDensity_OnCells(int d, double kappa)
    {
        var rng = new Xoshiro256PlusPlus(77 + d * 101 + (int)Math.Round(kappa * 7.0));
        double[] mu = Unit(rng, d);
        var law = PowerSpherical.Of(mu, kappa);
        double meanW = 0.5 * (d - 1) / (d - 1 + kappa);   // E[(1 − μᵀx)/2]
        double r = Math.Min(Math.PI / 2, 2.0 * Math.Sqrt(meanW));   // the bulk's angular scale
        Cells cells = d == 2 ? CircleCells(law, r) : SphereCells(rng, law, r);
        // On the sphere the antipode is a cone point of the density, |x + μ|^(2κ), inside a panel: there the rule holds to about 10⁻⁷.
        double held = d == 2 ? 1e-10 : 1e-6;
        Assert.True(Math.Abs(cells.Probabilities.Sum() - 1.0) <= held, $"the cells' quadrature sums to {cells.Probabilities.Sum():R}");

        const int n = 20_000;
        var counts = new int[cells.Probabilities.Length];
        var x = new double[d];
        for (int i = 0; i < n; i++)
        {
            law.Sample(rng, x);
            counts[cells.Locate(x)]++;
        }
        (double statistic, int freedom) = ChiSquare(counts, cells.Probabilities, n);
        double p = SpecialFunctions.RegularizedGammaQ(0.5 * freedom, 0.5 * statistic);
        output.WriteLine($"d = {d}, κ = {kappa}: χ² = {statistic:F1} on {freedom} degrees of freedom, p = {p:E2}; cells sum to {cells.Probabilities.Sum():R}");
        Assert.True(p >= CellFloor, $"χ² = {statistic:F1} on {freedom}: p = {p:E2}");
    }

    // ── The normalizer ──────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Through the tangent-normal decomposition dσ_(d−1) = (1 − t²)^((d−3)/2) dt dσ_(d−2), the density at cosine t times |S^(d−2)|(1 − t²)^((d−3)/2)
    /// is the cosine's density, ½·Beta(α, β) at (1 + t)/2: read at w = (1 − t)/2 as Beta(β, α) at w, so neither side loses w's digits.
    /// </summary>
    [Fact]
    public void Normalizer_MatchesTheBetaThroughTheDecomposition()
    {
        var rng = new Xoshiro256PlusPlus(5);
        double[] ws = [1e-12, 1e-6, 0.01, 0.3, 0.5, 0.75, 0.95];
        double[] kappas = [0.0, 1e-3, 0.5, 1.0, 7.5, 100.0, 1000.0, 1e4, 1e5];
        double worst = 0.0;
        for (int d = 2; d <= 40; d++)
            foreach (double kappa in kappas)
            {
                double[] mu = Unit(rng, d), e = Orthogonal(rng, mu);
                var law = PowerSpherical.Of(mu, kappa);
                double alpha = 0.5 * (d - 1) + kappa, beta = 0.5 * (d - 1);
                double logSphere = Math.Log(2.0) + beta * Math.Log(Math.PI) - SpecialFunctions.LogGamma(beta);   // |S^(d−2)|
                Beta cosine = Beta.Shapes(beta, alpha);
                foreach (double w in ws)
                {
                    double c = 2.0 * Math.Sqrt(w * (1.0 - w));
                    double[] x = [.. mu.Select((m, j) => (1.0 - 2.0 * w) * m + c * e[j])];
                    double left = law.LogPdf(x) + logSphere + 0.5 * (d - 3) * Math.Log(4.0 * w * (1.0 - w));
                    double right = -Math.Log(2.0) + cosine.LogPdf(w);
                    double scale = Math.Abs(SpecialFunctions.LogGamma(alpha)) + Math.Abs(SpecialFunctions.LogGamma(beta))
                                 + Math.Abs(SpecialFunctions.LogGamma(alpha + beta)) + Math.Abs(kappa * Math.Log(2.0 * (1.0 - w)))
                                 + Math.Abs((beta - 1.0) * Math.Log(w)) + Math.Abs((alpha - 1.0) * Math.Log(1.0 - w)) + 1.0;
                    double ratio = Math.Abs(left - right) / (64 * Eps * scale);
                    worst = Math.Max(worst, ratio);
                    Assert.True(ratio <= 1.0, $"d = {d}, κ = {kappa}, w = {w}: {left:R} against {right:R}");
                }
                if (kappa == 0.0)   // uniform: 1/|S^(d−1)|
                    Assert.Equal(-(Math.Log(2.0) + 0.5 * d * Math.Log(Math.PI) - SpecialFunctions.LogGamma(0.5 * d)), law.LogPdf(mu), 12);
            }
        output.WriteLine($"worst |left − right| / tolerance {worst:F3}");
    }

    [Fact]
    public void OnTheZeroSphere_ItKeepsTheSign()
    {
        var rng = new Xoshiro256PlusPlus(11);
        foreach (double sign in new[] { -1.0, 1.0 })
        {
            var law = PowerSpherical.Of([sign], 1000.0);
            Assert.Equal(0.0, law.LogPdf([sign]));
            Assert.Equal(double.NegativeInfinity, law.LogPdf([-sign]));
            for (int i = 0; i < 1000; i++) Assert.Equal(sign, law.Sample(rng)[0]);
            Assert.Equal(sign, PowerSpherical.Of([sign], 1e-9).Sample(rng)[0]);

            var uniform = PowerSpherical.Of([sign], 0.0);
            Assert.Equal(-Math.Log(2.0), uniform.LogPdf([1.0]), 15);
            Assert.Equal(-Math.Log(2.0), uniform.LogPdf([-1.0]), 15);
            int up = 0;
            for (int i = 0; i < 4000; i++) up += uniform.Sample(rng)[0] > 0 ? 1 : 0;
            Assert.True(Math.Abs(2.0 * up - 4000) / Math.Sqrt(4000) <= SignCut, $"{up} of 4000 up");
            Assert.Throws<InvalidOperationException>(() => law.CosineLaw);
        }
    }

    // ── The proposal's symmetry ─────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Density_IsSymmetricInMeanAndPoint_BitForBit()
    {
        var rng = new Xoshiro256PlusPlus(13);
        for (int c = 0; c < 4000; c++)
        {
            int d = 1 + rng.NextInt(12);
            double kappa = c % 5 == 0 ? 0.0 : Math.Pow(10.0, -3.0 + 8.0 * rng.NextDouble());
            double[] mu = Unit(rng, d);
            double[] x = c % 3 == 0 ? PowerSpherical.Of(mu, kappa).Sample(rng) : Unit(rng, d);
            Assert.Equal(PowerSpherical.Of(mu, kappa).LogPdf(x), PowerSpherical.Of(x, kappa).LogPdf(mu));
        }
    }

    /// <summary>
    /// A symmetric proposal under a uniform target is reversible, so with θ uniform and θ* ~ PS(θ, κ) the pair is exchangeable and θ*
    /// is uniform: θ*'s first coordinate against the uniform's cosine law, and three statistics odd under the swap, each with
    /// P(> 0) = ½ — h(θ) − h(θ*) for two functions h, and θ₁θ*₂ − θ*₁θ₂, which a proposal that turned one way would bias.
    /// </summary>
    [Theory]
    [InlineData(2, 2.0)]
    [InlineData(3, 20.0)]
    [InlineData(5, 1.0)]
    [InlineData(3, 1000.0)]
    [InlineData(4, 1e5)]
    public void OneStepFromUniform_LeavesItUniformAndThePairExchangeable(int d, double kappa)
    {
        var rng = new Xoshiro256PlusPlus(900 + d);
        const int n = 40_000;
        var first = new double[n];
        var signs = new int[3];
        var target = new double[d];
        for (int i = 0; i < n; i++)
        {
            double[] theta = Unit(rng, d);
            PowerSpherical.Of(theta, kappa).Sample(rng, target);
            first[i] = 0.5 * (1.0 + target[0]);
            double H(double[] v) => v[0] + 0.7 * v[1] * v[1];
            double K(double[] v) => v[0] * v[0] * v[0] - v[1] + 0.2 * v[d - 1];
            if (H(theta) > H(target)) signs[0]++;
            if (K(theta) > K(target)) signs[1]++;
            if (theta[0] * target[1] > target[0] * theta[1]) signs[2]++;
        }
        double stationary = Kolmogorov(first, Beta.Shapes(0.5 * (d - 1), 0.5 * (d - 1)).Cdf);
        double[] z = [.. signs.Select(s => Math.Abs(2.0 * s - n) / Math.Sqrt(n))];
        output.WriteLine($"d = {d}, κ = {kappa}: θ* uniform √n·D = {stationary:F2}; sign |z| = {string.Join(", ", z.Select(v => v.ToString("F2")))}");
        Assert.True(stationary <= KolmogorovCut, $"θ*'s first coordinate: √n·D = {stationary:F2}");
        Assert.True(z.All(v => v <= SignCut), $"sign |z| = {string.Join(", ", z)}");
    }

    // ── Refusals ────────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Of_RefusesWhatIsNotALaw()
    {
        Assert.Throws<ArgumentException>(() => PowerSpherical.Of([], 1.0));
        Assert.Throws<ArgumentException>(() => PowerSpherical.Of([0.6, 0.81], 1.0));
        Assert.Throws<ArgumentException>(() => PowerSpherical.Of([double.NaN, 1.0], 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => PowerSpherical.Of([0.6, 0.8], -1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => PowerSpherical.Of([0.6, 0.8], double.PositiveInfinity));
        Assert.Throws<ArgumentOutOfRangeException>(() => PowerSpherical.Of([0.6, 0.8], double.NaN));
        var law = PowerSpherical.Of([0.6, 0.8], 3.0);
        Assert.Throws<ArgumentException>(() => law.LogPdf([1.0, 0.0, 0.0]));
        Assert.Throws<ArgumentException>(() => law.LogPdf([1.0, 0.1]));
        Assert.Throws<ArgumentException>(() => law.Sample(new Xoshiro256PlusPlus(1), new double[3]));
        Assert.Throws<ArgumentNullException>(() => law.Sample(null!, new double[2]));
        // Rounding is not refused: a mean a few ulps off the sphere is kept as given.
        Assert.Equal(1.0 + 4e-16, PowerSpherical.Of([1.0 + 4e-16], 1.0).Mean[0]);
    }

    // ── Reading a law's draws ───────────────────────────────────────────────────────────────────────────────────────

    /// <summary>√n·D of the cosine and of the tangent along a fixed direction orthogonal to μ, over all draws and within each half split at
    /// the cosine's median (the tangent is independent of the cosine); at d = 2 the tangent is a sign, read as a count.</summary>
    private sealed record Laws(double Cosine, double Tangent, double TangentLow, double TangentHigh, double Sign, double NormError)
    {
        public double Worst => Math.Max(Math.Max(Cosine, Tangent), Math.Max(TangentLow, TangentHigh));

        public override string ToString() =>
            $"cosine {Cosine:F2}, tangent {Tangent:F2} (low half {TangentLow:F2}, high half {TangentHigh:F2}), sign |z| {Sign:F2}, |‖x‖² − 1| {NormError:E1}";
    }

    private static Laws Read(Xoshiro256PlusPlus rng, PowerSpherical law, int n)
    {
        int d = law.Dimension;
        double[] mu = law.Mean.ToArray(), e = Orthogonal(rng, mu);
        var w = new double[n];
        var s = new double[n];
        var x = new double[d];
        double normError = 0.0;
        for (int i = 0; i < n; i++)
        {
            law.Sample(rng, x);
            double norm = 0.0, gap = 0.0, t = 0.0;
            for (int j = 0; j < d; j++)
            {
                norm += x[j] * x[j];
                gap += (x[j] - mu[j]) * (x[j] - mu[j]);
                t += x[j] * mu[j];
            }
            Assert.True(x.All(double.IsFinite), "a draw has a non-finite coordinate");
            normError = Math.Max(normError, Math.Abs(norm - 1.0));
            w[i] = 0.25 * gap;   // (1 − μᵀx)/2 for unit x and μ
            double tangent = 0.0, along = 0.0;
            for (int j = 0; j < d; j++)
            {
                double v = x[j] - t * mu[j];
                tangent += v * v;
                along += v * e[j];
            }
            s[i] = tangent > 0.0 ? along / Math.Sqrt(tangent) : 0.0;
        }
        Beta cosine = Beta.Shapes(0.5 * (d - 1), 0.5 * (d - 1) + law.Concentration);
        double cosineStatistic = Kolmogorov((double[])w.Clone(), cosine.Cdf);
        if (d == 2)
        {
            int up = s.Count(v => v > 0.0);
            return new Laws(cosineStatistic, 0.0, 0.0, 0.0, Math.Abs(2.0 * up - n) / Math.Sqrt(n), normError);
        }
        Beta tangentLaw = Beta.Shapes(0.5 * (d - 2), 0.5 * (d - 2));
        double[] half = [.. s.Select(v => 0.5 * (1.0 + v))];
        double median = w.OrderBy(v => v).ElementAt(n / 2);
        double[] low = [.. half.Where((_, i) => w[i] < median)], high = [.. half.Where((_, i) => w[i] >= median)];
        return new Laws(cosineStatistic, Kolmogorov(half, tangentLaw.Cdf), Kolmogorov(low, tangentLaw.Cdf), Kolmogorov(high, tangentLaw.Cdf), 0.0, normError);
    }

    // √n·D, read at each distinct value's last rank so rounding's atoms count once.
    private static double Kolmogorov(double[] xs, Func<double, double> cdf)
    {
        Array.Sort(xs);
        int n = xs.Length, first = 0;
        double gap = 0.0;
        for (int i = 0; i < n; i++)
        {
            if (i + 1 < n && xs[i + 1] == xs[i]) continue;
            double f = cdf(xs[i]);
            gap = Math.Max(gap, Math.Max(Math.Abs(f - (i + 1.0) / n), Math.Abs(f - (double)first / n)));
            first = i + 1;
        }
        return Math.Sqrt(n) * gap;
    }

    // ── Cells ───────────────────────────────────────────────────────────────────────────────────────────────────────

    private sealed record Cells(double[] Probabilities, Func<double[], int> Locate);

    /// <summary>On the circle, 12 arcs over μ ± min(π, 5r), set off μ by 0.37 of an arc, and the rest of the circle as one more; each
    /// arc's probability by Gauss–Legendre on panels no wider than r/4.</summary>
    private static Cells CircleCells(PowerSpherical law, double r)
    {
        double[] mu = law.Mean.ToArray();
        double muAngle = Math.Atan2(mu[1], mu[0]);
        double window = Math.Min(Math.PI, 5.0 * r);
        const int arcs = 12;
        double width = 2.0 * window / arcs, start = muAngle - window + 0.37 * width;
        int count = window < Math.PI ? arcs + 1 : arcs;
        var edges = new double[count + 1];
        for (int k = 0; k <= arcs; k++) edges[k] = start + k * width;
        if (count > arcs) edges[count] = start + 2.0 * Math.PI;
        else edges[arcs] = start + 2.0 * Math.PI;   // the arcs cover the circle
        var point = new double[2];
        double Density(double phi)
        {
            point[0] = Math.Cos(phi);
            point[1] = Math.Sin(phi);
            return Math.Exp(law.LogPdf(point));
        }
        var probabilities = new double[count];
        // The density is |x + μ|^(2κ) near the antipode, a kink for κ = ½: an arc that holds it is split there.
        double antipode = start + (muAngle + Math.PI - start - 2.0 * Math.PI * Math.Floor((muAngle + Math.PI - start) / (2.0 * Math.PI)));
        for (int k = 0; k < count; k++)
            probabilities[k] = edges[k] < antipode && antipode < edges[k + 1]
                ? Integrate(Density, edges[k], antipode, r / 4.0) + Integrate(Density, antipode, edges[k + 1], r / 4.0)
                : Integrate(Density, edges[k], edges[k + 1], r / 4.0);
        int Locate(double[] x)
        {
            double phi = Math.Atan2(x[1], x[0]) - start;
            phi -= 2.0 * Math.PI * Math.Floor(phi / (2.0 * Math.PI));
            int k = (int)Math.Floor(phi / width);
            return Math.Min(k, count - 1);
        }
        return new Cells(probabilities, Locate);
    }

    /// <summary>On the sphere, cells about an axis e set 0.4r off μ: five polar bands over [0, min(π, 0.4r + 4r)] and the rest as a sixth,
    /// each cut into six azimuthal sectors turned 0.3 rad; each cell's probability by Gauss–Legendre in the polar angle and the azimuth
    /// against sin ϑ dϑ dφ.</summary>
    private static Cells SphereCells(Xoshiro256PlusPlus rng, PowerSpherical law, double r)
    {
        double[] mu = law.Mean.ToArray(), f0 = Orthogonal(rng, mu);
        double offset = 0.4 * r;
        double[] e = [.. mu.Select((m, j) => Math.Cos(offset) * m + Math.Sin(offset) * f0[j])];
        double[] f1 = Orthogonal(rng, e);
        double[] f2 = [e[1] * f1[2] - e[2] * f1[1], e[2] * f1[0] - e[0] * f1[2], e[0] * f1[1] - e[1] * f1[0]];
        double window = Math.Min(Math.PI, offset + 4.0 * r);
        const int bands = 5, sectors = 6;
        const double turn = 0.3;
        var polar = new List<double>();
        for (int k = 0; k <= bands; k++) polar.Add(window * k / bands);
        if (window < Math.PI) polar.Add(Math.PI);
        int rings = polar.Count - 1;
        double sector = 2.0 * Math.PI / sectors;
        var point = new double[3];
        double Density(double theta, double phi)
        {
            double st = Math.Sin(theta), ct = Math.Cos(theta), cp = Math.Cos(phi), sp = Math.Sin(phi);
            for (int j = 0; j < 3; j++) point[j] = ct * e[j] + st * (cp * f1[j] + sp * f2[j]);
            return Math.Exp(law.LogPdf(point));
        }
        // Past 12r from μ the density is below e^(−κ(12r)²/4) of its peak, about e^(−140) wherever 12r < π: there the panels widen to
        // an eighth of a radian.
        double far = 12.0 * r + offset;
        double Step(double theta) => theta < far ? r / 4.0 : 0.125;
        var probabilities = new double[rings * sectors];
        for (int b = 0; b < rings; b++)
            for (int a = 0; a < sectors; a++)
            {
                double phiLow = turn + a * sector, phiHigh = phiLow + sector;
                double Ring(double theta) => Math.Sin(theta) * Integrate(phi => Density(theta, phi), phiLow, phiHigh, Math.Max(r / 4.0, sector / 8.0));
                double low = polar[b], high = polar[b + 1];
                probabilities[b * sectors + a] = low < far && high > far
                    ? Integrate(Ring, low, far, Step(low)) + Integrate(Ring, far, high, Step(far))
                    : Integrate(Ring, low, high, Step(low));
            }
        int Locate(double[] x)
        {
            double along = 0.0, c1 = 0.0, c2 = 0.0;
            for (int j = 0; j < 3; j++)
            {
                along += x[j] * e[j];
                c1 += x[j] * f1[j];
                c2 += x[j] * f2[j];
            }
            double theta = Math.Atan2(Math.Sqrt(c1 * c1 + c2 * c2), along);
            double phi = Math.Atan2(c2, c1) - turn;
            phi -= 2.0 * Math.PI * Math.Floor(phi / (2.0 * Math.PI));
            int ring = 0;
            while (ring < rings - 1 && theta >= polar[ring + 1]) ring++;
            return ring * sectors + Math.Min((int)Math.Floor(phi / sector), sectors - 1);
        }
        return new Cells(probabilities, Locate);
    }

    // Composite 16-node Gauss–Legendre on panels no wider than the step.
    private static double Integrate(Func<double, double> f, double lower, double upper, double step)
    {
        GaussLegendre rule = GaussLegendre.Of(16);
        int panels = Math.Max(1, (int)Math.Ceiling((upper - lower) / step));
        double width = (upper - lower) / panels, sum = 0.0;
        for (int q = 0; q < panels; q++)
        {
            double mid = lower + (q + 0.5) * width, half = 0.5 * width;
            for (int i = 0; i < rule.Nodes.Length; i++) sum += half * rule.Weights[i] * f(mid + half * rule.Nodes[i]);
        }
        return sum;
    }

    // Pearson's χ², the cells expecting fewer than five pooled into one.
    private static (double Statistic, int Freedom) ChiSquare(int[] counts, double[] probabilities, int n)
    {
        double statistic = 0.0, pooledExpected = 0.0;
        int pooledObserved = 0, cells = 0;
        for (int k = 0; k < counts.Length; k++)
        {
            double expected = n * probabilities[k];
            if (expected < 5.0)
            {
                pooledExpected += expected;
                pooledObserved += counts[k];
                continue;
            }
            statistic += (counts[k] - expected) * (counts[k] - expected) / expected;
            cells++;
        }
        if (pooledExpected > 0.0 || pooledObserved > 0)
        {
            statistic += (pooledObserved - pooledExpected) * (pooledObserved - pooledExpected) / Math.Max(pooledExpected, 1e-300);
            cells++;
        }
        return (statistic, cells - 1);
    }

    // ── Points ──────────────────────────────────────────────────────────────────────────────────────────────────────

    private static double[] Unit(Xoshiro256PlusPlus rng, int d)
    {
        double[] v = [.. Enumerable.Range(0, d).Select(_ => Normal.SampleStandard(rng))];
        double norm = Math.Sqrt(v.Sum(c => c * c));
        return [.. v.Select(c => c / norm)];
    }

    private static double[] Axis(int d, int j, double sign)
    {
        var v = new double[d];
        v[j] = sign;
        return v;
    }

    // A unit vector with first coordinate c and the rest uniform on the sphere of radius √(1 − c²), read as √((1 − |c|)(1 + |c|)).
    private static double[] Tilted(Xoshiro256PlusPlus rng, int d, double c)
    {
        double rest = Math.Sqrt((1.0 - Math.Abs(c)) * (1.0 + Math.Abs(c)));
        double[] tail = Unit(rng, d - 1);
        return [c, .. tail.Select(v => rest * v)];
    }

    // A unit vector orthogonal to a unit v: a random one less its component along v, twice.
    private static double[] Orthogonal(Xoshiro256PlusPlus rng, double[] v)
    {
        double[] u = Unit(rng, v.Length);
        for (int pass = 0; pass < 2; pass++)
        {
            double dot = u.Zip(v, (a, b) => a * b).Sum();
            for (int j = 0; j < u.Length; j++) u[j] -= dot * v[j];
        }
        double norm = Math.Sqrt(u.Sum(c => c * c));
        return [.. u.Select(c => c / norm)];
    }
}
