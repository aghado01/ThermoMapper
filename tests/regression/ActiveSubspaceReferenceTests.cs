using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 0 gates (the active-subspace brief) on its reference functions (<see cref="ReferenceFunction"/>): each gradient against
/// fourth-order central differences of its value, and each closed-form C against <see cref="GradientQuadrature"/> of its gradient to
/// 10⁻¹² of C's largest entry, on chosen cases (ridges and their sums under a box and under a correlated Gaussian, polynomials, a
/// Legendre expansion, and Rumsey et al.'s benchmarks from arXiv 2307.11241, whose stated C is asserted) and on a seeded sweep whose
/// laws reach 10⁴ from the origin and widths down to 10⁻³. In the sweep the tolerance adds what rounding the arguments at the law's
/// far end costs the quadrature (<see cref="Conditioning"/>); both numbers are printed.
/// </summary>
public sealed class ActiveSubspaceReferenceTests(ITestOutputHelper output)
{
    private const double Tolerance = 1e-12;
    private const double Eps = 1.1102230246251565e-16;

    private static ReferenceBox Box(double[] lower, double[] upper) => new(lower, upper);

    private static readonly double[,] Correlated =
    {
        { 1.0, 0.6, 0.2 },
        { 0.6, 2.0, -0.3 },
        { 0.2, -0.3, 0.5 },
    };

    private static readonly double[,] Symmetric =
    {
        { 2.0, 0.5, -1.0 },
        { 0.5, 1.0, 0.3 },
        { -1.0, 0.3, 3.0 },
    };

    // Rumsey et al. (arXiv 2307.11241), §4.1: f(x) = x₁² + x₁x₂ + x₂³/9, with p − 2 inert inputs; §4.1.2's forms t₁, t₂ for
    // g(x) = f(t₁ᵀx, t₂ᵀx) on [0, 1]⁶; and the supplement's f₃(x) = x₁² + x₁x₂x₃ + x₂³/9, with p − 3 inert inputs.
    private static readonly double[] T1 = [4.0, 3.0, 2.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
    private static readonly double[] T2 = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 5.0, 2.0, 2.0];

    private static readonly Dictionary<string, Func<(ReferenceFunction F, ReferenceMeasure Law)>> Chosen = new()
    {
        ["sin ridge"] = () => (new RidgeSum("sin ridge", new RidgeTerm(ExponentialProfile.Sin(2.0 * Math.PI), [1.0, 0.5, -0.25], 0.3)),
                               ReferenceBox.Cube(3, 0.0, 1.0)),
        ["exp ridge"] = () => (new RidgeSum("exp ridge", new RidgeTerm(ExponentialProfile.Exp(0.8, 1.5), [0.6, -1.2, 0.3], 1.0)),
                               Box([-1.0, 0.0, 3.0], [2.0, 1.0, 4.0])),
        ["damped wave"] = () => (new RidgeSum("damped wave",
                                     new RidgeTerm(new ExponentialProfile(new Complex(1.0, 0.5), new Complex(-0.4, 3.0)), [0.57735, 0.57735, -0.57735])),
                                 ReferenceBox.Cube(3, 0.0, 1.0)),
        ["two ridges, four inputs"] = () => (new RidgeSum("two ridges, four inputs",
                                                 new RidgeTerm(ExponentialProfile.Sin(3.0), [1.0, -1.0, 0.5, 0.0]),
                                                 new RidgeTerm(ExponentialProfile.Exp(1.0, 0.5), [0.2, 0.4, 0.0, 1.0], 1.0)),
                                             Box([0.0, -1.0, 0.0, 0.5], [1.0, 1.0, 2.0, 1.5])),
        ["three ridges, five inputs"] = () => (new RidgeSum("three ridges, five inputs",
                                                   new RidgeTerm(ExponentialProfile.Cos(2.0), [0.5, 0.5, 0.5, 0.5, 0.0]),
                                                   new RidgeTerm(ExponentialProfile.Sin(1.0), [1.0, -1.0, 0.0, 0.0, 0.3]),
                                                   new RidgeTerm(ExponentialProfile.Exp(0.5), [0.0, 0.0, 1.0, -1.0, 1.0])),
                                               ReferenceBox.Cube(5, 0.0, 1.0)),
        ["two ridges, correlated Gaussian"] = () => (new RidgeSum("two ridges, correlated Gaussian",
                                                         new RidgeTerm(ExponentialProfile.Sin(1.2), [0.8, -0.6, 0.0]),
                                                         new RidgeTerm(ExponentialProfile.Exp(0.4, 0.3), [0.3, 0.3, 0.9], 0.5)),
                                                     new ReferenceGaussian([1.0, -2.0, 0.5], Correlated)),
        ["quadratic form"] = () => (PolynomialReference.Quadratic("quadratic form", Symmetric, [0.2, -0.1, 1.0]),
                                    Box([-1.0, 0.0, 1.0], [1.0, 2.0, 4.0])),
        ["cubic ridge"] = () => (PolynomialReference.Ridge("cubic ridge", [0.7, -0.2, 1.1], 0.4, [1.0, -2.0, 0.0, 1.0], [1.5, -2.0, 0.25]),
                                 Box([1.0, -3.0, 0.0], [2.0, -1.0, 0.5])),
        ["interacting polynomial"] = () => (new PolynomialReference("interacting polynomial",
                                                [new Monomial(3.0, [2, 1, 0]), new Monomial(-2.0, [0, 1, 3]), new Monomial(1.0, [1, 1, 1]), new Monomial(0.5, [0, 0, 4])],
                                                [0.5, 0.0, -0.2]),
                                            Box([0.0, -1.0, -0.5], [1.0, 1.0, 0.5])),
        ["Rumsey et al.'s benchmark"] = () => (Benchmark(2), ReferenceBox.Cube(2, 0.0, 1.0)),
        ["Rumsey et al.'s six-input composite"] = () => (PolynomialReference.OfLinearForms("Rumsey et al.'s six-input composite",
                                                             [new Monomial(1.0, [2, 0]), new Monomial(1.0, [1, 1]), new Monomial(1.0 / 9.0, [0, 3])],
                                                             [T1, T2], [0.0, 0.0]),
                                                         ReferenceBox.Cube(6, 0.0, 1.0)),
        ["Rumsey et al.'s three-way interaction"] = () => (Interaction(3), ReferenceBox.Cube(3, 0.0, 1.0)),
        ["Legendre expansion"] = () =>
        {
            ReferenceBox box = Box([0.0, -Math.PI, 2.0], [1.0, Math.PI, 5.0]);
            return (new LegendreExpansion("Legendre expansion", box,
                    [
                        new LegendreTerm(0.8, [1, 0, 0]), new LegendreTerm(-1.1, [0, 2, 0]), new LegendreTerm(0.5, [2, 1, 0]),
                        new LegendreTerm(0.3, [1, 1, 3]), new LegendreTerm(-0.7, [0, 3, 1]), new LegendreTerm(0.25, [4, 0, 2]),
                        new LegendreTerm(0.4, [3, 0, 0]),
                    ]),
                    box);
        },
    };

    private static PolynomialReference Benchmark(int p) => new("Rumsey et al.'s benchmark",
        [new Monomial(1.0, Exponents(p, (0, 2))), new Monomial(1.0, Exponents(p, (0, 1), (1, 1))), new Monomial(1.0 / 9.0, Exponents(p, (1, 3)))]);

    private static PolynomialReference Interaction(int p) => new("Rumsey et al.'s three-way interaction",
        [new Monomial(1.0, Exponents(p, (0, 2))), new Monomial(1.0, Exponents(p, (0, 1), (1, 1), (2, 1))), new Monomial(1.0 / 9.0, Exponents(p, (1, 3)))]);

    private static int[] Exponents(int p, params (int Input, int Power)[] powers)
    {
        var e = new int[p];
        foreach ((int input, int power) in powers) e[input] = power;
        return e;
    }

    public static TheoryData<string> ChosenNames()
    {
        var names = new TheoryData<string>();
        foreach (string name in Chosen.Keys) names.Add(name);
        return names;
    }

    [Theory]
    [MemberData(nameof(ChosenNames))]
    public void ClosedForm_MatchesQuadrature(string name)
    {
        (ReferenceFunction f, ReferenceMeasure law) = Chosen[name]();
        double[,] closed = f.C(law) ?? throw new InvalidOperationException($"{name} has no closed form under its law.");
        double error = RelativeError(closed, Quadrature(f, law));
        output.WriteLine($"{name}: max |C − quadrature| / max |C| = {error:E2}");
        Assert.True(error <= Tolerance, $"{name}: {error:E2}");
    }

    [Theory]
    [MemberData(nameof(ChosenNames))]
    public void Gradient_MatchesCentralDifferences(string name)
    {
        (ReferenceFunction f, ReferenceMeasure law) = Chosen[name]();
        var rng = new Xoshiro256PlusPlus(41);
        int p = f.Dimension;
        var g = new double[p];
        double worst = 0.0;
        for (int point = 0; point < 16; point++)
        {
            double[] x = Draw(rng, law);
            f.Gradient(x, g);
            double scale = g.Max(Math.Abs);
            for (int j = 0; j < p; j++)
            {
                double step = Math.Pow(2.0, Math.Round(Math.Log2(1e-3 * Spread(law, j))));
                double At(double offset)
                {
                    var y = (double[])x.Clone();
                    y[j] += offset;
                    return f.Value(y);
                }
                double difference = (-At(2 * step) + 8 * At(step) - 8 * At(-step) + At(-2 * step)) / (12 * step);
                worst = Math.Max(worst, Math.Abs(difference - g[j]) / scale);
            }
        }
        output.WriteLine($"{name}: worst |∂f − difference| / max |∇f| = {worst:E2}");
        Assert.True(worst <= 1e-8, $"{name}: {worst:E2}");
    }

    /// <summary>½(x − o)ᵀA(x − o) has gradient Ay, so C = A·E[yyᵀ]·A with E[yyᵀ] = ddᵀ + diag(h²/3), d the box's centre less the origin.</summary>
    [Fact]
    public void QuadraticForm_IsAEyyA()
    {
        (ReferenceFunction f, ReferenceMeasure law) = Chosen["quadratic form"]();
        var box = (ReferenceBox)law;
        double[] origin = [0.2, -0.1, 1.0];
        int p = 3;
        var second = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int k = 0; k < p; k++)
                second[i, k] = (box.Center[i] - origin[i]) * (box.Center[k] - origin[k]) + (i == k ? box.HalfWidths[i] * box.HalfWidths[i] / 3.0 : 0.0);
        double[,] expected = Product(Product(Symmetric, second), Symmetric);
        double error = RelativeError(f.C(law)!, expected);
        output.WriteLine($"C against A(ddᵀ + H²/3)A: {error:E2}");
        Assert.True(error <= 1e-14, $"{error:E2}");
    }

    /// <summary>
    /// Rumsey et al. (arXiv 2307.11241, §4.1) state the benchmark's C on [0, 1]² as [[120, 50], [50, 21]]/45, with first active
    /// direction ≈ (0.923, 0.385), and with inert inputs that block padded with zeros. The supplement's three-way interaction
    /// f₃, whose C the paper estimates by Monte Carlo, is likewise its three-input block padded; §4.1.2's composite
    /// f(t₁ᵀx, t₂ᵀx) has C of rank two on span(t₁, t₂). The paper's eigenvalues for the composite are a fit's, so the exact
    /// ones are printed beside them, not asserted.
    /// </summary>
    [Fact]
    public void RumseyBenchmarks_MatchThePaper()
    {
        double[,] stated = { { 120.0 / 45.0, 50.0 / 45.0 }, { 50.0 / 45.0, 21.0 / 45.0 } };
        double error = RelativeError(Benchmark(2).C(ReferenceBox.Cube(2, 0.0, 1.0))!, stated);
        (double top, double[] direction) = TopOfTwo(stated);
        double padding = Math.Max(RelativeError(Benchmark(8).C(ReferenceBox.Cube(8, 0.0, 1.0))!, Padded(Benchmark(2).C(ReferenceBox.Cube(2, 0.0, 1.0))!, 8)),
                                  RelativeError(Interaction(10).C(ReferenceBox.Cube(10, 0.0, 1.0))!, Padded(Interaction(3).C(ReferenceBox.Cube(3, 0.0, 1.0))!, 10)));

        (ReferenceFunction composite, ReferenceMeasure law) = Chosen["Rumsey et al.'s six-input composite"]();
        double[,] c = composite.C(law)!;
        double span = SpanResidual(c, [T1, T2]);
        double[][] q = Orthonormal([T1, T2]);
        var restricted = new double[2, 2];
        for (int a = 0; a < 2; a++)
            for (int b = 0; b < 2; b++)
                for (int i = 0; i < 6; i++)
                    for (int k = 0; k < 6; k++) restricted[a, b] += q[a][i] * c[i, k] * q[b][k];
        double first = TopOfTwo(restricted).Top, second = restricted[0, 0] + restricted[1, 1] - first;

        output.WriteLine($"benchmark against the stated C: {error:E2}; first direction ({direction[0]:F4}, {direction[1]:F4}), eigenvalue {top:F6}");
        output.WriteLine($"inert inputs padded with zeros: {padding:E2}");
        double[,] f3 = Interaction(3).C(ReferenceBox.Cube(3, 0.0, 1.0))!;
        output.WriteLine($"f₃'s C on [0, 1]³ (exact): [[{f3[0, 0]:F6}, {f3[0, 1]:F6}, {f3[0, 2]:F6}], [·, {f3[1, 1]:F6}, {f3[1, 2]:F6}], [·, ·, {f3[2, 2]:F6}]]");
        output.WriteLine($"composite: range in span(t₁, t₂) to {span:E2}; √(λ₂/λ₁) = {Math.Sqrt(second / first):F4} exact (the paper's fit: 0.125)");
        Assert.True(error <= 1e-15, $"{error:E2}");
        Assert.True(Math.Abs(direction[0] - 0.923) <= 5e-4 && Math.Abs(direction[1] - 0.385) <= 5e-4, $"({direction[0]}, {direction[1]})");
        Assert.True(padding == 0.0, $"{padding:E2}");
        Assert.True(span <= 1e-14, $"{span:E2}");
    }

    /// <summary>
    /// The cubic ridge's C, read from its monomials, is rank one on its direction, and equals the quadrature of the same ridge
    /// held as a profile and a direction: two representations that share no code.
    /// </summary>
    [Fact]
    public void PolynomialRidge_IsRankOneOnItsDirection()
    {
        (ReferenceFunction f, ReferenceMeasure law) = Chosen["cubic ridge"]();
        double[] w = [0.7, -0.2, 1.1];
        double[,] closed = f.C(law)!;
        double rankOne = SpanResidual(closed, [w]);
        var profile = new FunctionProfile("1 − 2t + t³", t => 1.0 - 2.0 * t + t * t * t, t => -2.0 + 3.0 * t * t);
        var asRidge = new RidgeSum("cubic ridge as a ridge", new RidgeTerm(profile, w, 0.4));
        double error = RelativeError(closed, GradientQuadrature.Uniform(asRidge, (ReferenceBox)law, 8, [1, 1, 1]));
        output.WriteLine($"rank-one residual {rankOne:E2}; against the ridge's quadrature {error:E2}");
        Assert.True(rankOne <= 1e-14, $"{rankOne:E2}");
        Assert.True(error <= Tolerance, $"{error:E2}");
    }

    /// <summary>A ridge with no closed form (tanh) still has its subspace in closed form: its quadrature C is rank one on w.</summary>
    [Fact]
    public void GeneralRidge_QuadratureIsRankOneOnItsDirection()
    {
        double[] w = [0.3, -0.8, 0.5];
        var f = new RidgeSum("tanh ridge", new RidgeTerm(new FunctionProfile("tanh 3t", t => Math.Tanh(3.0 * t), t => 3.0 / (Math.Cosh(3.0 * t) * Math.Cosh(3.0 * t))), w, -0.1));
        ReferenceBox box = Box([-1.0, 0.0, -0.5], [1.0, 1.0, 1.5]);
        Assert.Null(f.C(box));
        double residual = SpanResidual(GradientQuadrature.Uniform(f, box, 16, [2, 2, 2]), [w]);
        output.WriteLine($"rank-one residual {residual:E2}");
        Assert.True(residual <= 1e-14, $"{residual:E2}");
    }

    /// <summary>∫₀¹ψ′_mψ_n and ∫₀¹ψ′_mψ′_n against Gauss–Legendre, exact on their degrees, for m, n ≤ 12.</summary>
    [Fact]
    public void LegendreSlopeIntegrals_MatchQuadrature()
    {
        QuadratureRule rule = QuadratureRule.Uniform(0.0, 1.0, 16);
        double worst = 0.0;
        for (int m = 0; m <= 12; m++)
            for (int n = 0; n <= 12; n++)
            {
                double product = 0.0, gram = 0.0;
                for (int i = 0; i < rule.Nodes.Length; i++)
                {
                    var (vm, sm) = LegendreExpansion.Psi(m, rule.Nodes[i]);
                    var (vn, sn) = LegendreExpansion.Psi(n, rule.Nodes[i]);
                    product += rule.Weights[i] * sm * vn;
                    gram += rule.Weights[i] * sm * sn;
                }
                double scale = Math.Max(1.0, Math.Abs(LegendreExpansion.SlopeGram(Math.Max(m, n), Math.Max(m, n))));
                worst = Math.Max(worst, Math.Abs(product - LegendreExpansion.SlopeProduct(m, n)) / scale);
                worst = Math.Max(worst, Math.Abs(gram - LegendreExpansion.SlopeGram(m, n)) / scale);
            }
        output.WriteLine($"worst relative error {worst:E2}");
        Assert.True(worst <= 1e-13, $"{worst:E2}");
    }

    /// <summary>The closed forms are refused where the family has none: a polynomial under a Gaussian, a Legendre expansion
    /// off its own box (a general profile's refusal is <see cref="GeneralRidge_QuadratureIsRankOneOnItsDirection"/>'s).</summary>
    [Fact]
    public void ClosedForm_IsRefusedOffItsLaws()
    {
        (ReferenceFunction polynomial, _) = Chosen["interacting polynomial"]();
        (ReferenceFunction legendre, _) = Chosen["Legendre expansion"]();
        Assert.Null(polynomial.C(new ReferenceGaussian([0.0, 0.0, 0.0], Correlated)));
        Assert.Null(legendre.C(ReferenceBox.Cube(3, 0.0, 1.0)));
    }

    /// <summary>
    /// The seeded sweep: 100 random cases per family (ridge sums under a box, ridge sums under a Gaussian, polynomials, Legendre
    /// expansions), one to three inputs (two under the Gaussian), laws centred up to 10⁴ from the origin with widths from 10⁻³ to
    /// 10², rates scaled to the law's spread. Each closed form meets its quadrature within 10⁻¹² plus its <see cref="Conditioning"/>;
    /// where that conditioning is below 10⁻¹³ the worst error is printed apart, the flat gate's reading.
    /// </summary>
    [Fact]
    public void ClosedForms_SeededSweep()
    {
        var rng = new Xoshiro256PlusPlus(20261108);
        string[] families = ["ridges under a box", "ridges under a Gaussian", "polynomials", "Legendre expansions"];
        var worstError = new double[4];
        var worstRatio = new double[4];
        var worstConditioned = new double[4];
        var conditioned = new int[4];
        var failures = new List<string>();
        for (int i = 0; i < 400; i++)
        {
            int family = i % 4, p = 1 + rng.NextInt(family == 1 ? 2 : 3);   // the Gaussian's panelled rule is wide: two inputs at most
            bool far = rng.NextDouble() < 0.5;
            (ReferenceFunction f, ReferenceMeasure law) = family switch
            {
                0 => RandomRidges(rng, RandomBox(rng, p, far)),
                1 => RandomRidges(rng, RandomGaussian(rng, p, far)),
                2 => RandomPolynomial(rng, RandomBox(rng, p, far)),
                _ => RandomLegendre(rng, RandomBox(rng, p, far)),
            };
            double[,] closed = f.C(law)!;
            if (MaxAbs(closed) == 0.0) continue;
            double error = RelativeError(closed, Quadrature(f, law));
            double conditioning = Conditioning(f, law), tolerance = Tolerance + conditioning;
            worstError[family] = Math.Max(worstError[family], error);
            worstRatio[family] = Math.Max(worstRatio[family], error / tolerance);
            if (conditioning < 1e-13)
            {
                conditioned[family]++;
                worstConditioned[family] = Math.Max(worstConditioned[family], error);
            }
            if (!(error <= tolerance)) failures.Add($"case {i} ({families[family]}, p = {p}, far = {far}): {error:E2} against {tolerance:E2}");
        }
        for (int family = 0; family < 4; family++)
            output.WriteLine($"{families[family]}: worst relative error {worstError[family]:E2}, worst error/tolerance {worstRatio[family]:F3}; " +
                             $"{conditioned[family]} cases conditioned below 1e-13, worst {worstConditioned[family]:E2}");
        Assert.True(failures.Count == 0, string.Join(Environment.NewLine, failures));
    }

    // ── The quadrature each family needs ────────────────────────────────────────────────────────────────────────────

    private static double[,] Quadrature(ReferenceFunction f, ReferenceMeasure law)
    {
        int p = f.Dimension;
        switch (f, law)
        {
            case (RidgeSum ridges, ReferenceBox box):
            {
                // Each panel holds |μ_j|·(its half-width) ≤ 4 for every pair's exponent μ.
                var panels = new int[p];
                for (int j = 0; j < p; j++)
                {
                    double z = 2.0 * ridges.Ridges.Max(r => Rate(r).Magnitude * Math.Abs(r.Direction[j])) * box.HalfWidths[j];
                    panels[j] = Math.Max(1, (int)Math.Ceiling(z / 4.0));
                }
                return GradientQuadrature.Uniform(f, box, p >= 5 ? 14 : 20, panels);
            }
            case (RidgeSum ridges, ReferenceGaussian gaussian):
            {
                // In z = L⁻¹(x − m) a pair's exponent is μ = Lᵀ(λ_k w_k + λ_l w_l): the integrand's peak sits at Re μ, and each
                // panel of half-width 1/perUnit holds |μ|/perUnit ≤ 4.
                double real = 0.0, whole = 0.0;
                foreach (RidgeTerm r in ridges.Ridges)
                    for (int j = 0; j < p; j++)
                    {
                        double lw = 0.0;
                        for (int i = j; i < p; i++) lw += gaussian.Factor[i, j] * r.Direction[i];
                        real = Math.Max(real, Math.Abs(Rate(r).Real * lw));
                        whole = Math.Max(whole, Rate(r).Magnitude * Math.Abs(lw));
                    }
                double radius = 9.0 + 2.0 * real;
                int perUnit = Math.Max(1, (int)Math.Ceiling(2.0 * whole / 4.0));
                return GradientQuadrature.Gaussian(f, gaussian, radius, 16, (int)Math.Ceiling(radius) * perUnit);
            }
            case (PolynomialReference polynomial, ReferenceBox box):
                return GradientQuadrature.Uniform(f, box, polynomial.Terms.Max(t => t.Exponents.Max()) + 2, [.. Enumerable.Repeat(1, p)]);
            case (LegendreExpansion legendre, ReferenceBox box):
                return GradientQuadrature.Uniform(f, box, legendre.Terms.Max(t => t.Degrees.Max()) + 2, [.. Enumerable.Repeat(1, p)]);
            default:
                throw new ArgumentException($"No quadrature for {f.Name} under {law.GetType().Name}.");
        }
    }

    private static Complex Rate(RidgeTerm r) => ((ExponentialProfile)r.Profile).Rate;

    /// <summary>
    /// What rounding costs the quadrature, relative to C: it evaluates each argument at the law's points in raw units, where
    /// a ridge's s = wᵀx − c carries an error up to (p + 2)·ε·(Σ|w_j||x_j| + |c|), which moves g′ by |λ| times that and C by at most
    /// four times; a polynomial's y = x − o and a Legendre expansion's u = (x − a)/(b − a) carry ε|x| against a spread of h, moving a
    /// degree-n factor's slope by n (polynomial) or n² (Legendre) times that.
    /// </summary>
    private static double Conditioning(ReferenceFunction f, ReferenceMeasure law)
    {
        int p = f.Dimension;
        double Reach(int j) => law switch
        {
            ReferenceBox box => Math.Max(Math.Abs(box.Lower[j]), Math.Abs(box.Upper[j])),
            ReferenceGaussian g => Math.Abs(g.Center[j]) + 15.0 * Math.Sqrt(g.Covariance[j, j]),
            _ => throw new ArgumentException("An unknown law."),
        };
        switch (f)
        {
            case RidgeSum ridges:
                return ridges.Ridges.Max(r =>
                {
                    double reach = Math.Abs(r.Offset);
                    for (int j = 0; j < p; j++) reach += Math.Abs(r.Direction[j]) * Reach(j);
                    return 4.0 * Rate(r).Magnitude * (p + 2) * Eps * reach;
                });
            case PolynomialReference polynomial:
            {
                var box = (ReferenceBox)law;
                int degree = polynomial.Terms.Max(t => t.Exponents.Max());
                return 4.0 * degree * Eps * Enumerable.Range(0, p).Max(j => Reach(j) / box.HalfWidths[j]);
            }
            case LegendreExpansion legendre:
            {
                var box = (ReferenceBox)law;
                int degree = legendre.Terms.Max(t => t.Degrees.Max());
                return 4.0 * degree * degree * Eps * Enumerable.Range(0, p).Max(j => Reach(j) / (2.0 * box.HalfWidths[j]));
            }
            default:
                throw new ArgumentException($"No conditioning for {f.Name}.");
        }
    }

    // ── Random cases ────────────────────────────────────────────────────────────────────────────────────────────────

    private static ReferenceBox RandomBox(Xoshiro256PlusPlus rng, int p, bool far)
    {
        var lower = new double[p];
        var upper = new double[p];
        for (int j = 0; j < p; j++)
        {
            double center = far ? Sign(rng) * LogUniform(rng, 1.0, 4.0) : -2.0 + 4.0 * rng.NextDouble();
            double half = 0.5 * LogUniform(rng, -3.0, 2.0);
            (lower[j], upper[j]) = (center - half, center + half);
        }
        return new ReferenceBox(lower, upper);
    }

    private static ReferenceGaussian RandomGaussian(Xoshiro256PlusPlus rng, int p, bool far)
    {
        var mean = new double[p];
        var sd = new double[p];
        for (int j = 0; j < p; j++)
        {
            mean[j] = far ? Sign(rng) * LogUniform(rng, 1.0, 4.0) : -2.0 + 4.0 * rng.NextDouble();
            sd[j] = LogUniform(rng, -3.0, 2.0);
        }
        var covariance = new double[p, p];
        double rho = -0.9 + 1.8 * rng.NextDouble();
        for (int i = 0; i < p; i++)
            for (int k = 0; k < p; k++)
                covariance[i, k] = sd[i] * sd[k] * (i == k ? 1.0 : rho);
        return new ReferenceGaussian(mean, covariance);
    }

    // One to three ridges, random unit directions, rates scaled to the spread of wᵀx under the law, offsets near its centre.
    private static (ReferenceFunction, ReferenceMeasure) RandomRidges(Xoshiro256PlusPlus rng, ReferenceMeasure law)
    {
        int p = law.Dimension, count = 1 + rng.NextInt(3);
        var ridges = new RidgeTerm[count];
        for (int k = 0; k < count; k++)
        {
            double[] w = UnitVector(rng, p);
            double spread = law switch
            {
                ReferenceBox box => w.Select((v, j) => Math.Abs(v) * box.HalfWidths[j]).Sum(),
                ReferenceGaussian g => Math.Sqrt(Enumerable.Range(0, p).Sum(i => Enumerable.Range(0, p).Sum(j => w[i] * g.Covariance[i, j] * w[j]))),
                _ => throw new ArgumentException("An unknown law."),
            };
            double realTop = law is ReferenceBox ? 5.0 : 2.5, waveTop = law is ReferenceBox ? 8.0 : 5.0;
            Complex shape = rng.NextInt(3) switch
            {
                0 => new Complex(Sign(rng) * (0.2 + (realTop - 0.2) * rng.NextDouble()), 0.0),
                1 => new Complex(0.0, 0.2 + (waveTop - 0.2) * rng.NextDouble()),
                _ => new Complex(Sign(rng) * (0.1 + 0.4 * realTop * rng.NextDouble()), 0.2 + (waveTop - 0.2) * rng.NextDouble()),
            };
            Complex rate = shape / spread;
            double center = w.Select((v, j) => v * law.Center[j]).Sum();
            ridges[k] = new RidgeTerm(new ExponentialProfile(new Complex(Normal(rng), Normal(rng)), rate), w, center + (2.0 * rng.NextDouble() - 1.0) * spread);
        }
        return (new RidgeSum("random ridges", ridges), law);
    }

    private static (ReferenceFunction, ReferenceMeasure) RandomPolynomial(Xoshiro256PlusPlus rng, ReferenceBox box)
    {
        int p = box.Dimension, count = 1 + rng.NextInt(5);
        var terms = new Monomial[count];
        for (int t = 0; t < count; t++)
        {
            int[] e = [.. Enumerable.Range(0, p).Select(_ => rng.NextInt(4))];
            if (e.All(v => v == 0)) e[rng.NextInt(p)] = 1;
            terms[t] = new Monomial(Normal(rng), e);
        }
        double[] origin = [.. Enumerable.Range(0, p).Select(j => box.Center[j] + (2.0 * rng.NextDouble() - 1.0) * box.HalfWidths[j])];
        return (new PolynomialReference("random polynomial", terms, origin), box);
    }

    private static (ReferenceFunction, ReferenceMeasure) RandomLegendre(Xoshiro256PlusPlus rng, ReferenceBox box)
    {
        int p = box.Dimension, count = 1 + rng.NextInt(6);
        var terms = new LegendreTerm[count];
        for (int t = 0; t < count; t++)
        {
            int[] degrees = [.. Enumerable.Range(0, p).Select(_ => rng.NextInt(6))];
            if (degrees.All(v => v == 0)) degrees[rng.NextInt(p)] = 1 + rng.NextInt(5);
            terms[t] = new LegendreTerm(Normal(rng), degrees);
        }
        return (new LegendreExpansion("random Legendre expansion", box, terms), box);
    }

    // ── Helpers ─────────────────────────────────────────────────────────────────────────────────────────────────────

    private static double[] Draw(Xoshiro256PlusPlus rng, ReferenceMeasure law)
    {
        int p = law.Dimension;
        switch (law)
        {
            case ReferenceBox box:
                return [.. Enumerable.Range(0, p).Select(j => box.Lower[j] + (box.Upper[j] - box.Lower[j]) * rng.NextDouble())];
            case ReferenceGaussian g:
            {
                double[] z = [.. Enumerable.Range(0, p).Select(_ => Normal(rng))];
                return [.. Enumerable.Range(0, p).Select(i => g.Center[i] + Enumerable.Range(0, i + 1).Sum(j => g.Factor[i, j] * z[j]))];
            }
            default:
                throw new ArgumentException("An unknown law.");
        }
    }

    private static double Spread(ReferenceMeasure law, int j) => law switch
    {
        ReferenceBox box => 2.0 * box.HalfWidths[j],
        ReferenceGaussian g => Math.Sqrt(g.Covariance[j, j]),
        _ => throw new ArgumentException("An unknown law."),
    };

    // ‖C − PCP‖_max / ‖C‖_max, P the projector onto the directions' span.
    private static double SpanResidual(double[,] c, double[][] directions)
    {
        int p = c.GetLength(0);
        double[][] q = Orthonormal(directions);
        var projector = new double[p, p];
        foreach (double[] v in q)
            for (int i = 0; i < p; i++)
                for (int k = 0; k < p; k++) projector[i, k] += v[i] * v[k];
        double[,] projected = Product(Product(projector, c), projector);
        return RelativeError(c, projected);
    }

    // Gram–Schmidt, twice.
    private static double[][] Orthonormal(double[][] directions)
    {
        var q = new List<double[]>();
        foreach (double[] d in directions)
        {
            var v = (double[])d.Clone();
            for (int pass = 0; pass < 2; pass++)
                foreach (double[] u in q)
                {
                    double dot = v.Zip(u, (a, b) => a * b).Sum();
                    for (int i = 0; i < v.Length; i++) v[i] -= dot * u[i];
                }
            double norm = Math.Sqrt(v.Sum(a => a * a));
            q.Add([.. v.Select(a => a / norm)]);
        }
        return [.. q];
    }

    // The larger eigenvalue of a symmetric 2 × 2 and its unit eigenvector.
    private static (double Top, double[] Direction) TopOfTwo(double[,] m)
    {
        double a = m[0, 0], b = m[0, 1], c = m[1, 1];
        double top = 0.5 * (a + c) + Math.Sqrt(0.25 * (a - c) * (a - c) + b * b);
        double norm = Math.Sqrt(b * b + (top - a) * (top - a));
        return (top, [b / norm, (top - a) / norm]);
    }

    private static double[,] Padded(double[,] block, int p)
    {
        var padded = new double[p, p];
        for (int i = 0; i < block.GetLength(0); i++)
            for (int k = 0; k < block.GetLength(1); k++) padded[i, k] = block[i, k];
        return padded;
    }

    private static double[,] Product(double[,] a, double[,] b)
    {
        int n = a.GetLength(0), m = b.GetLength(1), inner = a.GetLength(1);
        var c = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < m; k++)
                for (int j = 0; j < inner; j++) c[i, k] += a[i, j] * b[j, k];
        return c;
    }

    private static double RelativeError(double[,] closed, double[,] other)
    {
        double worst = 0.0;
        for (int i = 0; i < closed.GetLength(0); i++)
            for (int k = 0; k < closed.GetLength(1); k++) worst = Math.Max(worst, Math.Abs(closed[i, k] - other[i, k]));
        return worst / MaxAbs(closed);
    }

    private static double MaxAbs(double[,] c)
    {
        double m = 0.0;
        foreach (double v in c) m = Math.Max(m, Math.Abs(v));
        return m;
    }

    private static double[] UnitVector(Xoshiro256PlusPlus rng, int p)
    {
        double[] w = [.. Enumerable.Range(0, p).Select(_ => Normal(rng))];
        double norm = Math.Sqrt(w.Sum(v => v * v));
        return [.. w.Select(v => v / norm)];
    }

    private static double Normal(Xoshiro256PlusPlus rng) =>
        Math.Sqrt(-2.0 * Math.Log(1.0 - rng.NextDouble())) * Math.Cos(2.0 * Math.PI * rng.NextDouble());

    private static double LogUniform(Xoshiro256PlusPlus rng, double lowerExponent, double upperExponent) =>
        Math.Pow(10.0, lowerExponent + (upperExponent - lowerExponent) * rng.NextDouble());

    private static double Sign(Xoshiro256PlusPlus rng) => rng.NextDouble() < 0.5 ? -1.0 : 1.0;
}
