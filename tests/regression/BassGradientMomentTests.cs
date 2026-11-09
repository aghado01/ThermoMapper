using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 1 gates (the active-subspace brief): the moment service's derivative moments (<see cref="DesignFactor.SlopeMean"/>,
/// <see cref="DesignFactor.SlopePair"/>, <see cref="BassMoments.GradientPair"/>) and C per draw on the exact route
/// (<see cref="BassMoments.GradientOuterProduct"/>), against references that share none of their integration:
/// <list type="bullet">
/// <item>each factor pair against one-variable quadrature of the factors' own values and derivatives, hinge degrees one to three and
/// Legendre degrees to ten, under a uniform law, a truncated normal–Student mixture and an empirical marginal whose atoms sit on the
/// knots; and a seeded sweep with knots off the support and mixture components far outside it;</item>
/// <item>each term pair against a tensor quadrature of the basis's own derivative columns (<see cref="BassBasis.CaseDerivative"/>),
/// categorical factors included;</item>
/// <item>C on fitted draws, hinge and Legendre, against the gradient readout (<see cref="BassFit.Gradient"/>) averaged over the
/// measure by a rule exact for the integrand, under the uniform law, a mixture and an empirical marginal, where the readout's
/// mean of one-sided limits at a knot is what the exact route must reproduce;</item>
/// <item>SG-4(b)'s identity E[x(1 − x)ψ′_mψ′_n] = n(n + 1)·δ_mn through the slope moments under the Jacobi weight, and on a
/// Legendre fit's draws against the coefficient readout's <see cref="PolynomialChaosSensitivity.Jacobi"/>.</item>
/// </list>
/// </summary>
public sealed class BassGradientMomentTests(ITestOutputHelper output)
{
    // ── One-variable laws and their quadrature ───────────────────────────────────────────────────────────────────────

    private static ContinuousMeasure UniformLaw(double lo, double hi) => new(Uniform.On(lo, hi), lo, hi);

    private static ContinuousMeasure MixtureLaw(double lo, double hi, params IContinuousDistribution[] parts)
        => new(Truncated.Of(Mixture.Of(parts, [.. Enumerable.Repeat(1.0, parts.Length)]), lo, hi), lo, hi);

    private static ContinuousMeasure EmpiricalLaw(double[] values)
    {
        Empirical e = Empirical.Of(values);
        return new(e, e.Lower, e.Upper);
    }

    // One variable's rule under its law: atoms as they are; a density by n-point Gauss–Legendre on panels split at the breaks and
    // no wider than width.
    internal static (double[] X, double[] W) Rule(ContinuousMeasure m, IEnumerable<double> breaks, int n, double width)
    {
        if (m.Measure is Empirical e) return (e.Values.ToArray(), [.. Enumerable.Repeat(1.0 / e.Count, e.Count)]);
        var density = (IContinuousDistribution)m.Measure;
        var cuts = new SortedSet<double> { m.Lower, m.Upper };
        foreach (double b in breaks)
            if (b > m.Lower && b < m.Upper) cuts.Add(b);
        double[] c = [.. cuts];
        GaussLegendre rule = GaussLegendre.Of(n);
        var x = new List<double>();
        var w = new List<double>();
        for (int i = 0; i + 1 < c.Length; i++)
        {
            int panels = Math.Max(1, (int)Math.Ceiling((c[i + 1] - c[i]) / width));
            double step = (c[i + 1] - c[i]) / panels;
            for (int q = 0; q < panels; q++)
            {
                double mid = c[i] + (q + 0.5) * step;
                for (int r = 0; r < rule.Count; r++)
                {
                    double at = mid + 0.5 * step * rule.Nodes[r];
                    x.Add(at);
                    w.Add(0.5 * step * rule.Weights[r] * Math.Exp(density.LogPdf(at)));
                }
            }
        }
        return ([.. x], [.. w]);
    }

    // ∫g dμ and ∫|g| dμ by the rule.
    private static (double Value, double Size) Integrate((double[] X, double[] W) rule, Func<double, double> g)
    {
        double value = 0.0, size = 0.0;
        for (int i = 0; i < rule.X.Length; i++)
        {
            double v = rule.W[i] * g(rule.X[i]);
            value += v;
            size += Math.Abs(v);
        }
        return (value, size);
    }

    private static readonly double[] Atoms = [0.0, 0.1, 0.1, 0.25, 0.3, 0.45, 0.6, 0.6, 0.6, 0.72, 0.8, 0.95, 1.0];

    // ── Factor pairs ─────────────────────────────────────────────────────────────────────────────────────────────────

    public static TheoryData<string> Laws() => new() { "uniform", "mixture", "empirical" };

    private static ContinuousMeasure Law(string name) => name switch
    {
        "uniform" => UniformLaw(0.1, 0.85),
        "mixture" => MixtureLaw(0.05, 0.95, Normal.MeanSd(0.3, 0.15), StudentT.LocationScale(2.0, 0.7, 0.2)),
        "empirical" => EmpiricalLaw(Atoms),
        _ => throw new ArgumentException(name),
    };

    /// <summary>E[φ′], E[φ′ψ] and E[φ′ψ′] for hinges of degree one to three, every sign pair, knots inside the support, at its ends
    /// and outside it, and on its atoms; at an atom on a knot the slope is the derivative's mean of one-sided limits.</summary>
    [Theory]
    [MemberData(nameof(Laws))]
    public void HingeSlopeMoments_MatchQuadrature(string law)
    {
        ContinuousMeasure m = Law(law);
        double[] knots = [-0.2, 0.0, 0.1, 0.3, 0.45, 0.6, 0.8, 1.0, 1.3];
        double worst = 0.0;
        string where = "";
        void Track(double gap, string what)
        {
            if (gap > worst) (worst, where) = (gap, what);
        }
        for (int q = 1; q <= 3; q++)
        {
            var factors = new List<HingeDesignFactor>();
            for (int k = 0; k < knots.Length; k++)
                foreach (int s in new[] { 1, -1 }) factors.Add(new HingeDesignFactor(0, s, k, knots[k], q));
            var rule = law == "mixture" ? Rule(m, knots, 20, 0.05) : Rule(m, knots, 2 * q + 2, 1.0);
            foreach (HingeDesignFactor a in factors)
            {
                (double mean, double ms) = Integrate(rule, a.Derivative);
                double got = a.SlopeMean(m, out _);
                Track(Gap(got, mean, ms), $"E[φ′], q {q}, ({a.Sign}, {a.KnotValue}): {got:R} against {mean:R}");
                foreach (HingeDesignFactor b in factors)
                {
                    (double sc, double scs) = Integrate(rule, x => a.Derivative(x) * b.Value(x));
                    (double ss, double sss) = Integrate(rule, x => a.Derivative(x) * b.Derivative(x));
                    (double gotC, double gotS) = a.SlopePair(b, m, out _);
                    Track(Gap(gotC, sc, scs), $"E[φ′ψ], q {q}, ({a.Sign}, {a.KnotValue}), ({b.Sign}, {b.KnotValue}): {gotC:R} against {sc:R}");
                    Track(Gap(gotS, ss, sss), $"E[φ′ψ′], q {q}, ({a.Sign}, {a.KnotValue}), ({b.Sign}, {b.KnotValue}): {gotS:R} against {ss:R}");
                }
            }
        }
        output.WriteLine($"{law}: worst gap over the integrand's size {worst:E2}, at {where}");
        Assert.True(worst <= 1e-12, $"{law}: {worst:E2}");
    }

    /// <summary>
    /// Legendre factors to degree ten: the slope moments against quadrature under each law, and under the uniform law on [0, 1]
    /// against the closed forms ∫ψ′_mψ_n and ∫ψ′_mψ′_n (<see cref="LegendreExpansion"/>). The factor expands ψ_α about each piece's
    /// centre, whose coefficients grow with α (<see cref="LegendreDesignFactor"/>), and a slope multiplies the kth by k, so the digits
    /// lost grow with the degree: the gate is 10⁻¹² to degree six and 10⁻¹⁰ to degree ten, both readings printed.
    /// </summary>
    [Theory]
    [MemberData(nameof(Laws))]
    public void LegendreSlopeMoments_MatchQuadrature(string law)
    {
        ContinuousMeasure m = Law(law);
        var rule = law == "mixture" ? Rule(m, [], 20, 0.05) : Rule(m, [], 12, 1.0);
        double[] worst = new double[2], closed = new double[2];   // degrees to six, and seven to ten
        ContinuousMeasure unit = UniformLaw(0.0, 1.0);
        for (int a = 1; a <= 10; a++)
        {
            var fa = new LegendreDesignFactor(0, a);
            (double mean, double ms) = Integrate(rule, fa.Derivative);
            worst[a <= 6 ? 0 : 1] = Math.Max(worst[a <= 6 ? 0 : 1], Gap(fa.SlopeMean(m, out _), mean, ms));
            for (int b = 1; b <= 10; b++)
            {
                int band = Math.Max(a, b) <= 6 ? 0 : 1;
                var fb = new LegendreDesignFactor(0, b);
                (double sc, double scs) = Integrate(rule, x => fa.Derivative(x) * fb.Value(x));
                (double ss, double sss) = Integrate(rule, x => fa.Derivative(x) * fb.Derivative(x));
                (double gotC, double gotS) = fa.SlopePair(fb, m, out _);
                worst[band] = Math.Max(worst[band], Math.Max(Gap(gotC, sc, scs), Gap(gotS, ss, sss)));
                (double c, double s) = fa.SlopePair(fb, unit, out _);
                double scale = Math.Max(1.0, LegendreExpansion.SlopeGram(Math.Max(a, b), Math.Max(a, b)));
                closed[band] = Math.Max(closed[band], Math.Max(Math.Abs(c - LegendreExpansion.SlopeProduct(a, b)), Math.Abs(s - LegendreExpansion.SlopeGram(a, b))) / scale);
            }
        }
        output.WriteLine($"{law}: worst gap to degree six {worst[0]:E2}, to ten {worst[1]:E2}; on [0, 1] against the closed forms {closed[0]:E2} and {closed[1]:E2}");
        Assert.True(Math.Max(worst[0], closed[0]) <= 1e-12, $"{law} to degree six: {worst[0]:E2}, {closed[0]:E2}");
        Assert.True(Math.Max(worst[1], closed[1]) <= 1e-10, $"{law} to degree ten: {worst[1]:E2}, {closed[1]:E2}");
    }

    /// <summary>
    /// The seeded sweep: 400 random hinge pairs, degree one to three, knots anywhere in [−1, 2], under a uniform law on a random
    /// interval, a truncated mixture whose components sit up to 40 scales outside the interval, or an empirical marginal on which
    /// some knots fall exactly.
    /// </summary>
    [Fact]
    public void SlopeMoments_SeededSweep_FarRangesIncluded()
    {
        var rng = new Xoshiro256PlusPlus(20261109);
        double worst = 0.0;
        int[] counts = new int[3];
        for (int i = 0; i < 400; i++)
        {
            int kind = i % 3, q = 1 + rng.NextInt(3);
            double lo = -0.5 + rng.NextDouble(), hi = lo + 0.05 + 1.5 * rng.NextDouble();
            ContinuousMeasure m;
            double width = 1.0;
            double[] atoms = [];
            switch (kind)
            {
                case 0:
                    m = UniformLaw(lo, hi);
                    break;
                case 1:
                {
                    double scale = 0.02 + 0.5 * rng.NextDouble(), side = rng.NextDouble() < 0.5 ? lo : hi;
                    double far = (rng.NextDouble() < 0.5 ? -1 : 1) * 40.0 * rng.NextDouble() * scale;
                    m = MixtureLaw(lo, hi, Normal.MeanSd(side + far, scale), StudentT.LocationScale(3.0, lo + (hi - lo) * rng.NextDouble(), 0.3 * scale + 0.01));
                    width = Math.Min(0.05, 0.25 * (0.3 * scale + 0.01));
                    break;
                }
                default:
                    atoms = [.. Enumerable.Range(0, 25).Select(_ => Math.Round(lo + (hi - lo) * rng.NextDouble(), 2))];
                    m = EmpiricalLaw(atoms);
                    break;
            }
            double Knot() => atoms.Length > 0 && rng.NextDouble() < 0.5 ? atoms[rng.NextInt(atoms.Length)] : -1.0 + 3.0 * rng.NextDouble();
            var a = new HingeDesignFactor(0, rng.NextDouble() < 0.5 ? 1 : -1, 0, Knot(), q);
            var b = new HingeDesignFactor(0, rng.NextDouble() < 0.5 ? 1 : -1, 1, Knot(), q);
            var rule = Rule(m, [a.KnotValue, b.KnotValue], kind == 1 ? 20 : 2 * q + 2, width);
            (double mean, double ms) = Integrate(rule, a.Derivative);
            (double sc, double scs) = Integrate(rule, x => a.Derivative(x) * b.Value(x));
            (double ss, double sss) = Integrate(rule, x => a.Derivative(x) * b.Derivative(x));
            (double gotC, double gotS) = a.SlopePair(b, m, out double error);
            double gap = Math.Max(Gap(a.SlopeMean(m, out _), mean, ms), Math.Max(Gap(gotC, sc, scs), Gap(gotS, ss, sss)));
            Assert.True(gap <= 1e-11 + error, $"case {i} (kind {kind}, q {q}, knots {a.KnotValue}, {b.KnotValue} on [{lo}, {hi}]): {gap:E2}");
            worst = Math.Max(worst, gap);
            counts[kind]++;
        }
        output.WriteLine($"{counts[0]} uniform, {counts[1]} far mixtures, {counts[2]} empirical: worst gap over the integrand's size {worst:E2}");
    }

    // The gap relative to the integrand's size ∫|g| dμ, with a floor: an integrand near zero is held to an absolute 10⁻¹⁶ at the 10⁻¹²
    // gate, rounding for factors of order one.
    private static double Gap(double got, double expected, double size) => Math.Abs(got - expected) / (size + 1e-4);

    // ── Term pairs ───────────────────────────────────────────────────────────────────────────────────────────────────

    private static readonly Dictionary<string, InputMeasure> FixtureMeasures = new()
    {
        ["x0"] = new MixtureMeasure([new NormalComponent(0.4, 0.3, 0.15), new StudentComponent(0.6, 2.0, 0.7, 0.2)])
        {
            Lower = 0.05, Upper = 0.95, Weighting = MixtureWeighting.ConditionalOnTruncation,
        },
        ["x1"] = new UniformMeasure(0.2, 0.9),
        ["z"] = new CategoricalMeasure(new Dictionary<string, double> { ["a"] = 0.0, ["b"] = 0.7, ["c"] = 0.3 }),
    };

    private static BassTerm[] FixtureTerms() =>
    [
        new([new HingeFactor(0, 1, 3)]),
        new([new HingeFactor(0, -1, 6), new HingeFactor(1, 1, 2)]),
        new([new HingeFactor(1, -1, 8)], [new CategoricalFactor(0, [0, 1])]),
        new(categorical: [new CategoricalFactor(0, [0])]),
        new([new HingeFactor(0, 1, 3), new HingeFactor(1, -1, 5)], [new CategoricalFactor(0, [1, 2])]),
        new([new HingeFactor(0, -1, 3)]),
        new([new HingeFactor(1, 1, 7)], [new CategoricalFactor(0, [1])]),
    ];

    /// <summary>E[∂ᵢB_a·∂ⱼB_b] for every pair of seven terms, hinge products with and without categorical factors, at degrees one to
    /// three, against a tensor quadrature of the basis's derivative columns over x0's mixture, x1's uniform law and z's levels.</summary>
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void GradientPair_MatchesTensorQuadratureOfTheBasisDerivatives(int degree)
    {
        BassData data = SensitivityFixture.Scalar();
        var transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        var measure = ProductMeasure.Resolve(transform, inputs, FixtureMeasures);
        var moments = BassMoments.On(inputs, degree, measure);
        var basis = new BassBasis(inputs, degree);
        var grid = Grid(inputs, measure, 2 * degree + 2, 0.05);
        BassTerm[] terms = FixtureTerms();
        double worst = 0.0;
        foreach (BassTerm a in terms)
            foreach (BassTerm b in terms)
            {
                TermPairGradientMoments g = moments.GradientPair(a, b);
                int[] rows = [.. a.Continuous.Select(h => h.Variable).Distinct().Order()];
                int[] columns = [.. b.Continuous.Select(h => h.Variable).Distinct().Order()];
                Assert.Equal(rows, g.Rows);
                Assert.Equal(columns, g.Columns);
                for (int r = 0; r < rows.Length; r++)
                    for (int s = 0; s < columns.Length; s++)
                    {
                        double[] da = Derivative(basis, a, rows[r], grid.Cases), db = Derivative(basis, b, columns[s], grid.Cases);
                        double value = 0.0, size = 0.0;
                        for (int i = 0; i < grid.W.Length; i++)
                        {
                            value += grid.W[i] * da[i] * db[i];
                            size += Math.Abs(grid.W[i] * da[i] * db[i]);
                        }
                        worst = Math.Max(worst, Math.Abs(g.Values[r, s] - value) / Math.Max(size, 1e-300));
                    }
            }
        output.WriteLine($"degree {degree}: worst gap over the integrand's size {worst:E2}; largest quadrature error {moments.LargestError:E2}");
        Assert.True(worst <= 1e-11, $"{worst:E2}");
    }

    private static double[] Derivative(BassBasis basis, BassTerm term, int variable, CaseInputs cases)
    {
        var column = new double[cases.Cases];
        basis.CaseDerivative(term, variable, cases, column);
        return column;
    }

    // The tensor rule over the fixture's x0, x1 and z: each continuous variable's rule split at its knots, the levels by their
    // probabilities.
    internal static (CaseInputs Cases, double[] W, double[][] X, int[] Z) Grid(BassInputs inputs, ProductMeasure measure, int n, double width)
    {
        var rules = new (double[] X, double[] W)[2];
        for (int k = 0; k < 2; k++)
        {
            ContinuousMeasure m = measure.Continuous(measure.OfContinuous(k))!;
            bool smooth = m.Measure is not Empirical && m.Measure is not Uniform;
            rules[k] = smooth ? Rule(m, inputs.ContinuousKnots(k).ToArray(), Math.Max(n, 20), width) : Rule(m, inputs.ContinuousKnots(k).ToArray(), n, 1.0);
        }
        IReadOnlyList<double> p = measure.Probabilities(measure.OfCategorical(0))!;
        int total = rules[0].X.Length * rules[1].X.Length * p.Count;
        var x0 = new double[total];
        var x1 = new double[total];
        var z = new int[total];
        var w = new double[total];
        int at = 0;
        for (int i = 0; i < rules[0].X.Length; i++)
            for (int j = 0; j < rules[1].X.Length; j++)
                for (int l = 0; l < p.Count; l++, at++)
                {
                    x0[at] = rules[0].X[i];
                    x1[at] = rules[1].X[j];
                    z[at] = l;
                    w[at] = rules[0].W[i] * rules[1].W[j] * p[l];
                }
        return (new CaseInputs([x0, x1], [z]), w, [x0, x1], z);
    }

    // ── C on fitted draws against the readout ────────────────────────────────────────────────────────────────────────

    public static TheoryData<int, string> FittedCases() => new()
    {
        { 1, "uniform" }, { 2, "uniform" }, { 1, "mixture" }, { 2, "mixture" }, { 1, "empirical" }, { 2, "empirical" },
    };

    /// <summary>
    /// Per draw of a hinge fit with a categorical input, C on the exact route against the gradient readout's outer product summed
    /// over a rule exact for the integrand (the knots split every panel; a smooth density takes narrow panels), with z's level
    /// probabilities. The fixture's box is [0, 1], so the readout's raw units are the transformed ones. Under the empirical marginal
    /// every atom is a knot, where the readout takes the mean of one-sided limits.
    /// </summary>
    [Theory]
    [MemberData(nameof(FittedCases))]
    public void ExactC_OnHingeDraws_IsTheReadoutAveragedOverTheMeasure(int degree, string law)
    {
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = SensitivityFixture.FitScalar(degree);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        var measures = new Dictionary<string, InputMeasure> { ["z"] = new CategoricalMeasure(new Dictionary<string, double> { ["a"] = 0.2, ["b"] = 0.5, ["c"] = 0.3 }) };
        if (law == "mixture") measures["x0"] = FixtureMeasures["x0"];
        if (law == "empirical") measures["x0"] = new EmpiricalMeasure();
        var measure = ProductMeasure.Resolve(fit.Transform, inputs, measures);
        var moments = BassMoments.On(inputs, degree, measure);
        var grid = Grid(inputs, measure, degree + 2, 0.05);
        string[] levels = ["a", "b", "c"];
        var cases = new PredictionCases([new ContinuousColumn("x0", grid.X[0]), new ContinuousColumn("x1", grid.X[1]),
                                         new CategoricalColumn("z", [.. grid.Z.Select(l => levels[l])])]);
        BassGradient readout = fit.Gradient(data, cases);
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        double worst = 0.0;
        for (int d = 0; d < draws.Length; d++)
        {
            BassTerm[] terms = [.. draws[d].Structure.TermIds.Select(id => fit.Terms[id])];
            double[,] exact = moments.GradientOuterProduct(terms, draws[d].Coefficients.AsSpan(1));
            worst = Math.Max(worst, Gap(exact, Average(readout, d, grid.W)));
        }
        output.WriteLine($"degree {degree}, {law}: {draws.Length} draws, worst max |ΔC| / max |C| {worst:E2}");
        Assert.True(worst <= 1e-12, $"{worst:E2}");
    }

    private static readonly Lazy<(BassData Data, BassFit Fit)> Legendre = new(() =>
    {
        BassData data = KhaosReadoutTests.Data(120, 20261109);
        return (data, Bass.Fit(KhaosReadoutTests.Config(6), data));
    });

    /// <summary>Per draw of a Legendre fit on three inputs, C on the exact route against the readout averaged by a Gauss–Legendre
    /// rule exact for the polynomials, under the uniform law and with x0's empirical marginal.</summary>
    [Theory]
    [InlineData("uniform")]
    [InlineData("empirical")]
    public void ExactC_OnLegendreDraws_IsTheReadoutAveragedOverTheMeasure(string law)
    {
        (BassData data, BassFit fit) = Legendre.Value;
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        var measures = new Dictionary<string, InputMeasure>();
        if (law == "empirical") measures["x0"] = new EmpiricalMeasure();
        var measure = ProductMeasure.Resolve(fit.Transform, inputs, measures);
        var moments = BassMoments.On(inputs, fit.Config.Structure.Degree, measure);
        var rules = Enumerable.Range(0, 3).Select(k => Rule(measure.Continuous(k)!, [], 8, 1.0)).ToArray();
        int total = rules.Aggregate(1, (acc, r) => acc * r.X.Length);
        var x = new double[3][];
        for (int k = 0; k < 3; k++) x[k] = new double[total];
        var w = new double[total];
        for (int i = 0; i < total; i++)
        {
            int rest = i;
            w[i] = 1.0;
            for (int k = 2; k >= 0; k--)
            {
                int at = rest % rules[k].X.Length;
                rest /= rules[k].X.Length;
                x[k][i] = rules[k].X[at];
                w[i] *= rules[k].W[at];
            }
        }
        BassGradient readout = fit.Gradient(data, new PredictionCases([.. Enumerable.Range(0, 3).Select(k => (PredictionColumn)new ContinuousColumn($"x{k}", x[k]))]));
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        double worst = 0.0;
        for (int d = 0; d < draws.Length; d++)
        {
            BassTerm[] terms = [.. draws[d].Structure.TermIds.Select(id => fit.Terms[id])];
            worst = Math.Max(worst, Gap(moments.GradientOuterProduct(terms, draws[d].Coefficients.AsSpan(1)), Average(readout, d, w)));
        }
        output.WriteLine($"{law}: {draws.Length} draws, worst max |ΔC| / max |C| {worst:E2}");
        Assert.True(worst <= 1e-12, $"{worst:E2}");
    }

    // Σᵢ wᵢ∇f∇fᵀ for draw d of the readout (one location).
    private static double[,] Average(BassGradient g, int d, double[] w)
    {
        int p = g.Partials.Count;
        var c = new double[p, p];
        for (int i = 0; i < w.Length; i++)
            for (int a = 0; a < p; a++)
            {
                double ga = w[i] * g.Partials[a][d, i, 0];
                for (int b = 0; b < p; b++) c[a, b] += ga * g.Partials[b][d, i, 0];
            }
        return c;
    }

    private static double Gap(double[,] exact, double[,] reference)
    {
        double worst = 0.0, scale = 0.0;
        for (int a = 0; a < exact.GetLength(0); a++)
            for (int b = 0; b < exact.GetLength(1); b++)
            {
                worst = Math.Max(worst, Math.Abs(exact[a, b] - reference[a, b]));
                scale = Math.Max(scale, Math.Abs(reference[a, b]));
            }
        return scale == 0.0 ? worst : worst / scale;
    }

    // ── SG-4(b) ──────────────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>The Jacobi weight's law on [0, 1], density 6x(1 − x), with exact interval moments: Beta(2, 2).</summary>
    private sealed class JacobiLaw : IIntervalMoments
    {
        public double IntervalMoments(double lower, double upper, double centre, double logScale, Span<double> moments)
        {
            double a = Math.Max(lower, 0.0), b = Math.Min(upper, 1.0);
            moments.Clear();
            if (!(b > a)) return 0.0;
            // 6x(1 − x) about c, with y = x − c: 6[c(1 − c) + (1 − 2c)y − y²].
            double c0 = 6.0 * centre * (1.0 - centre), c1 = 6.0 * (1.0 - 2.0 * centre), c2 = -6.0;
            double ya = a - centre, yb = b - centre, scale = Math.Exp(-logScale);
            for (int k = 0; k < moments.Length; k++)
            {
                double I(int power) => (Math.Pow(yb, power + 1) - Math.Pow(ya, power + 1)) / (power + 1);
                moments[k] = scale * (c0 * I(k) + c1 * I(k + 1) + c2 * I(k + 2));
            }
            return 0.0;
        }
    }

    /// <summary>
    /// SG-4(b): E_U[x(1 − x)ψ′_mψ′_n] = n(n + 1)·δ_mn, read through the slope moments under the Jacobi law (which carries the weight
    /// as 6x(1 − x)); and, on each draw of a Legendre fit, E[x_v(1 − x_v)(∂_v f)²] assembled from them, the other inputs' degrees
    /// matched by orthonormality, against <see cref="PolynomialChaosSensitivity.Jacobi"/>'s Σ b_α²α_v(α_v + 1).
    /// </summary>
    [Fact]
    public void Sg4b_JacobiIdentity_ThroughTheSlopeMoments()
    {
        var jacobi = new ContinuousMeasure(new JacobiLaw(), 0.0, 1.0);
        double[] factorGap = new double[2];   // degrees to six, and seven to ten
        for (int m = 1; m <= 10; m++)
            for (int n = 1; n <= 10; n++)
            {
                double got = new LegendreDesignFactor(0, m).SlopePair(new LegendreDesignFactor(0, n), jacobi, out _).SlopeSlope / 6.0;
                int band = Math.Max(m, n) <= 6 ? 0 : 1;
                factorGap[band] = Math.Max(factorGap[band], Math.Abs(got - (m == n ? n * (n + 1.0) : 0.0)) / Math.Max(m * (m + 1.0), n * (n + 1.0)));
            }

        (BassData data, BassFit fit) = Legendre.Value;
        double[,] coefficientRoute = BassSensitivity.FromCoefficients(fit, data).Jacobi;
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        double drawGap = 0.0;
        for (int d = 0; d < draws.Length; d++)
        {
            BassTerm[] terms = [.. draws[d].Structure.TermIds.Select(id => fit.Terms[id])];
            ReadOnlySpan<double> b = draws[d].Coefficients.AsSpan(1);
            for (int v = 0; v < 3; v++)
            {
                double sum = 0.0;
                for (int s = 0; s < terms.Length; s++)
                    for (int t = 0; t < terms.Length; t++)
                    {
                        int ds = DegreeOn(terms[s], v), dt = DegreeOn(terms[t], v);
                        if (ds == 0 || dt == 0 || !OthersMatch(terms[s], terms[t], v)) continue;
                        sum += b[s] * b[t] * new LegendreDesignFactor(v, ds).SlopePair(new LegendreDesignFactor(v, dt), jacobi, out _).SlopeSlope / 6.0;
                    }
                drawGap = Math.Max(drawGap, Math.Abs(sum - coefficientRoute[d, v]) / Math.Max(1.0, Math.Abs(coefficientRoute[d, v])));
            }
        }
        output.WriteLine($"factor pairs to degree six {factorGap[0]:E2}, to ten {factorGap[1]:E2}; {draws.Length} draws × 3 inputs against the coefficient readout: {drawGap:E2}");
        Assert.True(factorGap[0] <= 1e-12, $"to degree six: {factorGap[0]:E2}");
        Assert.True(factorGap[1] <= 1e-10, $"to degree ten: {factorGap[1]:E2}");
        Assert.True(drawGap <= 1e-12, $"{drawGap:E2}");
    }

    private static int DegreeOn(BassTerm term, int v)
    {
        foreach (LegendreFactor f in term.Legendre)
            if (f.Variable == v) return f.Degree;
        return 0;
    }

    private static bool OthersMatch(BassTerm s, BassTerm t, int v)
    {
        for (int u = 0; u < 3; u++)
            if (u != v && DegreeOn(s, u) != DegreeOn(t, u)) return false;
        return true;
    }
}
