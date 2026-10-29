using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's moment service (O5, <see cref="BassMoments"/>) against references that share none of its integration:
/// <list type="bullet">
/// <item>term means and cross moments at degrees 1–3 against a tensor quadrature of the basis's own columns
/// (<see cref="BassBasis.CaseColumn"/>) under a product of a truncated normal–Student mixture (weights conditional on
/// truncation, clipped inside the training range), a uniform on an inner interval, and level probabilities with a zero
/// level;</item>
/// <item>the empirical marginal against the average of the columns over the training cases, whose values are the knots,
/// so every atom sits on a breakpoint;</item>
/// <item>categorical moments against enumeration, the zero-probability level included;</item>
/// <item>each pair's covariance, the sum of its ANOVA components, against Cross − Product, and its totals against the
/// components that hold each variable.</item>
/// </list>
/// </summary>
public sealed class BassSensitivityMomentTests
{
    private static readonly BassData Data = SensitivityFixture.Scalar();

    private static readonly Dictionary<string, InputMeasure> Measures = new()
    {
        ["x0"] = new MixtureMeasure([new NormalComponent(0.4, 0.3, 0.15), new StudentComponent(0.6, 2.0, 0.7, 0.2)])
        {
            Lower = 0.05, Upper = 0.95, Weighting = MixtureWeighting.ConditionalOnTruncation,
        },
        ["x1"] = new UniformMeasure(0.2, 0.9),
        ["z"] = new CategoricalMeasure(new Dictionary<string, double> { ["a"] = 0.0, ["b"] = 0.7, ["c"] = 0.3 }),
    };

    private static BassTerm[] Terms() =>
    [
        new([new HingeFactor(0, 1, 3)]),
        new([new HingeFactor(0, -1, 6), new HingeFactor(1, 1, 2)]),
        new([new HingeFactor(1, -1, 8)], [new CategoricalFactor(0, [0, 1])]),
        new(categorical: [new CategoricalFactor(0, [0])]),
        new([new HingeFactor(0, 1, 3), new HingeFactor(1, -1, 5)], [new CategoricalFactor(0, [1, 2])]),
        new([new HingeFactor(0, -1, 3)]),
    ];

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void TermMoments_MatchTensorQuadratureOfTheBasisColumns(int degree)
    {
        BassInputs inputs = BassInputs.Prepare(Data, BassInputTransform.Fit(Data));
        var measure = ProductMeasure.Resolve(BassInputTransform.Fit(Data), inputs, Measures);
        var moments = BassMoments.On(inputs, degree, measure);
        var reference = new SensitivityReference(inputs, measure, degree);
        BassTerm[] terms = Terms();
        foreach (BassTerm a in terms)
        {
            double mean = reference.Expectation(a);
            Assert.True(Close(moments.Mean(a), mean), $"E[{a}] at q = {degree}: {moments.Mean(a)} against {mean}");
            foreach (BassTerm b in terms)
            {
                double cross = reference.Expectation(a, b);
                Assert.True(Close(moments.Cross(a, b), cross), $"E[{a} · {b}] at q = {degree}: {moments.Cross(a, b)} against {cross}");
            }
        }
        Assert.InRange(moments.LargestError, 0.0, 1e-12);   // the Student component's Gauss–Legendre estimate
        Assert.Equal(0.0, moments.Mean(terms[3]));         // the zero-probability level
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void EmpiricalMarginal_IsTheAverageOverTheTrainingValues(int degree)
    {
        var transform = BassInputTransform.Fit(Data);
        BassInputs inputs = BassInputs.Prepare(Data, transform);
        var measure = ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure> { ["x0"] = new EmpiricalMeasure() });
        var moments = BassMoments.On(inputs, degree, measure);
        var basis = new BassBasis(inputs, degree);
        int knots = inputs.ContinuousKnots(0).Length;
        var terms = new List<BassTerm>();
        for (int k = 0; k < knots; k++)
            foreach (int s in new[] { 1, -1 }) terms.Add(new BassTerm([new HingeFactor(0, s, k)]));
        double[] Column(BassTerm t)
        {
            var c = new double[inputs.Cases.Cases];
            basis.CaseColumn(t, inputs.Cases, c);
            return c;
        }
        foreach (BassTerm a in terms)
        {
            double[] ca = Column(a);
            Assert.Equal(ca.Average(), moments.Mean(a), 14);
            foreach (BassTerm b in terms)
            {
                double[] cb = Column(b);
                Assert.Equal(ca.Zip(cb, (u, v) => u * v).Average(), moments.Cross(a, b), 14);
            }
        }
    }

    [Fact]
    public void CategoricalMoments_AreSumsOfTheProbabilities_IntersectionsIncluded()
    {
        var transform = BassInputTransform.Fit(Data);
        BassInputs inputs = BassInputs.Prepare(Data, transform);
        var measure = ProductMeasure.Resolve(transform, inputs, Measures);
        var moments = BassMoments.On(inputs, 1, measure);
        double[] p = [0.0, 0.7, 0.3];
        int[][] subsets = [[0], [1], [0, 1], [1, 2], [0, 2]];
        foreach (int[] d1 in subsets)
        {
            var t1 = new BassTerm(categorical: [new CategoricalFactor(0, d1)]);
            Assert.Equal(d1.Sum(l => p[l]), moments.Mean(t1), 15);
            foreach (int[] d2 in subsets)
            {
                var t2 = new BassTerm(categorical: [new CategoricalFactor(0, d2)]);
                double both = d1.Intersect(d2).Sum(l => p[l]);
                TermPairMoments pair = moments.Pair(t1, t2);
                Assert.Equal(both, pair.Cross, 15);
                Assert.Equal(both - d1.Sum(l => p[l]) * d2.Sum(l => p[l]), pair.Covariance, 15);
                Assert.False(double.IsNaN(pair.Covariance));
            }
        }
        Assert.Equal(new[] { 0.0, 0.7, 0.3 }, measure.Probabilities(2)!);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(3)]
    public void PairStructure_SumsToTheCovariance_AndItsTotalsToTheirComponents(int degree)
    {
        var transform = BassInputTransform.Fit(Data);
        BassInputs inputs = BassInputs.Prepare(Data, transform);
        var moments = BassMoments.On(inputs, degree, ProductMeasure.Resolve(transform, inputs, Measures));
        BassTerm[] terms = Terms();
        foreach (BassTerm a in terms)
            foreach (BassTerm b in terms)
            {
                TermPairMoments pair = moments.Pair(a, b);
                double scale = Math.Abs(pair.Cross) + Math.Abs(pair.Product);
                Assert.True(Math.Abs(pair.Covariance - (pair.Cross - pair.Product)) <= 1e-13 * Math.Max(scale, 1e-300),
                    $"{a} · {b}: {pair.Covariance} against {pair.Cross - pair.Product}");
                int[] shared = [.. moments.Variables(a).Intersect(moments.Variables(b)).Order()];
                Assert.Equal((1 << shared.Length) - 1, pair.Components.Count);
                Assert.Equal(shared, pair.Totals.Select(t => t.Variable));
                foreach (VariableCovariance total in pair.Totals)
                {
                    double sum = pair.Components.Where(c => c.Variables.Contains(total.Variable)).Sum(c => c.Value);
                    Assert.True(Math.Abs(total.Value - sum) <= 1e-13 * Math.Max(scale, 1e-300), $"{a} · {b}, variable {total.Variable}");
                }
                Assert.Same(pair, moments.Pair(b, a));
            }
    }

    [Fact]
    public void Resolution_RefusesUnknownNamesAndWrongKinds()
    {
        var transform = BassInputTransform.Fit(Data);
        BassInputs inputs = BassInputs.Prepare(Data, transform);
        Assert.Throws<ArgumentException>(() => ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure> { ["w"] = new UniformMeasure() }));
        Assert.Throws<ArgumentException>(() => ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure> { ["z"] = new UniformMeasure() }));
        Assert.Throws<ArgumentException>(() => ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure> { ["x0"] = new CategoricalMeasure() }));
        Assert.Throws<ArgumentException>(() => ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure>
        {
            ["z"] = new CategoricalMeasure(new Dictionary<string, double> { ["d"] = 1.0 }),
        }));
        Assert.Throws<ArgumentException>(() => ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure>
        {
            ["x0"] = new MixtureMeasure([new NormalComponent(-1.0, 0.5, 0.1)]),
        }));
        Assert.Throws<ArgumentException>(() => ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure>
        {
            ["x0"] = new UniformMeasure(0.6, 0.4),
        }));
        // A component far outside the interval still has a truncated law (its log mass is finite), concentrated at the near end.
        var far = ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure>
        {
            ["x0"] = new MixtureMeasure([new NormalComponent(1.0, 5.0, 0.01)]) { Weighting = MixtureWeighting.ConditionalOnTruncation },
        });
        Assert.Equal(1.0, BassMoments.On(inputs, 1, far).Mean(new BassTerm([new HingeFactor(0, 1, 0)])), 4);
        var measure = ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure> { ["x1"] = new UniformMeasure(-0.5, 1.0) });
        Assert.True(measure.Continuous(1)!.Extrapolates);
        Assert.False(measure.Continuous(0)!.Extrapolates);
        Assert.IsType<UniformMeasure>(measure.Declared[0]);
        // A marginal sd of 0.2 at ν = 4 is a scale of 0.2·√(1/2).
        Assert.Equal(0.2 * Math.Sqrt(0.5), StudentComponent.FromSd(1.0, 4.0, 0.5, 0.2).Scale, 15);
        Assert.Throws<ArgumentOutOfRangeException>(() => StudentComponent.FromSd(1.0, 2.0, 0.5, 0.2));
    }

    private static bool Close(double actual, double expected)
        => Math.Abs(actual - expected) <= 1e-12 * Math.Max(Math.Abs(expected), 1e-3);
}

/// <summary>
/// The sensitivity tests' data and fits, and a brute-force reference: expectations under a <see cref="ProductMeasure"/> by a
/// tensor rule over each variable (composite Gauss–Legendre on the density, split at every knot; atoms for the empirical
/// measure; levels for a categorical one), evaluating terms with the basis's own columns.
/// </summary>
internal static class SensitivityFixture
{
    // x0, x1 on grids with ends 0 and 1, so the transform is the identity; z with levels a, b, c.
    public static BassData Scalar(int seed = 7)
    {
        var rng = new Maths.Rng.Xoshiro256PlusPlus(seed);
        const int n = 60;
        var x0 = new double[n];
        var x1 = new double[n];
        var z = new int[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = (i % 10) / 9.0;
            x1[i] = (i * 7 % 11) / 10.0;
            z[i] = i % 3;
            y[i] = 1 + 2 * Math.Max(x0[i] - 0.4, 0) + 1.5 * x1[i] * (z[i] == 1 ? 1 : 0) + 3 * x0[i] * x1[i] + 0.1 * (rng.NextDouble() - 0.5);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1), new CategoricalPredictor("z", z, ["a", "b", "c"])],
                            new ScalarResponse(y));
    }

    // A functional response over t on [0, 1] at 7 locations, x0 and x1 as above.
    public static BassData Functional(int seed = 9)
    {
        var rng = new Maths.Rng.Xoshiro256PlusPlus(seed);
        const int n = 40, m = 7;
        double[] t = [.. Enumerable.Range(0, m).Select(r => r / (m - 1.0))];
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            x0[i] = (i % 8) / 7.0;
            x1[i] = (i * 3 % 5) / 4.0;
            for (int r = 0; r < m; r++)
                y[i, r] = 1 + 2 * x0[i] * t[r] + Math.Max(x1[i] - 0.3, 0) + 0.1 * (rng.NextDouble() - 0.5);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)],
                            new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
    }

    public static BassConfig Config(int degree = 1, int maxBasis = 8)
        => new()
        {
            Structure = new BassStructureSpec { MaxInteraction = 3, MaxBasis = maxBasis, Degree = degree },
            Sampling = new BassSampling { Chains = 2, MasterSeed = 11, BurnIn = 300, Samples = 30, Thin = 2 },
        };

    private static readonly Lazy<BassFit> ScalarFit = new(() => Bass.Fit(Config(), Scalar()));
    private static readonly Lazy<BassFit> ScalarFit2 = new(() => Bass.Fit(Config(degree: 2), Scalar()));
    private static readonly Lazy<BassFit> FunctionalFit = new(() => Bass.Fit(Config(), Functional()));

    public static BassFit FitScalar(int degree = 1) => degree == 1 ? ScalarFit.Value : ScalarFit2.Value;

    public static BassFit FitFunctional() => FunctionalFit.Value;
}

/// <summary>Expectations of terms (and of any function of the case inputs) under a product measure by a tensor rule.</summary>
internal sealed class SensitivityReference
{
    private readonly BassInputs _inputs;
    private readonly ProductMeasure _measure;
    private readonly BassBasis _basis;

    public SensitivityReference(BassInputs inputs, ProductMeasure measure, int degree)
    {
        _inputs = inputs;
        _measure = measure;
        _basis = new BassBasis(inputs, degree);
    }

    // One variable's nodes and weights: atoms, levels, or composite Gauss–Legendre on the density over the support, split
    // at the variable's knots and into panels of width at most 0.05.
    public (double[] X, double[] W) Nodes(int variable)
    {
        if (_measure.Probabilities(variable) is { } p) return ([.. Enumerable.Range(0, p.Count).Select(l => (double)l)], [.. p]);
        ContinuousMeasure c = _measure.Continuous(variable)!;
        if (c.Measure is Empirical e) return (e.Values.ToArray(), [.. Enumerable.Repeat(1.0 / e.Count, e.Count)]);
        var density = (IContinuousDistribution)c.Measure;
        (bool _, int index) = _inputs.PredictorPosition(variable);
        var cuts = new SortedSet<double> { c.Lower, c.Upper };
        foreach (double k in _inputs.ContinuousKnots(index))
            if (k > c.Lower && k < c.Upper) cuts.Add(k);
        double[] b = [.. cuts];
        var x = new List<double>();
        var w = new List<double>();
        GaussLegendre rule = GaussLegendre.Of(20);
        for (int i = 0; i + 1 < b.Length; i++)
        {
            int panels = (int)Math.Ceiling((b[i + 1] - b[i]) / 0.05);
            double width = (b[i + 1] - b[i]) / panels;
            for (int q = 0; q < panels; q++)
            {
                double mid = b[i] + (q + 0.5) * width;
                for (int r = 0; r < rule.Count; r++)
                {
                    double at = mid + 0.5 * width * rule.Nodes[r];
                    x.Add(at);
                    w.Add(0.5 * width * rule.Weights[r] * Math.Exp(density.LogPdf(at)));
                }
            }
        }
        return ([.. x], [.. w]);
    }

    // E[Πₖ fₖ] over the variables the terms hold, the others fixed (the terms do not read them).
    public double Expectation(params BassTerm[] terms)
    {
        int[] vars = [.. Enumerable.Range(0, _measure.Count).Where(v => terms.Any(t => Holds(t, v)))];
        var nodes = vars.Select(Nodes).ToArray();
        int total = nodes.Aggregate(1, (acc, n) => acc * n.X.Length);
        var continuous = new double[_inputs.Cases.ContinuousCount][];
        var categorical = new int[_inputs.Cases.CategoricalCount][];
        for (int k = 0; k < continuous.Length; k++) continuous[k] = Enumerable.Repeat(0.5, total).ToArray();
        for (int k = 0; k < categorical.Length; k++) categorical[k] = new int[total];
        var weights = Enumerable.Repeat(1.0, total).ToArray();
        for (int i = 0; i < total; i++)
        {
            int rest = i;
            for (int k = vars.Length - 1; k >= 0; k--)
            {
                int at = rest % nodes[k].X.Length;
                rest /= nodes[k].X.Length;
                weights[i] *= nodes[k].W[at];
                (bool isContinuous, int index) = _inputs.PredictorPosition(vars[k]);
                if (isContinuous) continuous[index][i] = nodes[k].X[at];
                else categorical[index][i] = (int)nodes[k].X[at];
            }
        }
        var cases = new CaseInputs(continuous, categorical);
        var product = Enumerable.Repeat(1.0, total).ToArray();
        var column = new double[total];
        foreach (BassTerm t in terms)
        {
            _basis.CaseColumn(t, cases, column);
            for (int i = 0; i < total; i++) product[i] *= column[i];
        }
        double sum = 0.0;
        for (int i = 0; i < total; i++) sum += weights[i] * product[i];
        return sum;
    }

    private bool Holds(BassTerm t, int variable)
    {
        (bool isContinuous, int index) = _inputs.PredictorPosition(variable);
        return isContinuous ? t.Continuous.Any(h => h.Variable == index) : t.Categorical.Any(c => c.Variable == index);
    }
}
