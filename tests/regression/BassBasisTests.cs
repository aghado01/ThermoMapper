using System;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's terms and basis (stage 4, S2): terms are canonical and compare by value, and the table interns them; the
/// inputs' knots are the distinct scaled values; the case and output factors match R's construction by hand, its
/// normalization and its zero-constant guard included; support counts each part; and a direct functional design's
/// separable form matches an explicit nm × k design evaluated independently, Gram included. R parity on random terms
/// is the oracle test host's (<c>BassBasisParityTests</c>).
/// </summary>
public sealed class BassBasisTests
{
    private static BassData Grid(double[] x, int[]? codes = null, string[]? levels = null, double[,]? y = null, double[]? t = null)
    {
        var predictors = new System.Collections.Generic.List<BassPredictor> { new ContinuousPredictor("x", x) };
        if (codes is not null) predictors.Add(new CategoricalPredictor("z", codes, levels!));
        BassResponse response = y is null
            ? new ScalarResponse(x.Select((v, i) => Math.Sin(v) + i).ToArray())
            : new FunctionalResponse(y, [new ContinuousPredictor("t", t!)]);
        return new BassData(predictors, response);
    }

    private static BassBasis Basis(BassData data, int degree) => new(BassInputs.Prepare(data, BassInputTransform.Fit(data)), degree);

    private static double[] Case(BassBasis basis, BassTerm term)
    {
        basis.Validate(term);
        var column = new double[basis.Inputs.Cases.Cases];
        basis.CaseColumn(term, basis.Inputs.Cases, column);
        return column;
    }

    [Fact]
    public void Terms_AreCanonical_AndCompareByValue()
    {
        var a = new BassTerm([new HingeFactor(2, 1, 4), new HingeFactor(0, -1, 1)], [new CategoricalFactor(1, [3, 0])]);
        var b = new BassTerm([new HingeFactor(0, -1, 1), new HingeFactor(2, 1, 4)], [new CategoricalFactor(1, [0, 3])]);
        Assert.Equal(a, b);
        Assert.Equal(a.GetHashCode(), b.GetHashCode());
        Assert.Equal([0, 2], a.Continuous.Select(h => h.Variable));
        Assert.NotEqual(a, new BassTerm([new HingeFactor(0, -1, 1), new HingeFactor(2, -1, 4)], [new CategoricalFactor(1, [0, 3])]));
        Assert.NotEqual(new BassTerm([new HingeFactor(0, 1, 1)]), new BassTerm(functional: [new HingeFactor(0, 1, 1)]));

        Assert.Throws<ArgumentException>(() => new BassTerm());
        Assert.Throws<ArgumentException>(() => new BassTerm([new HingeFactor(0, 1, 1), new HingeFactor(0, -1, 2)]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new BassTerm([new HingeFactor(0, 0, 1)]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new BassTerm([new HingeFactor(0, 1, -1)]));
        Assert.Throws<ArgumentException>(() => new BassTerm(categorical: [new CategoricalFactor(0, [1]), new CategoricalFactor(0, [2])]));
        Assert.Throws<ArgumentException>(() => new CategoricalFactor(0, []));
        Assert.Throws<ArgumentException>(() => new CategoricalFactor(0, [1, 1]));
    }

    [Fact]
    public void Table_InternsEachDistinctTerm_Once_InOrderFirstSeen()
    {
        var table = new BassTermTable();
        var u = new BassTerm([new HingeFactor(0, 1, 2)]);
        var v = new BassTerm([new HingeFactor(0, -1, 2)]);
        Assert.Equal(0, table.Intern(u));
        Assert.Equal(1, table.Intern(v));
        Assert.Equal(0, table.Intern(new BassTerm([new HingeFactor(0, 1, 2)])));
        Assert.Equal(2, table.Count);
        Assert.Same(v, table[1]);
        Assert.True(table.TryGetId(v, out int id) && id == 1);
        Assert.False(table.TryGetId(new BassTerm([new HingeFactor(1, 1, 2)]), out _));
    }

    [Fact]
    public void Inputs_ScaleColumns_AndTakeDistinctScaledValues_AsKnots()
    {
        BassData data = Grid([3, 1, 5, 1, 3, 9], codes: [2, 0, 2, 2, 0, 0], levels: ["a", "unused", "c"]);
        BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
        Assert.Equal([0.25, 0.0, 0.5, 0.0, 0.25, 1.0], inputs.Cases.Continuous(0).ToArray());
        Assert.Equal([0.0, 0.25, 0.5, 1.0], inputs.ContinuousKnots(0).ToArray());
        Assert.Equal([1, 0, 1, 1, 0, 0], inputs.Cases.Categorical(0).ToArray());   // 'unused' leaves the training levels
        Assert.Equal(2, inputs.LevelCount(0));
        Assert.Equal((true, 0), inputs.PredictorPosition(0));
        Assert.Equal((false, 0), inputs.PredictorPosition(1));
        Assert.Null(inputs.Outputs);
    }

    [Fact]
    public void CaseFactors_MatchRsConstruction_ByHand()
    {
        BassBasis linear = Basis(Grid([0, 1, 2, 3, 4]), 1);   // scaled 0, .25, .5, .75, 1; knots the same
        Assert.Equal([0, 0, 0, 0.5, 1], Case(linear, new BassTerm([new HingeFactor(0, 1, 2)])));   // (x − .5)₊ / .5
        Assert.Equal([1, 0.5, 0, 0, 0], Case(linear, new BassTerm([new HingeFactor(0, -1, 2)])));  // (.5 − x)₊ / .5
        Assert.Equal([0, 0.25, 0.5, 0.75, 1], Case(linear, new BassTerm([new HingeFactor(0, 1, 0)])));    // (x − 0)₊ / 1, x itself

        BassBasis quadratic = Basis(Grid([0, 1, 2, 3, 4]), 2);
        Assert.Equal([1, 0.25, 0, 0, 0], Case(quadratic, new BassTerm([new HingeFactor(0, -1, 2)])));
        Assert.Equal([0, 0, 0, 0.25, 1], Case(quadratic, new BassTerm([new HingeFactor(0, 1, 2)])));

        // A hinge facing the end its knot sits at is zero on the training range, and R leaves its constant at 1.
        Assert.All(Case(linear, new BassTerm([new HingeFactor(0, 1, 4)])), v => Assert.Equal(0.0, v));
        Assert.All(Case(linear, new BassTerm([new HingeFactor(0, -1, 0)])), v => Assert.Equal(0.0, v));
    }

    [Fact]
    public void Interactions_MultiplyFactors_AndPeakAtOne()
    {
        var rng = new Xoshiro256PlusPlus(seed: 11);
        double[] x1 = Enumerable.Range(0, 40).Select(_ => rng.NextDouble()).ToArray();
        double[] x2 = Enumerable.Range(0, 40).Select(_ => rng.NextDouble()).Append(0.0).Append(1.0).ToArray();
        x1 = [.. x1, 1.0, 0.0];   // a corner each hinge pair peaks at
        int[] codes = Enumerable.Range(0, 42).Select(i => i % 3).ToArray();
        var data = new BassData(
            [new ContinuousPredictor("x1", x1), new ContinuousPredictor("x2", x2), new CategoricalPredictor("z", codes, ["a", "b", "c"])],
            new ScalarResponse(x1.Zip(x2, (a, b) => a * b).ToArray()));
        BassBasis basis = Basis(data, 1);
        ReadOnlySpan<double> s1 = basis.Inputs.Cases.Continuous(0), s2 = basis.Inputs.Cases.Continuous(1);
        double k1 = basis.Inputs.ContinuousKnots(0)[10], k2 = basis.Inputs.ContinuousKnots(1)[30];

        var term = new BassTerm([new HingeFactor(1, -1, 30), new HingeFactor(0, 1, 10)], [new CategoricalFactor(0, [0, 2])]);
        double[] column = Case(basis, term);
        double c = (1 - k1) * k2;
        for (int i = 0; i < column.Length; i++)
        {
            double expected = Math.Max(s1[i] - k1, 0) * Math.Max(k2 - s2[i], 0) / c * (codes[i] % 3 == 1 ? 0 : 1);
            Assert.Equal(expected, column[i], 15);
        }
        // The continuous part alone reaches one at the corner (x1 = 1, x2 = 0) and stays below it elsewhere.
        double[] continuousPart = Case(basis, new BassTerm([new HingeFactor(0, 1, 10), new HingeFactor(1, -1, 30)]));
        Assert.Equal(1.0, continuousPart.Max(), 15);
    }

    [Fact]
    public void CategoricalFactors_AreIndicators_AndAnUnseenLevelIsOutsideEverySubset()
    {
        BassData data = Grid([0, 1, 2, 3], codes: [0, 1, 2, 1], levels: ["a", "b", "c"]);
        BassBasis basis = Basis(data, 1);
        var term = new BassTerm(categorical: [new CategoricalFactor(0, [1, 2])]);
        Assert.Equal([0, 1, 1, 1], Case(basis, term));

        var newCases = new CaseInputs([[0.5, 0.5]], [[-1, 2]]);
        var column = new double[2];
        basis.CaseColumn(term, newCases, column);
        Assert.Equal([0, 1], column);

        Assert.Throws<ArgumentException>(() => basis.Validate(new BassTerm(categorical: [new CategoricalFactor(0, [0, 1, 2])])));
        Assert.Throws<ArgumentException>(() => basis.Validate(new BassTerm(categorical: [new CategoricalFactor(0, [3])])));
        Assert.Throws<ArgumentException>(() => basis.Validate(new BassTerm([new HingeFactor(0, 1, 4)])));
        Assert.Throws<ArgumentException>(() => basis.Validate(new BassTerm([new HingeFactor(1, 1, 0)])));
        Assert.Throws<ArgumentException>(() => basis.Validate(new BassTerm(functional: [new HingeFactor(0, 1, 0)])));
    }

    [Fact]
    public void Support_CountsEachPart()
    {
        var y = new double[5, 4];
        for (int i = 0; i < 5; i++) for (int r = 0; r < 4; r++) y[i, r] = i * r;
        BassData data = Grid([0, 1, 2, 3, 4], y: y, t: [0, 1, 2, 3]);
        BassBasis basis = Basis(data, 1);
        var term = new BassTerm([new HingeFactor(0, 1, 1)], functional: [new HingeFactor(0, -1, 2)]);
        TermSupport s = basis.Support(term);
        Assert.Equal((3, (0.25 + 0.5 + 0.75) / 0.75), (s.CaseNonzero, s.CaseMass));          // (x − .25)₊ / .75 at x = .5, .75, 1
        Assert.Equal((2, (2.0 / 3 + 1.0 / 3) / (2.0 / 3)), (s.OutputNonzero, s.OutputMass)); // (2/3 − t)₊ / (2/3) at t = 0, 1/3
        TermSupport functionalOnly = basis.Support(new BassTerm(functional: [new HingeFactor(0, 1, 0)]));
        Assert.Equal((5, 5.0), (functionalOnly.CaseNonzero, functionalOnly.CaseMass));
    }

    // f_m(x, t) straight from the definition, every factor at once: an evaluator independent of the separable split.
    private static double Direct(BassTerm term, BassInputs inputs, int i, int r, int q)
    {
        double value = 1.0, caseScale = 1.0, outputScale = 1.0;
        foreach (HingeFactor h in term.Continuous)
        {
            double t = inputs.ContinuousKnots(h.Variable)[h.Knot];
            value *= Math.Pow(Math.Max(h.Sign * (inputs.Cases.Continuous(h.Variable)[i] - t), 0.0), q);
            caseScale *= h.Sign > 0 ? 1 - t : t;
        }
        foreach (HingeFactor h in term.Functional)
        {
            double t = inputs.FunctionalKnots(h.Variable)[h.Knot];
            value *= Math.Pow(Math.Max(h.Sign * (inputs.Outputs!.Variable(h.Variable)[r] - t), 0.0), q);
            outputScale *= h.Sign > 0 ? 1 - t : t;
        }
        foreach (CategoricalFactor c in term.Categorical)
            if (!c.Levels.Contains(inputs.Cases.Categorical(c.Variable)[i])) value = 0.0;
        if (term.Continuous.Count > 0 && caseScale != 0) value /= Math.Pow(caseScale, q);
        if (term.Functional.Count > 0 && outputScale != 0) value /= Math.Pow(outputScale, q);
        return value;
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void FunctionalDesign_SeparableForm_MatchesTheExplicitDesign(int degree)
    {
        const int n = 6, m = 5;
        var rng = new Xoshiro256PlusPlus(seed: 23 + degree);
        double[] x = Enumerable.Range(0, n).Select(_ => rng.NextDouble()).ToArray();
        int[] codes = [0, 1, 2, 0, 1, 2];
        double[] t = [0.0, 0.1, 0.4, 0.7, 1.0];
        var y = new double[n, m];
        for (int i = 0; i < n; i++) for (int r = 0; r < m; r++) y[i, r] = rng.NextDouble();
        var data = new BassData([new ContinuousPredictor("x", x), new CategoricalPredictor("z", codes, ["a", "b", "c"])],
                                new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
        BassBasis basis = Basis(data, degree);
        BassTerm[] terms =
        [
            new([new HingeFactor(0, 1, 2)]),
            new(functional: [new HingeFactor(0, -1, 3)]),
            new([new HingeFactor(0, -1, 4)], [new CategoricalFactor(0, [1])], [new HingeFactor(0, 1, 1)]),
            new(categorical: [new CategoricalFactor(0, [0, 2])], functional: [new HingeFactor(0, 1, 2)]),
        ];
        var a = new double[terms.Length][];
        var h = new double[terms.Length][];
        for (int j = 0; j < terms.Length; j++)
        {
            basis.Validate(terms[j]);
            a[j] = new double[n];
            h[j] = new double[m];
            basis.CaseColumn(terms[j], basis.Inputs.Cases, a[j]);
            basis.OutputColumn(terms[j], basis.Inputs.Outputs!, h[j]);
        }

        var z = new double[n * m, terms.Length];
        for (int j = 0; j < terms.Length; j++)
            for (int i = 0; i < n; i++)
                for (int r = 0; r < m; r++)
                {
                    z[i * m + r, j] = Direct(terms[j], basis.Inputs, i, r, degree);
                    Assert.Equal(z[i * m + r, j], a[j][i] * h[j][r], 14);
                }

        double[,] gram = KhatriRao.Gram(a, h);
        for (int j = 0; j < terms.Length; j++)
            for (int k = 0; k < terms.Length; k++)
            {
                double explicitEntry = 0.0;
                for (int row = 0; row < n * m; row++) explicitEntry += z[row, j] * z[row, k];
                Assert.Equal(explicitEntry, gram[j, k], 12);
            }
    }
}
