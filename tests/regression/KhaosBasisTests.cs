using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// KHAOS's basis on the engine (the KHAOS brief's slice 1): the degree-split term space's size against counting and against
/// the brief's sum, at small sizes and far ranges; its uniform draw and the uniform composition; its refusals; the Legendre
/// part's column, derivatives and moments through the design-kernel contract; and the structure prior's uniform term measure.
/// </summary>
public sealed class KhaosBasisTests
{
    // Inputs on exactly [0, 1], so the range scaling is the identity and the columns read the values as given.
    internal static BassData UnitData(int n, int p, int seed, Func<double[], double>? f = null, double noise = 0.0)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[p][];
        for (int j = 0; j < p; j++)
        {
            x[j] = new double[n];
            for (int i = 0; i < n; i++) x[j][i] = rng.NextDouble();
            double min = x[j].Min(), max = x[j].Max();   // the ends land on 0 and 1 exactly
            for (int i = 0; i < n; i++) x[j][i] = (x[j][i] - min) / (max - min);
        }
        var y = new double[n];
        var point = new double[p];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < p; j++) point[j] = x[j][i];
            y[i] = (f?.Invoke(point) ?? point.Sum()) + noise * Normal.SampleStandard(rng);
        }
        return new BassData([.. Enumerable.Range(0, p).Select(j => new ContinuousPredictor($"x{j}", x[j]))], new ScalarResponse(y));
    }

    private static long Choose(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        long c = 1;
        for (int i = 1; i <= k; i++) c = c * (n - k + i) / i;
        return c;
    }

    // Every multi-index α ∈ {0, …, d}ᵖ with 1 ≤ |supp α| ≤ q and Σα ≤ d, counted directly: no composition, no subset.
    private static long CountMultiIndices(int p, int d, int q)
    {
        long count = 0;
        var alpha = new int[p];
        void Walk(int j, int sum, int support)
        {
            if (j == p)
            {
                if (support >= 1 && support <= q) count++;
                return;
            }
            for (int a = 0; a + sum <= d; a++)
            {
                alpha[j] = a;
                Walk(j + 1, sum + a, support + (a > 0 ? 1 : 0));
            }
        }
        Walk(0, 0, 0);
        return count;
    }

    [Theory]
    [InlineData(1, 1, 1)]
    [InlineData(3, 3, 2)]
    [InlineData(4, 3, 2)]
    [InlineData(5, 6, 3)]
    [InlineData(6, 4, 4)]
    [InlineData(4, 9, 4)]
    public void TermSpaceSize_MatchesCountingAndTheBriefsSum(int p, int d, int q)
    {
        var space = new LegendreTermSpace(p, d, q);
        List<BassTerm> terms = [.. space.Terms()];
        long brief = 0;
        for (int i = 1; i <= q; i++)
            for (int j = 1; j <= d; j++) brief += Choose(p, i) * Choose(j - 1, i - 1);
        Assert.Equal(brief, (long)space.Size);
        Assert.Equal(CountMultiIndices(p, d, q), (long)space.Size);
        Assert.Equal(terms.Count, (long)space.Size);
        Assert.Equal(terms.Count, terms.Distinct().Count());
        foreach (BassTerm t in terms) space.Validate(t);
        Assert.Equal(Math.Log(space.Size), space.LogSize, 1e-14);
    }

    [Fact]
    public void TermSpaceSize_HoldsAtFarRanges()
    {
        // The cost report's sizes, exact in longs; and past a long, the log size against the doubles' sum.
        var cost = new LegendreTermSpace(20, 15, 5);
        long brief = 0;
        for (int i = 1; i <= 5; i++) brief += Choose(20, i) * Choose(15, i);
        Assert.Equal(brief, (long)cost.Size);
        var wide = new LegendreTermSpace(400, 40, 12);
        double top = double.NegativeInfinity;
        var logs = new double[12];
        for (int i = 1; i <= 12; i++) top = Math.Max(top, logs[i - 1] = LegendreTermSpace.LogChoose(400, i) + LegendreTermSpace.LogChoose(40, i));
        double expected = top + Math.Log(logs.Sum(l => Math.Exp(l - top)));
        Assert.Equal(expected, wide.LogSize, 1e-12 * expected);
        Assert.Equal(-Math.Log(10.0), LegendreTermSpace.LogComposition(6, 3), 1e-14);
        Assert.Equal(0.0, LegendreTermSpace.LogComposition(7, 1));
        Assert.Equal(0.0, LegendreTermSpace.LogComposition(4, 4));
    }

    // Pearson's statistic against its degrees of freedom, a six-sigma band.
    private static void Uniform<T>(IReadOnlyCollection<T> support, IEnumerable<T> draws, int count, string what) where T : notnull
    {
        var tally = support.ToDictionary(s => s, _ => 0);
        foreach (T d in draws)
        {
            Assert.True(tally.ContainsKey(d), $"{what}: a draw outside the support");
            tally[d]++;
        }
        double expected = (double)count / support.Count, chi = tally.Values.Sum(c => (c - expected) * (c - expected) / expected);
        int df = support.Count - 1;
        Assert.True(chi < df + 6 * Math.Sqrt(2.0 * df), $"{what}: χ² {chi:F1} on {df} degrees of freedom");
    }

    [Theory]
    [InlineData(3, 3, 2)]
    [InlineData(4, 3, 2)]
    [InlineData(5, 4, 3)]
    public void TermSpaceDraw_IsUniform(int p, int d, int q)
    {
        var space = new LegendreTermSpace(p, d, q);
        List<BassTerm> terms = [.. space.Terms()];
        var rng = new Xoshiro256PlusPlus(20261107 + p);
        int draws = 2000 * terms.Count;
        Uniform(terms, Enumerable.Range(0, draws).Select(_ => space.Draw(rng)), draws, $"𝒜({p}, {d}, {q})");
    }

    [Fact]
    public void Composition_IsUniform_OverItsCuts()
    {
        var rng = new Xoshiro256PlusPlus(5);
        foreach ((int d, int q) in new[] { (6, 3), (5, 2), (8, 4) })
        {
            List<string> all = [];
            void Walk(string prefix, int left, int parts)
            {
                if (parts == 1) { all.Add(prefix + left); return; }
                for (int a = 1; a <= left - parts + 1; a++) Walk(prefix + a + ",", left - a, parts - 1);
            }
            Walk("", d, q);
            Assert.Equal(Math.Exp(-LegendreTermSpace.LogComposition(d, q)), all.Count, 1e-9);
            int draws = 4000 * all.Count;
            Uniform(all, Enumerable.Range(0, draws).Select(_ => string.Join(",", LegendreTermSpace.DrawComposition(rng, d, q))), draws, $"compositions of {d} into {q}");
        }
        Assert.Equal([5], LegendreTermSpace.DrawComposition(rng, 5, 1));
        Assert.Equal([1, 1, 1], LegendreTermSpace.DrawComposition(rng, 3, 3));
    }

    [Fact]
    public void TermSpace_RefusesTermsOutsideIt()
    {
        var space = new LegendreTermSpace(3, 4, 2);
        space.Validate(new BassTerm(legendre: [new(0, 2), new(2, 2)]));
        Assert.Throws<ArgumentException>(() => space.Validate(new BassTerm(legendre: [new(0, 1), new(1, 1), new(2, 1)])));   // order 3
        Assert.Throws<ArgumentException>(() => space.Validate(new BassTerm(legendre: [new(0, 3), new(1, 2)])));              // degree 5
        Assert.Throws<ArgumentException>(() => space.Validate(new BassTerm(legendre: [new(3, 1)])));                         // no variable 3
        Assert.Throws<ArgumentException>(() => space.Validate(new BassTerm([new HingeFactor(0, 1, 0)])));                   // a hinge
        Assert.Throws<ArgumentException>(() => new BassTerm(legendre: [new(0, 1), new(0, 2)]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new BassTerm(legendre: [new(0, 0)]));
        Assert.Throws<ArgumentOutOfRangeException>(() => new LegendreTermSpace(3, 1, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => new LegendreTermSpace(2, 5, 3));
    }

    [Fact]
    public void Terms_CompareByValue_AndAHingeTermStaysAsItWas()
    {
        var a = new BassTerm(legendre: [new(2, 1), new(0, 3)]);
        var b = new BassTerm(legendre: [new(0, 3), new(2, 1)]);
        Assert.Equal(a, b);
        Assert.Equal(a.GetHashCode(), b.GetHashCode());
        Assert.Equal(4, a.LegendreDegree);
        Assert.Equal("ψ3(x0) · ψ1(x2)", a.ToString());
        Assert.NotEqual(a, new BassTerm(legendre: [new(0, 3), new(2, 2)]));
        var hinge = new BassTerm([new HingeFactor(0, 1, 2)]);
        Assert.Empty(hinge.Legendre);
        Assert.Equal(0, hinge.LegendreDegree);
        Assert.NotEqual(hinge, new BassTerm([new HingeFactor(0, 1, 2)], legendre: [new(1, 1)]));
    }

    [Fact]
    public void Kernel_EvaluatesTheProduct_AndItsDerivatives()
    {
        BassData data = UnitData(37, 3, 3);
        BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
        var basis = new BassBasis(inputs, 1);
        var term = new BassTerm(legendre: [new(0, 3), new(2, 5)]);
        var column = new double[inputs.Cases.Cases];
        basis.CaseColumn(term, inputs.Cases, column);
        for (int i = 0; i < column.Length; i++)
        {
            double x0 = inputs.Cases.Continuous(0)[i], x2 = inputs.Cases.Continuous(2)[i];
            Assert.Equal(1.0 * ShiftedLegendre.Value(3, x0) * ShiftedLegendre.Value(5, x2), column[i]);
        }
        var support = basis.Support(term);
        Assert.Equal(column.Count(v => v != 0.0), support.CaseNonzero);

        // ∂/∂x₀ and ∂/∂x₂ against the closed forms, and ∂/∂x₁, a variable the term does not hold, zero.
        for (int v = 0; v < 3; v++)
        {
            var derivative = new double[column.Length];
            basis.CaseDerivative(term, v, inputs.Cases, derivative);
            for (int i = 0; i < column.Length; i++)
            {
                double x0 = inputs.Cases.Continuous(0)[i], x2 = inputs.Cases.Continuous(2)[i];
                double expected = v switch
                {
                    0 => ShiftedLegendre.Derivative(3, x0) * ShiftedLegendre.Value(5, x2),
                    2 => ShiftedLegendre.Value(3, x0) * ShiftedLegendre.Derivative(5, x2),
                    _ => 0.0,
                };
                Assert.Equal(expected, derivative[i], 1e-12 * Math.Max(1.0, Math.Abs(expected)));
            }
        }
        IDesignKernel<CaseInputs> legendre = basis.Kernels.Case[2];
        Assert.True(legendre.ReadsVariable(term, 2));
        Assert.False(legendre.ReadsVariable(term, 1));
        Assert.False(basis.Kernels.Case[0].ReadsVariable(term, 0));
        Assert.True(basis.Kernels.Case[0].ReadsVariable(new BassTerm([new HingeFactor(1, -1, 0)]), 1));
        Assert.Throws<ArgumentException>(() => basis.Validate(new BassTerm(legendre: [new(3, 1)])));
    }

    [Fact]
    public void Factors_AreOrthonormal_UnderTheUniformMeasure()
    {
        BassData data = UnitData(20, 2, 4);
        BassInputTransform transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        ProductMeasure measure = ProductMeasure.Resolve(transform, inputs);
        var basis = new BassBasis(inputs, 1);
        IDesignKernel<CaseInputs> legendre = basis.Kernels.Case[2];
        for (int a = 1; a <= 10; a++)
        {
            var fa = (LegendreDesignFactor)legendre.Factors(new BassTerm(legendre: [new(1, a)]), measure).Single();
            Assert.Equal(1.0, legendre.Scale(new BassTerm(legendre: [new(1, a)])));
            double ma = fa.Mean(measure, out double error);
            Assert.Equal(0.0, ma, 1e-11);
            Assert.Equal(0.0, error);
            Assert.Equal(ShiftedLegendre.Value(a, 0.3), fa.Value(0.3));
            Assert.Equal(ShiftedLegendre.Derivative(a, 0.3), fa.Derivative(0.3));
            for (int b = 1; b <= a; b++)
            {
                var fb = (LegendreDesignFactor)legendre.Factors(new BassTerm(legendre: [new(1, b)]), measure).Single();
                double mb = fb.Mean(measure, out _);
                (double cross, double covariance) = fa.Pair(fb, ma, mb, measure, out _);
                Assert.Equal(a == b ? 1.0 : 0.0, cross, 1e-10);
                Assert.Equal(a == b ? 1.0 : 0.0, covariance, 1e-10);
            }
        }
    }

    [Fact]
    public void StructurePrior_IsUniformOverTheSpace()
    {
        var space = new LegendreTermSpace(4, 3, 2);
        var prior = new BassStructurePrior(space, new BassCountPrior(4, 0.5));
        Assert.Same(space, prior.Legendre);
        Assert.True(prior.Has(BassPart.Legendre));
        Assert.False(prior.Has(BassPart.Continuous));
        Assert.Equal(2, prior.Cap(BassPart.Legendre));
        Assert.Equal(1, prior.MinimumOrder(BassPart.Legendre));
        Assert.Equal(4, prior.Variables(BassPart.Legendre));
        foreach (BassTerm t in space.Terms()) Assert.Equal(-Math.Log(30.0), prior.LogPrior(t), 1e-14);
        Assert.Throws<ArgumentException>(() => prior.LogPrior(new BassTerm([new HingeFactor(0, 1, 0)])));
        var rng = new Xoshiro256PlusPlus(9);
        for (int k = 0; k < 100; k++) prior.Validate(prior.DrawTerm(rng)!);

        var counts = new BassTermCounts(prior);
        var t1 = new BassTerm(legendre: [new(0, 1), new(3, 2)]);
        counts.Add(t1);
        counts.Add(new BassTerm(legendre: [new(3, 3)]));
        Assert.Equal(2, counts.Terms);
        Assert.Equal(1, counts.Appearances(BassPart.Legendre, 0));
        Assert.Equal(2, counts.Appearances(BassPart.Legendre, 3));
        Assert.Equal(1, counts.Order(BassPart.Legendre, 2));
        counts.Remove(t1);
        Assert.Equal(1, counts.Appearances(BassPart.Legendre, 3));

        // A hinge fit's prior refuses a Legendre term.
        BassData data = UnitData(12, 2, 2);
        BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
        var hinges = new BassStructurePrior(inputs, 2, 0, 0, CategoricalSubsetMeasure.SizeUniform, new BassCountPrior());
        Assert.Null(hinges.Legendre);
        Assert.Throws<ArgumentException>(() => hinges.Validate(t1));
    }
}
