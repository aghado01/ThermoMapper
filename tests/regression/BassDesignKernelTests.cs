using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The design-kernel contract (F-X3, O-X4): each of today's term types through <see cref="IDesignKernel{TPoints}"/>.
/// <list type="bullet">
/// <item>The kernels' parts compose a term's column, and the separable factors with the kernels' scale reproduce it, on the
/// scalar fixture's hinge, categorical and mixed terms and on a direct functional fit's terms, output side included.</item>
/// <item>Derivatives against central differences, per term type, at interior points and at the training inputs, which are
/// the knots: there a degree-one hinge's derivative jumps and the contract returns the mean of the one-sided limits, what a
/// central difference returns, and a higher degree's is zero. The factor derivative that sensitivity reads
/// (<see cref="BassMoments.FactorDerivative"/>) against the same differences.</item>
/// <item>A categorical part has no derivative; a term without a variable is constant in it.</item>
/// <item>A directional hinge, a kernel that reads several variables, implements the contract as it stands.</item>
/// </list>
/// </summary>
public sealed class BassDesignKernelTests
{
    private static readonly BassData Scalar = SensitivityFixture.Scalar();
    private static readonly BassData Functional = SensitivityFixture.Functional();

    private static BassTerm[] ScalarTerms() =>
    [
        new([new HingeFactor(0, 1, 3)]),
        new([new HingeFactor(0, -1, 6), new HingeFactor(1, 1, 2)]),
        new([new HingeFactor(1, -1, 8)], [new CategoricalFactor(0, [0, 1])]),
        new(categorical: [new CategoricalFactor(0, [0])]),
        new([new HingeFactor(0, 1, 3), new HingeFactor(1, -1, 5)], [new CategoricalFactor(0, [1, 2])]),
        new([new HingeFactor(0, 1, 0)]),   // a knot at the end its hinge faces away from: the guard c = 1 does not apply, c = 1 − 0
        new([new HingeFactor(1, -1, 0)]),  // a hinge facing the low end at the lowest knot: zero on the range, c = 0, left at 1
    ];

    private static BassTerm[] FunctionalTerms() =>
    [
        new([new HingeFactor(0, 1, 2)]),
        new(functional: [new HingeFactor(0, -1, 3)]),
        new([new HingeFactor(1, -1, 4)], functional: [new HingeFactor(0, 1, 1)]),
        new([new HingeFactor(0, 1, 2), new HingeFactor(1, 1, 1)], functional: [new HingeFactor(0, -1, 5)]),
    ];

    private static (BassBasis Basis, BassInputs Inputs, ProductMeasure Measure) Setup(BassData data, int degree, bool outputs = false)
    {
        BassInputTransform transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        return (new BassBasis(inputs, degree), inputs, ProductMeasure.Resolve(transform, inputs, integratedOutputs: outputs ? [.. transform.OutputVariables.Select(v => v.Name)] : null));
    }

    // The case inputs with continuous variable `variable` moved by `by`.
    private static CaseInputs Shifted(CaseInputs at, int variable, double by)
    {
        double[][] continuous = [.. Enumerable.Range(0, at.ContinuousCount).Select(v => at.Continuous(v).ToArray())];
        for (int i = 0; i < continuous[variable].Length; i++) continuous[variable][i] += by;
        return new CaseInputs(continuous, [.. Enumerable.Range(0, at.CategoricalCount).Select(v => at.Categorical(v).ToArray())]);
    }

    private static CaseInputs Interior(CaseInputs at)
    {
        double[][] continuous = [.. Enumerable.Range(0, at.ContinuousCount).Select(v => Enumerable.Range(0, at.Cases).Select(i => 0.05 + 0.0731 * ((i + 3 * v) % 13)).ToArray())];
        return new CaseInputs(continuous, [.. Enumerable.Range(0, at.CategoricalCount).Select(v => at.Categorical(v).ToArray())]);
    }

    private static double[] CentralDifference(BassBasis basis, BassTerm term, CaseInputs at, int variable, double h)
    {
        var up = new double[at.Cases];
        var down = new double[at.Cases];
        basis.CaseColumn(term, Shifted(at, variable, h), up);
        basis.CaseColumn(term, Shifted(at, variable, -h), down);
        return [.. up.Zip(down, (u, d) => (u - d) / (2.0 * h))];
    }

    [Fact]
    public void Kernels_AreTheTermTypes_AndReadTheirParts()
    {
        (BassBasis basis, _, _) = Setup(Scalar, 1);
        IDesignKernel<CaseInputs> hinge = basis.Kernels.Case[0], subsets = basis.Kernels.Case[1];
        Assert.True(hinge.Differentiable);
        Assert.False(subsets.Differentiable);
        var mixed = new BassTerm([new HingeFactor(0, 1, 3)], [new CategoricalFactor(0, [1])]);
        Assert.True(hinge.Reads(mixed) && subsets.Reads(mixed));
        Assert.False(subsets.Reads(new BassTerm([new HingeFactor(0, 1, 3)])));
        Assert.False(hinge.Reads(new BassTerm(categorical: [new CategoricalFactor(0, [1])])));
        Assert.True(basis.Kernels.Output[0].Differentiable);
        Assert.False(basis.Kernels.Output[0].Reads(mixed));

        // Each part supplies its Validate: the basis refuses what a kernel refuses, with its message.
        Assert.Throws<ArgumentException>(() => hinge.Validate(new BassTerm([new HingeFactor(5, 1, 0)])));
        Assert.Throws<ArgumentException>(() => hinge.Validate(new BassTerm([new HingeFactor(0, 1, 99)])));
        Assert.Throws<ArgumentException>(() => subsets.Validate(new BassTerm(categorical: [new CategoricalFactor(0, [0, 1, 2])])));
        Assert.Throws<ArgumentException>(() => basis.Kernels.Output[0].Validate(new BassTerm(functional: [new HingeFactor(0, 1, 0)])));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void ScalarTerms_FactorsWithTheScale_ReproduceTheColumn(int degree)
    {
        (BassBasis basis, BassInputs inputs, ProductMeasure measure) = Setup(Scalar, degree);
        var moments = BassMoments.On(inputs, degree, measure);
        foreach (BassTerm term in ScalarTerms())
        {
            var column = new double[inputs.Cases.Cases];
            basis.CaseColumn(term, inputs.Cases, column);
            IReadOnlyList<int> variables = moments.Variables(term);
            double kappa = moments.Scale(term);
            for (int i = 0; i < column.Length; i++)
            {
                double value = kappa;
                for (int k = 0; k < variables.Count; k++) value *= moments.FactorValue(term, k, At(inputs, measure, variables[k], i));
                Assert.True(Math.Abs(value - column[i]) <= 1e-13 * (1.0 + Math.Abs(column[i])), $"{term} at case {i}, q = {degree}: {value} against {column[i]}");
            }
        }
    }

    // The value variable `variable` of the measure takes at case `i`: a scaled continuous value, or a level index.
    private static double At(BassInputs inputs, ProductMeasure measure, int variable, int i)
    {
        for (int k = 0; k < inputs.Cases.ContinuousCount; k++)
            if (measure.OfContinuous(k) == variable) return inputs.Cases.Continuous(k)[i];
        for (int k = 0; k < inputs.Cases.CategoricalCount; k++)
            if (measure.OfCategorical(k) == variable) return inputs.Cases.Categorical(k)[i];
        throw new ArgumentException("An output variable has no case value.");
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void FunctionalTerms_FactorsWithTheScale_ReproduceBothSides(int degree)
    {
        (BassBasis basis, BassInputs inputs, ProductMeasure measure) = Setup(Functional, degree, outputs: true);
        var moments = BassMoments.On(inputs, degree, measure);
        OutputInputs outputs = inputs.Outputs!;
        foreach (BassTerm term in FunctionalTerms())
        {
            var a = new double[inputs.Cases.Cases];
            var h = new double[outputs.Locations];
            basis.CaseColumn(term, inputs.Cases, a);
            basis.OutputColumn(term, outputs, h);
            IReadOnlyList<int> variables = moments.Variables(term);
            double kappa = moments.Scale(term);
            for (int i = 0; i < a.Length; i++)
                for (int r = 0; r < h.Length; r++)
                {
                    double value = kappa;
                    for (int k = 0; k < variables.Count; k++)
                    {
                        int v = variables[k];
                        value *= moments.FactorValue(term, k, v < inputs.Cases.ContinuousCount ? inputs.Cases.Continuous(v)[i] : outputs.Variable(v - inputs.Cases.ContinuousCount)[r]);
                    }
                    double column = a[i] * h[r];
                    Assert.True(Math.Abs(value - column) <= 1e-13 * (1.0 + Math.Abs(column)), $"{term} at ({i}, {r}), q = {degree}: {value} against {column}");
                }
        }
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void CaseDerivatives_MatchCentralDifferences_AtInteriorPoints(int degree)
    {
        (BassBasis basis, BassInputs inputs, _) = Setup(Scalar, degree);
        CaseInputs at = Interior(inputs.Cases);
        foreach (BassTerm term in ScalarTerms())
            for (int v = 0; v < at.ContinuousCount; v++)
            {
                var derivative = new double[at.Cases];
                basis.CaseDerivative(term, v, at, derivative);
                double[] difference = CentralDifference(basis, term, at, v, 1e-6);
                for (int i = 0; i < derivative.Length; i++)
                    Assert.True(Math.Abs(derivative[i] - difference[i]) <= 1e-6 * (1.0 + Math.Abs(derivative[i])),
                                $"d{term}/dx{v} at case {i}, q = {degree}: {derivative[i]} against {difference[i]}");
            }
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void CaseDerivatives_AtTheKnots_AreTheMeanOfTheOneSidedLimits(int degree)
    {
        // The training inputs are the knots, so every case sits at a kink of the hinges at its own values.
        (BassBasis basis, BassInputs inputs, _) = Setup(Scalar, degree);
        CaseInputs at = inputs.Cases;
        bool atAKnot = false;
        foreach (BassTerm term in ScalarTerms())
            for (int v = 0; v < at.ContinuousCount; v++)
            {
                var derivative = new double[at.Cases];
                basis.CaseDerivative(term, v, at, derivative);
                double[] difference = CentralDifference(basis, term, at, v, 1e-7);
                for (int i = 0; i < derivative.Length; i++)
                {
                    Assert.True(Math.Abs(derivative[i] - difference[i]) <= 5e-5 * (1.0 + Math.Abs(derivative[i])),
                                $"d{term}/dx{v} at the knot of case {i}, q = {degree}: {derivative[i]} against {difference[i]}");
                    // Where the hinge on v has its knot at the case's value, a degree-one derivative jumps there.
                    foreach (HingeFactor h in term.Continuous)
                        if (h.Variable == v && at.Continuous(v)[i] == inputs.ContinuousKnots(v)[h.Knot]) atAKnot = true;
                }
            }
        Assert.True(atAKnot, "the fixture must put a hinge's knot under a case");
    }

    [Fact]
    public void ADegreeOneHinge_AtItsKnot_HasHalfItsSlope_AndAHigherDegreeNone()
    {
        foreach (int sign in new[] { 1, -1 })
        {
            Assert.Equal(0.5 * sign, Hinge.Derivative(1, sign, 0.4, 0.4));
            Assert.Equal(sign, Hinge.Derivative(1, sign, 0.4, 0.4 + sign * 0.1));
            Assert.Equal(0.0, Hinge.Derivative(1, sign, 0.4, 0.4 - sign * 0.1));
            foreach (int degree in new[] { 2, 3, 4 })
            {
                Assert.Equal(0.0, Hinge.Derivative(degree, sign, 0.4, 0.4));
                Assert.Equal(degree * sign * Math.Pow(0.1, degree - 1), Hinge.Derivative(degree, sign, 0.4, 0.4 + sign * 0.1), 12);
            }
        }
    }

    [Fact]
    public void ACategoricalPart_HasNoDerivative_AndATermWithoutTheVariableIsConstantInIt()
    {
        (BassBasis basis, BassInputs inputs, ProductMeasure measure) = Setup(Scalar, 1);
        var subsetOnly = new BassTerm(categorical: [new CategoricalFactor(0, [0])]);
        Assert.Throws<NotSupportedException>(() => basis.Kernels.Case[1].MultiplyDerivative(subsetOnly, 0, inputs.Cases, new double[inputs.Cases.Cases]));
        for (int v = 0; v < 2; v++)
        {
            var column = new double[inputs.Cases.Cases];
            Array.Fill(column, 7.0);
            basis.CaseDerivative(subsetOnly, v, inputs.Cases, column);
            Assert.All(column, x => Assert.Equal(0.0, x));
        }
        // A hinge on x0 alone is constant in x1.
        var hinge = new BassTerm([new HingeFactor(0, 1, 3)]);
        var d = new double[inputs.Cases.Cases];
        basis.CaseDerivative(hinge, 1, inputs.Cases, d);
        Assert.All(d, x => Assert.Equal(0.0, x));

        var moments = BassMoments.On(inputs, 1, measure);
        var mixed = new BassTerm([new HingeFactor(0, 1, 3)], [new CategoricalFactor(0, [1])]);
        Assert.Equal(2, moments.Variables(mixed).Count);
        Assert.Throws<NotSupportedException>(() => moments.FactorDerivative(mixed, 1, 1.0));   // the level factor, after the x0 hinge
        Assert.Throws<ArgumentOutOfRangeException>(() => basis.CaseDerivative(hinge, 2, inputs.Cases, d));
        Assert.Throws<ArgumentException>(() => basis.CaseDerivative(hinge, 0, inputs.Cases, new double[3]));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void FactorDerivatives_ThatSensitivityReads_MatchCentralDifferences_AtInteriorPointsAndKnots(int degree)
    {
        (BassBasis basis, BassInputs inputs, ProductMeasure measure) = Setup(Scalar, degree);
        var moments = BassMoments.On(inputs, degree, measure);
        BassTerm term = new([new HingeFactor(0, 1, 3), new HingeFactor(1, -1, 5)], [new CategoricalFactor(0, [1, 2])]);
        IReadOnlyList<int> variables = moments.Variables(term);
        for (int k = 0; k < 2; k++)   // the two hinges; the third factor is the level indicator
        {
            double knot = inputs.ContinuousKnots(variables[k])[term.Continuous[k].Knot];
            double[] values = [0.05, 0.2, 0.37, knot - 0.03, knot, knot + 0.03, 0.81, 0.97];
            foreach (double x in values)
            {
                const double h = 1e-7;
                double difference = (moments.FactorValue(term, k, x + h) - moments.FactorValue(term, k, x - h)) / (2.0 * h);
                Assert.True(Math.Abs(moments.FactorDerivative(term, k, x) - difference) <= 1e-5 * (1.0 + Math.Abs(difference)),
                            $"factor {k} at {x}, q = {degree}: {moments.FactorDerivative(term, k, x)} against {difference}");
            }
        }
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void OutputDerivatives_MatchCentralDifferences_AtInteriorPointsAndTheKnots(int degree)
    {
        (BassBasis basis, BassInputs inputs, _) = Setup(Functional, degree);
        OutputInputs training = inputs.Outputs!;
        double[] interior = [.. Enumerable.Range(0, training.Locations).Select(r => 0.04 + 0.131 * r)];
        foreach (double[] t in new[] { training.Variable(0).ToArray(), interior })
        {
            var at = new OutputInputs([t]);
            foreach (BassTerm term in FunctionalTerms())
            {
                var derivative = new double[at.Locations];
                basis.OutputDerivative(term, 0, at, derivative);
                const double h = 1e-7;
                var up = new double[at.Locations];
                var down = new double[at.Locations];
                basis.OutputColumn(term, new OutputInputs([[.. t.Select(x => x + h)]]), up);
                basis.OutputColumn(term, new OutputInputs([[.. t.Select(x => x - h)]]), down);
                for (int r = 0; r < derivative.Length; r++)
                {
                    double difference = (up[r] - down[r]) / (2.0 * h);
                    Assert.True(Math.Abs(derivative[r] - difference) <= 5e-5 * (1.0 + Math.Abs(derivative[r])),
                                $"d{term}/dt at {t[r]}, q = {degree}: {derivative[r]} against {difference}");
                }
            }
        }
    }

    [Fact]
    public void AFunctionalFits_Terms_GoThroughTheContract_OnBothSides()
    {
        BassFit fit = SensitivityFixture.FitFunctional();
        int degree = fit.Config.Structure.Degree;
        BassInputTransform transform = BassInputTransform.Fit(Functional);
        BassInputs inputs = BassInputs.Prepare(Functional, transform);
        var basis = new BassBasis(inputs, degree);
        var measure = ProductMeasure.Resolve(transform, inputs, integratedOutputs: [.. transform.OutputVariables.Select(v => v.Name)]);
        var moments = BassMoments.On(inputs, degree, measure);
        int checkedTerms = 0, withFunctional = 0;
        for (int id = 0; id < fit.Terms.Count; id++)
        {
            BassTerm term = fit.Terms[id];
            basis.Validate(term);
            var a = new double[inputs.Cases.Cases];
            var h = new double[inputs.Outputs!.Locations];
            basis.CaseColumn(term, inputs.Cases, a);
            basis.OutputColumn(term, inputs.Outputs, h);
            IReadOnlyList<int> variables = moments.Variables(term);
            double kappa = moments.Scale(term);
            for (int i = 0; i < a.Length; i += 7)
                for (int r = 0; r < h.Length; r++)
                {
                    double value = kappa;
                    for (int k = 0; k < variables.Count; k++)
                    {
                        int v = variables[k];
                        value *= moments.FactorValue(term, k, v < inputs.Cases.ContinuousCount ? inputs.Cases.Continuous(v)[i] : inputs.Outputs.Variable(v - inputs.Cases.ContinuousCount)[r]);
                    }
                    Assert.True(Math.Abs(value - a[i] * h[r]) <= 1e-12 * (1.0 + Math.Abs(a[i] * h[r])), $"{term} at ({i}, {r})");
                }
            if (term.Functional.Count > 0)
            {
                withFunctional++;
                var d = new double[inputs.Outputs.Locations];
                basis.OutputDerivative(term, 0, inputs.Outputs, d);
                var up = new double[d.Length];
                var down = new double[d.Length];
                const double eps = 1e-7;
                basis.OutputColumn(term, new OutputInputs([[.. inputs.Outputs.Variable(0).ToArray().Select(x => x + eps)]]), up);
                basis.OutputColumn(term, new OutputInputs([[.. inputs.Outputs.Variable(0).ToArray().Select(x => x - eps)]]), down);
                for (int r = 0; r < d.Length; r++)
                    Assert.True(Math.Abs(d[r] - (up[r] - down[r]) / (2.0 * eps)) <= 5e-5 * (1.0 + Math.Abs(d[r])), $"d{term}/dt at location {r}");
            }
            checkedTerms++;
        }
        Assert.True(checkedTerms > 0 && withFunctional > 0, "the functional fit should hold terms with an output part");
    }

    // (wᵀx − t)₊ on a case-side kernel: it reads several variables, so it is no product of one-variable factors and
    // declares none, but the contract's column evaluation and its partial derivatives hold as they stand.
    private sealed class DirectionalHinge(double[] w, double knot) : IDesignKernel<CaseInputs>
    {
        public bool Reads(BassTerm term) => true;
        public bool Differentiable => true;
        public void Validate(BassTerm term) { }

        public void Multiply(BassTerm term, CaseInputs at, Span<double> column)
        {
            for (int i = 0; i < column.Length; i++) column[i] *= Math.Max(Projection(at, i) - knot, 0.0);
        }

        public void MultiplyDerivative(BassTerm term, int variable, CaseInputs at, Span<double> column)
        {
            for (int i = 0; i < column.Length; i++) column[i] *= Projection(at, i) - knot > 0.0 ? w[variable] : 0.0;
        }

        public double Scale(BassTerm term) => 1.0;

        public IReadOnlyList<DesignFactor> Factors(BassTerm term, ProductMeasure measure) => [];

        private double Projection(CaseInputs at, int i)
        {
            double sum = 0.0;
            for (int v = 0; v < w.Length; v++) sum += w[v] * at.Continuous(v)[i];
            return sum;
        }
    }

    [Fact]
    public void ADirectionalHinge_FitsTheContractAsItStands()
    {
        (_, BassInputs inputs, _) = Setup(Scalar, 1);
        IDesignKernel<CaseInputs> kernel = new DirectionalHinge([0.6, -0.8], 0.1);
        CaseInputs at = Interior(inputs.Cases);
        var term = new BassTerm([new HingeFactor(0, 1, 0)]);
        for (int v = 0; v < 2; v++)
        {
            var derivative = new double[at.Cases];
            Array.Fill(derivative, 1.0);
            kernel.MultiplyDerivative(term, v, at, derivative);
            var up = new double[at.Cases];
            var down = new double[at.Cases];
            Array.Fill(up, 1.0);
            Array.Fill(down, 1.0);
            kernel.Multiply(term, Shifted(at, v, 1e-6), up);
            kernel.Multiply(term, Shifted(at, v, -1e-6), down);
            for (int i = 0; i < up.Length; i++) Assert.Equal((up[i] - down[i]) / 2e-6, derivative[i], 6);
        }
    }
}
