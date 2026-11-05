using System;
using System.Collections.Generic;
using Maths.Regression.Spline.Bass.Sensitivity;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's term types behind the design-kernel contract (<see cref="IDesignKernel{TPoints}"/>), bound to a fit's inputs
/// and hinge degree: on the case side the continuous hinge product and the categorical subsets, on the output side the
/// functional hinge product. <see cref="BassBasis"/> composes a term's column from them in this order, and
/// <see cref="BassMoments"/> its factors.
/// </summary>
public sealed class DesignKernels
{
    /// <summary>The kernels of <paramref name="inputs"/> at hinge degree <paramref name="degree"/>.</summary>
    public DesignKernels(BassInputs inputs, int degree)
    {
        ArgumentNullException.ThrowIfNull(inputs);
        if (degree < 1) throw new ArgumentOutOfRangeException(nameof(degree), "The degree is a positive integer.");
        Case = [new ContinuousHingeKernel(inputs, degree), new CategoricalSubsetKernel(inputs)];
        Output = [new FunctionalHingeKernel(inputs, degree)];
    }

    /// <summary>The case-side kernels, in the order a column applies them.</summary>
    public IReadOnlyList<IDesignKernel<CaseInputs>> Case { get; }

    /// <summary>The output-side kernels.</summary>
    public IReadOnlyList<IDesignKernel<OutputInputs>> Output { get; }
}

// A hinge product Πⱼ (sⱼ(xⱼ − tⱼ))₊^q / c, c its largest value over [0, 1] (R's makeBasis): the continuous and the
// functional part differ in the variables they read, and share the arithmetic.
internal abstract class HingeProductKernel<TPoints> : IDesignKernel<TPoints>
{
    protected HingeProductKernel(BassInputs inputs, int degree)
    {
        Inputs = inputs;
        Degree = degree;
    }

    protected BassInputs Inputs { get; }

    protected int Degree { get; }

    protected abstract IReadOnlyList<HingeFactor> Part(BassTerm term);

    protected abstract ReadOnlySpan<double> Values(TPoints at, int variable);

    protected abstract ReadOnlySpan<double> Knots(int variable);

    public bool Reads(BassTerm term) => Part(term).Count > 0;

    public bool Differentiable => true;

    public abstract void Validate(BassTerm term);

    public abstract IReadOnlyList<DesignFactor> Factors(BassTerm term, ProductMeasure measure);

    // The column times the part's values: each factor's hinge in one pass over the column (Hinge.Scale, four lanes at a time), the last
    // also dividing by c. The same operations per element, in the same order, as the scalar Apply: ((column·h₁)·h₂)…/c.
    public void Multiply(BassTerm term, TPoints at, Span<double> column)
    {
        IReadOnlyList<HingeFactor> factors = Part(term);
        if (factors.Count == 0) return;
        double peak = 1.0;
        foreach (HingeFactor f in factors)
        {
            double t = Knots(f.Variable)[f.Knot], s = f.Sign;
            peak *= (s + 1.0) / 2.0 - s * t;
        }
        double c = Hinge.Power(Degree, peak);
        bool divide = c != 0.0;   // a hinge facing its end of the range: zero there, so R leaves c at 1
        int last = factors.Count - 1;
        for (int k = 0; k <= last; k++)
        {
            HingeFactor f = factors[k];
            Hinge.Scale(column, Values(at, f.Variable), Degree, f.Sign, Knots(f.Variable)[f.Knot], divide && k == last, c);
        }
    }

    public void MultiplyDerivative(BassTerm term, int variable, TPoints at, Span<double> column)
    {
        if (variable < 0) throw new ArgumentOutOfRangeException(nameof(variable));
        Apply(term, variable, at, column);
    }

    public double Scale(BassTerm term) => Hinge.Normalizer(Degree, Peak(Part(term)));

    // Πⱼ ((sⱼ + 1)/2 − sⱼtⱼ) over the factors in variable order.
    private double Peak(IReadOnlyList<HingeFactor> factors)
    {
        double peak = 1.0;
        foreach (HingeFactor f in factors)
        {
            double s = f.Sign;
            peak *= (s + 1.0) / 2.0 - s * Knots(f.Variable)[f.Knot];
        }
        return peak;
    }

    // The column times the part's derivative in variable `differentiate` ≥ 0, scalar, element by element: the product rule has
    // one nonzero term, since a hinge product holds each variable once, and a part that lacks the variable is constant in it.
    // (The values, differentiate < 0, are Multiply's.)
    private void Apply(BassTerm term, int differentiate, TPoints at, Span<double> column)
    {
        IReadOnlyList<HingeFactor> factors = Part(term);
        if (differentiate >= 0)
        {
            bool reads = false;
            foreach (HingeFactor f in factors) reads |= f.Variable == differentiate;
            if (!reads)
            {
                column.Clear();
                return;
            }
        }
        if (factors.Count == 0) return;
        double peak = 1.0;
        foreach (HingeFactor f in factors)
        {
            double t = Knots(f.Variable)[f.Knot], s = f.Sign;
            peak *= (s + 1.0) / 2.0 - s * t;
            ReadOnlySpan<double> x = Values(at, f.Variable);
            if (f.Variable == differentiate)
                for (int i = 0; i < column.Length; i++) column[i] *= Hinge.Derivative(Degree, f.Sign, t, x[i]);
            else
                for (int i = 0; i < column.Length; i++) column[i] *= Hinge.Value(Degree, f.Sign, t, x[i]);
        }
        double c = Hinge.Power(Degree, peak);
        if (c == 0.0) return;   // a hinge facing its end of the range: zero there, so R leaves c at 1
        for (int i = 0; i < column.Length; i++) column[i] /= c;
    }
}

// The continuous predictors' hinge product, on the case side.
internal sealed class ContinuousHingeKernel(BassInputs inputs, int degree) : HingeProductKernel<CaseInputs>(inputs, degree)
{
    protected override IReadOnlyList<HingeFactor> Part(BassTerm term) => term.Continuous;

    protected override ReadOnlySpan<double> Values(CaseInputs at, int variable) => at.Continuous(variable);

    protected override ReadOnlySpan<double> Knots(int variable) => Inputs.ContinuousKnots(variable);

    public override void Validate(BassTerm term)
    {
        foreach (HingeFactor h in term.Continuous)
        {
            if (h.Variable >= Inputs.Cases.ContinuousCount) throw new ArgumentException($"Continuous variable {h.Variable} does not exist.", nameof(term));
            if (h.Knot >= Inputs.ContinuousKnots(h.Variable).Length) throw new ArgumentException($"Knot {h.Knot} of continuous variable {h.Variable} does not exist.", nameof(term));
        }
    }

    public override IReadOnlyList<DesignFactor> Factors(BassTerm term, ProductMeasure measure)
    {
        var factors = new List<DesignFactor>(term.Continuous.Count);
        foreach (HingeFactor h in term.Continuous)
        {
            if (h.Variable >= Inputs.Cases.ContinuousCount || h.Knot >= Inputs.ContinuousKnots(h.Variable).Length)
                throw new ArgumentException($"Term {term} does not index these inputs.", nameof(term));
            factors.Add(new HingeDesignFactor(measure.OfContinuous(h.Variable), h.Sign, h.Knot, Inputs.ContinuousKnots(h.Variable)[h.Knot], Degree));
        }
        return factors;
    }
}

// The output variables' hinge product, on the output side of a direct functional response.
internal sealed class FunctionalHingeKernel(BassInputs inputs, int degree) : HingeProductKernel<OutputInputs>(inputs, degree)
{
    protected override IReadOnlyList<HingeFactor> Part(BassTerm term) => term.Functional;

    protected override ReadOnlySpan<double> Values(OutputInputs at, int variable) => at.Variable(variable);

    protected override ReadOnlySpan<double> Knots(int variable) => Inputs.FunctionalKnots(variable);

    public override void Validate(BassTerm term)
    {
        if (term.Functional.Count == 0) return;
        OutputInputs outputs = Inputs.Outputs ?? throw new ArgumentException("A scalar response has no functional factors.", nameof(term));
        foreach (HingeFactor h in term.Functional)
        {
            if (h.Variable >= outputs.Count) throw new ArgumentException($"Output variable {h.Variable} does not exist.", nameof(term));
            if (h.Knot >= Inputs.FunctionalKnots(h.Variable).Length) throw new ArgumentException($"Knot {h.Knot} of output variable {h.Variable} does not exist.", nameof(term));
        }
    }

    public override IReadOnlyList<DesignFactor> Factors(BassTerm term, ProductMeasure measure)
    {
        var factors = new List<DesignFactor>(term.Functional.Count);
        foreach (HingeFactor h in term.Functional)
        {
            if (Inputs.Outputs is null || h.Variable >= Inputs.Outputs.Count || h.Knot >= Inputs.FunctionalKnots(h.Variable).Length)
                throw new ArgumentException($"Term {term} does not index these inputs.", nameof(term));
            int variable = measure.OfOutput(h.Variable);
            if (variable >= 0) factors.Add(new HingeDesignFactor(variable, h.Sign, h.Knot, Inputs.FunctionalKnots(h.Variable)[h.Knot], Degree));
        }
        return factors;
    }
}

// The categorical predictors' subset indicators Πₗ 1(zₗ ∈ Dₗ), on the case side: values 0 or 1, with no derivative.
internal sealed class CategoricalSubsetKernel(BassInputs inputs) : IDesignKernel<CaseInputs>
{
    public bool Reads(BassTerm term) => term.Categorical.Count > 0;

    public bool Differentiable => false;

    public void Validate(BassTerm term)
    {
        foreach (CategoricalFactor c in term.Categorical)
        {
            if (c.Variable >= inputs.Cases.CategoricalCount) throw new ArgumentException($"Categorical variable {c.Variable} does not exist.", nameof(term));
            int levels = inputs.LevelCount(c.Variable);
            if (c.Levels[^1] >= levels) throw new ArgumentException($"Categorical variable {c.Variable} has {levels} levels.", nameof(term));
            if (c.Levels.Count >= levels) throw new ArgumentException($"A subset of categorical variable {c.Variable} must leave a level out.", nameof(term));
        }
    }

    public void Multiply(BassTerm term, CaseInputs at, Span<double> column)
    {
        foreach (CategoricalFactor c in term.Categorical)
        {
            c.ZeroOutside(at.Categorical(c.Variable), column);
        }
    }

    public void MultiplyDerivative(BassTerm term, int variable, CaseInputs at, Span<double> column)
        => throw new NotSupportedException("A categorical part's variables are levels; it has no derivative.");

    public double Scale(BassTerm term) => 1.0;

    public IReadOnlyList<DesignFactor> Factors(BassTerm term, ProductMeasure measure)
    {
        var factors = new List<DesignFactor>(term.Categorical.Count);
        foreach (CategoricalFactor c in term.Categorical)
        {
            if (c.Variable >= inputs.Cases.CategoricalCount || c.Levels[^1] >= inputs.LevelCount(c.Variable))
                throw new ArgumentException($"Term {term} does not index these inputs.", nameof(term));
            factors.Add(new SubsetDesignFactor(measure.OfCategorical(c.Variable), c));
        }
        return factors;
    }
}
