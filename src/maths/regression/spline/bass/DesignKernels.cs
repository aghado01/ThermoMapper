using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Regression.Spline.Bass.Sensitivity;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's term types behind the design-kernel contract (<see cref="IDesignKernel{TPoints}"/>), bound to a fit's inputs
/// and hinge degree: on the case side the continuous hinge product, the categorical subsets and the Legendre product (E19), on
/// the output side the functional hinge product. <see cref="BassBasis"/> composes a term's column from them in this order, and
/// <see cref="BassMoments"/> its factors. A term without a part leaves the column as it is, so a hinge term's column is the
/// same bits with the Legendre kernel present.
/// </summary>
public sealed class DesignKernels
{
    private readonly ContinuousHingeKernel _continuous;

    /// <summary>
    /// The kernels of <paramref name="inputs"/> at hinge degree <paramref name="degree"/>, each hinge product normalized to a maximum
    /// of one over <paramref name="normalization"/>'s interval: the box, R's rule, by default.
    /// </summary>
    public DesignKernels(BassInputs inputs, int degree, HingeNormalization normalization = HingeNormalization.Box)
    {
        ArgumentNullException.ThrowIfNull(inputs);
        if (degree < 1) throw new ArgumentOutOfRangeException(nameof(degree), "The degree is a positive integer.");
        if (normalization is not (HingeNormalization.Box or HingeNormalization.TrainingRange))
            throw new ArgumentOutOfRangeException(nameof(normalization), $"No hinge normalization {normalization}.");
        Normalization = normalization;
        _continuous = new ContinuousHingeKernel(inputs, degree, normalization);
        Case = [_continuous, new CategoricalSubsetKernel(inputs), new LegendreProductKernel(inputs)];
        Output = [new FunctionalHingeKernel(inputs, degree, normalization)];
    }

    /// <summary>The interval each hinge product is normalized over.</summary>
    public HingeNormalization Normalization { get; }

    // The base of c for the term's continuous hinges, Πⱼ peakⱼ: what the prefix tables divide by, so their entries take the kernel's rule.
    internal double ContinuousPeak(BassTerm term) => _continuous.Peak(term.Continuous);

    /// <summary>The case-side kernels, in the order a column applies them.</summary>
    public IReadOnlyList<IDesignKernel<CaseInputs>> Case { get; }

    /// <summary>The output-side kernels.</summary>
    public IReadOnlyList<IDesignKernel<OutputInputs>> Output { get; }
}

// A hinge product Πⱼ (sⱼ(xⱼ − tⱼ))₊^q / c, c its largest value over each variable's normalization interval: the box [0, 1] (R's
// makeBasis), or the training range, the variable's knots' ends, which no declared domain moves. The continuous and the functional
// part differ in the variables they read, and share the arithmetic.
internal abstract class HingeProductKernel<TPoints> : IDesignKernel<TPoints>
{
    private readonly HingeNormalization _normalization;

    protected HingeProductKernel(BassInputs inputs, int degree, HingeNormalization normalization)
    {
        Inputs = inputs;
        Degree = degree;
        _normalization = normalization;
    }

    protected BassInputs Inputs { get; }

    protected int Degree { get; }

    protected abstract IReadOnlyList<HingeFactor> Part(BassTerm term);

    protected abstract ReadOnlySpan<double> Values(TPoints at, int variable);

    protected abstract ReadOnlySpan<double> Knots(int variable);

    public bool Reads(BassTerm term) => Part(term).Count > 0;

    public bool Differentiable => true;

    public bool ReadsVariable(BassTerm term, int variable)
    {
        foreach (HingeFactor f in Part(term))
            if (f.Variable == variable) return true;
        return false;
    }

    public abstract void Validate(BassTerm term);

    public abstract IReadOnlyList<DesignFactor> Factors(BassTerm term, ProductMeasure measure);

    // The column times the part's values: each factor's hinge in one pass over the column (Hinge.Scale, four lanes at a time), the last
    // also dividing by c. The same operations per element, in the same order, as the scalar Apply: ((column·h₁)·h₂)…/c.
    public void Multiply(BassTerm term, TPoints at, Span<double> column)
    {
        IReadOnlyList<HingeFactor> factors = Part(term);
        if (factors.Count == 0) return;
        double c = Hinge.Power(Degree, Peak(factors));
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

    // Πⱼ peakⱼ over the factors in order (Hinge.Peak), each over its variable's interval: on the box the same bits as R's
    // Πⱼ ((sⱼ + 1)/2 − sⱼtⱼ). Under the training range the evaluated column is the same function of the raw inputs whatever box
    // the transform reads, since each ratio (xⱼ − tⱼ)/(Uⱼ − tⱼ) cancels the box's affine map; c itself moves with the box.
    internal double Peak(IReadOnlyList<HingeFactor> factors)
    {
        double peak = 1.0;
        foreach (HingeFactor f in factors)
        {
            ReadOnlySpan<double> knots = Knots(f.Variable);
            peak *= _normalization == HingeNormalization.TrainingRange
                ? Hinge.Peak(f.Sign, knots[f.Knot], knots[0], knots[^1])
                : Hinge.Peak(f.Sign, knots[f.Knot], 0.0, 1.0);
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
        foreach (HingeFactor f in factors)
        {
            double t = Knots(f.Variable)[f.Knot];
            ReadOnlySpan<double> x = Values(at, f.Variable);
            if (f.Variable == differentiate)
                for (int i = 0; i < column.Length; i++) column[i] *= Hinge.Derivative(Degree, f.Sign, t, x[i]);
            else
                for (int i = 0; i < column.Length; i++) column[i] *= Hinge.Value(Degree, f.Sign, t, x[i]);
        }
        double c = Hinge.Power(Degree, Peak(factors));
        if (c == 0.0) return;   // a hinge facing its end of the range: zero there, so R leaves c at 1
        for (int i = 0; i < column.Length; i++) column[i] /= c;
    }
}

// The continuous predictors' hinge product, on the case side.
internal sealed class ContinuousHingeKernel(BassInputs inputs, int degree, HingeNormalization normalization)
    : HingeProductKernel<CaseInputs>(inputs, degree, normalization)
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
internal sealed class FunctionalHingeKernel(BassInputs inputs, int degree, HingeNormalization normalization)
    : HingeProductKernel<OutputInputs>(inputs, degree, normalization)
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

// The Legendre product Πⱼ ψ_{αⱼ}(xⱼ) over the continuous predictors' scaled values, on the case side (E19): no normalizer, the
// ψ being orthonormal under the uniform measure already. Each factor's values come in lanes (ShiftedLegendre.Multiply), the
// same bits as the scalar recurrence; a derivative takes the one factor's slope and the others' values.
internal sealed class LegendreProductKernel(BassInputs inputs) : IDesignKernel<CaseInputs>
{
    public bool Reads(BassTerm term) => term.Legendre.Count > 0;

    public bool Differentiable => true;

    public bool ReadsVariable(BassTerm term, int variable)
    {
        foreach (LegendreFactor f in term.Legendre)
            if (f.Variable == variable) return true;
        return false;
    }

    public void Validate(BassTerm term)
    {
        foreach (LegendreFactor f in term.Legendre)
            if (f.Variable >= inputs.Cases.ContinuousCount) throw new ArgumentException($"Continuous variable {f.Variable} does not exist.", nameof(term));
    }

    public void Multiply(BassTerm term, CaseInputs at, Span<double> column)
    {
        foreach (LegendreFactor f in term.Legendre) ShiftedLegendre.Multiply(f.Degree, at.Continuous(f.Variable), column);
    }

    public void MultiplyDerivative(BassTerm term, int variable, CaseInputs at, Span<double> column)
    {
        if (variable < 0) throw new ArgumentOutOfRangeException(nameof(variable));
        if (!ReadsVariable(term, variable))
        {
            column.Clear();
            return;
        }
        foreach (LegendreFactor f in term.Legendre)
        {
            ReadOnlySpan<double> x = at.Continuous(f.Variable);
            if (f.Variable == variable)
                for (int i = 0; i < column.Length; i++) column[i] *= ShiftedLegendre.Derivative(f.Degree, x[i]);
            else ShiftedLegendre.Multiply(f.Degree, x, column);
        }
    }

    public double Scale(BassTerm term) => 1.0;

    public IReadOnlyList<DesignFactor> Factors(BassTerm term, ProductMeasure measure)
    {
        var factors = new List<DesignFactor>(term.Legendre.Count);
        foreach (LegendreFactor f in term.Legendre)
        {
            if (f.Variable >= inputs.Cases.ContinuousCount) throw new ArgumentException($"Term {term} does not index these inputs.", nameof(term));
            factors.Add(new LegendreDesignFactor(measure.OfContinuous(f.Variable), f.Degree));
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
