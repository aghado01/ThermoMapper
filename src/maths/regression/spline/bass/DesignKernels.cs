using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Regression.Spline.Bass.Sensitivity;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's term types behind the design-kernel contract (<see cref="IDesignKernel{TPoints}"/>), bound to a fit's inputs
/// and hinge degree: on the case side the continuous hinge product, the categorical subsets, the Legendre product (E19) and the
/// directional hinge (E36), on the output side the functional hinge product. <see cref="BassBasis"/> composes a term's column from
/// them in this order, and <see cref="BassMoments"/> its factors. A term without a part leaves the column as it is, so a hinge
/// term's column is the same bits with the Legendre and directional kernels present.
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
        Case = [_continuous, new CategoricalSubsetKernel(inputs), new LegendreProductKernel(inputs), new DirectionalHingeKernel(inputs, degree, normalization)];
        Output = [new FunctionalHingeKernel(inputs, degree, normalization)];
    }

    /// <summary>The interval each hinge product is normalized over.</summary>
    public HingeNormalization Normalization { get; }

    /// <summary>Whether every part of <paramref name="term"/> is a product of one-variable factors the moment service integrates
    /// (<see cref="IDesignKernel{TPoints}.Separable"/>): false for a term holding a directional hinge.</summary>
    public bool Separable(BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(term);
        foreach (IDesignKernel<CaseInputs> kernel in Case)
            if (!kernel.Separable(term)) return false;
        foreach (IDesignKernel<OutputInputs> kernel in Output)
            if (!kernel.Separable(term)) return false;
        return true;
    }

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

// The directional hinge (θᵀz − t)₊^q / c on the case side (E36), z the standardized continuous predictors, each read from its scaled
// value as z = a + b·u (BassInputs.Standardizing), so prediction's cases read as training's did. One projection per case reads every
// variable of θ at once, so the part has no separable factors: Factors refuses, and the moment service with it. c is the hinge's
// largest value over each variable's normalization interval, the box's or the training range's, as for BASS's hinges: θᵀz is linear
// in u, so its largest value sits at the corner where each u_j is at its upper end for θ_j > 0 and at its lower end otherwise; where
// the hinge is zero over the whole box, c stays 1, R's guard. A derivative in u_v is Hinge.Derivative's on the projection times
// ∂(θᵀz)/∂u_v = θ_v·b_v, so at the knot at degree one it is the mean of the one-sided limits.
internal sealed class DirectionalHingeKernel(BassInputs inputs, int degree, HingeNormalization normalization) : IDesignKernel<CaseInputs>
{
    public bool Reads(BassTerm term) => term.Directional is not null;

    public bool Differentiable => true;

    public bool Separable(BassTerm term) => term.Directional is null;

    public bool ReadsVariable(BassTerm term, int variable) => term.Directional is { } d && d.Direction.WeightOn(variable) != 0.0;

    public void Validate(BassTerm term)
    {
        if (term.Directional is not { } d) return;
        foreach (int v in d.Direction.Variables)
        {
            if (v >= inputs.Cases.ContinuousCount) throw new ArgumentException($"Continuous variable {v} does not exist.", nameof(term));
            if (inputs.Standardizing(v) is null)
                throw new ArgumentException($"Continuous variable {v}'s scale carries no standardization, the metric a directional hinge reads; fit the transform from the data.", nameof(term));
        }
    }

    public void Multiply(BassTerm term, CaseInputs at, Span<double> column)
    {
        if (term.Directional is not { } d) return;
        double c = Hinge.Power(degree, Peak(d));
        bool divide = c != 0.0;
        double[] projection = Projection(d, at);
        for (int i = 0; i < column.Length; i++)
        {
            double v = column[i] * Hinge.Value(degree, 1, d.Knot, projection[i]);
            column[i] = divide ? v / c : v;
        }
    }

    public void MultiplyDerivative(BassTerm term, int variable, CaseInputs at, Span<double> column)
    {
        if (variable < 0) throw new ArgumentOutOfRangeException(nameof(variable));
        if (!ReadsVariable(term, variable))
        {
            column.Clear();
            return;
        }
        DirectionalHinge d = term.Directional!;
        double slope = d.Direction.WeightOn(variable) * inputs.Standardizing(variable)!.Value.Slope;
        double c = Hinge.Power(degree, Peak(d));
        bool divide = c != 0.0;
        double[] projection = Projection(d, at);
        for (int i = 0; i < column.Length; i++)
        {
            double v = column[i] * (Hinge.Derivative(degree, 1, d.Knot, projection[i]) * slope);
            column[i] = divide ? v / c : v;
        }
    }

    public double Scale(BassTerm term) => term.Directional is { } d ? Hinge.Normalizer(degree, Peak(d)) : 1.0;

    public IReadOnlyList<DesignFactor> Factors(BassTerm term, ProductMeasure measure)
        => term.Directional is null
            ? []
            : throw new NotSupportedException($"Term {term} holds a directional hinge, which reads several inputs through one projection and so has no "
                                              + "one-variable factors: the moment service, which Sobol indices, effects and E8's exact route read, cannot "
                                              + "integrate it. E8 reads its C on the sampled route.");

    // The hinge's largest value at degree one over the normalization box, θᵀz at the maximizing corner less t; 0 where that is below 0,
    // so that the power's guard leaves c at 1.
    internal double Peak(DirectionalHinge d)
    {
        double top = 0.0;
        IReadOnlyList<int> variables = d.Direction.Variables;
        IReadOnlyList<double> weights = d.Direction.Weights;
        for (int k = 0; k < variables.Count; k++)
        {
            int v = variables[k];
            double w = weights[k];
            (double a, double b) = inputs.Standardizing(v)!.Value;
            ReadOnlySpan<double> knots = inputs.ContinuousKnots(v);
            (double lower, double upper) = normalization == HingeNormalization.TrainingRange ? (knots[0], knots[^1]) : (0.0, 1.0);
            top += w * (a + b * (w > 0.0 ? upper : lower));
        }
        return Math.Max(top - d.Knot, 0.0);
    }

    private double[] Projection(DirectionalHinge d, CaseInputs at)
    {
        var projection = new double[at.Cases];
        Project(d.Direction, inputs, at, projection);
        return projection;
    }

    /// <summary>θᵀz at each case of <paramref name="at"/> into <paramref name="projection"/>, z = a + b·u from
    /// <paramref name="inputs"/>' standardization, accumulated over 𝒥 in ascending order: the column's projection, and the knot prior's
    /// (<see cref="DirectionalTermSpace"/>), so the knot's range reads the cases as the column does.</summary>
    internal static void Project(ProjectionDirection direction, BassInputs inputs, CaseInputs at, Span<double> projection)
    {
        projection.Clear();
        IReadOnlyList<int> variables = direction.Variables;
        IReadOnlyList<double> weights = direction.Weights;
        for (int k = 0; k < variables.Count; k++)
        {
            int v = variables[k];
            double w = weights[k];
            (double a, double b) = inputs.Standardizing(v)!.Value;
            ReadOnlySpan<double> u = at.Continuous(v);
            for (int i = 0; i < projection.Length; i++) projection[i] += w * (a + b * u[i]);
        }
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
