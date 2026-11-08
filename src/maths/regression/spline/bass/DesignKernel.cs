using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using Maths.Distributions;
using Maths.Regression.Spline.Bass.Sensitivity;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// The design-kernel contract (F-X3): what the design, the posterior layer and sensitivity need of one term type. A
/// term is a product of parts, each part of a type that sits behind this contract on its side of the Khatri–Rao column
/// (<see cref="CaseInputs"/> for the case side, <see cref="OutputInputs"/> for the output side), and a part supplies:
/// <list type="bullet">
/// <item>its <b>column evaluation</b>, <see cref="Multiply"/>, which scales the running product of the term's earlier
/// parts, so the term's column is the parts applied in order to a column of ones, and the term's support (its nonzero
/// count and mass, <see cref="BassBasis.Support"/>) is read off that product;</item>
/// <item>its <b>derivatives</b>, <see cref="MultiplyDerivative"/>, one variable at a time, which a differentiable part
/// supplies in place of its values (<see cref="BassBasis.CaseDerivative"/>); at a hinge's knot, where the derivative of
/// degree one jumps, the derivative is the mean of its one-sided limits, which is what a central difference returns;</item>
/// <item>its <b>moments</b>, as the separable <see cref="DesignFactor"/>s of its part with the normalizing
/// <see cref="Scale"/>, so that <see cref="BassMoments"/> serves every type's means, cross moments and ANOVA
/// structure by one rule: a term is κ·Πᵥ φᵥ(xᵥ).</item>
/// </list>
/// A variable belongs to one part of a term, so a derivative takes the differentiable part's derivative and the other
/// parts' values (the product rule has one nonzero term), and the part that reads it says so (<see cref="ReadsVariable"/>).
/// Today's members are the continuous hinge product, the categorical subsets, the functional hinge product and the Legendre
/// product of a polynomial-chaos fit (E19; <see cref="DesignKernels"/>).
/// <para>Later members join without changing the contract. A tensor-product B-spline part (ABMRS) is a case-side
/// kernel whose <see cref="Multiply"/> writes the product of one-variable B-splines, whose derivative differentiates one
/// factor, and whose factors are one-variable piecewise polynomials, which <see cref="PiecewisePolynomial"/> integrates
/// as it does a hinge. A directional hinge (learned projections) is a case-side kernel for (wᵀx − t)₊, differentiable in
/// every variable of w; it is not a product of one-variable factors, so it declares no <see cref="Factors"/> and
/// supplies its moments through a further interface the moment service tests for, which the active-subspace work writes.</para>
/// </summary>
/// <typeparam name="TPoints">The inputs the part is evaluated at: <see cref="CaseInputs"/> or <see cref="OutputInputs"/>.</typeparam>
public interface IDesignKernel<in TPoints>
{
    /// <summary>Whether <paramref name="term"/> has a part of this type.</summary>
    bool Reads(BassTerm term);

    /// <summary>Whether the part is a function of continuous variables and so has derivatives; the levels of a
    /// categorical part have none.</summary>
    bool Differentiable { get; }

    /// <summary>Refuses a part whose indices fall outside the inputs.</summary>
    void Validate(BassTerm term);

    /// <summary>Whether <paramref name="term"/>'s part of this type reads continuous variable <paramref name="variable"/> of this
    /// side, so that the term's partial derivative in it is this part's derivative times the other parts' values. The default
    /// is a differentiable part's <see cref="Reads"/>, for a part that reads every variable it is handed.</summary>
    bool ReadsVariable(BassTerm term, int variable) => Differentiable && Reads(term);

    /// <summary>Multiplies <paramref name="column"/>, one entry per row of <paramref name="at"/>, by the part's values; a term
    /// without the part leaves it as it is.</summary>
    void Multiply(BassTerm term, TPoints at, Span<double> column);

    /// <summary>Multiplies <paramref name="column"/> by the part's partial derivative in <paramref name="variable"/>, a
    /// continuous variable of this side by its position among its kind; zero for a part that does not read it. Only for a
    /// <see cref="Differentiable"/> part.</summary>
    void MultiplyDerivative(BassTerm term, int variable, TPoints at, Span<double> column);

    /// <summary>The constant κ the part's values carry, 1 for a part without one; the moment service multiplies the factors'
    /// means by it, where <see cref="Multiply"/> divides by its reciprocal.</summary>
    double Scale(BassTerm term);

    /// <summary>The part's separable factors, each on one variable of <paramref name="measure"/>; one on a variable the
    /// measure does not cover is left out.</summary>
    IReadOnlyList<DesignFactor> Factors(BassTerm term, ProductMeasure measure);
}

/// <summary>
/// One one-variable factor φ of a term, on a variable of a <see cref="ProductMeasure"/>: its value and derivative at a
/// scaled value, and its moments under the measure. A value type, so that <see cref="BassMoments"/> caches a factor's
/// moments, and a pair's, by the factor itself.
/// </summary>
public abstract record DesignFactor(int Variable)
{
    /// <summary>φ(x), unnormalized: the term's κ is separate.</summary>
    public abstract double Value(double x);

    /// <summary>φ′(x), the mean of the one-sided limits at a kink; refused for a factor on a level index.</summary>
    public abstract double Derivative(double x);

    /// <summary>E[φ] under <paramref name="measure"/>, with the quadrature error estimate of the moment (0 when exact).</summary>
    public abstract double Mean(ProductMeasure measure, out double error);

    /// <summary>E[φ·φ′] and Cov(φ, φ′) for this factor and <paramref name="other"/> on the same variable, given their means.</summary>
    public abstract (double Cross, double Covariance) Pair(DesignFactor other, double meanThis, double meanOther,
                                                           ProductMeasure measure, out double error);
}

/// <summary>A hinge (s(x − t))₊^q as a factor: the sign, the knot's index and value, and the degree.</summary>
public sealed record HingeDesignFactor(int Variable, int Sign, int Knot, double KnotValue, int Degree) : DesignFactor(Variable)
{
    /// <inheritdoc/>
    public override double Value(double x) => Hinge.Value(Degree, Sign, KnotValue, x);

    /// <inheritdoc/>
    public override double Derivative(double x) => Hinge.Derivative(Degree, Sign, KnotValue, x);

    /// <inheritdoc/>
    public override double Mean(ProductMeasure measure, out double error)
    {
        ContinuousMeasure m = measure.Continuous(Variable)!;
        (double mean, error) = PiecewisePolynomial.HingeProduct(Degree, [(Sign, KnotValue)]).Expectation(m.Measure, m.Lower, m.Upper);
        return mean;
    }

    /// <inheritdoc/>
    public override (double Cross, double Covariance) Pair(DesignFactor other, double meanThis, double meanOther,
                                                           ProductMeasure measure, out double error)
    {
        if (other is not HingeDesignFactor b || b.Variable != Variable) throw new ArgumentException("A hinge pairs with a hinge on its variable.", nameof(other));
        HingeDesignFactor a = this;
        double meanA = meanThis, meanB = meanOther;
        if (a.Sign > b.Sign || (a.Sign == b.Sign && a.Knot > b.Knot))
        {
            (a, b) = (b, a);
            (meanA, meanB) = (meanB, meanA);
        }
        ContinuousMeasure m = measure.Continuous(Variable)!;
        (int, double)[] hinges = [(a.Sign, a.KnotValue), (b.Sign, b.KnotValue)];
        (double cross, double ec) = PiecewisePolynomial.HingeProduct(Degree, hinges).Expectation(m.Measure, m.Lower, m.Upper);
        (double covariance, double ev) = PiecewisePolynomial.HingeProduct(Degree, hinges, [meanA, meanB]).Expectation(m.Measure, m.Lower, m.Upper);
        error = Math.Max(ec, ev);
        return (cross, covariance);
    }
}

/// <summary>
/// A standardized shifted Legendre polynomial ψ_α as a factor (E19): one piece over the line, whose expansion about any centre
/// <see cref="ShiftedLegendre.Coefficients"/> writes, so <see cref="PiecewisePolynomial"/> takes its moments under any measure
/// as it does a hinge's. Under the uniform measure on [0, 1] its mean is 0 for α ≥ 1 and its pairs are orthonormal; the
/// expansion's coefficients grow with α, so a moment past degree 15 or so loses digits to cancellation.
/// </summary>
public sealed record LegendreDesignFactor(int Variable, int Degree) : DesignFactor(Variable)
{
    /// <inheritdoc/>
    public override double Value(double x) => ShiftedLegendre.Value(Degree, x);

    /// <inheritdoc/>
    public override double Derivative(double x) => ShiftedLegendre.Derivative(Degree, x);

    /// <inheritdoc/>
    public override double Mean(ProductMeasure measure, out double error)
    {
        ContinuousMeasure m = measure.Continuous(Variable)!;
        (double mean, error) = Product(Degree, 0, 0.0, 0.0).Expectation(m.Measure, m.Lower, m.Upper);
        return mean;
    }

    /// <inheritdoc/>
    public override (double Cross, double Covariance) Pair(DesignFactor other, double meanThis, double meanOther,
                                                           ProductMeasure measure, out double error)
    {
        if (other is not LegendreDesignFactor b || b.Variable != Variable) throw new ArgumentException("A Legendre factor pairs with a Legendre factor on its variable.", nameof(other));
        ContinuousMeasure m = measure.Continuous(Variable)!;
        (double cross, double ec) = Product(Degree, b.Degree, 0.0, 0.0).Expectation(m.Measure, m.Lower, m.Upper);
        (double covariance, double ev) = Product(Degree, b.Degree, meanThis, meanOther).Expectation(m.Measure, m.Lower, m.Upper);
        error = Math.Max(ec, ev);
        return (cross, covariance);
    }

    // (ψ_a − m_a)(ψ_b − m_b) as one piece, b = 0 standing for the constant 1 with no offset: a single factor's mean.
    private static PiecewisePolynomial Product(int a, int b, double offsetA, double offsetB)
        => new([], a + b, (piece, centre, coefficients) =>
        {
            Span<double> ca = stackalloc double[a + 1], cb = stackalloc double[b + 1];
            ShiftedLegendre.Coefficients(a, centre, ca);
            ShiftedLegendre.Coefficients(b, centre, cb);
            ca[0] -= offsetA;
            cb[0] -= offsetB;
            for (int i = 0; i <= a; i++)
                for (int j = 0; j <= b; j++) coefficients[i + j] += ca[i] * cb[j];
            return true;
        });
}

/// <summary>The indicator 1(z ∈ D) of a categorical subset as a factor, on a level index.</summary>
public sealed record SubsetDesignFactor(int Variable, CategoricalFactor Subset) : DesignFactor(Variable)
{
    /// <inheritdoc/>
    public override double Value(double x) => Subset.Contains((int)x) ? 1.0 : 0.0;

    /// <inheritdoc/>
    public override double Derivative(double x) => throw new NotSupportedException("A level index has no derivative.");

    /// <inheritdoc/>
    public override double Mean(ProductMeasure measure, out double error)
    {
        error = 0.0;
        IReadOnlyList<double> p = measure.Probabilities(Variable)!;
        double sum = 0.0;
        foreach (int level in Subset.Levels) sum += p[level];
        return sum;
    }

    /// <inheritdoc/>
    public override (double Cross, double Covariance) Pair(DesignFactor other, double meanThis, double meanOther,
                                                           ProductMeasure measure, out double error)
    {
        if (other is not SubsetDesignFactor b || b.Variable != Variable) throw new ArgumentException("A subset pairs with a subset on its variable.", nameof(other));
        error = 0.0;
        IReadOnlyList<double> p = measure.Probabilities(Variable)!;
        double both = 0.0;
        foreach (int level in Subset.Levels)
            if (b.Subset.Contains(level)) both += p[level];
        return (both, both - meanThis * meanOther);
    }
}

/// <summary>R's hinge arithmetic: the value of (s(x − t))₊^q, its derivative, and the peak that normalizes a product.</summary>
internal static class Hinge
{
    // R's x^q: x at q = 1, x·x at q = 2, the library power above.
    public static double Power(int degree, double x) => degree switch
    {
        1 => x,
        2 => x * x,
        _ => Math.Pow(x, degree),
    };

    public static double Value(int degree, int sign, double knot, double x) => Power(degree, Math.Max(sign * (x - knot), 0.0));

    // column[i] *= Value(degree, sign, knot, x[i]), then column[i] /= divisor when `divide`: the scalar form's operations per element,
    // in its order, four lanes at a time for degrees one and two. The lanes' subtract, multiply, compare, select and divide are exact
    // per lane, so the bits are the scalar loop's on any hardware, and Math.Max's rules are kept: a NaN stays, a negative or −0 becomes +0.
    // The library power (degree three and up) and the tail stay scalar.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static void Scale(Span<double> column, ReadOnlySpan<double> x, int degree, int sign, double knot, bool divide, double divisor)
    {
        int n = column.Length, i = 0;
        if (x.Length != n) throw new ArgumentException("One value per row.", nameof(x));
        if (degree is 1 or 2 && n >= Vector256<double>.Count)
        {
            ref double c = ref MemoryMarshal.GetReference(column);
            ref double xr = ref MemoryMarshal.GetReference(x);
            Vector256<double> vk = Vector256.Create(knot), vs = Vector256.Create((double)sign), vd = Vector256.Create(divisor);
            Vector256<double> zero = Vector256<double>.Zero;
            for (; i + Vector256<double>.Count <= n; i += Vector256<double>.Count)
            {
                Vector256<double> u = Vector256.Multiply(vs, Vector256.Subtract(Vector256.LoadUnsafe(ref xr, (nuint)i), vk));
                Vector256<double> keep = Vector256.GreaterThan(u, zero) | ~Vector256.Equals(u, u);
                Vector256<double> h = Vector256.ConditionalSelect(keep, u, zero);
                if (degree == 2) h = Vector256.Multiply(h, h);
                Vector256<double> r = Vector256.Multiply(Vector256.LoadUnsafe(ref c, (nuint)i), h);
                if (divide) r = Vector256.Divide(r, vd);
                r.StoreUnsafe(ref c, (nuint)i);
            }
        }
        for (; i < n; i++)
        {
            double v = column[i] * Value(degree, sign, knot, x[i]);
            column[i] = divide ? v / divisor : v;
        }
    }

    // d/dx (s(x − t))₊^q = q·s·(s(x − t))₊^(q−1). At the knot, where degree one jumps from 0 to s, the mean of the
    // one-sided limits, s/2; from degree two up the hinge is differentiable there, with derivative 0.
    public static double Derivative(int degree, int sign, double knot, double x)
    {
        double u = sign * (x - knot);
        if (degree == 1) return u > 0.0 ? sign : u == 0.0 ? 0.5 * sign : 0.0;
        return u > 0.0 ? degree * sign * Power(degree - 1, u) : 0.0;
    }

    // 1/c with c the product's largest value on [0, 1], (Πⱼ ((sⱼ + 1)/2 − sⱼtⱼ))^q; 1 where c is 0, a hinge facing its end
    // of the range, which is zero there and which R leaves at c = 1.
    public static double Normalizer(int degree, double peak)
    {
        double c = Power(degree, peak);
        return c == 0.0 ? 1.0 : 1.0 / c;
    }
}
