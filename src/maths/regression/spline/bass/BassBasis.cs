using System;
using Maths.LinAlg;

namespace Maths.Regression.Spline.Bass;

/// <summary>A term's support on the training inputs, by part: its nonzero values and its mass (D01, D06).</summary>
public readonly record struct TermSupport(int CaseNonzero, double CaseMass, int OutputNonzero, double OutputMass);

/// <summary>
/// A fit's basis: the training inputs' knots and level counts, and the hinges' degree q. It checks terms against
/// them and evaluates each term as two factors, R's construction (<c>makeBasis</c>, <c>makeBasisCat</c>):
/// <list type="bullet">
/// <item><b>The case side</b>, a(x, z) = Πⱼ (sⱼ(xⱼ − tⱼ))₊^q / c · Πₗ 1(zₗ ∈ Dₗ), with
/// c = (Πⱼ ((sⱼ + 1)/2 − sⱼtⱼ))^q the hinge product's largest value over [0, 1], and c = 1 where that is 0 (a knot at
/// the end its hinge faces, which leaves the hinge zero on the training range). So a continuous part peaks at one, and
/// D06's column mass reads hinges normalized to a maximum of one.</item>
/// <item><b>The output side</b>, h(t), the functional part's hinges normalized the same way.</item>
/// </list>
/// An empty part is the constant 1. A scalar response's column is the case side; a direct functional response's
/// column is the Khatri–Rao product, a(xᵢ)·h(t_r) at case i and location r, which the design layer contracts without
/// forming (the hub's direct functional algebra). Powers follow R's <c>^</c>: x·x at q = 2, the library power above.
/// <para>Each term type is a <see cref="IDesignKernel{TPoints}"/> (<see cref="Kernels"/>); the basis applies them in order to
/// a column of ones, so a column, its support and its derivatives all come through the one contract.</para>
/// </summary>
public sealed class BassBasis
{
    /// <summary>The basis of <paramref name="inputs"/> at hinge degree <paramref name="degree"/>.</summary>
    public BassBasis(BassInputs inputs, int degree)
    {
        ArgumentNullException.ThrowIfNull(inputs);
        if (degree < 1) throw new ArgumentOutOfRangeException(nameof(degree), "The degree is a positive integer.");
        Inputs = inputs;
        Degree = degree;
        Kernels = new DesignKernels(inputs, degree);
    }

    /// <summary>The training inputs.</summary>
    public BassInputs Inputs { get; }

    /// <summary>The hinges' degree q.</summary>
    public int Degree { get; }

    /// <summary>The term types' kernels (F-X3): a term's columns, derivatives and factors are these applied in order.</summary>
    public DesignKernels Kernels { get; }

    /// <summary>
    /// Refuses a term whose indices fall outside the inputs: a variable or knot beyond its kind's, a level beyond its
    /// predictor's, a subset that holds every level (not a proper subset, so its indicator is the intercept's), or a
    /// functional factor when the response is scalar.
    /// </summary>
    public void Validate(BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(term);
        foreach (IDesignKernel<CaseInputs> kernel in Kernels.Case) kernel.Validate(term);
        foreach (IDesignKernel<OutputInputs> kernel in Kernels.Output) kernel.Validate(term);
    }

    /// <summary>The term's case-side factor at each case of <paramref name="at"/>, into <paramref name="column"/>.</summary>
    public void CaseColumn(BassTerm term, CaseInputs at, Span<double> column)
    {
        ArgumentNullException.ThrowIfNull(term);
        ArgumentNullException.ThrowIfNull(at);
        if (column.Length != at.Cases) throw new ArgumentException("The column needs one value per case.", nameof(column));
        column.Fill(1.0);
        foreach (IDesignKernel<CaseInputs> kernel in Kernels.Case) kernel.Multiply(term, at, column);
    }

    /// <summary>
    /// The partial derivative of the term's case-side factor in continuous predictor <paramref name="variable"/> (by its
    /// position among the continuous predictors) at each case of <paramref name="at"/>, into <paramref name="column"/>: the
    /// derivative of the part that reads the variable (<see cref="IDesignKernel{TPoints}.ReadsVariable"/>) times the other
    /// parts' values, zero for a term that does not hold the variable. At a hinge's knot the derivative is the mean of its
    /// one-sided limits (<see cref="IDesignKernel{TPoints}"/>).
    /// </summary>
    public void CaseDerivative(BassTerm term, int variable, CaseInputs at, Span<double> column)
    {
        ArgumentNullException.ThrowIfNull(term);
        ArgumentNullException.ThrowIfNull(at);
        if (column.Length != at.Cases) throw new ArgumentException("The column needs one value per case.", nameof(column));
        if ((uint)variable >= (uint)at.ContinuousCount) throw new ArgumentOutOfRangeException(nameof(variable), "No such continuous variable.");
        column.Fill(1.0);
        bool read = false;
        foreach (IDesignKernel<CaseInputs> kernel in Kernels.Case)
        {
            if (kernel.ReadsVariable(term, variable))
            {
                kernel.MultiplyDerivative(term, variable, at, column);
                read = true;
            }
            else kernel.Multiply(term, at, column);
        }
        if (!read) column.Clear();
    }

    /// <summary>The term's output-side factor at each location of <paramref name="at"/>, into <paramref name="column"/>.</summary>
    public void OutputColumn(BassTerm term, OutputInputs at, Span<double> column)
    {
        ArgumentNullException.ThrowIfNull(term);
        ArgumentNullException.ThrowIfNull(at);
        if (column.Length != at.Locations) throw new ArgumentException("The column needs one value per location.", nameof(column));
        column.Fill(1.0);
        foreach (IDesignKernel<OutputInputs> kernel in Kernels.Output) kernel.Multiply(term, at, column);
    }

    /// <summary>The partial derivative of the term's output-side factor in output variable <paramref name="variable"/> at each
    /// location of <paramref name="at"/>, into <paramref name="column"/> (<see cref="CaseDerivative"/>'s counterpart).</summary>
    public void OutputDerivative(BassTerm term, int variable, OutputInputs at, Span<double> column)
    {
        ArgumentNullException.ThrowIfNull(term);
        ArgumentNullException.ThrowIfNull(at);
        if (column.Length != at.Locations) throw new ArgumentException("The column needs one value per location.", nameof(column));
        if ((uint)variable >= (uint)at.Count) throw new ArgumentOutOfRangeException(nameof(variable), "No such output variable.");
        column.Fill(1.0);
        foreach (IDesignKernel<OutputInputs> kernel in Kernels.Output) kernel.MultiplyDerivative(term, variable, at, column);
    }

    /// <summary>The term's support on the training inputs: each part's nonzero count and mass.</summary>
    public TermSupport Support(BassTerm term)
    {
        Validate(term);
        var a = new double[Inputs.Cases.Cases];
        CaseColumn(term, Inputs.Cases, a);
        (int caseNonzero, double caseMass) = Count(a);
        if (Inputs.Outputs is not { } outputs) return new TermSupport(caseNonzero, caseMass, 0, 0.0);
        var h = new double[outputs.Locations];
        OutputColumn(term, outputs, h);
        (int outputNonzero, double outputMass) = Count(h);
        return new TermSupport(caseNonzero, caseMass, outputNonzero, outputMass);
    }

    private static (int Nonzero, double Mass) Count(double[] column)
    {
        int nonzero = 0;
        foreach (double v in column)
            if (v != 0.0) nonzero++;
        return (nonzero, PairwiseSum.Of(column));
    }
}
