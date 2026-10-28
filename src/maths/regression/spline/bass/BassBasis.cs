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
    }

    /// <summary>The training inputs.</summary>
    public BassInputs Inputs { get; }

    /// <summary>The hinges' degree q.</summary>
    public int Degree { get; }

    /// <summary>
    /// Refuses a term whose indices fall outside the inputs: a variable or knot beyond its kind's, a level beyond its
    /// predictor's, a subset that holds every level (not a proper subset, so its indicator is the intercept's), or a
    /// functional factor when the response is scalar.
    /// </summary>
    public void Validate(BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(term);
        CaseInputs cases = Inputs.Cases;
        foreach (HingeFactor h in term.Continuous)
        {
            if (h.Variable >= cases.ContinuousCount) throw new ArgumentException($"Continuous variable {h.Variable} does not exist.", nameof(term));
            if (h.Knot >= Inputs.ContinuousKnots(h.Variable).Length) throw new ArgumentException($"Knot {h.Knot} of continuous variable {h.Variable} does not exist.", nameof(term));
        }
        foreach (CategoricalFactor c in term.Categorical)
        {
            if (c.Variable >= cases.CategoricalCount) throw new ArgumentException($"Categorical variable {c.Variable} does not exist.", nameof(term));
            int levels = Inputs.LevelCount(c.Variable);
            if (c.Levels[^1] >= levels) throw new ArgumentException($"Categorical variable {c.Variable} has {levels} levels.", nameof(term));
            if (c.Levels.Count >= levels) throw new ArgumentException($"A subset of categorical variable {c.Variable} must leave a level out.", nameof(term));
        }
        if (term.Functional.Count == 0) return;
        OutputInputs outputs = Inputs.Outputs ?? throw new ArgumentException("A scalar response has no functional factors.", nameof(term));
        foreach (HingeFactor h in term.Functional)
        {
            if (h.Variable >= outputs.Count) throw new ArgumentException($"Output variable {h.Variable} does not exist.", nameof(term));
            if (h.Knot >= Inputs.FunctionalKnots(h.Variable).Length) throw new ArgumentException($"Knot {h.Knot} of output variable {h.Variable} does not exist.", nameof(term));
        }
    }

    /// <summary>The term's case-side factor at each case of <paramref name="at"/>, into <paramref name="column"/>.</summary>
    public void CaseColumn(BassTerm term, CaseInputs at, Span<double> column)
    {
        ArgumentNullException.ThrowIfNull(term);
        ArgumentNullException.ThrowIfNull(at);
        if (column.Length != at.Cases) throw new ArgumentException("The column needs one value per case.", nameof(column));
        if (term.Continuous.Count == 0) column.Fill(1.0);
        else Hinges(term.Continuous, at.Continuous, Inputs.ContinuousKnots, column);
        foreach (CategoricalFactor c in term.Categorical)
        {
            ReadOnlySpan<int> levels = at.Categorical(c.Variable);
            for (int i = 0; i < column.Length; i++)
                if (!c.Contains(levels[i])) column[i] = 0.0;
        }
    }

    /// <summary>The term's output-side factor at each location of <paramref name="at"/>, into <paramref name="column"/>.</summary>
    public void OutputColumn(BassTerm term, OutputInputs at, Span<double> column)
    {
        ArgumentNullException.ThrowIfNull(term);
        ArgumentNullException.ThrowIfNull(at);
        if (column.Length != at.Locations) throw new ArgumentException("The column needs one value per location.", nameof(column));
        if (term.Functional.Count == 0) column.Fill(1.0);
        else Hinges(term.Functional, at.Variable, Inputs.FunctionalKnots, column);
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

    private delegate ReadOnlySpan<double> Column(int variable);

    // Πⱼ (sⱼ(xⱼ − tⱼ))₊^q over the factors in variable order, divided by the product's largest value on [0, 1].
    private void Hinges(System.Collections.Generic.IReadOnlyList<HingeFactor> factors, Column values, Column knots, Span<double> column)
    {
        column.Fill(1.0);
        double peak = 1.0;
        foreach (HingeFactor f in factors)
        {
            double t = knots(f.Variable)[f.Knot], s = f.Sign;
            peak *= (s + 1.0) / 2.0 - s * t;
            ReadOnlySpan<double> x = values(f.Variable);
            for (int i = 0; i < column.Length; i++) column[i] *= Power(Math.Max(s * (x[i] - t), 0.0));
        }
        double c = Power(peak);
        if (c == 0.0) return;   // a hinge facing its end of the range: zero there, so R leaves c at 1
        for (int i = 0; i < column.Length; i++) column[i] /= c;
    }

    // R's x^q: x at q = 1, x·x at q = 2, the library power above.
    private double Power(double x) => Degree switch
    {
        1 => x,
        2 => x * x,
        _ => Math.Pow(x, Degree),
    };
}
