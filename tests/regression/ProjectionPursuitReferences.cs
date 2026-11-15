using System;
using System.Collections.Generic;
using System.Linq;

namespace Maths.Regression.Tests;

/// <summary>
/// Friedman's function f(x) = 10 sin(πx₁x₂) + 20(x₃ − ½)² + 10x₄ + 5x₅ (the BPPR paper's §3.1, after Friedman 1991), its inputs past the
/// fifth inert. ∂₄f : ∂₅f is 2 : 1 everywhere, so 10x₄ + 5x₅ is one ridge on (2, 1)/√5 and C has rank 4. Under the uniform law on
/// [0, 1]^p, with c = cos(πx₁x₂): E[x₂c] = 2/π², E[x₂²c²] = 1/6 − 1/(8π²), E[x₁x₂c²] = 1/8 − Cin(2π)/(8π²) with
/// Cin(z) = ∫₀^z (1 − cos t)/t dt; so C₁₁ = C₂₂ = 50π²/3 − 25/2, C₁₂ = 25π²/2 − 25·Cin(2π)/2, C₁₄ = C₂₄ = 200/π, C₁₅ = C₂₅ = 100/π,
/// C₃₃ = 400/3, C₄₄ = 100, C₄₅ = 50, C₅₅ = 25, and the rest 0. Other laws have no closed form here.
/// </summary>
internal sealed class FriedmanReference : ReferenceFunction
{
    public FriedmanReference(int dimension) : base("Friedman's function", dimension)
    {
        if (dimension < 5) throw new ArgumentOutOfRangeException(nameof(dimension), "Friedman's function reads five inputs.");
    }

    public override double Value(ReadOnlySpan<double> x) =>
        10.0 * Math.Sin(Math.PI * x[0] * x[1]) + 20.0 * (x[2] - 0.5) * (x[2] - 0.5) + 10.0 * x[3] + 5.0 * x[4];

    public override void Gradient(ReadOnlySpan<double> x, Span<double> gradient)
    {
        gradient.Clear();
        double c = 10.0 * Math.PI * Math.Cos(Math.PI * x[0] * x[1]);
        gradient[0] = c * x[1];
        gradient[1] = c * x[0];
        gradient[2] = 40.0 * (x[2] - 0.5);
        gradient[3] = 10.0;
        gradient[4] = 5.0;
    }

    public override double[,]? C(ReferenceMeasure measure)
    {
        if (measure.Dimension != Dimension) throw new ArgumentException("The law has other inputs than the function.");
        if (measure is not ReferenceBox box || box.Lower.Any(a => a != 0.0) || box.Upper.Any(b => b != 1.0)) return null;
        double pi2 = Math.PI * Math.PI;
        var c = new double[Dimension, Dimension];
        c[0, 0] = c[1, 1] = 50.0 * pi2 / 3.0 - 12.5;
        c[0, 1] = c[1, 0] = 12.5 * pi2 - 12.5 * CinTwoPi();
        c[0, 3] = c[3, 0] = c[1, 3] = c[3, 1] = 200.0 / Math.PI;
        c[0, 4] = c[4, 0] = c[1, 4] = c[4, 1] = 100.0 / Math.PI;
        c[2, 2] = 400.0 / 3.0;
        c[3, 3] = 100.0;
        c[3, 4] = c[4, 3] = 50.0;
        c[4, 4] = 25.0;
        return c;
    }

    /// <summary>Cin(2π) = ∫₀^(2π) (1 − cos t)/t dt = Σₖ (−1)^(k+1) (2π)^(2k)/(2k·(2k)!), summed in decimal: its terms reach about 16 against
    /// a sum of about 2.4, which would cost a double's sum a digit.</summary>
    public static double CinTwoPi()
    {
        const decimal TwoPi = 6.2831853071795864769252867666m;
        decimal z2 = TwoPi * TwoPi, term = z2 / 2m, sum = 0m;   // term = (2π)^(2k)/(2k)!
        for (int k = 1; term / (2 * k) > 1e-27m; k++)
        {
            sum += (k % 2 == 1 ? term : -term) / (2 * k);
            term *= z2 / ((2 * k + 1) * (2 * k + 2));
        }
        return (double)sum;
    }
}

/// <summary>
/// A reference for E36 (the BPPR brief, slice 0): a function from outside the ridge model, the law its C is under, and unit directions
/// spanning C's range, the truth that slice 2's recovery reads principal angles against. The directions are in the inputs' own units;
/// under a metric D, such as the fit's standardization (decision 3), the direction w reads as Dw renormalized.
/// </summary>
internal sealed record ProjectionPursuitReference(ReferenceFunction Function, ReferenceBox Law, double[][] Directions)
{
    public string Name => Function.Name;

    public override string ToString() => Name;
}

/// <summary>
/// E36's references: single indices with sparse, dense and box-skewed directions, a sum of two ridges sharing an input, and Friedman's
/// function as the paper runs it. None is drawn from BPPR's prior: the profiles are exponentials and waves, not one-sided splines.
/// </summary>
internal static class ProjectionPursuitReferences
{
    /// <summary>exp(0.6x₁ + 0.8x₂) on [0, 1]⁴, two inputs inert: two active coordinates, within a ridge's three.</summary>
    public static ProjectionPursuitReference SingleIndex() => new(
        new RidgeSum("single index", new RidgeTerm(ExponentialProfile.Exp(1.0), [0.6, 0.8, 0.0, 0.0])),
        ReferenceBox.Cube(4, 0.0, 1.0),
        [[0.6, 0.8, 0.0, 0.0]]);

    /// <summary>sin(π(x₁ − x₂ + ½x₃ + ½x₄)) on [0, 1]⁵, x₅ inert: four active coordinates, more than one ridge of BPPR's prior holds.</summary>
    public static ProjectionPursuitReference DenseSingleIndex()
    {
        double[] w = [1.0, -1.0, 0.5, 0.5, 0.0];
        return new(new RidgeSum("dense single index", new RidgeTerm(ExponentialProfile.Sin(Math.PI), w)),
                   ReferenceBox.Cube(5, 0.0, 1.0), [Normalized(w)]);
    }

    /// <summary>cos(2(x₁ + x₃)) on [0, 1] × [0, 4] × [−1, 1], x₂ inert: the sides differ, so standardizing turns (1, 0, 1)/√2 to
    /// (1, 0, 2)/√5.</summary>
    public static ProjectionPursuitReference SkewedSingleIndex()
    {
        double[] w = [1.0, 0.0, 1.0];
        return new(new RidgeSum("single index on a skewed box", new RidgeTerm(ExponentialProfile.Cos(2.0), w)),
                   new ReferenceBox([0.0, 0.0, -1.0], [1.0, 4.0, 1.0]), [Normalized(w)]);
    }

    /// <summary>sin(3(0.8x₁ + 0.6x₂)) + exp(1.2(0.6x₂ − 0.8x₃ − 0.4)) on [0, 1]⁵, two inputs inert: the directions share x₂, at cosine
    /// 0.36.</summary>
    public static ProjectionPursuitReference TwoRidges()
    {
        double[] w1 = [0.8, 0.6, 0.0, 0.0, 0.0], w2 = [0.0, 0.6, -0.8, 0.0, 0.0];
        return new(new RidgeSum("two ridges", new RidgeTerm(ExponentialProfile.Sin(3.0), w1), new RidgeTerm(ExponentialProfile.Exp(1.2), w2, 0.4)),
                   ReferenceBox.Cube(5, 0.0, 1.0), [w1, w2]);
    }

    /// <summary>Friedman's function on [0, 1]⁶, as the paper's §3.1 runs it: C's range is e₁, e₂, e₃ and (2e₄ + e₅)/√5.</summary>
    public static ProjectionPursuitReference Friedman() => new(
        new FriedmanReference(6),
        ReferenceBox.Cube(6, 0.0, 1.0),
        [[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], Normalized([0, 0, 0, 2, 1, 0])]);

    public static IReadOnlyList<ProjectionPursuitReference> All() => [SingleIndex(), DenseSingleIndex(), SkewedSingleIndex(), TwoRidges(), Friedman()];

    private static double[] Normalized(double[] w)
    {
        double norm = Math.Sqrt(w.Sum(v => v * v));
        return [.. w.Select(v => v / norm)];
    }
}
