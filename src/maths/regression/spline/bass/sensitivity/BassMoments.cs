using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass.Sensitivity;

/// <summary>
/// One subset's share of a term pair's covariance: Cov(B_a,u, B_b,u) for the terms' ANOVA components on the
/// <paramref name="Variables"/> u (ascending, in <see cref="ProductMeasure"/>'s numbering).
/// </summary>
public readonly record struct SubsetCovariance(int[] Variables, double Value);

/// <summary>One variable's share of a term pair's total-effect covariance, E[Cov(B_a, B_b | every other variable)].</summary>
public readonly record struct VariableCovariance(int Variable, double Value);

/// <summary>
/// Two terms' moments under a product measure and their ANOVA structure. A term is a product of one-variable factors,
/// B = κ·Πᵥ φᵥ(xᵥ), so its ANOVA component on a subset u is κ·Π_{v∈u}(φᵥ − E φᵥ)·Π_{v∉u} E φᵥ, nonzero only for u inside
/// its variables, and two terms' components on u covary only for u inside the variables they share:
/// Cov_u = κ_aκ_b·Π_{v∈u} Cov(φ_av, φ_bv)·Π_{v∈V_a∖u} E φ_av·Π_{v∈V_b∖u} E φ_bv.
/// </summary>
public sealed class TermPairMoments
{
    internal TermPairMoments(double cross, double product, SubsetCovariance[] components, VariableCovariance[] totals)
    {
        Cross = cross;
        Product = product;
        Components = Array.AsReadOnly(components);
        Totals = Array.AsReadOnly(totals);
        double covariance = 0.0;
        foreach (SubsetCovariance c in components) covariance += c.Value;
        Covariance = covariance;
    }

    /// <summary>E[B_a·B_b].</summary>
    public double Cross { get; }

    /// <summary>E[B_a]·E[B_b].</summary>
    public double Product { get; }

    /// <summary>Cov(B_a, B_b), the sum of <see cref="Components"/>: Cross − Product, without forming that difference.</summary>
    public double Covariance { get; }

    /// <summary>The covariance's share on every nonempty subset of the shared variables, smaller subsets first.</summary>
    public IReadOnlyList<SubsetCovariance> Components { get; }

    /// <summary>Each shared variable's total-effect share: the components on the subsets that hold it, summed.</summary>
    public IReadOnlyList<VariableCovariance> Totals { get; }
}

/// <summary>
/// BASS's moment service (D08, O5): its terms' means E[B], cross moments E[B_a B_b] and their ANOVA structure under a
/// <see cref="ProductMeasure"/>, for direct sensitivity and, unchanged, basis sensitivity's cross-component pairs. A term's
/// moment is a product over its variables of one-variable moments: a hinge's or hinge product's by
/// <see cref="PiecewisePolynomial"/>, split at the knots and clipped to the measure's support, at every degree; a
/// categorical factor's as a sum of the supplied level probabilities, an intersection's included.
/// <para>The factors come from the term types' design kernels (<see cref="IDesignKernel{TPoints}"/>), which supply each factor's
/// value, derivative (<see cref="FactorDerivative"/>) and moments; this service multiplies and caches them.</para>
/// <para>Terms are taken by value, so terms from different fits on the same inputs (a basis fit's components) meet here;
/// their knots and levels index <see cref="Inputs"/>. The normalization is the basis's: the case side's hinges divided by
/// their largest value on [0, 1], and the output side's likewise when the measure covers every output variable. Where it
/// leaves some out, a term's moments are over the factors it covers, and the caller holds the others at a coordinate, with
/// the output side's normalizer, which is one constant over the whole part.</para>
/// <para>Every factor, factor pair, term and term pair is computed once and cached, since draws repeat their structures.
/// One caller at a time.</para>
/// </summary>
public sealed class BassMoments
{
    private readonly Dictionary<DesignFactor, double> _means = [];
    private readonly Dictionary<(DesignFactor, DesignFactor), (double Cross, double Covariance)> _factorPairs = [];
    private readonly Dictionary<BassTerm, TermShape> _shapes = [];
    private readonly Dictionary<(BassTerm, BassTerm), TermPairMoments> _pairs = [];
    private readonly DesignKernels _kernels;

    private BassMoments(BassInputs inputs, int degree, ProductMeasure measure)
    {
        Inputs = inputs;
        Degree = degree;
        Measure = measure;
        _kernels = new DesignKernels(inputs, degree);
    }

    /// <summary>The service for terms over <paramref name="inputs"/> at hinge degree <paramref name="degree"/>.</summary>
    public static BassMoments On(BassInputs inputs, int degree, ProductMeasure measure)
    {
        ArgumentNullException.ThrowIfNull(inputs);
        ArgumentNullException.ThrowIfNull(measure);
        if (degree < 1) throw new ArgumentOutOfRangeException(nameof(degree), "The degree is a positive integer.");
        return new BassMoments(inputs, degree, measure);
    }

    /// <summary>The inputs the terms' knots and levels index.</summary>
    public BassInputs Inputs { get; }

    /// <summary>The hinges' degree q.</summary>
    public int Degree { get; }

    /// <summary>The product measure.</summary>
    public ProductMeasure Measure { get; }

    /// <summary>The largest quadrature error estimate among the one-variable moments so far: 0 unless a Student component's
    /// Gauss–Legendre rule was used.</summary>
    public double LargestError { get; private set; }

    /// <summary>E[B] for <paramref name="term"/>.</summary>
    public double Mean(BassTerm term)
    {
        TermShape s = Shape(term);
        double mean = s.Scale;
        foreach (double m in s.Means) mean *= m;
        return mean;
    }

    /// <summary>E[B_a·B_b].</summary>
    public double Cross(BassTerm a, BassTerm b) => Pair(a, b).Cross;

    /// <summary>The pair's moments and ANOVA structure.</summary>
    public TermPairMoments Pair(BassTerm a, BassTerm b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        if (_pairs.TryGetValue((a, b), out TermPairMoments? pair) || _pairs.TryGetValue((b, a), out pair)) return pair;
        pair = Build(Shape(a), Shape(b));
        _pairs[(a, b)] = pair;
        return pair;
    }

    /// <summary><paramref name="term"/>'s variables (ascending) in the measure's numbering.</summary>
    public IReadOnlyList<int> Variables(BassTerm term) => Shape(term).Variables;

    /// <summary>The scale κ that normalizes <paramref name="term"/>'s hinges, as the basis does.</summary>
    public double Scale(BassTerm term) => Shape(term).Scale;

    /// <summary>E[φᵥ] for each of <paramref name="term"/>'s factors, in <see cref="Variables"/>' order.</summary>
    public IReadOnlyList<double> FactorMeans(BassTerm term) => Shape(term).Means;

    /// <summary>
    /// <paramref name="term"/>'s factor on its <paramref name="position"/>-th variable at <paramref name="value"/>, unnormalized:
    /// (s(x − t))₊^q for a hinge, at a scaled value; 1(z ∈ D) for a categorical factor, at a level index.
    /// </summary>
    public double FactorValue(BassTerm term, int position, double value) => Shape(term).Factors[position].Value(value);

    /// <summary>
    /// The derivative of <see cref="FactorValue"/> in its variable (O-X4, through the design-kernel contract): q·s·(s(x − t))₊^(q−1)
    /// for a hinge, the mean of the one-sided limits at the knot, which is s/2 at q = 1 and 0 above. Refused for a categorical
    /// factor, whose variable is a level.
    /// </summary>
    public double FactorDerivative(BassTerm term, int position, double value) => Shape(term).Factors[position].Derivative(value);

    private TermPairMoments Build(TermShape a, TermShape b)
    {
        // Per variable of either term: E[φ_a φ_b] (a factor alone where the other term lacks the variable), and the
        // shared variables' covariances.
        var union = new List<int>();
        var cross = new List<double>();
        var shared = new List<int>();
        var covariance = new List<double>();
        var meanA = new List<double>();
        var meanB = new List<double>();
        int i = 0, j = 0;
        double restA = 1.0, restB = 1.0;   // Π of the unshared variables' means
        while (i < a.Variables.Length || j < b.Variables.Length)
        {
            int va = i < a.Variables.Length ? a.Variables[i] : int.MaxValue, vb = j < b.Variables.Length ? b.Variables[j] : int.MaxValue;
            if (va < vb)
            {
                union.Add(va);
                cross.Add(a.Means[i]);
                restA *= a.Means[i++];
            }
            else if (vb < va)
            {
                union.Add(vb);
                cross.Add(b.Means[j]);
                restB *= b.Means[j++];
            }
            else
            {
                (double c, double cov) = FactorPair(a.Factors[i], b.Factors[j], a.Means[i], b.Means[j]);
                union.Add(va);
                cross.Add(c);
                shared.Add(va);
                covariance.Add(cov);
                meanA.Add(a.Means[i++]);
                meanB.Add(b.Means[j++]);
            }
        }
        double scale = a.Scale * b.Scale;
        double crossMoment = scale, product = scale * restA * restB;
        foreach (double c in cross) crossMoment *= c;
        for (int k = 0; k < meanA.Count; k++) product *= meanA[k] * meanB[k];

        // The components: every nonempty subset u of the shared variables, Π_{u} Cov·Π_{shared∖u} E_a E_b·rest.
        int n = shared.Count;
        var components = new List<SubsetCovariance>((1 << n) - 1);
        for (int size = 1; size <= n; size++)
            for (int mask = 1; mask < 1 << n; mask++)
            {
                if (System.Numerics.BitOperations.PopCount((uint)mask) != size) continue;
                double value = scale * restA * restB;
                var subset = new int[size];
                int at = 0;
                for (int k = 0; k < n; k++)
                {
                    if ((mask >> k & 1) != 0)
                    {
                        value *= covariance[k];
                        subset[at++] = shared[k];
                    }
                    else value *= meanA[k] * meanB[k];
                }
                components.Add(new SubsetCovariance(subset, value));
            }
        // The totals: for shared v, Cov_v·Π over the union's other variables of E[φ_a φ_b].
        var totals = new VariableCovariance[n];
        for (int k = 0; k < n; k++)
        {
            double value = scale * covariance[k];
            for (int u = 0; u < union.Count; u++)
                if (union[u] != shared[k]) value *= cross[u];
            totals[k] = new VariableCovariance(shared[k], value);
        }
        return new TermPairMoments(crossMoment, product, [.. components], totals);
    }

    // (E[φ_a φ_b], Cov(φ_a, φ_b)) for two factors on one variable.
    private (double Cross, double Covariance) FactorPair(DesignFactor a, DesignFactor b, double meanA, double meanB)
    {
        if (_factorPairs.TryGetValue((a, b), out var cached)) return cached;
        (double cross, double covariance) = a.Pair(b, meanA, meanB, Measure, out double error);
        LargestError = Math.Max(LargestError, error);
        _factorPairs[(a, b)] = (cross, covariance);
        _factorPairs[(b, a)] = (cross, covariance);
        return (cross, covariance);
    }

    private TermShape Shape(BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(term);
        if (_shapes.TryGetValue(term, out TermShape? shape)) return shape;
        var factors = new List<DesignFactor>();
        double scale = 1.0;
        foreach (IDesignKernel<CaseInputs> kernel in _kernels.Case)
        {
            factors.AddRange(kernel.Factors(term, Measure));
            scale *= kernel.Scale(term);
        }
        foreach (IDesignKernel<OutputInputs> kernel in _kernels.Output)
        {
            if (!kernel.Reads(term)) continue;
            factors.AddRange(kernel.Factors(term, Measure));
            if (Measure.IncludesOutputs) scale *= kernel.Scale(term);
        }
        factors.Sort((x, y) => x.Variable.CompareTo(y.Variable));
        var variables = new int[factors.Count];
        var means = new double[factors.Count];
        var each = new DesignFactor[factors.Count];
        for (int k = 0; k < factors.Count; k++)
        {
            each[k] = factors[k];
            variables[k] = factors[k].Variable;
            means[k] = FactorMean(factors[k]);
        }
        shape = new TermShape(variables, each, means, scale);
        _shapes[term] = shape;
        return shape;
    }

    private double FactorMean(DesignFactor f)
    {
        if (_means.TryGetValue(f, out double mean)) return mean;
        mean = f.Mean(Measure, out double error);
        LargestError = Math.Max(LargestError, error);
        _means[f] = mean;
        return mean;
    }

    private sealed record TermShape(int[] Variables, DesignFactor[] Factors, double[] Means, double Scale);
}
