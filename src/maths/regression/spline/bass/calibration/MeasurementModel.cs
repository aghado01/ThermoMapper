using System;
using System.Collections.Generic;
using System.Text.Json.Serialization;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass.Calibration;

/// <summary>
/// The declared measurement-error model ε (D09, C-X1): a closed set, so robust models (GBMARS's seam) join as further
/// kinds with their own build. Gaussian with variance groups now.
/// </summary>
[JsonPolymorphic(TypeDiscriminatorPropertyName = "Kind")]
[JsonDerivedType(typeof(GaussianMeasurement), "Gaussian")]
public abstract record MeasurementErrorModel
{
    private protected MeasurementErrorModel() { }
}

/// <summary>
/// Gaussian measurement error with variance groups: ε_j ~ N(0, σ²_{g(j)}), each output j in exactly one group, declared
/// by output index (R's <c>s2.ind</c>). Under the multivariate formulation the errors are correlated across outputs by
/// its declared correlation.
/// </summary>
public sealed record GaussianMeasurement(IReadOnlyList<VarianceGroup> Groups) : MeasurementErrorModel;

/// <summary>A variance group: the outputs, by index, that share σ²_g, and σ²_g's prior.</summary>
public sealed record VarianceGroup(int[] Outputs, VariancePrior Prior);

/// <summary>A variance group's prior (D09): each group declares its own. Both kinds are proper, so either stands at a
/// b = 0 rung.</summary>
[JsonPolymorphic(TypeDiscriminatorPropertyName = "Kind")]
[JsonDerivedType(typeof(InverseGammaVariance), "InverseGamma")]
[JsonDerivedType(typeof(HalfCauchyDeviation), "HalfCauchy")]
public abstract record VariancePrior
{
    private protected VariancePrior() { }
}

/// <summary>
/// R's default: σ²_g ~ IG(ν/2, (ν/2)·ŝ²), from a prior estimate ŝ of the error's standard deviation held with ν degrees
/// of freedom (R's <c>sd.est</c> and <c>s2.df</c>). Needs ν &gt; 0: a group R would have switched to its improper kernel
/// (<c>s2.df</c> = 0) declares <see cref="HalfCauchyDeviation"/>.
/// </summary>
public sealed record InverseGammaVariance(double DegreesOfFreedom, double Scale) : VariancePrior;

/// <summary>
/// The heavy-tailed choice: σ_g ~ C⁺(0, ŝ) on the standard deviation (Gelman 2006), proper and scale-equivariant, in
/// place of R's (1 + σ²)⁻¹ kernel, which is improper and depends on y's units.
/// </summary>
public sealed record HalfCauchyDeviation(double Scale) : VariancePrior;

/// <summary>
/// A variance group's prior as a density on σ² (moves on log σ² add their Jacobian themselves), with exact draws for
/// starting states and a b = 0 rung. Every built prior here is proper; <see cref="IsProper"/> is what the propriety rule
/// reads.
/// </summary>
internal interface IVariancePrior
{
    /// <summary>ln p(σ²); −∞ for σ² ≤ 0.</summary>
    double LogDensity(double sigma2);

    /// <summary>A draw of σ².</summary>
    double Sample(Xoshiro256PlusPlus rng);

    /// <summary>Whether the density integrates to 1.</summary>
    bool IsProper { get; }
}

/// <summary>IG(ν/2, (ν/2)·ŝ²) on σ².</summary>
internal sealed class InverseGammaVariancePrior : IVariancePrior
{
    private readonly InverseGamma _law;

    public InverseGammaVariancePrior(double degreesOfFreedom, double scale)
        => _law = InverseGamma.ShapeScale(degreesOfFreedom / 2.0, degreesOfFreedom / 2.0 * scale * scale);

    public double LogDensity(double sigma2) => _law.LogPdf(sigma2);
    public double Sample(Xoshiro256PlusPlus rng) => _law.Sample(rng);
    public bool IsProper => true;
}

/// <summary>σ ~ C⁺(0, ŝ) as a density on σ²: p_σ(√σ²)/(2√σ²).</summary>
internal sealed class HalfCauchyVariancePrior : IVariancePrior
{
    private readonly HalfCauchy _law;

    public HalfCauchyVariancePrior(double scale) => _law = HalfCauchy.WithScale(scale);

    public double LogDensity(double sigma2)
    {
        if (!(sigma2 > 0.0) || double.IsPositiveInfinity(sigma2)) return double.NegativeInfinity;
        return _law.LogPdf(Math.Sqrt(sigma2)) - Math.Log(2.0) - 0.5 * Math.Log(sigma2);
    }

    public double Sample(Xoshiro256PlusPlus rng)
    {
        double sigma = _law.Sample(rng);
        return sigma * sigma;
    }

    public bool IsProper => true;
}

/// <summary>
/// The Gaussian measurement model bound to m outputs: each output's group and each group's prior, validated as a
/// partition of the outputs. It maps the group variances to the per-output variances the covariance formulations read.
/// </summary>
public sealed class MeasurementModel
{
    private readonly int[] _groupOf;
    private readonly int[][] _outputs;
    private readonly IVariancePrior[] _priors;

    /// <summary>The model <paramref name="declaration"/> declares, over <paramref name="outputs"/> outputs.</summary>
    public MeasurementModel(MeasurementErrorModel declaration, int outputs)
    {
        ArgumentNullException.ThrowIfNull(declaration);
        if (outputs < 1) throw new ArgumentOutOfRangeException(nameof(outputs), "Calibration needs at least one output.");
        if (declaration is not GaussianMeasurement gaussian) throw new ArgumentException("The measurement model is of no known kind.", nameof(declaration));
        Validate(gaussian);
        _groupOf = new int[outputs];
        Array.Fill(_groupOf, -1);
        int groups = gaussian.Groups.Count;
        _outputs = new int[groups][];
        _priors = new IVariancePrior[groups];
        for (int g = 0; g < groups; g++)
        {
            VarianceGroup group = gaussian.Groups[g];
            foreach (int j in group.Outputs)
            {
                if (j >= outputs) throw new ArgumentException($"Group {g} names output {j}, but there are {outputs} outputs.", nameof(declaration));
                _groupOf[j] = g;
            }
            _outputs[g] = (int[])group.Outputs.Clone();
            _priors[g] = Build(group.Prior);
        }
        int missing = Array.IndexOf(_groupOf, -1);
        if (missing >= 0) throw new ArgumentException($"Output {missing} is in no variance group; every output is in exactly one.", nameof(declaration));
    }

    // A model over explicit priors, for the propriety rule's checks against priors no declaration builds.
    internal MeasurementModel(int[] groupOf, IVariancePrior[] priors)
    {
        _groupOf = groupOf;
        _priors = priors;
        _outputs = new int[priors.Length][];
        for (int g = 0; g < priors.Length; g++)
        {
            var members = new List<int>();
            for (int j = 0; j < groupOf.Length; j++)
                if (groupOf[j] == g) members.Add(j);
            _outputs[g] = [.. members];
        }
    }

    /// <summary>The outputs m.</summary>
    public int Outputs => _groupOf.Length;

    /// <summary>The variance groups G.</summary>
    public int Groups => _priors.Length;

    /// <summary>Output <paramref name="output"/>'s group.</summary>
    public int GroupOf(int output) => _groupOf[output];

    /// <summary>Group <paramref name="group"/>'s outputs, in their declared order.</summary>
    public IReadOnlyList<int> OutputsOf(int group) => _outputs[group];

    /// <summary>Whether every group's prior is proper, as a b = 0 rung needs.</summary>
    public bool IsProper
    {
        get
        {
            foreach (IVariancePrior p in _priors)
                if (!p.IsProper) return false;
            return true;
        }
    }

    /// <summary>Each output's variance σ²_{g(j)} from the group variances <paramref name="variances"/>, into
    /// <paramref name="perOutput"/>.</summary>
    public void OutputVariances(ReadOnlySpan<double> variances, Span<double> perOutput)
    {
        if (variances.Length != Groups) throw new ArgumentException($"Give one variance per group ({Groups}).", nameof(variances));
        if (perOutput.Length != Outputs) throw new ArgumentException($"The result needs one value per output ({Outputs}).", nameof(perOutput));
        for (int j = 0; j < perOutput.Length; j++) perOutput[j] = variances[_groupOf[j]];
    }

    /// <summary>Σ_g ln p_g(σ²_g), the variances' log prior as a density on σ².</summary>
    public double LogPrior(ReadOnlySpan<double> variances)
    {
        if (variances.Length != Groups) throw new ArgumentException($"Give one variance per group ({Groups}).", nameof(variances));
        double sum = 0.0;
        for (int g = 0; g < variances.Length; g++) sum += _priors[g].LogDensity(variances[g]);
        return sum;
    }

    /// <summary>Group <paramref name="group"/>'s log prior at <paramref name="variance"/>.</summary>
    public double LogPrior(int group, double variance) => _priors[group].LogDensity(variance);

    /// <summary>An exact draw of every group's variance from its prior, into <paramref name="variances"/>.</summary>
    public void SamplePrior(Xoshiro256PlusPlus rng, Span<double> variances)
    {
        ArgumentNullException.ThrowIfNull(rng);
        for (int g = 0; g < variances.Length; g++) variances[g] = _priors[g].Sample(rng);
    }

    /// <summary>Checks a Gaussian declaration's own fields (the coverage of the outputs needs their count, which the
    /// constructor has): at least one group, no empty group, no output twice, and each prior's parameters.</summary>
    internal static void Validate(GaussianMeasurement declaration)
    {
        if (declaration.Groups is null || declaration.Groups.Count == 0) throw new ArgumentException("Declare at least one variance group.", nameof(declaration));
        var seen = new HashSet<int>();
        for (int g = 0; g < declaration.Groups.Count; g++)
        {
            VarianceGroup group = declaration.Groups[g] ?? throw new ArgumentException($"Variance group {g} is null.", nameof(declaration));
            if (group.Outputs is null || group.Outputs.Length == 0) throw new ArgumentException($"Variance group {g} names no output.", nameof(declaration));
            foreach (int j in group.Outputs)
            {
                if (j < 0) throw new ArgumentException($"Variance group {g} names output {j}; outputs are indexed from 0.", nameof(declaration));
                if (!seen.Add(j)) throw new ArgumentException($"Output {j} is in more than one variance group.", nameof(declaration));
            }
            _ = Build(group.Prior ?? throw new ArgumentException($"Variance group {g} declares no prior.", nameof(declaration)));
        }
    }

    private static IVariancePrior Build(VariancePrior prior) => prior switch
    {
        InverseGammaVariance ig when !(ig.DegreesOfFreedom > 0.0 && double.IsFinite(ig.DegreesOfFreedom))
            => throw new ArgumentException($"An inverse-gamma group needs finite degrees of freedom above 0, not {ig.DegreesOfFreedom}; " +
                                           "a group R would switch to its improper kernel (s2.df = 0) declares the half-Cauchy.", nameof(prior)),
        InverseGammaVariance ig when !(ig.Scale > 0.0 && double.IsFinite(ig.Scale))
            => throw new ArgumentException($"An inverse-gamma group needs a positive, finite scale ŝ, not {ig.Scale}.", nameof(prior)),
        InverseGammaVariance ig => new InverseGammaVariancePrior(ig.DegreesOfFreedom, ig.Scale),
        HalfCauchyDeviation hc when !(hc.Scale > 0.0 && double.IsFinite(hc.Scale))
            => throw new ArgumentException($"A half-Cauchy group needs a positive, finite scale ŝ, not {hc.Scale}.", nameof(prior)),
        HalfCauchyDeviation hc => new HalfCauchyVariancePrior(hc.Scale),
        _ => throw new ArgumentException("The variance prior is of no known kind.", nameof(prior)),
    };
}
