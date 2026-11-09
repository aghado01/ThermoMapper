using System;
using System.Linq;
using System.Text.Json.Serialization;
using Maths.Distributions;

namespace Maths.Regression.GPrior;

/// <summary>Which prior on several outputs' Σ an <see cref="OutputCovariancePrior"/> is; there, the parameters each kind takes.</summary>
public enum OutputCovariancePriorKind
{
    /// <summary>IW(K − 1 + 2c, 2c·S_Y), S_Y read from the response: E34's decision 2, the default.</summary>
    ScaledToResponse,

    /// <summary>IW(K − 1 + 2c, 2c·diag S_Y): each output's own scale alone.</summary>
    Standardized,

    /// <summary>IW(ν₀, Ψ₀) declared in the outputs' own units.</summary>
    InverseWishart,

    /// <summary>Huang and Wand's half-t marginal scales with conditional inverse-Wishart rates.</summary>
    HuangWand,
}

/// <summary>
/// The prior on several outputs' noise covariance Σ (E34), the K-output counterpart of <see cref="Sigma2Prior"/>, resolved per response to
/// the <see cref="Distributions.InverseWishart"/> the K-output block integrates Σ against (<see cref="OutputMarginal"/>).
/// <list type="bullet">
/// <item><see cref="ScaledToResponse"/>: decision 2, IW(K − 1 + 2c, 2c·S_Y), with S_Y the outputs' mean squares and cross-products about their
/// centre levels (scalar BASS's s² carried to K outputs), an empirical prior stated as one. At K = 1 it is IG(c, c·s²), the scalar default;
/// and since S_Y → LᵀS_YL when Y → YL, the structure posterior is invariant to any invertible rescaling or mixing of the outputs. Under a
/// diagonal S_Y each prior correlation has density ∝ (1 − ρ²)^(c − 1) (Huang and Wand's Property 3 at their ν = 2c), at c = 10⁻³ nearly all
/// at ±1, which the data overrule at moderate n/K; E34's slice 3 measures it at small n/K.</item>
/// <item><see cref="Standardized"/>: IW(K − 1 + 2c, 2c·diag S_Y), the default on outputs standardized to unit mean square: invariant to
/// rescaling the outputs, not to mixing them, and centred on uncorrelated outputs.</item>
/// <item><see cref="InverseWishart"/>: IW(ν₀, Ψ₀) in the outputs' own units, ν₀ &gt; K − 1, for an analysis that fixes it; free of the
/// data, as Geweke's test needs.</item>
/// <item><see cref="HuangWand"/>: half-t marginal standard deviations, represented by K gamma rates and their conditional inverse Wishart.</item>
/// </list>
/// A value, compared by value and serialized as it stands, Ψ₀ as rows; the constructor checks the parameters against the kind, so a prior
/// read back from JSON is checked as a built one is.
/// </summary>
public sealed record OutputCovariancePrior
{
    /// <summary>
    /// The prior of kind <paramref name="kind"/>: <see cref="OutputCovariancePriorKind.ScaledToResponse"/> and
    /// <see cref="OutputCovariancePriorKind.Standardized"/> take a shape c, positive and finite; <see cref="OutputCovariancePriorKind.InverseWishart"/>
    /// takes the degrees of freedom ν₀ and the scale Ψ₀ as rows, symmetric positive definite, ν₀ &gt; K − 1.
    /// Huang and Wand's prior takes positive ν and optional positive marginal scales.
    /// </summary>
    public OutputCovariancePrior(OutputCovariancePriorKind kind, double? shape = null, double? degreesOfFreedom = null, double[][]? scale = null,
                                 double[]? outputScales = null)
    {
        switch (kind)
        {
            case OutputCovariancePriorKind.ScaledToResponse or OutputCovariancePriorKind.Standardized:
                if (shape is not (> 0.0 and < double.PositiveInfinity)) throw new ArgumentOutOfRangeException(nameof(shape), "The shape must be positive and finite.");
                if (degreesOfFreedom is not null || scale is not null) throw new ArgumentException("A prior scaled to the response takes its shape alone.");
                break;
            case OutputCovariancePriorKind.InverseWishart:
                if (shape is not null) throw new ArgumentException("A declared inverse Wishart takes ν₀ and Ψ₀, not a shape.", nameof(shape));
                if (degreesOfFreedom is not { } nu || scale is null) throw new ArgumentException("A declared inverse Wishart takes ν₀ and Ψ₀.");
                _ = Distributions.InverseWishart.Of(nu, Matrix(scale));   // validates both
                break;
            case OutputCovariancePriorKind.HuangWand:
                if (shape is not null || scale is not null) throw new ArgumentException("Huang and Wand's prior takes ν and the marginal scales.");
                if (degreesOfFreedom is not (> 0 and < double.PositiveInfinity)) throw new ArgumentOutOfRangeException(nameof(degreesOfFreedom));
                if (outputScales is not null && (outputScales.Length == 0 || outputScales.Any(s => !(s > 0 && double.IsFinite(s)))))
                    throw new ArgumentException("The marginal scales must be positive and finite.", nameof(outputScales));
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(kind));
        }
        if (kind != OutputCovariancePriorKind.HuangWand && outputScales is not null)
            throw new ArgumentException("Marginal scales belong to Huang and Wand's prior.", nameof(outputScales));
        Kind = kind;
        Shape = shape;
        DegreesOfFreedom = degreesOfFreedom;
        Scale = scale is null ? null : [.. scale.Select(row => (double[])row.Clone())];
        OutputScales = outputScales is null ? null : (double[])outputScales.Clone();
    }

    /// <summary>Decision 2's IW(K − 1 + 2c, 2c·S_Y), c = <paramref name="shape"/>.</summary>
    public static OutputCovariancePrior ScaledToResponse(double shape = 1e-3) => new(OutputCovariancePriorKind.ScaledToResponse, shape);

    /// <summary>IW(K − 1 + 2c, 2c·diag S_Y), c = <paramref name="shape"/>.</summary>
    public static OutputCovariancePrior Standardized(double shape = 1e-3) => new(OutputCovariancePriorKind.Standardized, shape);

    /// <summary>Half-t(ν, Aₖ) marginal standard deviations, with uniform correlations at ν = 2.
    /// Null scales take each output's root mean square about its prior level.</summary>
    public static OutputCovariancePrior HuangWand(double nu = 2, double[]? scales = null)
        => new(OutputCovariancePriorKind.HuangWand, degreesOfFreedom: nu, outputScales: scales);

    /// <summary>IW(<paramref name="degreesOfFreedom"/>, <paramref name="scale"/>) in the outputs' own units.</summary>
    public static OutputCovariancePrior InverseWishart(double degreesOfFreedom, double[,] scale)
    {
        ArgumentNullException.ThrowIfNull(scale);
        int k = scale.GetLength(0);
        var rows = new double[k][];
        for (int i = 0; i < k; i++)
        {
            rows[i] = new double[scale.GetLength(1)];
            for (int j = 0; j < rows[i].Length; j++) rows[i][j] = scale[i, j];
        }
        return new(OutputCovariancePriorKind.InverseWishart, degreesOfFreedom: degreesOfFreedom, scale: rows);
    }

    /// <summary>Which prior this is.</summary>
    public OutputCovariancePriorKind Kind { get; }

    /// <summary>The shape c of the priors scaled to the response; null for a declared one.</summary>
    public double? Shape { get; }

    /// <summary>A declared inverse Wishart's ν₀ or Huang and Wand's marginal degrees of freedom ν.</summary>
    public double? DegreesOfFreedom { get; }

    /// <summary>A declared prior's Ψ₀, as rows; null otherwise. Not to be written.</summary>
    public double[][]? Scale { get; }

    /// <summary>Huang and Wand's declared marginal scales, or null for response root mean squares. Not to be written.</summary>
    [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
    public double[]? OutputScales { get; }

    /// <summary>Whether <see cref="For"/> reads the response's mean squares; read from the kind, so not serialized.</summary>
    [JsonIgnore]
    public bool ReadsResponse => Kind != OutputCovariancePriorKind.InverseWishart &&
                                 (Kind != OutputCovariancePriorKind.HuangWand || OutputScales is null);

    /// <summary>
    /// The prior on K outputs whose mean squares and cross-products about their levels are <paramref name="meanSquares"/> (K × K, positive
    /// definite for <see cref="OutputCovariancePriorKind.ScaledToResponse"/>, a positive diagonal for
    /// <see cref="OutputCovariancePriorKind.Standardized"/>); a declared prior refuses another K.
    /// For Huang and Wand's mixture this is the conditional inverse Wishart at the rates' prior means.
    /// </summary>
    public Distributions.InverseWishart For(double[,] meanSquares)
    {
        ArgumentNullException.ThrowIfNull(meanSquares);
        int k = meanSquares.GetLength(0);
        if (meanSquares.GetLength(1) != k) throw new ArgumentException("The mean squares are a square matrix.", nameof(meanSquares));
        if (ResolveRates(meanSquares) is { } rates) return rates.PriorAt(rates.Centre());
        if (Kind == OutputCovariancePriorKind.InverseWishart)
        {
            double[,] declared = Matrix(Scale!);
            if (declared.GetLength(0) != k) throw new ArgumentException($"The declared Ψ₀ is {declared.GetLength(0)} × {declared.GetLength(0)}, for {k} outputs.", nameof(meanSquares));
            return Distributions.InverseWishart.Of(DegreesOfFreedom!.Value, declared);
        }
        double c = Shape!.Value;
        var psi = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int d = 0; d < k; d++)
                psi[a, d] = Kind == OutputCovariancePriorKind.ScaledToResponse || a == d ? 2.0 * c * meanSquares[a, d] : 0.0;
        return Distributions.InverseWishart.Of(k - 1 + 2.0 * c, psi);
    }

    internal ResolvedOutputRates? ResolveRates(double[,] meanSquares)
    {
        if (Kind != OutputCovariancePriorKind.HuangWand) return null;
        int k = meanSquares.GetLength(0);
        if (OutputScales is not null && OutputScales.Length != k) throw new ArgumentException("One marginal scale per output is required.");
        double[] scales = OutputScales ?? [.. Enumerable.Range(0, k).Select(i => Math.Sqrt(meanSquares[i, i]))];
        return new ResolvedOutputRates(DegreesOfFreedom!.Value, scales);
    }

    /// <inheritdoc/>
    public bool Equals(OutputCovariancePrior? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        if (Kind != other.Kind || Shape != other.Shape || DegreesOfFreedom != other.DegreesOfFreedom) return false;
        if (!(OutputScales ?? []).AsSpan().SequenceEqual(other.OutputScales ?? [])) return false;
        if (Scale is null || other.Scale is null) return Scale is null && other.Scale is null;
        return Scale.Length == other.Scale.Length && Scale.Zip(other.Scale).All(rows => rows.First.AsSpan().SequenceEqual(rows.Second));
    }

    /// <inheritdoc/>
    public override int GetHashCode()
    {
        var hash = new HashCode();
        hash.Add(Kind);
        hash.Add(Shape);
        hash.Add(DegreesOfFreedom);
        if (OutputScales is not null)
            foreach (double s in OutputScales) hash.Add(s);
        if (Scale is not null)
            foreach (double[] row in Scale)
                foreach (double v in row) hash.Add(v);
        return hash.ToHashCode();
    }

    /// <summary>The spelling a log reads.</summary>
    public override string ToString() => Kind switch
    {
        OutputCovariancePriorKind.ScaledToResponse => FormattableString.Invariant($"ScaledToResponse({Shape:R})"),
        OutputCovariancePriorKind.Standardized => FormattableString.Invariant($"Standardized({Shape:R})"),
        OutputCovariancePriorKind.HuangWand => FormattableString.Invariant($"HuangWand({DegreesOfFreedom:R})"),
        _ => FormattableString.Invariant($"InverseWishart({DegreesOfFreedom:R}, {Scale!.Length} × {Scale.Length})"),
    };

    private static double[,] Matrix(double[][] rows)
    {
        int k = rows.Length;
        if (k == 0) throw new ArgumentException("Ψ₀ has no rows.", nameof(rows));
        var m = new double[k, k];
        for (int i = 0; i < k; i++)
        {
            double[] row = rows[i] ?? throw new ArgumentException($"Ψ₀'s row {i} is null.", nameof(rows));
            if (row.Length != k) throw new ArgumentException("Ψ₀ is square.", nameof(rows));
            for (int j = 0; j < k; j++) m[i, j] = row[j];
        }
        return m;
    }
}
