using System;
using System.Buffers;
using System.Runtime.CompilerServices;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A basis of directional hinges in place of BASS's axis hinges (E36's slice 2; decision 6, ridges only): each term one directional hinge
/// (θᵀz − t)₊^q / c over the standardized continuous predictors (<see cref="DirectionalHinge"/>), with Collins, Francom and Rumsey's priors
/// and moves. A direction's active count a is capped by <see cref="BassStructureSpec.MaxInteraction"/>, which resolves to the paper's
/// A = min(3, p) at its default; a birth draws a and 𝒥 by D01's order and variable law at <see cref="BassStructureSpec.OrderWeight"/> and
/// <see cref="BassStructureSpec.VariableWeight"/>, the paper's ω₀ and υ₀ (1 there, 5 by BASS's default). Null, the default, is BASS's hinges;
/// left out of the serialized form when null.
/// </summary>
public sealed record DirectionalBasis
{
    /// <summary>p₀ in (0, 1]: the prior probability that the knot lies above the least projected training case, so that the hinge is zero on
    /// some case. The knot is uniform on [U − (U − Q₀)/p₀, U] (the paper's eq. (UL)); 2/3, the paper's and <c>BayesPPR</c>'s default.</summary>
    public double LocalProbability { get; init; } = 2.0 / 3.0;

    /// <summary>About how many training cases the knot's upper end U leaves above it: U is the 1 − ActiveCases/n quantile of the projected
    /// training cases, R's type 7. Null takes min(20, n/10), <c>BayesPPR</c>'s <c>n_dat_min</c>.</summary>
    public double? ActiveCases { get; init; }

    /// <summary>κ, the change move's power spherical concentration on a direction's active coordinates; 1000, the paper's.</summary>
    public double Concentration { get; init; } = 1000.0;

    /// <summary>The first setting out of range, as its name and the reason; null when every one holds.</summary>
    public (string Setting, string Reason)? Problem()
    {
        if (!(LocalProbability > 0.0 && LocalProbability <= 1.0))
            return (nameof(LocalProbability), "must lie in (0, 1]: at 0 the knot's lower end is −∞");
        if (ActiveCases is { } active && !(active > 0.0 && double.IsFinite(active)))
            return (nameof(ActiveCases), "must be positive and finite, or null for min(20, n/10)");
        if (!(Concentration >= 0.0 && double.IsFinite(Concentration))) return (nameof(Concentration), "must be finite and not negative");
        return null;
    }
}

/// <summary>A knot's range along one direction, [<paramref name="Lower"/>, <paramref name="Upper"/>], on which its prior is uniform.</summary>
public readonly record struct KnotRange(double Lower, double Upper)
{
    /// <summary>U − L, whose reciprocal is the knot's prior density.</summary>
    public double Width => Upper - Lower;
}

/// <summary>
/// The term space of a directional basis (E36's slice 2) and its term measure's pieces, for <see cref="BassStructurePrior"/> and E36's
/// proposal (<see cref="DirectionalTermProposal"/>):
/// <list type="bullet">
/// <item><b>The direction</b> (decision 2, the user's): a uniform on {1, …, A}, 𝒥 uniform among the C(p, a) subsets, θ_𝒥 uniform on the
/// subsphere S^(a−1), whose density against surface measure is 1/|S^(a−1)| = Γ(a/2)/(2π^(a/2)); on S⁰ that is ½ on each sign.</item>
/// <item><b>The knot</b> (decision 4): t | θ uniform on [L, U], U = Q_q of the projected training cases θᵀzᵢ (R's type 7, q = 1 − n_active/n)
/// and L = U − (U − Q₀)/p₀ with Q₀ the least projection, the paper's eq. (UL). The range moves with θ, so a change draws its knot given the
/// proposed direction. Where U = Q₀, as ties can make it, the range has no width and the direction no knot: such a term lies outside the
/// support, as a structure outside F does, and the measure is not renormalized for it.</item>
/// </list>
/// Each direction's range is computed once and kept beside it, a pure cache: the value is the same computation whether read from it or not.
/// Immutable otherwise, so one serves every chain of a fit.
/// </summary>
public sealed class DirectionalTermSpace
{
    private const int SelectionLimit = 64;   // past this many order statistics, a sort rather than an insertion buffer

    private readonly BassInputs _inputs;
    private readonly double[] _logDirection;
    private readonly ConditionalWeakTable<ProjectionDirection, StrongBox<KnotRange?>> _ranges = new();

    /// <summary>The space over <paramref name="inputs"/>' continuous predictors, each carrying its standardization, with the active cap
    /// <paramref name="maxActive"/> (1 ≤ A ≤ p), the knot's upper quantile level <paramref name="knotLevel"/> in (0, 1) and p₀ =
    /// <paramref name="localProbability"/> in (0, 1].</summary>
    public DirectionalTermSpace(BassInputs inputs, int maxActive, double knotLevel, double localProbability)
    {
        ArgumentNullException.ThrowIfNull(inputs);
        int p = inputs.Cases.ContinuousCount;
        if (p < 1) throw new ArgumentException("A directional basis needs a continuous predictor.", nameof(inputs));
        if (inputs.Cases.Cases < 2) throw new ArgumentException("A directional basis needs two training cases for a knot's range.", nameof(inputs));
        for (int v = 0; v < p; v++)
            if (inputs.Standardizing(v) is null)
                throw new ArgumentException($"Continuous variable {v}'s scale carries no standardization, the metric a direction is uniform in.", nameof(inputs));
        if (maxActive < 1 || maxActive > p) throw new ArgumentOutOfRangeException(nameof(maxActive), $"A lies in [1, {p}].");
        if (!(knotLevel > 0.0 && knotLevel < 1.0)) throw new ArgumentOutOfRangeException(nameof(knotLevel), "The knot's upper quantile level lies in (0, 1).");
        if (!(localProbability > 0.0 && localProbability <= 1.0)) throw new ArgumentOutOfRangeException(nameof(localProbability), "p₀ lies in (0, 1].");
        _inputs = inputs;
        Variables = p;
        MaxActive = maxActive;
        KnotLevel = knotLevel;
        LocalProbability = localProbability;
        _logDirection = new double[maxActive + 1];
        for (int a = 1; a <= maxActive; a++)
            _logDirection[a] = SpecialFunctions.LogGamma(0.5 * a) - Math.Log(2.0) - 0.5 * a * Math.Log(Math.PI);
    }

    /// <summary>p, the continuous predictors a direction ranges over.</summary>
    public int Variables { get; }

    /// <summary>A, the cap on a direction's active variables.</summary>
    public int MaxActive { get; }

    /// <summary>q, the quantile level of the knot's upper end U.</summary>
    public double KnotLevel { get; }

    /// <summary>p₀, the prior probability that the knot lies above the least projected training case.</summary>
    public double LocalProbability { get; }

    /// <summary>ln of θ_𝒥's uniform density on S^(a−1) at <paramref name="active"/> = a: −ln|S^(a−1)| = ln Γ(a/2) − ln 2 − (a/2) ln π.</summary>
    public double LogDirection(int active)
        => active >= 1 && active <= MaxActive ? _logDirection[active] : throw new ArgumentOutOfRangeException(nameof(active), $"a lies in [1, {MaxActive}].");

    /// <summary>The knot's range along <paramref name="direction"/>, from the projected training cases; null where it has no width.</summary>
    public KnotRange? KnotRangeOf(ProjectionDirection direction)
    {
        ArgumentNullException.ThrowIfNull(direction);
        return _ranges.GetValue(direction, d => new StrongBox<KnotRange?>(Compute(d))).Value;
    }

    /// <summary>ln π(t | θ) for <paramref name="hinge"/>: −ln(U − L) with its knot in [L, U], −∞ outside it or where the range has no width.</summary>
    public double LogKnot(DirectionalHinge hinge)
    {
        ArgumentNullException.ThrowIfNull(hinge);
        if (KnotRangeOf(hinge.Direction) is not { } range) return double.NegativeInfinity;
        return hinge.Knot >= range.Lower && hinge.Knot <= range.Upper ? -Math.Log(range.Width) : double.NegativeInfinity;
    }

    /// <summary>Refuses a term outside the space: one holding any part but a directional hinge, a direction with more than A active variables,
    /// or a variable beyond p.</summary>
    public void Validate(BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(term);
        if (term.Directional is not { } hinge || term.Continuous.Count + term.Categorical.Count + term.Functional.Count + term.Legendre.Count > 0)
            throw new ArgumentException("A directional fit's terms hold one directional hinge alone.", nameof(term));
        ProjectionDirection d = hinge.Direction;
        if (d.Count > MaxActive) throw new ArgumentException($"The direction's {d.Count} active variables exceed the cap A = {MaxActive}.", nameof(term));
        if (d.Variables[^1] >= Variables) throw new ArgumentException($"Continuous variable {d.Variables[^1]} does not exist.", nameof(term));
    }

    /// <summary>θ_𝒥 uniform on the subsphere over <paramref name="variables"/>, distinct: a fair sign at a = 1, a normalized standard normal
    /// vector otherwise.</summary>
    public ProjectionDirection DrawDirection(Xoshiro256PlusPlus rng, ReadOnlySpan<int> variables)
    {
        ArgumentNullException.ThrowIfNull(rng);
        int a = variables.Length;
        if (a < 1 || a > MaxActive) throw new ArgumentOutOfRangeException(nameof(variables), $"a lies in [1, {MaxActive}].");
        if (a == 1) return new ProjectionDirection([variables[0]], [BassDraws.Sign(rng)]);
        var weights = new double[a];
        while (true)
        {
            double squares = 0.0;
            for (int k = 0; k < a; k++)
            {
                weights[k] = Normal.SampleStandard(rng);
                squares += weights[k] * weights[k];
            }
            if (!(squares > 0.0)) continue;
            double inverse = 1.0 / Math.Sqrt(squares);
            bool zero = false;
            for (int k = 0; k < a; k++)
            {
                weights[k] *= inverse;
                zero |= weights[k] == 0.0;
            }
            if (!zero) return new ProjectionDirection(variables.ToArray(), weights);
        }
    }

    /// <summary>A knot uniform on <paramref name="range"/>: U − (U − L)·r with r uniform on [0, 1), <c>BayesPPR</c>'s form.</summary>
    public static double DrawKnot(Xoshiro256PlusPlus rng, KnotRange range)
    {
        ArgumentNullException.ThrowIfNull(rng);
        return range.Upper - range.Width * rng.NextDouble();
    }

    // U as R's type-7 quantile at q, index 1 + (n − 1)q, lo its floor and h its fraction: (1 − h)·x₍lo₎ + h·x₍lo+1₎ where h > 0 and the two
    // differ, x₍lo₎ otherwise; Q₀ the least projection. The order statistics from lo up are kept in an ascending buffer, one comparison
    // for most cases; past SelectionLimit of them, a sort.
    private KnotRange? Compute(ProjectionDirection direction)
    {
        CaseInputs cases = _inputs.Cases;
        int n = cases.Cases;
        double[] rented = ArrayPool<double>.Shared.Rent(n);
        try
        {
            Span<double> projection = rented.AsSpan(0, n);
            DirectionalHingeKernel.Project(direction, _inputs, cases, projection);
            double index = 1.0 + (n - 1) * KnotLevel;
            double floor = Math.Floor(index);
            int lo = (int)floor;   // 1-based
            double h = index - floor;
            int count = n - lo + 1;   // the order statistics x₍lo₎ … x₍n₎
            double least, atLo, atHi;
            if (count > SelectionLimit)
            {
                projection.Sort();
                least = projection[0];
                atLo = projection[lo - 1];
                atHi = h > 0.0 ? projection[lo] : atLo;
            }
            else
            {
                Span<double> top = stackalloc double[count];
                int filled = 0;
                least = double.PositiveInfinity;
                foreach (double x in projection)
                {
                    if (x < least) least = x;
                    if (filled < count)
                    {
                        int at = filled++;
                        while (at > 0 && top[at - 1] > x) { top[at] = top[at - 1]; at--; }
                        top[at] = x;
                    }
                    else if (x > top[0])
                    {
                        int at = 0;
                        while (at + 1 < count && top[at + 1] < x) { top[at] = top[at + 1]; at++; }
                        top[at] = x;
                    }
                }
                atLo = top[0];
                atHi = h > 0.0 ? top[1] : atLo;
            }
            double upper = h > 0.0 && atHi != atLo ? (1.0 - h) * atLo + h * atHi : atLo;
            double width = (upper - least) / LocalProbability;
            if (!(width > 0.0)) return null;
            double lower = upper - width;
            return lower < upper ? new KnotRange(lower, upper) : null;
        }
        finally
        {
            ArrayPool<double>.Shared.Return(rented);
        }
    }
}
