using System;
using System.Text.Json.Serialization;
using Maths.Samplers.Ensemble;

namespace Maths.Samplers.FlatHistogram;

/// <summary>Where a flat-histogram run's weights start.</summary>
public enum FlatHistogramStart
{
    /// <summary>ln w = 0 over the window: SAMC discovers g from the prior.</summary>
    Flat,

    /// <summary>ln g from a tempering pilot's rungs by multiple-histogram reweighting (<see cref="FlatHistogram.Wham"/>), which the
    /// learning then corrects.</summary>
    Wham,
}

/// <summary>A window [ℓ_lo, ℓ_hi] of the untempered log likelihood.</summary>
public sealed record FlatHistogramWindow(double Low, double High);

/// <summary>
/// A flat-histogram run's configuration, values only (<see cref="MulticanonicalRun"/>): the chains, the learning and
/// production lengths, where the weights start, the window, b_lo, K and the gain. What is left null is read from the
/// tempering pilot the run is given: the window from ℓ₁'s quantiles at the pilot's hot end and at b = 1, with margins,
/// b_lo as the pilot's hottest rung, and K from the window's width over the posterior's ℓ₁ sd. Two options move the
/// window's top, both off by default: a margin above not held inside the pilot's draws, and the learning's own visits
/// raising it while it learns. Lengths are totals, which a run reaches in one call or in several, each resuming the last.
/// </summary>
public sealed record FlatHistogramSpec
{
    /// <summary>The default chain count, the ensemble's: Belardinelli and Pereyra (2016) find that walkers beyond a
    /// critical count no longer reduce the learning error.</summary>
    public const int DefaultChains = 4;

    /// <summary>The chains: they share the weights while learning and run independently in production.</summary>
    public int Chains { get; init; } = DefaultChains;

    /// <summary>The learning's iterations, each one step of every chain and one update; ≥ 0 (0 runs production on the
    /// starting weights).</summary>
    public required long LearningIterations { get; init; }

    /// <summary>Production draws per chain on the frozen weights, ≥ 1.</summary>
    public required int ProductionDraws { get; init; }

    /// <summary>Where the weights start.</summary>
    public FlatHistogramStart Start { get; init; } = FlatHistogramStart.Flat;

    /// <summary>The window; null reads it from the pilot.</summary>
    public FlatHistogramWindow? Window { get; init; }

    /// <summary>b_lo in [0, 1], ln w's slope below the window; null takes the pilot's hottest rung.</summary>
    public double? LowSlope { get; init; }

    /// <summary>K ≥ 2 nodes; null reads it from the pilot: ⌈width / sd⌉ + 1, sd the b = 1 rung's ℓ₁ sd, so each segment is
    /// at most a posterior sd wide.</summary>
    public int? Nodes { get; init; }

    /// <summary>The gain's scale t₀ ≥ 1; null is 10·K, inside Liang et al.'s suggested 2K to 100K.</summary>
    public double? GainScale { get; init; }

    /// <summary>The iteration the gain's count starts from, ≥ 0; null is 0 from flat and, from WHAM, the pilot's draws
    /// over the chain count, the iterations its draws are worth at one draw per chain visit.</summary>
    public long? GainStart { get; init; }

    /// <summary>q in (0, 0.5): a pilot window runs from ℓ₁'s q-quantile at the hot end to its (1 − q)-quantile at b = 1.</summary>
    public double WindowQuantile { get; init; } = 0.01;

    /// <summary>A pilot window's margins, in each end rung's ℓ₁ sd, ≥ 0, never past the pilot's extreme draws.</summary>
    public double WindowMargin { get; init; } = 1.0;

    /// <summary>A pilot window's margin above, in the b = 1 rung's ℓ₁ sd, ≥ 0: the top at that rung's (1 − q)-quantile plus
    /// this margin, not held inside the pilot's highest draw, so the window can reach above what the pilot drew, where the
    /// start's nodes hold no pilot draw and the learning fills them; past the highest ℓ₁ any state reaches they hold no mass,
    /// SAMC's limit there does not exist, and a WHAM start may be refused. Null takes <see cref="WindowMargin"/>, held
    /// there.</summary>
    public double? WindowTopMargin { get; init; }

    /// <summary>The learning iterations, from its start, during which its visits raise the window's top: a visit a node
    /// spacing or more above it adds nodes at that spacing up to the highest at or below the visit, ln w on the slope-1
    /// extension there, so the weights are unchanged as a function and only the update's reach grows. In
    /// [0, <see cref="LearningIterations"/>]; 0 keeps the window fixed. The learning is burn-in, so production's window is
    /// frozen with its weights (D10).</summary>
    public long ExtendTopUntil { get; init; }

    /// <summary>Whether a run needs a tempering pilot: a WHAM start, or a window, b_lo or K to read from one.</summary>
    [JsonIgnore]
    public bool NeedsPilot => Start == FlatHistogramStart.Wham || Window is null || LowSlope is null || Nodes is null;

    /// <summary>Refuse an ill-formed configuration (<see cref="ArgumentException"/>).</summary>
    public void Validate()
    {
        if (Chains < 1) throw new ArgumentException("A run needs at least one chain.", nameof(Chains));
        if (LearningIterations < 0) throw new ArgumentException("The learning's length must be non-negative.", nameof(LearningIterations));
        if (ProductionDraws < 1) throw new ArgumentException("Production needs at least one draw per chain.", nameof(ProductionDraws));
        if (Start is not (FlatHistogramStart.Flat or FlatHistogramStart.Wham)) throw new ArgumentException("Unknown start.", nameof(Start));
        if (Window is { } w && !(double.IsFinite(w.Low) && double.IsFinite(w.High) && w.High > w.Low))
            throw new ArgumentException("The window needs finite edges with ℓ_hi > ℓ_lo.", nameof(Window));
        if (LowSlope is double b && !(b >= 0.0 && b <= 1.0)) throw new ArgumentException("b_lo must lie in [0, 1].", nameof(LowSlope));
        if (Nodes is int k && k < 2) throw new ArgumentException("The window needs at least two nodes.", nameof(Nodes));
        if (GainScale is double t0 && !(t0 >= 1.0 && double.IsFinite(t0)))
            throw new ArgumentException("The gain's scale t₀ must be finite and at least 1.", nameof(GainScale));
        if (GainStart is long t && t < 0) throw new ArgumentException("The gain's start must be non-negative.", nameof(GainStart));
        if (!(WindowQuantile > 0.0 && WindowQuantile < 0.5)) throw new ArgumentException("The window's quantile must lie in (0, 0.5).", nameof(WindowQuantile));
        if (!(WindowMargin >= 0.0 && double.IsFinite(WindowMargin))) throw new ArgumentException("The window's margin must be finite and non-negative.", nameof(WindowMargin));
        if (WindowTopMargin is double top && !(top >= 0.0 && double.IsFinite(top)))
            throw new ArgumentException("The window's top margin must be finite and non-negative.", nameof(WindowTopMargin));
        if (ExtendTopUntil < 0 || ExtendTopUntil > LearningIterations)
            throw new ArgumentException("The top's extension runs within the learning: 0 ≤ ExtendTopUntil ≤ LearningIterations.", nameof(ExtendTopUntil));
    }

    // Every field, for a run's fingerprint.
    internal void AddTo(RunFingerprint fingerprint)
    {
        fingerprint.Add("flat-histogram").Add(Chains).Add((ulong)LearningIterations).Add(ProductionDraws).Add((int)Start)
                   .Add(Window is null ? 0 : 1).Add(Window?.Low ?? 0.0).Add(Window?.High ?? 0.0)
                   .Add(LowSlope is null ? 0 : 1).Add(LowSlope ?? 0.0).Add(Nodes ?? 0)
                   .Add(GainScale is null ? 0 : 1).Add(GainScale ?? 0.0).Add(GainStart is null ? 0 : 1).Add((ulong)(GainStart ?? 0))
                   .Add(WindowQuantile).Add(WindowMargin).Add(WindowTopMargin is null ? 0 : 1).Add(WindowTopMargin ?? 0.0).Add((ulong)ExtendTopUntil);
    }
}
