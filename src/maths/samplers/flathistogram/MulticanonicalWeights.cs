using System;
using System.Collections.Generic;

namespace Maths.Samplers.FlatHistogram;

/// <summary>
/// A multicanonical weight function w(ℓ) of the untempered log likelihood ℓ = ℓ₁ (Berg and Neuhaus 1992): ln w is
/// piecewise linear between K equally spaced nodes over a window [ℓ_lo, ℓ_hi], and continues outside it with the
/// Boltzmann slopes, b_lo below and 1 above, continuous at both edges. With w ≈ 1/g, g the prior's density of states,
/// ℓ's histogram under <c>prior · w(ℓ)</c> is flat over the window, while below it the chain sees the hot end's
/// tempered law and above it the posterior's. The field a <see cref="MulticanonicalTarget{TState}"/> reweights the
/// prior by, and the state <see cref="SamcLearning"/> learns. Immutable: the learning makes new weights.
/// </summary>
public sealed class MulticanonicalWeights
{
    private readonly double[] _logWeights;

    /// <param name="low">ℓ_lo, the window's low edge, finite.</param>
    /// <param name="high">ℓ_hi, above <paramref name="low"/>, finite.</param>
    /// <param name="lowSlope">b_lo in [0, 1], ln w's slope below the window: the hot end the weights reach down to.</param>
    /// <param name="logWeights">ln w at the K ≥ 2 nodes, finite; the array is copied. ln w is defined up to an additive
    /// constant, which changes no ratio a sampler or the reweighting reads.</param>
    public MulticanonicalWeights(double low, double high, double lowSlope, IReadOnlyList<double> logWeights)
    {
        ArgumentNullException.ThrowIfNull(logWeights);
        if (!double.IsFinite(low) || !double.IsFinite(high) || !(high > low))
            throw new ArgumentOutOfRangeException(nameof(high), "The window needs finite edges with ℓ_hi > ℓ_lo.");
        if (!(lowSlope >= 0.0 && lowSlope <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(lowSlope), "The low slope b_lo must lie in [0, 1].");
        if (logWeights.Count < 2) throw new ArgumentException("The window needs at least two nodes.", nameof(logWeights));
        _logWeights = new double[logWeights.Count];
        for (int k = 0; k < _logWeights.Length; k++)
        {
            if (!double.IsFinite(logWeights[k])) throw new ArgumentException($"ln w at node {k} is not finite.", nameof(logWeights));
            _logWeights[k] = logWeights[k];
        }
        Low = low;
        High = high;
        LowSlope = lowSlope;
        Spacing = (high - low) / (_logWeights.Length - 1);
        if (!(Spacing > 0.0)) throw new ArgumentException("The nodes are too many for the window's width.", nameof(logWeights));
    }

    /// <summary>Flat weights, ln w = 0 at every node: inside the window π_w is the prior; SAMC's start without a pilot.</summary>
    public static MulticanonicalWeights Flat(double low, double high, double lowSlope, int nodes)
    {
        if (nodes < 2) throw new ArgumentOutOfRangeException(nameof(nodes), "The window needs at least two nodes.");
        return new MulticanonicalWeights(low, high, lowSlope, new double[nodes]);
    }

    /// <summary>ℓ_lo.</summary>
    public double Low { get; }

    /// <summary>ℓ_hi.</summary>
    public double High { get; }

    /// <summary>b_lo, ln w's slope below the window; above it the slope is 1.</summary>
    public double LowSlope { get; }

    /// <summary>K, the number of nodes.</summary>
    public int Nodes => _logWeights.Length;

    /// <summary>The nodes' spacing Δ = (ℓ_hi − ℓ_lo)/(K − 1).</summary>
    public double Spacing { get; }

    /// <summary>ln w at the nodes, in order.</summary>
    public IReadOnlyList<double> LogWeights => _logWeights;

    /// <summary>Node k's ℓ: ℓ_lo + k·Δ, exactly ℓ_hi at the last.</summary>
    public double Node(int k)
    {
        if ((uint)k >= (uint)_logWeights.Length) throw new ArgumentOutOfRangeException(nameof(k));
        return k == _logWeights.Length - 1 ? High : Low + k * Spacing;
    }

    /// <summary>The same window and slope with other node values (copied).</summary>
    public MulticanonicalWeights WithLogWeights(IReadOnlyList<double> logWeights)
    {
        ArgumentNullException.ThrowIfNull(logWeights);
        if (logWeights.Count != _logWeights.Length)
            throw new ArgumentException($"The weights have {_logWeights.Length} nodes, not {logWeights.Count}.", nameof(logWeights));
        return new MulticanonicalWeights(Low, High, LowSlope, logWeights);
    }

    /// <summary>
    /// ln w(ℓ): the interpolant between the nodes inside the window, <c>ln w_0 + b_lo·(ℓ − ℓ_lo)</c> below it and
    /// <c>ln w_{K−1} + (ℓ − ℓ_hi)</c> above it. A zero likelihood, ℓ = −∞, has ln w = −∞ when b_lo &gt; 0 and ln w_0 at
    /// b_lo = 0, whose extension below is flat; NaN stays NaN.
    /// </summary>
    public double LogWeight(double level)
    {
        if (double.IsNaN(level)) return double.NaN;
        if (level < Low)
        {
            if (LowSlope == 0.0) return _logWeights[0];
            return double.IsNegativeInfinity(level) ? double.NegativeInfinity : _logWeights[0] + LowSlope * (level - Low);
        }
        if (level > High) return _logWeights[^1] + (level - High);
        Locate(level, out int segment, out double fraction);
        return (1.0 - fraction) * _logWeights[segment] + fraction * _logWeights[segment + 1];
    }

    /// <summary>
    /// Where <paramref name="level"/> falls: false outside the window (or NaN); inside it, the segment j in [0, K − 2]
    /// it lies in and its fraction u in [0, 1] of the way from node j to node j + 1, so the hat functions there are
    /// φ_j = 1 − u and φ_{j+1} = u. ℓ_hi falls in the last segment at u = 1.
    /// </summary>
    public bool TryLocate(double level, out int segment, out double fraction)
    {
        if (!(level >= Low && level <= High))
        {
            segment = -1;
            fraction = double.NaN;
            return false;
        }
        Locate(level, out segment, out fraction);
        return true;
    }

    /// <summary>
    /// The node occupancies a histogram flat over the window gives: each node's hat function's share of the window,
    /// 1/(K − 1) inside and half that at the two edge nodes, summing to 1. SAMC's desired sampling distribution.
    /// </summary>
    public double[] FlatOccupancy()
    {
        int k = _logWeights.Length;
        var pi = new double[k];
        double interior = 1.0 / (k - 1);
        for (int i = 0; i < k; i++) pi[i] = interior;
        pi[0] = pi[k - 1] = 0.5 * interior;
        return pi;
    }

    private void Locate(double level, out int segment, out double fraction)
    {
        int last = _logWeights.Length - 2;
        double position = (level - Low) / Spacing;
        segment = Math.Min((int)Math.Floor(position), last);
        if (segment < 0) segment = 0;
        fraction = Math.Clamp(position - segment, 0.0, 1.0);
    }
}
