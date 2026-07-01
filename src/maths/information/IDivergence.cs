// ============================================================================
// Maths.Information/Divergence.cs
// ============================================================================
// Divergence dispatch surface. Distinct from IMeasure<T> because divergences
// are not metrics — they are typically asymmetric (KL), need not satisfy the
// triangle inequality, and may not be bounded. Forcing them through
// IMeasure<T>.Distance(a, b) would lie about commutativity.
//
// The family this slot accommodates over time:
//   - KL (asymmetric f-divergence)
//   - Jeffreys (symmetric KL)
//   - Jensen-Shannon (bounded, symmetric, square root is a metric)
//   - Bhattacharyya, Hellinger
//   - Bregman family (Itakura-Saito, squared-Euclidean, etc.)
// ============================================================================

#nullable enable
using System;

namespace Maths.Information;

/// <summary>
/// Dispatch surface for directed divergences over a domain <typeparamref name="T"/>.
/// <para>
/// <see cref="Forward"/> is order-sensitive: D(a ‖ b) ≠ D(b ‖ a) in general.
/// <see cref="Symmetric"/> returns a symmetrized form (default = Jeffreys-style sum);
/// implementations override when a more efficient or canonical symmetrization exists
/// (e.g. Jensen-Shannon's mid-distribution path).
/// </para>
/// </summary>
public interface IDivergence<T>
{
    /// <summary>Directed divergence D(a ‖ b). Order matters.</summary>
    double Forward(T a, T b);

    /// <summary>Symmetrized form. Default implementation: Forward(a,b) + Forward(b,a).</summary>
    double Symmetric(T a, T b) => Forward(a, b) + Forward(b, a);
}

/// <summary>
/// Dispatch adapter for <see cref="KLDivergence"/> over dense fixed-support
/// distributions. Uses <see cref="ReadOnlyMemory{T}"/> so callers can pass
/// arrays, pooled buffers, or sliced Memory without copying.
/// </summary>
public readonly struct KlDivergence : IDivergence<ReadOnlyMemory<double>>
{
    private readonly double _eps;

    /// <summary>
    /// Construct with a custom q-floor. Pass <c>0</c> to opt out (require pre-smoothed q).
    /// </summary>
    public KlDivergence(double eps = 1e-12) => _eps = eps;

    /// <summary>Forward KL: KL(p ‖ q) = Σᵢ pᵢ · log(pᵢ / qᵢ), in nats.</summary>
    public double Forward(ReadOnlyMemory<double> p, ReadOnlyMemory<double> q)
        => KLDivergence.Forward(p.Span, q.Span, _eps);

    /// <summary>Symmetric KL (Jeffreys): KL(p ‖ q) + KL(q ‖ p), in nats.</summary>
    public double Symmetric(ReadOnlyMemory<double> p, ReadOnlyMemory<double> q)
        => KLDivergence.Symmetric(p.Span, q.Span, _eps);
}
