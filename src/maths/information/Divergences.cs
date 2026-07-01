// ============================================================================
// Maths.Information/Divergences.cs
// ============================================================================
// KL divergence primitives on dense fixed-support distributions.
//
// Lifted from SpcCore.SpcAnalysis (spc/thermo/kl.cs) and span-ified for
// general use. The thermo-specific wrappers (ComputeFisherInT, temperature
// sweep analysis) remain in SpcAnalysis; this class owns only the
// information-theoretic core.
//
// Zero conventions (consistent with thermo):
//   p[i] = 0 → term is 0 by convention (0 · log(0/q) = 0)
//   q[i] = 0, p[i] > 0 → q is floored at eps, not a hard error
//
// For sparse text distributions where many q[i] are structurally zero,
// pre-smooth q with Histogram.Normalize(alpha > 0) before calling Forward
// rather than relying on the eps floor — smoothing is semantically cleaner
// and makes the divergence more stable across vocabulary sizes.
// ============================================================================

#nullable enable
using System;
using System.Runtime.CompilerServices;

namespace Maths.Information;

public static class KLDivergence
{
    /// <summary>
    /// Forward KL divergence: KL(P ‖ Q) = Σᵢ pᵢ · log(pᵢ / qᵢ).
    /// Result is in nats. Always ≥ 0; equals 0 iff P = Q almost everywhere.
    /// <para>
    /// Explicit length guard: spans are walked in lockstep from <c>p.Length</c>;
    /// mismatched lengths would silently truncate rather than throw at indexing.
    /// </para>
    /// </summary>
    /// <param name="p">Source distribution P.</param>
    /// <param name="q">Target distribution Q.</param>
    /// <param name="eps">Floor applied to q[i] when p[i] &gt; 0 and q[i] ≈ 0.
    /// Prevents −∞ without requiring pre-smoothed Q. For sparse text use
    /// <see cref="T:Hashish.Histogram"/> smoothing instead.</param>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static double Forward(
        ReadOnlySpan<double> p,
        ReadOnlySpan<double> q,
        double eps = 1e-12)
    {
        if (q.Length != p.Length)
            throw new ArgumentException(
                $"KLDivergence requires equal-length distributions (p={p.Length}, q={q.Length}).");

        double sum = 0.0;
        for (int i = 0; i < p.Length; i++)
        {
            double pi = p[i];
            if (pi <= 0.0) continue;
            double qi = q[i] > eps ? q[i] : eps;
            sum += pi * Math.Log(pi / qi);
        }
        return sum;
    }

    /// <summary>
    /// Symmetric KL (Jeffreys divergence): KL(P ‖ Q) + KL(Q ‖ P).
    /// Symmetric and always ≥ 0; equals 0 iff P = Q almost everywhere.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Symmetric(
        ReadOnlySpan<double> p,
        ReadOnlySpan<double> q,
        double eps = 1e-12)
        => Forward(p, q, eps) + Forward(q, p, eps);
}
