// ============================================================================
// src/graphs/models/potts/BwdPottsCriticalEstimate.cs
// ============================================================================
// BWD's q-only estimate of the SP→PM pseudotransition temperature for the
// 1/K̂-normalized Gaussian-kernel Potts model:
//
//     T_ps(q) = e^(-1/2) / (4 · ln(1 + sqrt(q)))
//
// Provenance: Blatt, Wiseman and Domany, Neural Computation 9 (1997), eq. 4.2:
// T ≈ exp(−⟨d²⟩/2a²) / (4 log(1 + √q)), Wu's square-lattice T_c = J/ln(1 + √q)
// at four neighbours, with J ≈ exp(−⟨d²⟩/2a²)/K̂; e^(-1/2) takes ⟨d²⟩ = a².
// Valid ONLY when the coupling carries the 1/K̂ mean-degree prefactor (BWD1995:
// the careful 1/K̂ scaling "enables us to estimate T_c for any data sample") —
// without it the temperature scale is data-dependent and this q-only anchor is
// wrong; use a data-dependent estimate instead. q is read from configuration at
// runtime, never assumed: q=2→0.172, 10→0.106, 20→0.089, 100→0.063. It falls with
// q, as the random-cluster comparison inequality requires (Grimmett 2006, (3.22)).
// ============================================================================
using System;

namespace Graphs.Models.Potts;

public static class BwdPottsCriticalEstimate
{
    /// <summary>BWD's prefactor e^(−1/2)/4 ≈ 0.15163.</summary>
    public const double Prefactor = 0.15163266492815836;

    /// <summary>The q-only bracket's cold end, as a multiple of T_ps(q).</summary>
    public const double BracketColdFactor = 0.05;

    /// <summary>The q-only bracket's hot end, as a multiple of T_ps(q).</summary>
    public const double BracketHotFactor = 4.0;

    /// <summary>
    /// SP→PM upper-bound estimate T_ps(q) for the 1/K̂-normalized kernel —
    /// the physical anchor for q-only temperature brackets.
    /// </summary>
    public static double TpsUpperBound(int q)
    {
        if (q < 2)
            throw new ArgumentOutOfRangeException(nameof(q), q, "Potts q must be ≥ 2.");
        return Prefactor / Math.Log(1.0 + Math.Sqrt(q));
    }

    /// <summary>
    /// The q-only temperature bracket for the 1/K̂-normalized kernel, [0.05, 4]·T_ps(q). It reaches past the
    /// melts and most splits of the parity data; a weak cut can still split colder, so the cold end is not a guarantee.
    /// </summary>
    public static (double Lo, double Hi) Bracket(int q)
    {
        double tps = TpsUpperBound(q);
        return (BracketColdFactor * tps, BracketHotFactor * tps);
    }
}
