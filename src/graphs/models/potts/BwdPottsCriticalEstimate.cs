// ============================================================================
// src/graphs/models/potts/BwdPottsCriticalEstimate.cs
// ============================================================================
// BWD1996 q-only estimate of the SP→PM pseudotransition temperature for the
// 1/K̂-normalized Gaussian-kernel Potts model:
//
//     T_ps(q) = (e^(-1/2) / 4) · ln(1 + sqrt(q))
//
// Provenance: read off the PRL manually (the docling repair garbled the
// formula); recorded in .discussion/issues/spc-parity/domany-parity-plan.md
// § "Scheduling scope". Valid ONLY when the coupling carries the 1/K̂
// mean-degree prefactor (BWD1995: the careful 1/K̂ scaling "enables us to
// estimate T_c for any data sample") — without it the temperature scale is
// data-dependent and this q-only anchor is wrong; use a data-dependent
// estimate instead. q is read from configuration at runtime, never assumed:
// q=2→0.13, 10→0.22, 20→0.26, 100→0.36 — logarithmic in √q, so the anchor
// is flat and a sloppy q barely moves it (matches BWD's "influence of q on
// the assignment is very weak").
// ============================================================================
using System;

namespace Graphs.Models.Potts;

public static class BwdPottsCriticalEstimate
{
    /// <summary>The BWD1996 prefactor e^(−1/2)/4 ≈ 0.15163.</summary>
    public const double Prefactor = 0.15163266492815836;

    /// <summary>
    /// SP→PM upper-bound estimate T_ps(q) for the 1/K̂-normalized kernel —
    /// the physical anchor for q-only temperature brackets.
    /// </summary>
    public static double TpsUpperBound(int q)
    {
        if (q < 2)
            throw new ArgumentOutOfRangeException(nameof(q), q, "Potts q must be ≥ 2.");
        return Prefactor * Math.Log(1.0 + Math.Sqrt(q));
    }
}
