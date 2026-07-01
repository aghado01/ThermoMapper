using System;
using System.Runtime.CompilerServices;

namespace Graphs.Models.Potts;

/// <summary>
/// The Fortuin–Kasteleyn bond-activation probability of the Potts model:
/// <c>p = 1 − exp(−energy / T)</c> — the Edwards–Sokal / random-cluster kernel that
/// bridges the Potts couplings to FK bonds. Model vocabulary, shared by every
/// inference strategy: Swendsen–Wang feeds it the raw coupling <c>J_e</c> to place
/// bonds; PKWang feeds it the mean-field cumulative-energy ladder <c>Hcum</c> in
/// closed form. (The bare functional is the exponential CDF <c>1 − exp(−x)</c>; its
/// inverse — the quantile draws — lives in <c>Maths.Samplers.Iid.InverseTransform.Exponential</c>.)
/// </summary>
public static class FkKernel
{
    /// <summary>
    /// FK bond-activation probability <c>p = 1 − exp(−energy / temperature) ∈ [0,1]</c>.
    /// <paramref name="energy"/> is an edge coupling (SW) or cumulative-ladder value
    /// (PKWang); <paramref name="temperature"/> must be positive.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double BondProbability(double energy, double temperature)
        => 1.0 - Math.Exp(-energy / temperature);
}
