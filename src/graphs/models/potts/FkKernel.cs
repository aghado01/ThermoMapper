using System;
using System.Runtime.CompilerServices;

namespace Graphs.Models.Potts;

/// <summary>
/// The Fortuin–Kasteleyn bond-activation probability of the Potts model:
/// <c>p = 1 − exp(−energy / T)</c> — the Edwards–Sokal / random-cluster kernel that
/// bridges the Potts couplings to FK bonds. Model vocabulary, shared by every
/// inference strategy: Swendsen–Wang feeds it the raw coupling <c>J_e</c> to place
/// bonds; the cumulative-field solver feeds it the cumulative-energy ladder <c>Hcum</c> in
/// closed form. (The bare functional is the exponential CDF <c>1 − exp(−x)</c>; its
/// inverse — the quantile draws — lives in <c>Maths.Samplers.Iid.InverseTransform.Exponential</c>.)
/// </summary>
public static class FkKernel
{
    /// <summary>
    /// FK bond-activation probability <c>p = 1 − exp(−energy / temperature) ∈ [0,1]</c>.
    /// <paramref name="energy"/> is an edge coupling (SW) or cumulative-ladder value
    /// (the cumulative-field solver); <paramref name="temperature"/> must be positive.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double BondProbability(double energy, double temperature)
        => 1.0 - Math.Exp(-energy / temperature);

    /// <summary>
    /// The same probability from the logs of the energy and the temperature, <c>p = −expm1(−exp(ln E − ln T))</c>: exact
    /// when either lies past the double range (a coupling of e^−1250 at a temperature of e^−1240), and by expm1 it keeps
    /// its digits for a small ratio, where <c>1 − exp(−x)</c> rounds to 0 below about 1e-16. A zero energy
    /// (<c>ln E = −∞</c>) gives 0.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double LogBondProbability(double logEnergy, double logTemperature)
        => -ExpM1(-Math.Exp(logEnergy - logTemperature));

    /// <summary>
    /// e^x − 1 with full relative precision for small |x| (Kahan): .NET's <c>double.ExpM1</c> is <c>Math.Exp(x) − 1</c>, which
    /// returns 0 once e^x rounds to 1.
    /// </summary>
    private static double ExpM1(double x)
    {
        double u = Math.Exp(x);
        if (u == 1.0) return x;
        double um1 = u - 1.0;
        if (um1 == -1.0) return -1.0;
        return um1 * x / Math.Log(u);
    }
}
