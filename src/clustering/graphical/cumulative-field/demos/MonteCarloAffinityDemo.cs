using System;
using Clustering.Graphical.SPC.Runtime.Core;
using Clustering.Graphical.SPC.Runtime.Core.Solver;
using Graphs.Primitives;
using Maths.Rng;
using Maths.Samplers.Iid;

namespace Clustering.Graphical.SPC.Demos;

/// <summary>
/// <b>Demonstration, not production.</b> Wang 2020 as published — estimates the
/// per-edge affinity by averaging <paramref name="draws"/> Boltzmann energy
/// samples (an edge is bonded in a draw iff the budget did not reach its
/// <c>Hcum</c>). Converges to <see cref="PKWang.Solve"/> as <c>draws → ∞</c>;
/// the closed form <i>is</i> that limit with zero variance, so the production
/// field never needs this. Kept as the executable witness that Wang's Monte
/// Carlo was estimating a closed-form constant — the empirical face of the
/// scandal whose formal twin is sketched in <c>issues/lean/{spc-lemmas,
/// wang_2020_scandal}</c>. The energy draws are the generic iid
/// <see cref="InverseTransform.Exponential"/> sampler; only the per-edge
/// survival averaging is PKWang-specific.
/// </summary>
public static class WangMonteCarloDemo
{
    /// <summary>
    /// Monte-Carlo estimate of the closed-form affinity at temperature
    /// <paramref name="T"/> over <paramref name="draws"/> Boltzmann energy
    /// samples. Compare against <see cref="PKWang.Solve"/> to witness the
    /// zero-variance limit.
    /// </summary>
    public static Affinities Sample(PKWangContext context, double T, int draws, int? seed = null)
    {
        ArgumentNullException.ThrowIfNull(context);
        if (T <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(T), "Temperature must be positive.");
        if (draws < 1)
            throw new ArgumentOutOfRangeException(nameof(draws), "Draw count must be positive.");

        double[] energies = InverseTransform.Exponential(new Xoshiro256PlusPlus(seed), T, draws);

        double[] hcum = context.Hcum;
        var g = new double[hcum.Length];
        for (int e = 0; e < hcum.Length; e++)
        {
            double h = hcum[e];
            if (h <= 0.0) continue;
            int bonded = 0;
            for (int m = 0; m < draws; m++)
                if (energies[m] < h) bonded++;
            g[e] = bonded / (double)draws;
        }

        if (context.DirectedSymmetrize)
            EdgeFieldSymmetrization.Symmetrize(context.Graph, g, context.Mirror!, context.Rule);

        return new Affinities { Temperature = T, G = g, ReplicaIndex = 0 };
    }
}
