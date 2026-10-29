using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass.Presentation;

/// <summary>
/// One chain's traces (R's <c>plot.bass</c> trace panels; P09): M, σ², τ and λ of each retained draw in draw order, read
/// from the chain's <see cref="Maths.Samplers.Ensemble.DrawStore{TDraw}"/>, with each draw's sweep after burn-in, so the
/// draw axis is explicit under thinning.
/// </summary>
/// <param name="Chain">The chain's index.</param>
/// <param name="Iterations">Each draw's sweep after burn-in, a ladder's scan counting as one: (i·k + 1)·thin for draw i of a
/// store keeping one in k.</param>
/// <param name="Count">Each draw's M, its term count.</param>
/// <param name="Sigma2">Each draw's σ².</param>
/// <param name="Tau">Each draw's τ: the held value under fixed g.</param>
/// <param name="Lambda">Each draw's λ, the count prior's Poisson rate.</param>
public sealed record BassChainTrace(int Chain, long[] Iterations, int[] Count, double[] Sigma2, double[] Tau, double[] Lambda);

/// <summary>
/// A fit's traces (P09), one per chain in chain order: the data behind R's trace plots of M and σ², with τ and λ beside
/// them. <see cref="TauSampled"/> is false under fixed g, whose τ trace is constant.
/// </summary>
public sealed record BassTraces(IReadOnlyList<BassChainTrace> Chains, bool TauSampled, int Thin)
{
    internal static BassTraces Of(BassFit fit)
    {
        ArgumentNullException.ThrowIfNull(fit);
        int thin = fit.Config.Sampling.Thin;
        var chains = new BassChainTrace[fit.Draws.Count];
        for (int c = 0; c < chains.Length; c++)
        {
            var store = fit.Draws[c];
            int n = store.Count;
            var iterations = new long[n];
            var count = new int[n];
            var sigma2 = new double[n];
            var tau = new double[n];
            var lambda = new double[n];
            for (int i = 0; i < n; i++)
            {
                BassDraw draw = store[i];
                iterations[i] = ((long)i * store.Every + 1) * thin;
                count[i] = draw.Structure.TermIds.Length;
                sigma2[i] = draw.Sigma2;
                tau[i] = draw.Tau;
                lambda[i] = draw.Lambda;
            }
            chains[c] = new BassChainTrace(c, iterations, count, sigma2, tau, lambda);
        }
        return new BassTraces(Array.AsReadOnly(chains), fit.Config.Coefficients is SampledG, thin);
    }
}
