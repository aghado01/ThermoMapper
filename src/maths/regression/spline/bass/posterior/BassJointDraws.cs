using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers.Smc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// Joint draws across a basis fit's components (D07). The components are fitted independently on separate coefficient
/// responses, so their posteriors are independent and any pairing fixed in advance of the draws is valid:
/// <list type="bullet">
/// <item><b>Equal-weight draws</b> (MCMC) pair by index within each chain. When the components' counts differ, as under
/// the R̂ stop, each component's chain is thinned evenly to the shortest count before pairing (*(2026-10-29)*): the
/// j-th of L' kept from L is draw ⌊(2j + 1)L/(2L')⌋, which spreads them over the whole run and keeps every draw when
/// L' = L. The joint draws keep the chains, so PSIS reads each unit's relative efficiency from them.</item>
/// <item><b>Weighted draws</b> (SMC, stage 8) are resampled within each component to the smallest component's count,
/// then shuffled, so that the resampler's sorted ancestors do not pair one component's particle order with another's;
/// an equal-weight component beside them is thinned evenly over its pooled draws. The joint draws then weigh equally and
/// carry no chains.</item>
/// </list>
/// The components must share their chain count when every one is equal-weight.
/// </summary>
public sealed class BassJointDraws
{
    private readonly int[][] _draws;   // [component][joint] → index into the component's pooled draws

    private BassJointDraws(int[][] draws, int[]? chainLengths, bool resampled)
    {
        _draws = draws;
        ChainLengths = chainLengths is null ? null : Array.AsReadOnly(chainLengths);
        Resampled = resampled;
    }

    /// <summary>The joint draws S.</summary>
    public int Count => _draws[0].Length;

    /// <summary>The components r.</summary>
    public int Components => _draws.Length;

    /// <summary>The joint draws per chain, in order, when every component was equal-weight; null after resampling.</summary>
    public IReadOnlyList<int>? ChainLengths { get; }

    /// <summary>Whether weighted components were resampled.</summary>
    public bool Resampled { get; }

    /// <summary>Component <paramref name="component"/>'s draw in joint draw <paramref name="joint"/>, as an index into its
    /// draws pooled in chain order.</summary>
    public int Draw(int component, int joint) => _draws[component][joint];

    /// <summary>
    /// The joint draws of <paramref name="components"/>, a resampling stream per component derived from
    /// <paramref name="seed"/> by the seed tree; <paramref name="scheme"/> resamples weighted ones (multinomial by
    /// default, SMC's).
    /// </summary>
    public static BassJointDraws Pair(IReadOnlyList<BassFit> components, int seed, ResamplingScheme scheme = ResamplingScheme.Multinomial)
    {
        ArgumentNullException.ThrowIfNull(components);
        if (components.Count == 0) throw new ArgumentException("Joint draws need a component.", nameof(components));
        int r = components.Count;
        var pooled = new int[r];
        bool weighted = false;
        for (int k = 0; k < r; k++)
        {
            BassFit fit = components[k] ?? throw new ArgumentException($"Component {k} is null.", nameof(components));
            if (fit.Draws.Count == 0) throw new ArgumentException($"Component {k} holds no chain.", nameof(components));
            foreach (var chain in fit.Draws)
            {
                if (chain.Count == 0) throw new ArgumentException($"Component {k} has a chain with no draw.", nameof(components));
                pooled[k] += chain.Count;
            }
            if (fit.DrawWeights is { } w)
            {
                if (w.Length != pooled[k]) throw new ArgumentException($"Component {k}'s {w.Length} weights do not match its {pooled[k]} draws.", nameof(components));
                weighted = true;
            }
        }
        return weighted ? Resample(components, pooled, seed, scheme) : ByChain(components);
    }

    private static BassJointDraws ByChain(IReadOnlyList<BassFit> components)
    {
        int r = components.Count, chains = components[0].Draws.Count;
        for (int k = 1; k < r; k++)
            if (components[k].Draws.Count != chains)
                throw new ArgumentException($"Component {k} ran {components[k].Draws.Count} chains and component 0 {chains}; " +
                                            "equal-weight draws pair within chains, so the components must share their chain count.", nameof(components));
        var lengths = new int[chains];
        for (int c = 0; c < chains; c++)
        {
            lengths[c] = int.MaxValue;
            for (int k = 0; k < r; k++) lengths[c] = Math.Min(lengths[c], components[k].Draws[c].Count);
        }
        int total = 0;
        foreach (int length in lengths) total += length;
        var draws = new int[r][];
        for (int k = 0; k < r; k++)
        {
            draws[k] = new int[total];
            int start = 0, joint = 0;
            for (int c = 0; c < chains; c++)
            {
                int length = components[k].Draws[c].Count;
                for (int j = 0; j < lengths[c]; j++) draws[k][joint++] = start + Thinned(j, length, lengths[c]);
                start += length;
            }
        }
        return new BassJointDraws(draws, lengths, resampled: false);
    }

    private static BassJointDraws Resample(IReadOnlyList<BassFit> components, int[] pooled, int seed, ResamplingScheme scheme)
    {
        int r = components.Count, count = int.MaxValue;
        foreach (int p in pooled) count = Math.Min(count, p);
        int[] seeds = SeedTree.Derive(seed, r);
        var draws = new int[r][];
        for (int k = 0; k < r; k++)
        {
            draws[k] = new int[count];
            if (components[k].DrawWeights is not { } weights)
            {
                for (int j = 0; j < count; j++) draws[k][j] = Thinned(j, pooled[k], count);
                continue;
            }
            double total = 0.0;
            foreach (double w in weights)
            {
                if (!(w >= 0.0) || double.IsPositiveInfinity(w)) throw new ArgumentException($"Component {k}'s weights must be finite and not negative.", nameof(components));
                total += w;
            }
            if (!(total > 0.0)) throw new ArgumentException($"Component {k}'s weights are all zero.", nameof(components));
            var normalized = new double[weights.Length];
            for (int d = 0; d < normalized.Length; d++) normalized[d] = weights[d] / total;
            var rng = new Xoshiro256PlusPlus(seeds[k]);
            Resampling.Ancestors(scheme, normalized, rng, draws[k]);
            for (int j = count - 1; j > 0; j--)   // Fisher–Yates
            {
                int swap = rng.NextInt(j + 1);
                (draws[k][j], draws[k][swap]) = (draws[k][swap], draws[k][j]);
            }
        }
        return new BassJointDraws(draws, null, resampled: true);
    }

    // The j-th of `kept` draws spread evenly over `length`: ⌊(2j + 1)·length/(2·kept)⌋.
    private static int Thinned(int j, int length, int kept) => (int)((2L * j + 1) * length / (2L * kept));
}
