using System;
using Graphs.Primitives;
using Maths.Information;

namespace Graphs.Observables;

/// <summary>
/// Shannon entropy (nats) of the affinity field <see cref="Affinities.G"/> read as a distribution
/// over edges — a degree-0 (field→value) reduction of the universal currency measuring how evenly
/// bond activity is <i>spread</i> across the graph (high = dispersed, low = concentrated on a few
/// edges). The dispersion face of the per-edge field.
/// </summary>
/// <remarks>
/// <para><b>Distinct from <see cref="AffinityBinaryEntropySum"/></b> — same currency, different
/// question. This is the entropy of the <i>normalized</i> field (one scalar, activity dispersion);
/// <c>AffinityBinaryEntropySum</c> sums the per-edge <i>binary</i> entropy <c>H₂(G_e)</c> (bond
/// uncertainty / flickering). Both are kept as sibling channels.</para>
///
/// <para><b>Reads the minted currency only</b> — never a sampler's draw stream. The entropy is a
/// nonlinear reduction applied <i>once</i> to the already-reduced field (Swendsen–Wang's
/// <c>BondFormedCount/draws</c> or PKWang's closed form), so the non-commutative
/// accumulate-entropy-per-draw hazard is structurally inexpressible at this tier.</para>
/// </remarks>
public static class AffinityEntropy
{
    /// <summary>
    /// Dispersion entropy <c>H(G)</c> in nats; <see cref="Shannon.EntropyNats(System.ReadOnlySpan{double})"/>
    /// normalizes <see cref="Affinities.G"/> to a distribution before reducing.
    /// </summary>
    public static double EntropyNats(Affinities affinities)
    {
        ArgumentNullException.ThrowIfNull(affinities);
        return EntropyNats(affinities.G);
    }

    public static double EntropyNats(double[] edgeRates)
    {
        ArgumentNullException.ThrowIfNull(edgeRates);
        return Shannon.EntropyNats(edgeRates);
    }
}
