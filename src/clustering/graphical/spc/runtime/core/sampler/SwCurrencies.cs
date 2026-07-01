using System;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// Mints the per-edge currencies from a Swendsen–Wang <see cref="Accumulator"/> — the 1-form
/// reductions that collapse the SW edge counts to the model-agnostic currencies the graph-tier
/// observables consume. <c>BondFormedCount / draws → </c><see cref="Affinities"/> (bond-survival,
/// the universal currency); <c>SpinAgreementCount / draws → </c><see cref="Alignments"/>
/// (co-association, SW-native). A solver (PKWang) writes <see cref="Affinities"/> directly from the
/// closed-form survival kernel and never lands here; having no spin ensemble, it mints no
/// <see cref="Alignments"/> at all.
/// </summary>
/// <remarks>
/// These exist only where the SW pass materialized the per-edge arrays
/// (i.e. <see cref="ISwConfig.Affinities"/> or <see cref="ISwConfig.Alignments"/> was true);
/// on a <see cref="AccumulationSpec.None"/> accumulator the source arrays are null
/// and the mint throws. The currencies are the narrow waist between the sampler and the graph-tier
/// consumers — once minted, provenance (which sampler, how many draws) is erased.
/// </remarks>
public static class SwCurrencies
{
    /// <summary>
    /// Reduce the frozen-bond counts to the <see cref="Affinities"/> currency:
    /// <c>G[e] = BondFormedCount[e] / DrawCount</c>.
    /// </summary>
    public static Affinities ToAffinities(Accumulator accumulator)
    {
        ArgumentNullException.ThrowIfNull(accumulator);
        int[] counts = accumulator.BondFormedCount
            ?? throw new InvalidOperationException(
                "Accumulator carries no BondFormedCount — Affinities requires an edge-observable SW pass.");

        return new Affinities
        {
            Temperature  = accumulator.Temperature,
            G            = PerEdgeRate(counts, accumulator.DrawCount),
            ReplicaIndex = accumulator.ReplicaIndex,
        };
    }

    /// <summary>
    /// Reduce the spin-agreement counts to the <see cref="Alignments"/> currency:
    /// <c>G[e] = SpinAgreementCount[e] / DrawCount</c>.
    /// </summary>
    public static Alignments ToAlignments(Accumulator accumulator)
    {
        ArgumentNullException.ThrowIfNull(accumulator);
        int[] counts = accumulator.SpinAgreementCount
            ?? throw new InvalidOperationException(
                "Accumulator carries no SpinAgreementCount — Alignments requires an edge-observable SW pass.");

        return new Alignments
        {
            Temperature  = accumulator.Temperature,
            G            = PerEdgeRate(counts, accumulator.DrawCount),
            ReplicaIndex = accumulator.ReplicaIndex,
        };
    }

    /// <summary>
    /// Reduce the co-membership counts to the <see cref="CoMembership"/> currency:
    /// <c>G[e] = CoMembershipCount[e] / DrawCount</c> — fraction of draws where
    /// both endpoints fell in the same bond cluster (⟨n_ij⟩, improved estimator).
    /// </summary>
    public static CoMembership ToCoMembership(Accumulator accumulator)
    {
        ArgumentNullException.ThrowIfNull(accumulator);
        int[] counts = accumulator.CoMembershipCount
            ?? throw new InvalidOperationException(
                "Accumulator carries no CoMembershipCount — CoMembership requires an accumulation run with CoMembership=true.");

        return new CoMembership
        {
            Temperature  = accumulator.Temperature,
            Q            = accumulator.Q,
            G            = PerEdgeRate(counts, accumulator.DrawCount),
            ReplicaIndex = accumulator.ReplicaIndex,
        };
    }

    private static double[] PerEdgeRate(int[] counts, int draws)
    {
        if (draws <= 0)
            throw new InvalidOperationException(
                $"DrawCount must be positive to mint a currency; was {draws}.");

        double inv = 1.0 / draws;
        var g = new double[counts.Length];
        for (int e = 0; e < counts.Length; e++)
            g[e] = counts[e] * inv;
        return g;
    }
}
