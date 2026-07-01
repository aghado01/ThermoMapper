using System;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions.Strategies;

/// <summary>
/// Threshold-and-connected-components partition keyed on the BWD1996 eq-4
/// discriminant δ̄_ij = ((q−1)⟨n_ij⟩+1)/q. Uses the improved co-membership
/// estimator ⟨n_ij⟩ — P(same cluster), capturing transitive co-clustering via
/// multi-hop paths — rather than the raw bond-frequency or spin-agreement proxies.
/// </summary>
/// <remarks>
/// <para><b>Algorithm.</b> Transform each co-membership rate ⟨n_ij⟩ to the
/// normalized discriminant δ̄_ij = ((q−1)⟨n_ij⟩+1)/q, then threshold at
/// <see cref="Theta"/> and run connected components via union-find.</para>
///
/// <para><b>Why eq-4.</b> The raw ⟨n_ij⟩ ranges from 1/q (pure noise) to 1
/// (always co-clustered); δ̄_ij shifts the range to [1/q, 1] but normalizes
/// the noise floor to 1/q exactly, making the threshold interpretation
/// independent of q. At q=20, θ=0.5 thresholds on ⟨n_ij⟩ &gt; 9/19 ≈ 0.474.</para>
///
/// <para><b>Requires.</b> The <see cref="CoMembership"/> currency with Q≥2;
/// configure the sweep with
/// <see cref="Clustering.Graphical.SPC.Runtime.Core.Sampler.AccumulationSpec.Currencies"/>.</para>
/// </remarks>
public sealed class ThresholdCoMembership : IPartitionStrategy
{
    /// <summary>
    /// Cut threshold on δ̄_ij; edge included iff δ̄_ij &gt; Theta.
    /// Defaults to 0.5.
    /// </summary>
    public double Theta { get; init; } = 0.5;

    /// <summary>
    /// When <see langword="true"/>, also unions each node with its single
    /// highest-δ̄ neighbor regardless of <see cref="Theta"/> (Domany1999
    /// step 2 — peripheral capture). Default off.
    /// </summary>
    public bool PeripheralCapture { get; init; } = false;

    public Assignment Apply(CsrGraph graph, Affinities affinities, Alignments? alignments, CoMembership? coMembership = null)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(affinities);
        if (coMembership is null)
            throw new InvalidOperationException(
                $"{nameof(ThresholdCoMembership)} requires the CoMembership currency — " +
                "configure the sweep with AccumulationSpec.Currencies or pass --accumulation comembership.");
        if (coMembership.G.Length != graph.Targets.Length)
            throw new InvalidOperationException(
                $"CoMembership.G length ({coMembership.G.Length}) does not match " +
                $"CSR slot count ({graph.Targets.Length}).");
        if (Theta < 0.0 || Theta > 1.0)
            throw new InvalidOperationException($"Theta ({Theta}) must lie in [0, 1].");
        if (coMembership.Q < 2)
            throw new InvalidOperationException($"CoMembership.Q ({coMembership.Q}) must be ≥ 2.");

        // BWD1996 eq-4: δ̄_ij = ((q−1)·⟨n_ij⟩ + 1) / q
        double q = coMembership.Q;
        double[] g = coMembership.G;
        double[] delta = new double[g.Length];
        for (int e = 0; e < g.Length; e++)
            delta[e] = ((q - 1.0) * g[e] + 1.0) / q;

        return AffinityThreshold.Connect(graph, delta, Theta, PeripheralCapture);
    }
}
