using System;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// Public sampler boundary for one Potts run using a configuration spec.
/// </summary>
/// <remarks>
/// <para><b>Execution boundary.</b> This class is the public entrypoint for a
/// single Potts sampler execution. It constructs the correct internal model
/// specialization, performs optional burn-in and the requested number of
/// SWCycles, and returns the collected observables as a unified result.
/// </para>
///
/// <para><b>Implementation detail.</b> The internal stateful simulation remains
/// hidden inside the generic <see cref="PottsModel{TConfig}"/> engine and the
/// <see cref="ISwEngine"/> bridge. The public API does not expose the hot-
/// path cycle semantics or mutable session state.</para>
/// </remarks>
public static class SwRunner
{
    public static SwRunResult Run(SwRunSpec spec)
    {
        if (spec is null)
            throw new ArgumentNullException(nameof(spec));
        if (spec.Graph.NodeCount <= 0)
            throw new ArgumentException("Graph must be provided and non-empty.", nameof(spec));
        if (spec.Temperature <= 0)
            throw new ArgumentOutOfRangeException(nameof(spec.Temperature), "Temperature must be positive.");
        if (spec.Q < 2)
            throw new ArgumentOutOfRangeException(nameof(spec.Q), "Q must be at least 2.");
        if (spec.Budget.BurnIn < 0)
            throw new ArgumentOutOfRangeException(nameof(spec), "Budget.BurnIn cannot be negative.");
        if (spec.Budget.Cycles < 0)
            throw new ArgumentOutOfRangeException(nameof(spec), "Budget.Cycles cannot be negative.");

        var model = CreateModel(
            spec.Graph,
            spec.Temperature,
            spec.Q,
            spec.Accumulation,
            spec.Seed);

        if (spec.Budget.BurnIn > 0)
            model.BurnIn(spec.Budget.BurnIn);

        model.Run(spec.Budget.Cycles);

        var accumulator = model.GetCheckpoint() with { ReplicaIndex = spec.ReplicaIndex };

        return new SwRunResult
        {
            Accumulator = accumulator,
        };
    }

    private static ISwEngine CreateModel(
        CsrGraph graph,
        double temperature,
        int q,
        AccumulationSpec accumulation,
        int? seed)
        // The two per-edge currency gates select the monomorphized specialization;
        // the per-node landscapes and co-membership ride as runtime flags on the spec
        // (cheap O(N) or O(E) post-passes — they don't earn a JIT gate).
        => (accumulation.Affinities, accumulation.Alignments) switch
        {
            (false, false) => new SwendsenWang<NoCurrencies>(   graph, temperature, q, accumulation, seed),
            (true,  false) => new SwendsenWang<AffinitiesOnly>( graph, temperature, q, accumulation, seed),
            (false, true)  => new SwendsenWang<AlignmentsOnly>( graph, temperature, q, accumulation, seed),
            (true,  true)  => new SwendsenWang<BothCurrencies>( graph, temperature, q, accumulation, seed),
        };
}
