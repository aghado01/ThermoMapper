using System;
using System.Collections.Generic;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Dendrograms;

/// <summary>
/// The shared substrate of the two classical SPC T-stack resolvers: cut the
/// partition at <b>every</b> sweep temperature (not just the Blatt
/// phase-midpoints) from a rich sweep's pooled per-T co-membership columns,
/// yielding a <b>dense</b> <see cref="PartitionHierarchy"/> — one level per
/// grid temperature. The structure-side reading
/// SPC's <c>BlattPartitionStrategy</c> only sampled at phase representatives;
/// this completes it to the canonical dendrogram-across-T reading.
/// </summary>
/// <remarks>
/// <para><b>Schedule-agnostic.</b> The grid is whatever the frames carry —
/// ASCENDING and possibly NON-UNIFORM (the auto grid is log-spaced; BARS will
/// later feed a sparse+dense grid). Nothing here assumes a constant ΔT;
/// downstream persistence/lifetime reductions weight by the actual cell
/// widths.</para>
///
/// <para><b>The cut.</b> Each grid temperature's pooled eq-4 discriminant
/// column δ̄_ij(T) = ((q−1)⟨n_ij⟩+1)/q (from SPC's <c>SweepEdgeCurves</c>) is
/// thresholded at <paramref name="theta"/> and closed into connected
/// components — the identical threshold-and-connect step the chosen-T
/// SPC's <c>ThresholdCoMembership</c> cut applies, run at every T. Peripheral
/// capture (Domany step 2) is available per-level via the same shared
/// <see cref="AffinityThreshold"/> path.</para>
///
/// <para><b>Nesting.</b> Whether the dense stack is strictly nested is
/// <see cref="PartitionNesting.Holds"/>-checked and reported on the result —
/// the don't-warp gate that decides whether the stack bridges to a
/// <c>Clustering.Dendrograms.Dendrogram</c> (Track 1) or is resolved by
/// overlap-linked lineage persistence (Track 2, <c>LineagePersistence</c>).</para>
/// </remarks>
public static class DenseTStack
{
    /// <summary>
    /// Producer-agnostic core: cut a dense stack from explicit ascending
    /// temperatures and their grid-major per-edge δ̄ columns (CSR-slot indexed).
    /// SW supplies pooled sampled columns (SPC's <c>SweepDenseTStack</c>); a solver
    /// supplies its closed-form columns; BARS substitutes monotonized posterior
    /// columns upstream — all share this cut.
    /// </summary>
    public static PartitionHierarchy Build(
        CsrGraph graph,
        IReadOnlyList<double> temperatures,
        IReadOnlyList<double[]> deltaByGridPoint,
        double theta = 0.5,
        bool peripheralCapture = false)
    {
        ArgumentNullException.ThrowIfNull(graph);
        ArgumentNullException.ThrowIfNull(temperatures);
        ArgumentNullException.ThrowIfNull(deltaByGridPoint);
        if (temperatures.Count == 0)
            throw new ArgumentException("At least one grid temperature is required.", nameof(temperatures));
        if (deltaByGridPoint.Count != temperatures.Count)
            throw new ArgumentException(
                $"One δ̄ column per grid temperature: {deltaByGridPoint.Count} columns vs " +
                $"{temperatures.Count} temperatures.",
                nameof(deltaByGridPoint));
        for (int t = 1; t < temperatures.Count; t++)
            if (temperatures[t] <= temperatures[t - 1])
                throw new ArgumentException("Grid temperatures must be strictly ascending.", nameof(temperatures));

        int slots = graph.Targets.Length;
        var levels = new List<HierarchyLevel>(temperatures.Count);
        for (int t = 0; t < temperatures.Count; t++)
        {
            double[] column = deltaByGridPoint[t];
            if (column.Length != slots)
                throw new ArgumentException(
                    $"δ̄ column {t} length ({column.Length}) does not match CSR slot count ({slots}).",
                    nameof(deltaByGridPoint));

            Assignment partition = AffinityThreshold.Connect(graph, column, theta, peripheralCapture);
            levels.Add(new HierarchyLevel(
                Temperature: temperatures[t],
                Partition:   partition,
                Provenance:  $"dense T-stack cut: δ̄ > {theta} @ T={temperatures[t]:G6}" +
                             (peripheralCapture ? " (+ peripheral capture)" : "")));
        }

        bool nestingHolds = PartitionNesting.Holds(levels);
        return new PartitionHierarchy(levels, NestingHolds: nestingHolds);
    }
}
