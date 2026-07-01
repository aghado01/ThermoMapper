using System;
using Clustering.Graphical.SPC.Partitions;
using Clustering.Graphical.SPC.Partitions.Strategies;
using Clustering.Primitives;
using Graphs;
using Graphs.Models.Potts;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Core.Solver;

/// <summary>
/// Pan Ke Wang's 2020 direct (non-Monte-Carlo) SPC solver. Computes the
/// closed-form bond-activity affinity <c>G = 1 − exp(−Hcum/T)</c> — the
/// <c>M→∞</c> limit of the paper's inverse-transform draws — instead of
/// averaging M samples. Exact, zero-variance, and the clustering it induces is
/// thermal single-linkage (Lemma B).
/// </summary>
/// <remarks>
/// Credit: Wang et al. 2020 for the mean-field construction. This computes the
/// closed form the paper left on the table (it kept the inherited Monte-Carlo
/// wrapper). Reduction and prior-art live in
/// <c>.discussion/issues/spc-samplers/</c>, not here.
/// </remarks>
public static class PKWang
{
    /// <summary>
    /// Build the temperature-independent <c>Hcum</c> ladder once for the chosen
    /// field. <paramref name="rule"/> is consumed only by directed
    /// fields (LocalField); MeanField ignores it. The returned context is reused
    /// across an entire temperature sweep.
    /// </summary>
    public static PKWangContext Prepare(
        CsrGraph graph,
        EdgeWeightKind weightKind,
        Field field,
        SymmetrizationRule rule = SymmetrizationRule.Mutual)
    {
        if (graph.NodeCount <= 0)
            throw new ArgumentException("Graph must be non-empty.", nameof(graph));
        if (graph.Targets is null || graph.Weights is null || graph.RowPointers is null)
            throw new ArgumentException("Graph CSR arrays are not initialized.", nameof(graph));
        if (weightKind != EdgeWeightKind.Coupling)
            throw new ArgumentException(
                $"PKWang's survival kernel consumes couplings; got {weightKind}-weighted edges " +
                "(a Distance-weighted graph would silently feed distances to the kernel). " +
                "Feed a Coupling-weighted graph (GraphBuildResult.WeightKind == Coupling).",
                nameof(weightKind));

        return field switch
        {
            Field.Mean => Build<MeanField>(graph, rule),
            Field.Local => Build<LocalField>(graph, rule),
            _ => throw new ArgumentOutOfRangeException(nameof(field), field, "Unsupported field."),
        };
    }

    private static PKWangContext Build<TField>(CsrGraph graph, SymmetrizationRule rule)
        where TField : struct, IField
    {
        double[] hcum = TField.BuildHcum(graph);
        int[]? mirror = TField.DirectedSymmetrize ? graph.BuildReverseSlotMap() : null;
        return new PKWangContext(graph, hcum, TField.DirectedSymmetrize, rule, mirror);
    }

    /// <summary>
    /// The solver-native temperature bracket, read off its OWN cumulative-energy
    /// ladder. For the closed form an edge's cut sits at <c>T = Hcum/ln2</c>, so the
    /// temperatures that change anything span the range of <c>Hcum / ln2</c>: below
    /// the cold end every edge is active (one cluster), above the hot end every edge
    /// is cut (singletons). q-free and gauge-free — no borrowed Potts <c>T_ps(q)</c>
    /// (the solver has no susceptibility transition to anchor on; T is a
    /// single-linkage cut height, not a phase-transition temperature).
    /// </summary>
    /// <remarks>
    /// Robust to the tails: the absolute min Hcum is the single weakest (near-zero)
    /// edge, whose cut temperature is far below any real structure — using it would
    /// waste the whole log grid on an empty cold tail. So the endpoints are taken at
    /// the <paramref name="loQuantile"/>/<paramref name="hiQuantile"/> of the positive
    /// cumulative energies, then widened by <paramref name="coldPad"/>/<paramref name="hotPad"/>
    /// so the fully-merged and fully-fragmented ends are still reached.
    /// </remarks>
    public static (double Lo, double Hi) EstimateBracket(
        CsrGraph graph,
        EdgeWeightKind weightKind,
        Field field,
        SymmetrizationRule rule = SymmetrizationRule.Mutual,
        double loQuantile = 0.02,
        double hiQuantile = 0.98,
        double coldPad = 0.5,
        double hotPad = 2.0)
    {
        if (loQuantile < 0.0 || loQuantile >= hiQuantile || hiQuantile > 1.0)
            throw new ArgumentOutOfRangeException(nameof(loQuantile), "Require 0 ≤ loQuantile < hiQuantile ≤ 1.");
        if (coldPad <= 0.0 || coldPad >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(coldPad), "coldPad must lie in (0, 1).");
        if (hotPad <= 1.0)
            throw new ArgumentOutOfRangeException(nameof(hotPad), "hotPad must exceed 1.");

        double[] hcum = Prepare(graph, weightKind, field, rule).Hcum;
        int count = 0;
        for (int e = 0; e < hcum.Length; e++)
            if (hcum[e] > 0.0) count++;
        if (count == 0)
            throw new InvalidOperationException(
                "Graph carries no positive cumulative energy; cannot bracket the sweep.");

        var positive = new double[count];
        int k = 0;
        for (int e = 0; e < hcum.Length; e++)
            if (hcum[e] > 0.0) positive[k++] = hcum[e];
        Array.Sort(positive);

        double ln2 = Math.Log(2.0);
        double lo = positive[QuantileIndex(count, loQuantile)] / ln2;
        double hi = positive[QuantileIndex(count, hiQuantile)] / ln2;
        return (coldPad * lo, hotPad * hi);

        static int QuantileIndex(int n, double q) => Math.Clamp((int)(q * (n - 1)), 0, n - 1);
    }

    /// <summary>
    /// Apply the closed-form survival kernel at temperature <paramref name="T"/>,
    /// producing the per-edge affinity. No draws, no variance. For directed
    /// fields the two directions are reconciled per the context's
    /// <see cref="SymmetrizationRule"/>, written into the <c>j &gt; i</c> slot.
    /// </summary>
    public static Affinities Solve(PKWangContext context, double T, int replicaIndex = 0)
    {
        ArgumentNullException.ThrowIfNull(context);
        if (T <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(T), "Temperature must be positive.");

        double[] hcum = context.Hcum;
        var g = new double[hcum.Length];
        for (int e = 0; e < hcum.Length; e++)
        {
            double h = hcum[e];
            if (h > 0.0) g[e] = FkKernel.BondProbability(h, T);
        }

        if (context.DirectedSymmetrize)
            EdgeFieldSymmetrization.Symmetrize(context.Graph, g, context.Mirror!, context.Rule);

        return new Affinities { Temperature = T, G = g, ReplicaIndex = replicaIndex };
    }

    /// <summary>
    /// Solve then threshold-and-connect: union edges with <c>G &gt; theta</c>
    /// (default ½, i.e. the <c>Hcum &gt; T·ln2</c> cut) and densify the
    /// components. One temperature → one partition; the thermomapper backend
    /// sweeps <c>T</c> over this.
    /// </summary>
    public static Assignment Cluster(PKWangContext context, double T, double theta = 0.5)
    {
        Affinities affinity = Solve(context, T);
        return AffinityThreshold.Connect(context.Graph, affinity.G, theta);
    }

}
