using System;
using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Dendrograms;

/// <summary>
/// One cut's <see cref="ClosureGap"/>, in units where a false merge costs θ and a false split 1 − θ.
/// </summary>
/// <param name="Theta">The cut, on the co-membership scale.</param>
/// <param name="Pairs">The pairs scored: the graph's edges, or every pair.</param>
/// <param name="ClusterCount">The clusters of the closure.</param>
/// <param name="PairRuleLoss">The loss of deciding each pair alone, which no partition can beat.</param>
/// <param name="ClosureLoss">The loss of the closure's partition.</param>
/// <param name="FalseMerges">Pairs the closure joins against their own rule.</param>
/// <param name="FalseMergeCost">The part of the gap those pairs carry.</param>
/// <param name="FalseSplits">Pairs the closure parts against their own rule.</param>
/// <param name="FalseSplitCost">The part of the gap those pairs carry.</param>
public readonly record struct ClosureGapResult(
    double Theta,
    long Pairs,
    int ClusterCount,
    double PairRuleLoss,
    double ClosureLoss,
    long FalseMerges,
    double FalseMergeCost,
    long FalseSplits,
    double FalseSplitCost)
{
    /// <summary>ClosureLoss − PairRuleLoss: the most any partition could save over the closure.</summary>
    public double Gap => FalseMergeCost + FalseSplitCost;

    /// <summary>The gap as a share of the closure's loss; zero when the closure has no loss.</summary>
    public double GapShare => ClosureLoss > 0.0 ? Gap / ClosureLoss : 0.0;
}

/// <summary>
/// The price of closing a co-membership field into a partition, an observable of the flat single-linkage cut
/// <see cref="AffinityThreshold.Connect"/>. With costs a for a false merge and b for a false split (Binder's 1978 loss),
/// each pair alone links when ⟨n_ij⟩ &gt; θ = a/(a + b), so the θ-graph is the best decision when transitivity isn't
/// required, and <c>L(θ-graph) ≤ L(best partition) ≤ L(closure)</c>.
/// The gap <c>Σ |⟨n_ij⟩ − θ|</c>, over the pairs where the closure contradicts their own rule, bounds what any search
/// for a better partition could gain: at zero the closure is the exact optimum.
/// </summary>
/// <remarks>
/// <para><b>Two scopes.</b> Over the graph's edges only false merges can occur, and a pure chain reads zero, since
/// every edge on it passes. Over all pairs (<see cref="DenseCoMembership"/>) the chain's far pairs count.</para>
///
/// <para><b>Sampling noise.</b> Reading the rule and the costs off one estimate inflates the gap: a pair near θ whose
/// noise lands it on the wrong side always adds. Decide on one estimate and evaluate on another, independent one, and
/// that bias goes. The gap then measures the decisions made and can come out negative.</para>
///
/// <para><b>Integrated over θ.</b> The gap's integral over θ is <c>½ Σ (m_ij − ⟨n_ij⟩)²</c>, with m_ij the level at which
/// the closure first joins the pair: the single-linkage tree's distortion of the field, reduced to one number per
/// temperature.</para>
/// </remarks>
public static class ClosureGap
{
    /// <summary>The gap over the graph's edges, at <paramref name="theta"/> on the co-membership scale.</summary>
    /// <param name="graph">The graph the field lives on.</param>
    /// <param name="coMembership">The per-edge column deciding the closure and each pair's rule, by CSR slot.</param>
    /// <param name="theta">The cut, in [0, 1].</param>
    /// <param name="evaluate">An independent column to read the costs from; the deciding column when null.</param>
    public static ClosureGapResult OverEdges(CsrGraph graph, double[] coMembership, double theta, double[]? evaluate = null)
    {
        ArgumentNullException.ThrowIfNull(coMembership);
        double[] costs = evaluate ?? coMembership;
        RequireColumn(graph, coMembership, nameof(coMembership));
        RequireColumn(graph, costs, nameof(evaluate));

        Assignment closure = AffinityThreshold.Connect(graph, coMembership, theta);
        var tally = new Tally(theta);
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            bool joined = closure.Labels[edge.Source] == closure.Labels[edge.Target];
            tally.Add(joined, coMembership[edge.Slot], costs[edge.Slot]);
        }
        return tally.Result(closure.Count);
    }

    /// <summary>
    /// The gap over all pairs, at <paramref name="theta"/> on the co-membership scale. The closure runs on the graph's
    /// edges, as <see cref="AffinityThreshold.Connect"/> does, and every pair is scored against it.
    /// </summary>
    /// <param name="graph">The graph the closure runs on.</param>
    /// <param name="coMembership">The field deciding the closure and each pair's rule.</param>
    /// <param name="theta">The cut, in [0, 1].</param>
    /// <param name="evaluate">An independent field to read the costs from; the deciding field when null.</param>
    public static ClosureGapResult OverPairs(CsrGraph graph, DenseCoMembership coMembership, double theta, DenseCoMembership? evaluate = null)
    {
        ArgumentNullException.ThrowIfNull(coMembership);
        DenseCoMembership costs = evaluate ?? coMembership;
        if (costs.NodeCount != coMembership.NodeCount)
            throw new ArgumentException($"Node counts differ: {coMembership.NodeCount} and {costs.NodeCount}.", nameof(evaluate));

        Assignment closure = AffinityThreshold.Connect(graph, coMembership.EdgeColumn(graph), theta);
        int[] labels = closure.Labels;
        var tally = new Tally(theta);
        for (int i = 0; i < labels.Length; i++)
            for (int j = i + 1; j < labels.Length; j++)
                tally.Add(labels[i] == labels[j], coMembership.At(i, j), costs.At(i, j));
        return tally.Result(closure.Count);
    }

    /// <summary>The gap over the graph's edges, integrated over θ: <c>½ Σ_e (m_e − ⟨n_e⟩)²</c>.</summary>
    public static double IntegratedOverEdges(CsrGraph graph, double[] coMembership)
    {
        ArgumentNullException.ThrowIfNull(coMembership);
        RequireColumn(graph, coMembership, nameof(coMembership));

        var tree = new MergeTree(graph, coMembership);
        double total = 0.0;
        foreach (UndirectedEdge edge in graph.UndirectedEdges())
        {
            double d = tree.JoinLevel(edge.Source, edge.Target) - coMembership[edge.Slot];
            total += d * d;
        }
        return 0.5 * total;
    }

    /// <summary>The gap over all pairs, integrated over θ: <c>½ Σ_ij (m_ij − ⟨n_ij⟩)²</c>.</summary>
    public static double IntegratedOverPairs(CsrGraph graph, DenseCoMembership coMembership)
    {
        ArgumentNullException.ThrowIfNull(coMembership);

        var tree = new MergeTree(graph, coMembership.EdgeColumn(graph));
        double total = 0.0;
        for (int i = 0; i < graph.NodeCount; i++)
            for (int j = i + 1; j < graph.NodeCount; j++)
            {
                double d = tree.JoinLevel(i, j) - coMembership.At(i, j);
                total += d * d;
            }
        return 0.5 * total;
    }

    private static void RequireColumn(CsrGraph graph, double[] column, string name)
    {
        if (column.Length != graph.Targets.Length)
            throw new ArgumentException($"A column of {column.Length} does not match {graph.Targets.Length} CSR slots.", name);
    }

    private sealed class Tally(double theta)
    {
        private long _pairs, _falseMerges, _falseSplits;
        private double _pairRuleLoss, _closureLoss, _falseMergeCost, _falseSplitCost;

        public void Add(bool joined, double decide, double evaluate)
        {
            bool link = decide > theta;
            double costIfJoined = theta * (1.0 - evaluate);
            double costIfApart = (1.0 - theta) * evaluate;

            _pairs++;
            _pairRuleLoss += link ? costIfJoined : costIfApart;
            _closureLoss += joined ? costIfJoined : costIfApart;
            if (joined && !link)
            {
                _falseMerges++;
                _falseMergeCost += theta - evaluate;
            }
            else if (!joined && link)
            {
                _falseSplits++;
                _falseSplitCost += evaluate - theta;
            }
        }

        public ClosureGapResult Result(int clusterCount) => new(
            theta, _pairs, clusterCount, _pairRuleLoss, _closureLoss, _falseMerges, _falseMergeCost, _falseSplits, _falseSplitCost);
    }

    /// <summary>
    /// The single-linkage merge tree of a co-membership column, answering at which level the closure first joins a pair
    /// (the widest path's bottleneck). Union by size without path compression keeps each node's links descending toward
    /// the root, so a query walks both nodes up to where they meet, in O(log N) steps.
    /// </summary>
    private sealed class MergeTree
    {
        private readonly int[] _parent, _order;   // _order: when the node's link to its parent formed; roots never did
        private readonly double[] _level;

        public MergeTree(CsrGraph graph, double[] column)
        {
            int n = graph.NodeCount;
            _parent = new int[n];
            _order = new int[n];
            _level = new double[n];
            var size = new int[n];
            for (int i = 0; i < n; i++)
            {
                _parent[i] = i;
                _order[i] = int.MaxValue;
                size[i] = 1;
            }

            int edges = 0;
            foreach (UndirectedEdge _ in graph.UndirectedEdges()) edges++;
            var sources = new int[edges];
            var targets = new int[edges];
            var values = new double[edges];
            var byValue = new int[edges];
            int k = 0;
            foreach (UndirectedEdge edge in graph.UndirectedEdges())
            {
                sources[k] = edge.Source;
                targets[k] = edge.Target;
                values[k] = column[edge.Slot];
                byValue[k] = k;
                k++;
            }
            Array.Sort(byValue, (a, b) => values[b].CompareTo(values[a]));

            int order = 0;
            foreach (int e in byValue)
            {
                if (values[e] <= 0.0) break;    // no θ ≥ 0 passes it, so it joins nothing
                int a = Root(sources[e]), b = Root(targets[e]);
                if (a == b) continue;
                if (size[a] > size[b]) (a, b) = (b, a);
                _parent[a] = b;
                _order[a] = order++;
                _level[a] = values[e];
                size[b] += size[a];
            }
        }

        /// <summary>The level at which the closure first joins the pair; zero if no level does.</summary>
        public double JoinLevel(int i, int j)
        {
            double level = 0.0;
            while (i != j)
            {
                // The earlier link is the one lower in the tree, so taking it never climbs past the meeting point.
                int next = _order[i] < _order[j] ? i : j;
                if (_order[next] == int.MaxValue) return 0.0;   // two roots: never joined
                level = _level[next];
                if (next == i) i = _parent[i];
                else j = _parent[j];
            }
            return level;
        }

        private int Root(int i)
        {
            while (_parent[i] != i) i = _parent[i];
            return i;
        }
    }
}
