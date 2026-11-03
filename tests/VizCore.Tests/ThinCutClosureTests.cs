using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Rng;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The thin-cut closure against what can be counted by hand, and against brute force over every subset on small seeded
/// graphs from near-trees to dense ones.
/// </summary>
public sealed class ThinCutClosureTests
{
    private const int SweepGraphs = 300;

    private readonly record struct Edges(int NodeCount, int[] Sources, int[] Targets, double[] Values)
    {
        internal ThinCutClosure Closure(int maxCut) => new(NodeCount, Sources, Targets, Values, maxCut);
    }

    private static Edges Of(int nodeCount, params (int S, int T, double Value)[] edges) =>
        new(nodeCount, edges.Select(e => e.S).ToArray(), edges.Select(e => e.T).ToArray(), edges.Select(e => e.Value).ToArray());

    private static IEnumerable<(int S, int T, double Value)> Clique(int from, int count, double value)
    {
        for (int i = from; i < from + count; i++)
            for (int j = i + 1; j < from + count; j++)
                yield return (i, j, value);
    }

    private static IEnumerable<(int S, int T, double Value)> Cycle(int from, int count, double value)
    {
        for (int i = 0; i < count; i++)
            yield return (from + i, from + (i + 1) % count, value);
    }

    private static Edges RandomGraph(Xoshiro256PlusPlus rng, int maxNodes)
    {
        int n = 2 + rng.NextInt(maxNodes - 1);
        double density = rng.NextDouble();
        density *= density;
        var edges = new List<(int, int, double)>();
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (rng.NextDouble() < density)
                    edges.Add((i, j, rng.NextDouble()));
        return Of(n, edges.ToArray());
    }

    private static long Crossing(Edges edges, Func<int, bool> inside)
    {
        long weight = 0;
        for (int e = 0; e < edges.Sources.Length; e++)
            if (inside(edges.Sources[e]) != inside(edges.Targets[e]))
                weight += ThinCutClosure.EdgeWeight(edges.Values[e]);
        return weight;
    }

    // The lightest cut between a and b over every subset, and how many edges it holds.
    private static long BruteMinCut(Edges edges, int a, int b)
    {
        long best = long.MaxValue;
        for (int mask = 0; mask < 1 << edges.NodeCount; mask++)
            if ((mask >> a & 1) == 1 && (mask >> b & 1) == 0)
                best = Math.Min(best, Crossing(edges, v => (mask >> v & 1) == 1));
        return best;
    }

    private static bool Together(int[] labels, params int[] nodes) => nodes.All(v => labels[v] == labels[nodes[0]]);

    [Fact]
    public void ZeroCut_IsClosure()
    {
        Edges edges = Of(7, (0, 1, 0.9), (1, 2, 0.6), (3, 4, 0.7), (5, 6, 0.8), (6, 5, 0.1));
        int[] labels = edges.Closure(6).Partition(0, 1);
        Assert.True(Together(labels, 0, 1, 2));
        Assert.True(Together(labels, 3, 4));
        Assert.True(Together(labels, 5, 6));
        Assert.Equal(3, labels.Distinct().Count());

        var rng = new Xoshiro256PlusPlus(20261102);
        for (int g = 0; g < SweepGraphs; g++)
        {
            Edges random = RandomGraph(rng, 12);
            int[] set = Enumerable.Range(0, random.NodeCount).ToArray();
            int Find(int x) => set[x] == x ? x : set[x] = Find(set[x]);
            for (int e = 0; e < random.Sources.Length; e++) set[Find(random.Sources[e])] = Find(random.Targets[e]);

            int[] closure = random.Closure(3).Partition(0, 5);
            for (int i = 0; i < random.NodeCount; i++)
                for (int j = i + 1; j < random.NodeCount; j++)
                    Assert.Equal(Find(i) == Find(j), closure[i] == closure[j]);
        }
    }

    [Fact]
    public void Dumbbell_PartsOnlyWhenTheCutIsThinAndBothSidesLarge()
    {
        (int, int, double)[] bars = Clique(0, 6, 0.9).Concat(Clique(6, 6, 0.9)).ToArray();

        Edges one = Of(12, bars.Append((5, 6, 0.5)).ToArray());
        Assert.True(Together(one.Closure(6).Partition(0, 6), Enumerable.Range(0, 12).ToArray()));
        int[] parted = one.Closure(6).Partition(1, 6);
        Assert.True(Together(parted, 0, 1, 2, 3, 4, 5));
        Assert.True(Together(parted, 6, 7, 8, 9, 10, 11));
        Assert.NotEqual(parted[0], parted[6]);
        Assert.True(Together(one.Closure(6).Partition(1, 7), Enumerable.Range(0, 12).ToArray()));

        Edges three = Of(12, bars.Concat([(0, 6, 0.5), (1, 7, 0.5), (2, 8, 0.5)]).ToArray());
        Assert.True(Together(three.Closure(6).Partition(2, 6), Enumerable.Range(0, 12).ToArray()));
        Assert.NotEqual(three.Closure(6).Partition(3, 6)[0], three.Closure(6).Partition(3, 6)[6]);
    }

    [Fact]
    public void Pendants_StayWithTheirBlock()
    {
        // Node 12 hangs off one bar and 13 off the other, each through a cut of one edge with one point behind it.
        Edges edges = Of(14, Clique(0, 6, 0.9).Concat(Clique(6, 6, 0.9)).Concat([(5, 6, 0.5), (0, 12, 0.6), (11, 13, 0.6)]).ToArray());
        int[] labels = edges.Closure(6).Partition(1, 5);
        Assert.True(Together(labels, 0, 1, 2, 3, 4, 5, 12));
        Assert.True(Together(labels, 6, 7, 8, 9, 10, 11, 13));
        Assert.Equal(2, labels.Distinct().Count());

        // With no size rule both pendants stand alone.
        Assert.Equal(4, edges.Closure(6).Partition(1, 1).Distinct().Count());
    }

    [Fact]
    public void ChainsPartAtTheirWeakestLink()
    {
        // A bar, a node, a bar. Both links are one-edge cuts with two large sides, until the first is removed.
        Edges edges = Of(13, Clique(0, 6, 0.9).Concat(Clique(7, 6, 0.9)).Concat([(5, 6, 0.8), (6, 7, 0.55)]).ToArray());
        int[] labels = edges.Closure(6).Partition(1, 6);
        Assert.True(Together(labels, 0, 1, 2, 3, 4, 5, 6));
        Assert.True(Together(labels, 7, 8, 9, 10, 11, 12));
        Assert.NotEqual(labels[0], labels[7]);

        Edges swapped = Of(13, Clique(0, 6, 0.9).Concat(Clique(7, 6, 0.9)).Concat([(5, 6, 0.55), (6, 7, 0.8)]).ToArray());
        int[] other = swapped.Closure(6).Partition(1, 6);
        Assert.True(Together(other, 6, 7, 8, 9, 10, 11, 12));
        Assert.NotEqual(other[0], other[6]);
    }

    [Fact]
    public void AmongEqualCounts_TheWeakestCutIsTheOneRead()
    {
        // Two long cycles joined by two weak edges. Any two edges of a cycle also cut two, so counts alone can't tell.
        Edges edges = Of(40, Cycle(0, 20, 0.9).Concat(Cycle(20, 20, 0.9)).Concat([(3, 23, 0.5), (12, 31, 0.5)]).ToArray());
        int[] labels = edges.Closure(6).Partition(2, 15);
        Assert.True(Together(labels, Enumerable.Range(0, 20).ToArray()));
        Assert.True(Together(labels, Enumerable.Range(20, 20).ToArray()));
        Assert.NotEqual(labels[0], labels[20]);
    }

    [Fact]
    public void Sweep_EveryTreeEdgeIsARealCutAndPathsHoldTheMinimum()
    {
        var rng = new Xoshiro256PlusPlus(20261103);
        for (int g = 0; g < SweepGraphs; g++)
        {
            Edges edges = RandomGraph(rng, 9);
            int maxCut = rng.NextInt(5);
            ThinCutClosure closure = edges.Closure(maxCut);
            (int[] parent, int[] weight) = closure.Tree;

            int[] Path(int v)
            {
                var path = new List<int> { v };
                while (parent[v] != v) path.Add(v = parent[v]);
                return path.ToArray();
            }

            for (int v = 0; v < edges.NodeCount; v++)
            {
                if (parent[v] == v || weight[v] >= closure.Thick) continue;
                int child = v;
                Assert.Equal(weight[v], Crossing(edges, x => Path(x).Contains(child)));
            }

            for (int a = 0; a < edges.NodeCount; a++)
                for (int b = a + 1; b < edges.NodeCount; b++)
                {
                    int[] up = Path(a), down = Path(b);
                    long truth = BruteMinCut(edges, a, b);
                    if (up[^1] != down[^1])
                    {
                        Assert.Equal(0, truth);
                        continue;
                    }

                    int meet = up.First(down.Contains);
                    long least = long.MaxValue;
                    foreach (int x in up.TakeWhile(x => x != meet).Concat(down.TakeWhile(x => x != meet)))
                        least = Math.Min(least, weight[x]);
                    Assert.Equal(Math.Min(truth, closure.Thick), least);
                }
        }
    }

    [Fact]
    public void ClassMembersWithNoPathInsideTheirClass_StandApart()
    {
        // Three paths join 0 and 1, through 2, 3 and 4. No cut of two edges parts them, but each path leaves their class.
        Edges edges = Of(5, (0, 2, 0.9), (2, 1, 0.9), (0, 3, 0.9), (3, 1, 0.9), (0, 4, 0.9), (4, 1, 0.9));
        Assert.Single(edges.Closure(6).Partition(1, 1).Distinct());
        Assert.Equal(5, edges.Closure(6).Partition(2, 1).Distinct().Count());
    }

    [Fact]
    public void Sweep_NoSizeRule_GivesTheConnectedPartsOfTheEdgeConnectedClasses()
    {
        const int unit = 1 << 20;
        var rng = new Xoshiro256PlusPlus(20261104);
        for (int g = 0; g < SweepGraphs; g++)
        {
            Edges edges = RandomGraph(rng, 9);
            int n = edges.NodeCount;
            ThinCutClosure closure = edges.Closure(6);

            var counts = new long[n, n];
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    counts[a, b] = counts[b, a] = BruteMinCut(edges, a, b) / unit;

            for (int k = 0; k <= 3; k++)
            {
                int[] set = Enumerable.Range(0, n).ToArray();
                int Find(int x) => set[x] == x ? x : set[x] = Find(set[x]);
                for (int e = 0; e < edges.Sources.Length; e++)
                    if (counts[edges.Sources[e], edges.Targets[e]] > k)
                        set[Find(edges.Sources[e])] = Find(edges.Targets[e]);

                int[] labels = closure.Partition(k, 1);
                for (int a = 0; a < n; a++)
                    for (int b = a + 1; b < n; b++)
                        Assert.Equal(Find(a) == Find(b), labels[a] == labels[b]);
            }
        }
    }

    [Fact]
    public void Sweep_NoSizeRule_OnlyCoarsensAsEdgesAreAdded()
    {
        var rng = new Xoshiro256PlusPlus(20261105);
        for (int g = 0; g < SweepGraphs; g++)
        {
            Edges edges = RandomGraph(rng, 30);
            int keep = rng.NextInt(edges.Sources.Length + 1);
            Edges fewer = new(edges.NodeCount, edges.Sources[..keep], edges.Targets[..keep], edges.Values[..keep]);
            for (int k = 0; k <= 3; k++)
            {
                int[] finer = fewer.Closure(6).Partition(k, 1), coarser = edges.Closure(6).Partition(k, 1);
                for (int a = 0; a < edges.NodeCount; a++)
                    for (int b = a + 1; b < edges.NodeCount; b++)
                        if (finer[a] == finer[b])
                            Assert.Equal(coarser[a], coarser[b]);
            }
        }
    }

    [Fact]
    public void Sweep_ClustersAreConnectedAndNoOneEdgeCutRemains()
    {
        var rng = new Xoshiro256PlusPlus(20261106);
        for (int g = 0; g < SweepGraphs; g++)
        {
            Edges edges = RandomGraph(rng, 12);
            ThinCutClosure closure = edges.Closure(6);
            int m = 1 + rng.NextInt(4);

            foreach (int k in new[] { 1, 2, 3 })
            {
                int[] labels = closure.Partition(k, m);
                foreach (IGrouping<int, int> cluster in Enumerable.Range(0, edges.NodeCount).GroupBy(v => labels[v]))
                {
                    int[] members = cluster.ToArray();
                    var inside = new List<(int, int)>();
                    for (int e = 0; e < edges.Sources.Length; e++)
                        if (labels[edges.Sources[e]] == cluster.Key && labels[edges.Targets[e]] == cluster.Key)
                            inside.Add((edges.Sources[e], edges.Targets[e]));

                    var seen = new HashSet<int> { members[0] };
                    for (bool grew = true; grew;)
                    {
                        grew = false;
                        foreach ((int s, int t) in inside)
                            if (seen.Contains(s) != seen.Contains(t))
                                grew = seen.Add(s) | seen.Add(t);
                    }
                    Assert.Equal(members.Length, seen.Count);

                    if (k != 1) continue;

                    // A bridge is always a tree edge, so at one edge the rule misses no cut.
                    for (int mask = 1; mask < (1 << members.Length) - 1; mask++)
                    {
                        int chosen = System.Numerics.BitOperations.PopCount((uint)mask);
                        if (chosen < m || members.Length - chosen < m) continue;
                        var side = new HashSet<int>(members.Where((_, i) => (mask >> i & 1) == 1));
                        Assert.True(inside.Count(p => side.Contains(p.Item1) != side.Contains(p.Item2)) > 1);
                    }
                }
            }
        }
    }

    [Fact]
    public void Sweep_ThePartitionIgnoresNodeAndEdgeOrder()
    {
        var rng = new Xoshiro256PlusPlus(20261107);
        for (int g = 0; g < SweepGraphs; g++)
        {
            Edges edges = RandomGraph(rng, 30);
            int n = edges.NodeCount;
            int[] rename = Enumerable.Range(0, n).ToArray();
            int[] shuffle = Enumerable.Range(0, edges.Sources.Length).ToArray();
            for (int i = n - 1; i > 0; i--)
            {
                int j = rng.NextInt(i + 1);
                (rename[i], rename[j]) = (rename[j], rename[i]);
            }
            for (int i = shuffle.Length - 1; i > 0; i--)
            {
                int j = rng.NextInt(i + 1);
                (shuffle[i], shuffle[j]) = (shuffle[j], shuffle[i]);
            }

            Edges moved = new(n, shuffle.Select(e => rename[edges.Sources[e]]).ToArray(), shuffle.Select(e => rename[edges.Targets[e]]).ToArray(), shuffle.Select(e => edges.Values[e]).ToArray());
            int[] first = edges.Closure(6).Partition(2, 3), second = moved.Closure(6).Partition(2, 3);
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    Assert.Equal(first[a] == first[b], second[rename[a]] == second[rename[b]]);
        }
    }
}
