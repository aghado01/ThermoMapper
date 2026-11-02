using System;
using Clustering.Dendrograms;
using Maths.Rng;

namespace VizCore.Tests;

/// <summary>
/// A local search for the partition of least pair loss on a co-membership field: Binder's loss, for the probe in
/// <c>issues/spc/bayesian-clustering-design.md</c>, which asks how much of the closure gap a partition reaches. It
/// lives with the tests until that question says it should be anything more.
/// The pair loss is a constant minus <c>Σ_joined (⟨n_ij⟩ − θ)</c>, so this maximizes that sum, which is correlation
/// clustering. Seven starts, each polished by single-node moves (Gionis et al. 2007) and positive merges, and the best
/// score kept: the closure; agglomeration by the summed weight (Keuper et al. 2015) and by the mean; and four
/// sequential allocations in random order (Dahl et al. 2022). Finding is a lower bound on what a partition can reach,
/// so a shortfall only makes the probe conservative.
/// </summary>
internal static class PairLossSearch
{
    private const double Improvement = 1e-12;
    private const int PolishRounds = 200;

    internal readonly record struct Found(int[] Labels, double Score, string Start);

    /// <summary>The sum over joined pairs of ⟨n_ij⟩ − θ. Between two partitions, its rise is the pair loss's fall.</summary>
    internal static double Score(DenseCoMembership field, int[] labels, double theta)
    {
        double score = 0.0;
        for (int i = 0; i < labels.Length; i++)
            for (int j = i + 1; j < labels.Length; j++)
                if (labels[i] == labels[j]) score += field.At(i, j) - theta;
        return score;
    }

    /// <summary>The best partition found on <paramref name="field"/>, starting among others from <paramref name="closure"/>.</summary>
    internal static Found Best(DenseCoMembership field, double theta, int[] closure, int orders, int seed)
    {
        double[][] weights = Weights(field, theta);
        var rng = new Xoshiro256PlusPlus(seed);
        Found best = Polish(weights, (int[])closure.Clone(), "closure");
        Keep(ref best, Polish(weights, Linked(weights, byMean: false), "sum-linkage"));
        Keep(ref best, Polish(weights, Linked(weights, byMean: true), "mean-linkage"));
        for (int order = 0; order < orders; order++)
            Keep(ref best, Polish(weights, Sequential(weights, rng), "sequential"));
        return best;
    }

    /// <summary>One start polished alone, with no other start to cover for it.</summary>
    internal static Found Polished(DenseCoMembership field, double theta, int[] start)
        => Polish(Weights(field, theta), (int[])start.Clone(), "given");

    private static double[][] Weights(DenseCoMembership field, double theta)
    {
        int n = field.NodeCount;
        var weights = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double[] row = weights[i] = new double[n];
            for (int j = 0; j < n; j++)
                if (j != i) row[j] = field.At(i, j) - theta;
        }
        return weights;
    }

    private static void Keep(ref Found best, Found other)
    {
        if (other.Score > best.Score + Improvement) best = other;
    }

    private static Found Polish(double[][] weights, int[] labels, string start)
    {
        for (int round = 0; round < PolishRounds; round++)
        {
            bool moved = MoveNodes(weights, labels);
            bool merged = Agglomerate(weights, labels, byMean: false);
            if (!moved && !merged) break;
        }
        Compact(labels);

        double score = 0.0;
        for (int i = 0; i < labels.Length; i++)
        {
            double[] row = weights[i];
            for (int j = i + 1; j < labels.Length; j++)
                if (labels[i] == labels[j]) score += row[j];
        }
        return new Found(labels, score, start);
    }

    private static int[] Linked(double[][] weights, bool byMean)
    {
        var labels = new int[weights.Length];
        for (int i = 0; i < labels.Length; i++) labels[i] = i;
        Agglomerate(weights, labels, byMean);
        return labels;
    }

    /// <summary>Allocates the nodes in random order, each to the cluster it gains most from, or alone if none gains.</summary>
    private static int[] Sequential(double[][] weights, Xoshiro256PlusPlus rng)
    {
        int n = weights.Length;
        var order = new int[n];
        for (int i = 0; i < n; i++) order[i] = i;
        for (int i = n - 1; i > 0; i--)
        {
            int j = (int)(rng.NextDouble() * (i + 1));
            (order[i], order[j]) = (order[j], order[i]);
        }

        var labels = new int[n];
        var gain = new double[n];
        int clusters = 0;
        for (int placed = 0; placed < n; placed++)
        {
            int node = order[placed];
            double[] row = weights[node];
            Array.Clear(gain, 0, clusters);
            for (int p = 0; p < placed; p++)
                gain[labels[order[p]]] += row[order[p]];

            int target = clusters;
            double top = Improvement;
            for (int c = 0; c < clusters; c++)
                if (gain[c] > top)
                {
                    top = gain[c];
                    target = c;
                }
            labels[node] = target;
            if (target == clusters) clusters++;
        }
        return labels;
    }

    /// <summary>One or more sweeps of single-node moves, to whichever cluster gains or to a new one alone.</summary>
    private static bool MoveNodes(double[][] weights, int[] labels)
    {
        int n = labels.Length;
        Compact(labels);
        var sizes = new int[n];
        foreach (int label in labels) sizes[label]++;
        var free = new int[n];
        int freeCount = 0;
        for (int l = n - 1; l >= 0; l--)
            if (sizes[l] == 0) free[freeCount++] = l;

        var gain = new double[n];
        bool any = false, moved = true;
        for (int sweep = 0; moved && sweep < PolishRounds; sweep++)
        {
            moved = false;
            for (int node = 0; node < n; node++)
            {
                double[] row = weights[node];
                Array.Clear(gain);
                for (int j = 0; j < n; j++)
                    gain[labels[j]] += row[j];

                int here = labels[node];
                double top = gain[here] + Improvement;
                int target = here;
                for (int c = 0; c < n; c++)
                    if (sizes[c] > 0 && gain[c] > top)
                    {
                        top = gain[c];
                        target = c;
                    }
                if (top < 0.0) target = free[--freeCount];   // every cluster costs, and alone costs nothing
                if (target == here) continue;

                sizes[here]--;
                if (sizes[here] == 0) free[freeCount++] = here;
                sizes[target]++;
                labels[node] = target;
                moved = any = true;
            }
        }
        return any;
    }

    /// <summary>
    /// Merges clusters while some pair's total weight is positive, the pair with the largest total or the largest mean.
    /// Each cluster caches its best partner, rescanned only when a merge lowers that partner's key.
    /// </summary>
    private static bool Agglomerate(double[][] weights, int[] labels, bool byMean)
    {
        int k = Compact(labels);
        var between = new double[k][];
        for (int a = 0; a < k; a++) between[a] = new double[k];
        var sizes = new long[k];
        for (int i = 0; i < labels.Length; i++)
        {
            int a = labels[i];
            sizes[a]++;
            double[] row = weights[i];
            for (int j = i + 1; j < labels.Length; j++)
            {
                int b = labels[j];
                if (a == b) continue;
                between[a][b] += row[j];
                between[b][a] += row[j];
            }
        }

        var alive = new bool[k];
        var parent = new int[k];
        var partner = new int[k];
        var key = new double[k];
        for (int a = 0; a < k; a++)
        {
            alive[a] = true;
            parent[a] = a;
        }

        double Key(int a, int b) => byMean ? between[a][b] / (sizes[a] * sizes[b]) : between[a][b];

        void Refresh(int a)
        {
            partner[a] = -1;
            key[a] = double.NegativeInfinity;
            for (int b = 0; b < k; b++)
            {
                if (b == a || !alive[b]) continue;
                double value = Key(a, b);
                if (value > key[a])
                {
                    key[a] = value;
                    partner[a] = b;
                }
            }
        }

        for (int a = 0; a < k; a++) Refresh(a);

        bool merged = false;
        while (true)
        {
            int keep = -1;
            double top = Improvement;
            for (int c = 0; c < k; c++)
                if (alive[c] && key[c] > top)
                {
                    top = key[c];
                    keep = c;
                }
            if (keep < 0) break;

            int gone = partner[keep];
            alive[gone] = false;
            parent[gone] = keep;
            sizes[keep] += sizes[gone];
            merged = true;
            for (int c = 0; c < k; c++)
            {
                if (!alive[c] || c == keep) continue;
                between[keep][c] += between[gone][c];
                between[c][keep] = between[keep][c];
            }
            for (int c = 0; c < k; c++)
            {
                if (!alive[c] || c == keep) continue;
                double value = Key(c, keep);
                bool touched = partner[c] == keep || partner[c] == gone;
                if (value > key[c] || (touched && value == key[c]))
                {
                    key[c] = value;
                    partner[c] = keep;
                }
                else if (touched && key[c] > Improvement)
                {
                    Refresh(c);    // at or under zero the old key still bounds every partner, so no merge is missed
                }
            }
            Refresh(keep);
        }

        for (int i = 0; i < labels.Length; i++)
        {
            int root = labels[i];
            while (parent[root] != root) root = parent[root];
            labels[i] = root;
        }
        return merged;
    }

    /// <summary>Renames the labels to [0, k) by first appearance, and returns k.</summary>
    private static int Compact(int[] labels)
    {
        var names = new int[labels.Length];
        Array.Fill(names, -1);
        int k = 0;
        for (int i = 0; i < labels.Length; i++)
        {
            ref int name = ref names[labels[i]];
            if (name < 0) name = k++;
            labels[i] = name;
        }
        return k;
    }
}
