using System;
using System.Collections.Generic;
using Maths.Rng;

namespace Maths.Regression.Conformal;

/// <summary>How posterior conformal prediction weighs the units for a test case (<see cref="PosteriorConformal"/>).</summary>
public enum PosteriorWeighting
{
    /// <summary>
    /// Their randomized weights (Zhang and Candès, Theorem 1; <c>coverage_pcp</c>): L ~ Multinomial(m, π(x_test)) drawn for the case, and
    /// wᵢ ∝ Π_k π_k(xᵢ)^{L_k}, the test's own included: the posterior chance that unit i is the test unit given the draw. The rule covers at
    /// least 1 − α given L, so marginally too (<c>coverage_pcp_marginal</c>).
    /// </summary>
    Randomized,

    /// <summary>
    /// Their KL weights (Theorem 3; <c>pcp_cluster_gap</c>): wᵢ ∝ exp(−m·KL(π(x_test) ‖ π(xᵢ))), fixed given the inputs, the test's the largest.
    /// Under a mixture model of the scores the miscoverage given the inputs is at most α + Σᵢ qᵢ Σ_k |π_k(xᵢ) − π_k(x_test)| at the true
    /// memberships, half their constant (the Lean, through <c>refused_le_swap_tv</c>); marginal coverage is not exact, since the weights are not
    /// randomized (their Appendix on randomization).
    /// </summary>
    KullbackLeibler,
}

/// <summary>
/// The weights posterior conformal prediction reads for each test case (<see cref="PosteriorConformal.Weights"/>): <paramref name="Weights"/>,
/// per case, for the weighted rule's forms; each case's draw L (<paramref name="Draws"/>, null for the KL weights); the effective sample size
/// of each case's normalized weights over the units and the test, 1/Σ qᵢ²; each case's self-weight, the test's own normalized weight, which
/// is the chance its set is the whole space where it exceeds α; and, for the KL weights, each case's cluster gap
/// Σᵢ qᵢ Σ_k |π_k(xᵢ) − π_k(x_test)| read at the memberships given, the bound's excess over α where they are the true ones.
/// </summary>
public sealed record PosteriorWeights(ShiftWeights Weights, IReadOnlyList<int[]>? Draws, IReadOnlyList<double> EffectiveSampleSizes,
                                      IReadOnlyList<double> SelfWeights, IReadOnlyList<double>? ClusterGaps);

/// <summary>
/// Posterior conformal prediction (Zhang and Candès, arXiv 2409.19712; the conformal layer's slice 15): coverage near nominal per cluster of
/// a mixture model of the scores, beyond named strata. Each unit's scores are taken to follow Σ_k π_k(x) f_k, with memberships π(x) a
/// function of the inputs; a test case weighs the calibration (or training) units by how alike their memberships are to its own, and the
/// weighted rule (<see cref="ShiftWeights"/>, slice 14) ranks it under those weights, case by case (<see cref="ShiftWeights.PerCase"/>).
/// The memberships must be fixed given the inputs, or learned by an algorithm symmetric in the units and the test (their Proposition 2):
/// memberships fitted on a separate training split are both. The precision m trades how sharply the weights single out alike units against
/// their effective sample size.
/// </summary>
public static class PosteriorConformal
{
    /// <summary>
    /// The weights for each test case from the memberships at the units (<paramref name="unitMemberships"/>, a row each, in the units' order)
    /// and at the cases (<paramref name="caseMemberships"/>), each row nonnegative and finite with a positive sum, normalized here.
    /// </summary>
    /// <param name="unitMemberships">π at each unit, a row each in the units' order and a column per cluster.</param>
    /// <param name="caseMemberships">π at each test case, the same clusters.</param>
    /// <param name="precision">m, positive.</param>
    /// <param name="seed">The draws' seed under the randomized weights: case t draws from a stream of its own.</param>
    /// <param name="grids">Each case's grid, for the density forms, each case's weight spread over its columns; null for the closed forms.</param>
    /// <param name="weighting">The randomized weights, the default, or the KL weights.</param>
    public static PosteriorWeights Weights(double[,] unitMemberships, double[,] caseMemberships, int precision, int seed = 0,
                                           IReadOnlyList<double[]>? grids = null, PosteriorWeighting weighting = PosteriorWeighting.Randomized)
    {
        if (!Enum.IsDefined(weighting)) throw new ArgumentOutOfRangeException(nameof(weighting));
        ArgumentOutOfRangeException.ThrowIfLessThan(precision, 1);
        (double[][] units, double[][] tests) = Checked(unitMemberships, caseMemberships, grids);
        Func<int, int[]>? drawFor = weighting == PosteriorWeighting.Randomized ? t => Draw(tests[t], precision, Stream(seed, t)) : null;
        return Build(units, tests, precision, grids, drawFor);
    }

    /// <summary>The randomized weights at given draws, <paramref name="draws"/>[t] case t's L: the rule given L, which the exactness
    /// gate enumerates.</summary>
    internal static PosteriorWeights GivenDraws(double[,] unitMemberships, double[,] caseMemberships, IReadOnlyList<int[]> draws,
                                                IReadOnlyList<double[]>? grids = null)
    {
        (double[][] units, double[][] tests) = Checked(unitMemberships, caseMemberships, grids);
        if (draws.Count != tests.Length) throw new ArgumentException($"One draw per case is required: {draws.Count} draws for {tests.Length} cases.", nameof(draws));
        return Build(units, tests, 0, grids, t => draws[t]);
    }

    private static (double[][] Units, double[][] Tests) Checked(double[,] unitMemberships, double[,] caseMemberships, IReadOnlyList<double[]>? grids)
    {
        ArgumentNullException.ThrowIfNull(unitMemberships);
        ArgumentNullException.ThrowIfNull(caseMemberships);
        int cases = caseMemberships.GetLength(0), k = unitMemberships.GetLength(1);
        if (caseMemberships.GetLength(1) != k) throw new ArgumentException("The cases' memberships must cover the units' clusters.", nameof(caseMemberships));
        if (k == 0) throw new ArgumentException("At least one cluster is required.", nameof(unitMemberships));
        if (grids is not null && grids.Count != cases) throw new ArgumentException($"One grid per case is required: {grids.Count} grids for {cases} cases.", nameof(grids));
        return (Normalized(unitMemberships, nameof(unitMemberships)), Normalized(caseMemberships, nameof(caseMemberships)));
    }

    // Each case's weights: the multinomial kernel at its draw where drawFor is given, the KL kernel at precision m otherwise.
    private static PosteriorWeights Build(double[][] units, double[][] tests, int precision, IReadOnlyList<double[]>? grids, Func<int, int[]>? drawFor)
    {
        int n = units.Length, cases = tests.Length, k = tests.Length > 0 ? tests[0].Length : 0;
        var byCase = new double[cases][];
        var own = new double[cases];
        var ess = new double[cases];
        var self = new double[cases];
        int[][]? draws = drawFor is null ? null : new int[cases][];
        double[]? gaps = drawFor is null ? new double[cases] : null;
        var log = new double[n];
        for (int t = 0; t < cases; t++)
        {
            double ownLog;
            if (draws is not null)
            {
                int[] draw = drawFor!(t);
                if (draw.Length != k) throw new ArgumentException($"Case {t}'s draw counts {draw.Length} clusters for {k}.");
                draws[t] = draw;
                for (int i = 0; i < n; i++) log[i] = LogKernel(units[i], draw);
                ownLog = LogKernel(tests[t], draw);
                if (double.IsNegativeInfinity(ownLog)) throw new ArgumentException($"Case {t}'s draw counts a cluster its memberships give no chance.");
            }
            else
            {
                for (int i = 0; i < n; i++) log[i] = -precision * Kl(tests[t], units[i]);
                ownLog = 0.0;
            }
            double top = ownLog;
            foreach (double v in log) top = Math.Max(top, v);
            var w = new double[n];
            double total = Math.Exp(ownLog - top), squares = total * total;
            for (int i = 0; i < n; i++)
            {
                w[i] = Math.Exp(log[i] - top);
                total += w[i];
                squares += w[i] * w[i];
            }
            byCase[t] = w;
            own[t] = Math.Exp(ownLog - top);
            ess[t] = total * total / squares;
            self[t] = own[t] / total;
            if (gaps is not null)
            {
                double gap = 0.0;
                for (int i = 0; i < n; i++)
                {
                    double l1 = 0.0;
                    for (int c = 0; c < k; c++) l1 += Math.Abs(units[i][c] - tests[t][c]);
                    gap += w[i] / total * l1;
                }
                gaps[t] = gap;
            }
        }
        ShiftWeights weights = grids is null ? ShiftWeights.PerCase(byCase, own) : ShiftWeights.ByCase(byCase, own, grids);
        return new PosteriorWeights(weights, draws, ess, self, gaps);
    }

    /// <summary>ln Π_k π_k^{L_k}, the multinomial draw's likelihood under <paramref name="membership"/> up to its coefficient; −∞ where a drawn
    /// cluster has no membership.</summary>
    internal static double LogKernel(ReadOnlySpan<double> membership, int[] draw)
    {
        double log = 0.0;
        for (int c = 0; c < draw.Length; c++)
            if (draw[c] > 0) log += draw[c] * Math.Log(membership[c]);
        return log;
    }

    // KL(p ‖ q) = Σ p ln(p/q), +∞ where q is 0 and p is not.
    private static double Kl(double[] p, double[] q)
    {
        double kl = 0.0;
        for (int c = 0; c < p.Length; c++)
            if (p[c] > 0.0) kl += p[c] * (Math.Log(p[c]) - Math.Log(q[c]));
        return kl;
    }

    // m draws of a cluster from p, counted.
    internal static int[] Draw(double[] p, int m, Xoshiro256PlusPlus rng)
    {
        var counts = new int[p.Length];
        for (int d = 0; d < m; d++)
        {
            double u = rng.NextDouble(), cumulative = 0.0;
            int c = 0;
            for (; c < p.Length - 1; c++)
            {
                cumulative += p[c];
                if (u < cumulative) break;
            }
            while (p[c] == 0.0 && c > 0) c--;
            counts[c]++;
        }
        return counts;
    }

    private static double[][] Normalized(double[,] memberships, string name)
    {
        int rows = memberships.GetLength(0), k = memberships.GetLength(1);
        var result = new double[rows][];
        for (int i = 0; i < rows; i++)
        {
            double sum = 0.0;
            for (int c = 0; c < k; c++)
            {
                double v = memberships[i, c];
                if (!(v >= 0.0 && double.IsFinite(v))) throw new ArgumentOutOfRangeException(name, $"Row {i}'s membership {c} is {v}: each is finite and nonnegative.");
                sum += v;
            }
            if (!(sum > 0.0)) throw new ArgumentException($"Row {i}'s memberships sum to {sum}: a positive sum is required.", name);
            result[i] = new double[k];
            for (int c = 0; c < k; c++) result[i][c] = memberships[i, c] / sum;
        }
        return result;
    }

    // A stream per case from the seed: SplitMix64 over the seed and the case's index gives xoshiro's four words.
    internal static Xoshiro256PlusPlus Stream(int seed, int index)
    {
        ulong x = unchecked((ulong)(uint)seed * 0x9E3779B97F4A7C15UL ^ ((ulong)(uint)index + 1UL) * 0xD1B54A32D192ED03UL);
        ulong Next()
        {
            x = unchecked(x + 0x9E3779B97F4A7C15UL);
            ulong v = x;
            v = unchecked((v ^ (v >> 30)) * 0xBF58476D1CE4E5B9UL);
            v = unchecked((v ^ (v >> 27)) * 0x94D049BB133111EBUL);
            return v ^ (v >> 31);
        }
        ulong s0 = Next(), s1 = Next(), s2 = Next(), s3 = Next();
        return Xoshiro256PlusPlus.Restore(s0, s1, s2, s3 | 1UL);
    }
}
