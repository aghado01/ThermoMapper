#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Selector over the <see cref="DiagramMetrics"/> distance backends. Exact Hungarian matching is
/// O((n+m)³) in bar count — the SPRED block-scale wall (ISOLET brief, "H0 matching-cost gate");
/// the other two are the screening-scale alternatives.
/// </summary>
public enum DiagramDistanceKind
{
    /// <summary>Exact W_p: balanced diagonal-augmented Hungarian assignment, O((n+m)³).</summary>
    Wasserstein,

    /// <summary>Sliced W_p: deterministic evenly-spaced slices, O(L·(n+m)·log(n+m)). A distinct
    /// (strongly equivalent) metric, not a numerical approximation of W_p under L∞.</summary>
    SlicedWasserstein,

    /// <summary>Entropic W_p: log-domain Sinkhorn on the exact cost geometry,
    /// O(iters·(n+m)²); ε → 0 recovers the exact value.</summary>
    SinkhornWasserstein,
}

/// <summary>
/// Diagram distances on <see cref="Barcode"/> / <see cref="Bar"/>.
/// Unrelated to <c>Maths.Distance.Geodesic.Wasserstein1</c> (1-D PMF transport).
/// </summary>
public static class DiagramMetrics
{
    internal enum EssentialKind { InfiniteOnMismatch, FinitePenalty }

    /// <summary>
    /// Policy for essential (infinite-death) bars, matched separately from finite bars: with
    /// death = ∞ on both sides, the L∞ ground metric between two essentials reduces to |Δbirth|,
    /// so min-count essentials pair by birth (see <see cref="MatchEssentialBirths"/>) and the
    /// policy governs only the count surplus. <see cref="InfiniteOnMismatch"/> (the default via
    /// <c>default(EssentialPolicy)</c>) returns +∞ on any count mismatch; <see cref="FinitePenalty"/>
    /// charges perBar^p per surplus bar — a surplus essential has infinite persistence, so this
    /// deliberately caps a genuinely infinite transport term.
    /// </summary>
    public readonly record struct EssentialPolicy
    {
        internal EssentialKind Kind { get; init; }
        internal double PerBar { get; init; }

        public static EssentialPolicy InfiniteOnMismatch => default;

        public static EssentialPolicy FinitePenalty(double perBar)
        {
            if (!double.IsFinite(perBar) || perBar < 0.0)
                throw new ArgumentOutOfRangeException(nameof(perBar),
                    "Essential per-bar penalty is a distance scale: finite and >= 0 (zero deliberately " +
                    "disables the surplus charge; a negative value would reward essential-count mismatch).");
            return new() { Kind = EssentialKind.FinitePenalty, PerBar = perBar };
        }
    }

    /// <summary>
    /// Wasserstein distance W_p between two barcodes in the given homological dimension.
    /// L∞ ground metric; balanced (n+m) assignment (Kerber–Morozov–Nigmetov). Essential bars
    /// match by birth under the <see cref="EssentialPolicy"/>.
    /// </summary>
    public static double Wasserstein(
        Barcode a,
        Barcode b,
        int dimension,
        double p = 2.0,
        EssentialPolicy essential = default)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        ValidateP(p);

        var essA = CollectEssentialBirths(a, dimension);
        var essB = CollectEssentialBirths(b, dimension);

        if (essential.Kind == EssentialKind.InfiniteOnMismatch && essA.Count != essB.Count)
            return double.PositiveInfinity;

        var finiteA = CollectFinite(a, dimension);
        var finiteB = CollectFinite(b, dimension);

        double[,] cost = BuildCost(finiteA, finiteB, p, out double big);
        int s = finiteA.Count + finiteB.Count;
        double assignmentSum = s == 0 ? 0.0 : MinAssignment(cost, s);

        assignmentSum += MatchEssentialBirths(essA, essB, p);

        if (essential.Kind == EssentialKind.FinitePenalty)
        {
            int surplus = Math.Abs(essA.Count - essB.Count);
            if (surplus > 0)
                assignmentSum += surplus * Math.Pow(essential.PerBar, p);
        }

        return Math.Pow(assignmentSum, 1.0 / p);
    }

    /// <summary>
    /// Sliced Wasserstein distance between two barcodes, after Carrière–Cuturi–Oudot (2017),
    /// generalized from their p = 1 kernel to order p: each side is augmented with the diagonal
    /// projections of the other's points (mirroring the balanced assignment of
    /// <see cref="Wasserstein"/>), both are projected onto <paramref name="directions"/> evenly
    /// spaced lines, and each 1-D transport is the sorted matching. Returns
    /// ((1/L)·Σ_l W_p^p(θ_l))^(1/p). O(L·(n+m)·log(n+m)) — the cheap screening metric; a distinct
    /// (strongly equivalent) metric with L2 slice geometry, not a numerical approximation of the
    /// L∞-ground-metric W_p. Deterministic: fixed slices, no RNG. Essential bars contribute the
    /// same slice-independent birth-matched term and surplus policy as <see cref="Wasserstein"/> —
    /// an essential has no finite death to project, and its exact transport is already 1-D.
    /// </summary>
    public static double SlicedWasserstein(
        Barcode a,
        Barcode b,
        int dimension,
        double p = 2.0,
        EssentialPolicy essential = default,
        int directions = 50)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        ValidateP(p);
        if (directions < 1)
            throw new ArgumentOutOfRangeException(nameof(directions), "At least one slice direction is required.");

        var essA = CollectEssentialBirths(a, dimension);
        var essB = CollectEssentialBirths(b, dimension);
        if (essential.Kind == EssentialKind.InfiniteOnMismatch && essA.Count != essB.Count)
            return double.PositiveInfinity;

        var finiteA = CollectFinite(a, dimension);
        var finiteB = CollectFinite(b, dimension);
        int s = finiteA.Count + finiteB.Count;

        double meanPowerSum = 0.0;   // (1/L)·Σ_l W_p^p(θ_l)
        if (s > 0)
        {
            var projA = new double[s];
            var projB = new double[s];

            for (int l = 0; l < directions; l++)
            {
                double theta = Math.PI * (l + 0.5) / directions - Math.PI / 2.0;
                double cos = Math.Cos(theta), sin = Math.Sin(theta);
                double diagScale = 0.5 * (cos + sin);   // projection of ((b+d)/2, (b+d)/2)

                int idx = 0;
                for (int i = 0; i < finiteA.Count; i++) projA[idx++] = finiteA[i].Birth * cos + finiteA[i].Death * sin;
                for (int j = 0; j < finiteB.Count; j++) projA[idx++] = (finiteB[j].Birth + finiteB[j].Death) * diagScale;

                idx = 0;
                for (int j = 0; j < finiteB.Count; j++) projB[idx++] = finiteB[j].Birth * cos + finiteB[j].Death * sin;
                for (int i = 0; i < finiteA.Count; i++) projB[idx++] = (finiteA[i].Birth + finiteA[i].Death) * diagScale;

                Array.Sort(projA);
                Array.Sort(projB);

                double slice = 0.0;
                for (int i = 0; i < s; i++)
                    slice += Math.Pow(Math.Abs(projA[i] - projB[i]), p);
                meanPowerSum += slice;
            }
            meanPowerSum /= directions;
        }

        meanPowerSum += MatchEssentialBirths(essA, essB, p);   // slice-independent, enters once

        if (essential.Kind == EssentialKind.FinitePenalty)
        {
            int surplus = Math.Abs(essA.Count - essB.Count);
            if (surplus > 0)
                meanPowerSum += surplus * Math.Pow(essential.PerBar, p);
        }

        return Math.Pow(meanPowerSum, 1.0 / p);
    }

    /// <summary>
    /// Entropic (Sinkhorn) Wasserstein distance between two barcodes: log-domain Sinkhorn on the
    /// same diagonal-augmented balanced cost matrix as the exact <see cref="Wasserstein"/>, so the
    /// value converges to the exact one as <paramref name="epsilon"/> → 0 (the assignment LP is
    /// integral, and the entropic plan is LP-feasible, so the result is a near-upper bound).
    /// O(iters·(n+m)²) time and O((n+m)²) memory vs the exact O((n+m)³).
    /// <paramref name="epsilon"/> is dimensionless — the cost matrix is normalized by its largest
    /// finite entry before smoothing, so the same ε transfers across data scales.
    /// <para>Entropic bias: the self-distance is positive (smoothing smears mass onto
    /// near-diagonal escape cells whose cost is comparable to ε) and shrinks as ε → 0 — hold ε
    /// fixed when comparing values, as the SPRED objective does.</para>
    /// <para>Essential bars enter as the exact birth-matched term of <see cref="Wasserstein"/>,
    /// never smoothed — Betti-scale counts leave nothing worth relaxing.</para>
    /// </summary>
    public static double SinkhornWasserstein(
        Barcode a,
        Barcode b,
        int dimension,
        double p = 2.0,
        EssentialPolicy essential = default,
        double epsilon = 0.01,
        int maxIters = 500)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        ValidateP(p);
        if (!double.IsFinite(epsilon) || epsilon <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(epsilon), "Entropic regularization must be finite and > 0.");
        if (maxIters < 1)
            throw new ArgumentOutOfRangeException(nameof(maxIters), "At least one Sinkhorn iteration is required.");

        var essA = CollectEssentialBirths(a, dimension);
        var essB = CollectEssentialBirths(b, dimension);
        if (essential.Kind == EssentialKind.InfiniteOnMismatch && essA.Count != essB.Count)
            return double.PositiveInfinity;

        var finiteA = CollectFinite(a, dimension);
        var finiteB = CollectFinite(b, dimension);
        int s = finiteA.Count + finiteB.Count;

        double transportSum = 0.0;
        if (s > 0)
        {
            double[,] cost = BuildCost(finiteA, finiteB, p, out double big);
            transportSum = SinkhornAssignment(cost, s, big, epsilon, maxIters);
        }

        transportSum += MatchEssentialBirths(essA, essB, p);

        if (essential.Kind == EssentialKind.FinitePenalty)
        {
            int surplus = Math.Abs(essA.Count - essB.Count);
            if (surplus > 0)
                transportSum += surplus * Math.Pow(essential.PerBar, p);
        }

        return Math.Pow(transportSum, 1.0 / p);
    }

    /// <summary>Bottleneck distance — deferred; gate uses <see cref="Wasserstein"/>.</summary>
    public static double Bottleneck(
        Barcode a,
        Barcode b,
        int dimension,
        EssentialPolicy essential = default) =>
        throw new NotImplementedException("Bottleneck (Hopcroft–Karp threshold matching) is P1.");

    static void ValidateP(double p)
    {
        if (double.IsNaN(p) || p < 1.0)
            throw new ArgumentOutOfRangeException(nameof(p), "Wasserstein requires p >= 1.");
        if (double.IsPositiveInfinity(p))
            throw new ArgumentOutOfRangeException(nameof(p), "Use Bottleneck for p = infinity.");
    }

    static List<double> CollectEssentialBirths(Barcode barcode, int dimension)
    {
        var births = new List<double>();
        foreach (Bar bar in barcode.Bars)
            if (bar.Dimension == dimension && bar.IsInfinite)
                births.Add(bar.Birth);
        return births;
    }

    static List<Bar> CollectFinite(Barcode barcode, int dimension)
    {
        var list = new List<Bar>();
        foreach (Bar bar in barcode.Bars)
            if (bar.Dimension == dimension && !bar.IsInfinite)
                list.Add(bar);
        return list;
    }

    /// <summary>
    /// Matched-essential term: min(|A|,|B|) essential bars paired by birth — their only finite
    /// coordinate (death = ∞ on both sides collapses the L∞ ground metric to |Δbirth|). Solved
    /// with the same balanced <see cref="MinAssignment"/> as the finite matching, on a birth-only
    /// matrix whose surplus rows/columns stay at zero so the solver also chooses which surplus
    /// bars go unmatched — sorted-order pairing is exact for equal counts (1-D transport with
    /// convex ground distance is monotone) but cannot make that choice. Essential counts are
    /// Betti-number scale, so O(k³) is immaterial. The surplus itself is charged by the
    /// <see cref="EssentialPolicy"/> at the call site.
    /// </summary>
    static double MatchEssentialBirths(List<double> birthsA, List<double> birthsB, double p)
    {
        int n = birthsA.Count, m = birthsB.Count;
        if (n == 0 || m == 0)
            return 0.0;

        int k = Math.Max(n, m);
        var cost = new double[k, k];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                cost[i, j] = Math.Pow(Math.Abs(birthsA[i] - birthsB[j]), p);
        return MinAssignment(cost, k);
    }

    static double Dinf(in Bar p, in Bar q) =>
        Math.Max(Math.Abs(p.Birth - q.Birth), Math.Abs(p.Death - q.Death));

    static double Diag(in Bar p) => 0.5 * (p.Death - p.Birth);

    static double[,] BuildCost(IReadOnlyList<Bar> a, IReadOnlyList<Bar> b, double p, out double big)
    {
        int n = a.Count, m = b.Count, s = n + m;
        var c = new double[s, s];
        double maxReal = 0.0;

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                maxReal = Math.Max(maxReal, Dinf(a[i], b[j]));
        for (int i = 0; i < n; i++)
            maxReal = Math.Max(maxReal, Diag(a[i]));
        for (int j = 0; j < m; j++)
            maxReal = Math.Max(maxReal, Diag(b[j]));

        big = Math.Pow(maxReal + 1.0, p) * (s + 1) + 1.0;

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
                c[i, j] = Math.Pow(Dinf(a[i], b[j]), p);
            for (int k = 0; k < n; k++)
                c[i, m + k] = k == i ? Math.Pow(Diag(a[i]), p) : big;
        }

        for (int j = 0; j < m; j++)
        {
            for (int jj = 0; jj < m; jj++)
                c[n + j, jj] = jj == j ? Math.Pow(Diag(b[j]), p) : big;
            for (int k = 0; k < n; k++)
                c[n + j, m + k] = 0.0;
        }

        return c;
    }

    /// <summary>
    /// Log-domain Sinkhorn on a square cost matrix with unit row/column marginals — the entropic
    /// relaxation of <see cref="MinAssignment"/>. Entries are normalized by the largest admissible
    /// entry so <paramref name="epsilon"/> is dimensionless; forbidden cells (≥ <paramref name="big"/>)
    /// are masked as log-kernel −∞ inside both LSE sweeps and omitted from the transport sum, so
    /// the plan is confined to the declared diagonal-augmented support at every ε — a finite
    /// sentinel only underflows for small ε and would otherwise receive real mass. Returns the
    /// transport objective Σ π_ij·C_ij on the original scale.
    /// </summary>
    static double SinkhornAssignment(double[,] cost, int n, double big, double epsilon, int maxIters)
    {
        const double tol = 1e-9;

        double cMax = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (cost[i, j] < big && cost[i, j] > cMax)
                    cMax = cost[i, j];
        if (cMax <= 0.0)
            return 0.0;   // every admissible cell is free — the optimum is zero

        double eps = epsilon;                     // on the normalized scale C/cMax
        var f = new double[n];
        var g = new double[n];
        var scratch = new double[n];

        for (int iter = 0; iter < maxIters; iter++)
        {
            // f_i ← −ε·LSE_j((g_j − C'_ij)/ε);  g_j ← −ε·LSE_i((f_i − C'_ij)/ε)   (log a_i = log b_j = 0)
            // Forbidden cells enter as −∞ — outside the support at any ε, not merely expensive.
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    scratch[j] = cost[i, j] >= big
                        ? double.NegativeInfinity
                        : (g[j] - cost[i, j] / cMax) / eps;
                f[i] = -eps * LogSumExp(scratch);
            }

            double violation = 0.0;
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                    scratch[i] = cost[i, j] >= big
                        ? double.NegativeInfinity
                        : (f[i] - cost[i, j] / cMax) / eps;
                double lse = LogSumExp(scratch);
                violation = Math.Max(violation, Math.Abs(Math.Exp(g[j] / eps + lse) - 1.0));
                g[j] = -eps * lse;
            }

            if (violation < tol)
                break;
        }

        double transport = 0.0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (cost[i, j] >= big)                    // masked support — no mass at any ε
                    continue;
                double logPi = (f[i] + g[j] - cost[i, j] / cMax) / eps;
                if (logPi > -700.0)                       // exp underflow guard
                    transport += Math.Exp(logPi) * cost[i, j];
            }
        }
        return transport;
    }

    static double LogSumExp(ReadOnlySpan<double> values)
    {
        double max = double.NegativeInfinity;
        for (int i = 0; i < values.Length; i++)
            if (values[i] > max) max = values[i];
        if (double.IsNegativeInfinity(max))
            return double.NegativeInfinity;

        double sum = 0.0;
        for (int i = 0; i < values.Length; i++)
            sum += Math.Exp(values[i] - max);
        return max + Math.Log(sum);
    }

    /// <summary>Potential-based Hungarian (Kuhn–Munkres), O(n³). Square n×n, min total cost.</summary>
    static double MinAssignment(double[,] a, int n)
    {
        const double inf = double.PositiveInfinity;
        var u = new double[n + 1];
        var v = new double[n + 1];
        var p = new int[n + 1];
        var way = new int[n + 1];

        for (int i = 1; i <= n; i++)
        {
            p[0] = i;
            int j0 = 0;
            var minv = new double[n + 1];
            var used = new bool[n + 1];
            Array.Fill(minv, inf);

            do
            {
                used[j0] = true;
                int i0 = p[j0];
                int j1 = -1;
                double delta = inf;

                for (int j = 1; j <= n; j++)
                {
                    if (used[j]) continue;
                    double cur = a[i0 - 1, j - 1] - u[i0] - v[j];
                    if (cur < minv[j])
                    {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta)
                    {
                        delta = minv[j];
                        j1 = j;
                    }
                }

                for (int j = 0; j <= n; j++)
                {
                    if (used[j])
                    {
                        u[p[j]] += delta;
                        v[j] -= delta;
                    }
                    else
                    {
                        minv[j] -= delta;
                    }
                }

                j0 = j1;
            } while (p[j0] != 0);

            do
            {
                int j1 = way[j0];
                p[j0] = p[j1];
                j0 = j1;
            } while (j0 != 0);
        }

        double res = 0.0;
        for (int j = 1; j <= n; j++)
            res += a[p[j] - 1, j - 1];
        return res;
    }
}
