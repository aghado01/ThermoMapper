using System;
using System.Collections.Generic;
using System.Linq;

namespace Maths.Samplers.FlatHistogram;

/// <summary>
/// A WHAM estimate of the density of states from a tempering pilot: each pooled draw's log mass under ĝ and the rungs'
/// log normalizers f_i = ln Z(b_i), up to one constant, at the fixed point of the multiple-histogram equations.
/// </summary>
/// <param name="Levels">The pooled draws' distinct ℓ₁, ascending.</param>
/// <param name="LogMasses">ĝ's log mass at each: <c>ln m − ln Σ_k ñ_k e^{b_k·ℓ − f_k}</c>, m the draws at that ℓ₁ each
/// counted 1/g of its rung, and ñ_k = n_k/g_k rung k's effective draws.</param>
/// <param name="LogNormalizers">f_i per rung, the hottest at 0.</param>
/// <param name="Iterations">Newton's iterations.</param>
/// <param name="Residual">The fixed point's largest residual, max_i |Σ_j e^{b_i·ℓ_j − f_i}·ĝ_j − 1|.</param>
/// <param name="EffectiveDraws">ñ_i per rung, the weights its histogram entered with.</param>
public sealed record WhamEstimate(double[] Levels, double[] LogMasses, double[] LogNormalizers, int Iterations, double Residual,
                                  double[] EffectiveDraws);

/// <summary>
/// The warm start: multiple-histogram reweighting (Ferrenberg and Swendsen 1989; WHAM, in its binless form) of a
/// tempering pilot's rung draws into ln g, then SAMC's own limit on that density of states
/// (<see cref="SamcLearning.FixedPoint"/>), so the learning starts where it would converge were the pilot exact and
/// corrects the pilot rather than discovering g from flat. It also reads the run's window, b_lo and K from the pilot.
/// <para>The rungs' laws are <c>prior·e^{b_i·ℓ}/Z(b_i)</c>, so the draws pooled over every rung estimate g as the measure
/// <c>ĝ_j ∝ (1/g_{i(j)}) / Σ_k ñ_k e^{b_k·ℓ_j − f_k}</c> with <c>e^{f_i} = Σ_j ĝ_j e^{b_i·ℓ_j}</c>, the free energies'
/// fixed point: Ferrenberg and Swendsen's equations (3) and (4) in their list form. Each rung enters by its effective
/// draws ñ_i = n_i/g_i, g_i = 1 + 2τ_i its statistical inefficiency (<see cref="TemperingPilot.Inefficiencies"/>), their
/// error-minimizing combination; any positive weights would be consistent, these lower the error where rungs mix at
/// different rates. It is the stationary point of the convex
/// <c>F(f) = Σ_j (1/g_{i(j)}) ln Σ_k ñ_k e^{b_k·ℓ_j − f_k} + Σ_i ñ_i f_i</c> (the MBAR objective), which Newton's method
/// solves from the path-sampling estimate of f, the derivative-guided iteration they recommend. Where the pilot drew
/// nothing, as between two phases it does not cross, ĝ has no mass and the start's weights are interpolated; where it drew
/// too little for SAMC's limit on ĝ to exist, the start ties the thin nodes to their neighbours
/// (<see cref="WarmStart(TemperingPilot, MulticanonicalWeights, double, out int[])"/>).</para>
/// </summary>
public static class Wham
{
    /// <summary>The density of states the pilot's rungs estimate.</summary>
    public static WhamEstimate Estimate(TemperingPilot pilot)
    {
        ArgumentNullException.ThrowIfNull(pilot);
        int rungs = pilot.Rungs;
        int draws = checked((int)pilot.Draws);
        var pooled = new double[draws];
        var drawWeight = new double[draws];   // 1/g of the draw's rung
        var counts = new double[rungs];       // ñ_i = n_i/g_i
        var b = new double[rungs];
        int at = 0;
        for (int i = 0; i < rungs; i++)
        {
            ReadOnlySpan<double> l = pilot.Levels(i);
            double g = pilot.Inefficiencies[i];
            l.CopyTo(pooled.AsSpan(at));
            drawWeight.AsSpan(at, l.Length).Fill(1.0 / g);
            at += l.Length;
            counts[i] = l.Length / g;
            b[i] = pilot.InverseTemperatures[i];
        }
        var logCounts = Array.ConvertAll(counts, Math.Log);

        // Draws at one ℓ₁ share their mass, so a lattice energy's draws collapse to its levels, each weighing its draws.
        Array.Sort(pooled, drawWeight);
        var levelList = new List<double>();
        var multiplicity = new List<double>();
        for (int j = 0; j < draws; j++)
        {
            if (levelList.Count > 0 && levelList[^1] == pooled[j]) multiplicity[^1] += drawWeight[j];
            else
            {
                levelList.Add(pooled[j]);
                multiplicity.Add(drawWeight[j]);
            }
        }
        double[] levels = [.. levelList], weight = [.. multiplicity];
        int total = levels.Length;

        // The start: path sampling's ln Z(b_i) − ln Z(b_hot), the trapezoid of the rungs' mean ℓ₁ in b.
        var f = new double[rungs];
        var order = new int[rungs];
        for (int i = 0; i < rungs; i++) order[i] = i;
        Array.Sort(order, (p, q) => b[p].CompareTo(b[q]));
        for (int r = 1; r < rungs; r++)
        {
            int lo = order[r - 1], hi = order[r];
            f[hi] = f[lo] + 0.5 * (b[hi] - b[lo]) * (Mean(pilot.Levels(lo)) + Mean(pilot.Levels(hi)));
        }

        // Newton on F with the hottest rung pinned at f = 0.
        int pin = pilot.HotRung;
        var free = new List<int>(rungs - 1);
        for (int i = 0; i < rungs; i++) if (i != pin) free.Add(i);
        int n = free.Count;
        var occupancy = new double[rungs];       // Σ_j P_ji
        var cross = new double[rungs * rungs];   // Σ_j P_ji P_jk
        var logDenominators = new double[total];
        double value = Objective(f, b, logCounts, counts, levels, weight, occupancy, cross, logDenominators);
        int iterations = 0;
        for (; iterations < 100 && n > 0; iterations++)
        {
            var gradient = new double[n];
            double norm = 0.0;
            for (int a = 0; a < n; a++)
            {
                int i = free[a];
                gradient[a] = counts[i] - occupancy[i];
                norm = Math.Max(norm, Math.Abs(gradient[a]) / counts[i]);
            }
            if (norm < 1e-10) break;
            var hessian = new double[n * n];
            for (int a = 0; a < n; a++)
                for (int c = 0; c < n; c++)
                    hessian[a * n + c] = (a == c ? occupancy[free[a]] : 0.0) - cross[free[a] * rungs + free[c]];
            double[] direction = SamcLearning.SolveSpd(hessian, gradient, n);
            double slope = 0.0;
            for (int a = 0; a < n; a++) slope -= gradient[a] * direction[a];
            double t = 1.0, trialValue = double.NaN;
            var trial = (double[])f.Clone();
            bool descended = false;
            for (int halving = 0; halving < 40; halving++, t *= 0.5)
            {
                for (int a = 0; a < n; a++) trial[free[a]] = f[free[a]] - t * direction[a];
                trialValue = Objective(trial, b, logCounts, counts, levels, weight, occupancy, cross, logDenominators);
                if (trialValue <= value + 1e-4 * t * slope) { descended = true; break; }
            }
            if (!descended)
            {
                Objective(f, b, logCounts, counts, levels, weight, occupancy, cross, logDenominators);
                break;
            }
            Array.Copy(trial, f, rungs);
            bool stalled = value - trialValue <= 1e-14 * Math.Max(1.0, Math.Abs(value));
            value = trialValue;
            if (stalled) break;
        }

        double residual = 0.0;
        for (int i = 0; i < rungs; i++) residual = Math.Max(residual, Math.Abs(occupancy[i] / counts[i] - 1.0));
        var logMasses = new double[total];
        for (int j = 0; j < total; j++) logMasses[j] = Math.Log(weight[j]) - logDenominators[j];
        return new WhamEstimate(levels, logMasses, f, iterations, residual, counts);
    }

    /// <summary>Fewer of the pilot's draws than this under a node's hat (<see cref="NodeDraws"/>), and the warm start ties
    /// the node where SAMC's fixed point refuses ĝ as it stands: 10, the least of the thresholds measured (1, 2, 5 and 10)
    /// at which every refused pilot read, on fixture A, its gap pilots, the double well and Friedman, was accepted (the warm
    /// start's chip, 2026-10-31). A pilot too small for two nodes to hold ten draws is refused.</summary>
    public const double ThinDraws = 10.0;

    /// <summary>The warm-start weights over <paramref name="shape"/>'s window: SAMC's limit on WHAM's ĝ, or where the
    /// fixed point refuses ĝ, its limit with the nodes under <see cref="ThinDraws"/> pilot draws tied.</summary>
    public static MulticanonicalWeights WarmStart(TemperingPilot pilot, MulticanonicalWeights shape) => WarmStart(pilot, shape, ThinDraws, out _);

    /// <summary>
    /// The warm-start weights over <paramref name="shape"/>'s window: SAMC's limit on WHAM's ĝ
    /// (<see cref="SamcLearning.FixedPoint(MulticanonicalWeights, ReadOnlySpan{double}, ReadOnlySpan{double}, out bool[])"/>),
    /// returned as it stands wherever the fixed point accepts ĝ. ĝ is a discrete measure, its atoms the pilot's draws, and
    /// where a node's hat sees it only through a few atoms far from the node, beside an empty segment or the window's edge,
    /// no finite limit exists and the fixed point refuses it. Then the nodes holding fewer than
    /// <paramref name="thinDraws"/> of the pilot's draws are tied to their neighbours, their ln w on the line through the
    /// nearest held nodes and past the outermost at the extensions' slopes
    /// (<see cref="SamcLearning.FixedPoint(MulticanonicalWeights, ReadOnlySpan{double}, ReadOnlySpan{double}, ReadOnlySpan{bool})"/>),
    /// and the limit on that coarser interpolant is the start, which SAMC's learning then corrects at every node. Where even
    /// that is refused, or fewer than two nodes hold the draws, the start is refused (<see cref="InvalidOperationException"/>),
    /// saying which, and where flat occupancy leaves the hull of ĝ's atoms' hats (<see cref="SamcLearning.Infeasibility"/>).
    /// </summary>
    /// <param name="pilot">The tempering pilot.</param>
    /// <param name="shape">The window, slope and node count; its values are not read.</param>
    /// <param name="thinDraws">The draws under a node's hat below which a refused ĝ's node is tied; positive.</param>
    /// <param name="tied">The nodes tied, ascending; none where the fixed point accepted ĝ as it stands.</param>
    public static MulticanonicalWeights WarmStart(TemperingPilot pilot, MulticanonicalWeights shape, double thinDraws, out int[] tied)
    {
        ArgumentNullException.ThrowIfNull(shape);
        if (!(thinDraws > 0.0) || double.IsInfinity(thinDraws))
            throw new ArgumentOutOfRangeException(nameof(thinDraws), "The thin nodes' draws must be positive and finite.");
        WhamEstimate estimate = Estimate(pilot);
        MulticanonicalWeights? start = SamcLearning.Solve(shape, estimate.Levels, estimate.LogMasses, out _, out double departure);
        if (start is not null)
        {
            tied = [];
            return start;
        }
        double[] draws = NodeDraws(pilot, shape);
        var thin = new bool[shape.Nodes];
        for (int k = 0; k < thin.Length; k++) thin[k] = draws[k] < thinDraws;
        tied = [.. Enumerable.Range(0, thin.Length).Where(k => thin[k])];
        string refusal = SamcLearning.Refusal(shape, estimate.Levels, estimate.LogMasses, departure);
        if (thin.Length - tied.Length < 2)
            throw new InvalidOperationException($"{refusal}; and fewer than two nodes hold {thinDraws:G} of the pilot's draws to tie the rest to.");
        start = SamcLearning.SolveTied(shape, estimate.Levels, estimate.LogMasses, thin, out double tiedDeparture);
        return start ?? throw new InvalidOperationException(
            $"{refusal}; and with the nodes under {thinDraws:G} pilot draws tied [{string.Join(", ", tied)}], a merged node's occupancy departs " +
            $"from flat by {tiedDeparture:G3}.");
    }

    /// <summary>The pilot's draws under each node's hat over <paramref name="shape"/>'s window, every rung's: Σ φ_k(ℓ₁) over
    /// the draws inside the window. Where a node holds none, ĝ has no mass there.</summary>
    public static double[] NodeDraws(TemperingPilot pilot, MulticanonicalWeights shape)
    {
        ArgumentNullException.ThrowIfNull(pilot);
        ArgumentNullException.ThrowIfNull(shape);
        var counts = new double[shape.Nodes];
        for (int i = 0; i < pilot.Rungs; i++)
            foreach (double level in pilot.Levels(i))
            {
                if (!shape.TryLocate(level, out int segment, out double fraction)) continue;
                counts[segment] += 1.0 - fraction;
                counts[segment + 1] += fraction;
            }
        return counts;
    }

    /// <summary>
    /// The window from the pilot: from ℓ₁'s <paramref name="quantile"/>-quantile at the hottest rung less
    /// <paramref name="margin"/> of its sd, to its (1 − q)-quantile at b = 1 plus <paramref name="margin"/> of its sd, each
    /// held inside the pilot's extreme draws, past which the extensions carry the hot end's and the posterior's own laws.
    /// With <paramref name="topMargin"/>, the top is the (1 − q)-quantile at b = 1 plus that many of its sd, not held inside
    /// the pilot's highest draw, so the window can reach above what the pilot drew.
    /// </summary>
    public static FlatHistogramWindow Window(TemperingPilot pilot, double quantile, double margin, double? topMargin = null)
    {
        ArgumentNullException.ThrowIfNull(pilot);
        if (!(quantile > 0.0 && quantile < 0.5)) throw new ArgumentOutOfRangeException(nameof(quantile));
        if (!(margin >= 0.0)) throw new ArgumentOutOfRangeException(nameof(margin));
        if (topMargin is double t && !(t >= 0.0 && double.IsFinite(t))) throw new ArgumentOutOfRangeException(nameof(topMargin));
        double min = double.PositiveInfinity, max = double.NegativeInfinity;
        for (int i = 0; i < pilot.Rungs; i++)
            foreach (double l in pilot.Levels(i))
            {
                min = Math.Min(min, l);
                max = Math.Max(max, l);
            }
        ReadOnlySpan<double> hot = pilot.Levels(pilot.HotRung), cold = pilot.Levels(pilot.ColdRung);
        double low = Math.Max(Quantile(hot, quantile) - margin * StandardDeviation(hot), min);
        double high = topMargin is double above
            ? Quantile(cold, 1.0 - quantile) + above * StandardDeviation(cold)
            : Math.Min(Quantile(cold, 1.0 - quantile) + margin * StandardDeviation(cold), max);
        if (!(high > low)) throw new ArgumentException($"The pilot gives an empty window [{low}, {high}].", nameof(pilot));
        return new FlatHistogramWindow(low, high);
    }

    /// <summary>K from the pilot: ⌈width / sd₁⌉ + 1 nodes, sd₁ the b = 1 rung's ℓ₁ sd, so no segment is wider than
    /// the posterior's spread; at least 2.</summary>
    public static int Nodes(TemperingPilot pilot, FlatHistogramWindow window)
    {
        ArgumentNullException.ThrowIfNull(pilot);
        ArgumentNullException.ThrowIfNull(window);
        double sd = StandardDeviation(pilot.Levels(pilot.ColdRung));
        if (!(sd > 0.0)) throw new ArgumentException("The pilot's b = 1 rung has no spread in ℓ₁ to size the nodes by.", nameof(pilot));
        return Math.Max(2, checked((int)Math.Ceiling((window.High - window.Low) / sd)) + 1);
    }

    // F(f), with each rung's Σ_j m_j P_ji, the pairs' Σ_j m_j P_ji P_jk, and each level's ln Σ_k ñ_k e^{b_k ℓ_j − f_k}; a
    // level stands for its draws, m_j their summed weights 1/g.
    private static double Objective(double[] f, double[] b, double[] logCounts, double[] counts, double[] levels, double[] weight,
                                    double[] occupancy, double[] cross, double[] logDenominators)
    {
        int rungs = f.Length;
        Array.Clear(occupancy);
        Array.Clear(cross);
        Span<double> terms = stackalloc double[rungs];
        double value = 0.0;
        for (int j = 0; j < levels.Length; j++)
        {
            double max = double.NegativeInfinity;
            for (int k = 0; k < rungs; k++)
            {
                terms[k] = logCounts[k] + b[k] * levels[j] - f[k];
                max = Math.Max(max, terms[k]);
            }
            double sum = 0.0;
            for (int k = 0; k < rungs; k++) sum += terms[k] = Math.Exp(terms[k] - max);
            double logD = max + Math.Log(sum);
            double m = weight[j];
            logDenominators[j] = logD;
            value += m * logD;
            for (int i = 0; i < rungs; i++)
            {
                double p = terms[i] / sum;
                occupancy[i] += m * p;
                for (int k = 0; k < rungs; k++) cross[i * rungs + k] += m * p * terms[k] / sum;
            }
        }
        for (int i = 0; i < rungs; i++) value += counts[i] * f[i];
        return value;
    }

    private static double Mean(ReadOnlySpan<double> values)
    {
        double s = 0.0;
        foreach (double v in values) s += v;
        return s / values.Length;
    }

    private static double StandardDeviation(ReadOnlySpan<double> values)
    {
        double mean = Mean(values), ss = 0.0;
        foreach (double v in values) ss += (v - mean) * (v - mean);
        return Math.Sqrt(ss / (values.Length - 1));
    }

    // The empirical quantile, linear between order statistics (type 7).
    private static double Quantile(ReadOnlySpan<double> values, double p)
    {
        double[] sorted = values.ToArray();
        Array.Sort(sorted);
        double h = (sorted.Length - 1) * p;
        int lo = (int)Math.Floor(h);
        int hi = Math.Min(lo + 1, sorted.Length - 1);
        return sorted[lo] + (h - lo) * (sorted[hi] - sorted[lo]);
    }
}
