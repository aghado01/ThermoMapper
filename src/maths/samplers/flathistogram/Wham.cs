using System;
using System.Collections.Generic;

namespace Maths.Samplers.FlatHistogram;

/// <summary>
/// A WHAM estimate of the density of states from a tempering pilot: each pooled draw's log mass under ĝ and the rungs'
/// log normalizers f_i = ln Z(b_i), up to one constant, at the fixed point of the multiple-histogram equations.
/// </summary>
/// <param name="Levels">The pooled draws' distinct ℓ₁, ascending.</param>
/// <param name="LogMasses">ĝ's log mass at each: <c>ln m − ln Σ_k n_k e^{b_k·ℓ − f_k}</c>, m the draws at that ℓ₁.</param>
/// <param name="LogNormalizers">f_i per rung, the hottest at 0.</param>
/// <param name="Iterations">Newton's iterations.</param>
/// <param name="Residual">The fixed point's largest residual, max_i |Σ_j e^{b_i·ℓ_j − f_i}·ĝ_j − 1|.</param>
public sealed record WhamEstimate(double[] Levels, double[] LogMasses, double[] LogNormalizers, int Iterations, double Residual);

/// <summary>
/// The warm start: multiple-histogram reweighting (Ferrenberg and Swendsen 1989; WHAM, in its binless form) of a
/// tempering pilot's rung draws into ln g, then SAMC's own limit on that density of states
/// (<see cref="SamcLearning.FixedPoint"/>), so the learning starts where it would converge were the pilot exact and
/// corrects the pilot rather than discovering g from flat. It also reads the run's window, b_lo and K from the pilot.
/// <para>The rungs' laws are <c>prior·e^{b_i·ℓ}/Z(b_i)</c>, so the draws pooled over every rung estimate g as the measure
/// <c>ĝ_j ∝ 1/Σ_k n_k e^{b_k·ℓ_j − f_k}</c> with <c>e^{f_i} = Σ_j ĝ_j e^{b_i·ℓ_j}</c>, the free energies' fixed point.
/// It is the stationary point of the convex <c>F(f) = Σ_j ln Σ_k n_k e^{b_k·ℓ_j − f_k} + Σ_i n_i f_i</c> (the MBAR
/// objective), which Newton's method solves from the path-sampling estimate of f. Raw draw counts weigh the rungs:
/// the estimate is consistent whatever their autocorrelation, which its error, not its value, depends on. Where the
/// pilot drew nothing, as between two phases it does not cross, ĝ has no mass and the start's weights are interpolated.</para>
/// </summary>
public static class Wham
{
    /// <summary>The density of states the pilot's rungs estimate.</summary>
    public static WhamEstimate Estimate(TemperingPilot pilot)
    {
        ArgumentNullException.ThrowIfNull(pilot);
        int rungs = pilot.Rungs;
        var pooled = new double[checked((int)pilot.Draws)];
        var counts = new double[rungs];
        var b = new double[rungs];
        int at = 0;
        for (int i = 0; i < rungs; i++)
        {
            ReadOnlySpan<double> l = pilot.Levels(i);
            l.CopyTo(pooled.AsSpan(at));
            at += l.Length;
            counts[i] = l.Length;
            b[i] = pilot.InverseTemperatures[i];
        }
        var logCounts = Array.ConvertAll(counts, Math.Log);

        // Draws at one ℓ₁ share their mass, so a lattice energy's draws collapse to its levels.
        Array.Sort(pooled);
        var levelList = new List<double>();
        var multiplicity = new List<double>();
        foreach (double l in pooled)
        {
            if (levelList.Count > 0 && levelList[^1] == l) multiplicity[^1]++;
            else
            {
                levelList.Add(l);
                multiplicity.Add(1.0);
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
        return new WhamEstimate(levels, logMasses, f, iterations, residual);
    }

    /// <summary>The warm-start weights over <paramref name="shape"/>'s window: SAMC's limit on WHAM's ĝ.</summary>
    public static MulticanonicalWeights WarmStart(TemperingPilot pilot, MulticanonicalWeights shape)
    {
        ArgumentNullException.ThrowIfNull(shape);
        WhamEstimate estimate = Estimate(pilot);
        return SamcLearning.FixedPoint(shape, estimate.Levels, estimate.LogMasses, out _);
    }

    /// <summary>
    /// The window from the pilot: from ℓ₁'s <paramref name="quantile"/>-quantile at the hottest rung less
    /// <paramref name="margin"/> of its sd, to its (1 − q)-quantile at b = 1 plus <paramref name="margin"/> of its sd, each
    /// held inside the pilot's extreme draws, past which the extensions carry the hot end's and the posterior's own laws.
    /// </summary>
    public static FlatHistogramWindow Window(TemperingPilot pilot, double quantile, double margin)
    {
        ArgumentNullException.ThrowIfNull(pilot);
        if (!(quantile > 0.0 && quantile < 0.5)) throw new ArgumentOutOfRangeException(nameof(quantile));
        if (!(margin >= 0.0)) throw new ArgumentOutOfRangeException(nameof(margin));
        double min = double.PositiveInfinity, max = double.NegativeInfinity;
        for (int i = 0; i < pilot.Rungs; i++)
            foreach (double l in pilot.Levels(i))
            {
                min = Math.Min(min, l);
                max = Math.Max(max, l);
            }
        ReadOnlySpan<double> hot = pilot.Levels(pilot.HotRung), cold = pilot.Levels(pilot.ColdRung);
        double low = Math.Max(Quantile(hot, quantile) - margin * StandardDeviation(hot), min);
        double high = Math.Min(Quantile(cold, 1.0 - quantile) + margin * StandardDeviation(cold), max);
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

    // F(f), with each rung's Σ_j P_ji, the pairs' Σ_j P_ji P_jk, and each level's ln Σ_k n_k e^{b_k ℓ_j − f_k}; a level
    // stands for its m_j draws.
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
