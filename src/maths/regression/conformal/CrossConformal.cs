using System;
using System.Collections.Generic;

namespace Maths.Regression.Conformal;

/// <summary>Row i contains the fit excluding training unit i: its prediction at that unit and at every test input.</summary>
public sealed record LeaveOneOutReadings(IReadOnlyList<double> Response, IReadOnlyList<double> HeldOut, double[,] Predictions);

/// <summary>A jackknife+ interval, including its attainable finite-sample lower bound (<c>jackknifePlus_coverage</c>).
/// The bound assumes exchangeable units and leave-out fits symmetric under relabeling, including their randomness.</summary>
public sealed record JackknifePlusInterval(int Case, int Units, int Level, double Lower, double Upper)
{
    public bool WholeSpace => Level == 0;
    public bool Empty => Lower > Upper;
    public double Width => Empty ? 0.0 : Upper - Lower;
    public bool Contains(double y) => y >= Lower && y <= Upper;
    public double Coverage => CrossConformal.Coverage(Units, Level);
}

public sealed record JackknifePlusIntervals(IReadOnlyList<JackknifePlusInterval> Intervals, double Alpha);

/// <summary>Reduce leave-one-out readings to jackknife+ intervals or n-fold pairwise sets.
/// The guarantees (<c>crossFit_coverage</c>, <c>jackknifePlus_coverage</c>) require exchangeability of all n + 1 units.
/// Each fit must read its training set symmetrically, possibly under a relabeling of its randomness. An order-reading sampler
/// can use one uniform random order shared by all leave-out fits, or randomness attached to the units.</summary>
public static class CrossConformal
{
    internal static double Coverage(int units, int level)
        => level == 0 ? 1.0 : Math.Max(0.0, 1.0 - (2.0 * level - 1) / units);

    /// <summary>Closed ends from the k-th smallest μ₋ᵢ(x) − |yᵢ − μ₋ᵢ(xᵢ)| and k-th largest μ₋ᵢ(x) + |yᵢ − μ₋ᵢ(xᵢ)|,
    /// k = ⌊α(n + 1)⌋ (<c>JackRefuses</c>). At k = 0 return the whole line. Rows must come from actual leave-one-out fits.</summary>
    public static JackknifePlusIntervals JackknifePlus(LeaveOneOutReadings readings, double alpha = 0.1)
    {
        ArgumentNullException.ThrowIfNull(readings);
        ArgumentNullException.ThrowIfNull(readings.Response);
        ArgumentNullException.ThrowIfNull(readings.HeldOut);
        ArgumentNullException.ThrowIfNull(readings.Predictions);
        SplitConformal.CheckAlpha(alpha);
        int n = readings.Response.Count, m = checked(n + 1), k = (int)Math.Floor(alpha * m);
        if (readings.HeldOut.Count != n || readings.Predictions.GetLength(0) != n)
            throw new ArgumentException("One held-out response and one prediction row per leave-one-out fit are required.", nameof(readings));
        var residuals = new double[n];
        for (int i = 0; i < n; i++)
        {
            if (!double.IsFinite(readings.Response[i]) || !double.IsFinite(readings.HeldOut[i]))
                throw new ArgumentException("Responses and predictions must be finite.", nameof(readings));
            residuals[i] = Math.Abs(readings.Response[i] - readings.HeldOut[i]);
        }
        foreach (double v in readings.Predictions)
            if (!double.IsFinite(v)) throw new ArgumentException("Predictions must be finite.", nameof(readings));
        int cases = readings.Predictions.GetLength(1);
        var intervals = new JackknifePlusInterval[cases];
        var lower = new double[n];
        var upper = new double[n];
        for (int t = 0; t < cases; t++)
        {
            if (k == 0)
            {
                intervals[t] = new(t, m, k, double.NegativeInfinity, double.PositiveInfinity);
                continue;
            }
            for (int i = 0; i < n; i++)
            {
                lower[i] = readings.Predictions[i, t] - residuals[i];
                upper[i] = readings.Predictions[i, t] + residuals[i];
            }
            Array.Sort(lower);
            Array.Sort(upper);
            intervals[t] = new(t, m, k, lower[k - 1], upper[n - k]);
        }
        return new(intervals, alpha);
    }

    /// <summary>The n-fold pairwise rule (<c>AcceptsPairwise</c>): candidate score [i,g] and held-out score [i] use the same
    /// fit excluding unit i. Higher scores conform better. A strict candidate &lt; held-out comparison counts as strange, and the
    /// candidate is accepted when the count is below n + 1 − k. Ties count for neither side.</summary>
    public static ConformalPrediction Pairwise(IReadOnlyList<double> heldOut, double[,] candidateScores,
                                               IReadOnlyList<double[]> grids, double alpha = 0.1)
    {
        ArgumentNullException.ThrowIfNull(heldOut);
        ArgumentNullException.ThrowIfNull(candidateScores);
        ArgumentNullException.ThrowIfNull(grids);
        SplitConformal.CheckAlpha(alpha);
        if (grids.Count == 0) throw new ArgumentException("At least one test case is required.", nameof(grids));
        if (candidateScores.GetLength(0) != heldOut.Count)
            throw new ArgumentException("One candidate-score row per leave-one-out fit is required.", nameof(candidateScores));
        foreach (double v in heldOut)
            if (double.IsNaN(v)) throw new ArgumentException("A held-out score is NaN.", nameof(heldOut));
        foreach (double v in candidateScores)
            if (double.IsNaN(v)) throw new ArgumentException("A candidate score is NaN.", nameof(candidateScores));
        int columns = 0;
        var copies = new double[grids.Count][];
        for (int t = 0; t < grids.Count; t++)
        {
            double[] grid = grids[t] ?? throw new ArgumentException("A grid is null.", nameof(grids));
            if (grid.Length == 0) throw new ArgumentException("A grid is empty.", nameof(grids));
            for (int g = 0; g < grid.Length; g++)
                if (!double.IsFinite(grid[g]) || (g > 0 && !(grid[g] > grid[g - 1])))
                    throw new ArgumentException("Grid values must be finite and strictly ascending.", nameof(grids));
            copies[t] = (double[])grid.Clone();
            columns = checked(columns + grid.Length);
        }
        if (columns != candidateScores.GetLength(1))
            throw new ArgumentException("The grids must cover all candidate-score columns.", nameof(candidateScores));
        int m = checked(heldOut.Count + 1), k = (int)Math.Floor(alpha * m);
        var sets = new ConformalSet[grids.Count];
        int column = 0;
        for (int t = 0; t < sets.Length; t++)
        {
            double[] grid = copies[t];
            int first = column;
            column += grid.Length;
            if (k == 0)
            {
                sets[t] = new(t, 0, m, k, grid, null, null, null, null, null, Form: ConformalForm.CrossFitted);
                continue;
            }
            var below = new int[grid.Length];
            var tied = new int[grid.Length];
            var p = new double[grid.Length];
            for (int i = 0; i < heldOut.Count; i++)
                for (int g = 0; g < grid.Length; g++)
                {
                    double candidate = candidateScores[i, first + g];
                    if (heldOut[i] < candidate) below[g]++;
                    else if (heldOut[i] == candidate) tied[g]++;
                }
            for (int g = 0; g < grid.Length; g++) p[g] = (below[g] + tied[g] + 1) / (double)m;
            sets[t] = new(t, 0, m, k, grid, p, below, tied, null, null, Form: ConformalForm.CrossFitted);
        }
        return new(sets, alpha, ConformalForm.CrossFitted);
    }

    /// <summary>Pairwise sets scored by each leave-one-out fit's predictive density. Each held-out seam holds one response on
    /// one singleton grid. Each candidate seam holds the same test grids. Draws align within each fit and can differ across fits.</summary>
    public static ConformalPrediction Density(IReadOnlyList<CandidateLogDensity> heldOut,
                                              IReadOnlyList<CandidateLogDensity> candidates, double alpha = 0.1)
    {
        ArgumentNullException.ThrowIfNull(heldOut);
        ArgumentNullException.ThrowIfNull(candidates);
        SplitConformal.CheckAlpha(alpha);
        if (heldOut.Count == 0 || candidates.Count != heldOut.Count)
            throw new ArgumentException("A nonempty set of paired leave-one-out seams is required.", nameof(heldOut));
        return Density(heldOut.Count, i => (heldOut[i], candidates[i]), alpha);
    }

    /// <summary>Reduce one leave-one-out fit at a time. The callback's matrices can be released after each call, retaining
    /// only the fit-by-candidate scores.</summary>
    public static ConformalPrediction Density(int fits, Func<int, (CandidateLogDensity HeldOut, CandidateLogDensity Candidates)> read,
                                              double alpha = 0.1)
    {
        ArgumentNullException.ThrowIfNull(read);
        ArgumentOutOfRangeException.ThrowIfLessThan(fits, 1);
        SplitConformal.CheckAlpha(alpha);
        var own = new double[fits];
        double[,]? scores = null;
        IReadOnlyList<double[]>? grids = null;
        double[]? buffer = null;
        for (int i = 0; i < fits; i++)
        {
            (CandidateLogDensity at, CandidateLogDensity test) = read(i);
            SplitConformal.CheckCalibration(at);
            if (at.Cases != 1) throw new ArgumentException("Each fit must score its one held-out unit.", nameof(read));
            new SplitDraws(at).Check(test);
            if (at.StratumOf[0] != 0 || Array.Exists(test.StratumOf, h => h != 0))
                throw new ArgumentException("Supply one exchangeable group in stratum zero.", nameof(read));
            grids ??= test.Grids;
            scores ??= new double[fits, test.Values.GetLength(1)];
            buffer ??= new double[test.Values.GetLength(1)];
            if (test.Cases != grids.Count) throw new ArgumentException("Every fit must score the same test cases.", nameof(read));
            for (int t = 0; t < grids.Count; t++)
                if (!ConformalBayes.SameSequence(grids[t], test.Grids[t]))
                    throw new ArgumentException("Every fit must score the same test grids.", nameof(read));
            (int[] rows, double[]? weights) = PredictiveChecks.Rows(at.Weights, at.Draws);
            Span<double> value = own.AsSpan(i, 1);
            SplitConformal.LogPredictive(at.Values, rows, weights, 0, 1, value);
            SplitConformal.LogPredictive(test.Values, rows, weights, 0, buffer.Length, buffer);
            for (int g = 0; g < buffer.Length; g++) scores[i, g] = buffer[g];
        }
        return Pairwise(own, scores!, grids!, alpha);
    }
}
