using System;
using System.Collections.Generic;

namespace Maths.Regression.Conformal;

/// <summary>
/// Highest-density sets with their threshold chosen by risk control (<see cref="RiskControlledSets.Density"/>): each case's grid values
/// whose density score is at most <paramref name="Threshold"/>, the chosen threshold's index among the calibration scores
/// (<paramref name="Choice"/>), and the confidence β and level α it was chosen at.
/// </summary>
public sealed record HighDensitySets(IReadOnlyList<double[]> Grids, IReadOnlyList<bool[]> Accepts, double Threshold, RiskThreshold Choice,
                                     double Alpha, double Confidence)
{
    /// <summary>The accepted share of case <paramref name="t"/>'s grid cells, read as <see cref="ConformalSet.Width"/> reads them.</summary>
    public double Width(int t)
    {
        double[] grid = Grids[t];
        if (grid.Length == 1) return Accepts[t][0] ? 1.0 : 0.0;
        double width = 0.0;
        for (int g = 0; g < grid.Length; g++)
            if (Accepts[t][g]) width += ConformalBayes.Cell(grid, g);
        return width;
    }

    /// <summary>Whether case <paramref name="t"/>'s set holds the grid value nearest <paramref name="y"/>.</summary>
    public bool Holds(int t, double y)
    {
        double[] grid = Grids[t];
        int best = 0;
        for (int g = 1; g < grid.Length; g++) if (Math.Abs(grid[g] - y) < Math.Abs(grid[best] - y)) best = g;
        return Accepts[t][best];
    }
}

/// <summary>
/// Bayesian conformal prediction as a decision risk problem (Wu, Lohmanova, Kaski and Caprio, arXiv 2602.03331; the conformal layer's
/// slice 17): the split form's density score s = −ln p̂(y | x), from a fit to training data apart from the calibration units, and the set
/// {y : s(x, y) ≤ λ*}, with λ* the least threshold whose Gaffke bound on the calibration miscoverage, b*_β of ℓᵢ(λ) = 1{sᵢ &gt; λ}, is at most
/// α (<see cref="RiskControl.HighProbability"/>). The losses fall as λ rises, so the set's miscoverage exceeds α with chance at most 1 − β
/// over calibration sets of independent units (<c>fail_le_fixed</c>, through Vlassis and Thomas's calibration): coverage 1 − α with
/// confidence β, a training-conditional guarantee where CB-1's is marginal. Their Bayesian quadrature estimates the expected set size to
/// choose among feasible thresholds; the least feasible one is the smallest set, so the guarantee needs no quadrature.
/// </summary>
public static class RiskControlledSets
{
    /// <summary>
    /// The sets for <paramref name="candidates"/>' cases from <paramref name="calibration"/>, each calibration unit on a grid of its observed
    /// response, both seams the same fit's draws: thresholds at every calibration score, ascending, and +∞ beyond them.
    /// </summary>
    /// <param name="calibration">The calibration units' seam, each on a grid of its observed response.</param>
    /// <param name="candidates">The test cases' seam on their grids, the same draws.</param>
    /// <param name="alpha">α, the miscoverage held.</param>
    /// <param name="confidence">β, the chance over calibration sets that the miscoverage is at most α.</param>
    public static HighDensitySets Density(CandidateLogDensity calibration, CandidateLogDensity candidates, double alpha = 0.1, double confidence = 0.9)
    {
        SplitConformal.CheckCalibration(calibration);
        new SplitDraws(calibration).Check(candidates);
        SplitConformal.CheckAlpha(alpha);
        if (!(confidence > 0.0 && confidence < 1.0)) throw new ArgumentOutOfRangeException(nameof(confidence), "β lies in (0, 1).");
        (int[] rows, double[]? logWeights) = PredictiveChecks.Rows(calibration.Weights, calibration.Draws);
        int n = calibration.Cases;
        var logDensity = new double[n];
        SplitConformal.LogPredictive(calibration.Values, rows, logWeights, 0, n, logDensity);
        var scores = new double[n];
        for (int i = 0; i < n; i++) scores[i] = -logDensity[i];
        double[] thresholds = [.. scores, double.PositiveInfinity];
        Array.Sort(thresholds);
        var losses = new double[n, thresholds.Length];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < thresholds.Length; j++) losses[i, j] = scores[i] > thresholds[j] ? 1.0 : 0.0;
        RiskThreshold choice = RiskControl.HighProbability(losses, 1.0, alpha, confidence);
        double threshold = thresholds[choice.Index];

        int columns = candidates.Values.GetLength(1);
        var candidateDensity = new double[columns];
        SplitConformal.LogPredictive(candidates.Values, rows, logWeights, 0, columns, candidateDensity);
        var accepts = new bool[candidates.Cases][];
        for (int t = 0; t < candidates.Cases; t++)
        {
            int offset = candidates.Offset(t), g = candidates.Grids[t].Length;
            accepts[t] = new bool[g];
            for (int j = 0; j < g; j++) accepts[t][j] = -candidateDensity[offset + j] <= threshold;
        }
        return new HighDensitySets(candidates.Grids, accepts, threshold, choice, alpha, confidence);
    }
}
