using System;
using System.Collections.Generic;
using Maths.Rng;

namespace Maths.Regression.Conformal;

/// <summary>Each held-out calibration group's one score for <see cref="GroupConformal.Density"/>, which decides what its sets cover.</summary>
public enum GroupScore
{
    /// <summary>The group's least-conforming case: a new group's cases are all held at once, with chance at least 1 − α.</summary>
    WorstCase,

    /// <summary>One case drawn per group, uniformly, on a seeded stream: a case drawn at random from a new group is held with chance at least
    /// 1 − α, and so is each of its cases where they are exchangeable within the group.</summary>
    OneCase,
}

/// <summary>
/// The split density form with whole groups as its units (the conformal layer's follow-up for new groups, agreed with the GBASS chip for
/// E21's slice 5): for a test case of a group no training unit holds, its effect integrated from its prior. A new group's cases are
/// correlated, so a calibration on cases does not cover them; held-out groups, each with one score, are exchangeable units if the groups
/// are drawn alike, and CB-1's rule on them covers: m = G + 1 units, G the calibration groups, k = ⌊α·m⌋, a candidate accepted where
/// (B + T + 1)/m &gt; k/m, B and T the groups scoring below and level with it (CB-1). The calibration groups must be held out of training and
/// scored as new groups are, their effects integrated; the fit is independent of them and of the test groups. Drawn alike includes their
/// sizes: the worst case of more cases is lower, so under <see cref="GroupScore.WorstCase"/> a new group larger than the calibration's
/// typical one is held less often (6 cases against groups of 3 to 8: 0.897 and 0.890 at 20,000 replications, 0.901 drawn alike, α 0.1).
/// <see cref="GroupScore.OneCase"/> reads one drawn case, whose law a group's size does not move where sizes are drawn apart from the
/// effects. Sets need G ≥ ⌈1/α⌉ − 1 to be
/// other than the whole space. The idea of one score per held-out group is Dunn, Wasserman and Ramdas's (two-layer conformal), as the
/// agreement names it; their paper is not on the shelf, and the construction here rests on CB-1 alone.
/// </summary>
public static class GroupConformal
{
    /// <summary>
    /// The sets for <paramref name="candidates"/>' cases, each of a new group, from <paramref name="calibration"/>'s cases, each on a grid of its
    /// observed response and of the group <paramref name="groupOf"/> names; both seams the same fit's draws, one stratum.
    /// </summary>
    /// <param name="calibration">The held-out groups' cases, each at its observed response, scored as new groups.</param>
    /// <param name="groupOf">Each calibration case's group.</param>
    /// <param name="candidates">The new groups' test cases on their grids.</param>
    /// <param name="alpha">α.</param>
    /// <param name="score">How each group is scored (<see cref="GroupScore"/>).</param>
    /// <param name="seed">The draws' seed under <see cref="GroupScore.OneCase"/>.</param>
    public static ConformalPrediction Density(CandidateLogDensity calibration, IReadOnlyList<int> groupOf, CandidateLogDensity candidates,
                                              double alpha = 0.1, GroupScore score = GroupScore.WorstCase, int seed = 0)
    {
        ArgumentNullException.ThrowIfNull(groupOf);
        SplitConformal.CheckAlpha(alpha);
        SplitConformal.CheckCalibration(calibration);
        // The units are groups, which span cells by design: the cell rule binds case units only.
        (new SplitDraws(calibration) with { StrataCells = null }).Check(candidates);
        if (!Enum.IsDefined(score)) throw new ArgumentOutOfRangeException(nameof(score));
        if (groupOf.Count != calibration.Cases)
            throw new ArgumentException($"One group per calibration case is required: {groupOf.Count} labels for {calibration.Cases} cases.", nameof(groupOf));
        if (Array.Exists(calibration.StratumOf, h => h != 0) || Array.Exists(candidates.StratumOf, h => h != 0))
            throw new ArgumentException("Groups are the units here, one exchangeable population of them: strata are not read.", nameof(candidates));
        if (calibration.NewGroup is { } calibrationNew && Array.Exists(calibrationNew, f => !f))
            throw new ArgumentException("A calibration case's group is known to the fit: the calibration groups are held out and scored as new groups.", nameof(calibration));
        if (candidates.NewGroup is { } candidateNew && Array.Exists(candidateNew, f => !f))
            throw new ArgumentException("A test case's group is known to the fit: such cases take the case-level split form within their cell.", nameof(candidates));

        (int[] rows, double[]? logWeights) = PredictiveChecks.Rows(calibration.Weights, calibration.Draws);
        var conformity = new double[calibration.Cases];
        SplitConformal.LogPredictive(calibration.Values, rows, logWeights, 0, calibration.Cases, conformity);
        double[] groupScores = Scores(conformity, groupOf, score, seed);
        var fit = new SplitCalibration(groupScores, null, alpha);
        int m = fit.Units(0), k = fit.Level(0);
        var sets = new ConformalSet[candidates.Cases];
        for (int t = 0; t < candidates.Cases; t++)
        {
            double[] grid = candidates.Grids[t];
            if (k == 0)
            {
                sets[t] = new ConformalSet(t, 0, m, k, grid, null, null, null, null, null, Form: ConformalForm.Split);
                continue;
            }
            int g = grid.Length;
            var own = new double[g];
            SplitConformal.LogPredictive(candidates.Values, rows, logWeights, candidates.Offset(t), g, own);
            var below = new int[g];
            var tied = new int[g];
            var p = new double[g];
            for (int j = 0; j < g; j++)
            {
                (below[j], tied[j]) = fit.Rank(0, own[j]);
                p[j] = (below[j] + tied[j] + 1) / (double)m;
            }
            sets[t] = new ConformalSet(t, 0, m, k, grid, p, below, tied, own, null, Form: ConformalForm.Split);
        }
        return new ConformalPrediction(sets, alpha, ConformalForm.Split);
    }

    // Each group's score, in the groups' ascending order of label: its least conformity, or one case's drawn on the group's own stream.
    internal static double[] Scores(double[] conformity, IReadOnlyList<int> groupOf, GroupScore score, int seed)
    {
        var members = new SortedDictionary<int, List<int>>();
        for (int i = 0; i < conformity.Length; i++)
        {
            if (!members.TryGetValue(groupOf[i], out List<int>? list)) members.Add(groupOf[i], list = []);
            list.Add(i);
        }
        var scores = new double[members.Count];
        int at = 0;
        foreach ((int label, List<int> cases) in members)
        {
            if (score == GroupScore.WorstCase)
            {
                double least = double.PositiveInfinity;
                foreach (int i in cases) least = Math.Min(least, conformity[i]);
                scores[at++] = least;
            }
            else
            {
                Xoshiro256PlusPlus rng = PosteriorConformal.Stream(seed, label);
                scores[at++] = conformity[cases[rng.NextInt(cases.Count)]];
            }
        }
        return scores;
    }
}
