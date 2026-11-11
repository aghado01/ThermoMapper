using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Maths.Regression.Conformal;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A random-effect fit's cells (E21's slice 5, the candidate seam as agreed with the conformal layer): a cell is one level of every grouping
/// jointly. Cases are exchangeable only within a cell, since a case's term reads its groups' effects, so every conformal stratum lies inside
/// one, in the full step and the split form alike (<see cref="ConformalBayes"/>). The training cases' cells are numbered by first appearance, 0 to
/// <paramref name="TrainingCells"/> − 1. A candidate in a cell no training case holds, a new group's among them, is numbered after them in
/// order of first appearance; candidates share a new cell when they name the same levels, and a candidate whose grouping column is left
/// out is a cell of its own. The cells are not a training density's default strata: a fit's strata are the caller's to give.
/// </summary>
/// <param name="Training">Each training case's cell.</param>
/// <param name="Candidates">Each candidate's cell; null when no cases were given.</param>
/// <param name="NewGroup">Each candidate: whether some grouping's level is one no training case took, or its column is left out, so its
/// term integrates that group's effect from its prior; null when no cases were given.</param>
/// <param name="TrainingCells">The cells the training cases hold.</param>
public sealed record BassGroupCells(int[] Training, int[]? Candidates, bool[]? NewGroup, int TrainingCells)
{
    /// <summary>The cells of <paramref name="model"/>'s training cases and, when given, of <paramref name="cases"/>, whose groups are read as
    /// prediction reads them (<see cref="PredictionCases"/>, <paramref name="unseenLevels"/>).</summary>
    internal static BassGroupCells Of(BassModel model, PredictionCases? cases, UnseenLevel unseenLevels)
    {
        RandomEffectsBlock block = model.RandomEffects ?? throw new NotSupportedException("A fit without random effects reads no group: its cases share one cell.");
        var numbers = new Dictionary<string, int>(StringComparer.Ordinal);
        var training = new int[block.Cases];
        var key = new StringBuilder();
        for (int i = 0; i < training.Length; i++)
        {
            key.Clear();
            for (int k = 0; k < block.Count; k++) AppendKey(key, "known", block.Layout(k).GroupOf(i).ToString(System.Globalization.CultureInfo.InvariantCulture));
            training[i] = Number(numbers, key.ToString());
        }
        int trainingCells = numbers.Count;
        if (cases is null) return new BassGroupCells(training, null, null, trainingCells);

        int n = cases.Cases;
        int[][] groups = [.. Enumerable.Range(0, block.Count).Select(k => cases.GroupsOf(block.Layout(k), block.Groupings[k], unseenLevels))];
        string[]?[] labels = [.. Enumerable.Range(0, block.Count).Select(k => cases.GroupLabelsOf(block.Layout(k), block.Groupings[k]))];
        var candidates = new int[n];
        var fresh = new bool[n];
        for (int t = 0; t < n; t++)
        {
            key.Clear();
            for (int k = 0; k < block.Count; k++)
            {
                int g = groups[k][t];
                if (g >= 0) AppendKey(key, "known", g.ToString(System.Globalization.CultureInfo.InvariantCulture));
                else
                {
                    fresh[t] = true;
                    // A new level is keyed by its name, so its cases share a cell; a grouping left out keys the case alone.
                    if (labels[k] is { } named) AppendKey(key, "new", named[t]);
                    else AppendKey(key, "case", t.ToString(System.Globalization.CultureInfo.InvariantCulture));
                }
            }
            candidates[t] = Number(numbers, key.ToString());
        }
        return new BassGroupCells(training, candidates, fresh, trainingCells);
    }

    // Length framing keeps arbitrary level names (including separators and type prefixes) distinct.
    private static void AppendKey(StringBuilder key, string kind, string value)
        => key.Append(kind).Append(':').Append(value.Length).Append(':').Append(value).Append('|');

    // The cell of `key`, numbered on its first appearance.
    private static int Number(Dictionary<string, int> numbers, string key)
    {
        if (numbers.TryGetValue(key, out int cell)) return cell;
        numbers.Add(key, numbers.Count);
        return numbers.Count - 1;
    }

    /// <summary>Refuses a stratum of <paramref name="stratumOf"/> whose members' cells differ: cases are exchangeable only within a cell.</summary>
    internal static void CheckStrata(IReadOnlyList<int> stratumOf, Func<int, IEnumerable<int>> cellsOf, string parameter)
    {
        var cellOf = new Dictionary<int, int>();
        for (int u = 0; u < stratumOf.Count; u++)
            foreach (int cell in cellsOf(u))
            {
                if (!cellOf.TryGetValue(stratumOf[u], out int held)) cellOf.Add(stratumOf[u], cell);
                else if (held != cell)
                    throw new ArgumentException($"Stratum {stratumOf[u]} spans cells {held} and {cell}: under random effects cases are exchangeable only "
                                                + "within a cell, one level of every grouping jointly, so each stratum lies inside one (BassFit.GroupCells).", parameter);
            }
    }
}
