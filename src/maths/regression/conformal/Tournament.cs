using System;
using System.Collections.Generic;

namespace Maths.Regression.Conformal;

public static partial class CrossConformal
{
    /// <summary>
    /// The tournament correction (Bhattacharyya, Zhang and Barber 2026, arXiv 2605.29200, <i>Approximating full conformal prediction:
    /// distribution free guarantees via the tournament correction</i>; the conformal layer's slice 13): full data at 1 − 2α whatever the
    /// fits' quality. Each training unit i plays one game against the candidate, both scored by one fit to the rest, D₋ᵢ, with the pair
    /// {zᵢ, candidate} handed over as values: <paramref name="unitScores"/>[i, c] is zᵢ's conformity there and
    /// <paramref name="candidateScores"/>[i, c] the candidate's, at candidate column c, higher conforming better. The candidate wins game i
    /// when it is stranger by more than δ = <paramref name="inflation"/>, its score plus δ below zᵢ's (<c>strangerThanBy</c>; δ = 0 is the
    /// plain rule, <c>strangerThan</c>, ties winning nothing), and it is accepted when its wins plus k fall below m = n + 1, at
    /// k = ⌊α·m⌋ (<c>tournament_count_iff</c>, <c>AcceptsPairwiseBy</c>).
    /// <list type="bullet">
    /// <item><b>Coverage</b> at least 1 − (2k − 1)/m, so at least 1 − 2α (<c>tournament_coverage</c>, Theorem 1), when each game's fit
    /// receives the pair as values, not marking which is the candidate (<c>tournFit_comp</c>), and its fit to the rest is symmetric in its
    /// training units given its randomness, or reads them in index order with the left-out place closed up and a seed that does not move
    /// with i (<c>tournament_orderReading</c>): a chain run on each D₋ᵢ from one seed qualifies, mixed or not. A fit that hands the pair
    /// over as (training, test) breaks it.</item>
    /// <item><b>The inflated set</b> (δ = 2ε) holds the plain one, so it keeps that bound; under their stability, each ordered pair's score
    /// straying from a reference by more than ε with chance at most ν, it covers at least <see cref="TournamentStableCoverage"/>
    /// (<c>coverage_tournament_stable</c>, Theorem 2).</item>
    /// </list>
    /// Each set is read as the cross-fitted form's (<see cref="ConformalForm.CrossFitted"/>): <c>Below</c> and <c>Tied</c> count the games
    /// whose zᵢ scores below and equal to the candidate's score plus δ, so the candidate's wins are the rest, and π = (Below + Tied + 1)/m.
    /// </summary>
    public static ConformalPrediction Tournament(double[,] unitScores, double[,] candidateScores, IReadOnlyList<double[]> grids,
                                                 double alpha = 0.1, double inflation = 0.0)
    {
        ArgumentNullException.ThrowIfNull(unitScores);
        ArgumentNullException.ThrowIfNull(candidateScores);
        ArgumentNullException.ThrowIfNull(grids);
        SplitConformal.CheckAlpha(alpha);
        if (!(inflation >= 0.0 && double.IsFinite(inflation))) throw new ArgumentOutOfRangeException(nameof(inflation), "δ must be non-negative and finite.");
        if (grids.Count == 0) throw new ArgumentException("At least one test case is required.", nameof(grids));
        int n = unitScores.GetLength(0);
        if (n == 0) throw new ArgumentException("At least one game is required.", nameof(unitScores));
        if (candidateScores.GetLength(0) != n || candidateScores.GetLength(1) != unitScores.GetLength(1))
            throw new ArgumentException("The candidate's scores must have the units' layout: one row per game, one column per candidate.", nameof(candidateScores));
        foreach (double v in unitScores)
            if (double.IsNaN(v)) throw new ArgumentException("A unit's score is NaN.", nameof(unitScores));
        foreach (double v in candidateScores)
            if (double.IsNaN(v)) throw new ArgumentException("A candidate's score is NaN.", nameof(candidateScores));
        double[][] copies = CopyGrids(grids, unitScores.GetLength(1));
        int m = checked(n + 1), k = (int)Math.Floor(alpha * m);
        var sets = new ConformalSet[copies.Length];
        int column = 0;
        for (int t = 0; t < sets.Length; t++)
        {
            double[] grid = copies[t];
            int first = column;
            column += grid.Length;
            if (k == 0)
            {
                sets[t] = new(t, 0, m, k, grid, null, null, null, null, null, Inflated: inflation > 0.0, Form: ConformalForm.CrossFitted);
                continue;
            }
            var below = new int[grid.Length];
            var tied = new int[grid.Length];
            var p = new double[grid.Length];
            for (int i = 0; i < n; i++)
                for (int g = 0; g < grid.Length; g++)
                {
                    double unit = unitScores[i, first + g], raised = candidateScores[i, first + g] + inflation;
                    if (unit < raised) below[g]++;
                    else if (unit == raised) tied[g]++;
                }
            for (int g = 0; g < grid.Length; g++) p[g] = (below[g] + tied[g] + 1) / (double)m;
            sets[t] = new(t, 0, m, k, grid, p, below, tied, null, null, Inflated: inflation > 0.0, Form: ConformalForm.CrossFitted);
        }
        return new(sets, alpha, ConformalForm.CrossFitted);
    }

    /// <summary>
    /// The tournament from n leave-one-out fits' draws, read one fit at a time (<paramref name="read"/>): fit i's held-out seam holds
    /// ln p(zᵢ | θ_s) on its one observed value, and its candidate seam ln p(y | θ_s) on the test grids, over the same draws of a fit to
    /// D₋ᵢ. Game i's fit to the pair reweights those draws by the pair's likelihoods, w̃_s ∝ W_s·p(zᵢ | θ_s)·p(y | θ_s), the two terms
    /// summed first, so the weight is one function of the pair as values; each member then scores ln Σ_s w̃_s p(· | θ_s), its in-sample
    /// predictive under the fit to D₋ᵢ plus the pair, by one function with the roles swapped (<see cref="PairScore"/>). Where no draw gives
    /// the pair positive density, each member reads its own predictive under the fit to the rest. The draws need not mix and the fits need
    /// not be exact: the guarantee needs the fits to the rest to be run alike for every i, from one seed (see the overload on scores). The
    /// seams are checked as the split form's: one exchangeable group in stratum 0.
    /// </summary>
    public static ConformalPrediction Tournament(int fits, Func<int, (CandidateLogDensity HeldOut, CandidateLogDensity Candidates)> read,
                                                 double alpha = 0.1, double inflation = 0.0)
    {
        ArgumentNullException.ThrowIfNull(read);
        ArgumentOutOfRangeException.ThrowIfLessThan(fits, 1);
        SplitConformal.CheckAlpha(alpha);
        double[,]? units = null, candidates = null;
        IReadOnlyList<double[]>? grids = null;
        double[] held = [], pair = [];
        for (int i = 0; i < fits; i++)
        {
            (CandidateLogDensity at, CandidateLogDensity test) = read(i);
            SplitConformal.CheckCalibration(at);
            if (at.Cases != 1) throw new ArgumentException("Each fit must score its one held-out unit.", nameof(read));
            new SplitDraws(at).Check(test);
            if (at.StratumOf[0] != 0 || Array.Exists(test.StratumOf, h => h != 0))
                throw new ArgumentException("Supply one exchangeable group in stratum zero.", nameof(read));
            grids ??= test.Grids;
            int columns = test.Values.GetLength(1);
            units ??= new double[fits, columns];
            candidates ??= new double[fits, columns];
            if (test.Cases != grids.Count || columns != units.GetLength(1)) throw new ArgumentException("Every fit must score the same test cases.", nameof(read));
            for (int t = 0; t < grids.Count; t++)
                if (!ConformalBayes.SameSequence(grids[t], test.Grids[t]))
                    throw new ArgumentException("Every fit must score the same test grids.", nameof(read));
            (int[] rows, double[]? weights) = PredictiveChecks.Rows(at.Weights, at.Draws);
            int s = rows.Length;
            double[] logWeights = weights ?? new double[s];
            if (held.Length != s)
            {
                held = new double[s];
                pair = new double[s];
            }
            for (int d = 0; d < s; d++) held[d] = at.Values[rows[d], 0];
            for (int c = 0; c < columns; c++)
            {
                for (int d = 0; d < s; d++) pair[d] = test.Values[rows[d], c];
                units[i, c] = PairScore(logWeights, held, pair);
                candidates[i, c] = PairScore(logWeights, pair, held);
            }
        }
        return Tournament(units!, candidates!, grids!, alpha, inflation);
    }

    /// <summary>
    /// The conformity of <paramref name="scored"/> under one fit to the rest, <paramref name="logWeights"/> over its draws, plus the pair
    /// {<paramref name="scored"/>, <paramref name="other"/>}, each given by its log density per draw: ln Σ_s w̃_s e^(scored_s) with
    /// ln w̃_s = ln W_s + (scored_s + other_s) normalized, the pair's two terms summed before the weight is added, so swapping the arguments
    /// swaps the two members' scores bit for bit and the fit cannot tell which is the candidate. Where no draw gives the pair positive density,
    /// the scored member's own predictive under the fit to the rest, ln Σ_s W_s e^(scored_s) − ln Σ_s W_s, again one function of the pair.
    /// </summary>
    internal static double PairScore(ReadOnlySpan<double> logWeights, ReadOnlySpan<double> scored, ReadOnlySpan<double> other)
    {
        int s = scored.Length;
        double top = double.NegativeInfinity;
        for (int d = 0; d < s; d++) top = Math.Max(top, logWeights[d] + (scored[d] + other[d]));
        if (double.IsNegativeInfinity(top)) return LogSum(logWeights, scored) - LogSum(logWeights, default);
        double total = 0.0, kept = 0.0, best = double.NegativeInfinity;
        for (int d = 0; d < s; d++) best = Math.Max(best, logWeights[d] + (scored[d] + other[d]) + scored[d]);
        for (int d = 0; d < s; d++)
        {
            double l = logWeights[d] + (scored[d] + other[d]);
            total += Math.Exp(l - top);
            kept += Math.Exp(l + scored[d] - best);
        }
        return best + Math.Log(kept) - (top + Math.Log(total));
    }

    // ln Σ_s e^(logWeights_s + terms_s), or of the weights alone where the terms are not given.
    private static double LogSum(ReadOnlySpan<double> logWeights, ReadOnlySpan<double> terms)
    {
        double top = double.NegativeInfinity;
        for (int d = 0; d < logWeights.Length; d++) top = Math.Max(top, logWeights[d] + (terms.IsEmpty ? 0.0 : terms[d]));
        if (double.IsNegativeInfinity(top)) return top;
        double sum = 0.0;
        for (int d = 0; d < logWeights.Length; d++) sum += Math.Exp(logWeights[d] + (terms.IsEmpty ? 0.0 : terms[d]) - top);
        return top + Math.Log(sum);
    }

    /// <summary>
    /// The inflated tournament set's coverage under their stability (<c>coverage_tournament_stable</c>, <c>card_refused_by_le</c>): with m
    /// units at level k, each ordered pair's score straying from a reference score of the data set by more than ε with chance at most ν under
    /// an exchangeable law, the 2ε-inflated set covers at least 1 − min over a ≥ 0 of [(k + a)/m + 2(m − 1)ν/(a + 1)], the Markov count's
    /// threshold optimized. Their form fixes it near 2√ν·m and reads 1 − α − 3√ν: 0.673 here against their 0.600 at m 20, k 2 (α 0.1),
    /// ν 0.01, the excess over α 0.227 against 0.300. The reference need be neither symmetric nor a valid score. The step from a data set to the law is the averaging every
    /// CB statement uses, not formalized in the Lean.
    /// </summary>
    public static double TournamentStableCoverage(int units, int level, double strayChance)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(units, 1);
        if (level < 0 || level > units) throw new ArgumentOutOfRangeException(nameof(level), "The level lies in 0 … m.");
        if (!(strayChance >= 0.0 && strayChance <= 1.0)) throw new ArgumentOutOfRangeException(nameof(strayChance), "ν lies in [0, 1].");
        double best = double.PositiveInfinity;
        for (int a = 0; a <= units; a++)
            best = Math.Min(best, (level + a) / (double)units + 2.0 * (units - 1) * strayChance / (a + 1));
        return Math.Max(0.0, 1.0 - best);
    }

    // Each grid finite and strictly ascending, and the grids covering the scores' columns.
    private static double[][] CopyGrids(IReadOnlyList<double[]> grids, int columns)
    {
        int total = 0;
        var copies = new double[grids.Count][];
        for (int t = 0; t < grids.Count; t++)
        {
            double[] grid = grids[t] ?? throw new ArgumentException("A grid is null.", nameof(grids));
            if (grid.Length == 0) throw new ArgumentException("A grid is empty.", nameof(grids));
            for (int g = 0; g < grid.Length; g++)
                if (!double.IsFinite(grid[g]) || (g > 0 && !(grid[g] > grid[g - 1])))
                    throw new ArgumentException("Grid values must be finite and strictly ascending.", nameof(grids));
            copies[t] = (double[])grid.Clone();
            total = checked(total + grid.Length);
        }
        if (total != columns) throw new ArgumentException("The grids must cover all score columns.", nameof(grids));
        return copies;
    }
}
