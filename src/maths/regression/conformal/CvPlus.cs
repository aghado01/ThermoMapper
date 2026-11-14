using System;
using System.Collections.Generic;

namespace Maths.Regression.Conformal;

/// <summary>Assumptions supplied by the caller for the CV+ bounds. The layer validates dimensions and folds,
/// while the caller establishes the data law and the fits' behavior.</summary>
[Flags]
public enum CvPlusAssumptions
{
    None = 0,

    /// <summary>Each fit reads neither its held-out fold nor the test unit, including preprocessing and tuning.</summary>
    ExcludesFoldAndTest = 1,

    /// <summary>Each fold with the test unit is exchangeable given the rest and the fitting randomness.</summary>
    FoldTestExchangeability = 2,

    /// <summary>The joint law of n + m units and fitting randomness is invariant under fold and slot permutations.
    /// I.i.d. units with independent fitting randomness supply this extended law.</summary>
    ExtendedFoldExchangeability = 4,

    /// <summary>Fits are symmetric under the fold and slot permutations, up to the same relabeling of their randomness.</summary>
    SymmetricFit = 8,

    MergedRanks = ExcludesFoldAndTest | FoldTestExchangeability,
    FoldTournament = ExtendedFoldExchangeability | SymmetricFit,
    Both = MergedRanks | FoldTournament,
}

/// <summary>Responses and their own fold's held-out predictions, with unit-to-fold labels starting at zero.
/// Predictions has one row per fold fit and one column per test input. Every fold has the same positive size.</summary>
public sealed record FoldReadings(IReadOnlyList<double> Response, IReadOnlyList<double> HeldOut,
                                  IReadOnlyList<int> FoldOf, double[,] Predictions);

/// <summary>A conditional lower bound with its theorem, required assumptions and whether those assumptions were supplied.</summary>
public sealed record CvPlusBound(double Value, string Theorem, CvPlusAssumptions RequiredAssumptions, bool Applies);

/// <summary>What shift weights for K-fold CV+ state (slice 14's follow-up 2), which decides the bounds it reports.</summary>
public enum CvPlusWeighting
{
    /// <summary>
    /// Likelihood ratios w(z) travelling with the units, the test unit's its own (JAWS's setting): w(x) under covariate shift. Three bounds,
    /// the least miss reported:
    /// <list type="bullet">
    /// <item>one fold's, min(1, α·W_tot/(W_k + w_t)) missed at most, read at the fold of largest weight (<c>cvPlus_oneFold_likelihoodRatio</c>),
    /// which serves two folds;</item>
    /// <item>across folds, 2α + W_max/W, W_max bounding every fold's weight, the test's fold-mates the interval never reads included, so it
    /// needs a bound on the ratio (<c>coverage_across_likelihoodRatio</c>), which serves many folds of few units;</item>
    /// <item>merged, the least of min(1, 2c(w_max)), c(γ) = (αW + (K − 1)γ)/(W + (K − 1)γ), W the bag's weight with the test point's and
    /// w_max its largest weight (<c>cvPlus_merged_likelihoodRatio</c>), and the refined 2L(W + (K − 1)e)/N², L = αW + (K − 1)γ,
    /// N = W + (K − 1)γ, e = Σw²/W, at γ the largest weight at most αW where 2L ≤ N, 0 where no weight is (no unit can then be refused)
    /// (<c>cvPlus_merged_likelihoodRatio_refined</c>); it serves many units.</item>
    /// </list>
    /// </summary>
    LikelihoodRatio,

    /// <summary>
    /// Weights fixed by position, the test position's the largest (Barber, Candès, Ramdas and Tibshirani 2023): the merged bound 2c,
    /// c = (α(W + w_t) + (K − 1)·w_t)/(W + K·w_t), given each fold's weighted p-value super-uniform, which exchangeable data supply
    /// (<c>cvPlus_merged_weighted</c>); at equal weights slice 9's merged bound.
    /// </summary>
    Position,
}

/// <summary>Shift weights for K-fold CV+ (<see cref="CrossConformal.CvPlus"/>): <paramref name="Weights"/>, one per training unit (or per
/// test case) and one per test case; what they state, <paramref name="Weighting"/>; and, for likelihood ratios, a bound w̄ on the ratio,
/// <paramref name="RatioBound"/>, which the across-folds bound reads as W_max ≤ max(the folds' weights, fold size × w̄).</summary>
public sealed record CvPlusWeights(ShiftWeights Weights, CvPlusWeighting Weighting, double? RatioBound = null);

/// <summary>The finite CV+ bounds. Coverage selects the largest applicable bound, or null if none applies.
/// At level zero the interval is the whole line and coverage is one without further assumptions.</summary>
public sealed class CvPlusBounds
{
    private readonly bool _weighted;

    // A weighted case's bounds (slice 14's follow-up 2), each a share of the weight, with no level shortcut.
    private CvPlusBounds(int trainingUnits, int folds, int foldSize, int level, CvPlusAssumptions assumptions,
                         CvPlusBound merged, CvPlusBound oneFold, CvPlusBound tournament)
    {
        TrainingUnits = trainingUnits;
        Folds = folds;
        FoldSize = foldSize;
        Level = level;
        Assumptions = assumptions;
        Merged = merged;
        OneFold = oneFold;
        Tournament = tournament;
        _weighted = true;
    }

    /// <summary>A weighted case's bounds from its fold weights W_k, its unit weights and its test weight w_t.</summary>
    internal static CvPlusBounds Weighted(int trainingUnits, int folds, int foldSize, int level, CvPlusAssumptions assumptions, double alpha,
                                          double[] foldWeights, IReadOnlyList<double> unitWeights, double testWeight, CvPlusWeights weights)
    {
        double training = 0.0, heaviest = 0.0, largestUnit = 0.0;
        foreach (double w in foldWeights)
        {
            training += w;
            heaviest = Math.Max(heaviest, w);
        }
        for (int i = 0; i < unitWeights.Count; i++) largestUnit = Math.Max(largestUnit, unitWeights[i]);
        double total = training + testWeight;
        bool ranks = (assumptions & CvPlusAssumptions.MergedRanks) == CvPlusAssumptions.MergedRanks;
        bool tournament = (assumptions & CvPlusAssumptions.FoldTournament) == CvPlusAssumptions.FoldTournament;
        CvPlusBound none(string theorem, CvPlusAssumptions required) => new(double.NaN, theorem, required, false);
        if (weights.Weighting == CvPlusWeighting.Position)
        {
            double c = (alpha * total + (folds - 1) * testWeight) / (training + folds * testWeight);
            var merged = new CvPlusBound(Math.Max(0.0, 1.0 - 2.0 * c), "cvPlus_merged_weighted", CvPlusAssumptions.MergedRanks,
                                         ranks && testWeight >= largestUnit);
            return new(trainingUnits, folds, foldSize, level, assumptions, merged,
                       none("cvPlus_oneFold_likelihoodRatio", CvPlusAssumptions.MergedRanks), none("coverage_across_likelihoodRatio", CvPlusAssumptions.FoldTournament));
        }
        double share = heaviest + testWeight > 0.0 ? Math.Min(1.0, alpha * total / (heaviest + testWeight)) : 1.0;
        var oneFold = new CvPlusBound(1.0 - share, "cvPlus_oneFold_likelihoodRatio", CvPlusAssumptions.MergedRanks, ranks);
        CvPlusBound across = weights.RatioBound is double bound && total > 0.0
            ? new(Math.Max(0.0, 1.0 - 2.0 * alpha - Math.Max(heaviest, foldSize * bound) / total), "coverage_across_likelihoodRatio",
                  CvPlusAssumptions.FoldTournament, tournament)
            : none("coverage_across_likelihoodRatio", CvPlusAssumptions.FoldTournament);
        var mergedRatio = new CvPlusBound(1.0 - MergedMiss(alpha, folds, unitWeights, testWeight, total), "cvPlus_merged_likelihoodRatio",
                                          CvPlusAssumptions.MergedRanks, ranks);
        return new(trainingUnits, folds, foldSize, level, assumptions, mergedRatio, oneFold, across);
    }

    // The merged bound's miss under likelihood ratios: the least of min(1, 2c(w_max)) and the refined form, over the bag's weights, the test
    // point's included, W their total.
    private static double MergedMiss(double alpha, int folds, IReadOnlyList<double> unitWeights, double testWeight, double total)
    {
        if (!(total > 0.0)) return 0.0;
        double largest = testWeight, light = testWeight <= alpha * total ? testWeight : double.NegativeInfinity, squares = testWeight * testWeight;
        for (int i = 0; i < unitWeights.Count; i++)
        {
            double w = unitWeights[i];
            largest = Math.Max(largest, w);
            if (w <= alpha * total) light = Math.Max(light, w);
            squares += w * w;
        }
        double extra = folds - 1;
        double simple = Math.Min(1.0, 2.0 * (alpha * total + extra * largest) / (total + extra * largest));
        if (double.IsNegativeInfinity(light)) return 0.0;   // no weight at most αW: no unit can be refused
        double l = alpha * total + extra * light, n = total + extra * light, mean = squares / total;
        double refined = 2.0 * l <= n ? 2.0 * l * (total + extra * mean) / (n * n) : 1.0;
        return Math.Min(simple, refined);
    }

    internal CvPlusBounds(int trainingUnits, int folds, int foldSize, int level, CvPlusAssumptions assumptions)
    {
        TrainingUnits = trainingUnits;
        Folds = folds;
        FoldSize = foldSize;
        Level = level;
        Assumptions = assumptions;
        Merged = new(Math.Max(0, 1 - 2.0 * (level + (double)folds - 1) / (trainingUnits + (double)folds)),
                     "cvPlus_merged", CvPlusAssumptions.MergedRanks,
                     (assumptions & CvPlusAssumptions.MergedRanks) == CvPlusAssumptions.MergedRanks);
        OneFold = new(Math.Max(0, 1 - level / (foldSize + 1.0)), "cvPlus_oneFold", CvPlusAssumptions.MergedRanks,
                      (assumptions & CvPlusAssumptions.MergedRanks) == CvPlusAssumptions.MergedRanks);
        Tournament = new(Math.Max(0, 1 - Math.Max(0, 2.0 * level + foldSize - 2) / (trainingUnits + (double)foldSize)),
                         "cvPlus_tournament", CvPlusAssumptions.FoldTournament,
                         (assumptions & CvPlusAssumptions.FoldTournament) == CvPlusAssumptions.FoldTournament);
    }

    public int TrainingUnits { get; }
    public int Folds { get; }
    public int FoldSize { get; }
    public int Level { get; }
    public CvPlusAssumptions Assumptions { get; }
    public CvPlusBound Merged { get; }
    public CvPlusBound OneFold { get; }
    public CvPlusBound Tournament { get; }

    /// <summary>The largest applicable bound, null if none applies. Under weights each bound is a share of the weight, read per case.</summary>
    public double? Coverage
    {
        get
        {
            if (Level == 0 && !_weighted) return 1.0;
            double? coverage = Merged.Applies ? Merged.Value : null;
            if (OneFold.Applies) coverage = Math.Max(coverage ?? 0, OneFold.Value);
            if (Tournament.Applies) coverage = Math.Max(coverage ?? 0, Tournament.Value);
            return coverage;
        }
    }
}

/// <summary>One closed CV+ interval, with the assumption-specific coverage bounds: a weighted case's own (slice 14's follow-up 2).</summary>
public sealed record CvPlusInterval(int Case, double Lower, double Upper, CvPlusBounds Bounds)
{
    /// <summary>The weighted rule's α where the interval is weighted CV+'s; null for the counts' interval, equal weights' included.</summary>
    public double? WeightedAlpha { get; init; }

    /// <summary>Whether the interval is the whole line: level 0 for the counts' interval; both ends unbounded for the weighted one.</summary>
    public bool WholeSpace => WeightedAlpha is null ? Bounds.Level == 0 : double.IsNegativeInfinity(Lower) && double.IsPositiveInfinity(Upper);
    public bool Empty => Lower > Upper;
    public double Width => Empty ? 0.0 : Upper - Lower;
    public bool Contains(double y) => y >= Lower && y <= Upper;
    public double? Coverage => Bounds.Coverage;
}

/// <summary>The CV+ intervals, with the counts' bounds; under weights each interval carries its own (<see cref="CvPlusInterval.Bounds"/>).</summary>
public sealed record CvPlusIntervals(IReadOnlyList<CvPlusInterval> Intervals, double Alpha, CvPlusBounds Bounds);

public static partial class CrossConformal
{
    /// <summary>CV+ closed endpoints (<c>CVRefuses</c>), at ℓ = ⌊α(n + 1)⌋. Unit i uses fold FoldOf[i]'s prediction
    /// and its held-out absolute residual. Folds are nonempty and equal-sized. At ℓ = 0 the interval is the whole line.
    /// Each conditional bound is enabled only by its stated assumptions. Fits and their preprocessing come from the caller.
    /// With <paramref name="weights"/> (slice 14's follow-up 2) a case whose weights are not all one value takes weighted CV+, the
    /// jackknife+'s weighted ends on the folds' predictions (<c>JackRefusesW</c>), and its own bounds, which
    /// <see cref="CvPlusWeights.Weighting"/> decides.</summary>
    public static CvPlusIntervals CvPlus(FoldReadings readings, CvPlusAssumptions assumptions, double alpha = 0.1, CvPlusWeights? weights = null)
    {
        ArgumentNullException.ThrowIfNull(readings);
        ArgumentNullException.ThrowIfNull(readings.Response);
        ArgumentNullException.ThrowIfNull(readings.HeldOut);
        ArgumentNullException.ThrowIfNull(readings.FoldOf);
        ArgumentNullException.ThrowIfNull(readings.Predictions);
        SplitConformal.CheckAlpha(alpha);
        if ((assumptions & ~CvPlusAssumptions.Both) != 0) throw new ArgumentOutOfRangeException(nameof(assumptions));
        int n = readings.Response.Count, folds = readings.Predictions.GetLength(0);
        if (n == 0 || folds == 0 || n % folds != 0 || readings.FoldOf.Count != n)
            throw new ArgumentException("Nonempty, equally sized folds with one label per response are required.", nameof(readings));
        int units = checked(n + 1);
        int size = n / folds;
        var counts = new int[folds];
        foreach (int fold in readings.FoldOf)
        {
            if ((uint)fold >= (uint)folds) throw new ArgumentException("Fold labels must index the prediction rows.", nameof(readings));
            counts[fold]++;
        }
        foreach (int count in counts)
            if (count != size) throw new ArgumentException("Every fold must have the same positive size.", nameof(readings));
        double[] residuals = AbsoluteResiduals(readings.Response, readings.HeldOut);
        foreach (double prediction in readings.Predictions)
            if (!double.IsFinite(prediction)) throw new ArgumentException("Predictions must be finite.", nameof(readings));
        int level = (int)Math.Floor(alpha * units);
        var bounds = new CvPlusBounds(n, folds, size, level, assumptions);
        int cases = readings.Predictions.GetLength(1);
        if (weights is not null)
        {
            ArgumentNullException.ThrowIfNull(weights.Weights, nameof(weights));
            if (!Enum.IsDefined(weights.Weighting)) throw new ArgumentOutOfRangeException(nameof(weights));
            if (weights.RatioBound is double bound && !(bound > 0.0 && double.IsFinite(bound)))
                throw new ArgumentOutOfRangeException(nameof(weights), "The ratio's bound is positive and finite.");
            if (weights.RatioBound is not null && weights.Weighting != CvPlusWeighting.LikelihoodRatio)
                throw new ArgumentException("A bound on the ratio serves likelihood-ratio weights.", nameof(weights));
            CheckWeights(weights.Weights, n, cases, cases, nameof(weights));
        }
        var intervals = new CvPlusInterval[cases];
        var lower = new double[n];
        var upper = new double[n];
        for (int t = 0; t < cases; t++)
        {
            bool weighed = weights is not null && !Equal(weights.Weights.UnitsFor(t), weights.Weights.Candidates, t, 1);
            if (level == 0 && !weighed)
            {
                intervals[t] = new(t, double.NegativeInfinity, double.PositiveInfinity, bounds);
                continue;
            }
            for (int i = 0; i < n; i++)
            {
                double prediction = readings.Predictions[readings.FoldOf[i], t];
                lower[i] = prediction - residuals[i];
                upper[i] = prediction + residuals[i];
            }
            if (weighed)
            {
                IReadOnlyList<double> unitWeights = weights!.Weights.UnitsFor(t);
                double testWeight = weights.Weights.Candidates[t];
                var foldWeights = new double[folds];
                for (int i = 0; i < n; i++) foldWeights[readings.FoldOf[i]] += unitWeights[i];
                GridInterval w = WeightedEnds(lower, upper, unitWeights, testWeight, alpha);
                intervals[t] = new(t, w.Lower, w.Upper, CvPlusBounds.Weighted(n, folds, size, level, assumptions, alpha, foldWeights, unitWeights, testWeight, weights))
                {
                    WeightedAlpha = alpha,
                };
                continue;
            }
            GridInterval ends = OrderedEnds(lower, upper, level);
            intervals[t] = new(t, ends.Lower, ends.Upper, bounds);
        }
        return new(intervals, alpha, bounds);
    }
}
