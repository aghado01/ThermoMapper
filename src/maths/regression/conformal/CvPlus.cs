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

/// <summary>The finite CV+ bounds. Coverage selects the largest applicable bound, or null if none applies.
/// At level zero the interval is the whole line and coverage is one without further assumptions.</summary>
public sealed class CvPlusBounds
{
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

    public double? Coverage
    {
        get
        {
            if (Level == 0) return 1.0;
            double? coverage = Merged.Applies ? Merged.Value : null;
            if (OneFold.Applies) coverage = Math.Max(coverage ?? 0, OneFold.Value);
            if (Tournament.Applies) coverage = Math.Max(coverage ?? 0, Tournament.Value);
            return coverage;
        }
    }
}

/// <summary>One closed CV+ interval, with the assumption-specific coverage bounds.</summary>
public sealed record CvPlusInterval(int Case, double Lower, double Upper, CvPlusBounds Bounds)
{
    public bool WholeSpace => Bounds.Level == 0;
    public bool Empty => Lower > Upper;
    public double Width => Empty ? 0.0 : Upper - Lower;
    public bool Contains(double y) => y >= Lower && y <= Upper;
    public double? Coverage => Bounds.Coverage;
}

public sealed record CvPlusIntervals(IReadOnlyList<CvPlusInterval> Intervals, double Alpha, CvPlusBounds Bounds);

public static partial class CrossConformal
{
    /// <summary>CV+ closed endpoints (<c>CVRefuses</c>), at ℓ = ⌊α(n + 1)⌋. Unit i uses fold FoldOf[i]'s prediction
    /// and its held-out absolute residual. Folds are nonempty and equal-sized. At ℓ = 0 the interval is the whole line.
    /// Each conditional bound is enabled only by its stated assumptions. Fits and their preprocessing come from the caller.</summary>
    public static CvPlusIntervals CvPlus(FoldReadings readings, CvPlusAssumptions assumptions, double alpha = 0.1)
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
        var intervals = new CvPlusInterval[cases];
        var lower = new double[n];
        var upper = new double[n];
        for (int t = 0; t < cases; t++)
        {
            if (level == 0)
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
            GridInterval ends = OrderedEnds(lower, upper, level);
            intervals[t] = new(t, ends.Lower, ends.Upper, bounds);
        }
        return new(intervals, alpha, bounds);
    }
}
