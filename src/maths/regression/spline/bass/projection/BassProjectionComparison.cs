using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// What a projected refit is compared with its native fit on (E8, slice 7c): held-out cases, cross-fitting's folds (the user's
/// decision 5, opt-in), how each fold's frame is read, and the refit's readout for its consistency with the frame. A declarative
/// record.
/// </summary>
public sealed record ProjectionComparisonSpec
{
    /// <summary>Held-out cases in the original units; none when null.</summary>
    public PredictionCases? TestCases { get; init; }

    /// <summary>The held-out responses, cases × output locations (one column for a scalar response), with <see cref="TestCases"/>.</summary>
    public double[,]? TestResponse { get; init; }

    /// <summary>Cross-fitting's folds K: 0, the default, for none; otherwise at least 2. It costs K native fits and K refits.</summary>
    public int Folds { get; init; }

    /// <summary>The seed of the folds' assignment.</summary>
    public int FoldSeed { get; init; } = 1;

    /// <summary>How each fold's frame is read (the law and the metric), as the compared analysis was; the default spec when null.</summary>
    public ActiveSubspaceSpec? Analysis { get; init; }

    /// <summary>The refit's readout points, for its consistency with the frame.</summary>
    public int Points { get; init; } = 4096;

    /// <summary>The readout points' seed.</summary>
    public int Seed { get; init; } = 1;
}

/// <summary>One fit's side of the comparison (slice 7c).</summary>
/// <param name="BasisFunctions">The mean number of basis functions per draw, under the draws' weights.</param>
/// <param name="Loo">PSIS-LOO on the fit's own data (<see cref="PredictiveChecks.PsisLoo"/>), in-sample: see the comparison's caveat.</param>
/// <param name="Evidence">ln p(y) when the run was tempered to an exact reference rung (<see cref="BassFit.Evidence"/>); null otherwise.</param>
/// <param name="HeldOutRmse">The posterior mean's RMSE at the held-out cases; null without them.</param>
/// <param name="HeldOutExplained">100·(1 − RMSE/sd) of the held-out responses, the percentage of their sd accounted for (2307.11241's
/// measure, its definition not stated there); null without them.</param>
public sealed record ProjectionFitReport(double BasisFunctions, PredictiveCheck Loo, BassEvidence? Evidence, double? HeldOutRmse, double? HeldOutExplained);

/// <summary>
/// Cross-fitting (slice 7c, decision 5): over K disjoint folds, a native fit on the other folds gives each fold's frame and refit,
/// and both are scored on the fold held out, so no score reads a frame its cases helped estimate. Every fold's fit scales each input by
/// the full data's box, so the frames share their coordinates.
/// </summary>
/// <param name="Folds">Each fold's cases, ascending.</param>
/// <param name="Frames">Each fold's frame, read from the native fit on the other folds.</param>
/// <param name="NativeRmse">Each fold's native RMSE on its held-out cases.</param>
/// <param name="RefitRmse">Each fold's refit RMSE on its held-out cases.</param>
/// <param name="NativePooled">The native RMSE over every held-out case.</param>
/// <param name="RefitPooled">The refit RMSE over every held-out case.</param>
/// <param name="NativeExplained">100·(1 − pooled RMSE/sd of the responses), native.</param>
/// <param name="RefitExplained">The same for the refit.</param>
public sealed record CrossFitReport(IReadOnlyList<int[]> Folds, IReadOnlyList<BassProjection> Frames, double[] NativeRmse, double[] RefitRmse,
                                   double NativePooled, double RefitPooled, double NativeExplained, double RefitExplained);

/// <summary>
/// A projected refit beside its native fit (E8, slice 7c): each fit's size, in-sample criteria and held-out error; the refit's own
/// subspace, mapped to the original inputs, against the frame; the native draws' spread about the frame, which the refit, taking the
/// frame as a plug-in, does not carry; and cross-fitting when asked. <see cref="Caveat"/> is in every report.
/// </summary>
/// <param name="Native">The native fit's side.</param>
/// <param name="Refit">The refit's side.</param>
/// <param name="RefitToFrame">For k = 1 to d, the largest principal angle between the refit's mapped E[C]'s leading k eigenvectors and
/// the frame's leading k, in radians: 0 at k = d whenever the refit uses every active variable; below d it reads whether the refit
/// orders the directions as the frame does, which a small eigengap leaves arbitrary.</param>
/// <param name="FrameSpread">Each native draw's d principal angles between its leading d eigenvectors and the frame, ascending.</param>
/// <param name="CrossFitted">Cross-fitting's scores; null unless asked.</param>
public sealed record BassProjectionComparison(ProjectionFitReport Native, ProjectionFitReport Refit, IReadOnlyList<double> RefitToFrame,
                                              double[][] FrameSpread, CrossFitReport? CrossFitted)
{
    /// <summary>The data-reuse caveat (the user's decision 5): in every report.</summary>
    public string Caveat => "The frame was read from the same data the refit is fitted and scored on, so the in-sample criteria (PSIS-LOO and the " +
                            "evidence) favour the refit; only the cross-fitted scores read each fold's frame from the other folds alone.";
}

/// <summary>The comparison's computation (slice 7c).</summary>
internal static class BassProjectionComparisons
{
    public static BassProjectionComparison Of(BassProjectedFit refit, BassFit native, BassData data, ActiveSubspaceAnalysis analysis, ProjectionComparisonSpec spec)
    {
        ArgumentNullException.ThrowIfNull(native);
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(analysis);
        ArgumentNullException.ThrowIfNull(spec);
        BassProjection projection = refit.Projection;
        int p = projection.Inputs.Count, d = projection.Dimension;
        for (int i = 0; i < p; i++)
            for (int j = 0; j < d; j++)
                if (analysis.Posterior.MeanCVectors[i, j] != projection.Frame[i, j])
                    throw new ArgumentException("The analysis is not the one the projection's frame was read from.", nameof(analysis));
        if ((spec.TestCases is null) != (spec.TestResponse is null))
            throw new ArgumentException("Held-out cases and their responses come together.", nameof(spec));
        if (spec.Folds is 1 or < 0) throw new ArgumentOutOfRangeException(nameof(spec), "Cross-fitting needs at least two folds; 0 asks for none.");

        ProjectionFitReport nativeSide = Side(native, data, spec.TestCases is null ? null : native.Predict(data, spec.TestCases).Mean(), spec.TestResponse);
        ProjectionFitReport refitSide = Side(refit.Fit, refit.Data, spec.TestCases is null ? null : refit.Predict(spec.TestCases).Mean(), spec.TestResponse);

        ProjectedActiveSubspace read = refit.ActiveSubspace(new ActiveSubspaceSpec { Points = spec.Points, Seed = spec.Seed });
        double[] toFrame = [.. Enumerable.Range(1, d).Select(k => SubspacePosterior.PrincipalAngles(Leading(read.InputSpace.MeanCVectors, k), Leading(projection.Frame, k)).Max())];
        double[][] spread = [.. analysis.Posterior.Draws.Select(draw =>
            SubspacePosterior.PrincipalAngles(Leading(draw.Vectors, Math.Min(d, draw.Vectors.GetLength(1))), projection.Frame))];
        CrossFitReport? cross = spec.Folds == 0 ? null : CrossFit(refit, native, data, spec);
        return new BassProjectionComparison(nativeSide, refitSide, Array.AsReadOnly(toFrame), spread, cross);
    }

    private static ProjectionFitReport Side(BassFit fit, BassData data, double[,]? mean, double[,]? response)
    {
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        IReadOnlyList<double>? w = fit.DrawWeights;
        double total = w?.Sum() ?? pooled.Length;
        double terms = Enumerable.Range(0, pooled.Length).Sum(s => (w?[s] ?? 1.0) * pooled[s].Structure.TermIds.Length) / total;
        PredictiveCheck loo = PredictiveChecks.PsisLoo(fit.PointwiseLogDensity(data));
        double? rmse = null, explained = null;
        if (mean is not null) (rmse, explained) = Score(mean, response!);
        return new ProjectionFitReport(terms, loo, fit.Evidence, rmse, explained);
    }

    private static (double Rmse, double Explained) Score(double[,] mean, double[,] response)
    {
        if (mean.GetLength(0) != response.GetLength(0) || mean.GetLength(1) != response.GetLength(1))
            throw new ArgumentException($"The held-out responses are {mean.GetLength(0)} × {mean.GetLength(1)}, one row per case and one column per location.");
        double sse = 0.0, sum = 0.0;
        int count = response.Length;
        foreach (double y in response) sum += y;
        double average = sum / count, ss = 0.0;
        for (int r = 0; r < response.GetLength(0); r++)
            for (int l = 0; l < response.GetLength(1); l++)
            {
                double e = mean[r, l] - response[r, l], c = response[r, l] - average;
                sse += e * e;
                ss += c * c;
            }
        double rmse = Math.Sqrt(sse / count), sd = Math.Sqrt(ss / (count - 1));
        return (rmse, 100.0 * (1.0 - rmse / sd));
    }

    // K folds: the cases shuffled by the seed and cut into K parts of sizes differing by at most one; each fold's frame from a native
    // fit on the other folds, its refit on them, and both scored on it. Every fold's data take the full data's boxes as their domains.
    private static CrossFitReport CrossFit(BassProjectedFit refit, BassFit native, BassData data, ProjectionComparisonSpec spec)
    {
        if (data.Groupings.Count > 0 || data.Weights is not null)
            throw new ArgumentException("Cross-fitting reads data without groupings or case weights; subsetting those is not built.", nameof(data));
        int n = data.Response.Cases, k = spec.Folds;
        if (k > n) throw new ArgumentOutOfRangeException(nameof(spec), $"There are {n} cases, so at most {n} folds.");
        int[] order = [.. Enumerable.Range(0, n)];
        var rng = new Xoshiro256PlusPlus(spec.FoldSeed);
        for (int i = n - 1; i > 0; i--)
        {
            int j = rng.NextInt(i + 1);
            (order[i], order[j]) = (order[j], order[i]);
        }
        int[][] folds = [.. Enumerable.Range(0, k).Select(f => order[(int)((long)f * n / k)..(int)((long)(f + 1) * n / k)].Order().ToArray())];
        var boxes = new Dictionary<string, InputDomain>(StringComparer.Ordinal);
        foreach (RangeScale s in native.Transform.Predictors.OfType<RangeScale>()) boxes[s.Name] = new InputDomain(s.Min, s.Max);

        BassProjection projection = refit.Projection;
        var frames = new BassProjection[k];
        double[] nativeRmse = new double[k], refitRmse = new double[k];
        double nativeSse = 0.0, refitSse = 0.0, sum = 0.0, square = 0.0;
        int count = 0;
        for (int f = 0; f < k; f++)
        {
            int[] held = folds[f];
            int[] train = [.. Enumerable.Range(0, n).Except(held)];
            BassData trainData = Subset(data, train, boxes);
            BassFit foldNative = Bass.Fit(native.Config, trainData);
            ActiveSubspaceAnalysis foldAnalysis = foldNative.ActiveSubspace(trainData, spec.Analysis ?? new ActiveSubspaceSpec());
            frames[f] = BassProjection.Of(foldNative.Transform, foldAnalysis, projection.Dimension, projection.Names);
            BassProjectedFit foldRefit = BassProjectedFit.Of(frames[f], trainData, refit.Fit.Config);
            PredictionCases cases = Cases(data, held);
            double[,] response = Responses(data, held);
            double[,] nativeMean = foldNative.Predict(trainData, cases).Mean(), refitMean = foldRefit.Predict(cases).Mean();
            double ns = 0.0, rs = 0.0;
            for (int r = 0; r < held.Length; r++)
                for (int l = 0; l < response.GetLength(1); l++)
                {
                    double y = response[r, l], en = nativeMean[r, l] - y, er = refitMean[r, l] - y;
                    ns += en * en;
                    rs += er * er;
                    sum += y;
                    square += y * y;
                    count++;
                }
            nativeRmse[f] = Math.Sqrt(ns / response.Length);
            refitRmse[f] = Math.Sqrt(rs / response.Length);
            nativeSse += ns;
            refitSse += rs;
        }
        double sd = Math.Sqrt((square - sum * sum / count) / (count - 1));
        double nativePooled = Math.Sqrt(nativeSse / count), refitPooled = Math.Sqrt(refitSse / count);
        return new CrossFitReport(Array.AsReadOnly(folds), Array.AsReadOnly(frames), nativeRmse, refitRmse, nativePooled, refitPooled,
                                  100.0 * (1.0 - nativePooled / sd), 100.0 * (1.0 - refitPooled / sd));
    }

    // The cases' rows of the data, each continuous predictor taking its box as its declared domain.
    private static BassData Subset(BassData data, int[] rows, IReadOnlyDictionary<string, InputDomain> boxes)
    {
        var predictors = new List<BassPredictor>();
        foreach (BassPredictor predictor in data.Predictors)
            predictors.Add(predictor switch
            {
                ContinuousPredictor c => new ContinuousPredictor(c.Name, Pick(c.Values.Span, rows), boxes[c.Name]),
                CategoricalPredictor z => new CategoricalPredictor(z.Name, Pick(z.Codes.Span, rows), z.Levels),
                _ => throw new ArgumentException($"A predictor of kind {predictor.GetType().Name} is not subset here.", nameof(data)),
            });
        BassResponse response = data.Response switch
        {
            ScalarResponse s => new ScalarResponse(Pick(s.Y.Span, rows)),
            FunctionalResponse fr => new FunctionalResponse(Rows(fr, rows), fr.OutputVariables),
            _ => throw new ArgumentException($"A response of kind {data.Response.GetType().Name} is not subset here.", nameof(data)),
        };
        return new BassData(predictors, response);
    }

    private static PredictionCases Cases(BassData data, int[] rows)
    {
        var columns = new List<PredictionColumn>();
        foreach (BassPredictor predictor in data.Predictors)
            columns.Add(predictor switch
            {
                ContinuousPredictor c => new ContinuousColumn(c.Name, Pick(c.Values.Span, rows)),
                CategoricalPredictor z => new CategoricalColumn(z.Name, [.. Pick(z.Codes.Span, rows).Select(code => z.Levels[code])]),
                _ => throw new ArgumentException($"A predictor of kind {predictor.GetType().Name} is not read here.", nameof(data)),
            });
        return new PredictionCases(columns);
    }

    private static double[,] Responses(BassData data, int[] rows)
    {
        int m = data.Response.Locations;
        var y = new double[rows.Length, m];
        for (int r = 0; r < rows.Length; r++)
            for (int l = 0; l < m; l++)
                y[r, l] = data.Response switch
                {
                    ScalarResponse s => s.Y.Span[rows[r]],
                    FunctionalResponse fr => fr[rows[r], l],
                    _ => throw new ArgumentException("An unknown response."),
                };
        return y;
    }

    private static double[,] Rows(FunctionalResponse response, int[] rows)
    {
        int m = response.Locations;
        var y = new double[rows.Length, m];
        for (int r = 0; r < rows.Length; r++)
            for (int l = 0; l < m; l++) y[r, l] = response[rows[r], l];
        return y;
    }

    private static T[] Pick<T>(ReadOnlySpan<T> values, int[] rows)
    {
        var picked = new T[rows.Length];
        for (int r = 0; r < rows.Length; r++) picked[r] = values[rows[r]];
        return picked;
    }

    private static double[,] Leading(double[,] vectors, int k)
    {
        int p = vectors.GetLength(0);
        var basis = new double[p, k];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < k; j++) basis[i, j] = vectors[i, j];
        return basis;
    }
}
