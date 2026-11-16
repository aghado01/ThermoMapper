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

    /// <summary>The held-out cases' weights, one per case, at the data's meaning (slice 8b, D2): required with <see cref="TestCases"/>
    /// when the data carry case weights, refused otherwise. Replication weights count each case wᵢ times in the RMSE, the sd and the
    /// density, as duplicated rows would; precision weights enter the density as σ²/wᵢ alone, and each case counts once in the RMSE.</summary>
    public IReadOnlyList<double>? TestWeights { get; init; }

    /// <summary>Cross-fitting's folds K: 0, the default, for none; otherwise at least 2. It costs K native fits and K refits.</summary>
    public int Folds { get; init; }

    /// <summary>The seed of the folds' assignment.</summary>
    public int FoldSeed { get; init; } = 1;

    /// <summary>Cross-fitting's unit (slice 8b, D1): null, the default, assigns cases to folds, so a grouping's levels spread across
    /// them; the name of one of the data's groupings, which random effects read, assigns its levels whole, so each fold scores groups no training fold saw, the
    /// leave-groups-out question. A held-out case in a level no training fold took is predicted at the effect's prior mean, zero, and its
    /// density integrates the effect.</summary>
    public string? FoldGrouping { get; init; }

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
public sealed record ProjectionFitReport(double BasisFunctions, PredictiveCheck Loo, BassEvidence? Evidence, double? HeldOutRmse, double? HeldOutExplained)
{
    /// <summary>The held-out log predictive density Σᵢ cᵢ·ln Σ_s w̃_s p(yᵢ | xᵢ, θ_s) (slice 8b, C), through the candidate seam
    /// (<see cref="BassFit.CandidateLogDensity"/>): every noise law, the case weights and random effects as the seam reads them, a new
    /// group's effect integrated from its prior; cᵢ = wᵢ under replication weights, else 1. Null without held-out cases, and for a
    /// response that is not scalar: the seam serves one observation of a scalar response (C1).</summary>
    public double? HeldOutLogDensity { get; init; }
}

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
                                   double NativePooled, double RefitPooled, double NativeExplained, double RefitExplained)
{
    /// <summary>The grouping whose levels were assigned whole to folds (slice 8b, D1); null when cases were.</summary>
    public string? FoldGrouping { get; init; }

    /// <summary>Each fold's native held-out log predictive density (slice 8b, C), as <see cref="ProjectionFitReport.HeldOutLogDensity"/>;
    /// null for a response that is not scalar.</summary>
    public double[]? NativeLogDensity { get; init; }

    /// <summary>Each fold's refit held-out log predictive density; null for a response that is not scalar.</summary>
    public double[]? RefitLogDensity { get; init; }

    /// <summary>The native densities summed over the folds, every case held out once.</summary>
    public double? NativeLogDensityTotal { get; init; }

    /// <summary>The refit densities summed over the folds.</summary>
    public double? RefitLogDensityTotal { get; init; }
}

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
        if (spec.FoldGrouping is not null && spec.Folds == 0) throw new ArgumentException("A fold grouping assigns cross-fitting's folds; ask for folds.", nameof(spec));
        double[]? testWeights = null;
        if (spec.TestCases is not null)
        {
            if (data.Weights is not null && spec.TestWeights is null)
                throw new ArgumentException($"The data carry case weights ({data.Weights.Meaning}), so each held-out case needs its own.", nameof(spec));
            if (data.Weights is null && spec.TestWeights is not null)
                throw new ArgumentException("The data carry no case weights, so a held-out case's would state nothing: give none.", nameof(spec));
            if (spec.TestWeights is not null) testWeights = CaseWeightsOf(spec.TestWeights, spec.TestResponse!.GetLength(0));
        }
        else if (spec.TestWeights is not null) throw new ArgumentException("Held-out weights come with held-out cases.", nameof(spec));

        PredictionCases? projected = spec.TestCases is null ? null : refit.Projection.ProjectCases(spec.TestCases);
        BassPredictOptions? options = Unseen(data);
        ProjectionFitReport nativeSide = Side(native, data, spec.TestCases, spec.TestCases is null ? null : native.Predict(data, spec.TestCases, options).Mean(),
                                              spec.TestResponse, testWeights);
        ProjectionFitReport refitSide = Side(refit.Fit, refit.Data, projected, spec.TestCases is null ? null : refit.Predict(spec.TestCases, options).Mean(),
                                             spec.TestResponse, testWeights);

        ProjectedActiveSubspace read = refit.ActiveSubspace(new ActiveSubspaceSpec { Points = spec.Points, Seed = spec.Seed });
        double[] toFrame = [.. Enumerable.Range(1, d).Select(k => SubspacePosterior.PrincipalAngles(Leading(read.InputSpace.MeanCVectors, k), Leading(projection.Frame, k)).Max())];
        double[][] spread = [.. analysis.Posterior.Draws.Select(draw =>
            SubspacePosterior.PrincipalAngles(Leading(draw.Vectors, Math.Min(d, draw.Vectors.GetLength(1))), projection.Frame))];
        CrossFitReport? cross = spec.Folds == 0 ? null : CrossFit(refit, native, data, spec);
        return new BassProjectionComparison(nativeSide, refitSide, Array.AsReadOnly(toFrame), spread, cross);
    }

    // A projected basis fit beside its native basis fit (slice 8c): each side's basis functions per joint draw summed over the components,
    // PSIS-LOO on the coefficients, and the reconstruction's held-out error; no evidence and no density (C1); the refit's subspace against
    // the frame and the native draws' spread about it. Cross-fitting would decompose each fold's training outputs again (E2), which needs
    // the decomposition's spec, and held-out weights have no meaning for a basis fit's components; both are refused.
    public static BassProjectionComparison Of(BassProjectedBasisFit refit, BassBasisFit native, ActiveSubspaceAnalysis analysis, ProjectionComparisonSpec spec)
    {
        ArgumentNullException.ThrowIfNull(native);
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
        if (spec.Folds != 0 || spec.FoldGrouping is not null)
            throw new ArgumentException("Cross-fitting a basis fit would decompose each fold's training outputs again, which needs the decomposition's own spec; "
                                        + "it is not built.", nameof(spec));
        if (spec.TestWeights is not null) throw new ArgumentException("A basis fit's components carry no case weights, so a held-out case's would state nothing.", nameof(spec));
        ProjectionFitReport nativeSide = BasisSide(native, spec.TestCases is null ? null : native.Predict(spec.TestCases).Mean(), spec.TestResponse);
        ProjectionFitReport refitSide = BasisSide(refit.Fit, spec.TestCases is null ? null : refit.Predict(spec.TestCases).Mean(), spec.TestResponse);
        ProjectedActiveSubspace read = refit.ActiveSubspace(new ActiveSubspaceSpec { Points = spec.Points, Seed = spec.Seed });
        double[] toFrame = [.. Enumerable.Range(1, d).Select(k => SubspacePosterior.PrincipalAngles(Leading(read.InputSpace.MeanCVectors, k), Leading(projection.Frame, k)).Max())];
        double[][] spread = [.. analysis.Posterior.Draws.Select(draw =>
            SubspacePosterior.PrincipalAngles(Leading(draw.Vectors, Math.Min(d, draw.Vectors.GetLength(1))), projection.Frame))];
        return new BassProjectionComparison(nativeSide, refitSide, Array.AsReadOnly(toFrame), spread, null);
    }

    // A basis fit's side: Σ_k of each joint draw's component draw's basis functions, averaged over the joint draws; PSIS-LOO on the
    // coefficients (BassBasisFit.PointwiseLogDensity's default space); the reconstruction's held-out RMSE.
    private static ProjectionFitReport BasisSide(BassBasisFit fit, double[,]? mean, double[,]? response)
    {
        BassDraw[][] pooled = [.. fit.Components.Select(c => c.Draws.SelectMany(chain => chain).ToArray())];
        double terms = 0.0;
        for (int b = 0; b < fit.JointDraws.Count; b++)
            for (int k = 0; k < pooled.Length; k++)
                terms += pooled[k][fit.JointDraws.Draw(k, b)].Structure.TermCount(fit.Components[k].Terms);
        terms /= fit.JointDraws.Count;
        PredictiveCheck loo = PredictiveChecks.PsisLoo(fit.PointwiseLogDensity());
        double? rmse = null, explained = null;
        if (mean is not null) (rmse, explained) = Score(mean, response!, null);
        return new ProjectionFitReport(terms, loo, null, rmse, explained);
    }

    // One fit's side; cases are the held-out cases in this fit's inputs, for the density.
    private static ProjectionFitReport Side(BassFit fit, BassData data, PredictionCases? cases, double[,]? mean, double[,]? response, double[]? weights)
    {
        BassDraw[] pooled = [.. fit.Draws.SelectMany(c => c)];
        IReadOnlyList<double>? w = fit.DrawWeights;
        int exponent = WeightScale.Exponent(w);
        double Weight(int s) => w is null ? 1.0 : WeightScale.Scaled(w[s], exponent);
        double total = w is null ? pooled.Length : Enumerable.Range(0, pooled.Length).Sum(Weight);
        double terms = Enumerable.Range(0, pooled.Length).Sum(s => Weight(s) * pooled[s].Structure.TermCount(fit.Terms)) / total;
        PredictiveCheck loo = PredictiveChecks.PsisLoo(fit.PointwiseLogDensity(data));
        double? rmse = null, explained = null, density = null;
        if (mean is not null)
        {
            (rmse, explained) = Score(mean, response!, Counts(data, weights));
            density = HeldOutLogDensity(fit, data, cases!, response!, weights);
        }
        return new ProjectionFitReport(terms, loo, fit.Evidence, rmse, explained) { HeldOutLogDensity = density };
    }

    // The RMSE and the percentage of the responses' sd, each case's row counted c_r times (1 unless replication weights count it).
    private static (double Rmse, double Explained) Score(double[,] mean, double[,] response, double[]? counts)
    {
        if (mean.GetLength(0) != response.GetLength(0) || mean.GetLength(1) != response.GetLength(1))
            throw new ArgumentException($"The held-out responses are {mean.GetLength(0)} × {mean.GetLength(1)}, one row per case and one column per location.");
        double sse = 0.0, sum = 0.0, count = 0.0;
        for (int r = 0; r < response.GetLength(0); r++)
        {
            double c = counts?[r] ?? 1.0;
            for (int l = 0; l < response.GetLength(1); l++)
            {
                sum += c * response[r, l];
                count += c;
            }
        }
        double average = sum / count, ss = 0.0;
        for (int r = 0; r < response.GetLength(0); r++)
        {
            double weight = counts?[r] ?? 1.0;
            for (int l = 0; l < response.GetLength(1); l++)
            {
                double e = mean[r, l] - response[r, l], c = response[r, l] - average;
                sse += weight * e * e;
                ss += weight * c * c;
            }
        }
        double rmse = Math.Sqrt(sse / count), sd = Math.Sqrt(ss / (count - 1));
        return (rmse, 100.0 * (1.0 - rmse / sd));
    }

    // Each held-out case's count (slice 8b, D2): its weight under replication weights, as duplicated rows count; null, every case
    // once, otherwise.
    private static double[]? Counts(BassData data, double[]? weights)
        => data.Weights?.Meaning == WeightMeaning.Replication ? weights : null;

    // Σᵢ cᵢ·ln Σ_s w̃_s p(yᵢ | xᵢ, θ_s) over held-out cases (slice 8b, C) through the candidate seam, each case's response its one
    // candidate. A weighted fit's seam reads each case's weight at its meaning: a precision as σ²/wᵢ; under replication one observation,
    // weight 1, counted wᵢ times here. A noise read by quadrature and a new group under a noise take the seam's quadrature; a new group's
    // effect is integrated from its prior. Null for a response that is not scalar (C1).
    private static double? HeldOutLogDensity(BassFit fit, BassData data, PredictionCases cases, double[,] response, double[]? weights)
    {
        if (data.Response is not ScalarResponse) return null;
        BassModel model = fit.ModelFor(data);
        int n = response.GetLength(0);
        double[]? counts = Counts(data, weights);
        bool reads = model.Observation.Law is null && !model.Working.IsUnweighted;
        IReadOnlyList<double>? caseWeights = !reads ? null : counts is not null ? Enumerable.Repeat(1.0, n).ToArray() : weights;
        UnseenLevel unseen = Unseen(data)?.UnseenLevels ?? UnseenLevel.Refuse;
        int[]? strata = model.RandomEffects is null ? null : fit.GroupCells(data, cases, unseen).Candidates;
        double[][] grids = [.. Enumerable.Range(0, n).Select(r => new[] { response[r, 0] })];
        CandidateLogDensity density = fit.CandidateLogDensity(data, cases, grids, strata, unseen, quadrature: true, caseWeights: caseWeights);
        int draws = density.Draws;
        var logWeight = new double[draws];
        for (int s = 0; s < draws; s++) logWeight[s] = density.Weights is { } w ? Math.Log(w[s]) : 0.0;
        double normalizer = LogSumExp(logWeight, null, 0), total = 0.0;
        for (int t = 0; t < n; t++)
            total += (counts?[t] ?? 1.0) * (LogSumExp(logWeight, density.Values, density.Offset(t)) - normalizer);
        return total;
    }

    // ln Σ_s exp(a_s + v_s), v the column of values or 0.
    private static double LogSumExp(double[] a, double[,]? values, int column)
    {
        double top = double.NegativeInfinity;
        for (int s = 0; s < a.Length; s++) top = Math.Max(top, a[s] + (values?[s, column] ?? 0.0));
        if (double.IsNegativeInfinity(top)) return top;
        double sum = 0.0;
        for (int s = 0; s < a.Length; s++) sum += Math.Exp(a[s] + (values?[s, column] ?? 0.0) - top);
        return top + Math.Log(sum);
    }

    // Data with groupings predict a level no training case took as a new group, the effect at its prior mean, zero (slice 8b, D1);
    // this also places a predictor's unseen level outside every subset.
    private static BassPredictOptions? Unseen(BassData data)
        => data.Groupings.Count > 0 ? new BassPredictOptions { UnseenLevels = UnseenLevel.OutsideEverySubset } : null;

    private static double[] CaseWeightsOf(IReadOnlyList<double> weights, int cases)
    {
        if (weights.Count != cases) throw new ArgumentException($"One held-out weight per case: {weights.Count} for {cases} cases.", "spec");
        var w = new double[cases];
        for (int r = 0; r < cases; r++)
        {
            w[r] = weights[r];
            if (!(w[r] > 0.0 && double.IsFinite(w[r]))) throw new ArgumentException($"Held-out case {r}'s weight is {w[r]:R}; a weight is positive and finite.", "spec");
        }
        return w;
    }

    // K folds: by default the cases shuffled by the seed and cut into K parts of sizes differing by at most one; under a fold grouping its
    // observed levels shuffled and cut so, each level's cases in one fold (slice 8b, D1). Each fold's frame from a native fit on the other
    // folds, its refit on them, and both scored on it. Every fold's data take the full data's boxes as their domains, and carry the
    // groupings' columns and the weights' rows; the held-out cases carry their group columns into prediction (slice 8b).
    private static CrossFitReport CrossFit(BassProjectedFit refit, BassFit native, BassData data, ProjectionComparisonSpec spec)
    {
        int n = data.Response.Cases, k = spec.Folds;
        int[][] folds = spec.FoldGrouping is null ? CaseFolds(n, k, spec.FoldSeed) : GroupFolds(data, spec.FoldGrouping, k, spec.FoldSeed);
        var boxes = new Dictionary<string, InputDomain>(StringComparer.Ordinal);
        foreach (RangeScale s in native.Transform.Predictors.OfType<RangeScale>()) boxes[s.Name] = new InputDomain(s.Min, s.Max);
        BassPredictOptions? options = Unseen(data);
        bool scalar = data.Response is ScalarResponse;

        BassProjection projection = refit.Projection;
        var frames = new BassProjection[k];
        double[] nativeRmse = new double[k], refitRmse = new double[k];
        double[]? nativeDensity = scalar ? new double[k] : null, refitDensity = scalar ? new double[k] : null;
        double nativeSse = 0.0, refitSse = 0.0, sum = 0.0, square = 0.0, count = 0.0;
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
            double[]? weights = data.Weights is null ? null : Pick(data.Weights.Values.Span, held), counts = Counts(data, weights);
            double[,] nativeMean = foldNative.Predict(trainData, cases, options).Mean(), refitMean = foldRefit.Predict(cases, options).Mean();
            double ns = 0.0, rs = 0.0, foldCount = 0.0;
            for (int r = 0; r < held.Length; r++)
            {
                double c = counts?[r] ?? 1.0;
                for (int l = 0; l < response.GetLength(1); l++)
                {
                    double y = response[r, l], en = nativeMean[r, l] - y, er = refitMean[r, l] - y;
                    ns += c * en * en;
                    rs += c * er * er;
                    sum += c * y;
                    square += c * y * y;
                    foldCount += c;
                }
            }
            count += foldCount;
            nativeRmse[f] = Math.Sqrt(ns / foldCount);
            refitRmse[f] = Math.Sqrt(rs / foldCount);
            nativeSse += ns;
            refitSse += rs;
            if (scalar)
            {
                nativeDensity![f] = HeldOutLogDensity(foldNative, trainData, cases, response, weights)!.Value;
                refitDensity![f] = HeldOutLogDensity(foldRefit.Fit, foldRefit.Data, foldRefit.Projection.ProjectCases(cases), response, weights)!.Value;
            }
        }
        double sd = Math.Sqrt((square - sum * sum / count) / (count - 1));
        double nativePooled = Math.Sqrt(nativeSse / count), refitPooled = Math.Sqrt(refitSse / count);
        return new CrossFitReport(Array.AsReadOnly(folds), Array.AsReadOnly(frames), nativeRmse, refitRmse, nativePooled, refitPooled,
                                  100.0 * (1.0 - nativePooled / sd), 100.0 * (1.0 - refitPooled / sd))
        {
            FoldGrouping = spec.FoldGrouping,
            NativeLogDensity = nativeDensity,
            RefitLogDensity = refitDensity,
            NativeLogDensityTotal = nativeDensity?.Sum(),
            RefitLogDensityTotal = refitDensity?.Sum(),
        };
    }

    // The cases shuffled by the seed and cut into K parts of sizes differing by at most one (slice 7c).
    private static int[][] CaseFolds(int n, int k, int seed)
    {
        if (k > n) throw new ArgumentOutOfRangeException("spec", $"There are {n} cases, so at most {n} folds.");
        int[] order = Shuffled(n, seed);
        return [.. Enumerable.Range(0, k).Select(f => order[(int)((long)f * n / k)..(int)((long)(f + 1) * n / k)].Order().ToArray())];
    }

    // A grouping's observed levels, in code order, shuffled by the seed and cut into K parts as the cases are; each fold the cases whose
    // level is in its part, ascending (slice 8b, D1).
    private static int[][] GroupFolds(BassData data, string name, int k, int seed)
    {
        CategoricalPredictor grouping = data.Groupings.FirstOrDefault(g => g.Name == name)
            ?? throw new ArgumentException($"The data hold no grouping named '{name}'.", "spec");
        ReadOnlySpan<int> codes = grouping.Codes.Span;
        int[] observed = [.. grouping.Codes.ToArray().Distinct().Order()];
        if (k > observed.Length) throw new ArgumentOutOfRangeException("spec", $"The grouping '{name}' has {observed.Length} observed levels, so at most {observed.Length} folds.");
        int[] order = Shuffled(observed.Length, seed);
        var foldOf = new Dictionary<int, int>();
        for (int f = 0; f < k; f++)
            for (int at = (int)((long)f * observed.Length / k); at < (int)((long)(f + 1) * observed.Length / k); at++) foldOf[observed[order[at]]] = f;
        var folds = Enumerable.Range(0, k).Select(_ => new List<int>()).ToArray();
        for (int i = 0; i < codes.Length; i++) folds[foldOf[codes[i]]].Add(i);
        return [.. folds.Select(fold => fold.ToArray())];
    }

    // 0 … n − 1 in Fisher–Yates order from the seed.
    private static int[] Shuffled(int n, int seed)
    {
        int[] order = [.. Enumerable.Range(0, n)];
        var rng = new Xoshiro256PlusPlus(seed);
        for (int i = n - 1; i > 0; i--)
        {
            int j = rng.NextInt(i + 1);
            (order[i], order[j]) = (order[j], order[i]);
        }
        return order;
    }

    // The cases' rows of the data, each continuous predictor taking its box as its declared domain, with the groupings' columns and the
    // case weights' rows (slice 8b); a grouping keeps its levels, so a level no training row takes reads as a new group.
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
        CategoricalPredictor[] groupings = [.. data.Groupings.Select(g => new CategoricalPredictor(g.Name, Pick(g.Codes.Span, rows), g.Levels))];
        CaseWeights? weights = data.Weights is { } w ? new CaseWeights(Pick(w.Values.Span, rows), w.Meaning, w.PriorGram) : null;
        return new BassData(predictors, response, groupings, weights);
    }

    // The cases' rows as prediction cases: the predictors' columns, then the groupings' (slice 8b).
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
        foreach (CategoricalPredictor g in data.Groupings)
            columns.Add(new CategoricalColumn(g.Name, [.. Pick(g.Codes.Span, rows).Select(code => g.Levels[code])]));
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
