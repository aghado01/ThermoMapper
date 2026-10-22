using System;
using Maths.Regression.Spline.Bars;
using Maths.Samplers.Tests;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The regression suite's comparative assertions, each with the seeds its test uses, rerun across seeds by
/// <see cref="SeedSweep"/> (step 3.5). An assertion appears as written in its test, as Lower &lt; Upper. Set
/// <c>SEED_SWEEP</c> to the number of seeds per axis.
/// </summary>
public sealed class ComparativeSweep
{
    private readonly ITestOutputHelper _out;
    public ComparativeSweep(ITestOutputHelper output) => _out = output;

    internal static readonly SweptClaim[] Claims =
    [
        new("AdaptiveKernelTests.AdaptiveBurn_InEnsemble", DataSeed: 41, ChainSeeds: [5], (d, c) =>
        {
            var r = AdaptiveKernelTests.AdaptiveAgainstFixedKernel(d, c[0]);
            return [new("peak-location bulk ESS: too-tight fixed kernel < adaptive",
                        r.FixedBad.Peak.LocationConvergence.BulkEss.Value, r.Adaptive.Peak.LocationConvergence.BulkEss.Value)];
        }),
        new("BapsTests.AdaptiveBaps_BeatsGlobal", DataSeed: 123, ChainSeeds: [5], (d, c) =>
        {
            var r = BapsTests.AdaptiveAgainstGlobal(d, c[0]);
            return [new("MSE: adaptive λ < global λ", r.MseAdaptive, r.MseGlobal),
                    new("local λ: at the bump < on the flat ends", r.Central, r.Outer)];
        }),
        new("EbicTests.HigherGamma_SelectsFewerKnots", DataSeed: 17, ChainSeeds: [3], (d, c) =>
        {
            var r = EbicTests.FlatAndStrict(d, c[0]);
            return [new("mean knots: γ = 1 < γ = 0", r.Strict.MeanKnots, r.Flat.MeanKnots)];
        }),
        new("RobustnessTests.Robust_DownweightsOutliers", DataSeed: 44, ChainSeeds: [9], (d, c) =>
        {
            var r = RobustnessTests.AgainstLeastSquares(new StudentTWeights(4.0), d, c[0]);
            return [new("MSE: Student-t weights < least squares", r.MseRobust, r.MseLs)];
        }),
        new("RobustnessTests.ContaminatedNormal_RejectsOutliers", DataSeed: 44, ChainSeeds: [9], (d, c) =>
        {
            var r = RobustnessTests.AgainstLeastSquares(new ContaminatedNormalWeights(0.1, 12.0), d, c[0]);
            return [new("MSE: contaminated-normal weights < least squares", r.MseRobust, r.MseLs)];
        }),
        new("GibbsLossTests.GibbsLoss_RejectsOutliers", DataSeed: 44, ChainSeeds: [9], (d, c) =>
        {
            var r = GibbsLossTests.AgainstLeastSquares(d, c[0]);
            return [new("MSE: Gibbs loss < least squares", r.MseRobust, r.MseLs)];
        }),
        new("HierarchicalPSplineTests.Hierarchical_BorrowsStrength", DataSeed: 2025, ChainSeeds: [7], (d, c) =>
        {
            var r = HierarchicalPSplineTests.HierarchicalAgainstIndependent(d, c[0]);
            return [new("per-replicate MSE: hierarchical < independent REML", r.HierarchicalMse, r.IndependentMse)];
        }),
        new("AnisotropicTensorTests.Anisotropic_BeatsIsotropic", DataSeed: 808, ChainSeeds: [], (d, _) =>
        {
            var r = AnisotropicTensorTests.AnisotropicAgainstIsotropic(d);
            return [new("MSE: anisotropic < isotropic", r.AnisoMse, r.IsoMse),
                    new("smoothing: λy < λx", r.LambdaY, r.LambdaX)];
        }),
        new("ExactChangepointTests.Cauchy_IsRobust", DataSeed: 123, ChainSeeds: [], (d, _) =>
        {
            var r = ExactChangepointTests.GaussianAndCauchy(d);
            return [new("MAP segments: Cauchy < Gaussian", r.CauchyFit.MapSegments, r.GaussianFit.MapSegments)];
        }),
        new("WeightedNormalMarginalTests.Posterior_PrefersTrueModel", DataSeed: 12345, ChainSeeds: [], (d, _) =>
        {
            var r = WeightedNormalMarginalTests.TrueAndBloated(d);
            return [new("log posterior: nine spurious knots < the true two", r.BloatScore, r.TrueScore)];
        }),
        new("PenalizedSplineTests.RemlEvidence_SelectsInteriorSmoothing", DataSeed: 73, ChainSeeds: [], (d, _) =>
        {
            var r = PenalizedSplineTests.RemlAgainstExtremes(d);
            return [new("REML evidence: better grid end < argmax", Math.Max(r.Evidence[0], r.Evidence[^1]), r.Evidence[r.Best]),
                    new("MSE: REML λ* < λ → 0", r.MseStar, r.MseUnder),
                    new("MSE: REML λ* < λ → ∞", r.MseStar, r.MseOver)];
        }),
        new("TensorPSplineTests.TensorPSpline_RecoversSurface", DataSeed: 202, ChainSeeds: [], (d, _) =>
        {
            var r = TensorPSplineTests.RemlAgainstExtremes(d);
            return [new("REML argmax: 0 < steps from the nearer grid end", 0.0, Math.Min(r.Best, r.Lambdas.Length - 1 - r.Best)),
                    new("MSE: REML λ* < λ → 0", r.MseStar, r.MseUnder),
                    new("MSE: REML λ* < λ → ∞", r.MseStar, r.MseOver)];
        }),
        new("PenaltyTests.DifferencePenalty_SmoothsOverRichBasis", DataSeed: 31, ChainSeeds: [], (d, _) =>
        {
            var r = PenaltyTests.PenalizedAgainstOverfit(d);
            return [new("MSE: best of five penalties < unpenalized", r.MseBest, r.MseOver)];
        }),
    ];

    [Fact]
    public void ComparativeAssertions_AcrossSeeds()
    {
        if (SeedSweep.Seeds is not int seeds) return;   // heavy; set SEED_SWEEP to run
        SeedSweep.Run("Maths.Regression.Tests", Claims, seeds, _out);
    }
}
