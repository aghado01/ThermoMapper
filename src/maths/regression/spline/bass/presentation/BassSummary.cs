using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bass.Presentation;

/// <summary>
/// The data's shape as a fit holds it (R's <c>print</c>; P12, P13): the cases and output locations, each predictor's
/// training range or levels (<see cref="RangeScale"/>, <see cref="TrainingLevels"/>), and a direct functional fit's output
/// variables with their ranges; none for a scalar or basis fit.
/// </summary>
public sealed record BassDataShape(int Cases, int Locations, IReadOnlyList<PredictorScale> Predictors, IReadOnlyList<RangeScale> OutputVariables)
{
    /// <summary>The categorical predictors, R's <c>p.cat</c>.</summary>
    public int CategoricalPredictors => Predictors.Count(p => p is TrainingLevels);

    internal static BassDataShape Of(BassInputTransform transform, int cases, int locations)
        => new(cases, locations, transform.Predictors, transform.OutputVariables);
}

/// <summary>A run's lengths: its chains, the burn-in sweeps, the draws asked for per chain and the thinning, and the draws
/// each chain kept, which the R̂ stop can shorten (<see cref="BassDiagnostics.SamplesUsed"/>). An SMC fit reads no sweep
/// count: <see cref="Chains"/> is then its replicates, <see cref="DrawsPerChain"/> each replicate's population,
/// <see cref="BurnIn"/>, <see cref="Samples"/> and <see cref="Thin"/> null, and <see cref="Smc"/> its specification.</summary>
public sealed record BassRunLengths(int Chains, int? BurnIn, int? Samples, int? Thin, IReadOnlyList<int> DrawsPerChain, SmcSpec? Smc = null)
{
    /// <summary>The draws kept over every chain.</summary>
    public int Draws => DrawsPerChain.Sum();
}

/// <summary>
/// A fit's metadata (R's <c>print.bass</c>; P12): the data's shape, the config and what it resolved to on the data with its
/// auto-picks, the run's identity and lengths, and its account: convergence, the moves' acceptance and refusals, the
/// support tally and the tempering readout (<see cref="BassDiagnostics"/>), with the evidence when the run reports one. The
/// records are the fit's own, so the metadata round-trips what the fit holds.
/// </summary>
public sealed record BassMetadata(
    BassDataShape Data, BassConfig Config, BassResolved Resolved, RunIdentity Identity, BassRunLengths Lengths,
    BassDiagnostics Diagnostics, BassEvidence? Evidence)
{
    internal static BassMetadata Of(BassFit fit)
    {
        ArgumentNullException.ThrowIfNull(fit);
        BassSampling sampling = fit.Config.Sampling;
        bool sweeps = sampling.Smc is null;
        var lengths = new BassRunLengths(fit.Identity.Chains, sweeps ? sampling.BurnIn : null, sweeps ? sampling.Samples : null,
                                         sweeps ? sampling.Thin : null, Array.AsReadOnly(fit.Draws.Select(c => c.Count).ToArray()),
                                         sampling.Smc);
        return new BassMetadata(BassDataShape.Of(fit.Transform, fit.Resolved.Cases, fit.Resolved.Locations), fit.Config, fit.Resolved,
                                fit.Identity, lengths, fit.Diagnostics, fit.Evidence);
    }
}

/// <summary>
/// The posterior summaries of a fit's draws pooled over its chains under its draw weights (R's <c>summary.bass</c>; P14):
/// M's, whose <see cref="PosteriorSummary.Min"/> and <see cref="PosteriorSummary.Max"/> are R's range of the basis count;
/// σ²'s, whose mean is R's posterior mean error variance; τ's when it is sampled (null under fixed g); and λ's.
/// </summary>
public sealed record BassPosteriorSummary(PosteriorSummary Count, PosteriorSummary Sigma2, PosteriorSummary? Tau, PosteriorSummary Lambda)
{
    internal static BassPosteriorSummary Of(BassFit fit, IReadOnlyList<double> probabilities)
    {
        var draws = fit.Draws.SelectMany(c => c).ToArray();
        if (fit.DrawWeights is { } w && w.Length != draws.Length)
            throw new ArgumentException($"The fit has {draws.Length} draws and {w.Length} weights.", nameof(fit));
        PosteriorSummary Summarized(Func<BassDraw, double> value) => PosteriorSummary.Over(Array.ConvertAll(draws, d => value(d)), fit.DrawWeights, probabilities);
        return new BassPosteriorSummary(Summarized(d => d.Structure.TermIds.Length), Summarized(d => d.Sigma2),
                                        fit.Config.Coefficients is SampledG ? Summarized(d => d.Tau) : null, Summarized(d => d.Lambda));
    }
}

/// <summary>A fit's summary (R's <c>summary.bass</c>; P14): its metadata and its posterior summaries.</summary>
public sealed record BassSummary(BassMetadata Metadata, BassPosteriorSummary Posterior)
{
    internal static BassSummary Of(BassFit fit, IReadOnlyList<double>? probabilities)
        => new(BassMetadata.Of(fit), BassPosteriorSummary.Of(fit, PosteriorSummary.Checked(probabilities)));
}

/// <summary>A decomposition's truncation residuals, as scalars (D07): their kind and RMS, each set's RMS, and the
/// basis-stability ratio (<see cref="TruncationResiduals"/>).</summary>
public sealed record BassTruncationSummary(TruncationResidual Kind, double Rms, double? InSampleRms, double? LeaveOneOutRms, double? StabilityRatio);

/// <summary>
/// A basis fit's metadata (R's <c>print.bassBasis</c>; P13): the data's shape, n cases at m locations over the shared
/// predictors; the decomposition's components r, its PCA account (null for a supplied basis) and its truncation residuals;
/// the master seed and the joint draws, with their chain lengths (null after resampling); and each component fit's own
/// metadata, in component order.
/// </summary>
public sealed record BassBasisMetadata(
    BassDataShape Data, int Components, BassPcaAccount? Pca, BassTruncationSummary Truncation, int MasterSeed,
    int JointDraws, IReadOnlyList<int>? JointChainLengths, bool Resampled, IReadOnlyList<BassMetadata> ComponentFits)
{
    internal static BassBasisMetadata Of(BassBasisFit fit)
    {
        ArgumentNullException.ThrowIfNull(fit);
        BassDecomposition decomposition = fit.Decomposition;
        TruncationResiduals t = decomposition.Truncation;
        return new BassBasisMetadata(
            BassDataShape.Of(fit.Components[0].Transform, decomposition.Cases, decomposition.Locations),
            decomposition.Rank, decomposition.Pca, new BassTruncationSummary(t.Kind, t.Rms, t.InSampleRms, t.LeaveOneOutRms, t.StabilityRatio),
            fit.MasterSeed, fit.JointDraws.Count, fit.JointDraws.ChainLengths, fit.JointDraws.Resampled,
            Array.AsReadOnly(fit.Components.Select(BassMetadata.Of).ToArray()));
    }
}

/// <summary>A basis fit's summary (R's <c>summary.bassBasis</c>; P15): its metadata and each component's posterior
/// summaries over its own draws, in component order.</summary>
public sealed record BassBasisSummary(BassBasisMetadata Metadata, IReadOnlyList<BassPosteriorSummary> Components)
{
    internal static BassBasisSummary Of(BassBasisFit fit, IReadOnlyList<double>? probabilities)
    {
        IReadOnlyList<double> p = PosteriorSummary.Checked(probabilities);
        return new(BassBasisMetadata.Of(fit), Array.AsReadOnly(fit.Components.Select(c => BassPosteriorSummary.Of(c, p)).ToArray()));
    }
}
