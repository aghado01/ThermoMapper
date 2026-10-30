using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.ExceptionServices;
using System.Threading.Tasks;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS, Bayesian adaptive spline surfaces (Francom and Sansó 2020): the entry points (the contract shapes). A fit runs
/// <see cref="BassSampling.Chains"/> chains of D04's sweep on the chain ensemble, each started from an exact draw of the
/// prior on F, and keeps every retained draw; or, with <see cref="BassSampling.Smc"/>, <see cref="SmcSpec.Replicates"/>
/// SMC runs from the prior to the posterior (stage 8), and keeps each one's final population with its weights.
/// </summary>
public static class Bass
{
    /// <summary>
    /// Fits <paramref name="config"/> to <paramref name="data"/>: <see cref="BassSampling.BurnIn"/> sweeps per chain, then
    /// <see cref="BassSampling.Samples"/> retained draws <see cref="BassSampling.Thin"/> sweeps apart, or fewer when the
    /// R̂ stop is met. <paramref name="run"/>'s options never change the result. A failure the chains share, such as a
    /// start the prior sampler refuses, is thrown as itself rather than wrapped by the fan-out.
    /// </summary>
    public static BassFit Fit(BassConfig config, BassData data, BassRunOptions? run = null) => Run(config, data, run, null);

    /// <summary>
    /// Continues <paramref name="checkpoint"/>'s run (a fit's <see cref="BassFit.Checkpoint"/>) for
    /// <see cref="BassSampling.Samples"/> more draws per chain, or fewer when the R̂ stop is met, exactly as the
    /// uninterrupted run would have gone on; the fit covers every draw since burn-in. An SMC fit's checkpoint continues
    /// each replicate for up to <see cref="SmcSpec.MaxStages"/> more stages towards b = 1, as one uninterrupted run with
    /// the larger budget would; a run already at b = 1 has nothing left, and gives its fit again. <paramref name="config"/>
    /// and <paramref name="data"/> must be those the checkpoint was taken on: one whose fingerprint differs is refused.
    /// </summary>
    public static BassFit Resume(EnsembleCheckpoint checkpoint, BassConfig config, BassData data, BassRunOptions? run = null)
    {
        ArgumentNullException.ThrowIfNull(checkpoint);
        return Run(config, data, run, checkpoint);
    }

    private static BassFit Run(BassConfig config, BassData data, BassRunOptions? run, EnsembleCheckpoint? resumeFrom)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(data);
        BassRunOptions options = run ?? new BassRunOptions();
        if (options.MaxDegreeOfParallelism is 0 or < -1)
            throw new ArgumentOutOfRangeException(nameof(run), "MaxDegreeOfParallelism is −1, no bound, or at least 1.");
        if (options.ProgressInterval < 1) throw new ArgumentOutOfRangeException(nameof(run), "ProgressInterval must be at least 1.");
        BassModel model = BassModel.Create(config, data);
        BassSampling sampling = config.Sampling;

        string fingerprint = Fingerprint(config, data);
        if (resumeFrom is not null && resumeFrom.Identity.Fingerprint != fingerprint)
            throw new ArgumentException(
                $"The checkpoint was taken on another config or data: its fingerprint is {resumeFrom.Identity.Fingerprint}, theirs {fingerprint}.",
                "checkpoint");
        if (sampling.Smc is not null)
        {
            try
            {
                return BassSmcFit.Run(model, fingerprint, options, resumeFrom);
            }
            catch (AggregateException e) when (Shared(e) is { } failure)
            {
                ExceptionDispatchInfo.Throw(failure);
                throw;
            }
        }
        long planned = sampling.Chains * ((resumeFrom is null ? sampling.BurnIn : 0L) + (long)sampling.Samples * sampling.Thin);
        var ensemble = new BassEnsemble(model, fingerprint, options, planned, (resumeFrom?.DrawsPerChain ?? 0) + sampling.Samples);
        var parallel = new ParallelOptions { MaxDegreeOfParallelism = options.MaxDegreeOfParallelism, CancellationToken = options.Cancellation };
        EnsembleRun summary;
        IReadOnlyList<IChainRun<BassState>> chains;
        try
        {
            (summary, chains) = new ChainEnsemble().Run(ensemble, sampling.Chains, sampling.MasterSeed, sampling.BurnIn, sampling.Samples,
                sampling.RHatTarget, sampling.BatchSize, sampling.MinEssPerChain, options.Checkpoint, resumeFrom, parallel,
                stoppingScope: model.Resolved.StoppingGroups);
        }
        catch (AggregateException e) when (Shared(e) is { } failure)
        {
            ExceptionDispatchInfo.Throw(failure);
            throw;
        }
        ensemble.ReportEnd();
        return ensemble.Assemble(summary, chains);
    }

    // The chains' one kind of failure, which the parallel fan-out wraps: the first, when every inner exception has its type.
    private static Exception? Shared(AggregateException e)
    {
        IReadOnlyList<Exception> inner = e.Flatten().InnerExceptions;
        return inner.Count > 0 && inner.All(x => x.GetType() == inner[0].GetType()) ? inner[0] : null;
    }

    /// <summary>
    /// The run's fingerprint, which a checkpoint must match to resume: FNV-1a (<see cref="RunFingerprint"/>) over the
    /// config's canonical JSON (<see cref="BassConfig.ToJson"/>) and the data's exact bits, predictors in order with
    /// their kinds and names, then the response and its output variables. The config is all values, so the
    /// fingerprint covers every setting, not only the components' types.
    /// </summary>
    public static string Fingerprint(BassConfig config, BassData data)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(data);
        RunFingerprint f = new RunFingerprint().Add("BASS").Add(config.ToJson()).Add(data.Predictors.Count);
        foreach (BassPredictor p in data.Predictors)
            switch (p)
            {
                case ContinuousPredictor c:
                    f.Add("continuous").Add(c.Name).Add(c.Values.Span);
                    break;
                case CategoricalPredictor k:
                    f.Add("categorical").Add(k.Name).Add(k.Levels.Count);
                    foreach (string level in k.Levels) f.Add(level);
                    foreach (int code in k.Codes.Span) f.Add(code);
                    break;
            }
        BassResponse response = data.Response;
        f.Add(response is FunctionalResponse ? "functional" : "scalar").Add(response.Cases).Add(response.Locations).Add(response.Values);
        if (response is FunctionalResponse functional)
        {
            f.Add(functional.OutputVariables.Count);
            foreach (ContinuousPredictor v in functional.OutputVariables) f.Add(v.Name).Add(v.Values.Span);
        }
        return f.ToString();
    }
}
