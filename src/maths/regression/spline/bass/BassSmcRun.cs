using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Maths.Regression.GPrior;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Smc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's fit by SMC (stage 8; D05): <see cref="SmcSpec.Replicates"/> independent <see cref="SmcRun{TState}"/>s, each
/// from exact draws of the prior on F at b = 0 to the posterior at b = 1 along D05's likelihood path, each next b placed by
/// <see cref="SmcSampler.NextInverseTemperature"/> (the loop of <see cref="SmcSampler.RunTempered"/>, its stage budget a
/// return rather than an exception), mutated by <see cref="BassSmcKernel"/>, then assembled into a <see cref="BassFit"/>.
/// <list type="bullet">
/// <item><b>Streams.</b> The master seed derives each replicate's seed (<see cref="SeedTree"/>), which fans out into its
/// run's seed and its composition's; the composition draws each particle on its own stream, so the draws depend neither
/// on the parallelism nor on the order of composition.</item>
/// <item><b>Parallelism.</b> Replicates run at once within the run's bound, and each one's particles share what is
/// left of it; no result depends on either.</item>
/// <item><b>Draws.</b> Each final particle is composed as a retained draw is (σ², the coefficients and λ; D04), with its
/// normalized weight over the replicates' count, so the fit's weights sum to 1 and each replicate holds 1/R
/// (<see cref="BassFit.DrawWeights"/>). The terms merge in replicate order.</item>
/// <item><b>The budget.</b> A replicate ends at b = 1 or after <see cref="SmcSpec.MaxStages"/> stages. A fit short of
/// b = 1 is flagged (<see cref="BassSmcReadout.ReachedPosterior"/>) and reports no evidence or residual check, and its
/// checkpoint continues it.</item>
/// <item><b>Restart state.</b> Each replicate's <see cref="SmcCheckpoint{TState}"/>, its kernel's account and its prior
/// attempts fill the ensemble checkpoint's engine slot, its moments and trace empty, since a population keeps neither.
/// A resume runs up to <see cref="SmcSpec.MaxStages"/> more stages, exactly as one uninterrupted run with the larger
/// budget would: the stages' streams are derived by position.</item>
/// </list>
/// </summary>
internal static class BassSmcFit
{
    /// <summary>The fit of <paramref name="model"/>, whose config sets <see cref="BassSampling.Smc"/>, or the continuation of
    /// <paramref name="resumeFrom"/>, which must be an SMC fit's checkpoint on the same config and data.</summary>
    public static BassFit Run(BassModel model, string fingerprint, BassRunOptions options, EnsembleCheckpoint? resumeFrom)
    {
        BassSampling sampling = model.Config.Sampling;
        SmcSpec spec = sampling.Smc ?? throw new ArgumentException("The config runs MCMC.", nameof(model));
        int replicates = spec.Replicates;
        BassSmcState[]? restart = null;
        if (resumeFrom is not null)
        {
            if (resumeFrom.Chains.Count != replicates || resumeFrom.ChainSeeds.Length != replicates)
                throw new ArgumentException($"The checkpoint holds {resumeFrom.Chains.Count} runs; the config's SMC runs {replicates}.", "checkpoint");
            restart = new BassSmcState[replicates];
            for (int c = 0; c < replicates; c++)
                restart[c] = resumeFrom.Chains[c].Engine as BassSmcState
                             ?? throw new ArgumentException("Not a BASS SMC fit's checkpoint: a chain's restart state is MCMC's.", "checkpoint");
        }
        int[] seeds = resumeFrom is null ? SeedTree.Derive(sampling.MasterSeed, replicates) : (int[])resumeFrom.ChainSeeds.Clone();

        int bound = options.MaxDegreeOfParallelism == -1 ? Environment.ProcessorCount : options.MaxDegreeOfParallelism;
        int outer = Math.Min(replicates, bound), inner = Math.Max(1, bound / outer);
        long planned = (long)replicates * spec.Particles * spec.MutationSteps * spec.MaxStages;
        var progress = new Progress(options, planned);
        var runs = new Replicate[replicates];
        Parallel.For(0, replicates, new ParallelOptions { MaxDegreeOfParallelism = outer, CancellationToken = options.Cancellation },
                     c => runs[c] = RunOne(model, spec, seeds[c], restart?[c], inner, options, progress));
        progress.End();
        return Assemble(model, fingerprint, spec, seeds, runs, options.Checkpoint);
    }

    // One replicate: started, or resumed with its account, and tempered until b = 1 or its budget; then each final particle
    // composed on its own stream into the replicate's own table.
    private static Replicate RunOne(BassModel model, SmcSpec spec, int seed, BassSmcState? restart, int parallelism, BassRunOptions options,
                                    Progress progress)
    {
        int[] streams = SeedTree.Derive(seed, 2);   // the run's, the composition's
        var target = new BassParticleTarget(model.Target);
        var particles = new BassParticles(model, restart?.PriorAttempts ?? 0);
        var kernel = new BassSmcKernel(model, restart?.Account, options.Cancellation);
        int budget = (restart?.Run.StagesDone ?? 0) + spec.MaxStages;
        var smc = new SmcOptions
        {
            Particles = spec.Particles, Resampling = spec.Resampling, ResampleBelow = spec.ResampleBelow, MutationSteps = spec.MutationSteps,
            MaxStages = budget, MaxDegreeOfParallelism = parallelism, TraceStages = spec.TraceStages,
        };
        SmcRun<BassParticle> run = restart is null
            ? SmcRun<BassParticle>.Start(new SmcStage<BassParticle>(target, new Temperature(0.0, TemperingPath.Likelihood)), particles, kernel, smc, streams[0])
            : SmcRun<BassParticle>.Resume(Retarget(restart.Run, target), particles, kernel, smc);
        while (run.Stage.Temperature.InverseTemperature < 1.0 && run.StagesDone < budget)
        {
            options.Cancellation.ThrowIfCancellationRequested();
            double next = SmcSampler.NextInverseTemperature(run, spec.TargetConditionalEss);
            run.Advance(new SmcStage<BassParticle>(target, new Temperature(next, TemperingPath.Likelihood)));
            progress.Add((long)spec.Particles * spec.MutationSteps);
        }
        SmcResult<BassParticle> result = run.Result();
        BassSmcState? state = options.Checkpoint ? new BassSmcState(run.Checkpoint(), kernel.Account.Copy(), particles.PriorAttempts) : null;

        int[] composition = SeedTree.Derive(streams[1], spec.Particles);
        var terms = new BassTermTable();
        var draws = new BassDraw[spec.Particles];
        for (int i = 0; i < draws.Length; i++)
        {
            BassParticle p = result.Particles[i];
            draws[i] = p.Moves.Compose(p.State, new Xoshiro256PlusPlus(composition[i]), terms);
        }
        return new Replicate(seed, result, terms, draws, kernel.Account.Copy(), particles.PriorAttempts, state);
    }

    // The checkpoint's run on this model's target: the stage's target replaced and each particle evaluated under it, which
    // gives the evaluations the run was saved with, since an evaluation is a function of the state.
    private static SmcCheckpoint<BassParticle> Retarget(SmcCheckpoint<BassParticle> checkpoint, BassParticleTarget target)
    {
        var evaluations = new ITargetEvaluation[checkpoint.Particles.Length];
        for (int i = 0; i < evaluations.Length; i++) evaluations[i] = target.Evaluate(checkpoint.Particles[i]);
        return checkpoint with { Stage = new SmcStage<BassParticle>(target, checkpoint.Stage.Temperature), Evaluations = evaluations };
    }

    private static BassFit Assemble(BassModel model, string fingerprint, SmcSpec spec, int[] seeds, Replicate[] runs, bool checkpoint)
    {
        BassConfig config = model.Config;
        int r = runs.Length, n = spec.Particles;
        var terms = new BassTermTable();
        var stores = new DrawStore<BassDraw>[r];
        var weights = new double[r * n];
        for (int c = 0; c < r; c++)
        {
            BassTermTable own = runs[c].Terms;
            var map = new int[own.Count];
            for (int id = 0; id < map.Length; id++) map[id] = terms.Intern(own[id]);
            stores[c] = new DrawStore<BassDraw>();
            for (int i = 0; i < n; i++)
            {
                stores[c].Add(BassEnsemble.Remap(runs[c].Draws[i], map));
                weights[c * n + i] = runs[c].Result.Weights[i] / r;
            }
        }

        // Each draw's monitored functionals, in the chain ensemble's order: M, τ in the g mode, the log posterior at b = 1,
        // σ², and the fitted value at each monitored point.
        var inputs = new BassPointInputs(model);
        IReadOnlyList<string> names = inputs.FunctionalNames;
        var pointValues = new Dictionary<int, double[]>();
        var functionals = new double[r][][];
        for (int c = 0; c < r; c++)
        {
            functionals[c] = new double[n][];
            for (int i = 0; i < n; i++)
            {
                BassDraw d = stores[c][i];
                var f = new double[names.Count];
                int k = 0;
                f[k++] = d.Structure.TermIds.Length;
                if (model.Target.SamplesTau) f[k++] = d.Tau;
                f[k++] = Temperature.LogPosterior(model.Target.Evaluate(runs[c].Result.Particles[i].State));
                f[k++] = d.Sigma2;
                for (int p = k; p < f.Length; p++) f[p] = d.Coefficients[0];
                int[] ids = d.Structure.TermIds;
                for (int j = 0; j < ids.Length; j++)
                {
                    if (!pointValues.TryGetValue(ids[j], out double[]? values)) pointValues.Add(ids[j], values = inputs.Values(terms[ids[j]]));
                    double a = d.Coefficients[j + 1];
                    for (int p = 0; p < values.Length; p++) f[k + p] += a * values[p];
                }
                functionals[c][i] = f;
            }
        }
        int level = model.Target.SamplesTau ? 2 : 1, first = level + 2;
        var spread = new BassSmcSpread[names.Count];
        for (int q = 0; q < spread.Length; q++) spread[q] = Spread(names[q], functionals, runs, q);
        var replicateLevel = new double[r];
        var replicateFitted = new double[r][];
        for (int c = 0; c < r; c++)
        {
            replicateLevel[c] = spread[level].ReplicateMeans[c];
            replicateFitted[c] = new double[names.Count - first];
            for (int p = 0; p < replicateFitted[c].Length; p++) replicateFitted[c][p] = spread[first + p].ReplicateMeans[c];
        }

        var accounts = new BassSmcReplicate[r];
        var account = new BassSmcAccount(model.Prior);
        double massAtCap = 0.0;
        for (int c = 0; c < r; c++)
        {
            SmcResult<BassParticle> result = runs[c].Result;
            double atCap = 0.0;
            for (int i = 0; i < n; i++)
                if (stores[c][i].Structure.TermIds.Length == model.MaxBasis) atCap += result.Weights[i];
            massAtCap += atCap / r;
            SmcStageReport last = result.Stages[^1];
            accounts[c] = new BassSmcReplicate(runs[c].Seed, result.Stages, last.Temperature.InverseTemperature, result.LogNormalizerRatio,
                                               result.LogNormalizerError, SmcWeights.EffectiveSampleSize(result.Weights), last.Families,
                                               runs[c].PriorAttempts, runs[c].Account, atCap);
            account.MergeFrom(runs[c].Account);
        }
        var readout = new BassSmcReadout(accounts, spread);
        bool posterior = readout.ReachedPosterior;

        BassSampling sampling = config.Sampling;
        ConvergenceDiagnostic none = new(RHatDiagnostic.Unavailable, RHatDiagnostic.Unavailable, RHatDiagnostic.Unavailable,
                                         EssDiagnostic.Unavailable, EssDiagnostic.Unavailable);
        var fitted = new BassMonitoredValue[inputs.Points.Count];
        for (int p = 0; p < fitted.Length; p++) fitted[p] = new BassMonitoredValue(inputs.Points[p], none);
        var verdict = new ConvergenceVerdict(ConvergenceOutcome.Undetermined, ChainEnsemble.ReferenceRHatTarget, sampling.MinEssPerChain,
            [.. new[] { BassDiagnostics.PredictiveGroup, BassDiagnostics.LevelGroup, BassDiagnostics.StructuralGroup }
                .Select(g => new GroupVerdict(g, ConvergenceOutcome.Undetermined, false, []))]);
        var diagnostics = new BassDiagnostics(none, model.Target.SamplesTau ? none : null, none, none, fitted, account.Moves, account.AcceptanceRate,
                                              account.Kinds, account.Refusals, account.EmptyBirths, account.Support, massAtCap, n, false, verdict,
                                              BassModeSeparation.Of(replicateLevel, replicateFitted), [], null, readout);

        bool evidence = posterior && sampling.ReportEvidence && model.Target.Sigma2 is not null && model.Target.Intercept != GPriorCentre.Flat;
        var identity = new RunIdentity(sampling.MasterSeed, r, fingerprint);
        EnsembleCheckpoint? restart = checkpoint
            ? new EnsembleCheckpoint(identity, (int[])seeds.Clone(), n,
                                     [.. runs.Select(x => new ChainCheckpoint(x.State!, [], new ChainTrace(0), 0, 0))])
            : null;
        return new BassFit(config, model.Resolved, model.Transform, identity, terms, stores, weights, diagnostics,
                           evidence ? BassEvidence.OfSmc(model, [.. runs.Select(x => x.Result)]) : null,
                           posterior ? BassResidualDependence.Of(model, terms, stores, weights: weights) : null, restart);
    }

    // Quantity q's between-run reading: each replicate's weighted mean, the pooled posterior's mean and sd (each replicate
    // weighing 1/R), and the SD of the replicates' means over that sd.
    private static BassSmcSpread Spread(string name, double[][][] functionals, Replicate[] runs, int q)
    {
        int r = runs.Length;
        var means = new double[r];
        double pooled = 0.0;
        for (int c = 0; c < r; c++)
        {
            double[] w = runs[c].Result.Weights;
            double m = 0.0, total = 0.0;
            for (int i = 0; i < w.Length; i++)
            {
                m += w[i] * functionals[c][i][q];
                total += w[i];
            }
            means[c] = m / total;
            pooled += means[c] / r;
        }
        double variance = 0.0;
        for (int c = 0; c < r; c++)
        {
            double[] w = runs[c].Result.Weights;
            double total = 0.0, squares = 0.0;
            for (int i = 0; i < w.Length; i++)
            {
                double e = functionals[c][i][q] - pooled;
                squares += w[i] * e * e;
                total += w[i];
            }
            variance += squares / total / r;
        }
        double sd = Math.Sqrt(variance), between = double.NaN;
        if (r >= 2)
        {
            double s = 0.0;
            foreach (double m in means) s += (m - pooled) * (m - pooled);
            between = Math.Sqrt(s / (r - 1));
        }
        double spread = r < 2 ? double.NaN : sd > 0.0 ? between / sd : between > 0.0 ? double.PositiveInfinity : 0.0;
        return new BassSmcSpread(name, means, pooled, sd, spread);
    }

    // A replicate's run, its composed draws in its own table, its kernel's account, its prior attempts and, when asked, its
    // restart state.
    private sealed record Replicate(int Seed, SmcResult<BassParticle> Result, BassTermTable Terms, BassDraw[] Draws, BassSmcAccount Account,
                                    long PriorAttempts, BassSmcState? State);

    // Sweeps over every replicate, reported each time they pass a multiple of the interval, and once at the end.
    private sealed class Progress(BassRunOptions options, long planned)
    {
        private long _sweeps;

        public void Add(long sweeps)
        {
            if (options.Progress is not { } progress) return;
            long before = Interlocked.Add(ref _sweeps, sweeps) - sweeps, after = before + sweeps;
            if (after / options.ProgressInterval > before / options.ProgressInterval) progress.Report(new BassProgress(after, planned));
        }

        public void End()
        {
            if (options.Progress is not { } progress) return;
            long done = Interlocked.Read(ref _sweeps);
            if (done == 0 || done % options.ProgressInterval != 0) progress.Report(new BassProgress(done, planned));
        }
    }
}

/// <summary>An SMC replicate's restart state: its run's checkpoint, its kernel's account and its prior draws' attempts.</summary>
internal sealed record BassSmcState(SmcCheckpoint<BassParticle> Run, BassSmcAccount Account, long PriorAttempts) : IChainState;
