using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Maths.Rng;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// The kernel-agnostic ensemble shell: fans out C independent shared-nothing replicas (one RNG stream each from a
/// <see cref="SeedTree"/>), burns and then samples each in rounds, reduces every draw to the model's functionals,
/// and pools them into cross-chain R̂ (the consensus map) + ESS — with an R̂ target able to stop sampling early.
/// Parallelism is across replicas, never data. This is the substance the engine family shares: the diagnostics are
/// on reported functionals, so trans-dimensionality / rotational non-identifiability / label symmetry are all
/// handled by the one discipline "reduce to invariant functionals first, R̂ across replicas second". The model
/// (<see cref="IEnsembleModel{TDraw}"/>) supplies the kernel and the reduction; this owns everything else.
/// </summary>
public sealed class ChainEnsemble
{
    /// <summary>
    /// The default effective sample size per chain every monitored functional needs before an R̂ target may stop
    /// sampling: Vehtari et al. (2021) ask for at least 100 per chain, in the bulk and in the tails, before trusting R̂.
    /// </summary>
    public const double MinEssPerChain = 100.0;

    /// <param name="model">Supplies each replica's kernel + start policy and the draw→functional reduction.</param>
    /// <param name="chains">Number of independent replicas (R̂ needs ≥ 2; an R̂ target requires it).</param>
    /// <param name="masterSeed">The single integer the whole fan-out is reproducible from.</param>
    /// <param name="burn">Burn-in steps per replica (no accumulation).</param>
    /// <param name="samples">Post-burn samples per replica to collect (an early R̂ stop may use fewer).</param>
    /// <param name="rHatTarget">When &gt; 0, stop at the end of the first round in which every functional's R̂ meets
    /// this target (<see cref="RHatDiagnostic.MeetsTarget"/>: rank-normalized for the monitored functionals, classic
    /// for the rest), at least one is an estimate, and every monitored functional's bulk and tail ESS reach
    /// <paramref name="minEssPerChain"/> per chain. A frozen or unavailable reading never passes; constant functionals
    /// pass but cannot certify a stop, so an all-constant set runs to <paramref name="samples"/>. It needs at least two
    /// chains and one monitored functional.</param>
    /// <param name="batchSize">Samples per round between R̂ checks (≤ 0 = a single round of <paramref name="samples"/>).</param>
    /// <param name="minEssPerChain">The ESS gate of an R̂ stop, per chain; <see cref="MinEssPerChain"/> by default.</param>
    /// <returns>The pooled <see cref="EnsembleRun"/> plus the per-chain handles, so the client can pool its own sinks.</returns>
    public (EnsembleRun Summary, IReadOnlyList<IChainRun<TDraw>> Chains) Run<TDraw>(
        IEnsembleModel<TDraw> model, int chains, int masterSeed, int burn, int samples,
        double rHatTarget = 0.0, int batchSize = 0, double minEssPerChain = MinEssPerChain)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
        if (rHatTarget > 0.0 && chains < 2)
            throw new ArgumentOutOfRangeException(nameof(chains), "An R̂ stopping target needs at least two chains.");
        if (samples < 1) throw new ArgumentOutOfRangeException(nameof(samples));
        if (!(minEssPerChain >= 0.0)) throw new ArgumentOutOfRangeException(nameof(minEssPerChain));

        int dim = model.FunctionalDim;
        int monitored = model.MonitoredDim;
        if (monitored < 0 || monitored > dim)
            throw new ArgumentOutOfRangeException(nameof(model), "MonitoredDim must be in [0, FunctionalDim].");
        if (rHatTarget > 0.0 && monitored < 1)
            throw new ArgumentException("An R̂ stopping target needs a monitored functional for its ESS gate.", nameof(model));

        // One SeedTree expansion, split into slices: kernel streams take children [0, chains) and
        // the decoupled readout streams (e.g. BARS's per-sample coefficient draws) take
        // [chains, 2·chains) — the ParallelTempering pattern. The former additive master salt
        // aliased across runs: master m+salt's kernel streams were master m's readout streams.
        int[] all = SeedTree.Derive(masterSeed, 2 * chains);
        int[] seeds = all[..chains];
        int[] readoutSeeds = all[chains..];

        var runs = new IChainRun<TDraw>[chains];
        var moments = new RunningMoments[chains][];   // [ci][f] streaming moments of functional f
        var traces = new ChainTrace[chains];          // [ci] the monitored functionals' post-burn draws
        var burnAccepted = new long[chains];
        var burnAttempts = new long[chains];

        Parallel.For(0, chains, ci =>
        {
            IChainRun<TDraw> run = model.StartChain(seeds[ci], readoutSeeds[ci]);
            run.Burn(burn);
            runs[ci] = run;
            burnAccepted[ci] = run.Chain.Accepted;
            burnAttempts[ci] = run.Chain.Attempts;
            moments[ci] = new RunningMoments[dim];
            traces[ci] = new ChainTrace(monitored, samples);
        });

        int effBatch = batchSize <= 0 ? samples : batchSize;
        int collected = 0;
        while (collected < samples)
        {
            int step = Math.Min(effBatch, samples - collected);
            Parallel.For(0, chains, ci =>
            {
                IChain<TDraw> chain = runs[ci].Chain;
                RunningMoments[] cm = moments[ci];
                ChainTrace trace = traces[ci];
                var fn = new double[dim];
                for (int s = 0; s < step; s++)
                {
                    TDraw draw = chain.Step();
                    runs[ci].Accumulate(draw, fn);
                    for (int f = 0; f < dim; f++) cm[f].Add(fn[f]);
                    trace.Add(fn);
                }
            });
            collected += step;

            if (rHatTarget > 0.0 && collected < samples && MayStop(moments, traces, dim, monitored, rHatTarget, minEssPerChain))
                break;
        }

        var mean = new double[dim];
        var variance = new double[dim];
        var rhat = new RHatDiagnostic[dim];
        var chainMoments = new RunningMoments[chains];
        for (int f = 0; f < dim; f++)
        {
            double s = 0.0;
            var pooled = new RunningMoments();
            for (int ci = 0; ci < chains; ci++)
            {
                chainMoments[ci] = moments[ci][f];
                s += moments[ci][f].Mean;
                pooled = RunningMoments.Merge(pooled, moments[ci][f]);
            }
            mean[f] = s / chains;   // every chain holds `collected` draws
            variance[f] = pooled.PopulationVariance;
            rhat[f] = ChainDiagnostics.RHat(chainMoments);
        }

        var diagnostics = new ConvergenceDiagnostic[monitored];
        for (int f = 0; f < monitored; f++) diagnostics[f] = ChainDiagnostics.Diagnose(traces, f);

        // The sampling phase only: the burn-in's counts, adaptation included, are subtracted.
        long accA = 0, accT = 0;
        for (int ci = 0; ci < chains; ci++)
        {
            accA += runs[ci].Chain.Accepted - burnAccepted[ci];
            accT += runs[ci].Chain.Attempts - burnAttempts[ci];
        }
        double acceptance = accT > 0 ? (double)accA / accT : 0.0;

        var summary = new EnsembleRun(mean, variance, rhat, diagnostics, traces, collected, acceptance, seeds);
        return (summary, runs);
    }

    // Whether sampling may stop: every functional's R̂ meets the target (rank-normalized for the monitored ones,
    // classic for the rest), at least one is an estimate, and every monitored functional's bulk and tail ESS reach
    // minEssPerChain per chain. A frozen or unavailable reading blocks the stop; constant functionals pass but cannot
    // certify it on their own.
    private static bool MayStop(RunningMoments[][] moments, ChainTrace[] traces, int dim, int monitored, double target,
                                double minEssPerChain)
    {
        int c = moments.Length;
        bool anyEstimate = false;
        for (int f = 0; f < monitored; f++)
        {
            ConvergenceDiagnostic d = ChainDiagnostics.Diagnose(traces, f);
            if (!d.RHat.MeetsTarget(target) || !d.Ess.MeetsMinimum(minEssPerChain * c)) return false;
            anyEstimate |= d.RHat.State == DiagnosticState.Estimated;
        }
        var chainMoments = new RunningMoments[c];
        for (int f = monitored; f < dim; f++)
        {
            for (int ci = 0; ci < c; ci++) chainMoments[ci] = moments[ci][f];
            RHatDiagnostic r = ChainDiagnostics.RHat(chainMoments);
            if (!r.MeetsTarget(target)) return false;
            anyEstimate |= r.State == DiagnosticState.Estimated;
        }
        return anyEstimate;
    }
}
