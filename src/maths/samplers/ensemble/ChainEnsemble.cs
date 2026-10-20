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
    /// Draws per chain before an R̂ target may stop sampling. Vehtari et al. (2021) ask for an effective sample size
    /// of at least 100 per chain before trusting R̂, and the draw count bounds the effective sample size from above.
    /// </summary>
    public const int MinDrawsForStop = 100;

    /// <param name="model">Supplies each replica's kernel + start policy and the draw→functional reduction.</param>
    /// <param name="chains">Number of independent replicas (R̂ needs ≥ 2; an R̂ target requires it).</param>
    /// <param name="masterSeed">The single integer the whole fan-out is reproducible from.</param>
    /// <param name="burn">Burn-in steps per replica (no accumulation).</param>
    /// <param name="samples">Post-burn samples per replica to collect (an early R̂ stop may use fewer).</param>
    /// <param name="rHatTarget">When &gt; 0, stop once every functional meets this target
    /// (<see cref="RHatDiagnostic.MeetsTarget"/>) and at least one is an estimate, checked each round from
    /// <see cref="MinDrawsForStop"/> draws per chain on. A frozen or unavailable R̂ never meets it; an all-constant
    /// set runs to <paramref name="samples"/>.</param>
    /// <param name="batchSize">Samples per round between R̂ checks (≤ 0 = a single round of <paramref name="samples"/>).</param>
    /// <returns>The pooled <see cref="EnsembleRun"/> plus the per-chain handles, so the client can pool its own sinks.</returns>
    public (EnsembleRun Summary, IReadOnlyList<IChainRun<TDraw>> Chains) Run<TDraw>(
        IEnsembleModel<TDraw> model, int chains, int masterSeed, int burn, int samples,
        double rHatTarget = 0.0, int batchSize = 0)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
        if (rHatTarget > 0.0 && chains < 2)
            throw new ArgumentOutOfRangeException(nameof(chains), "An R̂ stopping target needs at least two chains.");
        if (samples < 1) throw new ArgumentOutOfRangeException(nameof(samples));

        int dim = model.FunctionalDim;
        int essDim = model.EssDim;
        if (essDim < 0 || essDim > dim) throw new ArgumentOutOfRangeException(nameof(model), "EssDim must be in [0, FunctionalDim].");

        // One SeedTree expansion, split into slices: kernel streams take children [0, chains) and
        // the decoupled readout streams (e.g. BARS's per-sample coefficient draws) take
        // [chains, 2·chains) — the ParallelTempering pattern. The former additive master salt
        // aliased across runs: master m+salt's kernel streams were master m's readout streams.
        int[] all = SeedTree.Derive(masterSeed, 2 * chains);
        int[] seeds = all[..chains];
        int[] readoutSeeds = all[chains..];

        var runs = new IChainRun<TDraw>[chains];
        var moments = new RunningMoments[chains][];   // [ci][f] streaming moments of functional f
        var seqs = new List<double>[chains][];        // [ci][f<essDim] the post-burn sequence (for ESS)

        Parallel.For(0, chains, ci =>
        {
            IChainRun<TDraw> run = model.StartChain(seeds[ci], readoutSeeds[ci]);
            run.Burn(burn);
            runs[ci] = run;
            moments[ci] = new RunningMoments[dim];
            var s = new List<double>[essDim];
            for (int f = 0; f < essDim; f++) s[f] = new List<double>(samples);
            seqs[ci] = s;
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
                List<double>[] cseq = seqs[ci];
                var fn = new double[dim];
                for (int s = 0; s < step; s++)
                {
                    TDraw draw = chain.Step();
                    runs[ci].Accumulate(draw, fn);
                    for (int f = 0; f < dim; f++) cm[f].Add(fn[f]);
                    for (int f = 0; f < essDim; f++) cseq[f].Add(fn[f]);
                }
            });
            collected += step;

            if (rHatTarget > 0.0 && collected >= Math.Max(effBatch, MinDrawsForStop)
                && AllMeetTarget(moments, dim, rHatTarget))
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

        var ess = new double[essDim];
        for (int f = 0; f < essDim; f++)
        {
            var chainSeq = new double[chains][];
            for (int ci = 0; ci < chains; ci++) chainSeq[ci] = seqs[ci][f].ToArray();
            ess[f] = ChainDiagnostics.Ess(chainSeq);
        }

        long accA = 0, accT = 0;
        for (int ci = 0; ci < chains; ci++) { accA += runs[ci].Chain.Accepted; accT += runs[ci].Chain.Attempts; }
        double acceptance = accT > 0 ? (double)accA / accT : 0.0;

        var summary = new EnsembleRun(mean, variance, rhat, ess, collected, acceptance, seeds);
        return (summary, runs);
    }

    // Whether every functional's R̂ meets the target and at least one of them is an actual estimate. A frozen or
    // unavailable functional blocks the stop; constant functionals don't block it but can't certify it on their own.
    private static bool AllMeetTarget(RunningMoments[][] moments, int dim, double target)
    {
        int c = moments.Length;
        var chainMoments = new RunningMoments[c];
        bool anyEstimate = false;
        for (int f = 0; f < dim; f++)
        {
            for (int ci = 0; ci < c; ci++) chainMoments[ci] = moments[ci][f];
            RHatDiagnostic r = ChainDiagnostics.RHat(chainMoments);
            if (!r.MeetsTarget(target)) return false;
            anyEstimate |= r.State == RHatState.Estimated;
        }
        return anyEstimate;
    }
}
