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
    /// <param name="model">Supplies each replica's kernel + start policy and the draw→functional reduction.</param>
    /// <param name="chains">Number of independent replicas (R̂/ESS need ≥ 2 to be defined).</param>
    /// <param name="masterSeed">The single integer the whole fan-out is reproducible from.</param>
    /// <param name="burn">Burn-in steps per replica (no accumulation).</param>
    /// <param name="samples">Post-burn samples per replica to collect (an early R̂ stop may use fewer).</param>
    /// <param name="rHatTarget">When &gt; 0, stop once every functional's R̂ is ≤ this (checked each round; a
    /// degenerate-variance functional, R̂ NaN/∞, is treated as converged so it cannot block the stop).</param>
    /// <param name="batchSize">Samples per round between R̂ checks (≤ 0 = a single round of <paramref name="samples"/>).</param>
    /// <returns>The pooled <see cref="EnsembleRun"/> plus the per-chain handles, so the client can pool its own sinks.</returns>
    public (EnsembleRun Summary, IReadOnlyList<IChainRun<TDraw>> Chains) Run<TDraw>(
        IEnsembleModel<TDraw> model, int chains, int masterSeed, int burn, int samples,
        double rHatTarget = 0.0, int batchSize = 0)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
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
        var sums = new double[chains][];          // [ci][f] running Σ of functional f
        var sumSqs = new double[chains][];        // [ci][f] running Σ² of functional f
        var seqs = new List<double>[chains][];    // [ci][f<essDim] the post-burn sequence (for ESS)

        Parallel.For(0, chains, ci =>
        {
            IChainRun<TDraw> run = model.StartChain(seeds[ci], readoutSeeds[ci]);
            run.Burn(burn);
            runs[ci] = run;
            sums[ci] = new double[dim];
            sumSqs[ci] = new double[dim];
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
                double[] cs = sums[ci], css = sumSqs[ci];
                List<double>[] cseq = seqs[ci];
                var fn = new double[dim];
                for (int s = 0; s < step; s++)
                {
                    TDraw draw = chain.Step();
                    runs[ci].Accumulate(draw, fn);
                    for (int f = 0; f < dim; f++) { double v = fn[f]; cs[f] += v; css[f] += v * v; }
                    for (int f = 0; f < essDim; f++) cseq[f].Add(fn[f]);
                }
            });
            collected += step;

            if (rHatTarget > 0.0 && collected >= effBatch && MaxClampedRHat(sums, sumSqs, dim, collected) <= rHatTarget)
                break;
        }

        double denom = chains * (double)collected;
        var mean = new double[dim];
        var rhat = new double[dim];
        var chainSums = new double[chains];
        var chainSqs = new double[chains];
        for (int f = 0; f < dim; f++)
        {
            double s = 0.0;
            for (int ci = 0; ci < chains; ci++) { double v = sums[ci][f]; s += v; chainSums[ci] = v; chainSqs[ci] = sumSqs[ci][f]; }
            mean[f] = s / denom;
            rhat[f] = ChainDiagnostics.RHat(chainSums, chainSqs, collected);
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

        var summary = new EnsembleRun(mean, rhat, ess, collected, acceptance, seeds);
        return (summary, runs);
    }

    // Max R̂ over all functionals, treating a non-finite (degenerate-variance) functional as converged (→ 0) so it
    // cannot stall the early stop — the same clamp the client applies to the per-functional consensus map.
    private static double MaxClampedRHat(double[][] sums, double[][] sumSqs, int dim, int n)
    {
        int c = sums.Length;
        var cs = new double[c];
        var css = new double[c];
        double max = 0.0;
        for (int f = 0; f < dim; f++)
        {
            for (int ci = 0; ci < c; ci++) { cs[ci] = sums[ci][f]; css[ci] = sumSqs[ci][f]; }
            double r = ChainDiagnostics.RHat(cs, css, n);
            double rv = double.IsNaN(r) || double.IsInfinity(r) ? 0.0 : r;
            if (rv > max) max = rv;
        }
        return max;
    }
}
