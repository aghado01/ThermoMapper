using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Maths.Rng;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// The kernel-agnostic ensemble shell: fans out C independent shared-nothing replicas (one RNG stream each from a
/// <see cref="SeedTree"/>), burns and then samples each in rounds, reduces every draw to the model's functionals,
/// and pools them into cross-chain R̂ (the consensus map) + ESS — with an R̂ target able to stop sampling early, and a
/// convergence verdict (<see cref="ConvergenceVerdict"/>) on every run.
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

    /// <summary>
    /// The R̂ target a verdict reads when the run sets none: 1.01, Vehtari et al.'s (2021) threshold for rank-normalized
    /// R̂ (D10).
    /// </summary>
    public const double ReferenceRHatTarget = 1.01;

    /// <param name="model">Supplies each replica's kernel + start policy and the draw→functional reduction.</param>
    /// <param name="chains">Number of independent replicas (R̂ needs ≥ 2; an R̂ target requires it).</param>
    /// <param name="masterSeed">The single integer the whole fan-out is reproducible from.</param>
    /// <param name="burn">Burn-in steps per replica (no accumulation).</param>
    /// <param name="samples">Post-burn samples per replica to collect (an early R̂ stop may use fewer).</param>
    /// <param name="rHatTarget">When &gt; 0, stop at the end of the first round whose verdict reads
    /// <see cref="ConvergenceOutcome.Converged"/> against this target: in every gating group, every functional's R̂ meets
    /// it (<see cref="RHatDiagnostic.MeetsTarget"/>: rank-normalized for the monitored functionals, classic for the rest)
    /// and every monitored functional's bulk and tail ESS reach <paramref name="minEssPerChain"/> per chain, and some
    /// gating group holds an estimate. A frozen or unavailable reading never passes; constant functionals pass but cannot
    /// certify a stop, so an all-constant scope runs to <paramref name="samples"/>. It needs at least two chains and a
    /// monitored functional in a gating group. The end-of-run verdict reads the same rule against this target, or
    /// <see cref="ReferenceRHatTarget"/> when it is 0.</param>
    /// <param name="batchSize">Samples per round between R̂ checks (≤ 0 = a single round of <paramref name="samples"/>).</param>
    /// <param name="minEssPerChain">The ESS gate of an R̂ stop, per chain; <see cref="MinEssPerChain"/> by default.</param>
    /// <param name="checkpoint">When true, <see cref="EnsembleRun.Checkpoint"/> holds the run's restart state at its end;
    /// every chain run must then be an <see cref="IResumableChainRun{TDraw}"/>.</param>
    /// <param name="resumeFrom">When given, continue that checkpoint's run instead of starting one: its chains are
    /// rebuilt by <see cref="IEnsembleModel{TDraw}.RestoreChain"/>, <paramref name="burn"/> is not run again, and
    /// <paramref name="samples"/> more draws per chain are collected. The checkpoint's identity must match
    /// <paramref name="masterSeed"/>, <paramref name="chains"/> and the model's <see cref="IEnsembleModel{TDraw}.Fingerprint"/>.
    /// A run resumed from a checkpoint continues exactly as the uninterrupted run would have.</param>
    /// <param name="parallel">Bounds the fan-out over chains (<see cref="ParallelOptions.MaxDegreeOfParallelism"/>) and
    /// carries its cancellation; null runs every chain at once. Chains share nothing, so neither changes a result.</param>
    /// <param name="stoppingScope">The names of the model's groups (<see cref="IEnsembleModel{TDraw}.Groups"/>) that gate
    /// an R̂ stop, and whose conjunction is the verdict's outcome; null for every group. The verdict reports every group
    /// whatever the scope. An unknown name, or an empty scope, is refused.</param>
    /// <returns>The pooled <see cref="EnsembleRun"/> plus the per-chain handles, so the client can pool its own sinks.</returns>
    public (EnsembleRun Summary, IReadOnlyList<IChainRun<TDraw>> Chains) Run<TDraw>(
        IEnsembleModel<TDraw> model, int chains, int masterSeed, int burn, int samples,
        double rHatTarget = 0.0, int batchSize = 0, double minEssPerChain = MinEssPerChain,
        bool checkpoint = false, EnsembleCheckpoint? resumeFrom = null, ParallelOptions? parallel = null,
        IReadOnlyCollection<string>? stoppingScope = null)
    {
        ArgumentNullException.ThrowIfNull(model);
        parallel ??= new ParallelOptions();
        if (chains < 1) throw new ArgumentOutOfRangeException(nameof(chains));
        if (rHatTarget > 0.0 && chains < 2)
            throw new ArgumentOutOfRangeException(nameof(chains), "An R̂ stopping target needs at least two chains.");
        if (samples < 1) throw new ArgumentOutOfRangeException(nameof(samples));
        if (!(minEssPerChain >= 0.0)) throw new ArgumentOutOfRangeException(nameof(minEssPerChain));

        int dim = model.FunctionalDim;
        int monitored = model.MonitoredDim;
        if (monitored < 0 || monitored > dim)
            throw new ArgumentOutOfRangeException(nameof(model), "MonitoredDim must be in [0, FunctionalDim].");
        (IReadOnlyList<string> names, IReadOnlyList<FunctionalGroup> groups) = Declarations(model, dim);
        bool[] gates = Gates(groups, stoppingScope);
        if (rHatTarget > 0.0 && !GatesMonitored(groups, gates, monitored))
            throw new ArgumentException("An R̂ stopping target needs a monitored functional in a gating group for its ESS gate.",
                                        monitored < 1 ? nameof(model) : nameof(stoppingScope));
        double verdictTarget = rHatTarget > 0.0 ? rHatTarget : ReferenceRHatTarget;

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
        int collected = 0;

        if (resumeFrom is null)
        {
            Parallel.For(0, chains, parallel, ci =>
            {
                IChainRun<TDraw> run = model.StartChain(seeds[ci], readoutSeeds[ci]);
                run.Burn(burn);
                runs[ci] = run;
                burnAccepted[ci] = run.Chain.Accepted;
                burnAttempts[ci] = run.Chain.Attempts;
                moments[ci] = new RunningMoments[dim];
                traces[ci] = new ChainTrace(monitored, samples);
            });
        }
        else
        {
            var identity = new RunIdentity(masterSeed, chains, model.Fingerprint);
            if (resumeFrom.Identity != identity || resumeFrom.Chains.Count != chains)
                throw new ArgumentException(
                    $"The checkpoint belongs to another run ({resumeFrom.Identity}, not {identity}).", nameof(resumeFrom));
            seeds = (int[])resumeFrom.ChainSeeds.Clone();
            collected = resumeFrom.DrawsPerChain;
            for (int ci = 0; ci < chains; ci++)   // in order, so a model that cannot restore says so directly
            {
                ChainCheckpoint cp = resumeFrom.Chains[ci];
                runs[ci] = model.RestoreChain(cp.Engine);
                burnAccepted[ci] = cp.BurnAccepted;
                burnAttempts[ci] = cp.BurnAttempts;
                moments[ci] = (RunningMoments[])cp.Moments.Clone();
                traces[ci] = cp.Trace.Copy();
            }
        }

        int total = collected + samples;
        int effBatch = batchSize <= 0 ? samples : batchSize;
        while (collected < total)
        {
            int step = Math.Min(effBatch, total - collected);
            Parallel.For(0, chains, parallel, ci =>
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

            if (rHatTarget > 0.0 && collected < total)
            {
                (RHatDiagnostic[] classic, ConvergenceDiagnostic[] ranked) = Readings(moments, traces, dim, monitored);
                ConvergenceVerdict round = ConvergenceVerdict.Of(classic, ranked, chains, names, groups, gates, rHatTarget, minEssPerChain);
                if (round.Outcome == ConvergenceOutcome.Converged) break;
            }
        }

        var mean = new double[dim];
        var variance = new double[dim];
        for (int f = 0; f < dim; f++)
        {
            double s = 0.0;
            var pooled = new RunningMoments();
            for (int ci = 0; ci < chains; ci++)
            {
                s += moments[ci][f].Mean;
                pooled = RunningMoments.Merge(pooled, moments[ci][f]);
            }
            mean[f] = s / chains;   // every chain holds `collected` draws
            variance[f] = pooled.PopulationVariance;
        }

        (RHatDiagnostic[] rhat, ConvergenceDiagnostic[] diagnostics) = Readings(moments, traces, dim, monitored);
        ConvergenceVerdict verdict = ConvergenceVerdict.Of(rhat, diagnostics, chains, names, groups, gates, verdictTarget, minEssPerChain);

        // The sampling phase only: the burn-in's counts, adaptation included, are subtracted.
        long accA = 0, accT = 0;
        for (int ci = 0; ci < chains; ci++)
        {
            accA += runs[ci].Chain.Accepted - burnAccepted[ci];
            accT += runs[ci].Chain.Attempts - burnAttempts[ci];
        }
        double acceptance = accT > 0 ? (double)accA / accT : 0.0;

        EnsembleCheckpoint? restart = null;
        if (checkpoint)
        {
            var parts = new ChainCheckpoint[chains];
            for (int ci = 0; ci < chains; ci++)
            {
                if (runs[ci] is not IResumableChainRun<TDraw> resumable)
                    throw new InvalidOperationException($"{runs[ci].GetType().Name} cannot save its restart state.");
                parts[ci] = new ChainCheckpoint(resumable.SaveState(), (RunningMoments[])moments[ci].Clone(), traces[ci].Copy(),
                                                burnAccepted[ci], burnAttempts[ci]);
            }
            restart = new EnsembleCheckpoint(new RunIdentity(masterSeed, chains, model.Fingerprint), (int[])seeds.Clone(),
                                             collected, parts);
        }

        var summary = new EnsembleRun(mean, variance, rhat, diagnostics, traces, collected, acceptance, seeds, verdict, restart);
        return (summary, runs);
    }

    // Every functional's classic R̂ from the chains' streaming moments, and the monitored ones' rank-normalized diagnostics
    // from their traces: what both the stop between rounds and the end-of-run verdict read.
    private static (RHatDiagnostic[] Classic, ConvergenceDiagnostic[] Monitored) Readings(
        RunningMoments[][] moments, ChainTrace[] traces, int dim, int monitored)
    {
        int c = moments.Length;
        var classic = new RHatDiagnostic[dim];
        var chainMoments = new RunningMoments[c];
        for (int f = 0; f < dim; f++)
        {
            for (int ci = 0; ci < c; ci++) chainMoments[ci] = moments[ci][f];
            classic[f] = ChainDiagnostics.RHat(chainMoments);
        }
        var ranked = new ConvergenceDiagnostic[monitored];
        for (int f = 0; f < monitored; f++) ranked[f] = ChainDiagnostics.Diagnose(traces, f);
        return (classic, ranked);
    }

    // The model's functional names and groups, checked: one name per functional, and groups with distinct names that
    // between them hold every functional exactly once.
    private static (IReadOnlyList<string> Names, IReadOnlyList<FunctionalGroup> Groups) Declarations<TDraw>(
        IEnsembleModel<TDraw> model, int dim)
    {
        IReadOnlyList<string> names = model.FunctionalNames ?? throw new ArgumentException("FunctionalNames is missing.", nameof(model));
        if (names.Count != dim)
            throw new ArgumentException($"FunctionalNames holds {names.Count} names for {dim} functionals.", nameof(model));
        for (int f = 0; f < dim; f++)
            if (names[f] is null) throw new ArgumentException($"FunctionalNames has no name for functional {f}.", nameof(model));

        IReadOnlyList<FunctionalGroup> groups = model.Groups ?? throw new ArgumentException("Groups is missing.", nameof(model));
        var seen = new HashSet<string>(StringComparer.Ordinal);
        var owner = new int[dim];
        Array.Fill(owner, -1);
        for (int g = 0; g < groups.Count; g++)
        {
            FunctionalGroup group = groups[g] ?? throw new ArgumentException($"Group {g} is missing.", nameof(model));
            if (string.IsNullOrWhiteSpace(group.Name)) throw new ArgumentException($"Group {g} has no name.", nameof(model));
            if (!seen.Add(group.Name)) throw new ArgumentException($"Two groups are named \"{group.Name}\".", nameof(model));
            if (group.Functionals is not { Count: > 0 })
                throw new ArgumentException($"Group \"{group.Name}\" holds no functional.", nameof(model));
            foreach (int f in group.Functionals)
            {
                if (f < 0 || f >= dim)
                    throw new ArgumentException($"Group \"{group.Name}\" holds functional {f}, outside [0, {dim}).", nameof(model));
                if (owner[f] >= 0)
                    throw new ArgumentException(
                        $"Functional {f} ({names[f]}) belongs to two groups, \"{groups[owner[f]].Name}\" and \"{group.Name}\".", nameof(model));
                owner[f] = g;
            }
        }
        for (int f = 0; f < dim; f++)
            if (owner[f] < 0) throw new ArgumentException($"Functional {f} ({names[f]}) belongs to no group.", nameof(model));
        return (names, groups);
    }

    // Which groups gate the stop: every group when the scope is null, otherwise the named ones, each of which must exist.
    private static bool[] Gates(IReadOnlyList<FunctionalGroup> groups, IReadOnlyCollection<string>? scope)
    {
        var gates = new bool[groups.Count];
        if (scope is null)
        {
            Array.Fill(gates, true);
            return gates;
        }
        if (scope.Count == 0)
            throw new ArgumentException("A stopping scope names at least one group; null takes every group.", "stoppingScope");
        foreach (string name in scope)
        {
            int g = 0;
            while (g < groups.Count && groups[g].Name != name) g++;
            if (g == groups.Count)
            {
                var known = new List<string>(groups.Count);
                foreach (FunctionalGroup group in groups) known.Add(group.Name);
                throw new ArgumentException(
                    $"The stopping scope names \"{name}\", which is not one of the model's groups ({string.Join(", ", known)}).", "stoppingScope");
            }
            gates[g] = true;
        }
        return gates;
    }

    // Whether a gating group holds a monitored functional, which an R̂ stop's ESS gate needs.
    private static bool GatesMonitored(IReadOnlyList<FunctionalGroup> groups, bool[] gates, int monitored)
    {
        for (int g = 0; g < groups.Count; g++)
            if (gates[g])
                foreach (int f in groups[g].Functionals)
                    if (f < monitored) return true;
        return false;
    }
}
