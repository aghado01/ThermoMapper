using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace Maths.Samplers.FlatHistogram;

/// <summary>The weights at one iteration of the learning, for its trace.</summary>
public sealed record LearningSnapshot(long Iteration, MulticanonicalWeights Weights);

/// <summary>
/// A multicanonical run's restart state: its identity, the learning's state, diagnostics and trace, and either the
/// chains' restart states (while learning, or at its end before production) or the ensemble shell's checkpoint (once
/// production has begun), whose chains carry their production records.
/// </summary>
public sealed record MulticanonicalCheckpoint(
    RunIdentity Identity, SamcLearningState Learning, FlatHistogramDiagnostics LearningDiagnostics,
    IReadOnlyList<LearningSnapshot> WeightTrace, IReadOnlyList<IChainState>? Chains, EnsembleCheckpoint? Production);

/// <summary>
/// A multicanonical run's outcome, so far.
/// </summary>
/// <param name="Weights">The weights: frozen once learning is done, the learning's current ones otherwise.</param>
/// <param name="Learning">The learning's state.</param>
/// <param name="LearningDiagnostics">The learning's visits: occupancy, outside counts, round trips.</param>
/// <param name="WeightTrace">The learning's weights, at its start and every hundredth of its length, and at its end.</param>
/// <param name="Production">The ensemble shell's summary of production, D10's verdict on π_w's functionals among it;
/// null before production.</param>
/// <param name="Draws">The production draws, for the reweighting; null before production.</param>
/// <param name="ProductionDiagnostics">Production's visits; null before production.</param>
/// <param name="Chains">The live chains.</param>
/// <param name="Checkpoint">The run's restart state, when asked for.</param>
/// <param name="Complete">Whether the run reached its configured lengths.</param>
public sealed record MulticanonicalResult<TState>(
    MulticanonicalWeights Weights, SamcLearningState Learning, FlatHistogramDiagnostics LearningDiagnostics,
    IReadOnlyList<LearningSnapshot> WeightTrace, EnsembleRun? Production, MulticanonicalDraws? Draws,
    FlatHistogramDiagnostics? ProductionDiagnostics, IReadOnlyList<IMulticanonicalChain<TState>> Chains,
    MulticanonicalCheckpoint? Checkpoint, bool Complete);

/// <summary>
/// The flat-histogram run: SAMC learns the multicanonical weights during burn-in, the chains sharing one weight
/// function (<see cref="SamcLearning"/>), and production then runs the chains independently on the frozen weights,
/// over the ensemble shell (<see cref="ChainEnsemble"/>), so D10's verdict reads π_w's functionals and the shell's
/// checkpoint carries production. A learning iteration steps every chain once, in parallel at fixed weights, then
/// applies one update from their visits in chain order and reweighs each chain's cached evaluation; the result does not
/// depend on the parallelism bound. The model supplies the chains on the wrapped target
/// (<see cref="IMulticanonicalModel{TState}"/>); the seeds are the shell's, from <see cref="SeedTree"/>.
/// </summary>
public static class MulticanonicalRun
{
    /// <summary>
    /// Run <paramref name="spec"/> on <paramref name="target"/>, or continue <paramref name="resumeFrom"/>'s run, taking
    /// at most <paramref name="iterations"/> iterations (a learning update or a production draw per chain each) before
    /// returning; null runs to the spec's lengths.
    /// </summary>
    /// <param name="spec">The configuration.</param>
    /// <param name="target">The target π_w wraps.</param>
    /// <param name="model">The chains and their functionals.</param>
    /// <param name="masterSeed">The seed the chains' streams derive from.</param>
    /// <param name="resumeFrom">A checkpoint of this run (same spec, seed and model fingerprint) to continue.</param>
    /// <param name="iterations">At most this many iterations in this call, ≥ 1; null for all that remain.</param>
    /// <param name="checkpoint">Whether to return the run's restart state; every chain must then save its state.</param>
    /// <param name="parallel">Bounds the fan-out over chains and carries cancellation; neither changes a result.</param>
    public static MulticanonicalResult<TState> Run<TState>(
        FlatHistogramSpec spec, ISamplingTarget<TState> target, IMulticanonicalModel<TState> model, int masterSeed,
        MulticanonicalCheckpoint? resumeFrom = null, long? iterations = null,
        bool checkpoint = false, ParallelOptions? parallel = null)
    {
        ArgumentNullException.ThrowIfNull(spec);
        ArgumentNullException.ThrowIfNull(target);
        ArgumentNullException.ThrowIfNull(model);
        spec.Validate();
        if (iterations is long budget0 && budget0 < 1) throw new ArgumentOutOfRangeException(nameof(iterations), "A call takes at least one iteration.");
        parallel ??= new ParallelOptions();
        int chains = spec.Chains;
        var identity = new RunIdentity(masterSeed, chains, Fingerprint(spec, model));

        // The start: resolved from the spec and the pilot, or restored.
        SamcLearning learning;
        HistogramAccumulator learnHistogram;
        List<LearningSnapshot> trace;
        IMulticanonicalChain<TState>[] live = new IMulticanonicalChain<TState>[chains];
        MulticanonicalTarget<TState> wrapped;
        EnsembleCheckpoint? production = null;
        int[] seeds = SeedTree.Derive(masterSeed, 2 * chains);
        if (resumeFrom is null)
        {
            (MulticanonicalWeights start, double t0, long gainStart) = Resolve(spec);
            learning = new SamcLearning(start, t0, gainStart);
            wrapped = new MulticanonicalTarget<TState>(target, start);
            learnHistogram = new HistogramAccumulator(start, chains);
            trace = [new LearningSnapshot(0, start)];
            Parallel.For(0, chains, parallel, c => live[c] = Started(model, wrapped, seeds[c], seeds[chains + c]));
        }
        else
        {
            if (resumeFrom.Identity != identity)
                throw new ArgumentException($"The checkpoint belongs to another run ({resumeFrom.Identity}, not {identity}).", nameof(resumeFrom));
            learning = new SamcLearning(resumeFrom.Learning);
            wrapped = new MulticanonicalTarget<TState>(target, learning.Weights);
            learnHistogram = new HistogramAccumulator(learning.Weights, resumeFrom.LearningDiagnostics);
            trace = [.. resumeFrom.WeightTrace];
            production = resumeFrom.Production;
            if (production is null)
            {
                IReadOnlyList<IChainState> states = resumeFrom.Chains
                    ?? throw new ArgumentException("The checkpoint holds neither chains nor production.", nameof(resumeFrom));
                if (states.Count != chains) throw new ArgumentException("The checkpoint holds another chain count.", nameof(resumeFrom));
                for (int c = 0; c < chains; c++) live[c] = Restored(model, wrapped, states[c]);   // in order, so a refusal is direct
            }
        }

        long remaining = iterations ?? long.MaxValue;
        long traceEvery = Math.Max(1, spec.LearningIterations / 100);

        // Learning: every chain steps at the current weights, then one update from their visits, in chain order.
        if (production is null && learning.Iterations < spec.LearningIterations)
        {
            var levels = new double[chains];
            long steps = Math.Min(spec.LearningIterations - learning.Iterations, remaining);
            for (long i = 0; i < steps; i++)
            {
                Parallel.For(0, chains, parallel, c => live[c].Step());
                for (int c = 0; c < chains; c++)
                {
                    levels[c] = live[c].CurrentEvaluation.UntemperedLogLikelihood;
                    learnHistogram.Add(c, levels[c]);
                }
                wrapped.Install(learning.Update(levels));
                for (int c = 0; c < chains; c++) live[c].Reweigh();
                if (learning.Iterations % traceEvery == 0 || learning.Iterations == spec.LearningIterations)
                    trace.Add(new LearningSnapshot(learning.Iterations, learning.Weights));
            }
            remaining -= steps;
        }

        bool learned = learning.Iterations >= spec.LearningIterations;
        int produced = production?.DrawsPerChain ?? 0;
        if (produced >= spec.ProductionDraws) throw new InvalidOperationException("The checkpoint's run is complete.");
        if (!learned || remaining <= 0)
        {
            MulticanonicalCheckpoint? held = checkpoint
                ? new MulticanonicalCheckpoint(identity, learning.Save(), learnHistogram.Snapshot(), [.. trace], Saved(live), null)
                : null;
            return new MulticanonicalResult<TState>(learning.Weights, learning.Save(), learnHistogram.Snapshot(), [.. trace],
                null, null, null, live, held, Complete: false);
        }

        // Production: independent chains on the frozen weights, over the shell, burn-in already done.
        int samples = (int)Math.Min(spec.ProductionDraws - produced, remaining);
        var shellModel = new ProductionModel<TState>(model, wrapped, identity.Fingerprint, production is null ? live : null, seeds[..chains]);
        (EnsembleRun summary, IReadOnlyList<IChainRun<TState>> runs) = new ChainEnsemble().Run(
            shellModel, chains, masterSeed, burn: 0, samples, checkpoint: checkpoint, resumeFrom: production, parallel: parallel);

        var levelsOut = new double[chains][];
        var functionals = new double[chains][];
        var diagnostics = new FlatHistogramDiagnostics[chains];
        var handles = new IMulticanonicalChain<TState>[chains];
        for (int c = 0; c < chains; c++)
        {
            var run = (ProductionRun<TState>)runs[c];
            levelsOut[c] = [.. run.Levels];
            functionals[c] = [.. run.Functionals];
            diagnostics[c] = run.Histogram.Snapshot();
            handles[c] = run.Chain;
        }
        var draws = new MulticanonicalDraws(learning.Weights, levelsOut, functionals, model.FunctionalDim);
        MulticanonicalCheckpoint? restart = checkpoint
            ? new MulticanonicalCheckpoint(identity, learning.Save(), learnHistogram.Snapshot(), [.. trace], null, summary.Checkpoint)
            : null;
        return new MulticanonicalResult<TState>(learning.Weights, learning.Save(), learnHistogram.Snapshot(), [.. trace], summary, draws,
            FlatHistogramDiagnostics.Merge(diagnostics), handles, restart, summary.SamplesUsed >= spec.ProductionDraws);
    }

    // The starting weights, t₀ and the gain's start, from the spec.
    private static (MulticanonicalWeights Start, double GainScale, long GainStart) Resolve(FlatHistogramSpec spec)
    {
        if (spec.NeedsPilot)
            throw new ArgumentException("A run without a pilot needs a flat start with its window, b_lo and K given.", nameof(spec));
        var shape = MulticanonicalWeights.Flat(spec.Window!.Low, spec.Window.High, spec.LowSlope!.Value, spec.Nodes!.Value);
        return (shape, spec.GainScale ?? 10.0 * shape.Nodes, spec.GainStart ?? 0);
    }

    private static string Fingerprint<TState>(FlatHistogramSpec spec, IMulticanonicalModel<TState> model)
    {
        var fingerprint = new RunFingerprint();
        spec.AddTo(fingerprint);
        fingerprint.Add(model.Fingerprint ?? string.Empty);
        return fingerprint.ToString();
    }

    private static IMulticanonicalChain<TState> Started<TState>(IMulticanonicalModel<TState> model, MulticanonicalTarget<TState> target,
                                                                int chainSeed, int readoutSeed)
        => model.StartChain(target, chainSeed, readoutSeed)
           ?? throw new InvalidOperationException($"{model.GetType().Name} started no chain.");

    private static IMulticanonicalChain<TState> Restored<TState>(IMulticanonicalModel<TState> model, MulticanonicalTarget<TState> target,
                                                                 IChainState state)
    {
        IMulticanonicalChain<TState> chain = model.RestoreChain(target, state)
            ?? throw new InvalidOperationException($"{model.GetType().Name} restored no chain.");
        chain.Reweigh();   // the target's weights are the checkpoint's
        return chain;
    }

    private static IChainState[] Saved<TState>(IMulticanonicalChain<TState>[] chains)
    {
        var states = new IChainState[chains.Length];
        for (int c = 0; c < chains.Length; c++) states[c] = chains[c].SaveState();
        return states;
    }

    // One production chain's restart state: the model's chain and the run's records of it.
    private sealed record ProductionState(IChainState Chain, double[] Levels, double[] Functionals, FlatHistogramDiagnostics Histogram) : IChainState;

    // The shell's view of the production phase: each chain the learning left, found by its seed, or restored.
    private sealed class ProductionModel<TState> : IEnsembleModel<TState>
    {
        private readonly IMulticanonicalModel<TState> _model;
        private readonly MulticanonicalTarget<TState> _target;
        private readonly Dictionary<int, Queue<IMulticanonicalChain<TState>>>? _learned;
        private readonly object _gate = new();

        public ProductionModel(IMulticanonicalModel<TState> model, MulticanonicalTarget<TState> target, string fingerprint,
                               IMulticanonicalChain<TState>[]? learned, int[] seeds)
        {
            _model = model;
            _target = target;
            Fingerprint = fingerprint;
            if (learned is null) return;
            _learned = [];
            for (int c = 0; c < learned.Length; c++)
            {
                if (!_learned.TryGetValue(seeds[c], out Queue<IMulticanonicalChain<TState>>? queue)) _learned[seeds[c]] = queue = new();
                queue.Enqueue(learned[c]);
            }
        }

        public int FunctionalDim => _model.FunctionalDim;
        public int MonitoredDim => _model.MonitoredDim;
        public IReadOnlyList<string> FunctionalNames => _model.FunctionalNames;
        public IReadOnlyList<FunctionalGroup> Groups => _model.Groups;
        public string Fingerprint { get; }

        public IChainRun<TState> StartChain(int chainSeed, int readoutSeed)
        {
            IMulticanonicalChain<TState> chain;
            lock (_gate)
            {
                // Two chains of one seed learned identically, so either serves.
                if (_learned is null || !_learned.TryGetValue(chainSeed, out Queue<IMulticanonicalChain<TState>>? queue) || queue.Count == 0)
                    throw new InvalidOperationException("No learned chain holds this seed.");
                chain = queue.Dequeue();
            }
            return new ProductionRun<TState>(chain, _target.Weights, _model.FunctionalDim);
        }

        public IChainRun<TState> RestoreChain(IChainState state)
        {
            if (state is not ProductionState saved) throw new ArgumentException("Not a production chain's state.", nameof(state));
            IMulticanonicalChain<TState> chain = Restored(_model, _target, saved.Chain);
            return new ProductionRun<TState>(chain, _target.Weights, _model.FunctionalDim, saved);
        }
    }

    // One production chain: the model's chain plus the run's records of its draws.
    private sealed class ProductionRun<TState> : IResumableChainRun<TState>
    {
        private readonly int _dim;

        public ProductionRun(IMulticanonicalChain<TState> chain, MulticanonicalWeights weights, int dim, ProductionState? saved = null)
        {
            Chain = chain;
            _dim = dim;
            Histogram = saved is null ? new HistogramAccumulator(weights, 1) : new HistogramAccumulator(weights, saved.Histogram);
            if (saved is null) return;
            Levels.AddRange(saved.Levels);
            Functionals.AddRange(saved.Functionals);
        }

        public IMulticanonicalChain<TState> Chain { get; }
        IChain<TState> IChainRun<TState>.Chain => Chain;
        public List<double> Levels { get; } = [];
        public List<double> Functionals { get; } = [];
        public HistogramAccumulator Histogram { get; }

        public void Burn(int steps)
        {
            for (int i = 0; i < steps; i++) Chain.Step();
        }

        public void Accumulate(in TState draw, Span<double> functionals)
        {
            Chain.Accumulate(draw, functionals);
            double level = Chain.CurrentEvaluation.UntemperedLogLikelihood;
            Levels.Add(level);
            for (int f = 0; f < _dim; f++) Functionals.Add(functionals[f]);
            Histogram.Add(0, level);
        }

        public IChainState SaveState() => new ProductionState(Chain.SaveState(), [.. Levels], [.. Functionals], Histogram.Snapshot());
    }
}
