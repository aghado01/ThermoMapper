using System;
using System.Collections.Generic;
using System.Threading;
using Maths.Rng;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's client of the ensemble shell (<see cref="ChainEnsemble"/>; stage 4's S7): one <see cref="BassChain"/> per
/// chain, started from an exact draw of the prior on F (<see cref="BassPriorSampler"/>), each retained state composed
/// into a <see cref="BassDraw"/>, kept in the chain's <see cref="DrawStore{TDraw}"/>, and reduced to the functionals the
/// shell monitors: M, τ in the g mode, and σ².
/// <list type="bullet">
/// <item><b>Streams.</b> The shell's chain seed fans out into three: the start's draws, the moves and τ's block. Its
/// readout seed drives composition, so the trajectory depends neither on how often draws are retained nor on how many
/// attempts the start took.</item>
/// <item><b>Lengths.</b> Burn-in runs <see cref="BassSampling.BurnIn"/> sweeps; each retained draw follows
/// <see cref="BassSampling.Thin"/> more.</item>
/// <item><b>Terms.</b> Each chain interns its draws' terms in its own <see cref="BassTermTable"/>; <see cref="Assemble"/>
/// merges them in chain order.</item>
/// </list>
/// The model is shared read-only; each chain owns its design, streams and stores.
/// </summary>
internal sealed class BassEnsemble : IEnsembleModel<BassState>
{
    private readonly BassModel _model;
    private readonly BassRunOptions _options;
    private readonly long _plannedSweeps;
    private long _sweeps;

    public BassEnsemble(BassModel model, string fingerprint, BassRunOptions options, long plannedSweeps)
    {
        _model = model;
        Fingerprint = fingerprint;
        _options = options;
        _plannedSweeps = plannedSweeps;
    }

    /// <summary>M, τ in the g mode, and σ², each monitored.</summary>
    public int FunctionalDim => MonitoredDim;

    /// <inheritdoc/>
    public int MonitoredDim => _model.Target.SamplesTau ? 3 : 2;

    /// <summary><see cref="Bass.Fingerprint"/> of the config and data.</summary>
    public string Fingerprint { get; }

    /// <summary>Sweeps run so far over every chain, in this call.</summary>
    public long Sweeps => Interlocked.Read(ref _sweeps);

    /// <inheritdoc/>
    public IChainRun<BassState> StartChain(int chainSeed, int readoutSeed)
    {
        int[] seeds = SeedTree.Derive(chainSeed, 3);
        var sampler = new BassPriorSampler(_model);
        BassState start = sampler.DrawState(new Xoshiro256PlusPlus(seeds[0])).State;
        var chain = new BassChain(_model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(seeds[1]), new Xoshiro256PlusPlus(seeds[2]));
        return new BassChainRun(this, chain, new Xoshiro256PlusPlus(readoutSeed), new BassTermTable(), new DrawStore<BassDraw>());
    }

    /// <summary>The final progress report, unless the last sweep made one.</summary>
    public void ReportEnd()
    {
        if (_options.Progress is not { } progress) return;
        long done = Sweeps;
        if (done == 0 || done % _options.ProgressInterval != 0) progress.Report(new BassProgress(done, _plannedSweeps));
    }

    /// <summary>The fit from the shell's summary and its chains: the chains' tables merged in chain order, each draw's ids
    /// remapped into the merged table and re-sorted with its coefficients.</summary>
    public BassFit Assemble(EnsembleRun summary, IReadOnlyList<IChainRun<BassState>> chains)
    {
        var runs = new BassChainRun[chains.Count];
        for (int c = 0; c < runs.Length; c++) runs[c] = (BassChainRun)chains[c];
        var terms = new BassTermTable();
        var draws = new DrawStore<BassDraw>[runs.Length];
        for (int c = 0; c < runs.Length; c++)
        {
            BassTermTable own = runs[c].Terms;
            var map = new int[own.Count];
            for (int id = 0; id < map.Length; id++) map[id] = terms.Intern(own[id]);
            draws[c] = new DrawStore<BassDraw>();
            foreach (BassDraw d in runs[c].Draws) draws[c].Add(Remap(d, map));
        }

        BassConfig config = _model.Config;
        var diagnostics = new BassDiagnostics(summary.SamplesUsed, summary.AcceptanceRate);
        return new BassFit(config, _model.Resolved, _model.Transform,
                           new RunIdentity(config.Sampling.MasterSeed, config.Sampling.Chains, Fingerprint),
                           terms, draws, null, diagnostics, summary.Checkpoint);
    }

    // A draw's ids through `map`, sorted, with its coefficients (the intercept, then one per id) following them.
    private static BassDraw Remap(BassDraw draw, int[] map)
    {
        int[] ids = draw.Structure.TermIds;
        var mapped = new int[ids.Length];
        bool same = true;
        for (int k = 0; k < ids.Length; k++) same &= (mapped[k] = map[ids[k]]) == ids[k];
        if (same) return draw;
        var order = new int[ids.Length];
        for (int k = 0; k < order.Length; k++) order[k] = k;
        Array.Sort((int[])mapped.Clone(), order);
        var sorted = new int[ids.Length];
        var coefficients = new double[draw.Coefficients.Length];
        coefficients[0] = draw.Coefficients[0];
        for (int k = 0; k < order.Length; k++)
        {
            sorted[k] = mapped[order[k]];
            coefficients[k + 1] = draw.Coefficients[order[k] + 1];
        }
        return draw with { Structure = new BassStructure(sorted), Coefficients = coefficients };
    }

    // One sweep of a chain: the moves, then τ. Cancellation is checked before it, and progress counted after.
    private void Sweep(BassChainRun run)
    {
        _options.Cancellation.ThrowIfCancellationRequested();
        run.Chain.Step();
        long done = Interlocked.Increment(ref _sweeps);
        if (_options.Progress is { } progress && done % _options.ProgressInterval == 0)
            progress.Report(new BassProgress(done, _plannedSweeps));
    }

    // The retained draw at `state`, kept, and its functionals: M, τ (g mode), σ².
    private void Accumulate(BassChainRun run, BassState state, Span<double> functionals)
    {
        BassDraw draw = run.Chain.Structure.Compose(state, run.Composition, run.Terms);
        run.Draws.Add(draw);
        int f = 0;
        functionals[f++] = state.Count;
        if (_model.Target.SamplesTau) functionals[f++] = draw.Tau;
        functionals[f] = draw.Sigma2;
    }

    /// <summary>One chain: its <see cref="BassChain"/>, composition stream, term table and retained draws. As the
    /// shell's <see cref="IChain{TDraw}"/>, a step is <see cref="BassSampling.Thin"/> sweeps, and its acceptance counts
    /// the structure moves only.</summary>
    private sealed class BassChainRun(BassEnsemble owner, BassChain chain, Xoshiro256PlusPlus composition, BassTermTable terms,
                                      DrawStore<BassDraw> draws) : IChainRun<BassState>, IChain<BassState>
    {
        public BassChain Chain { get; } = chain;
        public Xoshiro256PlusPlus Composition { get; } = composition;
        public BassTermTable Terms { get; } = terms;
        public DrawStore<BassDraw> Draws { get; } = draws;

        IChain<BassState> IChainRun<BassState>.Chain => this;
        public long Accepted => Chain.Moves.Accepted;
        public long Attempts => Chain.Moves.Attempts;

        public BassState Step()
        {
            for (int t = 0; t < owner._model.Config.Sampling.Thin; t++) owner.Sweep(this);
            return Chain.Current;
        }

        public void Burn(int steps)
        {
            for (int t = 0; t < steps; t++) owner.Sweep(this);
        }

        public void Accumulate(in BassState draw, Span<double> functionals) => owner.Accumulate(this, draw, functionals);
    }
}
