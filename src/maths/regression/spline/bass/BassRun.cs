using System;
using System.Collections.Generic;
using System.Threading;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's client of the ensemble shell (<see cref="ChainEnsemble"/>; stage 4's S7): one <see cref="BassChain"/> per
/// chain, started from an exact draw of the prior on F (<see cref="BassPriorSampler"/>), each retained state composed
/// into a <see cref="BassDraw"/>, kept in the chain's <see cref="DrawStore{TDraw}"/>, and reduced to the functionals the
/// shell monitors: M, τ in the g mode, σ², and the draw's fitted value at each monitored point.
/// <list type="bullet">
/// <item><b>Streams.</b> The shell's chain seed fans out into three: the start's draws, the moves and τ's block. Its
/// readout seed drives composition, so the trajectory depends neither on how often draws are retained nor on how many
/// attempts the start took.</item>
/// <item><b>Lengths.</b> Burn-in runs <see cref="BassSampling.BurnIn"/> sweeps; each retained draw follows
/// <see cref="BassSampling.Thin"/> more. At the end of burn-in a chain's tallies restart, so the account covers the
/// sampling phase.</item>
/// <item><b>Terms.</b> Each chain interns its draws' terms in its own <see cref="BassTermTable"/>; <see cref="Assemble"/>
/// merges them in chain order.</item>
/// <item><b>Restart state.</b> Each chain saves its chain's snapshot (the state, the design as it is, the three
/// streams, the counters and tallies), its composition stream, its term table, its draws, its start and its move counts at
/// the end of burn-in, so a restored chain continues bit for bit.</item>
/// </list>
/// The model is shared read-only; each chain owns its design, streams and stores.
/// </summary>
internal sealed class BassEnsemble : IEnsembleModel<BassState>
{
    private readonly BassModel _model;
    private readonly BassRunOptions _options;
    private readonly long _plannedSweeps;
    private readonly int _plannedDraws;
    private readonly IReadOnlyList<BassMonitoredPoint> _points;
    private readonly CaseInputs _pointCases;         // the monitored cases' inputs
    private readonly OutputInputs? _pointLocation;   // the monitored location's, for a direct functional response
    private long _sweeps;

    /// <summary>The client for <paramref name="model"/>'s run: <paramref name="plannedSweeps"/> over every chain, ending
    /// with <paramref name="plannedDraws"/> per chain unless the R̂ stop ends it sooner.</summary>
    public BassEnsemble(BassModel model, string fingerprint, BassRunOptions options, long plannedSweeps, int plannedDraws)
    {
        _model = model;
        Fingerprint = fingerprint;
        _options = options;
        _plannedSweeps = plannedSweeps;
        _plannedDraws = plannedDraws;
        BassInputs inputs = model.Inputs;
        _points = BassDiagnostics.DefaultPoints(model.Resolved, inputs.Outputs is not null);

        CaseInputs cases = inputs.Cases;
        var continuous = new double[cases.ContinuousCount][];
        for (int v = 0; v < continuous.Length; v++)
        {
            continuous[v] = new double[_points.Count];
            for (int p = 0; p < _points.Count; p++) continuous[v][p] = cases.Continuous(v)[_points[p].Case];
        }
        var categorical = new int[cases.CategoricalCount][];
        for (int v = 0; v < categorical.Length; v++)
        {
            categorical[v] = new int[_points.Count];
            for (int p = 0; p < _points.Count; p++) categorical[v][p] = cases.Categorical(v)[_points[p].Case];
        }
        _pointCases = new CaseInputs(continuous, categorical);
        if (inputs.Outputs is { } outputs && _points[0].Location is int location)
        {
            var at = new double[outputs.Count][];
            for (int v = 0; v < at.Length; v++) at[v] = [outputs.Variable(v)[location]];
            _pointLocation = new OutputInputs(at);
        }
    }

    /// <summary>Every functional is monitored.</summary>
    public int FunctionalDim => MonitoredDim;

    /// <summary>M, τ in the g mode, σ², and the fitted value at each monitored point, in that order.</summary>
    public int MonitoredDim => (_model.Target.SamplesTau ? 3 : 2) + _points.Count;

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
        return new BassChainRun(this, chain, new Xoshiro256PlusPlus(readoutSeed), new BassTermTable(), new DrawStore<BassDraw>(),
                                new BassStart(start.Terms, _model.Target.FixedTau ?? start.Tau, sampler.Attempts));
    }

    /// <summary>Rebuilds a chain from the restart state its run saved (<see cref="BassChainState"/>).</summary>
    public IChainRun<BassState> RestoreChain(IChainState state)
    {
        if (state is not BassChainState s) throw new ArgumentException("Not a BASS chain's restart state.", nameof(state));
        BassChain chain = BassChain.Restore(_model, s.Chain);
        (ulong c0, ulong c1, ulong c2, ulong c3) = s.Composition;
        var terms = new BassTermTable();
        foreach (BassTerm term in s.Terms) terms.Intern(term);
        return new BassChainRun(this, chain, Xoshiro256PlusPlus.Restore(c0, c1, c2, c3), terms, s.Draws.Copy(), s.Start)
        {
            BurnMoves = s.BurnMoves,
        };
    }

    /// <summary>The final progress report, unless the last sweep made one.</summary>
    public void ReportEnd()
    {
        if (_options.Progress is not { } progress) return;
        long done = Sweeps;
        if (done == 0 || done % _options.ProgressInterval != 0) progress.Report(new BassProgress(done, _plannedSweeps));
    }

    /// <summary>
    /// The fit from the shell's summary and its chains: the chains' tables merged in chain order, each draw's ids remapped
    /// into the merged table and re-sorted with its coefficients, and the diagnostics, the chains' accounts summed.
    /// </summary>
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
        return new BassFit(config, _model.Resolved, _model.Transform,
                           new RunIdentity(config.Sampling.MasterSeed, config.Sampling.Chains, Fingerprint),
                           terms, draws, null, Diagnostics(summary, runs), null, summary.Checkpoint);
    }

    private BassDiagnostics Diagnostics(EnsembleRun summary, BassChainRun[] runs)
    {
        var accounts = new BassChainAccount[runs.Length];
        var refusals = new Tally<BassRefusal>();
        var support = new BassSupportTally(_model.Prior);
        long emptyBirths = 0, atCap = 0, retained = 0;
        MoveStat[]? moves = null;
        for (int c = 0; c < runs.Length; c++)
        {
            BassChainRun run = runs[c];
            BassStructureMoves structure = run.Chain.Structure;
            IReadOnlyList<MoveStat> now = run.Chain.Moves.MoveStats();
            var sampling = new MoveStat[now.Count];
            for (int i = 0; i < now.Count; i++)
                sampling[i] = new MoveStat(now[i].Key, now[i].Attempts - run.BurnMoves[i].Attempts, now[i].Accepted - run.BurnMoves[i].Accepted);
            moves ??= new MoveStat[now.Count];
            for (int i = 0; i < now.Count; i++)
                moves[i] = new MoveStat(sampling[i].Key, moves[i].Attempts + sampling[i].Attempts, moves[i].Accepted + sampling[i].Accepted);

            long chainAtCap = 0;
            foreach (BassDraw d in run.Draws)
                if (d.Structure.TermIds.Length == _model.MaxBasis) chainAtCap++;
            atCap += chainAtCap;
            retained += run.Draws.Count;
            accounts[c] = new BassChainAccount(run.Start, sampling, structure.Refusals.Copy(), structure.EmptyBirths, structure.Support.Copy(),
                                               (double)chainAtCap / run.Draws.Count);
            refusals.MergeFrom(structure.Refusals);
            emptyBirths += structure.EmptyBirths;
            support.MergeFrom(structure.Support);
        }

        int f = 0;
        ConvergenceDiagnostic count = summary.Monitored[f++];
        ConvergenceDiagnostic? tau = _model.Target.SamplesTau ? summary.Monitored[f++] : null;
        ConvergenceDiagnostic sigma2 = summary.Monitored[f++];
        var fitted = new BassMonitoredValue[_points.Count];
        for (int p = 0; p < fitted.Length; p++) fitted[p] = new BassMonitoredValue(_points[p], summary.Monitored[f++]);
        return new BassDiagnostics(count, tau, sigma2, fitted, moves!, summary.AcceptanceRate, refusals, emptyBirths, support,
                                   (double)atCap / retained, summary.SamplesUsed, summary.SamplesUsed < _plannedDraws, accounts);
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

    // The retained draw at `state`, kept, and its functionals: M, τ (g mode), σ², then the fitted value
    // f = a₀ + Σⱼ aⱼBⱼ at each monitored point.
    private void Accumulate(BassChainRun run, BassState state, Span<double> functionals)
    {
        BassDraw draw = run.Chain.Structure.Compose(state, run.Composition, run.Terms);
        run.Draws.Add(draw);
        int f = 0;
        functionals[f++] = state.Count;
        if (_model.Target.SamplesTau) functionals[f++] = draw.Tau;
        functionals[f++] = draw.Sigma2;
        Span<double> fitted = functionals[f..];
        fitted.Fill(draw.Coefficients[0]);
        int[] ids = draw.Structure.TermIds;
        for (int k = 0; k < ids.Length; k++)
        {
            double[] values = PointValues(run, ids[k]);
            double a = draw.Coefficients[k + 1];
            for (int p = 0; p < fitted.Length; p++) fitted[p] += a * values[p];
        }
    }

    // The chain's term `id` at each monitored point: its case side there, times its output side at the monitored
    // location; cached by the chain's own ids.
    private double[] PointValues(BassChainRun run, int id)
    {
        List<double[]?> cache = run.PointValues;
        while (cache.Count <= id) cache.Add(null);
        if (cache[id] is { } cached) return cached;
        BassTerm term = run.Terms[id];
        var values = new double[_points.Count];
        _model.Basis.CaseColumn(term, _pointCases, values);
        if (_pointLocation is { } at)
        {
            Span<double> h = stackalloc double[1];
            _model.Basis.OutputColumn(term, at, h);
            for (int p = 0; p < values.Length; p++) values[p] *= h[0];
        }
        return cache[id] = values;
    }

    /// <summary>One chain's restart state: its <see cref="BassChain"/>'s snapshot, the composition stream's words, its
    /// term table's terms in id order, its draws, its start and its move counts at the end of burn-in.</summary>
    private sealed record BassChainState(
        BassChainSnapshot Chain, (ulong, ulong, ulong, ulong) Composition, BassTerm[] Terms, DrawStore<BassDraw> Draws, BassStart Start,
        IReadOnlyList<MoveStat> BurnMoves) : IChainState;

    /// <summary>One chain: its <see cref="BassChain"/>, composition stream, term table, retained draws and start, and its
    /// move counts at the end of burn-in. As the shell's <see cref="IChain{TDraw}"/>, a step is
    /// <see cref="BassSampling.Thin"/> sweeps, and its acceptance counts the structure moves only.</summary>
    private sealed class BassChainRun(BassEnsemble owner, BassChain chain, Xoshiro256PlusPlus composition, BassTermTable terms,
                                      DrawStore<BassDraw> draws, BassStart start) : IResumableChainRun<BassState>, IChain<BassState>
    {
        public BassChain Chain { get; } = chain;
        public Xoshiro256PlusPlus Composition { get; } = composition;
        public BassTermTable Terms { get; } = terms;
        public DrawStore<BassDraw> Draws { get; } = draws;
        public BassStart Start { get; } = start;
        public List<double[]?> PointValues { get; } = [];
        public IReadOnlyList<MoveStat> BurnMoves { get; set; } = [];

        IChain<BassState> IChainRun<BassState>.Chain => this;
        public long Accepted => Chain.Moves.Accepted;
        public long Attempts => Chain.Moves.Attempts;

        public BassState Step()
        {
            for (int t = 0; t < owner._model.Config.Sampling.Thin; t++) owner.Sweep(this);
            return Chain.Current;
        }

        // The burn-in's sweeps; then the tallies restart and the move counts are kept, so the account covers sampling.
        public void Burn(int steps)
        {
            for (int t = 0; t < steps; t++) owner.Sweep(this);
            Chain.Structure.ResetTallies();
            BurnMoves = Chain.Moves.MoveStats();
        }

        public void Accumulate(in BassState draw, Span<double> functionals) => owner.Accumulate(this, draw, functionals);

        // A snapshot of everything above; the live chain is unchanged and continues as if none were taken.
        public IChainState SaveState()
        {
            var terms = new BassTerm[Terms.Count];
            for (int id = 0; id < terms.Length; id++) terms[id] = Terms[id];
            return new BassChainState(Chain.Snapshot(), Composition.SaveState(), terms, Draws.Copy(), Start, BurnMoves);
        }
    }
}
