using System;
using System.Collections.Generic;
using System.Threading;
using Maths.Regression.GPrior;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's client of the ensemble shell (<see cref="ChainEnsemble"/>; stage 4's S7 and S8): one <see cref="BassChain"/>
/// per chain, or with <see cref="BassSampling.Tempering"/> one <see cref="BassLadder"/>, started from an exact draw of
/// the prior on F (<see cref="BassPriorSampler"/>), each retained state composed into a <see cref="BassDraw"/>, kept in
/// the chain's <see cref="DrawStore{TDraw}"/>, and reduced to the functionals the shell monitors: M, τ when sampled, the
/// log posterior, σ², and the draw's fitted value at each monitored point, in D10's three groups (structural, level,
/// predictive). A ladder's draws are its b = 1 replica's.
/// <list type="bullet">
/// <item><b>Streams.</b> The shell's chain seed fans out into four: the start's draws, the moves, τ's block and a
/// ladder's own seed, of which an untempered chain uses the first three and a ladder the first and the last. Its readout
/// seed drives composition, so the trajectory depends neither on how often draws are retained nor on how many attempts
/// the start took.</item>
/// <item><b>Lengths.</b> Burn-in runs <see cref="BassSampling.BurnIn"/> sweeps, a ladder's scan counting as one; each
/// retained draw follows <see cref="BassSampling.Thin"/> more. An adaptive ladder re-spaces its rungs after rounds of 1,
/// 2, 4, … sweeps while the burn-in lasts (Syed et al.'s Algorithm 4, as BARS's ladders do) and runs the rest on the last
/// ladder, which sampling keeps (D10). A self-sizing ladder sizes itself after every round but the last, the longest,
/// which equalizes the count it settled on. At the end of burn-in a chain's tallies restart, so the account covers the
/// sampling phase, and a ladder that reports an evidence starts tracing its rungs.</item>
/// <item><b>Terms.</b> Each chain interns its draws' terms in its own <see cref="BassTermTable"/>; <see cref="Assemble"/>
/// merges them in chain order.</item>
/// <item><b>Restart state.</b> Each chain saves its chain's snapshot (the state, the design as it is, the three
/// streams, the counters and tallies) or its ladder's (every replica's, the swap and reference streams, the ladder, its
/// statistics and traces), its composition stream, its term table, its draws, its start and its move counts at the end of
/// burn-in, so a restored chain continues bit for bit.</item>
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
    private readonly BassPointInputs _pointInputs;
    private readonly TemperingSpec? _tempering;
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
        _tempering = model.Config.Sampling.Tempering;
        _pointInputs = new BassPointInputs(model);
        _points = _pointInputs.Points;
    }

    /// <summary>Every functional is monitored.</summary>
    public int FunctionalDim => MonitoredDim;

    /// <summary>M, τ when sampled, the log posterior, σ², and the fitted value at each monitored point, in that order: the
    /// structural group, then the level group, then the predictive group.</summary>
    public int MonitoredDim => (_model.Target.SamplesTau ? 4 : 3) + _points.Count;

    /// <summary>"M", "tau", "log_posterior", "sigma2", and "fitted[i]" at case i, or "fitted[i,r]" at case i and output
    /// location r, in the functionals' order.</summary>
    public IReadOnlyList<string> FunctionalNames => _pointInputs.FunctionalNames;

    /// <summary>D10's groups: <see cref="BassDiagnostics.PredictiveGroup"/> (σ² and the fitted values),
    /// <see cref="BassDiagnostics.LevelGroup"/> (the log posterior) and <see cref="BassDiagnostics.StructuralGroup"/> (M, and τ
    /// when it is sampled).</summary>
    public IReadOnlyList<FunctionalGroup> Groups
    {
        get
        {
            int level = _model.Target.SamplesTau ? 2 : 1;
            var predictive = new int[1 + _points.Count];
            for (int k = 0; k < predictive.Length; k++) predictive[k] = level + 1 + k;
            return
            [
                new FunctionalGroup(BassDiagnostics.PredictiveGroup, predictive),
                new FunctionalGroup(BassDiagnostics.LevelGroup, [level]),
                new FunctionalGroup(BassDiagnostics.StructuralGroup, _model.Target.SamplesTau ? [0, 1] : [0]),
            ];
        }
    }

    /// <summary><see cref="Bass.Fingerprint"/> of the config and data.</summary>
    public string Fingerprint { get; }

    /// <summary>Sweeps run so far over every chain, in this call.</summary>
    public long Sweeps => Interlocked.Read(ref _sweeps);

    /// <summary>
    /// Whether the run reports its evidence (D03, D05): asked for (<see cref="BassSampling.ReportEvidence"/>), tempered
    /// down to the exact reference rung, under proper priors: neither the flat intercept, the Jeffreys σ² prior nor an
    /// improper τ prior.
    /// </summary>
    public bool ReportsEvidence => _model.Config.Sampling.ReportEvidence && _tempering is { ReferenceRung: true }
                                   && _model.Target.Sigma2 is not null && _model.Target.Intercept != GPriorCentre.Flat && _model.Target.Scale.IsProper;

    /// <inheritdoc/>
    public IChainRun<BassState> StartChain(int chainSeed, int readoutSeed)
    {
        int[] seeds = SeedTree.Derive(chainSeed, 4);
        var sampler = new BassPriorSampler(_model);
        BassState start = sampler.DrawState(new Xoshiro256PlusPlus(seeds[0])).State;
        BassChain? chain = null;
        BassLadder? ladder = null;
        if (_tempering is { } spec) ladder = new BassLadder(_model, spec, start.Terms, start.LogTau, seeds[3]);
        else chain = new BassChain(_model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(seeds[1]), new Xoshiro256PlusPlus(seeds[2]));
        return new BassChainRun(this, chain, ladder, new Xoshiro256PlusPlus(readoutSeed), new BassTermTable(), new DrawStore<BassDraw>(),
                                new BassStart(start.Terms, _model.Target.HeldTau ?? start.Tau, sampler.Attempts));
    }

    /// <summary>Rebuilds a chain from the restart state its run saved (<see cref="BassChainState"/>).</summary>
    public IChainRun<BassState> RestoreChain(IChainState state)
    {
        if (state is not BassChainState s) throw new ArgumentException("Not a BASS chain's restart state.", nameof(state));
        if ((s.Ladder is null) != (_tempering is null))
            throw new ArgumentException("The chain ran " + (_tempering is null ? "tempered" : "untempered") + ".", nameof(state));
        BassChain? chain = s.Chain is { } c ? BassChain.Restore(_model, c) : null;
        BassLadder? ladder = s.Ladder is { } l ? BassLadder.Restore(_model, _tempering!, l) : null;
        (ulong c0, ulong c1, ulong c2, ulong c3) = s.Composition;
        var terms = new BassTermTable();
        foreach (BassTerm term in s.Terms) terms.Intern(term);
        return new BassChainRun(this, chain, ladder, Xoshiro256PlusPlus.Restore(c0, c1, c2, c3), terms, s.Draws.Copy(), s.Start)
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
        BassLadder[]? ladders = _tempering is null ? null : Array.ConvertAll(runs, r => r.Ladder!);
        BassDiagnostics diagnostics = Diagnostics(summary, runs, ladders);
        return new BassFit(config, _model.Resolved, _model.Transform,
                           new RunIdentity(config.Sampling.MasterSeed, config.Sampling.Chains, Fingerprint),
                           terms, draws, null, diagnostics,
                           ReportsEvidence ? BassEvidence.Of(_model, ladders!, diagnostics.Tempering!) : null,
                           BassResidualDependence.Of(_model, terms, draws), summary.Checkpoint);
    }

    private BassDiagnostics Diagnostics(EnsembleRun summary, BassChainRun[] runs, BassLadder[]? ladders)
    {
        var accounts = new BassChainAccount[runs.Length];
        var refusals = new Tally<BassRefusal>();
        var support = new BassSupportTally(_model.Prior);
        var kinds = new BassMoveTally();
        long emptyBirths = 0, atCap = 0, retained = 0;
        MoveStat[]? moves = null;
        for (int c = 0; c < runs.Length; c++)
        {
            BassChainRun run = runs[c];
            BassStructureMoves structure = run.Structure;
            IReadOnlyList<MoveStat> now = run.Moves.MoveStats();
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
            accounts[c] = new BassChainAccount(run.Start, sampling, structure.Kinds.Copy(), structure.Refusals.Copy(), structure.EmptyBirths,
                                               structure.Support.Copy(), (double)chainAtCap / run.Draws.Count);
            kinds.MergeFrom(structure.Kinds);
            refusals.MergeFrom(structure.Refusals);
            emptyBirths += structure.EmptyBirths;
            support.MergeFrom(structure.Support);
        }

        int f = 0;
        ConvergenceDiagnostic count = summary.Monitored[f++];
        ConvergenceDiagnostic? tau = _model.Target.SamplesTau ? summary.Monitored[f++] : null;
        int level = f;
        ConvergenceDiagnostic logPosterior = summary.Monitored[f++];
        ConvergenceDiagnostic sigma2 = summary.Monitored[f++];
        int first = f;
        var fitted = new BassMonitoredValue[_points.Count];
        for (int p = 0; p < fitted.Length; p++) fitted[p] = new BassMonitoredValue(_points[p], summary.Monitored[f++]);
        return new BassDiagnostics(count, tau, logPosterior, sigma2, fitted, moves!, summary.AcceptanceRate, kinds, refusals, emptyBirths,
                                   support, (double)atCap / retained, summary.SamplesUsed, summary.SamplesUsed < _plannedDraws,
                                   summary.Verdict, Separation(summary.Traces, level, first), accounts,
                                   ladders is null ? null : BassTemperingReadout.Of(ladders));
    }

    // The mode-separation readings from each chain's trace: its mean log posterior (functional `level`), and its mean fitted
    // values (the functionals from `first` on).
    private BassModeSeparation Separation(IReadOnlyList<ChainTrace> traces, int level, int first)
    {
        var logPosterior = new double[traces.Count];
        var fitted = new double[traces.Count][];
        for (int c = 0; c < traces.Count; c++)
        {
            logPosterior[c] = Mean(traces[c].Draws(level));
            fitted[c] = new double[_points.Count];
            for (int p = 0; p < _points.Count; p++) fitted[c][p] = Mean(traces[c].Draws(first + p));
        }
        return BassModeSeparation.Of(logPosterior, fitted);

        static double Mean(double[] draws)
        {
            double sum = 0.0;
            foreach (double d in draws) sum += d;
            return sum / draws.Length;
        }
    }

    // A draw's ids through `map`, sorted, with its coefficients (the intercept, then one per id) following them.
    internal static BassDraw Remap(BassDraw draw, int[] map)
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

    // One sweep of a chain: the moves, then τ; a ladder's scan of every replica, then the swaps. Cancellation is checked
    // before it, and progress counted after.
    private void Sweep(BassChainRun run)
    {
        _options.Cancellation.ThrowIfCancellationRequested();
        if (run.Ladder is { } ladder) ladder.Step();
        else run.Chain!.Step();
        long done = Interlocked.Increment(ref _sweeps);
        if (_options.Progress is { } progress && done % _options.ProgressInterval == 0)
            progress.Report(new BassProgress(done, _plannedSweeps));
    }

    // The retained draw at `state`, kept, and its functionals: M, τ (when sampled), the log posterior, σ², then the fitted value
    // f = a₀ + Σⱼ aⱼBⱼ at each monitored point. The log posterior is the collapsed ln π(structure, u | y) at b = 1, the
    // state's prior and log marginal, which is the chain's own density there (a ladder's b = 1 replica's).
    private void Accumulate(BassChainRun run, BassState state, Span<double> functionals)
    {
        BassDraw draw = run.Structure.Compose(state, run.Composition, run.Terms);
        run.Draws.Add(draw);
        int f = 0;
        functionals[f++] = state.Count;
        if (_model.Target.SamplesTau) functionals[f++] = draw.Tau;
        functionals[f++] = Temperature.LogPosterior(_model.Target.Evaluate(state));
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

    // The chain's term `id` at each monitored point: its case side at the point's case, times its output side at the
    // point's location; cached by the chain's own ids.
    private double[] PointValues(BassChainRun run, int id)
    {
        List<double[]?> cache = run.PointValues;
        while (cache.Count <= id) cache.Add(null);
        if (cache[id] is { } cached) return cached;
        return cache[id] = _pointInputs.Values(run.Terms[id]);
    }

    /// <summary>One chain's restart state: its <see cref="BassChain"/>'s snapshot or its <see cref="BassLadder"/>'s, the
    /// composition stream's words, its term table's terms in id order, its draws, its start and its move counts at the end
    /// of burn-in.</summary>
    private sealed record BassChainState(
        BassChainSnapshot? Chain, BassLadderSnapshot? Ladder, (ulong, ulong, ulong, ulong) Composition, BassTerm[] Terms,
        DrawStore<BassDraw> Draws, BassStart Start, IReadOnlyList<MoveStat> BurnMoves) : IChainState;

    /// <summary>One chain: its <see cref="BassChain"/> or its <see cref="BassLadder"/>, composition stream, term table,
    /// retained draws and start, and its move counts at the end of burn-in. As the shell's <see cref="IChain{TDraw}"/>, a
    /// step is <see cref="BassSampling.Thin"/> sweeps, and its acceptance counts the structure moves only, a ladder's b = 1
    /// replica's.</summary>
    private sealed class BassChainRun(BassEnsemble owner, BassChain? chain, BassLadder? ladder, Xoshiro256PlusPlus composition,
                                      BassTermTable terms, DrawStore<BassDraw> draws, BassStart start) : IResumableChainRun<BassState>, IChain<BassState>
    {
        public BassChain? Chain { get; } = chain;
        public BassLadder? Ladder { get; } = ladder;
        public Xoshiro256PlusPlus Composition { get; } = composition;
        public BassTermTable Terms { get; } = terms;
        public DrawStore<BassDraw> Draws { get; } = draws;
        public BassStart Start { get; } = start;
        public List<double[]?> PointValues { get; } = [];
        public IReadOnlyList<MoveStat> BurnMoves { get; set; } = [];

        // The sampling chain: the lone chain's, or the ladder's b = 1 replica, with its structure moves.
        public ReversibleJumpChain<BassState> Moves => Ladder?.Moves ?? Chain!.Moves;
        public BassStructureMoves Structure => Ladder?.Structure ?? Chain!.Structure;
        public BassState Current => Ladder?.Current ?? Chain!.Current;

        IChain<BassState> IChainRun<BassState>.Chain => this;
        public long Accepted => Moves.Accepted;
        public long Attempts => Moves.Attempts;

        public BassState Step()
        {
            for (int t = 0; t < owner._model.Config.Sampling.Thin; t++) owner.Sweep(this);
            return Current;
        }

        // The burn-in's sweeps, an adaptive ladder re-spaced after each round of 1, 2, 4, … while the burn-in lasts, a
        // self-sizing one sized after every round but the last, which equalizes the count it settled on; then the ladder
        // freezes, the tallies restart and the move counts are kept, so the account covers sampling.
        public void Burn(int steps)
        {
            int done = 0;
            if (Ladder is { Spec.Adaptive: true } adaptive)
                for (int round = 1; done + round <= steps; round *= 2)
                {
                    for (int k = 0; k < round; k++) owner.Sweep(this);
                    done += round;
                    if (adaptive.SizesItself && done + 2L * round <= steps) adaptive.Resize();
                    else adaptive.Adapt();
                }
            for (; done < steps; done++) owner.Sweep(this);
            if (Ladder is { } ladder) ladder.EndBurnIn(owner.ReportsEvidence);
            else Chain!.Structure.ResetTallies();
            BurnMoves = Moves.MoveStats();
        }

        public void Accumulate(in BassState draw, Span<double> functionals) => owner.Accumulate(this, draw, functionals);

        // A snapshot of everything above; the live chain is unchanged and continues as if none were taken.
        public IChainState SaveState()
        {
            var terms = new BassTerm[Terms.Count];
            for (int id = 0; id < terms.Length; id++) terms[id] = Terms[id];
            return new BassChainState(Chain?.Snapshot(), Ladder?.Snapshot(), Composition.SaveState(), terms, Draws.Copy(), Start, BurnMoves);
        }
    }
}

/// <summary>
/// The monitored points' inputs (<see cref="BassResolved.MonitoredPoints"/>): each point's case and, for a direct
/// functional response, its output location, so a term's value at every point is one evaluation; and the monitored
/// functionals' names, which the chain ensemble and SMC's between-run reading share.
/// </summary>
internal sealed class BassPointInputs
{
    private readonly BassModel _model;
    private readonly CaseInputs _cases;          // the monitored cases' inputs
    private readonly OutputInputs? _location;    // each monitored point's location, for a direct functional response

    public BassPointInputs(BassModel model)
    {
        _model = model;
        Points = model.Resolved.MonitoredPoints;
        BassInputs inputs = model.Inputs;
        CaseInputs cases = inputs.Cases;
        var continuous = new double[cases.ContinuousCount][];
        for (int v = 0; v < continuous.Length; v++)
        {
            continuous[v] = new double[Points.Count];
            for (int p = 0; p < Points.Count; p++) continuous[v][p] = cases.Continuous(v)[Points[p].Case];
        }
        var categorical = new int[cases.CategoricalCount][];
        for (int v = 0; v < categorical.Length; v++)
        {
            categorical[v] = new int[Points.Count];
            for (int p = 0; p < Points.Count; p++) categorical[v][p] = cases.Categorical(v)[Points[p].Case];
        }
        _cases = new CaseInputs(continuous, categorical);
        if (inputs.Outputs is { } outputs)
        {
            var at = new double[outputs.Count][];
            for (int v = 0; v < at.Length; v++)
            {
                at[v] = new double[Points.Count];
                for (int p = 0; p < Points.Count; p++) at[v][p] = outputs.Variable(v)[Points[p].Location!.Value];
            }
            _location = new OutputInputs(at);
        }
    }

    /// <summary>The monitored points.</summary>
    public IReadOnlyList<BassMonitoredPoint> Points { get; }

    /// <summary>"M", "tau" when τ is sampled, "log_posterior", "sigma2", and "fitted[i]" at case i, or "fitted[i,r]" at case i
    /// and output location r.</summary>
    public IReadOnlyList<string> FunctionalNames
    {
        get
        {
            var names = new List<string>(4 + Points.Count) { "M" };
            if (_model.Target.SamplesTau) names.Add("tau");
            names.Add("log_posterior");
            names.Add("sigma2");
            foreach (BassMonitoredPoint p in Points)
                names.Add(p.Location is int r ? $"fitted[{p.Case},{r}]" : $"fitted[{p.Case}]");
            return names;
        }
    }

    /// <summary><paramref name="term"/> at each monitored point: its case side at the point's case, times its output side at
    /// the point's location.</summary>
    public double[] Values(BassTerm term)
    {
        var values = new double[Points.Count];
        _model.Basis.CaseColumn(term, _cases, values);
        if (_location is { } at)
        {
            var h = new double[Points.Count];
            _model.Basis.OutputColumn(term, at, h);
            for (int p = 0; p < values.Length; p++) values[p] *= h[p];
        }
        return values;
    }
}
