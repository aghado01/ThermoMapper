using System;
using System.Collections.Generic;
using System.Threading;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// BASS's client of the flat-histogram run (N17 Phase 2; <see cref="MulticanonicalRun"/>): one <see cref="BassChain"/> per
/// chain on the <see cref="MulticanonicalTarget{TState}"/> over the collapsed target, each started from an exact draw of
/// the prior on F, each production draw composed into a <see cref="BassDraw"/>, kept in the chain's
/// <see cref="DrawStore{TDraw}"/>, and reduced to the functionals the ensemble monitors: M, τ when sampled, the log
/// posterior, σ², and the fitted value at each monitored point, in D10's three groups, as an MCMC fit's are. D10's verdict
/// then reads them under π_w, the chains' own target.
/// <list type="bullet">
/// <item><b>Streams.</b> The run's chain seed fans out into three, the start's draws, the moves and τ's block; its readout
/// seed drives composition, as an MCMC fit's does.</item>
/// <item><b>The account.</b> A chain's tallies restart when it takes its first production step, its step count then
/// equal to the learning's length, so its account covers production, as an MCMC chain's covers sampling.</item>
/// </list>
/// The model is shared read-only; each chain owns its design, streams and stores.
/// </summary>
internal sealed class BassMulticanonicalModel : IMulticanonicalModel<BassState>
{
    private readonly BassPointInputs _pointInputs;

    /// <summary>The client for <paramref name="model"/>'s run, whose learning takes <paramref name="learningIterations"/>
    /// steps per chain, its sweeps counted into <paramref name="sweeps"/>.</summary>
    public BassMulticanonicalModel(BassModel model, string fingerprint, long learningIterations, BassSweepCount sweeps)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(fingerprint);
        ArgumentNullException.ThrowIfNull(sweeps);
        Model = model;
        Fingerprint = fingerprint;
        LearningIterations = learningIterations;
        Sweeps = sweeps;
        _pointInputs = new BassPointInputs(model);
    }

    /// <summary>The model.</summary>
    public BassModel Model { get; }

    /// <summary>The learning's length, per chain.</summary>
    public long LearningIterations { get; }

    /// <summary>The run's sweep count, which reports progress and carries cancellation.</summary>
    public BassSweepCount Sweeps { get; }

    /// <summary>The monitored points' inputs.</summary>
    public BassPointInputs PointInputs => _pointInputs;

    /// <summary>Every functional is monitored.</summary>
    public int FunctionalDim => MonitoredDim;

    /// <summary>M, τ when sampled, the log posterior, σ², and the fitted value at each monitored point.</summary>
    public int MonitoredDim => (Model.Target.SamplesTau ? 4 : 3) + _pointInputs.Points.Count;

    /// <summary>The functionals' names, an MCMC fit's.</summary>
    public IReadOnlyList<string> FunctionalNames => _pointInputs.FunctionalNames;

    /// <summary>D10's groups, an MCMC fit's: predictive (σ² and the fitted values), level (the log posterior) and
    /// structural (M, and τ when it is sampled).</summary>
    public IReadOnlyList<FunctionalGroup> Groups
    {
        get
        {
            int level = LevelIndex;
            var predictive = new int[1 + _pointInputs.Points.Count];
            for (int k = 0; k < predictive.Length; k++) predictive[k] = level + 1 + k;
            return
            [
                new FunctionalGroup(BassDiagnostics.PredictiveGroup, predictive),
                new FunctionalGroup(BassDiagnostics.LevelGroup, [level]),
                new FunctionalGroup(BassDiagnostics.StructuralGroup, Model.Target.SamplesTau ? [0, 1] : [0]),
            ];
        }
    }

    /// <summary>The log posterior's functional.</summary>
    public int LevelIndex => Model.Target.SamplesTau ? 2 : 1;

    /// <summary><see cref="Bass.Fingerprint"/> of the config and data.</summary>
    public string Fingerprint { get; }

    /// <inheritdoc/>
    public IMulticanonicalChain<BassState> StartChain(MulticanonicalTarget<BassState> target, int chainSeed, int readoutSeed)
    {
        ArgumentNullException.ThrowIfNull(target);
        int[] seeds = SeedTree.Derive(chainSeed, 3);   // the start's draws, the moves, τ's block
        var sampler = new BassPriorSampler(Model);
        BassState start = sampler.DrawState(new Xoshiro256PlusPlus(seeds[0])).State;
        var chain = new BassChain(Model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(seeds[1]), new Xoshiro256PlusPlus(seeds[2]), target);
        return new BassMulticanonicalChain(this, chain, target, new Xoshiro256PlusPlus(readoutSeed), new BassTermTable(), new DrawStore<BassDraw>(),
                                           new BassStart(start.Terms, Model.Target.HeldTau ?? start.Tau, sampler.Attempts), null);
    }

    /// <inheritdoc/>
    public IMulticanonicalChain<BassState> RestoreChain(MulticanonicalTarget<BassState> target, IChainState state)
    {
        ArgumentNullException.ThrowIfNull(target);
        if (state is not BassMulticanonicalChainState s) throw new ArgumentException("Not a BASS flat-histogram chain's restart state.", nameof(state));
        BassChain chain = BassChain.Restore(Model, s.Chain, target);
        (ulong c0, ulong c1, ulong c2, ulong c3) = s.Composition;
        var terms = new BassTermTable();
        foreach (BassTerm term in s.Terms) terms.Intern(term);
        return new BassMulticanonicalChain(this, chain, target, Xoshiro256PlusPlus.Restore(c0, c1, c2, c3), terms, s.Draws.Copy(), s.Start,
                                           s.LearningMoves);
    }
}

/// <summary>
/// One BASS chain on the multicanonical target (<see cref="IMulticanonicalChain{TState}"/>): D04's sweep, the structure
/// moves and τ's slice step, each reading π_w through the target, so a move accepts on the inner prior and proposal
/// ratios with ln w(ℓ₁) in the likelihood's place, and τ's slice targets p(u)·w(ℓ₁(S, u)). <see cref="Reweigh"/> re-reads
/// ln w from the current evaluation's cached ℓ₁, with no new marginal. A production draw is composed and kept, as an MCMC
/// fit's retained draw is, with its full state (D04).
/// </summary>
internal sealed class BassMulticanonicalChain : IMulticanonicalChain<BassState>
{
    private readonly BassMulticanonicalModel _owner;
    private readonly MulticanonicalTarget<BassState> _target;
    private readonly List<double[]?> _pointValues = [];

    public BassMulticanonicalChain(BassMulticanonicalModel owner, BassChain chain, MulticanonicalTarget<BassState> target, Xoshiro256PlusPlus composition,
                                   BassTermTable terms, DrawStore<BassDraw> draws, BassStart start, IReadOnlyList<MoveStat>? learningMoves)
    {
        _owner = owner;
        Chain = chain;
        _target = target;
        Composition = composition;
        Terms = terms;
        Draws = draws;
        Start = start;
        LearningMoves = learningMoves;
        if (Chain.Moves.CurrentEvaluation is not MulticanonicalEvaluation)
            throw new ArgumentException("The chain does not sample a multicanonical target.", nameof(chain));
    }

    /// <summary>D04's sweep on π_w.</summary>
    public BassChain Chain { get; }

    /// <summary>The composition stream.</summary>
    public Xoshiro256PlusPlus Composition { get; }

    /// <summary>The chain's own term table, which its draws' ids index.</summary>
    public BassTermTable Terms { get; }

    /// <summary>The production draws, in order.</summary>
    public DrawStore<BassDraw> Draws { get; }

    /// <summary>The chain's start.</summary>
    public BassStart Start { get; }

    /// <summary>The move counts when production began, from which its account counts; null while learning.</summary>
    public IReadOnlyList<MoveStat>? LearningMoves { get; private set; }

    /// <inheritdoc/>
    public long Accepted => Chain.Moves.Accepted;

    /// <inheritdoc/>
    public long Attempts => Chain.Moves.Attempts;

    /// <inheritdoc/>
    public MulticanonicalEvaluation CurrentEvaluation => (MulticanonicalEvaluation)Chain.Moves.CurrentEvaluation;

    /// <summary>One sweep on π_w. The first production sweep restarts the tallies first.</summary>
    public BassState Step()
    {
        _owner.Sweeps.Cancellation.ThrowIfCancellationRequested();
        if (LearningMoves is null && Chain.Scan.Steps == _owner.LearningIterations)
        {
            Chain.Structure.ResetTallies();
            LearningMoves = Chain.Moves.MoveStats();
        }
        BassState state = Chain.Step();
        _owner.Sweeps.Add(1);
        return state;
    }

    /// <inheritdoc/>
    public void Reweigh() => Chain.Moves.SetState(Chain.Moves.Current, CurrentEvaluation.Reweighed(_target.Weights));

    /// <summary>
    /// The production draw at <paramref name="draw"/>, composed and kept, and its functionals: M, τ (when sampled), the log
    /// posterior at b = 1 from the cached evaluation (the inner prior and ℓ₁), σ², then the fitted value at each monitored
    /// point.
    /// </summary>
    public void Accumulate(in BassState draw, Span<double> functionals)
    {
        BassDraw composed = Chain.Structure.Compose(draw, Composition, Terms);
        Draws.Add(composed);
        MulticanonicalEvaluation evaluation = CurrentEvaluation;
        int f = 0;
        functionals[f++] = draw.Count;
        if (_owner.Model.Target.SamplesTau) functionals[f++] = composed.Tau;
        functionals[f++] = evaluation.LogPrior + evaluation.UntemperedLogLikelihood;
        functionals[f++] = composed.Sigma2;
        Span<double> fitted = functionals[f..];
        fitted.Fill(composed.Coefficients[0]);
        int[] ids = composed.Structure.TermIds;
        for (int k = 0; k < ids.Length; k++)
        {
            double[] values = PointValues(ids[k]);
            double a = composed.Coefficients[k + 1];
            for (int p = 0; p < fitted.Length; p++) fitted[p] += a * values[p];
        }
    }

    /// <inheritdoc/>
    public IChainState SaveState()
    {
        var terms = new BassTerm[Terms.Count];
        for (int id = 0; id < terms.Length; id++) terms[id] = Terms[id];
        return new BassMulticanonicalChainState(Chain.Snapshot(), Composition.SaveState(), terms, Draws.Copy(), Start, LearningMoves);
    }

    // The chain's term `id` at each monitored point, cached by the chain's own ids.
    private double[] PointValues(int id)
    {
        while (_pointValues.Count <= id) _pointValues.Add(null);
        return _pointValues[id] ??= _owner.PointInputs.Values(Terms[id]);
    }
}

/// <summary>A BASS flat-histogram chain's restart state: its <see cref="BassChain"/>'s snapshot, the composition stream's
/// words, its term table's terms in id order, its production draws, its start and its move counts when production
/// began (null while learning).</summary>
internal sealed record BassMulticanonicalChainState(
    BassChainSnapshot Chain, (ulong, ulong, ulong, ulong) Composition, BassTerm[] Terms, DrawStore<BassDraw> Draws, BassStart Start,
    IReadOnlyList<MoveStat>? LearningMoves) : IChainState;

/// <summary>
/// A run's sweeps over every chain, counted from any thread: each time they pass a multiple of the options' interval the
/// progress hears of it, and once at the end; it carries the run's cancellation.
/// </summary>
internal sealed class BassSweepCount(BassRunOptions options, long planned)
{
    private long _sweeps;

    /// <summary>The run's cancellation.</summary>
    public CancellationToken Cancellation => options.Cancellation;

    /// <summary>The sweeps so far.</summary>
    public long Sweeps => Interlocked.Read(ref _sweeps);

    /// <summary>Counts <paramref name="sweeps"/> more.</summary>
    public void Add(long sweeps)
    {
        long after = Interlocked.Add(ref _sweeps, sweeps), before = after - sweeps;
        if (options.Progress is { } progress && after / options.ProgressInterval > before / options.ProgressInterval)
            progress.Report(new BassProgress(after, planned));
    }

    /// <summary>The final report, unless the last count made one.</summary>
    public void End()
    {
        if (options.Progress is not { } progress) return;
        long done = Sweeps;
        if (done == 0 || done % options.ProgressInterval != 0) progress.Report(new BassProgress(done, planned));
    }
}
