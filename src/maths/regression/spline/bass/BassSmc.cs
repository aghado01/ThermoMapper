using System;
using System.Collections.Generic;
using System.Threading;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// One SMC particle (stage 8): a state of D04's chain with its own structure moves, so its design, counts and tallies
/// follow it alone. The run holds each particle in one slot; resampling's further offspring are forks
/// (<see cref="BassParticles.Fork"/>). Not thread-safe: one mover at a time.
/// </summary>
public sealed class BassParticle
{
    internal BassParticle(BassStructureMoves moves, BassState state)
    {
        Moves = moves;
        State = state;
    }

    /// <summary>The particle's state.</summary>
    public BassState State { get; internal set; }

    /// <summary>Its structure moves, over its own design.</summary>
    public BassStructureMoves Moves { get; }
}

/// <summary>
/// BASS's particles for SMC (stage 8; D05), the sampler of the sequence's b = 0 end. A draw is an exact draw of the prior
/// on F, τ from its prior in the g mode (<see cref="BassPriorSampler"/>), loaded into the particle's own design. A fork
/// copies the design as the ladder's hand-over does, the factor and the fit copied and the columns shared, and gives the
/// fork its own copy of the term array, so no two particles share a structure's identity (<see cref="BassLadder"/>'s
/// clone-safe replicas). Draws may run concurrently, each on its own sampler; <see cref="PriorAttempts"/> counts the
/// collections they drew, a sum that does not depend on their order.
/// </summary>
public sealed class BassParticles : ISmcParticles<BassParticle>
{
    private readonly BassModel _model;
    private long _attempts;

    /// <summary>Particles of <paramref name="model"/>, their prior attempts counted from <paramref name="priorAttempts"/>
    /// (a resumed run's).</summary>
    public BassParticles(BassModel model, long priorAttempts = 0)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentOutOfRangeException.ThrowIfNegative(priorAttempts);
        _model = model;
        _attempts = priorAttempts;
    }

    /// <summary>The collections the draws took from the prior, accepted ones included: the attempts per particle are the
    /// prior mass F holds, inverted.</summary>
    public long PriorAttempts => Interlocked.Read(ref _attempts);

    /// <inheritdoc/>
    public BassParticle Draw(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        var sampler = new BassPriorSampler(_model);
        BassState start = sampler.DrawState(rng).State;
        Interlocked.Add(ref _attempts, sampler.Attempts);
        var moves = new BassStructureMoves(_model, start.Terms, start.LogTau);
        return new BassParticle(moves, moves.Start);
    }

    /// <inheritdoc/>
    public BassParticle Fork(BassParticle particle)
    {
        ArgumentNullException.ThrowIfNull(particle);
        BassState state = particle.State;
        BassDesignSnapshot design = particle.Moves.DesignFor(state);
        var copy = new BassState((BassTerm[])state.TermArray.Clone(), state.Fit, state.LogTau);
        var moves = new BassStructureMoves(_model, copy, new BassMovesSnapshot(design, new Tally<BassRefusal>(), 0, new BassSupportTally(_model.Prior),
                                                                               new BassMoveTally()));
        return new BassParticle(moves, copy);
    }
}

/// <summary>The model's collapsed target (<see cref="BassTarget"/>) read through a particle's state: the target an SMC run
/// tempers along D05's likelihood path, in the target's form.</summary>
public sealed class BassParticleTarget : ISamplingTarget<BassParticle>
{
    /// <summary>The target <paramref name="target"/> read through particles.</summary>
    public BassParticleTarget(BassTarget target)
    {
        ArgumentNullException.ThrowIfNull(target);
        Target = target;
    }

    /// <summary>The collapsed target.</summary>
    public BassTarget Target { get; }

    /// <inheritdoc/>
    public ITargetEvaluation Evaluate(BassParticle particle)
    {
        ArgumentNullException.ThrowIfNull(particle);
        return Target.Evaluate(particle.State);
    }
}

/// <summary>
/// BASS's SMC mutation (stage 8; D04, D05): each move is one sweep of D04's chain at the stage's b, the structure move on
/// D01's schedule, then τ's slice step in the g mode. The structure move accepts on ℓ_b(new) − ℓ_b(old) plus the
/// untempered prior and proposal ratios, and τ's slice targets p(τ)·p_b, so every sweep leaves the stage's π_b invariant.
/// <list type="bullet">
/// <item><b>Streams.</b> A sweep runs on the stream SMC hands its particle, through a reversible-jump chain over the
/// particle's own moves, so it reads no other particle's state.</item>
/// <item><b>No adaptation</b> (D10's policy): the kernel reads nothing from the population, so the default
/// <see cref="ISmcKernel{TState}.Prepare"/> stands.</item>
/// <item><b>The account.</b> After each sweep the particle's tallies, its accepted candidate committed first, join the
/// kernel's <see cref="Account"/> and restart, so a particle that resampling drops loses none of them. The account holds
/// integer counts, whose sums do not depend on the order particles move in.</item>
/// </list>
/// Its moves count as accepted when the structure move was; τ's step always moves.
/// </summary>
public sealed class BassSmcKernel : ISmcKernel<BassParticle>
{
    private readonly BassModel _model;
    private readonly BassTauBlock? _tau;
    private readonly CancellationToken _cancellation;
    private readonly Lock _gate = new();

    /// <summary>The kernel of <paramref name="model"/>, its account continued from <paramref name="account"/> (a resumed
    /// run's) or empty; <paramref name="cancellation"/> stops it before a sweep.</summary>
    public BassSmcKernel(BassModel model, BassSmcAccount? account = null, CancellationToken cancellation = default)
    {
        ArgumentNullException.ThrowIfNull(model);
        _model = model;
        _tau = model.Target.SamplesTau ? new BassTauBlock(model.Target) : null;
        _cancellation = cancellation;
        Account = account?.Copy() ?? new BassSmcAccount(model.Prior);
    }

    /// <summary>The sweeps' account so far; read it between stages, not while particles move.</summary>
    public BassSmcAccount Account { get; }

    /// <inheritdoc/>
    public SmcMove<BassParticle> Move(BassParticle particle, ITargetEvaluation evaluation, SmcStage<BassParticle> stage, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(particle);
        ArgumentNullException.ThrowIfNull(rng);
        if (stage.Target is not BassParticleTarget { Target: var target } || !ReferenceEquals(target, _model.Target))
            throw new ArgumentException("The stage's target is not this kernel's model's.", nameof(stage));
        _cancellation.ThrowIfCancellationRequested();
        BassStructureMoves moves = particle.Moves;
        var chain = new ReversibleJumpChain<BassState>(moves.Moves, target, particle.State, rng, stage.Temperature);
        chain.Step();
        bool accepted = chain.Accepted > 0;
        _tau?.Update(chain, rng);
        particle.State = chain.Current;
        moves.DesignAt(particle.State);   // commits an accepted candidate, which tallies its acceptance
        lock (_gate) Account.Add(chain.MoveStats(), moves);
        moves.ResetTallies();
        return new SmcMove<BassParticle>(particle, chain.CurrentEvaluation, accepted);
    }
}

/// <summary>
/// An SMC kernel's account over its run (<see cref="BassSmcKernel"/>): attempts and acceptances per move type, and the
/// moves' tallies by kind, refusals by reason, empty births and support by order, summed over every particle's sweeps.
/// </summary>
public sealed class BassSmcAccount
{
    private readonly List<MoveStat> _moves = [];

    /// <summary>An empty account over <paramref name="prior"/>'s parts and orders.</summary>
    public BassSmcAccount(BassStructurePrior prior)
    {
        ArgumentNullException.ThrowIfNull(prior);
        Support = new BassSupportTally(prior);
    }

    private BassSmcAccount(BassSmcAccount source)
    {
        _moves.AddRange(source._moves);
        Kinds = source.Kinds.Copy();
        Refusals = source.Refusals.Copy();
        EmptyBirths = source.EmptyBirths;
        Support = source.Support.Copy();
        Sweeps = source.Sweeps;
    }

    /// <summary>Attempts and acceptances per move type, in the palette's order; empty before the first sweep.</summary>
    public IReadOnlyList<MoveStat> Moves => _moves;

    /// <summary>The moves by kind (<see cref="BassMoveTally"/>).</summary>
    public BassMoveTally Kinds { get; } = new();

    /// <summary>Refused candidates by reason.</summary>
    public Tally<BassRefusal> Refusals { get; } = new();

    /// <summary>Births whose draw came out the intercept: no move.</summary>
    public long EmptyBirths { get; private set; }

    /// <summary>Candidate terms' support outcomes by their orders.</summary>
    public BassSupportTally Support { get; }

    /// <summary>The sweeps counted.</summary>
    public long Sweeps { get; private set; }

    /// <summary>The structure moves' acceptance over their attempts; NaN before the first attempt.</summary>
    public double AcceptanceRate
    {
        get
        {
            long attempts = 0, accepted = 0;
            foreach (MoveStat m in _moves)
            {
                attempts += m.Attempts;
                accepted += m.Accepted;
            }
            return attempts > 0 ? (double)accepted / attempts : double.NaN;
        }
    }

    /// <summary>An independent copy.</summary>
    public BassSmcAccount Copy() => new(this);

    /// <summary>Adds <paramref name="other"/>'s counts, which must be over the same palette.</summary>
    public void MergeFrom(BassSmcAccount other)
    {
        ArgumentNullException.ThrowIfNull(other);
        AddMoves(other._moves);
        Kinds.MergeFrom(other.Kinds);
        Refusals.MergeFrom(other.Refusals);
        EmptyBirths += other.EmptyBirths;
        Support.MergeFrom(other.Support);
        Sweeps += other.Sweeps;
    }

    // One sweep: its chain's move counts, and the moves' tallies since they last restarted.
    internal void Add(IReadOnlyList<MoveStat> moves, BassStructureMoves tallies)
    {
        AddMoves(moves);
        Kinds.MergeFrom(tallies.Kinds);
        Refusals.MergeFrom(tallies.Refusals);
        EmptyBirths += tallies.EmptyBirths;
        Support.MergeFrom(tallies.Support);
        Sweeps++;
    }

    private void AddMoves(IReadOnlyList<MoveStat> moves)
    {
        if (moves.Count == 0) return;
        if (_moves.Count == 0)
            foreach (MoveStat m in moves) _moves.Add(new MoveStat(m.Key, 0, 0));
        if (_moves.Count != moves.Count) throw new ArgumentException("The sweeps ran different palettes.", nameof(moves));
        for (int i = 0; i < moves.Count; i++)
        {
            if (moves[i].Key != _moves[i].Key) throw new ArgumentException("The sweeps ran different palettes.", nameof(moves));
            _moves[i] = new MoveStat(_moves[i].Key, _moves[i].Attempts + moves[i].Attempts, _moves[i].Accepted + moves[i].Accepted);
        }
    }
}
