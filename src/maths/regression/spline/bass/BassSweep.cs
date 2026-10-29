using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// One chain's structure moves (D01, D04): its own design and counts, and the birth, death and change moves over them,
/// as the reversible-jump palette a <see cref="ReversibleJumpChain{TState}"/> runs on D01's move schedule.
/// <list type="bullet">
/// <item><b>The schedule.</b> Equal over the moves available at M (birth only at 0; death and change at M_max; all three
/// between), or the prior-cancelling schedule b_M = C·min(1, p(M+1)/p(M)), d_M = C·min(1, p(M−1)/p(M)), the rest change,
/// so b_M p(M) = d_{M+1} p(M+1) and the count prior cancels from births and deaths. The chain owns the selection ratio,
/// so these are the moves' weights.</item>
/// <item><b>The proposals.</b> A birth draws its term from D01's proposal (an empty draw is no move); a death and a change
/// pick a term uniformly. Each candidate carries its fit from the design, or the design's refusal, which the target
/// scores off the support, and its within-move ratio from <see cref="BassTermProposal"/>.</item>
/// <item><b>The design follows the chain.</b> A proposal leaves the design's candidate pending; when the chain next
/// shows a state holding that candidate's structure it was accepted, and the design commits it. A replica handed a
/// state by a swap adopts its giver's design, and a resize's new replica its source's (<see cref="BassLadder"/>); any
/// other structure is loaded whole. A chain restored from its snapshot gets its design back as it was.</item>
/// </list>
/// Refusals are tallied by reason, the target's exact-fit refusal included, and each birth's and change's candidate term
/// by its orders with its support outcome (D01's tilt). Not thread-safe: one per chain.
/// </summary>
public sealed class BassStructureMoves
{
    private enum Kind { Birth, Death, Change }

    private readonly BassModel _model;
    private readonly BassDesign _design;
    private BassTermCounts _counts;
    private BassState _committed;
    private (BassState State, Kind Kind, BassTerm? Added, BassTerm? Removed)? _candidate;

    /// <summary>The moves of a chain started at <paramref name="terms"/>, which must lie in F, with u =
    /// <paramref name="logTau"/> (the held τ's log under fixed g).</summary>
    public BassStructureMoves(BassModel model, IReadOnlyList<BassTerm> terms, double logTau)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(terms);
        _model = model;
        _design = model.NewDesign();
        if (_design.Load(terms) is { } refusal) throw new ArgumentException($"The start lies outside F ({refusal}).", nameof(terms));
        BassTerm[] loaded = [.. _design.Terms];
        Start = new BassState(loaded, Fit(loaded), model.Target.SamplesTau ? logTau : Math.Log(model.Target.FixedTau!.Value));
        if (model.Target.Refusal(Start.Fit) is { } exact) throw new ArgumentException($"The start lies outside F ({exact}).", nameof(terms));
        _counts = CountsOf(loaded);
        _committed = Start;
        Support = new BassSupportTally(model.Prior);
        Moves = [new Birth(this), new Death(this), new Change(this)];
    }

    // The moves restored at `current` from `snapshot` (Snapshot): the design as it was, the counts from the structure,
    // the tallies; they continue bit for bit.
    internal BassStructureMoves(BassModel model, BassState current, BassMovesSnapshot snapshot)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(current);
        ArgumentNullException.ThrowIfNull(snapshot);
        _model = model;
        _design = model.NewDesign();
        _design.Restore(snapshot.Design);
        if (!_design.Terms.SequenceEqual(current.Terms)) throw new ArgumentException("The snapshot's design holds another structure.", nameof(snapshot));
        Start = current;
        _counts = CountsOf(current.TermArray);
        _committed = current;
        Refusals.Assign(snapshot.Refusals);
        EmptyBirths = snapshot.EmptyBirths;
        Support = snapshot.Support.Copy();
        Moves = [new Birth(this), new Death(this), new Change(this)];
    }

    /// <summary>The state these moves started from: the chain's start, or the state they were restored at.</summary>
    public BassState Start { get; }

    /// <summary>The birth, death and change moves, keyed "birth", "death" and "change".</summary>
    public IReadOnlyList<IRjMove<BassState>> Moves { get; }

    /// <summary>Refused candidates by reason: the design's (cap, support, numerical, variance inflation) and the
    /// target's exact fits.</summary>
    public Tally<BassRefusal> Refusals { get; } = new();

    /// <summary>Births whose draw left every part empty, the intercept: no move (D01).</summary>
    public long EmptyBirths { get; private set; }

    /// <summary>Each birth's and change's candidate term by its orders, refused for support or not (D01's tilt).</summary>
    public BassSupportTally Support { get; }

    /// <summary>Clears the tallies (<see cref="Refusals"/>, <see cref="EmptyBirths"/>, <see cref="Support"/>), so they
    /// count from here: a run's sampling phase.</summary>
    public void ResetTallies()
    {
        Refusals.Assign(new Tally<BassRefusal>());
        EmptyBirths = 0;
        Support.Clear();
    }

    // The restart state at `current`, which the design is first synced to: the design's own, and the tallies.
    internal BassMovesSnapshot Snapshot(BassState current)
    {
        Sync(current);
        return new BassMovesSnapshot(_design.Snapshot(), Refusals.Copy(), EmptyBirths, Support.Copy());
    }

    /// <summary>The chain's design, holding <paramref name="state"/>'s structure.</summary>
    public BassDesign DesignAt(BassState state)
    {
        Sync(state);
        return _design;
    }

    // Whether the design can show `state`'s structure without a load: the committed one, or the pending candidate.
    internal bool Holds(BassState state) => state.SharesStructure(_committed) || _candidate is { } c && state.SharesStructure(c.State);

    // An exact copy of the design at `state`, which it is first synced to, for another replica to adopt (BassLadder): the
    // factor and the fit copied, the columns shared.
    internal BassDesignSnapshot DesignFor(BassState state)
    {
        Sync(state);
        return _design.Snapshot(shareColumns: true);
    }

    // Takes over `state`, which a swap handed to these moves, with the design its giver held for it (DesignFor): the
    // design as the giver's, bit for bit, so the structure is neither refactored, re-admitted nor re-evaluated here.
    internal void Adopt(BassState state, BassDesignSnapshot design)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(design);
        if (!design.Terms.AsSpan().SequenceEqual(state.TermArray)) throw new ArgumentException("The design holds another structure.", nameof(design));
        _design.Restore(design);
        _candidate = null;
        _counts = CountsOf(state.TermArray);
        _committed = state;
    }

    /// <summary>
    /// The retained draw at <paramref name="state"/> (D04's composition): σ² from its collapsed conditional, the
    /// coefficients given σ², and λ | M ~ Gamma(h₁ + M, rate h₂ + 1), each exact given (structure, τ). The terms are
    /// interned in <paramref name="terms"/>, and the coefficients move to the uncentred basis in the ids' order.
    /// </summary>
    public BassDraw Compose(BassState state, Xoshiro256PlusPlus rng, BassTermTable terms)
    {
        ArgumentNullException.ThrowIfNull(rng);
        ArgumentNullException.ThrowIfNull(terms);
        Sync(state);
        BassTarget target = _model.Target;
        double tau = target.SamplesTau ? state.Tau : target.FixedTau!.Value;
        CollapsedGPrior block = CollapsedGPrior.FromCentredLeastSquares(_design.Factor, _design.Coefficients.ToArray(), _design.ResidualSquares,
                                                                         _design.FittedSquares, _design.Observations, _design.ResponseMean, target.Intercept);
        var centred = new double[block.Columns];   // (α, β_c) in the design's order
        double sigma2 = block.At(tau).Sample(target.Sigma2, rng, centred);
        BassCountPrior count = _model.Prior.Count;
        double lambda = Gamma.ShapeRate(count.Shape + state.Count, count.Rate + 1.0).Sample(rng);

        int m = state.Count;
        var ids = new int[m];
        for (int j = 0; j < m; j++) ids[j] = terms.Intern(state.Terms[j]);
        int[] order = [.. Enumerable.Range(0, m).OrderBy(j => ids[j])];
        var coefficients = new double[m + 1];
        double intercept = centred[0];
        for (int k = 0; k < m; k++)
        {
            int j = order[k];
            coefficients[k + 1] = centred[j + 1];
            intercept -= centred[j + 1] * ColumnMean(state.Terms[j]);
        }
        coefficients[0] = intercept;
        return new BassDraw(new BassStructure([.. order.Select(j => ids[j])]), tau, sigma2, coefficients, lambda);
    }

    // The mean of the term's column over all N observations, the design's centring: a case column's mean, or āh̄.
    private double ColumnMean(BassTerm term)
    {
        BassBasis basis = _model.Basis;
        var a = new double[basis.Inputs.Cases.Cases];
        basis.CaseColumn(term, basis.Inputs.Cases, a);
        double mean = PairwiseSum.Of(a) / a.Length;
        if (basis.Inputs.Outputs is not { } outputs) return mean;
        var h = new double[outputs.Locations];
        basis.OutputColumn(term, outputs, h);
        return mean * (PairwiseSum.Of(h) / h.Length);
    }

    // The design and counts follow the chain to `state`: a τ update keeps the structure; the pending candidate's
    // structure means the chain accepted it; anything else is loaded whole.
    private void Sync(BassState state)
    {
        ArgumentNullException.ThrowIfNull(state);
        if (state.SharesStructure(_committed))
        {
            _committed = state;
            return;
        }
        if (_candidate is { } c && state.SharesStructure(c.State))
        {
            _design.Accept();
            if (c.Removed is { } removed) _counts.Remove(removed);
            if (c.Added is { } added) _counts.Add(added);
            _candidate = null;
            _committed = state;
            return;
        }
        _candidate = null;
        if (_design.Load(state.Terms) is { } refusal) throw new InvalidOperationException($"The chain holds a structure outside F ({refusal}).");
        if (!_design.Terms.SequenceEqual(state.Terms)) throw new InvalidOperationException("The design reordered a loaded structure.");
        _counts = CountsOf(state.TermArray);
        _committed = state;
    }

    private BassTermCounts CountsOf(IEnumerable<BassTerm> terms)
    {
        var counts = new BassTermCounts(_model.Prior);
        foreach (BassTerm t in terms) counts.Add(t);
        return counts;
    }

    private BassStructureFit Fit(BassTerm[] terms) => BassStructureFit.Of(_design, LogTermPrior(terms));

    private double LogTermPrior(BassTerm[] terms)
    {
        double sum = 0.0;
        foreach (BassTerm t in terms) sum += _model.Prior.LogPrior(t);
        return sum;
    }

    // The candidate for a proposal the design has weighed: admitted, it is pending until the chain shows it; refused,
    // by the design or the target, it is tallied and scored off the support.
    private BassState Candidate(BassState current, BassTerm[] terms, DesignProposal proposal, Kind kind, BassTerm? added, BassTerm? removed)
    {
        BassStructureFit fit = BassStructureFit.Of(proposal, proposal.Admitted ? LogTermPrior(terms) : double.NaN);
        if (added is not null) Support.Add(added, proposal.Refusal == BassRefusal.Support);
        var candidate = new BassState(terms, fit, current.LogTau);
        if (_model.Target.Refusal(fit) is { } reason)
        {
            Refusals.Add(reason);
            _candidate = null;
        }
        else _candidate = (candidate, kind, added, removed);
        return candidate;
    }

    // D01's schedule: each move's weight at M, normalized by the chain over the palette.
    private double Weight(Kind kind, int m)
    {
        int cap = _model.MaxBasis;
        bool birth = m < cap, other = m > 0;
        switch (_model.Config.Sampling.Moves)
        {
            case PriorCancelling { C: var c }:
                double b = birth ? c * Math.Min(1.0, Math.Exp(_model.Prior.LogCountRatio(m))) : 0.0;
                double d = other ? c * Math.Min(1.0, Math.Exp(-_model.Prior.LogCountRatio(m - 1))) : 0.0;
                // The change keeps 1 − b − d even at M = 0, where it proposes nothing, so the weights sum to one and
                // the prior cancels at the boundary too (BARS's DMGK schedule).
                return kind switch { Kind.Birth => b, Kind.Death => d, _ => Math.Max(0.0, 1.0 - b - d) };
            default:
                return kind == Kind.Birth ? (birth ? 1.0 : 0.0) : (other ? 1.0 : 0.0);
        }
    }

    private sealed class Birth(BassStructureMoves owner) : IRjMove<BassState>
    {
        public string Key => "birth";
        public string ReverseKey => "death";
        public double Weight(BassState state) => owner.Weight(Kind.Birth, state.Count);

        public Proposal<BassState>? Propose(BassState current, Xoshiro256PlusPlus rng)
        {
            owner.Sync(current);
            if (owner._model.Proposal.DrawBirth(rng, owner._counts) is not { } term)
            {
                owner.EmptyBirths++;
                owner._candidate = null;
                return null;
            }
            double logRatio = owner._model.Proposal.LogBirthRatio(term, owner._counts);
            DesignProposal proposal = owner._design.ProposeBirth(term);
            BassTerm[] terms = [.. current.TermArray, term];
            return new Proposal<BassState>(owner.Candidate(current, terms, proposal, Kind.Birth, term, null), logRatio, 0.0);
        }
    }

    private sealed class Death(BassStructureMoves owner) : IRjMove<BassState>
    {
        public string Key => "death";
        public string ReverseKey => "birth";
        public double Weight(BassState state) => owner.Weight(Kind.Death, state.Count);

        public Proposal<BassState>? Propose(BassState current, Xoshiro256PlusPlus rng)
        {
            owner.Sync(current);
            if (current.Count == 0) return null;
            int index = BassDraws.Below(rng, current.Count);
            BassTerm term = current.TermArray[index];
            double logRatio = owner._model.Proposal.LogDeathRatio(term, owner._counts);
            DesignProposal proposal = owner._design.ProposeDeath(index);
            BassTerm[] terms = [.. current.TermArray[..index], .. current.TermArray[(index + 1)..]];
            return new Proposal<BassState>(owner.Candidate(current, terms, proposal, Kind.Death, null, term), logRatio, 0.0);
        }
    }

    private sealed class Change(BassStructureMoves owner) : IRjMove<BassState>
    {
        public string Key => "change";
        public string ReverseKey => "change";
        public double Weight(BassState state) => owner.Weight(Kind.Change, state.Count);

        public Proposal<BassState>? Propose(BassState current, Xoshiro256PlusPlus rng)
        {
            owner.Sync(current);
            if (current.Count == 0) return null;
            int index = BassDraws.Below(rng, current.Count);
            BassTerm old = current.TermArray[index];
            BassTerm term = owner._model.Proposal.DrawChange(rng, old, out double logRatio);
            owner._candidate = null;
            if (term.Equals(old))   // the same state, whose term has its support
            {
                owner.Support.Add(term, false);
                return new Proposal<BassState>(current, 0.0, 0.0);
            }
            DesignProposal proposal = owner._design.ProposeChange(index, term);
            // The design moves a changed term to the end, so the candidate's order follows it.
            BassTerm[] terms = [.. current.TermArray[..index], .. current.TermArray[(index + 1)..], term];
            return new Proposal<BassState>(owner.Candidate(current, terms, proposal, Kind.Change, term, old), logRatio, 0.0);
        }
    }
}

/// <summary>
/// τ's block (D04): a slice step on u = ln τ (<see cref="SliceSampler"/>) against the chain's target at the chain's own
/// temperature, each evaluation O(1) from the structure's cached fit. The target's density is already on u, so there is
/// no Jacobian to add. A chain block, so a lone chain runs it beside its moves in a <see cref="BlockScan{TState}"/> and a
/// tempering replica between its moves and the swaps.
/// </summary>
public sealed class BassTauBlock : IBlock<ReversibleJumpChain<BassState>>
{
    /// <summary>The slice's initial width in u.</summary>
    public const double Width = 1.0;

    /// <summary>The stepping-out limit.</summary>
    public const int MaxSteps = 32;

    private readonly ISamplingTarget<BassState> _target;

    /// <summary>τ's block on the chain's target <paramref name="target"/>.</summary>
    public BassTauBlock(ISamplingTarget<BassState> target)
    {
        ArgumentNullException.ThrowIfNull(target);
        _target = target;
    }

    /// <inheritdoc/>
    public string Name => "tau";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        BassState current = state.Current;
        Temperature temperature = state.Temperature;
        double u = current.LogTau, logDensity = state.CurrentLogDensity;
        SliceSampler.Step(ref u, ref logDensity, Width, MaxSteps,
                          v => temperature.LogDensity(_target.Evaluate(current.WithLogTau(v))), rng);
        BassState next = current.WithLogTau(u);
        state.SetState(next, _target.Evaluate(next));
        return true;
    }
}

/// <summary>
/// One untempered chain (D04's sweep): each step is a structure move on D01's schedule, then τ's slice step in the g
/// mode, run as a <see cref="BlockScan{TState}"/> over the reversible-jump chain. The moves draw from the chain's stream
/// and τ's block from the scan's. <see cref="Compose"/> gives a retained draw at the current state.
/// </summary>
public sealed class BassChain
{
    /// <summary>A chain of <paramref name="model"/> from <paramref name="start"/> and u = <paramref name="logTau"/>,
    /// sampling <paramref name="target"/> (the model's by default).</summary>
    public BassChain(BassModel model, IReadOnlyList<BassTerm> start, double logTau, Xoshiro256PlusPlus moveRng,
                     Xoshiro256PlusPlus blockRng, ISamplingTarget<BassState>? target = null)
    {
        ArgumentNullException.ThrowIfNull(model);
        Structure = new BassStructureMoves(model, start, logTau);
        ISamplingTarget<BassState> sampled = target ?? model.Target;
        Moves = new ReversibleJumpChain<BassState>(Structure.Moves, sampled, Structure.Start, moveRng);
        Scan = new BlockScan<ReversibleJumpChain<BassState>>(Blocks(model, sampled), Moves, blockRng, chain => chain);
    }

    private BassChain(BassStructureMoves structure, ReversibleJumpChain<BassState> moves, BlockScan<ReversibleJumpChain<BassState>> scan)
    {
        Structure = structure;
        Moves = moves;
        Scan = scan;
    }

    // The scan's blocks: the moves, then τ in the g mode.
    private static List<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassModel model, ISamplingTarget<BassState> target)
    {
        List<IBlock<ReversibleJumpChain<BassState>>> blocks = [new MoveStep<BassState>()];
        if (model.Target.SamplesTau) blocks.Add(new BassTauBlock(target));
        return blocks;
    }

    /// <summary>
    /// The chain's restart state, the chain unchanged: the reversible-jump chain's snapshot (state, temperature, move
    /// stream, counters), the moves' (the design as it is, the tallies) and the scan's (τ's stream, step count and block
    /// counters).
    /// </summary>
    internal BassChainSnapshot Snapshot()
    {
        BlockScanSnapshot<ReversibleJumpChain<BassState>> scan = Scan.Snapshot();
        return new BassChainSnapshot(Moves.Snapshot(), Structure.Snapshot(Current), (scan.Rng0, scan.Rng1, scan.Rng2, scan.Rng3),
                                     scan.Steps, scan.BlockAttempts, scan.BlockAccepts);
    }

    /// <summary>A chain of <paramref name="model"/> that continues bit for bit where <paramref name="snapshot"/>'s left
    /// off, sampling <paramref name="target"/> (the model's by default), which must be the one it ran with.</summary>
    internal static BassChain Restore(BassModel model, BassChainSnapshot snapshot, ISamplingTarget<BassState>? target = null)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(snapshot);
        var structure = new BassStructureMoves(model, snapshot.Chain.Current, snapshot.Structure);
        ISamplingTarget<BassState> sampled = target ?? model.Target;
        var moves = ReversibleJumpChain<BassState>.Restore(structure.Moves, sampled, snapshot.Chain);
        (ulong r0, ulong r1, ulong r2, ulong r3) = snapshot.BlockRng;
        var scan = BlockScan<ReversibleJumpChain<BassState>>.Restore(Blocks(model, sampled),
            new BlockScanSnapshot<ReversibleJumpChain<BassState>>(moves, r0, r1, r2, r3, snapshot.Steps, snapshot.BlockAttempts,
                                                                  snapshot.BlockAccepts), chain => chain);
        return new BassChain(structure, moves, scan);
    }

    /// <summary>The chain's structure moves, design and tallies.</summary>
    public BassStructureMoves Structure { get; }

    /// <summary>The reversible-jump chain.</summary>
    public ReversibleJumpChain<BassState> Moves { get; }

    /// <summary>The scan: the moves, then τ.</summary>
    public BlockScan<ReversibleJumpChain<BassState>> Scan { get; }

    /// <summary>The current state.</summary>
    public BassState Current => Moves.Current;

    /// <summary>One sweep; the state after it.</summary>
    public BassState Step()
    {
        Scan.Step();
        return Moves.Current;
    }

    /// <summary>The retained draw at the current state (<see cref="BassStructureMoves.Compose"/>).</summary>
    public BassDraw Compose(Xoshiro256PlusPlus rng, BassTermTable terms) => Structure.Compose(Current, rng, terms);
}

/// <summary>A chain's structure moves' restart state: the design's, and the tallies.</summary>
internal sealed record BassMovesSnapshot(BassDesignSnapshot Design, Tally<BassRefusal> Refusals, long EmptyBirths, BassSupportTally Support);

/// <summary>A <see cref="BassChain"/>'s restart state: the reversible-jump chain's snapshot, the moves', and the scan's
/// stream (τ's block), step count and block counters.</summary>
internal sealed record BassChainSnapshot(
    ChainSnapshot<BassState> Chain, BassMovesSnapshot Structure, (ulong, ulong, ulong, ulong) BlockRng, long Steps,
    long[] BlockAttempts, long[] BlockAccepts);
