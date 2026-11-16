using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// E36's structure proposals over a directional basis (slice 2; decisions 4 and 5), the paper's moves with the corrections slice 0 read in
/// <c>BayesPPR</c> and Lean's ratios (<c>DirectionMoves.lean</c>):
/// <list type="bullet">
/// <item><b>The birth.</b> a by D01's order law, I(a) ∝ w₁ + (terms of active count a); 𝒥 by D01's variable law, uniform at a = 1 and
/// Wallenius' with z_u ∝ w₂ + (terms using u) otherwise, which is the paper's law corrected (its eq. (mhb) divides by Wallenius' weight at
/// a = 1, where its Algorithm 2 draws 𝒥 uniformly). θ and the knot are drawn from their priors given a, 𝒥 and θ, so their densities cancel against
/// the prior's and the birth accepts on L′/L · λ-part · π(a)π(𝒥|a)/(q(a)q(𝒥|a)) (<c>birth_green_ratio</c>, <c>prior_draw_cancel</c>).</item>
/// <item><b>The death.</b> A term uniformly; its reverse is the dying term's birth under the counts without it, the supplement's death ratio
/// corrected to the reduced counts, as <c>BayesPPR</c> takes them.</item>
/// <item><b>The change.</b> a and 𝒥 kept; θ_𝒥 drawn from the power spherical law centred on the current θ_𝒥 with concentration κ, or at
/// a = 1 its sign from a fair coin (the user's decision, 2026-11-14: <c>BayesPPR</c>'s form, since the law on S⁰ keeps the sign); then the
/// knot from its prior given the proposed direction, whose range moves with it. The law is symmetric, so with the knot's prior density
/// cancelling the move accepts on the collapsed likelihood ratio alone (<c>change_ratio_powerSpherical</c>).</item>
/// <item><b>Categorical gates</b> (slice 5, decision 7): beside a categorical part a birth draws both parts' orders by D01's law, each from 0,
/// the categorical variables by D01's variable law and their subsets from the prior; both empty is an empty birth. The knot's range is
/// read over the term's cells. A change on a term holding both parts takes the direction or a subset with probability ½ each: the direction's
/// as above, the knot from its prior over the cells; a subset's D01's redraw from the subset measure, then the knot from its prior over the
/// new cells, so both accept on the likelihood ratio alone (the subset measure's and the knot's ratios cancel the prior's).</item>
/// </list>
/// A drawn direction or subset whose knot range has none (no width, or cells too small or empty) proposes nothing: the birth is a failed one
/// and the change a null move, a proposal outside the support turned down. Immutable, with no scratch state, so one serves every chain of a
/// fit.
/// </summary>
public sealed class DirectionalTermProposal : IBassProposalKernel
{
    private readonly BassStructurePrior.Part _part;
    private readonly BassStructurePrior.Part? _gate;   // the categorical part, beside the direction (slice 5); null without one
    private readonly DirectionalTermSpace _space;

    /// <summary>Proposals over <paramref name="prior"/>'s directional term space with D01's pseudo-counts w₁ and w₂, both positive, and the
    /// change's concentration κ ≥ 0.</summary>
    public DirectionalTermProposal(BassStructurePrior prior, double orderWeight, double variableWeight, double concentration)
    {
        ArgumentNullException.ThrowIfNull(prior);
        if (prior.Directional is not { } space || prior[BassPart.Directional] is not { } part)
            throw new ArgumentException("The prior has no directional part.", nameof(prior));
        if (!(orderWeight > 0.0 && double.IsFinite(orderWeight)))
            throw new ArgumentOutOfRangeException(nameof(orderWeight), "w₁ must be positive and finite, or a death could have no reverse birth.");
        if (!(variableWeight > 0.0 && double.IsFinite(variableWeight)))
            throw new ArgumentOutOfRangeException(nameof(variableWeight), "w₂ must be positive and finite, or a death could have no reverse birth.");
        if (!(concentration >= 0.0 && double.IsFinite(concentration)))
            throw new ArgumentOutOfRangeException(nameof(concentration), "κ must be finite and not negative.");
        Prior = prior;
        _part = part;
        _gate = prior[BassPart.Categorical];
        _space = space;
        OrderWeight = orderWeight;
        VariableWeight = variableWeight;
        Concentration = concentration;
    }

    /// <summary>The proposals <paramref name="config"/> sets over <paramref name="prior"/>'s directional term space.</summary>
    public static DirectionalTermProposal From(BassConfig config, BassStructurePrior prior)
    {
        ArgumentNullException.ThrowIfNull(config);
        DirectionalBasis basis = config.Structure.Directional ?? throw new ArgumentException("The config declares no directional basis.", nameof(config));
        return new(prior, config.Structure.OrderWeight, config.Structure.VariableWeight, basis.Concentration);
    }

    /// <inheritdoc/>
    public BassStructurePrior Prior { get; }

    /// <summary>w₁, the pseudo-count on active counts, the paper's ω₀.</summary>
    public double OrderWeight { get; }

    /// <summary>w₂, the pseudo-count on variables, the paper's υ₀.</summary>
    public double VariableWeight { get; }

    /// <summary>κ, the change's power spherical concentration.</summary>
    public double Concentration { get; }

    /// <inheritdoc/>
    public BassTerm? DrawBirth(Xoshiro256PlusPlus rng, BassTermCounts counts)
    {
        ArgumentNullException.ThrowIfNull(rng);
        Check(counts);
        int active = BassTermProposal.DrawOrderOf(rng, OrderWeight, _part, counts.Orders(BassPart.Directional));
        int k = _gate is null ? 0 : BassTermProposal.DrawOrderOf(rng, OrderWeight, _gate, counts.Orders(BassPart.Categorical));
        if (active == 0 && k == 0) return null;
        ProjectionDirection? direction = null;
        if (active > 0)
        {
            int[] variables = BassTermProposal.DrawVariablesOf(rng, VariableWeight, _part, counts.AppearanceCounts(BassPart.Directional), active);
            direction = _space.DrawDirection(rng, variables);
        }
        var cells = new CategoricalFactor[k];
        if (k > 0)
        {
            int[] variables = BassTermProposal.DrawVariablesOf(rng, VariableWeight, _gate!, counts.AppearanceCounts(BassPart.Categorical), k);
            for (int f = 0; f < k; f++) cells[f] = _gate!.DrawSubset(rng, variables[f], Prior.Subsets);
        }
        return BassStructurePrior.WithKnot(rng, _space, direction, cells);
    }

    /// <inheritdoc/>
    public double LogBirth(BassTerm term, BassTermCounts counts) => LogBirth(term, counts, 0);

    /// <inheritdoc/>
    public double LogBirthRatio(BassTerm born, BassTermCounts counts) => -Math.Log(counts.Terms + 1.0) - LogBirth(born, counts, 0);

    /// <inheritdoc/>
    public double LogDeathRatio(BassTerm dying, BassTermCounts counts) => LogBirth(dying, counts, 1) + Math.Log(counts.Terms);

    /// <summary>
    /// A change of <paramref name="term"/>: θ_𝒥 from PS(θ_𝒥, κ), or a fair sign at a = 1, then the knot from its prior given the new direction
    /// and the term's cells; or, beside a categorical part, one of its subsets redrawn (<see cref="ProposeSubset"/>), each with probability ½
    /// where the term holds both. The within-move ratio is the knot's prior ratio, and a subset's measure's, which the prior's ratio cancels,
    /// so the move accepts on the likelihood ratio.
    /// </summary>
    public BassTermMove ProposeChange(Xoshiro256PlusPlus rng, BassTerm term, BassTermCounts counts)
    {
        ArgumentNullException.ThrowIfNull(rng);
        Check(counts);
        Prior.Validate(term);
        if (term.Directional is not { } hinge || (term.Categorical.Count > 0 && rng.NextDouble() < 0.5)) return ProposeSubset(rng, term);
        IReadOnlyList<CategoricalFactor> cells = term.Categorical;
        ProjectionDirection from = hinge.Direction;
        ProjectionDirection to;
        if (from.Count == 1) to = new ProjectionDirection([from.Variables[0]], [BassDraws.Sign(rng)]);
        else
        {
            var weights = new double[from.Count];
            PowerSpherical.Of([.. from.Weights], Concentration).Sample(rng, weights);
            foreach (double w in weights)
                if (w == 0.0) return new BassTermMove(BassMoveKind.DirectionChange, null, double.NaN);   // off the open subsphere
            to = new ProjectionDirection(from.Variables, weights);
        }
        if (_space.KnotRangeOf(to, cells) is not { } range) return new BassTermMove(BassMoveKind.DirectionChange, null, double.NaN);
        var moved = new DirectionalHinge(to, DirectionalTermSpace.DrawKnot(rng, range));
        return new BassTermMove(BassMoveKind.DirectionChange, new BassTerm(categorical: cells, directional: moved),
                                _space.LogKnot(hinge, cells) - _space.LogKnot(moved, cells));
    }

    /// <summary>
    /// A subset of <paramref name="term"/>'s categorical part redrawn (slice 5): one factor uniformly, its subset from the prior's measure,
    /// D01's redraw; then, where the term bears a direction, the knot from its prior over the new cells, a null move where they give no range.
    /// The within-move ratio is ln π_v(D) − ln π_v(D′) + ln π(t | θ, cells) − ln π(t* | θ, cells′), which the prior's ratio cancels.
    /// </summary>
    public BassTermMove ProposeSubset(Xoshiro256PlusPlus rng, BassTerm term)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (_gate is null || term.Categorical.Count == 0) throw new ArgumentException("The term holds no categorical factor.", nameof(term));
        CategoricalFactor[] cells = [.. term.Categorical];
        int pick = BassDraws.Below(rng, cells.Length);
        CategoricalFactor old = cells[pick];
        cells[pick] = _gate.DrawSubset(rng, old.Variable, Prior.Subsets);
        double[] log = _gate.LogSubset![old.Variable];
        double ratio = log[old.Levels.Count] - log[cells[pick].Levels.Count];
        if (term.Directional is not { } hinge) return new BassTermMove(BassMoveKind.Redraw, new BassTerm(categorical: cells), ratio);
        if (_space.KnotRangeOf(hinge.Direction, cells) is not { } range) return new BassTermMove(BassMoveKind.Redraw, null, double.NaN);
        var moved = new DirectionalHinge(hinge.Direction, DirectionalTermSpace.DrawKnot(rng, range));
        return new BassTermMove(BassMoveKind.Redraw, new BassTerm(categorical: cells, directional: moved),
                                ratio + _space.LogKnot(hinge, term.Categorical) - _space.LogKnot(moved, cells));
    }

    // ln q_b(term) under the counts less `shift` copies of the term: 0 for a birth into them, 1 for the reverse of a death from them. The
    // order's and the variables' laws are D01's; θ's and the knot's are their priors.
    private double LogBirth(BassTerm term, BassTermCounts counts, int shift)
    {
        Check(counts);
        Prior.Validate(term);
        int terms = counts.Terms - shift;
        if (terms < 0) throw new ArgumentException("The counts hold no such term.", nameof(term));
        double lq = 0.0;
        foreach ((BassPart kind, BassStructurePrior.Part? part) in new[] { (BassPart.Directional, (BassStructurePrior.Part?)_part), (BassPart.Categorical, _gate) })
        {
            if (part is null) continue;
            int k = BassStructurePrior.Order(term, kind);
            int count = counts.Orders(kind)[k] - shift;
            if (count < 0) throw new ArgumentException("The counts hold no such term.", nameof(term));
            lq += BassTermProposal.LogOrderOf(OrderWeight, part, count, terms);
            if (k == 0) continue;
            lq += BassTermProposal.LogVariablesOf(VariableWeight, term, kind, part, counts.AppearanceCounts(kind), shift);
            lq += part.LogFactors(term);
        }
        return lq;
    }

    private void Check(BassTermCounts counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        if (!ReferenceEquals(counts.Prior, Prior)) throw new ArgumentException("The counts are over another prior's term space.", nameof(counts));
    }
}
