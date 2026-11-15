using System;
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
/// </list>
/// A drawn direction whose knot range has no width proposes nothing: the birth is a failed one and the change a null move, a proposal
/// outside the support turned down. Immutable, with no scratch state, so one serves every chain of a fit.
/// </summary>
public sealed class DirectionalTermProposal : IBassProposalKernel
{
    private readonly BassStructurePrior.Part _part;
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
        int[] variables = BassTermProposal.DrawVariablesOf(rng, VariableWeight, _part, counts.AppearanceCounts(BassPart.Directional), active);
        ProjectionDirection direction = _space.DrawDirection(rng, variables);
        return _space.KnotRangeOf(direction) is { } range
            ? new BassTerm(directional: new DirectionalHinge(direction, DirectionalTermSpace.DrawKnot(rng, range)))
            : null;
    }

    /// <inheritdoc/>
    public double LogBirth(BassTerm term, BassTermCounts counts) => LogBirth(term, counts, 0);

    /// <inheritdoc/>
    public double LogBirthRatio(BassTerm born, BassTermCounts counts) => -Math.Log(counts.Terms + 1.0) - LogBirth(born, counts, 0);

    /// <inheritdoc/>
    public double LogDeathRatio(BassTerm dying, BassTermCounts counts) => LogBirth(dying, counts, 1) + Math.Log(counts.Terms);

    /// <summary>
    /// A change of <paramref name="term"/>: θ_𝒥 from PS(θ_𝒥, κ), or a fair sign at a = 1, then the knot from its prior given the new direction.
    /// The within-move ratio is ln π(t | θ) − ln π(t* | θ*), which the prior's ratio cancels, so the move accepts on the likelihood ratio.
    /// </summary>
    public BassTermMove ProposeChange(Xoshiro256PlusPlus rng, BassTerm term, BassTermCounts counts)
    {
        ArgumentNullException.ThrowIfNull(rng);
        Check(counts);
        Prior.Validate(term);
        DirectionalHinge hinge = term.Directional!;
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
        if (_space.KnotRangeOf(to) is not { } range) return new BassTermMove(BassMoveKind.DirectionChange, null, double.NaN);
        var moved = new DirectionalHinge(to, DirectionalTermSpace.DrawKnot(rng, range));
        return new BassTermMove(BassMoveKind.DirectionChange, new BassTerm(directional: moved), _space.LogKnot(hinge) - _space.LogKnot(moved));
    }

    // ln q_b(term) under the counts less `shift` copies of the term: 0 for a birth into them, 1 for the reverse of a death from them. The
    // order's and the variables' laws are D01's; θ's and the knot's are their priors.
    private double LogBirth(BassTerm term, BassTermCounts counts, int shift)
    {
        Check(counts);
        Prior.Validate(term);
        int terms = counts.Terms - shift;
        if (terms < 0) throw new ArgumentException("The counts hold no such term.", nameof(term));
        int active = term.Directional!.Direction.Count;
        int count = counts.Orders(BassPart.Directional)[active] - shift;
        if (count < 0) throw new ArgumentException("The counts hold no such term.", nameof(term));
        return BassTermProposal.LogOrderOf(OrderWeight, _part, count, terms)
             + BassTermProposal.LogVariablesOf(VariableWeight, term, BassPart.Directional, _part, counts.AppearanceCounts(BassPart.Directional), shift)
             + _part.LogFactors(term);
    }

    private void Check(BassTermCounts counts)
    {
        ArgumentNullException.ThrowIfNull(counts);
        if (!ReferenceEquals(counts.Prior, Prior)) throw new ArgumentException("The counts are over another prior's term space.", nameof(counts));
    }
}
