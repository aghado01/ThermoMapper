using System;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

// Random effects under GBASS's noises (E21's slice 4; decision 7 as the user amended it on 2026-11-10). u is a latent on every rung,
// which the reading takes as the offset Z̃u beside σβv, with σ² on the record: each rung's kernel then tempers the likelihood given u
// and the latents, as the noise's blocks do, so ladders, SMC and the reference rung are exact. At b = 1 the structure step may integrate
// u and draw it afresh before any block reads it, the partially collapsed order of the β·v move's fix: an opt-in since 2026-11-10
// (BassSampling.CollapseEffectsInMoves), off by default, on slice 4's measurement.

/// <summary>
/// u's block (the coefficients' slot, after a): u from its conditional given a, σ², the latents and the scales, at the rung's power b,
/// N(K·Z̃ᵀ(bW)e, σ²K) in contrasts with K = (S⁻² + Z̃ᵀ(bW)Z̃)⁻¹ and e = y − σβv − f, through the weighted whitening at bW
/// (<see cref="WeightedWhitening"/>; at b = 0 its prior). The chain follows the offset at the new u, and the sweep's fit is handed on with
/// Z̃u in it, so the blocks after read the residual net of u and σ²'s and interweaving's priors its energy.
/// </summary>
public sealed class BassGroupEffectsBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private readonly BassModel _model;
    private readonly VarianceMeanMixture _mixture;

    internal BassGroupEffectsBlock(BassModel model)
    {
        ArgumentNullException.ThrowIfNull(model);
        _model = model;
        _mixture = model.Observation as VarianceMeanMixture ?? throw new ArgumentException("u is a latent under a noise.", nameof(model));
    }

    /// <inheritdoc/>
    public string Name => "group effects";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        BassStructureMoves moves = BassStructureMoves.Of(state);
        BassSweepFit fit = moves.SweepFit is { } drawn && ReferenceEquals(drawn.State, state.Current)
            ? drawn
            : throw new InvalidOperationException("u reads this sweep's coefficients: the coefficients' block runs first, at the current state.");
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        RandomEffectsLatent current = working.Effects ?? throw new InvalidOperationException("The reading carries no random effects.");
        ObservationLatents latents = working.Latents ?? throw new InvalidOperationException("The working response holds no latents.");
        double sigma2 = latents.Sigma2 ?? fit.Sigma2, power = BassCoefficientsBlock.Power(state);
        ReadOnlySpan<double> values = working.Values, weights = working.Weights;
        int n = values.Length;
        var old = new double[n];
        current.AddTo(old);
        // f, the coefficients' part of the fit; e = y − σβv − f, the reading's values with the old Z̃u put back.
        var terms = new double[n];
        var residual = new double[n];
        var tempered = new double[n];
        for (int i = 0; i < n; i++)
        {
            terms[i] = fit.Fitted[i] - old[i];
            residual[i] = values[i] + old[i] - terms[i];
            tempered[i] = power * weights[i];
        }
        (double[][] effects, _) = new WeightedWhitening(current.Block, tempered).DrawEffects(residual, sigma2, rng);
        var latent = new RandomEffectsLatent(current.Block, effects);
        if (BassStructureMoves.Follow(state, _mixture.At(working, latent)) is not null) return false;
        latent.AddTo(terms);
        moves.Hand(fit with { State = state.Current, Fitted = terms, Groups = latent.Groups, GroupPrior = latent.Prior });
        return true;
    }
}

/// <summary>
/// A drawn scale's step under a noise (τ's slot): u is a latent, so s_k reads only u_k's prior N(0, σ²s_k²) on its q_k contrasts and its
/// own half-Cauchy. A Metropolis step on t = ln s_k, with the Jacobian eᵗ, accepts with
/// <para>ln A = ln p(s′) − ln p(s) + t − q_k·t − ‖u_k‖²/(2σ²)·(1/s′² − 1/s²) + ln q(−t | x′) − ln q(t | x),</para>
/// <para>t ~ N(0, h(x)²), h = 2.4/√C, C = 2q_k, the information q_k effects carry about ln s, plus the prior's curvature in ln s, at least 1.
/// No refit: the offset Z̃u does not move. Its law at every rung is u's prior's, which tempering leaves alone.</para>
/// </summary>
public sealed class BassLatentGroupScaleBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private const double Step = 2.4;
    private readonly int _group;
    private readonly string _name;
    private readonly BassGroupScalePrior _prior;
    private readonly VarianceMeanMixture _mixture;

    internal BassLatentGroupScaleBlock(BassModel model, int group)
    {
        ArgumentNullException.ThrowIfNull(model);
        RandomEffectsBlock block = model.RandomEffects ?? throw new ArgumentException("The model holds no random effects.", nameof(model));
        if ((uint)group >= (uint)block.Count) throw new ArgumentOutOfRangeException(nameof(group));
        _group = group;
        _prior = block.Specifications[group].ScalePrior ?? throw new ArgumentException("The scale is held.", nameof(group));
        _mixture = model.Observation as VarianceMeanMixture ?? throw new ArgumentException("u is a latent under a noise.", nameof(model));
        IBassRandomEffect spec = block.Specifications[group];
        _name = spec is BassRandomSlope slope ? $"group scale {spec.Grouping} on {slope.Input}" : $"group scale {spec.Grouping}";
    }

    /// <inheritdoc/>
    public string Name => _name;

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        WorkingResponse working = BassStructureMoves.Of(state).DesignAt(state.Current).Working;
        RandomEffectsLatent latent = working.Effects ?? throw new InvalidOperationException("The reading carries no random effects.");
        double sigma2 = working.Latents?.Sigma2 ?? throw new InvalidOperationException("Random effects under a noise keep σ² on the record.");
        int q = latent.Block.Layout(_group).Dimension;
        double s = latent.Block.Scales.Span[_group], h = Scale(s, q), t = h * Normal.SampleStandard(rng), moved = s * Math.Exp(t);
        if (!(moved > 0.0) || !double.IsFinite(moved)) return false;
        double reverse = Scale(moved, q), squares = latent.Squares(_group);
        double logA = _prior.LogDensity(moved) - _prior.LogDensity(s) + t - q * t - squares / (2.0 * sigma2) * (1.0 / (moved * moved) - 1.0 / (s * s))
                      + Math.Log(h) - Math.Log(reverse) - t * t / (2.0 * reverse * reverse) + t * t / (2.0 * h * h);
        if (!(Math.Log(1.0 - rng.NextDouble()) < logA)) return false;
        return BassStructureMoves.Follow(state, _mixture.At(working, latent.At(latent.Block.At(_group, moved)))) is null;
    }

    private double Scale(double s, int q)
    {
        double c = 2.0 * q + _prior.CurvatureInLog(s);
        return Step / Math.Sqrt(double.IsFinite(c) ? Math.Max(c, 1.0) : 1.0);
    }
}

/// <summary>
/// A drawn scale's non-centred step under a noise (E21's slice 6; decision 7, the user's call of 2026-11-13): the ancillary half of Yu and
/// Meng's interweaving, beside the centred step (<see cref="BassLatentGroupScaleBlock"/>), when <see cref="BassSampling.InterweaveGroupScales"/>
/// asks for it. The centred step moves s_k only through u_k, and where the data pin u_k weakly that is a funnel. This step holds
/// ũ_k = u_k/(σs_k), so t = ln(s′/s) scales u_k by eᵗ and moves the offset Z̃u, and accepts with
/// <para>ln A = ℓ_b(s′) − ℓ_b(s) + ln p(s′) − ln p(s) + t + ln q(−t | x′) − ln q(t | x).</para>
/// <para>ℓ_b is the rung's tempered target at the reading, with a integrated, as τ's step reads it before the coefficients' slot draws a
/// (the partially collapsed order). ũ's standard normal law is free of s and cancels. Under a noise the target holds no scale prior, which
/// a whitening alone adds (<see cref="BassTarget"/>), so the step adds ln p(s′) − ln p(s) and the Jacobian. The proposal is t ~ N(0, h(x)²)
/// with h = 2.4/√C, where C = b·Σᵢ wᵢ(z̃ᵢu_{k,g(i)})²/σ², the information the data carry about ln s with ũ held, plus the prior's
/// curvature in ln s, at least 1. A scalar response's; several outputs' U is E34's.</para>
/// </summary>
public sealed class BassGroupScaleInterweavingBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private const double Step = 2.4;
    private readonly int _group;
    private readonly string _name;
    private readonly BassGroupScalePrior _prior;
    private readonly VarianceMeanMixture _mixture;

    internal BassGroupScaleInterweavingBlock(BassModel model, int group)
    {
        ArgumentNullException.ThrowIfNull(model);
        RandomEffectsBlock block = model.RandomEffects ?? throw new ArgumentException("The model holds no random effects.", nameof(model));
        if ((uint)group >= (uint)block.Count) throw new ArgumentOutOfRangeException(nameof(group));
        if (model.Data.Response is MultiOutputResponse) throw new ArgumentException("Several outputs' U has steps of its own (E34).", nameof(model));
        _group = group;
        _prior = block.Specifications[group].ScalePrior ?? throw new ArgumentException("The scale is held.", nameof(group));
        _mixture = model.Observation as VarianceMeanMixture ?? throw new ArgumentException("u is a latent under a noise.", nameof(model));
        IBassRandomEffect spec = block.Specifications[group];
        _name = (spec is BassRandomSlope slope ? $"group scale {spec.Grouping} on {slope.Input}" : $"group scale {spec.Grouping}") + " interwoven";
    }

    /// <inheritdoc/>
    public string Name => _name;

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        BassStructureMoves moves = BassStructureMoves.Of(state);
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        RandomEffectsLatent latent = working.Effects ?? throw new InvalidOperationException("The reading carries no random effects.");
        double sigma2 = working.Latents?.Sigma2 ?? throw new InvalidOperationException("Random effects under a noise keep σ² on the record.");
        double power = BassCoefficientsBlock.Power(state);
        double s = latent.Block.Scales.Span[_group], h = Scale(latent, working.Weights, sigma2, power);
        double t = h * Normal.SampleStandard(rng), moved = s * Math.Exp(t);
        if (!(moved > 0.0) || !double.IsFinite(moved)) return false;
        var effects = (double[][])latent.Effects.Clone();
        double ratio = moved / s;
        effects[_group] = [.. latent.Effects[_group].Select(u => u * ratio)];
        var scaled = new RandomEffectsLatent(latent.Block.At(_group, moved), effects);
        double reverse = Scale(scaled, working.Weights, sigma2, power);
        BassState candidate = moves.ProposeFollow(state.Current, _mixture.At(working, scaled));
        var target = state.Target as BassTarget ?? throw new InvalidOperationException("The chain's target is BASS's.");
        if (target.Refusal(candidate.Fit) is not null) return false;
        ITargetEvaluation evaluation = state.Target.Evaluate(candidate);
        double logA = state.Temperature.LogDensity(evaluation) - state.CurrentLogDensity
                      + _prior.LogDensity(moved) - _prior.LogDensity(s) + t
                      + Math.Log(h) - Math.Log(reverse) - t * t / (2.0 * reverse * reverse) + t * t / (2.0 * h * h);
        if (!(Math.Log(1.0 - rng.NextDouble()) < logA)) return false;
        state.SetState(candidate, evaluation);
        return true;
    }

    // h = 2.4/√C: C = b·Σᵢ wᵢ(z̃ᵢu_{k,g(i)})²/σ² plus the prior's curvature in ln s at the latent's s, at least 1.
    private double Scale(RandomEffectsLatent latent, ReadOnlySpan<double> weights, double sigma2, double power)
    {
        RandomEffectGrouping layout = latent.Block.Layout(_group);
        double[] u = latent.Effects[_group];
        double information = 0.0;
        for (int i = 0; i < weights.Length; i++)
        {
            double shift = layout.Multiplier(i) * u[layout.GroupOf(i)];
            information += weights[i] * shift * shift;
        }
        double c = power * information / sigma2 + _prior.CurvatureInLog(latent.Block.Scales.Span[_group]);
        return Step / Math.Sqrt(double.IsFinite(c) ? Math.Max(c, 1.0) : 1.0);
    }
}

/// <summary>
/// The structure step at b = 1 under a noise with random effects, when <see cref="BassSampling.CollapseEffectsInMoves"/> asks for it
/// (decision 7, an opt-in since 2026-11-10; off by default, form B's move step in its place): the chain follows the u-integrated reading at the
/// current weights (<see cref="VarianceMeanMixture.Integrated"/>), takes one reversible-jump move on that marginal, draws a and then u from
/// their conditionals (<see cref="BassStructureMoves.DrawSweepFit"/>), and follows the reading with u as its offset again. A step on the
/// u-marginal followed by u's draw before any block reads it: the partially collapsed order, so the rung's kernel leaves its target over
/// (structure, latents, u) invariant, and swaps read that shared state. At b ≠ 1 it does nothing; there the moves read the likelihood
/// given u. A plain chain takes it in place of its move step; a replica or a particle after its own moves.
/// </summary>
public sealed class BassCollapsedStructureBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private readonly VarianceMeanMixture _mixture;

    internal BassCollapsedStructureBlock(BassModel model)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (model.RandomEffects is null) throw new ArgumentException("The model holds no random effects.", nameof(model));
        _mixture = model.Observation as VarianceMeanMixture ?? throw new ArgumentException("u is a latent under a noise.", nameof(model));
    }

    /// <inheritdoc/>
    public string Name => "collapsed moves";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        if (BassCoefficientsBlock.Power(state) != 1.0) return false;
        BassStructureMoves moves = BassStructureMoves.Of(state);
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        RandomEffectsLatent current = working.Effects ?? throw new InvalidOperationException("The reading carries no random effects.");
        WorkingResponse integrated = _mixture.Integrated(working);
        if (BassStructureMoves.Follow(state, integrated) is not null)
        {
            // Outside F on the integrated reading: the chain stays on the reading given u.
            if (BassStructureMoves.Follow(state, working) is not null) throw new InvalidOperationException("The chain could not return to its reading.");
            return false;
        }
        BassState before = state.Current;
        BassState after = state.Step();
        BassSweepFit fit = moves.DrawSweepFit(after, rng, state.Target as BassTarget, 1.0);
        BassGroupDraw[] groups = fit.Groups ?? throw new InvalidOperationException("The integrated reading's fit draws u.");
        var latent = new RandomEffectsLatent(current.Block, [.. current.Block.Ordered(groups).Select(g => g.Effects)]);
        if (BassStructureMoves.Follow(state, _mixture.At(integrated, latent)) is not null)
            throw new InvalidOperationException("The structure lies in F on the integrated reading but not on the reading given u.");
        moves.Hand(fit with { State = state.Current });
        return !ReferenceEquals(before, after);
    }
}
