using System;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

// Random intercepts for several outputs under a noise (E34's slice 4, stage 4; decisions 6 and 7): u_g ∈ ℝᴷ ~ N(0, s²Σ) is a latent on every
// rung, read as the offset Z̃u beside vᵢγ (form B, E21's decision 7 as amended). Σ stays integrated in the structure step, so u's prior enters
// its integral as pseudo-rows, IW(ν₀ + Σ_k q_k, Ψ₀ + Σ_k U_kᵀU_k/s_k²) beside γ's (BassTarget.WithPseudoRows): s's step reads the target.

/// <summary>
/// Several outputs' group effects (the coefficients' slot, after Σ and A): U from its conditional given A, Σ, the latents and the scales at
/// the rung's power b, the matrix normal MN(K·Z̃ᵀ(bW)E, K, Σ) in contrasts with K = (S⁻² + Z̃ᵀ(bW)Z̃)⁻¹ and E = Y − vγᵀ − F, through the weighted
/// whitening at bW (<see cref="WeightedWhitening"/>; at b = 0 its prior). Each output's mean is the scalar block's on its residual; the
/// spread is the standardized draw once per output, mixed across the outputs by Σ's factor. The chain follows the offset at the new U, and
/// the sweep's fit is handed on with Z̃U in it, so the blocks after read the residual net of it.
/// </summary>
public sealed class BassOutputGroupEffectsBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private readonly VarianceMeanMixture _mixture;

    internal BassOutputGroupEffectsBlock(BassModel model)
    {
        ArgumentNullException.ThrowIfNull(model);
        _mixture = model.Observation as VarianceMeanMixture ?? throw new ArgumentException("U is a latent under a noise.", nameof(model));
    }

    /// <inheritdoc/>
    public string Name => "group effects";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        BassStructureMoves moves = BassStructureMoves.Of(state);
        BassOutputSweepFit fit = moves.OutputSweepFit is { } drawn && ReferenceEquals(drawn.State, state.Current)
            ? drawn
            : throw new InvalidOperationException("U reads this sweep's Σ and coefficients: the coefficients' block runs first, at the current state.");
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        RandomEffectsLatent current = working.Effects ?? throw new InvalidOperationException("The reading carries no random effects.");
        ObservationLatents latents = working.Latents ?? throw new InvalidOperationException("The working response holds no latents.");
        int k = fit.Sigma.GetLength(0), n = latents.Values.Length;
        double power = BassCoefficientsBlock.Power(state);
        ReadOnlySpan<double> y = working.Response.Values, v = latents.Values, weights = working.Weights;
        // F, the coefficients' part of the fit (the sweep's fit carries the old Z̃U); E = Y − vγᵀ − F.
        double[] terms = (double[])fit.Fitted.Clone();
        var old = new double[n * k];
        current.AddTo(old);
        for (int i = 0; i < terms.Length; i++) terms[i] -= old[i];
        var tempered = new double[n];
        for (int i = 0; i < n; i++) tempered[i] = power * weights[i];
        var whitening = new WeightedWhitening(current.Block, tempered);
        double[]? gamma = latents.OutputSkew;
        var means = new double[k][][];
        var noises = new double[k][][];
        var residual = new double[n];
        for (int l = 0; l < k; l++)
        {
            for (int i = 0; i < n; i++) residual[i] = y[i * k + l] - v[i] * (gamma?[l] ?? 0.0) - terms[i * k + l];
            means[l] = whitening.MeanEffects(residual);
            noises[l] = whitening.DrawEffects(new double[n], 1.0, rng).Effects;
        }
        double[,] root = MatrixTPredictive.Root(fit.Sigma);
        var effects = new double[current.Block.Count][][];
        for (int g = 0; g < effects.Length; g++)
        {
            effects[g] = new double[k][];
            for (int l = 0; l < k; l++)
            {
                int levels = means[l][g].Length;
                var u = new double[levels];
                for (int h = 0; h < levels; h++)
                {
                    double mixed = 0.0;
                    for (int q = 0; q <= l; q++) mixed += root[l, q] * noises[q][g][h];
                    u[h] = means[l][g][h] + mixed;
                }
                effects[g][l] = u;
            }
        }
        RandomEffectsLatent latent = RandomEffectsLatent.OfOutputs(current.Block, effects);
        if (BassStructureMoves.Follow(state, _mixture.At(working, latent)) is not null) return false;
        latent.AddTo(terms);
        moves.HandOutputs(fit with { State = state.Current, Fitted = terms });
        return true;
    }
}

/// <summary>
/// Several outputs' drawn group scale under a noise (τ's slot): U is a latent and Σ is integrated, so s_k reads U_k's prior through the
/// target, whose Σ integral holds U_k's q_k pseudo-rows U_kᵀU_k/s_k² (the matrix t of U given s, <see cref="BassTarget.WithPseudoRows"/>),
/// and its own prior. A Metropolis step on t = ln s_k, with the Jacobian eᵗ, accepts with
/// <para>ln A = Δ ln π_b + ln p(s′) − ln p(s) + t + ln q(−t | x′) − ln q(t | x),</para>
/// <para>t ~ N(0, h(x)²), h = 2.4/√C, C = 2K·q_k, the information K·q_k effects carry about ln s, plus the prior's curvature in ln s, at least
/// 1. The offset Z̃U does not move, nor do the weights; the reading takes the block at s′. At K = 1 with σ² on the record this is E21's
/// <see cref="BassLatentGroupScaleBlock"/>; here σ² is Σ's, integrated.</para>
/// </summary>
public sealed class BassOutputGroupScaleBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private const double Step = 2.4;
    private readonly int _group, _outputs;
    private readonly string _name;
    private readonly BassGroupScalePrior _prior;
    private readonly VarianceMeanMixture _mixture;

    internal BassOutputGroupScaleBlock(BassModel model, int group)
    {
        ArgumentNullException.ThrowIfNull(model);
        RandomEffectsBlock block = model.RandomEffects ?? throw new ArgumentException("The model holds no random effects.", nameof(model));
        if ((uint)group >= (uint)block.Count) throw new ArgumentOutOfRangeException(nameof(group));
        _group = group;
        _prior = block.Specifications[group].ScalePrior ?? throw new ArgumentException("The scale is held.", nameof(group));
        _mixture = model.Observation as VarianceMeanMixture ?? throw new ArgumentException("U is a latent under a noise.", nameof(model));
        _outputs = _mixture.Outputs;
        _name = $"group scale {block.Specifications[group].Grouping}";
    }

    /// <inheritdoc/>
    public string Name => _name;

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        if (Propose(state, Normal.SampleStandard(rng)) is not { } move) return false;
        if (!(Math.Log(1.0 - rng.NextDouble()) < move.LogAcceptance)) return false;
        state.SetState(move.Candidate, move.Evaluation);
        return true;
    }

    /// <summary>A proposed step: the candidate state, its evaluation, and ln A.</summary>
    internal readonly record struct Move(BassState Candidate, ITargetEvaluation Evaluation, double LogAcceptance);

    // The step of z standard deviations of h from the current state, or null when it leaves the scale's support or F.
    internal Move? Propose(ReversibleJumpChain<BassState> state, double z)
    {
        BassStructureMoves moves = BassStructureMoves.Of(state);
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        RandomEffectsLatent latent = working.Effects ?? throw new InvalidOperationException("The reading carries no random effects.");
        int q = latent.Block.Layout(_group).Dimension;
        double s = latent.Block.Scales.Span[_group], h = Scale(s, q), t = h * z, moved = s * Math.Exp(t);
        if (!(moved > 0.0) || !double.IsFinite(moved)) return null;
        double reverse = Scale(moved, q);
        BassState candidate = moves.ProposeFollow(state.Current, _mixture.At(working, latent.At(latent.Block.At(_group, moved))));
        if ((state.Target as BassTarget ?? throw new InvalidOperationException("The chain's target is BASS's.")).Refusal(candidate.Fit) is not null) return null;
        ITargetEvaluation evaluation = state.Target.Evaluate(candidate);
        double target = state.Temperature.LogDensity(evaluation) - state.CurrentLogDensity;
        double logA = target + _prior.LogDensity(moved) - _prior.LogDensity(s) + t
                      + Math.Log(h) - Math.Log(reverse) - t * t / (2.0 * reverse * reverse) + t * t / (2.0 * h * h);
        return new Move(candidate, evaluation, logA);
    }

    private double Scale(double s, int q)
    {
        double c = 2.0 * _outputs * q + _prior.CurvatureInLog(s);
        return Step / Math.Sqrt(double.IsFinite(c) ? Math.Max(c, 1.0) : 1.0);
    }
}
