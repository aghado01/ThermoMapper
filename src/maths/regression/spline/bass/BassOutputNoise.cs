using System;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// Several outputs' coefficients block's draw in a sweep under a noise (E34's decision 6): Σ and the centred coefficients A from their
/// conditional at the rung's power given the structure, τ and the reading (its weights 1/(c·vᵢ) and offset vᵢγ), with the fit f at every
/// case and output, case by case and output by output. Redrawn each sweep before it is read, so not state.
/// </summary>
internal sealed record BassOutputSweepFit(BassState State, double[,] Sigma, double[,] Centred, double[] Fitted)
{
    /// <summary>The structure the draw was made at.</summary>
    public BassTerm[] Structure => State.TermArray;
}

/// <summary>
/// Several outputs' latent scales (E34's decision 6, the latent slot): each case's vᵢ, shared by its K outputs, by its law's step under its
/// conditional given the residual rᵢ = yᵢ − fᵢ of this sweep's Σ and coefficients, ∝ p(v)·v^(−bK/2)·exp(−b·rᵢᵀΣ⁻¹rᵢ/(2cv) − b·γᵀΣ⁻¹γ·v/(2c))
/// (<see cref="IMixingLaw.StepConditionalForm"/>); then the chain follows the reading at the new v. At K = 1 it is one response's latent
/// block, rᵢ²/σ² in place of rᵢᵀΣ⁻¹rᵢ. A draw that puts the structure outside F is refused, as one response's is.
/// </summary>
public sealed class BassOutputLatentScaleBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private readonly VarianceMeanMixture _model;

    /// <summary>The block of <paramref name="model"/>'s latents for several outputs.</summary>
    public BassOutputLatentScaleBlock(VarianceMeanMixture model)
    {
        ArgumentNullException.ThrowIfNull(model);
        _model = model;
    }

    /// <inheritdoc/>
    public string Name => "latent scales";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        BassStructureMoves moves = BassStructureMoves.Of(state);
        BassOutputSweepFit fit = moves.OutputSweepFit is { } drawn && ReferenceEquals(drawn.State, state.Current)
            ? drawn
            : throw new InvalidOperationException("The latent scales read this sweep's Σ and coefficients: the coefficients' block runs first, at the current state.");
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        ObservationLatents latents = working.Latents ?? throw new InvalidOperationException("The working response holds no latents.");
        ReadOnlySpan<double> y = working.Response.Values, current = latents.Values;
        int k = fit.Sigma.GetLength(0), n = current.Length;
        double[,] root = MatrixTPredictive.Root(fit.Sigma);
        double c = _model.VarianceScale, power = BassCoefficientsBlock.Power(state);
        double drift = latents.OutputSkew is { } gamma ? Energy(root, gamma) / c : 0.0;
        IMixingLaw mixing = _model.MixingAt(latents);
        var v = new double[n];
        var r = new double[k];
        for (int i = 0; i < n; i++)
        {
            for (int l = 0; l < k; l++) r[l] = y[i * k + l] - fit.Fitted[i * k + l];
            v[i] = mixing.StepConditionalForm(current[i], Energy(root, r) / c, drift, k, power, rng);
        }
        WorkingResponse moved = _model.At(working, new ObservationLatents(v, 0.0, null, latents.Gamma) { OutputSkew = latents.OutputSkew });
        return BassStructureMoves.Follow(state, moved) is null;
    }

    /// <summary>xᵀΣ⁻¹x = ‖L⁻¹x‖² for Σ = LLᵀ, by forward substitution.</summary>
    internal static double Energy(double[,] root, ReadOnlySpan<double> x)
    {
        int k = x.Length;
        Span<double> u = stackalloc double[k];
        double sum = 0.0;
        for (int a = 0; a < k; a++)
        {
            double e = x[a];
            for (int q = 0; q < a; q++) e -= root[a, q] * u[q];
            u[a] = e / root[a, a];
            sum += u[a] * u[a];
        }
        return sum;
    }
}

/// <summary>
/// Several outputs' noise shape (E34's decision 6, the scale slot): the skew γ ∈ ℝᴷ where it is drawn, and the Normal–Wald's tail where it
/// is, each from its conditional. With rᵢ = yᵢ − fᵢ of this sweep's Σ and coefficients, γ | Σ ~ N(0, s²Σ) a priori and
/// rᵢ = vᵢγ + √(c·vᵢ)·L_Σzᵢ, γ | · ~ N(m, Σ/P) with P = bΣvᵢ/c + 1/s² and m = (b/c)Σrᵢ/P, at the rung's power b; at K = 1 it is one
/// response's β = γ/σ draw at a zero prior mean. The tail reads v alone, as one response's does. The chain follows the reading at the new
/// γ, its weights shared; a draw that puts the structure outside F is refused.
/// </summary>
public sealed class BassOutputShapeBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private readonly VarianceMeanMixture _model;

    /// <summary>The block of <paramref name="model"/>'s skew and tail for several outputs.</summary>
    public BassOutputShapeBlock(VarianceMeanMixture model)
    {
        ArgumentNullException.ThrowIfNull(model);
        _model = model;
    }

    /// <inheritdoc/>
    public string Name => "noise shape";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        BassStructureMoves moves = BassStructureMoves.Of(state);
        BassOutputSweepFit fit = moves.OutputSweepFit is { } drawn && ReferenceEquals(drawn.Structure, state.Current.TermArray)
            ? drawn
            : throw new InvalidOperationException("The noise's shape reads this sweep's Σ and coefficients: the coefficients' block runs first, at the current structure.");
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        ObservationLatents latents = working.Latents ?? throw new InvalidOperationException("The working response holds no latents.");
        ReadOnlySpan<double> y = working.Response.Values, v = latents.Values;
        int k = fit.Sigma.GetLength(0), n = v.Length;
        double sumV = 0.0;
        foreach (double vi in v) sumV += vi;
        double[]? gamma = latents.OutputSkew;
        if (_model.Noise.SkewPrior is { IsHeld: false } skewPrior)
        {
            double c = _model.VarianceScale, power = BassCoefficientsBlock.Power(state);
            double precision = power * sumV / c + 1.0 / (skewPrior.Sd * skewPrior.Sd);
            double[,] root = MatrixTPredictive.Root(fit.Sigma);
            var sumR = new double[k];
            for (int i = 0; i < n; i++)
                for (int l = 0; l < k; l++) sumR[l] += y[i * k + l] - fit.Fitted[i * k + l];
            var z = new double[k];
            for (int l = 0; l < k; l++) z[l] = Normal.SampleStandard(rng);
            gamma = new double[k];
            double spread = 1.0 / Math.Sqrt(precision);
            for (int l = 0; l < k; l++)
            {
                double e = 0.0;
                for (int q = 0; q <= l; q++) e += root[l, q] * z[q];
                gamma[l] = power * sumR[l] / c / precision + spread * e;
            }
        }
        double? tail = latents.Gamma;
        if (_model.Noise.GammaPrior is { IsHeld: false } gammaPrior)
        {
            double inverse = 1.0 / (gammaPrior.Sd * gammaPrior.Sd), precision = sumV + inverse;
            tail = Normal.MeanSd((v.Length + gammaPrior.Mean * inverse) / precision, 1.0 / Math.Sqrt(precision)).SampleOn(rng, 0.0, double.PositiveInfinity);
        }
        WorkingResponse moved = _model.At(working, latents.WithOutputShape(gamma, tail), weightsMoved: false);
        return BassStructureMoves.Follow(state, moved) is null;
    }
}

/// <summary>
/// Several outputs' skew rescaling, the γ·v move (E34's decision 6; one response's β·v move carried to K outputs, the latent slot after v's
/// block, under a drawn skew): (γ, v) → (γ/κ, κ·v) along the ridge vᵢγ holds, the offset standing and the weights 1/(c·vᵢ) moving to W/κ,
/// refitted once. On the sweep's joint with Σ and the coefficients integrated, π ∝ p(structure, τ)·p(γ)·Πp_v(vᵢ)·ℓ_b, the target reading γ's
/// prior (its pseudo-row in Σ's posterior and its multivariate t), the map's Jacobian is κ^(N−K), so with t = ln κ it accepts with
/// <para>ln A = Δ ln π_b + Σᵢ [ln p_v(κ·vᵢ) − ln p_v(vᵢ)] + (N − K)·t + ln q(−t | x′) − ln q(t | x),</para>
/// <para>declining a move that leaves F. It is a step on the coefficients' marginal, so after an accepted move Σ and the coefficients are
/// drawn afresh before a block reads them (van Dyk and Park 2008), and its proposal reads the state alone: t ~ N(0, h(x)²), h = 2.4/√C(x),
/// C the likelihood's b·N·K/2 less the second derivative of v's prior along t, at least 1.</para>
/// </summary>
public sealed class BassOutputSkewRescaleBlock : IBlock<ReversibleJumpChain<BassState>>
{
    private const double Step = 2.4;
    private const double Difference = 1e-2;
    private readonly VarianceMeanMixture _model;

    internal BassOutputSkewRescaleBlock(VarianceMeanMixture model)
    {
        ArgumentNullException.ThrowIfNull(model);
        if (model.Noise.SkewPrior is not { IsHeld: false }) throw new ArgumentException("The skew's rescaling needs a drawn skew.", nameof(model));
        _model = model;
    }

    /// <inheritdoc/>
    public string Name => "skew rescale";

    /// <inheritdoc/>
    public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(state);
        ArgumentNullException.ThrowIfNull(rng);
        BassStructureMoves moves = BassStructureMoves.Of(state);
        WorkingResponse working = moves.DesignAt(state.Current).Working;
        ObservationLatents latents = working.Latents ?? throw new InvalidOperationException("The working response holds no latents.");
        double[] gamma = latents.OutputSkew ?? throw new InvalidOperationException("A drawn skew's latents hold γ.");
        double power = BassCoefficientsBlock.Power(state);
        IMixingLaw mixing = _model.MixingAt(latents);
        double[] v = latents.Values.ToArray();
        int k = gamma.Length;
        double h = Scale(v, k, power, mixing), t = h * Normal.SampleStandard(rng);
        if (Propose(state, working, latents, v, gamma, mixing, power, h, t) is not { } move) return false;
        if (!(Math.Log(1.0 - rng.NextDouble()) < move.LogAcceptance)) return false;
        state.SetState(move.Candidate, move.Evaluation);
        moves.DrawOutputSweepFit(state.Current, rng, state.Target as BassTarget, power);   // a step on the coefficients' marginal
        return true;
    }

    /// <summary>A proposed move: the candidate state, its evaluation, and ln A.</summary>
    internal readonly record struct Move(BassState Candidate, ITargetEvaluation Evaluation, double LogAcceptance);

    // h = Step/√C, C = b·N·K/2 − ∂²/∂t² Σᵢ ln p_v(vᵢ·eᵗ) at t = 0, at least 1.
    internal static double Scale(double[] v, int k, double power, IMixingLaw mixing)
    {
        double up = Math.Exp(Difference), down = Math.Exp(-Difference), centre = 0.0, above = 0.0, below = 0.0;
        foreach (double vi in v)
        {
            centre += mixing.LogPrior(vi);
            above += mixing.LogPrior(vi * up);
            below += mixing.LogPrior(vi * down);
        }
        double curvature = 0.5 * power * v.Length * k - (above - 2.0 * centre + below) / (Difference * Difference);
        return Step / Math.Sqrt(double.IsFinite(curvature) ? Math.Max(curvature, 1.0) : 1.0);
    }

    // The move at t, or null when it leaves F.
    internal Move? Propose(ReversibleJumpChain<BassState> state, WorkingResponse working, ObservationLatents latents, double[] v, double[] gamma,
                           IMixingLaw mixing, double power, double h, double t)
    {
        double kappa = Math.Exp(t), inverse = Math.Exp(-t);
        int n = v.Length, k = gamma.Length;
        var moved = new double[n];
        double latentRatio = 0.0;
        for (int i = 0; i < n; i++)
        {
            moved[i] = kappa * v[i];
            latentRatio += mixing.LogPrior(moved[i]) - mixing.LogPrior(v[i]);
        }
        var skew = new double[k];
        for (int l = 0; l < k; l++) skew[l] = gamma[l] * inverse;
        double reverse = Scale(moved, k, power, mixing);
        double proposalRatio = Math.Log(h) - Math.Log(reverse) - t * t / (2.0 * reverse * reverse) + t * t / (2.0 * h * h);
        WorkingResponse next = _model.At(working, new ObservationLatents(moved, 0.0, null, latents.Gamma) { OutputSkew = skew });
        BassState candidate = BassStructureMoves.Of(state).ProposeFollow(state.Current, next);
        if ((state.Target as BassTarget ?? throw new InvalidOperationException("The chain's target is BASS's.")).Refusal(candidate.Fit) is not null) return null;
        ITargetEvaluation evaluation = state.Target.Evaluate(candidate);
        double target = state.Temperature.LogDensity(evaluation) - state.CurrentLogDensity;
        return new Move(candidate, evaluation, target + latentRatio + (n - k) * t + proposalRatio);
    }
}
