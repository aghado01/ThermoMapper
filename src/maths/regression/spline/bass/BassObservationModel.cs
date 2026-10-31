using System;
using System.Collections.Generic;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// What an observation model hands the Gaussian core (F-X1): the working response y − offset and the observation weights
/// W, the precisions the core's normal equations carry (ZᵀWZ, ZᵀWỹ). Gaussian BASS hands the response itself and unit
/// weights, so its marginal and draws are the core's own. A variance-mean mixture (generalized-hyperbolic noise) hands
/// y − β√w·v with weights 1/(w·v) at the current latent state, where a scale mixture hands weights alone.
/// <para>The core fits the working response. Weights are carried by the contract and refused by today's core
/// (<see cref="BassDesign"/>), which holds the unweighted Gram: the weighted Gram, centring and marginal arrive with the
/// extension that supplies weights.</para>
/// </summary>
public sealed class WorkingResponse
{
    private readonly double[]? _shifted;
    private readonly double[]? _offset, _weights;

    /// <summary>The working response of <paramref name="response"/> less <paramref name="offset"/> (one value per
    /// observation, case by case and location by location; none is zero), with <paramref name="weights"/> (positive, one per
    /// observation; none is unit). Both are copied.</summary>
    public WorkingResponse(BassResponse response, ReadOnlySpan<double> offset = default, ReadOnlySpan<double> weights = default)
    {
        ArgumentNullException.ThrowIfNull(response);
        int n = response.Observations;
        Response = response;
        if (!offset.IsEmpty)
        {
            if (offset.Length != n) throw new ArgumentException($"The offset has {offset.Length} values for {n} observations.", nameof(offset));
            _offset = offset.ToArray();
            _shifted = new double[n];
            ReadOnlySpan<double> y = response.Values;
            for (int i = 0; i < n; i++)
            {
                if (!double.IsFinite(_offset[i])) throw new ArgumentException("The offset has a non-finite value.", nameof(offset));
                _shifted[i] = y[i] - _offset[i];
            }
        }
        if (!weights.IsEmpty)
        {
            if (weights.Length != n) throw new ArgumentException($"The weights hold {weights.Length} values for {n} observations.", nameof(weights));
            _weights = weights.ToArray();
            foreach (double w in _weights)
                if (!(w > 0.0 && double.IsFinite(w))) throw new ArgumentException("A weight is positive and finite.", nameof(weights));
        }
    }

    /// <summary>The response the working response derives from.</summary>
    public BassResponse Response { get; }

    /// <summary>y − offset, case by case and location by location; the response itself when there is no offset.</summary>
    public ReadOnlySpan<double> Values => _shifted ?? Response.Values;

    /// <summary>The offset; empty when there is none, which is the zero offset.</summary>
    public ReadOnlySpan<double> Offset => _offset;

    /// <summary>The weights W; empty when there are none, which is unit weights.</summary>
    public ReadOnlySpan<double> Weights => _weights;

    /// <summary>Whether the weights are unit.</summary>
    public bool IsUnweighted => _weights is null;
}

/// <summary>
/// The seam between the noise and the Gaussian core (F-X1, F-X2): the observation model supplies the core its
/// <see cref="WorkingResponse"/> and, for a noise with latent variables, the blocks that update them in the sweep's slots
/// (<see cref="BassSweepSlot"/>). Gaussian BASS is the one member (<see cref="GaussianObservation"/>); generalized-hyperbolic
/// noise (GBMARS) joins by supplying an offset and weights at its latent state, and the blocks of its latent scales.
/// </summary>
public interface IBassObservationModel
{
    /// <summary>The working response and weights the core fits for <paramref name="response"/>.</summary>
    WorkingResponse Working(BassResponse response);

    /// <summary>The blocks that fill <paramref name="slot"/> of the sweep, in the order they run; none leaves the slot empty,
    /// which draws nothing and consumes no random number. Asked of <see cref="BassSweepSlot.Coefficients"/>,
    /// <see cref="BassSweepSlot.Latent"/> and <see cref="BassSweepSlot.Scale"/> only: the structure's slot and τ's are BASS's.</summary>
    IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot);
}

/// <summary>Independent Gaussian noise: the working response is the response, with zero offset and unit weights, so the
/// core's marginal and draws are BASS's own, and every latent slot is empty.</summary>
public sealed class GaussianObservation : IBassObservationModel
{
    /// <summary>The one Gaussian observation model.</summary>
    public static GaussianObservation Instance { get; } = new();

    private GaussianObservation() { }

    /// <inheritdoc/>
    public WorkingResponse Working(BassResponse response) => new(response);

    /// <inheritdoc/>
    public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot) => [];
}

/// <summary>
/// The slots of the sweep (F-X2), in D04's general order, the robust and generalized-hyperbolic one: the structure given the
/// latents, with the coefficients integrated out; τ; the coefficients a, drawn every sweep just before the latents that read
/// the residuals; the latent scales; then σ² when it is not integrated, and the skew. Gaussian BASS fills the first two and
/// draws a and σ² at retention by composition (D04), so its other slots are empty.
/// </summary>
public enum BassSweepSlot
{
    /// <summary>The reversible-jump structure moves.</summary>
    Structure,

    /// <summary>τ's slice step, when τ is sampled.</summary>
    Tau,

    /// <summary>The coefficients a, drawn each sweep when a latent block reads their residuals.</summary>
    Coefficients,

    /// <summary>The per-observation latent scales.</summary>
    Latent,

    /// <summary>σ² when it is not integrated out, and the skew.</summary>
    Scale,
}

/// <summary>The sweep assembled from its slots, in <see cref="BassSweepSlot"/> order: the one place the order is written, for
/// the chain, the tempering replicas and the SMC kernel.</summary>
public static class BassSweepSlots
{
    /// <summary>The slots in the order the sweep runs them.</summary>
    public static IReadOnlyList<BassSweepSlot> Order { get; } =
        [BassSweepSlot.Structure, BassSweepSlot.Tau, BassSweepSlot.Coefficients, BassSweepSlot.Latent, BassSweepSlot.Scale];

    /// <summary>
    /// The blocks of a sweep of <paramref name="model"/> from slot <paramref name="first"/> on, each block's target
    /// <paramref name="target"/>: the structure's move step, τ's slice step when τ is sampled, then the observation model's
    /// blocks for the slots after τ. A chain starts at the structure slot; a tempering replica or an SMC particle, whose moves
    /// the engine runs itself, at τ's.
    /// </summary>
    public static List<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassModel model, ISamplingTarget<BassState> target,
                                                                      BassSweepSlot first = BassSweepSlot.Structure)
    {
        ArgumentNullException.ThrowIfNull(model);
        ArgumentNullException.ThrowIfNull(target);
        var blocks = new List<IBlock<ReversibleJumpChain<BassState>>>();
        foreach (BassSweepSlot slot in Order)
        {
            if (slot < first) continue;
            switch (slot)
            {
                case BassSweepSlot.Structure:
                    blocks.Add(new MoveStep<BassState>());
                    break;
                case BassSweepSlot.Tau:
                    if (model.Target.SamplesTau) blocks.Add(new BassTauBlock(target));
                    break;
                default:
                    blocks.AddRange(model.Observation.Blocks(slot));
                    break;
            }
        }
        return blocks;
    }
}
