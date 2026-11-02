using System;
using System.Collections.Generic;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>What a weight states about its observation.</summary>
public enum WeightMeaning
{
    /// <summary>yᵢ ~ N(fᵢ, σ²/wᵢ): one observation at a known relative precision, what a latent scale or Pólya–Gamma hands the
    /// core. σ²'s shape counts the rows, and the likelihood carries ½Σ ln wᵢ.</summary>
    Precision,

    /// <summary>yᵢ counted wᵢ times: the likelihood is Πᵢ N(yᵢ | fᵢ, σ²)^wᵢ, so σ²'s shape counts Σwᵢ and integer weights are
    /// duplicated rows.</summary>
    Replication,
}

/// <summary>Which columns the coefficient prior reads under weights (D03's robust row).</summary>
public enum WeightedPriorGram
{
    /// <summary>The columns as they stand, N rows, free of the weights: the g-form's τG_c on the columns centred over N rows, and
    /// the ridge's D from them. The level and s² read the working response unweighted, so a latent scale's conditional stays
    /// in closed form, at a decomposition per structure (<see cref="PriorPencil"/>). D03's default under robust noise.</summary>
    Unweighted,

    /// <summary>The weighted columns: τZ_cᵀWZ_c and its D, with the level and s² weighted. The closed forms of the unweighted core
    /// hold on the weighted statistics, at O(1) per (τ, b), and integer replication weights are duplicated rows. A latent scale's
    /// conditional gains |Z_cᵀWZ_c|^½, since the prior moves with W.</summary>
    Weighted,
}

/// <summary>
/// What an observation model hands the Gaussian core (F-X1, F-X6): the working response y − offset, the observation
/// weights W, the precisions the core's normal equations carry (ZᵀWZ, ZᵀWỹ), and, for a direct functional response, the
/// inner product Q on the output index that a residual correlated along it declares (covariance σ²(I ⊗ R), Q = R⁻¹).
/// Gaussian BASS hands the response itself, unit weights and Q = I, so its marginal and draws are the core's own. A
/// variance-mean mixture (generalized-hyperbolic noise) hands y − β√w·v with weights 1/(w·v) at the current latent state,
/// where a scale mixture hands weights alone.
/// <para>The core fits the working response under its weights (<see cref="BassDesign"/>), read as <see cref="Meaning"/> says,
/// with the coefficient prior on <see cref="PriorGram"/>'s columns. A Q other than the identity is carried by the contract and
/// refused by today's core, which holds the independent Gram: the centring under Q and the marginal's Σᵢ yᵢᵀQyᵢ and log|R|
/// terms arrive with the extension that supplies them.</para>
/// </summary>
public sealed class WorkingResponse
{
    private readonly double[]? _shifted;
    private readonly double[]? _offset, _weights;

    /// <summary>The working response of <paramref name="response"/> less <paramref name="offset"/> (one value per
    /// observation, case by case and location by location; none is zero), with <paramref name="weights"/> (positive, one per
    /// observation; none is unit), read as <paramref name="meaning"/> says with the coefficient prior on
    /// <paramref name="priorGram"/>'s columns, and <paramref name="output"/>, the inner product on the output index (none is the
    /// identity; a direct functional response's only). The arrays are copied.</summary>
    public WorkingResponse(BassResponse response, ReadOnlySpan<double> offset = default, ReadOnlySpan<double> weights = default,
                           IInnerProduct? output = null, WeightMeaning meaning = WeightMeaning.Precision,
                           WeightedPriorGram priorGram = WeightedPriorGram.Unweighted)
    {
        ArgumentNullException.ThrowIfNull(response);
        if (meaning is not (WeightMeaning.Precision or WeightMeaning.Replication)) throw new ArgumentOutOfRangeException(nameof(meaning));
        if (priorGram is not (WeightedPriorGram.Unweighted or WeightedPriorGram.Weighted)) throw new ArgumentOutOfRangeException(nameof(priorGram));
        int n = response.Observations;
        Response = response;
        Meaning = meaning;
        PriorGram = priorGram;
        Weighting = GPriorWeights.Unit(n);
        Output = output ?? IdentityInnerProduct.Instance;
        if (!Output.IsIdentity)
        {
            if (response is not FunctionalResponse) throw new ArgumentException("A scalar response has no output index.", nameof(output));
            if (Output.Dimension != response.Locations)
                throw new ArgumentException($"The inner product acts on {Output.Dimension} output locations, the response has {response.Locations}.", nameof(output));
        }
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
            var logs = new double[n];
            for (int i = 0; i < n; i++)
            {
                double w = _weights[i];
                if (!(w > 0.0 && double.IsFinite(w))) throw new ArgumentException("A weight is positive and finite.", nameof(weights));
                logs[i] = Math.Log(w);
            }
            double mass = PairwiseSum.Of(_weights);
            Weighting = meaning == WeightMeaning.Replication ? new GPriorWeights(mass, mass, 0.0) : new GPriorWeights(n, mass, PairwiseSum.Of(logs));
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

    /// <summary>What the weights state; unit weights state nothing.</summary>
    public WeightMeaning Meaning { get; }

    /// <summary>Which columns the coefficient prior reads under weights; unit weights make them one.</summary>
    public WeightedPriorGram PriorGram { get; }

    /// <summary>The count, the mass and Σ ln wᵢ the marginal reads beyond the weighted statistics.</summary>
    public GPriorWeights Weighting { get; }

    /// <summary>Whether a coefficient prior apart from the likelihood's Gram is needed: weights, with the prior free of them.</summary>
    public bool PriorStandsApart => _weights is not null && PriorGram == WeightedPriorGram.Unweighted;

    /// <summary>ȳ, the working response's mean under W over all N observations: Σwᵢyᵢ/Σwᵢ, the mean at unit weights.</summary>
    public double ResponseMean => _weights is null
        ? PairwiseSum.Of(Values) / Response.Observations
        : PairwiseSum.Dot(_weights, Values) / Weighting.Mass;

    /// <summary>The weights' own array, not to be written; null at unit weights.</summary>
    internal double[]? WeightArray => _weights;

    /// <summary>A column's mean under W, the design's centring, from its case side <paramref name="a"/> and its output side
    /// <paramref name="h"/> for a direct response: āh̄ at unit weights.</summary>
    internal double ColumnMean(double[] a, double[]? h)
    {
        if (_weights is null)
        {
            double mean = PairwiseSum.Of(a) / a.Length;
            return h is null ? mean : mean * (PairwiseSum.Of(h) / h.Length);
        }
        if (h is null) return PairwiseSum.Dot(_weights, a) / Weighting.Mass;
        var byCase = new double[a.Length];
        for (int i = 0; i < a.Length; i++) byCase[i] = PairwiseSum.Dot(_weights.AsSpan(i * h.Length, h.Length), h);
        return PairwiseSum.Dot(a, byCase) / Weighting.Mass;
    }

    /// <summary>Q, the inner product on the output index; the identity for independent noise.</summary>
    public IInnerProduct Output { get; }
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
