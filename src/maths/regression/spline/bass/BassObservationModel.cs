using System;

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
/// The seam between the noise and the Gaussian core (F-X1): the observation model supplies the core its
/// <see cref="WorkingResponse"/> at the current latent state. Gaussian BASS is the one member
/// (<see cref="GaussianObservation"/>); generalized-hyperbolic noise (GBMARS) joins by supplying an offset and weights,
/// and its latent scales through the sweep's slots.
/// </summary>
public interface IBassObservationModel
{
    /// <summary>The working response and weights the core fits for <paramref name="response"/>.</summary>
    WorkingResponse Working(BassResponse response);
}

/// <summary>Independent Gaussian noise: the working response is the response, with zero offset and unit weights, so the
/// core's marginal and draws are BASS's own.</summary>
public sealed class GaussianObservation : IBassObservationModel
{
    /// <summary>The one Gaussian observation model.</summary>
    public static GaussianObservation Instance { get; } = new();

    private GaussianObservation() { }

    /// <inheritdoc/>
    public WorkingResponse Working(BassResponse response) => new(response);
}
