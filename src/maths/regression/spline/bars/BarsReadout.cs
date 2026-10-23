using System;
using System.Collections.Generic;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// One BARS chain's readout statistics, mergeable across chains (<see cref="IMergeable{TSelf}"/>): the conditional
/// variance field behind <c>FitSd</c>, the per-draw peak location, height and significant-peak count, the
/// peak-location sample behind its credible interval, every evaluated configuration's admission outcome, and every
/// draw's exact peak readout, from which λ and π are read. Each draw folds once into every field, so their total
/// weights agree; nonlinear reductions (standard deviations, quantiles, modes) read the merged readout once. The peak
/// location's ordered draws, which its diagnostics read, are a <see cref="ChainTrace"/> kept beside this.
/// </summary>
internal sealed class BarsReadout : IMergeable<BarsReadout>
{
    public BarsReadout(int gridSize)
    {
        ConditionalVariance = new WeightedSums(gridSize);
    }

    /// <summary>Σ of each draw's conditional variance of the fit at each grid point (folded only when the model has one).</summary>
    public WeightedSums ConditionalVariance { get; }

    /// <summary>Per draw: the argmax location, the height there, and the significant-peak count.</summary>
    public WeightedSums Peak { get; } = new(3);

    /// <summary>The argmax locations, for quantiles.</summary>
    public DrawSample PeakSample { get; } = new();

    /// <summary>Every evaluated configuration's admission outcome, the chain's proposals and its ladder's included.</summary>
    public Tally<DesignRejection> Evaluations { get; } = new();

    /// <summary>Every draw's exact peak readout, in draw order; a merge appends the other readout's draws.</summary>
    public List<PeakDraw> Draws { get; } = new();

    /// <summary>An independent copy (a merge into an empty readout, which adds nothing to any value).</summary>
    public BarsReadout Copy()
    {
        var copy = new BarsReadout(ConditionalVariance.Length);
        copy.MergeFrom(this);
        return copy;
    }

    public void MergeFrom(BarsReadout other)
    {
        ArgumentNullException.ThrowIfNull(other);
        ConditionalVariance.MergeFrom(other.ConditionalVariance);
        Peak.MergeFrom(other.Peak);
        PeakSample.MergeFrom(other.PeakSample);
        Evaluations.MergeFrom(other.Evaluations);
        Draws.AddRange(other.Draws);
    }
}
