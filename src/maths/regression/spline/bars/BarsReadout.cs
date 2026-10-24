using System;
using System.Collections.Generic;
using Maths.Samplers.Ensemble;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// One BARS chain's readout statistics, mergeable across chains (<see cref="IMergeable{TSelf}"/>): the conditional
/// variance field behind <c>FitSd</c>, every evaluated configuration's admission outcome, and every draw's retained
/// state (<see cref="BarsDraw"/>), in draw order. The peak readouts (the argmax summary, its trace and sample, the peak
/// count, λ and π) are read from the retained draws after the run.
/// </summary>
internal sealed class BarsReadout : IMergeable<BarsReadout>
{
    public BarsReadout(int gridSize)
    {
        ConditionalVariance = new WeightedSums(gridSize);
    }

    /// <summary>Σ of each draw's conditional variance of the fit at each grid point (folded only when the model has one).</summary>
    public WeightedSums ConditionalVariance { get; }

    /// <summary>Every evaluated configuration's admission outcome, the chain's proposals and its ladder's included.</summary>
    public Tally<DesignRejection> Evaluations { get; } = new();

    /// <summary>Every draw's retained state, in draw order; a merge appends the other readout's draws. A draw's arrays
    /// are never written after it is retained, so a copy shares them.</summary>
    public DrawStore<BarsDraw> Draws { get; } = new();

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
        Evaluations.MergeFrom(other.Evaluations);
        Draws.AppendFrom(other.Draws);
    }
}
