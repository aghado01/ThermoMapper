using System;
using System.Collections.Generic;

namespace Maths.Samplers.FlatHistogram;

/// <summary>
/// A multicanonical run's production draws, per chain in order: each draw's untempered log likelihood ℓ₁ and its
/// functionals, under the frozen weights. What the reweighting reads (<see cref="Reweighting"/>): draw s carries the
/// importance weight <c>e^{b·ℓ₁}/w(ℓ₁)</c> toward the tempered posterior at b, so one run serves every b. A model's own
/// sinks align with these by chain and draw order.
/// </summary>
public sealed class MulticanonicalDraws
{
    private readonly double[][] _levels;
    private readonly double[][] _functionals;   // per chain, [draw · FunctionalDim + f]

    /// <param name="weights">The frozen weights the draws were taken under.</param>
    /// <param name="levels">Per chain, each draw's ℓ₁; every chain holds the same number of draws.</param>
    /// <param name="functionals">Per chain, the draws' functionals, row-major, <paramref name="functionalDim"/> per draw.</param>
    /// <param name="functionalDim">The functionals per draw.</param>
    public MulticanonicalDraws(MulticanonicalWeights weights, IReadOnlyList<double[]> levels, IReadOnlyList<double[]> functionals, int functionalDim)
    {
        ArgumentNullException.ThrowIfNull(weights);
        ArgumentNullException.ThrowIfNull(levels);
        ArgumentNullException.ThrowIfNull(functionals);
        if (levels.Count == 0) throw new ArgumentException("At least one chain is required.", nameof(levels));
        if (functionals.Count != levels.Count) throw new ArgumentException("One functional array per chain is required.", nameof(functionals));
        if (functionalDim < 0) throw new ArgumentOutOfRangeException(nameof(functionalDim));
        int draws = levels[0].Length;
        _levels = new double[levels.Count][];
        _functionals = new double[levels.Count][];
        for (int c = 0; c < levels.Count; c++)
        {
            if (levels[c].Length != draws) throw new ArgumentException("Every chain must hold the same number of draws.", nameof(levels));
            if (functionals[c].Length != draws * functionalDim)
                throw new ArgumentException($"Chain {c} holds {functionals[c].Length} functional values, not {draws * functionalDim}.", nameof(functionals));
            _levels[c] = (double[])levels[c].Clone();
            _functionals[c] = (double[])functionals[c].Clone();
        }
        Weights = weights;
        Draws = draws;
        FunctionalDim = functionalDim;
    }

    /// <summary>The frozen weights.</summary>
    public MulticanonicalWeights Weights { get; }

    /// <summary>The chains.</summary>
    public int Chains => _levels.Length;

    /// <summary>Draws per chain.</summary>
    public int Draws { get; }

    /// <summary>Functionals per draw.</summary>
    public int FunctionalDim { get; }

    /// <summary>Chain <paramref name="chain"/>'s ℓ₁, draw by draw.</summary>
    public ReadOnlySpan<double> Levels(int chain) => _levels[chain];

    /// <summary>Draw <paramref name="draw"/> of chain <paramref name="chain"/>'s functionals.</summary>
    public ReadOnlySpan<double> Functionals(int chain, int draw)
    {
        if ((uint)draw >= (uint)Draws) throw new ArgumentOutOfRangeException(nameof(draw));
        return new ReadOnlySpan<double>(_functionals[chain], draw * FunctionalDim, FunctionalDim);
    }

    /// <summary>Functional <paramref name="functional"/> of every draw of chain <paramref name="chain"/>, in order.</summary>
    public double[] Functional(int chain, int functional)
    {
        if ((uint)functional >= (uint)FunctionalDim) throw new ArgumentOutOfRangeException(nameof(functional));
        var values = new double[Draws];
        double[] row = _functionals[chain];
        for (int s = 0; s < Draws; s++) values[s] = row[s * FunctionalDim + functional];
        return values;
    }
}
