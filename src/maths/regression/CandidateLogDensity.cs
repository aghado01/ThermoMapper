using System;
using System.Collections.Generic;

namespace Maths.Regression;

/// <summary>
/// The candidate density seam (E16): log p(y | x_t, draw) for each retained draw and each value y on a test case t's grid,
/// beside a training <see cref="PointwiseLogDensity"/> over the same draws in the same rows. An engine joins the conformal
/// layer (<see cref="Conformal.ConformalBayes"/>) by producing it at new cases; the layer knows no engine. Each candidate is one
/// more unit of the training partition's kind, a value of the response at case t, and its column is read as that unit's term
/// would be: <see cref="Latents"/> reports how each latent was treated, which the layer requires to match the training
/// density's (CB-4). The columns run case by case, each case's grid in ascending order (<see cref="Offset"/>); each case
/// sits in a stratum of the training partition (<see cref="UnitPartition.StratumOf"/>), or in one no training unit holds.
/// </summary>
public sealed record CandidateLogDensity
{
    private readonly IReadOnlyList<double>? _weights;

    /// <param name="values">The draws × candidates matrix of log terms, case by case; finite or −∞ (a value outside the
    /// draw's support).</param>
    /// <param name="grids">Each case's candidate values, finite and strictly ascending, at least one.</param>
    /// <param name="latents">How each latent was treated in the candidates' terms.</param>
    /// <param name="stratumOf">Each case's stratum, numbered from 0; null puts every case in stratum 0.</param>
    /// <param name="chainLengths">The draws per chain, in row order: the training density's.</param>
    public CandidateLogDensity(double[,] values, IReadOnlyList<double[]> grids, IReadOnlyList<LatentReport> latents,
                               int[]? stratumOf = null, IReadOnlyList<int>? chainLengths = null)
    {
        ArgumentNullException.ThrowIfNull(values);
        ArgumentNullException.ThrowIfNull(grids);
        ArgumentNullException.ThrowIfNull(latents);
        if (grids.Count == 0) throw new ArgumentException("At least one test case is required.", nameof(grids));
        var offsets = new int[grids.Count + 1];
        var copies = new double[grids.Count][];
        for (int t = 0; t < grids.Count; t++)
        {
            double[] grid = grids[t] ?? throw new ArgumentException($"Case {t}'s grid is null.", nameof(grids));
            if (grid.Length == 0) throw new ArgumentException($"Case {t}'s grid holds no value.", nameof(grids));
            for (int g = 0; g < grid.Length; g++)
            {
                if (!double.IsFinite(grid[g])) throw new ArgumentException($"Case {t}'s grid holds a value that is not finite.", nameof(grids));
                if (g > 0 && !(grid[g] > grid[g - 1])) throw new ArgumentException($"Case {t}'s grid must be strictly ascending.", nameof(grids));
            }
            copies[t] = (double[])grid.Clone();
            offsets[t + 1] = checked(offsets[t] + grid.Length);
        }
        if (values.GetLength(1) != offsets[^1])
            throw new ArgumentException($"The matrix has {values.GetLength(1)} columns; the grids hold {offsets[^1]} candidates.", nameof(values));
        if (values.GetLength(0) < 1) throw new ArgumentException("At least one draw is required.", nameof(values));
        foreach (double v in values)
            if (double.IsNaN(v) || double.IsPositiveInfinity(v)) throw new ArgumentException("A candidate's log term is NaN or +∞.", nameof(values));
        if (stratumOf is not null)
        {
            if (stratumOf.Length != grids.Count) throw new ArgumentException("One stratum per case is required.", nameof(stratumOf));
            foreach (int s in stratumOf)
                if (s < 0) throw new ArgumentOutOfRangeException(nameof(stratumOf), "Stratum numbers start at 0.");
        }
        Values = values;
        Grids = copies;
        Latents = latents;
        StratumOf = stratumOf is null ? new int[grids.Count] : (int[])stratumOf.Clone();
        ChainLengths = chainLengths;
        _offsets = offsets;
    }

    private readonly int[] _offsets;

    /// <summary>The draws × candidates matrix of log terms.</summary>
    public double[,] Values { get; }

    /// <summary>Each case's candidate values, ascending.</summary>
    public IReadOnlyList<double[]> Grids { get; }

    /// <summary>How each latent was treated in the candidates' terms.</summary>
    public IReadOnlyList<LatentReport> Latents { get; }

    /// <summary>Each case's stratum.</summary>
    public int[] StratumOf { get; }

    /// <summary>The draws per chain, in row order; null for independent draws. They must sum to <see cref="Draws"/>.</summary>
    public IReadOnlyList<int>? ChainLengths
    {
        get;
        init => field = PointwiseLogDensity.CheckChains(value, Values);
    }

    /// <summary>The draws' weights in row order, as <see cref="PointwiseLogDensity.Weights"/>: the training density's.</summary>
    public IReadOnlyList<double>? Weights
    {
        get => _weights;
        init => _weights = PointwiseLogDensity.CheckWeights(value, Values, ChainLengths);
    }

    /// <summary>The number of draws S.</summary>
    public int Draws => Values.GetLength(0);

    /// <summary>The number of test cases.</summary>
    public int Cases => Grids.Count;

    /// <summary>Case <paramref name="t"/>'s first column; its grid's values follow in order.</summary>
    public int Offset(int t) => _offsets[t];
}
