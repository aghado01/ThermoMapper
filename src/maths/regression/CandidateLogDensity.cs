using System;
using System.Collections.Generic;

namespace Maths.Regression;

/// <summary>The forms of the conformal layer a candidate seam serves (<see cref="CandidateLogDensity.Forms"/>), as flags.</summary>
[Flags]
public enum ConformalForms
{
    /// <summary>The full step, which ranks a candidate among the training units it is exchangeable with, every one scored in sample.</summary>
    Full = 1,

    /// <summary>The split form, and the inductive streams that score as it does: units ranked among calibration units scored by a fixed fit.</summary>
    Split = 2,

    /// <summary>Both, the default.</summary>
    Both = Full | Split,
}

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

    /// <summary>Each draw's log CDF at each candidate, ln F(y | x_t, θ_s), in <see cref="Values"/>' layout and with the latents treated
    /// as <see cref="Latents"/> reports, both tails kept as <see cref="PointwiseLogDensity.LogCdf"/> keeps them; null where the engine
    /// gave none. The conformal predictive distribution reads it.</summary>
    public double[,]? LogCdf
    {
        get;
        init => field = PointwiseLogDensity.CheckLogCdf(value, Values);
    }

    /// <summary>
    /// The forms of the conformal layer this seam serves, both by default. An engine whose candidates break a form's premise names the
    /// forms that remain, with the reason (<see cref="FormsReason"/>), and the layer refuses the seam in the form left out, quoting it.
    /// </summary>
    public ConformalForms Forms
    {
        get;
        init => field = value is ConformalForms.Full or ConformalForms.Split or ConformalForms.Both
            ? value
            : throw new ArgumentOutOfRangeException(nameof(Forms), "A seam serves the full step, the split form or both.");
    } = ConformalForms.Both;

    /// <summary>Why a form is left out of <see cref="Forms"/>; null where the seam serves both.</summary>
    public string? FormsReason { get; init; }

    /// <summary>Refuses the seam in <paramref name="form"/> when it does not serve it, quoting the engine's reason.</summary>
    internal void CheckServes(ConformalForms form, string parameter)
    {
        if ((Forms & form) != 0) return;
        string name = form == ConformalForms.Full ? "the full step" : "the split form";
        throw new ArgumentException($"This seam does not serve {name}: {FormsReason ?? "its engine gave no reason."}", parameter);
    }

    /// <summary>The number of draws S.</summary>
    public int Draws => Values.GetLength(0);

    /// <summary>The number of test cases.</summary>
    public int Cases => Grids.Count;

    /// <summary>Case <paramref name="t"/>'s first column; its grid's values follow in order.</summary>
    public int Offset(int t) => _offsets[t];
}

/// <summary>
/// ln F(e) under a noise at a draw, each half on its own tail, for every engine's candidate seam: above the median ln(1 − S(e)) from the
/// survival, which a law integrates on the upper side. The mixtures' tables end at the 10⁻¹² tails, and beyond the upper one their
/// <c>LogCdf</c> integrates the lower side, losing 1 − F's digits and reading a few units of 10⁻¹⁵ above 0. The closed forms' own
/// <c>LogCdf</c> takes this path already.
/// </summary>
internal static class NoiseLogCdf
{
    private const double Ln2 = 0.69314718055994530942;

    public static double Of(Maths.Distributions.IContinuousDistribution noise, double e)
    {
        double lower = noise.LogCdf(e);
        return lower > -Ln2 ? GPrior.LogFunctions.Log1mExp(noise.LogSurvival(e)) : lower;
    }

    /// <summary>The same split for a law given by its two tails' logs, ln F(e) and ln S(e), each read where it keeps its digits.</summary>
    public static double FromTails(double logCdf, double logSurvival) => logCdf > -Ln2 ? GPrior.LogFunctions.Log1mExp(logSurvival) : logCdf;
}

/// <summary>What every engine's candidate seam checks alike: one grid per case, and units of one observation where a CDF is read.</summary>
internal static class CandidateSeam
{
    /// <summary>Each case's first column, the grids' candidates counted case by case; one grid per case of <paramref name="cases"/>.</summary>
    public static int[] Offsets(IReadOnlyList<double[]> grids, int cases)
    {
        ArgumentNullException.ThrowIfNull(grids);
        if (grids.Count != cases) throw new ArgumentException($"One grid per case is required: {grids.Count} grids for {cases} cases.", nameof(grids));
        var offsets = new int[cases + 1];
        for (int t = 0; t < cases; t++)
            offsets[t + 1] = checked(offsets[t] + (grids[t] ?? throw new ArgumentException($"Case {t}'s grid is null.", nameof(grids))).Length);
        return offsets;
    }

    /// <summary>Refuses a partition with a unit of more than one observation: a CDF orders one value.</summary>
    public static void OneObservationEach(UnitPartition units, string parameter)
    {
        for (int u = 0; u < units.Units; u++)
            if (units.Members[u].Length != 1)
                throw new ArgumentException($"Unit {u} holds {units.Members[u].Length} observations; a CDF orders one value.", parameter);
    }
}
