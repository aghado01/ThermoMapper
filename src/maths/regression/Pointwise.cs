using System;
using System.Collections.Generic;

namespace Maths.Regression;

/// <summary>
/// How observations group for a pointwise density (OQ17). The unit map sends each observation to a unit: what a
/// leave-one-out or add-one-in step removes or adds. The stratum map sends each unit to an exchangeability stratum, which
/// only the consumers read: coverage then holds within each stratum (Mondrian conformal prediction, Vovk; grouped
/// conformal Bayes, Fong &amp; Holmes 2021). The caller chooses both; the engine is built for no one of them. SPC × BARS
/// takes a replicate's observation at a temperature as the unit and the temperature as the stratum, functional data take
/// a curve, and plain regression an observation with one stratum.
/// </summary>
public sealed record UnitPartition
{
    /// <param name="unitOf">Each observation's unit; units are numbered 0..U−1 and each has at least one observation.</param>
    /// <param name="stratumOf">Each unit's stratum, numbered from 0; null puts every unit in one stratum.</param>
    public UnitPartition(int[] unitOf, int[]? stratumOf = null)
    {
        ArgumentNullException.ThrowIfNull(unitOf);
        if (unitOf.Length == 0) throw new ArgumentException("A partition needs at least one observation.", nameof(unitOf));
        int units = 0;
        foreach (int u in unitOf)
        {
            if (u < 0) throw new ArgumentOutOfRangeException(nameof(unitOf), "Unit numbers start at 0.");
            units = Math.Max(units, u + 1);
        }
        var members = new List<int>[units];
        for (int k = 0; k < units; k++) members[k] = new List<int>();
        for (int i = 0; i < unitOf.Length; i++) members[unitOf[i]].Add(i);
        for (int k = 0; k < units; k++)
            if (members[k].Count == 0) throw new ArgumentException($"Unit {k} has no observation; number the units 0..U−1.", nameof(unitOf));
        if (stratumOf is not null)
        {
            if (stratumOf.Length != units) throw new ArgumentException("One stratum per unit is required.", nameof(stratumOf));
            foreach (int s in stratumOf)
                if (s < 0) throw new ArgumentOutOfRangeException(nameof(stratumOf), "Stratum numbers start at 0.");
        }
        UnitOf = (int[])unitOf.Clone();
        StratumOf = stratumOf is null ? new int[units] : (int[])stratumOf.Clone();
        Units = units;
        var arrays = new int[units][];
        for (int k = 0; k < units; k++) arrays[k] = members[k].ToArray();
        Members = arrays;
    }

    /// <summary>Each observation's unit.</summary>
    public int[] UnitOf { get; }

    /// <summary>Each unit's stratum.</summary>
    public int[] StratumOf { get; }

    /// <summary>The number of units U.</summary>
    public int Units { get; }

    /// <summary>The number of observations.</summary>
    public int Observations => UnitOf.Length;

    /// <summary>Each unit's observations, ascending.</summary>
    public IReadOnlyList<int[]> Members { get; }

    /// <summary>Each observation its own unit, all in one stratum: the default.</summary>
    public static UnitPartition EachObservation(int observations)
    {
        if (observations < 1) throw new ArgumentOutOfRangeException(nameof(observations));
        var unitOf = new int[observations];
        for (int i = 0; i < observations; i++) unitOf[i] = i;
        return new UnitPartition(unitOf);
    }

    /// <summary>Consecutive blocks of <paramref name="size"/> observations as units, all in one stratum: whole curves when
    /// replicates are stored one after another.</summary>
    public static UnitPartition ConsecutiveBlocks(int blocks, int size)
    {
        if (blocks < 1) throw new ArgumentOutOfRangeException(nameof(blocks));
        if (size < 1) throw new ArgumentOutOfRangeException(nameof(size));
        var unitOf = new int[blocks * size];
        for (int i = 0; i < unitOf.Length; i++) unitOf[i] = i / size;
        return new UnitPartition(unitOf);
    }
}

/// <summary>Where a latent variable lives: in one observation's term, shared by a group of observations, or by all.</summary>
public enum LatentScope
{
    /// <summary>One observation depends on it, as on a robust weight wᵢ.</summary>
    Observation,

    /// <summary>A group of observations depends on it, as on a replicate's curve β_r or a segment's level.</summary>
    Group,

    /// <summary>Every observation depends on it, as on the coefficients of a single curve or σ².</summary>
    Global,
}

/// <summary>How a pointwise density treated a latent.</summary>
public enum LatentTreatment
{
    /// <summary>Integrated out: everything that depends on it lay inside the unit.</summary>
    Integrated,

    /// <summary>Conditioned on: something that depends on it lay outside the unit.</summary>
    Conditioned,

    /// <summary>Integrated in some units or draws and conditioned in others.</summary>
    Mixed,
}

/// <summary>A latent variable's name, its scope, and how a pointwise density treated it.</summary>
public sealed record LatentReport(string Name, LatentScope Scope, LatentTreatment Treatment);

/// <summary>
/// log p(y_unit | draw) for each retained draw and each unit (OQ17): rows are draws, pooled over chains in chain order, and
/// columns units. A latent is integrated out of a unit's term when everything that depends on it lies inside the unit,
/// and conditioned on otherwise, so the choice of unit is the choice of predictive question (a new curve, or a new point
/// on a seen curve); <see cref="Latents"/> says which latents were treated which way, so a consumer can check the question
/// it answers.
/// </summary>
public sealed record PointwiseLogDensity(double[,] Values, UnitPartition Units, IReadOnlyList<LatentReport> Latents)
{
    /// <summary>The number of draws S.</summary>
    public int Draws => Values.GetLength(0);
}

/// <summary>Tracks whether a latent was integrated, conditioned on, or each, across units and draws.</summary>
internal sealed class TreatmentTally
{
    private bool _integrated, _conditioned;

    public void Integrated() => _integrated = true;

    public void Conditioned() => _conditioned = true;

    public LatentTreatment Result => _integrated && _conditioned ? LatentTreatment.Mixed
                                   : _integrated ? LatentTreatment.Integrated : LatentTreatment.Conditioned;
}
