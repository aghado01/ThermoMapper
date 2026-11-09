using System;
using System.Collections.Generic;
using Maths.Rng;

namespace Maths.Regression.Conformal;

/// <summary>One arrival's smoothed p-value (E35, M-1): its conformity <see cref="Score"/>, the scores so far strictly below it
/// (<see cref="Below"/>) and equal to it with itself (<see cref="Tied"/>), the tie-breaking draw θ, and p = (Below + θ·Tied)/n with n the
/// arrivals so far, itself included (<see cref="Count"/>).</summary>
public readonly record struct SmoothedPValue(int Count, double Score, int Below, int Tied, double Theta, double P);

/// <summary>
/// The smoothed conformal p-values of a stream (E35, M-1; Vovk, Nouretdinov and Gammerman 2003, Theorem 1; Vovk et al. 2021,
/// Proposition 1). Each new conformity score, higher conforming better, is ranked among every score so far, itself included:
/// pₙ = (Lₙ + θₙTₙ)/n, with Lₙ the scores strictly below it, Tₙ those equal to it with itself, and θₙ uniform on (0, 1], so p never
/// reads 0. A small p is a nonconforming newest unit.
/// <para><b>The guarantee</b> (<c>seqRank_uniform</c>, <c>seqRank_law</c>): if the stream's units are exchangeable and each score is one
/// fixed function of its unit, p₁, p₂, … are i.i.d. uniform. That is the inductive form, a fit made on a training set proper that
/// is not in the stream (CB-2's conditioning; <see cref="StreamScores"/>). The exchangeability is the caller's premise; the rest the
/// type enforces:</para>
/// <list type="bullet">
/// <item>every score is ranked against the whole stream so far, never a fixed calibration set, whose p-values for several units are
/// each uniform but dependent;</item>
/// <item>ties are broken by a fresh θₙ from the stream's own seeded generator, independent of the scores, never by CB-1's rule, which
/// counts ties for the unit and is not uniform;</item>
/// <item>each score is taken once and never revised, so no score kept from an earlier fit mixes with later ones.</item>
/// </list>
/// Each arrival costs O(log n) (<see cref="OrderStatisticTree{T}"/>). Not thread-safe: one per stream.
/// </summary>
public sealed class SmoothedPValues
{
    private readonly OrderStatisticTree<double> _scores = new();
    private readonly Xoshiro256PlusPlus _theta;

    /// <param name="seed">The tie-breaking draws' seed; the draws are independent of the scores by construction.</param>
    public SmoothedPValues(int seed) => _theta = new Xoshiro256PlusPlus(seed);

    /// <summary>The arrivals so far.</summary>
    public int Count => _scores.Count;

    /// <summary>Ranks <paramref name="score"/> among every score so far, itself included, and gives its smoothed p-value.</summary>
    /// <param name="score">The newest unit's conformity score: finite or ±∞, never NaN.</param>
    public SmoothedPValue Add(double score)
    {
        if (double.IsNaN(score)) throw new ArgumentException("A conformity score is NaN; a NaN does not order.", nameof(score));
        (int below, int tied) = _scores.Add(score);
        tied++;
        int n = _scores.Count;
        double theta = 1.0 - _theta.NextDouble();
        return new SmoothedPValue(n, score, below, tied, theta, (below + theta * tied) / n);
    }

    /// <summary>Each of <paramref name="scores"/> in turn (<see cref="Add"/>).</summary>
    public SmoothedPValue[] AddRange(IEnumerable<double> scores)
    {
        ArgumentNullException.ThrowIfNull(scores);
        var steps = new List<SmoothedPValue>();
        foreach (double score in scores) steps.Add(Add(score));
        return [.. steps];
    }
}

/// <summary>
/// A fixed fit's conformity scores at a stream's units, each higher conforming better (E35's inductive form). The fit is made on a
/// training set proper outside the stream, so every unit is scored by one fixed function (M-1's premise), as the split form scores its
/// calibration units (<see cref="SplitConformal"/>).
/// </summary>
public static class StreamScores
{
    /// <summary>The fit's predictive density ln p̂(y | x) = ln Σ_s W̄_s p(y | x, θ_s) at each unit's response, through the candidate seam:
    /// each case's grid holds its one observed response. The split form's density score, bit for bit.</summary>
    public static double[] Density(CandidateLogDensity units)
    {
        ArgumentNullException.ThrowIfNull(units);
        CheckResponses(units);
        (int[] rows, double[]? logWeights) = PredictiveChecks.Rows(units.Weights, units.Draws);
        var scores = new double[units.Cases];
        SplitConformal.LogPredictive(units.Values, rows, logWeights, 0, units.Cases, scores);
        return scores;
    }

    /// <summary>The standardized residual −|y − μ(x)|/ρ(x) at each unit (Lei et al.'s ρ; ρ = 1 without a scale), the split form's.</summary>
    public static double[] Residuals(ResidualReadings readings, IReadOnlyList<double> response)
    {
        ArgumentNullException.ThrowIfNull(readings);
        ArgumentNullException.ThrowIfNull(response);
        if (readings.Center.Count != response.Count || (readings.Scale is { } s && s.Count != response.Count))
            throw new ArgumentException("One centre, scale and response per unit is required.", nameof(readings));
        var scores = new double[response.Count];
        for (int i = 0; i < scores.Length; i++)
        {
            double scale = readings.Scale?[i] ?? 1.0;
            if (!(scale > 0.0 && double.IsFinite(scale))) throw new ArgumentException($"Unit {i}'s scale must be positive and finite.", nameof(readings));
            scores[i] = SplitConformal.ResidualScore(readings.Center[i], scale, response[i]);
        }
        return scores;
    }

    /// <summary>The PIT's centrality min(F̂, 1 − F̂) at each unit's response, F̂(y | x) = Σ_s W̄_s F(y | x, θ_s) read by both tails from the
    /// seam's log CDF, as the split form's predictive distribution reads it: a unit far in either tail of the fit's predictive conforms
    /// least. It orders units as the PIT's distance from ½ does, reversed.</summary>
    public static double[] Centrality(CandidateLogDensity units)
    {
        ArgumentNullException.ThrowIfNull(units);
        CheckResponses(units);
        double[,] logCdf = units.LogCdf ?? throw new ArgumentException("The units carry no log CDF; the PIT reads the fit's predictive CDF.", nameof(units));
        (int[] rows, double[]? logWeights) = PredictiveChecks.Rows(units.Weights, units.Draws);
        var cdf = new CdfScore[units.Cases];
        CdfScores.Fixed(logCdf, rows, CdfScores.Weights(rows, logWeights), 0, units.Cases, cdf);
        var scores = new double[units.Cases];
        for (int i = 0; i < scores.Length; i++) scores[i] = Math.Min(cdf[i].Lower, cdf[i].Upper);
        return scores;
    }

    private static void CheckResponses(CandidateLogDensity units)
    {
        for (int i = 0; i < units.Cases; i++)
            if (units.Grids[i].Length != 1)
                throw new ArgumentException($"Unit {i}'s grid holds {units.Grids[i].Length} values: each unit is read at its observed response.", nameof(units));
    }
}
