using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.LinAlg;

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

    /// <summary>Integrated out against its posterior given the observations outside the unit, so the term is a
    /// leave-unit-out predictive (<see cref="PointwiseContent.CollapsedPredictive"/>).</summary>
    IntegratedGivenOthers,
}

/// <summary>What each retained draw contributes to a pointwise density: OQ17's two draw contents. Each engine picks its
/// default and keeps the other selectable where it has one.</summary>
public enum PointwiseContent
{
    /// <summary>The draw's retained parameters: log p(y_unit | draw), the plain likelihood, with only the latents private
    /// to the unit integrated out. WAIC and PSIS-LOO both read it.</summary>
    RetainedParameters,

    /// <summary>
    /// The collapsed predictive. The draw's structure is kept, and its coefficients and noise scale are integrated out
    /// against their posterior given the observations outside the unit, so each term is the leave-unit-out predictive
    /// log p(y_unit | y_rest, structure). PSIS-LOO reads it as the conditional term of a non-factorized model (Bürkner,
    /// Gabry &amp; Vehtari 2021). WAIC does not read it, since its terms must condition on the unit's own data.
    /// Its importance ratios vary only with the structure. Where the structure is discrete, as BARS's knots are, the chain
    /// revisits a few values, and the Pareto fit reads that lumpiness as a heavy tail. So k̂ runs pessimistic while the
    /// estimate holds: on a BARS fit the elpd agreed with the retained content's to 0.006 nats per unit, with 30 units
    /// above k̂ = 0.5 against 5.
    /// </summary>
    CollapsedPredictive,
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

    /// <summary>Whether the terms are leave-unit-out predictives (some latent was <see cref="LatentTreatment.IntegratedGivenOthers"/>),
    /// as the collapsed content's are.</summary>
    public bool LeavesUnitsOut
    {
        get
        {
            foreach (LatentReport latent in Latents)
                if (latent.Treatment == LatentTreatment.IntegratedGivenOthers) return true;
            return false;
        }
    }
}

/// <summary>
/// The leave-unit-out predictive of a Gaussian model (OQ17's collapsed content), shared by the engines that integrate
/// their coefficients out in closed form. Take y | σ² ~ N(·, σ²Q⁻¹), with any location parameters already integrated
/// into the precision Q, so Q may be a projected precision. A unit U of m observations then has block Q_UU, and (Qr)_U for
/// the residual r. With S = rᵀQr over n_eff effective observations:
/// <list type="bullet">
///   <item>the leave-out residual is ẽ = Q_UU⁻¹(Qr)_U, and the conditional is y_U | y_rest, σ² ~ N(y_U − ẽ, σ²Q_UU⁻¹);</item>
///   <item>σ² | y_rest ~ IG(a + (n_eff − m)/2, b + S₋U/2), with S₋U = S − (Qr)_Uᵀẽ;</item>
///   <item>so ln p(y_U | y_rest) = lnΓ(a + n_eff/2) − lnΓ(a + (n_eff − m)/2) − (m/2) ln 2π + ½ ln|Q_UU|
///   + (a + (n_eff − m)/2) ln(b + S₋U/2) − (a + n_eff/2) ln(b + S/2),</item>
/// </list>
/// a multivariate Student t. a = b = 0 is the Jeffreys prior on σ². <see cref="Gaussian"/> is the known-scale form,
/// ln p = −(m/2) ln 2π + ½ ln|Q_UU| − ½(Qr)_Uᵀẽ.
/// </summary>
internal static class LeaveOutPredictive
{
    public static double StudentT(double[,] qUU, ReadOnlySpan<double> qr, double s, double effectiveCount,
                                  double priorShape = 0.0, double priorScale = 0.0)
    {
        int m = qr.Length;
        (double logDet, double quadratic) = Reduce(qUU, qr);
        double leaveOut = s - quadratic;
        double shapeRest = priorShape + 0.5 * (effectiveCount - m), shapeAll = priorShape + 0.5 * effectiveCount;
        if (!(shapeRest > 0.0)) throw new ArgumentException("The unit leaves too few observations to integrate the scale.", nameof(qr));
        if (!(leaveOut > 0.0)) throw new ArgumentException("The observations outside the unit leave no residual scale.", nameof(qr));
        return SpecialFunctions.LogGamma(shapeAll) - SpecialFunctions.LogGamma(shapeRest) - 0.5 * m * Math.Log(2.0 * Math.PI)
             + 0.5 * logDet + shapeRest * Math.Log(priorScale + 0.5 * leaveOut) - shapeAll * Math.Log(priorScale + 0.5 * s);
    }

    public static double Gaussian(double[,] qUU, ReadOnlySpan<double> qr)
    {
        (double logDet, double quadratic) = Reduce(qUU, qr);
        return -0.5 * qr.Length * Math.Log(2.0 * Math.PI) + 0.5 * logDet - 0.5 * quadratic;
    }

    // ln|Q_UU| and (Qr)_Uᵀ Q_UU⁻¹ (Qr)_U, through the block's Cholesky factor.
    private static (double LogDet, double Quadratic) Reduce(double[,] qUU, ReadOnlySpan<double> qr)
    {
        int m = qr.Length;
        if (qUU.GetLength(0) != m || qUU.GetLength(1) != m) throw new ArgumentException("The block must be m × m.", nameof(qUU));
        if (m == 1)
        {
            double q = qUU[0, 0];
            if (!(q > 0.0)) throw new InvalidOperationException("The unit's precision block is not positive definite.");
            return (Math.Log(q), qr[0] * qr[0] / q);
        }
        var factor = new DenseCholesky(m);
        factor.Decompose(qUU);
        if (!factor.Status.PositiveDefinite) throw new InvalidOperationException("The unit's precision block is not positive definite.");
        double[] leaveOut = factor.Solve(qr);
        double quadratic = 0.0;
        for (int a = 0; a < m; a++) quadratic += qr[a] * leaveOut[a];
        return (factor.LogDet, quadratic);
    }
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
