using System;
using System.Globalization;

namespace Maths.Samplers.Ensemble;

/// <summary>What a Gelman–Rubin R̂ reading is: an estimate, or one of the cases with no estimate.</summary>
public enum RHatState
{
    /// <summary>A finite estimate from positive within-chain variance.</summary>
    Estimated,

    /// <summary>Every draw of every chain is the same value: nothing to converge, and no estimate.</summary>
    Constant,

    /// <summary>No within-chain variance, but the chains sit at different values: not converged.</summary>
    Frozen,

    /// <summary>Fewer than two chains or two draws per chain, or a non-finite draw: no assessment possible.</summary>
    Unavailable,
}

/// <summary>
/// One functional's R̂ with its <see cref="RHatState"/>. <see cref="Value"/> is the estimate when
/// <see cref="RHatState.Estimated"/>, +∞ when <see cref="RHatState.Frozen"/>, and NaN otherwise. Only an estimate
/// at or below a target, or a constant functional, satisfies a stopping rule; frozen and unavailable readings never do.
/// </summary>
public readonly record struct RHatDiagnostic(RHatState State, double Value) : IFormattable
{
    public static RHatDiagnostic Unavailable => new(RHatState.Unavailable, double.NaN);
    public static RHatDiagnostic Constant => new(RHatState.Constant, double.NaN);
    public static RHatDiagnostic Frozen => new(RHatState.Frozen, double.PositiveInfinity);
    public static RHatDiagnostic Estimate(double value) => new(RHatState.Estimated, value);

    /// <summary>
    /// Whether this functional satisfies an R̂ stopping <paramref name="target"/>: an estimate at or below it, or a
    /// constant functional (exempt — it has nothing to converge).
    /// </summary>
    public bool MeetsTarget(double target)
        => State == RHatState.Constant || (State == RHatState.Estimated && Value <= target);

    /// <summary>
    /// The least favourable reading of a set — the summary a consensus map reports. Order: constant, then estimates
    /// by value, then frozen, then unavailable. An empty set is unavailable.
    /// </summary>
    public static RHatDiagnostic Worst(ReadOnlySpan<RHatDiagnostic> readings)
    {
        if (readings.IsEmpty) return Unavailable;
        RHatDiagnostic worst = readings[0];
        for (int i = 1; i < readings.Length; i++)
            if (WorseThan(readings[i], worst)) worst = readings[i];
        return worst;
    }

    private static bool WorseThan(RHatDiagnostic a, RHatDiagnostic b)
    {
        int ra = Rank(a.State), rb = Rank(b.State);
        if (ra != rb) return ra > rb;
        return a.State == RHatState.Estimated && a.Value > b.Value;
    }

    private static int Rank(RHatState state) => state switch
    {
        RHatState.Constant => 0,
        RHatState.Estimated => 1,
        RHatState.Frozen => 2,
        _ => 3,
    };

    public override string ToString() => ToString(null, CultureInfo.InvariantCulture);

    public string ToString(string? format, IFormatProvider? formatProvider) => State switch
    {
        RHatState.Estimated => Value.ToString(format, formatProvider),
        RHatState.Constant => "constant",
        RHatState.Frozen => "frozen",
        _ => "unavailable",
    };
}
