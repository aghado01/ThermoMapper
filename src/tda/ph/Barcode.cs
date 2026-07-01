#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// One bar in a persistence barcode: a topological feature born at
/// <see cref="Birth"/> and dying at <see cref="Death"/> along the sweep axis.
/// <para><see cref="Dimension"/>: 0 = connected component, 1 = loop, 2 = void.</para>
/// <para><see cref="Generator"/>: optional opaque identifier for the filtration
/// element that gave rise to this bar; reserved for viewer back-references and
/// (co)cycle representatives.</para>
/// <para><see cref="Cocycle"/>: when requested from <see cref="PersistentCohomology"/>,
/// indices of filtration simplices forming a representative cocycle (Z/2).</para>
/// <para><see cref="Cycle"/>: when requested from <see cref="PersistentInvolutedHomology"/>,
/// indices of filtration simplices forming a representative homology cycle (Z/2).
/// Use <see cref="BarCycleEdges"/> to obtain LMP-consumable vertex pairs for H1.</para>
/// </summary>
public readonly record struct Bar(
    double Birth,
    double Death,
    int Dimension,
    int? Generator = null,
    int[]? Cocycle = null,
    int[]? Cycle = null,
    IntervalEnd BirthEnd = IntervalEnd.Closed,
    IntervalEnd DeathEnd = IntervalEnd.Closed)
{
    public double Persistence =>
        double.IsPositiveInfinity(Death) ? double.PositiveInfinity : Death - Birth;

    public bool IsInfinite => double.IsPositiveInfinity(Death);
}

/// <summary>
/// Persistence diagram output currency — the narrow waist both PH paths emit and
/// every consumer reads (viz, diagnostics, diagram distances, SPRED cost).
/// </summary>
public sealed class Barcode
{
    public IReadOnlyList<Bar> Bars { get; }
    public string AxisLabel { get; }

    public Barcode(IReadOnlyList<Bar> bars, string axisLabel = "")
    {
        ArgumentNullException.ThrowIfNull(bars);
        Bars = bars;
        AxisLabel = axisLabel;
    }
}
