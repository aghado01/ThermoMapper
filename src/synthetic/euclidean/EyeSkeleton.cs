using System;
using System.Collections.Generic;

namespace Synthetic.Euclidean;

/// <summary>Cross-section shape of a stroke's tube, sampled in the normal plane.</summary>
public enum CrossSectionShape
{
    /// <summary>Hollow tube surface (the boundary ellipse) — the canonical thin ring.</summary>
    Shell,

    /// <summary>Filled disk — a solid "cylindrical donut".</summary>
    Solid,

    /// <summary>Flat band: wide in-plane, near-zero out-of-plane — the 2-D-toy-as-ribbon look.</summary>
    Ribbon,
}

/// <summary>
/// One stroke of the eye: a geodesic arc (a geodesic circle of radius
/// <see cref="MajorRadius"/> about the skeleton center, swept over
/// [<see cref="ArcStart"/>, <see cref="ArcEnd"/>]) carrying a tube of
/// cross-section radius <see cref="MinorRadius"/>. Coordinate-free — the
/// realizer turns it into points on whichever manifold it is handed, so the
/// same stroke serves the Euclidean and (later) Poincaré-ball eyes.
/// </summary>
public sealed class EyeStroke
{
    /// <summary>Geodesic radius of the swept circle about the skeleton center.</summary>
    public required double MajorRadius { get; init; }

    /// <summary>Cross-section (tube) radius.</summary>
    public required double MinorRadius { get; init; }

    public double ArcStart { get; init; } = 0.0;
    public double ArcEnd { get; init; } = 2.0 * Math.PI;

    /// <summary>Tube z-extent: 0 = a flat ribbon in-plane, 1 = a full round tube.</summary>
    public double ZThickness { get; init; } = 0.3;

    /// <summary>Stroke center offset relative to the skeleton center (asymmetry).</summary>
    public double[] Offset { get; init; } = { 0.0, 0.0, 0.0 };

    /// <summary>Density bias toward the tube's inner edge (0 = uniform, 1 = inner).
    /// Honored by Shell/Solid cross-sections; ignored by Ribbon.</summary>
    public double DensityGradientStrength { get; init; } = 0.0;

    /// <summary>Cross-section shape (Shell hollow / Solid filled / Ribbon flat band).</summary>
    public CrossSectionShape CrossSection { get; init; } = CrossSectionShape.Shell;

    /// <summary>Out-of-plane half-thickness for <see cref="CrossSectionShape.Ribbon"/>.</summary>
    public double RibbonThickness { get; init; } = 0.04;

    /// <summary>
    /// Continuous cross-section taper along the arc, in [0,1]: 0 = uniform tube,
    /// 1 = sin(π·t) almond (zero width at the arc ends, max in the middle) — the
    /// continuously-varying brow / eye-bag irregularity. Meaningful only for
    /// partial arcs; leave 0 on the full central ring.
    /// </summary>
    public double EndTaper { get; init; } = 0.0;

    // ── Ambient frame ────────────────────────────────────────────────────────
    // The geodesic circle is swept in the 2-plane span(e[SweepAxisA], e[SweepAxisB]);
    // the tube's out-of-plane offset goes along e[AxialAxis]. All three must be
    // distinct and < the skeleton's ambient dimension. Defaults reproduce the
    // canonical 3-D eye (sweep in XY, thickness in Z). In 4-D, giving strokes
    // distinct sweep planes places their rings in independent 2-planes — circles
    // are unlinked in codimension ≥ 3, removing the ambient false-proximity that
    // would otherwise bridge them.

    public int SweepAxisA { get; init; } = 0;
    public int SweepAxisB { get; init; } = 1;
    public int AxialAxis { get; init; } = 2;

    public required int PointCount { get; init; }
    public required int Label { get; init; }
}

/// <summary>
/// A solid geodesic ball at the skeleton center — the pupil that fits in the
/// iris void, dilatable up toward the inner bore.
/// </summary>
public sealed class EyePupil
{
    public required double Radius { get; init; }
    public required int PointCount { get; init; }
    public required int Label { get; init; }
}

/// <summary>
/// Intrinsic, coordinate-free description of the eye — a center point and a
/// set of geodesic-arc <see cref="EyeStroke"/>s (plus an optional
/// <see cref="Pupil"/>) — plus the realization knobs a curved realizer honors
/// (a flat realizer leaves the curvature-only ones inert). The skeleton is the
/// shared field both the Euclidean and hyperbolic generators reduce to points;
/// only the realizing manifold differs.
/// </summary>
public sealed class EyeSkeleton
{
    /// <summary>Skeleton center. Default origin — for the Poincaré ball this is
    /// the ball center, where the warp is radially symmetric and most readable.</summary>
    public double[] Center { get; init; } = { 0.0, 0.0, 0.0 };

    /// <summary>
    /// Geodesic shell cap: structure is confined within this geodesic radius of
    /// <see cref="Center"/>, keeping a curved realization inside an inner shell
    /// (away from the boundary blow-up). In the flat realizer this is an ordinary
    /// Euclidean-radius clamp. Default ∞ — no cap.
    /// </summary>
    public double MaxGeodesicRadius { get; init; } = double.PositiveInfinity;

    /// <summary>
    /// Realization-fidelity dial in [0,1]: 1 = geometry-faithful (full curvature),
    /// 0 = cosmetically Euclidean (curvature cancelled so a curved render matches
    /// the flat eye). Active only in a curved realizer; inert when the manifold is
    /// flat. The pedagogical control for the distortion demo.
    /// </summary>
    public double WarpStrength { get; init; } = 1.0;

    public required IReadOnlyList<EyeStroke> Strokes { get; init; }

    /// <summary>Optional central pupil (a solid ball in the iris void). Null = absent.</summary>
    public EyePupil? Pupil { get; init; }
}
