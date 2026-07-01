using System;
using System.Collections.Generic;
using Synthetic;

namespace Synthetic.Euclidean;

// ── Enums ────────────────────────────────────────────────────────────────────

public enum EllipsoidPlacement
{
    /// <summary>Ellipsoid near the crescent's open face (original default).
    /// Stresses proximity rules at the open gap.</summary>
    NearOpenFace,

    /// <summary>Ellipsoid semi-major axis aligned with Z (perpendicular to the
    /// crescent plane), center placed at the crescent elbow. The canonical stress
    /// test for manifold-aware metrics: Euclidean distance fuses the clusters at
    /// the intersection while geodesic / graph-based distances correctly separate them.</summary>
    OrthogonalElbowIntersect,

    /// <summary>Ellipsoid centered at the upper arc tip (angle = π − arcHalfAngle).
    /// Semi-major axis aligned outward along the radial direction at that tip.</summary>
    IntersectUpperTip,

    /// <summary>Ellipsoid centered at the lower arc tip (angle = π + arcHalfAngle).
    /// Mirror of IntersectUpperTip about the XZ plane.</summary>
    IntersectLowerTip,

    /// <summary>Fully manual: use ellipsoidCenter and ellipsoidEulerXYZ directly.</summary>
    Manual,
}

/// <summary>Controls the radial density profile used when sampling the ellipsoid cluster.</summary>
public enum EllipsoidShellMode
{
    /// <summary>Uniform volume density: cbrt-CDF radial scaling. Hard boundary.</summary>
    Solid,

    /// <summary>True multivariate Gaussian: raw standard-normal vector without radial scaling.</summary>
    Gaussian,

    /// <summary>Surface shell only: all points on the ellipsoid boundary (r = 1 before Cholesky warp).</summary>
    Hollow,

    /// <summary>Thick shell: inner radius 0.75, outer radius 1.0 (uniform in that band).</summary>
    Annular,
}

// ── Generator ────────────────────────────────────────────────────────────────

/// <summary>
/// Two clusters with fundamentally different topologies:
///
///   Cluster 0 — a crescent (half-annulus arc) in 3D ambient space.
///     Non-convex; no single Gaussian can represent it faithfully.
///     Exposes false-bridge formation in Euclidean k-NN graphs.
///
///   Cluster 1 — a compact anisotropic ellipsoid. Convex; a single
///     full-covariance Gaussian fits it well.
///
/// v2 adds full positional and orientational control. See <see cref="EllipsoidPlacement"/>
/// for preset arrangements targeting specific stress scenarios.
/// Reference: custom crescent-and-ellipsoid topology stress test benchmark.
/// </summary>
public static class CrescentAndEllipsoid
{
    public static SyntheticDataset Generate(
        // --- Crescent ---
        int crescentPoints = 5000,           // bumped from 1500 → plausible-data scale
        double crescentRadius = 3.0,
        double crescentWidth = 0.40,
        double arcHalfAngle = 2.04,
        // --- Ellipsoid ---
        int ellipsoidPoints = 2500,          // bumped from 700 → 7500 total
        double[]? ellipsoidAxes = null,
        EllipsoidShellMode ellipsoidShellMode = EllipsoidShellMode.Hollow,
        // --- Relative placement ---
        EllipsoidPlacement placement = EllipsoidPlacement.NearOpenFace,
        double[]? ellipsoidCenter = null,
        double[]? ellipsoidEulerXYZ = null,
        double intersectDepth = 0.0,
        double intersectRadialShift = 0.0,
        double gapScale = 1.0,
        int seed = 42)
    {
        ellipsoidAxes ??= new[] { 3.0, 0.5, 0.5 };

        if (ellipsoidAxes.Length != 3)
            throw new ArgumentException("ellipsoidAxes must have 3 elements.");

        var rng = new Random(seed);
        int n = crescentPoints + ellipsoidPoints;
        var features = new double[n][];
        var labels = new int[n];

        // ------------------------------------------------------------------
        // Cluster 0: crescent arc
        // ------------------------------------------------------------------
        double arcStart = Math.PI - arcHalfAngle;
        double arcEnd = Math.PI + arcHalfAngle;

        double elbowX = crescentRadius * Math.Cos(Math.PI);
        double elbowY = crescentRadius * Math.Sin(Math.PI);
        double elbowZ = 0.0;

        int spineRes = Math.Max(64, crescentPoints / 4);
        var spineSamples = new double[spineRes][];
        for (int s = 0; s < spineRes; s++)
        {
            double t = arcStart + (arcEnd - arcStart) * s / (spineRes - 1);
            spineSamples[s] = new[] {
                crescentRadius * Math.Cos(t),
                crescentRadius * Math.Sin(t),
                0.0
            };
        }

        for (int p = 0; p < crescentPoints; p++)
        {
            double u = SyntheticData.SampleStandardNormal(rng) * 0.35;
            u = Math.Max(-0.5, Math.Min(0.5, u));
            double t = Math.PI + u * (arcEnd - arcStart);

            double cos = Math.Cos(t), sin = Math.Sin(t);
            double nx = cos, ny = sin;
            double txArc = -sin, tyArc = cos;

            double normalizedAngle = u * 2.0;
            double taper = Math.Cos(normalizedAngle * Math.PI / 2.0);
            double localWidth = crescentWidth * taper;

            double nr = localWidth * SyntheticData.SampleStandardNormal(rng);
            double nt = crescentWidth * 0.12 * SyntheticData.SampleStandardNormal(rng);
            double nz = localWidth * SyntheticData.SampleStandardNormal(rng);

            features[p] = new[] {
                crescentRadius * cos + nr * nx + nt * txArc,
                crescentRadius * sin + nr * ny + nt * tyArc,
                nz
            };
            labels[p] = 0;
        }

        // ------------------------------------------------------------------
        // Cluster 1: anisotropic ellipsoid
        // ------------------------------------------------------------------
        double[] center = ResolveCenter(
            placement, ellipsoidCenter, ellipsoidAxes,
            crescentRadius, crescentWidth, arcHalfAngle, gapScale,
            elbowX, elbowY, elbowZ,
            intersectDepth, intersectRadialShift);

        double[,] rot = ResolveRotation(
            placement, ellipsoidEulerXYZ, arcHalfAngle, rng);

        double[,] cov = SyntheticData.BuildCovariance(rot, ellipsoidAxes);
        double[,] L = SyntheticData.CholeskyLower(cov, 3);

        for (int p = 0; p < ellipsoidPoints; p++)
        {
            double[] z = {
                SyntheticData.SampleStandardNormal(rng),
                SyntheticData.SampleStandardNormal(rng),
                SyntheticData.SampleStandardNormal(rng)
            };
            double norm = Math.Sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
            if (norm < 1e-12) norm = 1.0;

            double r = ellipsoidShellMode switch
            {
                EllipsoidShellMode.Gaussian => norm,
                EllipsoidShellMode.Solid => Math.Pow(rng.NextDouble(), 1.0 / 3.0),
                EllipsoidShellMode.Hollow => 1.0,
                EllipsoidShellMode.Annular => 0.75 + 0.25 * rng.NextDouble(),
                _ => Math.Pow(rng.NextDouble(), 1.0 / 3.0),
            };

            z[0] = z[0] / norm * r;
            z[1] = z[1] / norm * r;
            z[2] = z[2] / norm * r;

            double[] pt = SyntheticData.MultiplyMatrixVector(L, z);
            features[crescentPoints + p] = new[] {
                center[0] + pt[0],
                center[1] + pt[1],
                center[2] + pt[2]
            };
            labels[crescentPoints + p] = 1;
        }

        var arcGeom = new ArcGeometry
        {
            SpineSamples = spineSamples,
            Radius = crescentRadius,
            NoiseScale = crescentWidth,
        };

        var ellipsoidGeom = new EllipsoidGeometry
        {
            Center = center,
            Covariance = cov,
        };

        return new SyntheticDataset
        {
            Features = features,
            Labels = labels,
            ClusterCount = 2,
            LabelsByLevel = new[] { labels },
            Parameters = new Dictionary<string, object>
            {
                ["generator"] = "CrescentAndEllipsoid",
                ["crescentPoints"] = crescentPoints,
                ["ellipsoidPoints"] = ellipsoidPoints,
                ["crescentRadius"] = crescentRadius,
                ["crescentWidth"] = crescentWidth,
                ["arcHalfAngle"] = arcHalfAngle,
                ["ellipsoidAxes"] = ellipsoidAxes,
                ["placement"] = placement.ToString(),
                ["ellipsoidCenter"] = center,
                ["gapScale"] = gapScale,
                ["intersectDepth"] = intersectDepth,
                ["intersectRadialShift"] = intersectRadialShift,
                ["ellipsoidShellMode"] = ellipsoidShellMode.ToString(),
                ["seed"] = seed,
            },
            Metadata = new SyntheticDatasetMeta(
                GeneratorName: nameof(CrescentAndEllipsoid),
                GeometryClass: "Euclidean",
                TopologyTag: "non-convex",
                HierarchyTag: "none",
                GTNumClusters: 2,
                AmbientDimensionality: 3,
                LiteratureReference: "custom crescent-and-ellipsoid topology stress test benchmark"),
            ClusterGeometries = new ClusterGeometry[] { arcGeom, ellipsoidGeom },
        };
    }

    // ── Placement helpers ─────────────────────────────────────────────────────

    private static double[] ResolveCenter(
        EllipsoidPlacement placement,
        double[]? manualCenter,
        double[] axes,
        double crescentRadius,
        double crescentWidth,
        double arcHalfAngle,
        double gapScale,
        double elbowX, double elbowY, double elbowZ,
        double intersectDepth,
        double intersectRadialShift)
    {
        switch (placement)
        {
            case EllipsoidPlacement.Manual:
                if (manualCenter is null || manualCenter.Length != 3)
                    throw new ArgumentException(
                        "ellipsoidCenter must be 3-element when placement == Manual.");
                return (double[])manualCenter.Clone();

            case EllipsoidPlacement.OrthogonalElbowIntersect:
                return new[] {
                    elbowX + intersectRadialShift,
                    elbowY,
                    elbowZ + intersectDepth
                };

            case EllipsoidPlacement.IntersectUpperTip:
            case EllipsoidPlacement.IntersectLowerTip:
                {
                    double ySgn = placement == EllipsoidPlacement.IntersectUpperTip ? 1.0 : -1.0;
                    double radX = -Math.Cos(arcHalfAngle);
                    double radY = ySgn * Math.Sin(arcHalfAngle);
                    double spine = crescentRadius + intersectRadialShift;
                    return new[] {
                        radX * spine,
                        radY * spine,
                        elbowZ + intersectDepth
                    };
                }

            case EllipsoidPlacement.NearOpenFace:
            default:
                return new[] {
                    0.0 + intersectRadialShift,
                    0.0,
                    intersectDepth
                };
        }
    }

    private static double[,] ResolveRotation(
        EllipsoidPlacement placement,
        double[]? eulerXYZ,
        double arcHalfAngle,
        Random rng)
    {
        if (eulerXYZ is not null)
        {
            if (eulerXYZ.Length != 3)
                throw new ArgumentException("ellipsoidEulerXYZ must be 3-element.");
            return SyntheticData.EulerToRotation(eulerXYZ[0], eulerXYZ[1], eulerXYZ[2]);
        }

        switch (placement)
        {
            case EllipsoidPlacement.OrthogonalElbowIntersect:
                return SyntheticData.EulerToRotation(0.0, Math.PI / 2.0, 0.0);

            case EllipsoidPlacement.IntersectUpperTip:
                return SyntheticData.EulerToRotation(0.0, 0.0, Math.PI - arcHalfAngle);

            case EllipsoidPlacement.IntersectLowerTip:
                return SyntheticData.EulerToRotation(0.0, 0.0, Math.PI + arcHalfAngle);

            case EllipsoidPlacement.NearOpenFace:
            default:
                return SyntheticData.EulerToRotation(0.0, 0.0, Math.PI / 2.0);
        }
    }
}


