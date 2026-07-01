using System;
using System.Collections.Generic;
using Synthetic;

namespace Synthetic.Euclidean;

// ── Enums and geometry DTO ────────────────────────────────────────────────────

public enum TubeCrossSection
{
    Ribbon,                // flat rectangular: uniform in N×B rectangle; sharpest topology signal
    GaussianIsotropic,     // soft boundary, round cross-section
    GaussianAnisotropic,   // elliptical cross-section; major axis twists with strip
    UniformDisk,           // hard boundary, uniform density
    Annular,               // hollow tube; maximally stresses global topology
}

public enum MobiusPlacement
{
    /// <summary>Near θ≈0/2π seam — maximum Euclidean false-bridge stress.</summary>
    NearSeam,
    /// <summary>Ellipsoid at center-cross zone with semi-major axis orthogonal to the Möbius plane (along Z).</summary>
    CenterCrossOrtho,
    /// <summary>Near the peripheral turning elbow (θ≈¾π).</summary>
    PeripheralElbow,
    /// <summary>Ellipsoid at center-cross zone with semi-major axis coplanar with the Möbius plane.</summary>
    CenterCrossCoPlanar,
    /// <summary>Fully manual; use ellipsoidCenter and ellipsoidEulerXYZ directly.</summary>
    Manual,
}

public enum MobiusSpineShape
{
    /// <summary>Classic circular spine — the standard Möbius strip topology.</summary>
    Circle,
    /// <summary>Lissajous figure-eight spine (x=R·sin θ, y=R/2·sin 2θ): two loops connected at the origin.</summary>
    FigureEight,
}

public sealed class MobiusTubeGeometry : ClusterGeometry
{
    public double[][] SpineSamples { get; set; } = Array.Empty<double[]>();
    /// <summary>Per-spine-sample [T, N, B] frame triplet. Shape: M × 3 × 3.</summary>
    public double[][][] LocalFrames { get; set; } = Array.Empty<double[][]>();
    public double SpineRadius { get; set; }
    public double HalfWidth { get; set; }
    public double HalfThickness { get; set; }
    public double TwistCount { get; set; } = 1.0;
    public TubeCrossSection CrossSection { get; set; }
    public double RadialBias { get; set; } = 1.0;
}

// ── Generator ────────────────────────────────────────────────────────────────

/// <summary>
/// Two clusters with fundamentally different topologies:
///
///   Cluster 0 — a solid Möbius tube with non-orientable topology. In 3D it
///     self-intersects at the gluing seam (θ=0/2π), creating genuine false-neighbor
///     opportunities for Euclidean proximity graphs. Optional 4D embedding resolves
///     the self-intersection cleanly.
///
///   Cluster 1 — a compact anisotropic ellipsoid. Convex; a single full-covariance
///     Gaussian fits it faithfully.
/// Reference: custom Möbius tube + ellipsoid topology stress test benchmark.
/// </summary>
public static class MobiusAndEllipsoid
{
    public static SyntheticDataset Generate(
        // --- Möbius Tube ---
        int mobiusPoints = 5000,             // bumped from 1500 → plausible-data scale
        double spineRadius = 2.5,
        double halfWidth = 1.1,
        double halfThickness = 0.12,
        double noiseSigma = 0.06,
        double twistCount = 1.0,
        TubeCrossSection crossSection = TubeCrossSection.Annular,
        double radialBias = 1.0,
        MobiusSpineShape spineShape = MobiusSpineShape.FigureEight,
        double splayFactor = 0.7,
        // --- Ellipsoid ---
        int ellipsoidPoints = 2500,          // bumped from 700 → 7500 total
        double[]? ellipsoidAxes = null,
        EllipsoidShellMode ellipsoidShellMode = EllipsoidShellMode.Solid,
        // --- Placement ---
        MobiusPlacement placement = MobiusPlacement.CenterCrossOrtho,
        double[]? ellipsoidCenter = null,
        double[]? ellipsoidEulerXYZ = null,
        double intersectDepth = 0.0,
        double intersectRadialShift = 0.0,
        double gapScale = 1.0,
        // --- Ambient ---
        int dimensions = 3,
        int seed = 42)
    {
        ellipsoidAxes ??= new[] { 3.0, 0.5, 0.5 };
        if (ellipsoidAxes.Length != 3)
            throw new ArgumentException("ellipsoidAxes must have 3 elements.");
        if (dimensions < 3 || dimensions > 4)
            throw new ArgumentException("dimensions must be 3 or 4.");

        var rng = new Random(seed);
        int n = mobiusPoints + ellipsoidPoints;
        var features = new double[n][];
        var labels = new int[n];

        // ====================== Möbius Tube ======================
        var spineSamples = new List<double[]>();
        var localFrames = new List<double[][]>();

        double effBias = Math.Max(radialBias, 0.05);
        const double kDensityMax = 1.4;

        int p = 0;
        while (p < mobiusPoints)
        {
            double theta = 2.0 * Math.PI * rng.NextDouble();

            double cx, cy, sht;
            double[] T_frame, N_tw, B_tw;

            if (spineShape == MobiusSpineShape.Circle)
            {
                double thetaWeight = 1.0 + 0.4 * Math.Abs(Math.Cos(theta));
                if (rng.NextDouble() * kDensityMax > thetaWeight) continue;

                double ht = theta * twistCount * 0.5;
                double cht = Math.Cos(ht);
                sht = Math.Sin(ht);
                cx = spineRadius * Math.Cos(theta);
                cy = spineRadius * Math.Sin(theta);
                T_frame = new[] { -Math.Sin(theta), Math.Cos(theta), 0.0 };
                N_tw = new[] { cht * Math.Cos(theta), cht * Math.Sin(theta), sht };
                B_tw = new[] { -sht * Math.Cos(theta), -sht * Math.Sin(theta), cht };
            }
            else // FigureEight
            {
                double txr = spineRadius * Math.Cos(theta);
                double tyr = spineRadius * Math.Cos(2.0 * theta);
                double tLen = Math.Sqrt(txr * txr + tyr * tyr) + 1e-12;
                cx = spineRadius * Math.Sin(theta);
                cy = spineRadius / 2.0 * Math.Sin(2.0 * theta);
                double Nfx = -tyr / tLen, Nfy = txr / tLen;
                double ht = theta * twistCount * 0.5;
                double cht = Math.Cos(ht);
                sht = Math.Sin(ht);
                T_frame = new[] { txr / tLen, tyr / tLen, 0.0 };
                N_tw = new[] { cht * Nfx, cht * Nfy, sht };
                B_tw = new[] { -sht * Nfx, -sht * Nfy, cht };
            }

            double effHalfWidth = (spineShape == MobiusSpineShape.FigureEight)
                ? halfWidth * effBias * (1.0 - splayFactor * (1.0 - Math.Abs(Math.Sin(theta))))
                : halfWidth * effBias;
            double effHalfThick = halfThickness * effBias;

            double px, py, pz, w4;
            if (crossSection == TubeCrossSection.Ribbon)
            {
                double uN = (rng.NextDouble() * 2.0 - 1.0) * effHalfWidth;
                double uB = (rng.NextDouble() * 2.0 - 1.0) * effHalfThick;
                px = cx + uN * N_tw[0] + uB * B_tw[0];
                py = cy + uN * N_tw[1] + uB * B_tw[1];
                pz = uN * N_tw[2] + uB * B_tw[2];
                w4 = uN * sht * 1.2;
            }
            else
            {
                double rN = SampleCrossSection(rng, effHalfWidth, crossSection);
                double rB = SampleCrossSection(rng, effHalfThick, crossSection);
                double phi = 2.0 * Math.PI * rng.NextDouble();
                double cosPhi = Math.Cos(phi);
                double sinPhi = Math.Sin(phi);
                px = cx + rN * cosPhi * N_tw[0] + rB * sinPhi * B_tw[0];
                py = cy + rN * cosPhi * N_tw[1] + rB * sinPhi * B_tw[1];
                pz = rN * cosPhi * N_tw[2] + rB * sinPhi * B_tw[2];
                w4 = rN * sinPhi * sht * 1.2;
            }

            px += noiseSigma * SyntheticData.SampleStandardNormal(rng);
            py += noiseSigma * SyntheticData.SampleStandardNormal(rng);
            pz += noiseSigma * SyntheticData.SampleStandardNormal(rng);

            var point = new double[dimensions];
            point[0] = px; point[1] = py; point[2] = pz;

            if (dimensions == 4)
                point[3] = w4;

            features[p] = point;
            labels[p] = 0;
            p++;
        }

        // Spine samples on a uniform θ grid
        int spineRes = Math.Max(80, mobiusPoints / 5);
        for (int s = 0; s < spineRes; s++)
        {
            double spTheta = 2.0 * Math.PI * s / spineRes;
            double spHalf = spTheta * twistCount * 0.5;
            double spCht = Math.Cos(spHalf), spSht = Math.Sin(spHalf);
            if (spineShape == MobiusSpineShape.Circle)
            {
                spineSamples.Add(new[] { spineRadius * Math.Cos(spTheta), spineRadius * Math.Sin(spTheta), 0.0 });
                localFrames.Add(new[] {
                    new[] { -Math.Sin(spTheta), Math.Cos(spTheta), 0.0 },
                    new[] { spCht * Math.Cos(spTheta), spCht * Math.Sin(spTheta), spSht },
                    new[] { -spSht * Math.Cos(spTheta), -spSht * Math.Sin(spTheta), spCht },
                });
            }
            else // FigureEight
            {
                double spTxr = spineRadius * Math.Cos(spTheta);
                double spTyr = spineRadius * Math.Cos(2.0 * spTheta);
                double spTLen = Math.Sqrt(spTxr * spTxr + spTyr * spTyr) + 1e-12;
                double spTx = spTxr / spTLen, spTy = spTyr / spTLen;
                double spNx = -spTy, spNy = spTx;
                spineSamples.Add(new[] {
                    spineRadius * Math.Sin(spTheta),
                    spineRadius / 2.0 * Math.Sin(2.0 * spTheta),
                    0.0
                });
                localFrames.Add(new[] {
                    new[] { spTx, spTy, 0.0 },
                    new[] { spCht * spNx, spCht * spNy, spSht },
                    new[] { -spSht * spNx, -spSht * spNy, spCht },
                });
            }
        }

        // ====================== Ellipsoid ======================
        double[] center = ResolveCenter(
            placement, spineRadius, halfWidth,
            intersectDepth, intersectRadialShift, gapScale, ellipsoidCenter, spineShape, rng);

        double[,] rot = ResolveRotation(placement, ellipsoidEulerXYZ, rng);
        double[,] cov = SyntheticData.BuildCovariance(rot, ellipsoidAxes);
        double[,] L = SyntheticData.CholeskyLower(cov, 3);

        for (int ep = 0; ep < ellipsoidPoints; ep++)
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
            for (int d = 0; d < 3; d++) z[d] = z[d] / norm * r;

            double[] pt = SyntheticData.MultiplyMatrixVector(L, z);
            var point = new double[dimensions];
            for (int d = 0; d < 3; d++) point[d] = center[d] + pt[d];
            if (dimensions == 4) point[3] = ellipsoidAxes[2] * SyntheticData.SampleStandardNormal(rng);

            features[mobiusPoints + ep] = point;
            labels[mobiusPoints + ep] = 1;
        }

        var mobiusGeom = new MobiusTubeGeometry
        {
            SpineSamples = spineSamples.ToArray(),
            LocalFrames = localFrames.ToArray(),
            SpineRadius = spineRadius,
            HalfWidth = halfWidth,
            HalfThickness = halfThickness,
            TwistCount = twistCount,
            CrossSection = crossSection,
            RadialBias = radialBias,
        };

        var ellipsoidGeom = new EllipsoidGeometry { Center = center, Covariance = cov };

        return new SyntheticDataset
        {
            Features = features,
            Labels = labels,
            ClusterCount = 2,
            LabelsByLevel = new[] { labels },
            Parameters = new Dictionary<string, object>
            {
                ["generator"] = "MobiusAndEllipsoid",
                ["dimensions"] = dimensions,
                ["mobiusPoints"] = mobiusPoints,
                ["ellipsoidPoints"] = ellipsoidPoints,
                ["spineRadius"] = spineRadius,
                ["halfWidth"] = halfWidth,
                ["halfThickness"] = halfThickness,
                ["noiseSigma"] = noiseSigma,
                ["twistCount"] = twistCount,
                ["crossSection"] = crossSection.ToString(),
                ["radialBias"] = radialBias,
                ["spineShape"] = spineShape.ToString(),
                ["splayFactor"] = splayFactor,
                ["placement"] = placement.ToString(),
                ["intersectDepth"] = intersectDepth,
                ["intersectRadialShift"] = intersectRadialShift,
                ["gapScale"] = gapScale,
                ["ellipsoidAxes"] = ellipsoidAxes,
                ["ellipsoidShellMode"] = ellipsoidShellMode.ToString(),
                ["seed"] = seed,
            },
            Metadata = new SyntheticDatasetMeta(
                GeneratorName: nameof(MobiusAndEllipsoid),
                GeometryClass: "Euclidean",
                TopologyTag: "topological",
                HierarchyTag: "none",
                GTNumClusters: 2,
                AmbientDimensionality: dimensions,
                LiteratureReference: "custom Möbius tube + ellipsoid topology stress test benchmark"),
            ClusterGeometries = new ClusterGeometry[] { mobiusGeom, ellipsoidGeom },
        };
    }

    /// <summary>
    /// Projects 4D points to 3D by rotating in the x₃–x₄ plane.
    /// Sweep <paramref name="angle4"/> from 0 → 2π for an animated 4D rotation.
    /// Points with fewer than 4 dimensions are passed through unchanged.
    /// </summary>
    public static double[][] Project4DTo3D(double[][] pts4D, double angle4 = 0.0)
    {
        double c = Math.Cos(angle4), s = Math.Sin(angle4);
        var proj = new double[pts4D.Length][];
        for (int i = 0; i < pts4D.Length; i++)
        {
            var pt = pts4D[i];
            if (pt.Length < 4) { proj[i] = (double[])pt.Clone(); continue; }
            proj[i] = new[] { pt[0], pt[1], pt[2] * c + pt[3] * s };
        }
        return proj;
    }

    // ── Placement helpers ─────────────────────────────────────────────────────

    private static double[] ResolveCenter(
        MobiusPlacement placement, double R, double w,
        double intersectDepth, double intersectRadialShift, double gapScale,
        double[]? manualCenter, MobiusSpineShape spineShape, Random rng)
    {
        if (placement == MobiusPlacement.Manual)
        {
            if (manualCenter is null || manualCenter.Length != 3)
                throw new ArgumentException(
                    "Manual placement requires a 3-element ellipsoidCenter.");
            return (double[])manualCenter.Clone();
        }

        bool isFigureEight = spineShape == MobiusSpineShape.FigureEight;

        switch (placement)
        {
            case MobiusPlacement.NearSeam:
                return new[] { R + intersectRadialShift, 0.0, intersectDepth };

            case MobiusPlacement.CenterCrossOrtho:
                if (isFigureEight)
                    return new[] { intersectRadialShift, 0.0, intersectDepth };
                return new[] { -R + intersectRadialShift, 0.0, intersectDepth };

            case MobiusPlacement.PeripheralElbow:
                if (isFigureEight)
                {
                    double sq2 = Math.Sqrt(2.0);
                    return new[] {
                        R / sq2 + intersectRadialShift,
                        R / 2.0,
                        intersectDepth * 0.8
                    };
                }
                double alpha = Math.PI * 0.75;
                return new[] {
                    R * Math.Cos(alpha) + intersectRadialShift,
                    R * Math.Sin(alpha),
                    intersectDepth * 0.8
                };

            case MobiusPlacement.CenterCrossCoPlanar:
                return new[] { intersectRadialShift, 0.0, intersectDepth };

            default:
                return new[] { R * 0.2, -R * 0.9, 0.8 };
        }
    }

    private static double[,] ResolveRotation(
        MobiusPlacement placement, double[]? eulerXYZ, Random rng)
    {
        if (eulerXYZ is not null)
        {
            if (eulerXYZ.Length != 3)
                throw new ArgumentException("ellipsoidEulerXYZ must have 3 elements.");
            return SyntheticData.EulerToRotation(eulerXYZ[0], eulerXYZ[1], eulerXYZ[2]);
        }

        return placement switch
        {
            MobiusPlacement.NearSeam
                => SyntheticData.EulerToRotation(0.0, Math.PI / 2.0, 0.0),
            MobiusPlacement.CenterCrossOrtho
                => SyntheticData.EulerToRotation(0.0, Math.PI / 2.0, 0.0),
            MobiusPlacement.PeripheralElbow
                => SyntheticData.EulerToRotation(0.0, 0.0, Math.PI * 0.75),
            MobiusPlacement.CenterCrossCoPlanar
                => SyntheticData.EulerToRotation(0.0, 0.0, Math.PI / 4.0),
            _ => SyntheticData.EulerToRotation(0.0, Math.PI / 2.0, 0.0),
        };
    }

    private static double SampleCrossSection(
        Random rng, double scale, TubeCrossSection mode)
    {
        return mode switch
        {
            TubeCrossSection.GaussianIsotropic or TubeCrossSection.GaussianAnisotropic
                => Math.Abs(SyntheticData.SampleStandardNormal(rng)) * scale,
            TubeCrossSection.UniformDisk
                => scale * Math.Sqrt(rng.NextDouble()),
            TubeCrossSection.Annular
                => scale * (0.4 + 0.6 * rng.NextDouble()),
            _ => scale * Math.Sqrt(rng.NextDouble()),
        };
    }
}


