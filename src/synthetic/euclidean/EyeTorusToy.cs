using System;
using System.Collections.Generic;
using Maths.Rng;
using Maths.Samplers;
using Synthetic;

namespace Synthetic.Euclidean;

/// <summary>
/// "Eye" toy in arbitrary ambient dimension (default 3-D): a central full torus
/// (iris ring) plus upper/lower half-tori (asymmetric arcs) over a noise
/// blanket — a sampler that lifts the canonical Blatt–Wiseman–Domany "Toy"
/// (BWD 1996, PRL 76:3251, Fig. 1) into D-space. The intrinsic structure is an
/// <see cref="EyeSkeleton"/>; <see cref="SampleLocal"/> samples it into
/// tangent-space (eye-local) vectors of dimension <c>Center.Length</c>, which
/// this generator realizes on flat space and <c>HyperbolicEyeTorus</c> realizes
/// on the Poincaré ball. In 4-D the strokes occupy independent 2-planes, so the
/// rings are unlinked (no ambient false-bridges) — a controlled test of
/// ambient-vs-intrinsic clustering. <see cref="FlattenToPlane"/> recovers the
/// 2-D plate by orthographic projection onto the first two axes.
/// </summary>
/// <remarks>
/// Two-level ground truth via <c>LabelsByLevel</c>: level 0 = signal vs.
/// background (0/1); level 1 = structure groups then background. Cross-section
/// profiles (Shell/Solid/Ribbon) and the brow/eye-bag taper live in the
/// tangent-space layer and carry to the hyperbolic realizer. The colored-noise
/// FFT blanket is 3-D-only; for other dimensions the background is uniform
/// (N-D colored noise is a separate lift). Default D=3 with default stroke
/// frames reproduces the canonical eye exactly.
/// </remarks>
public static class EyeTorusToy
{
    public sealed class EyeTorusToyConfig
    {
        // Core geometry (BWD-figure scale, radius ~1–5).
        public double CentralMajorR { get; set; } = 2.5;
        public double CentralMinorR { get; set; } = 0.6;
        public double UpperMajorR { get; set; } = 2.0;   // smaller
        public double UpperMinorR { get; set; } = 0.4;   // thinner
        public double LowerMajorR { get; set; } = 3.8;   // larger
        public double LowerMinorR { get; set; } = 0.9;   // fatter

        public double ZThickness { get; set; } = 0.3;    // 0 = flat rings, 1 = full circular tubes

        // Angular spans for the arcs (radians; central is the full 2π).
        public double UpperArcStart { get; set; } = -Math.PI * 0.6;
        public double UpperArcEnd { get; set; } = Math.PI * 0.6;
        public double LowerArcStart { get; set; } = Math.PI * 0.7;
        public double LowerArcEnd { get; set; } = Math.PI * 2.3;
        public double[] LowerOffset { get; set; } = { 0.0, -0.5, 0.0 }; // asymmetry

        // Sampling.
        public int CentralPoints { get; set; } = 4000;
        public int UpperPoints { get; set; } = 2500;
        public int LowerPoints { get; set; } = 3500;
        public double StructureNoiseSigma { get; set; } = 0.08; // isotropic per-point jitter

        // Cross-section profiles (tangent-space; carry to the hyperbolic realizer).
        public CrossSectionShape CrossSection { get; set; } = CrossSectionShape.Shell;     // all strokes
        public CrossSectionShape? CentralCrossSection { get; set; }                        // override central (e.g. Solid donut)
        public double RibbonThickness { get; set; } = 0.04;
        public double HalfArcTaper { get; set; } = 0.0;  // brow/eye-bag almond taper on the half-tori

        // Pupil (solid ball in the iris void; dilation 0 → point, 1 → fills the bore R_central − r_central).
        public int PupilPoints { get; set; }             // 0 = no pupil
        public double PupilDilation { get; set; } = 0.6;

        // Background blanket (colored-noise density field; 3-D only).
        public double BackgroundDensityRatio { get; set; } = 0.12; // bg points / structure points
        public double SpectralExponent { get; set; } = SpectralNoiseField.Pink; // β: 0 white, 1 pink, 2 brown
        public double BackgroundContrast { get; set; } = 1.0;  // log-density gain; 0 = uniform background
        public int NoiseGridSize { get; set; } = 32;           // power of two

        // Density gradient within each tube, biased toward its inner edge.
        public double DensityGradientStrength { get; set; } = 0.4; // 0 = uniform, 1 = fully inner-biased

        public double GlobalScale { get; set; } = 1.0;

        // Ambient dimension. 3 = canonical eye; 4 = rings in independent 2-planes (unlinked).
        public int Dimension { get; set; } = 3;

        // Skeleton realization knobs (hyperbolic-friendly; see EyeSkeleton).
        public double[] Center { get; set; } = { 0.0, 0.0, 0.0 };       // origin-centering (padded/truncated to Dimension)
        public double MaxGeodesicRadius { get; set; } = double.PositiveInfinity; // ρ_max shell cap
        public double WarpStrength { get; set; } = 1.0;                 // faithful↔cosmetic dial

        public int Seed { get; set; } = 42;
    }

    /// <summary>
    /// Build the intrinsic, coordinate-free <see cref="EyeSkeleton"/> from a
    /// config — the shared structure both the Euclidean and hyperbolic realizers
    /// consume. Ambient dimension comes from <see cref="EyeTorusToyConfig.Dimension"/>;
    /// stroke frames default to the canonical XY/Z eye in 3-D and fan into
    /// independent 2-planes in 4-D+.
    /// </summary>
    public static EyeSkeleton BuildSkeleton(EyeTorusToyConfig cfg)
    {
        ArgumentNullException.ThrowIfNull(cfg);
        int dim = cfg.Dimension;
        if (dim < 3) throw new ArgumentException("EyeTorusToy requires Dimension >= 3.", nameof(cfg));

        // Stroke frames: (sweepA, sweepB, axial). 3-D = coplanar canonical eye;
        // 4-D+ = a fan of 2-planes sharing axis 0 so the rings are unlinked.
        (int a, int b, int ax) centralF, upperF, lowerF;
        if (dim == 3)
        {
            centralF = upperF = lowerF = (0, 1, 2);
        }
        else
        {
            centralF = (0, 1, 2);
            upperF = (0, 2, 3);
            lowerF = (0, 3, 1);
        }

        EyePupil? pupil = null;
        if (cfg.PupilPoints > 0)
        {
            double bore = Math.Max(0.0, cfg.CentralMajorR - cfg.CentralMinorR);
            double radius = Math.Clamp(cfg.PupilDilation, 0.0, 0.98) * bore;
            pupil = new EyePupil { Radius = radius, PointCount = cfg.PupilPoints, Label = 3 };
        }

        return new EyeSkeleton
        {
            Center = Resize(cfg.Center, dim),
            MaxGeodesicRadius = cfg.MaxGeodesicRadius,
            WarpStrength = cfg.WarpStrength,
            Pupil = pupil,
            Strokes = new[]
            {
                new EyeStroke
                {
                    MajorRadius = cfg.CentralMajorR, MinorRadius = cfg.CentralMinorR,
                    ArcStart = 0.0, ArcEnd = 2.0 * Math.PI, ZThickness = cfg.ZThickness,
                    DensityGradientStrength = cfg.DensityGradientStrength,
                    CrossSection = cfg.CentralCrossSection ?? cfg.CrossSection,
                    RibbonThickness = cfg.RibbonThickness, EndTaper = 0.0,
                    SweepAxisA = centralF.a, SweepAxisB = centralF.b, AxialAxis = centralF.ax,
                    Offset = new double[dim],
                    PointCount = cfg.CentralPoints, Label = 0,
                },
                new EyeStroke
                {
                    MajorRadius = cfg.UpperMajorR, MinorRadius = cfg.UpperMinorR,
                    ArcStart = cfg.UpperArcStart, ArcEnd = cfg.UpperArcEnd, ZThickness = cfg.ZThickness,
                    DensityGradientStrength = cfg.DensityGradientStrength,
                    CrossSection = cfg.CrossSection, RibbonThickness = cfg.RibbonThickness,
                    EndTaper = cfg.HalfArcTaper,
                    SweepAxisA = upperF.a, SweepAxisB = upperF.b, AxialAxis = upperF.ax,
                    Offset = new double[dim],
                    PointCount = cfg.UpperPoints, Label = 1,
                },
                new EyeStroke
                {
                    MajorRadius = cfg.LowerMajorR, MinorRadius = cfg.LowerMinorR,
                    ArcStart = cfg.LowerArcStart, ArcEnd = cfg.LowerArcEnd, ZThickness = cfg.ZThickness,
                    Offset = Resize(cfg.LowerOffset, dim),
                    DensityGradientStrength = cfg.DensityGradientStrength,
                    CrossSection = cfg.CrossSection, RibbonThickness = cfg.RibbonThickness,
                    EndTaper = cfg.HalfArcTaper,
                    SweepAxisA = lowerF.a, SweepAxisB = lowerF.b, AxialAxis = lowerF.ax,
                    PointCount = cfg.LowerPoints, Label = 2,
                },
            },
        };
    }

    /// <summary>
    /// Sample the skeleton into eye-local tangent vectors of dimension
    /// <c>skeleton.Center.Length</c> (relative to the center, before any manifold
    /// placement) plus labels — the shared field a realizer pushes through
    /// <c>exp_center</c>. Each stroke is realized in its own ambient frame.
    /// </summary>
    public static (double[][] Local, int[] Fine, int[] Coarse) SampleLocal(
        EyeSkeleton skeleton, double structureNoiseSigma, Random rng)
    {
        ArgumentNullException.ThrowIfNull(skeleton);
        ArgumentNullException.ThrowIfNull(rng);
        int dim = skeleton.Center.Length;

        int strokePoints = 0;
        foreach (var s in skeleton.Strokes) strokePoints += s.PointCount;
        int pupilPoints = skeleton.Pupil?.PointCount ?? 0;
        int total = strokePoints + pupilPoints;

        var local = new double[total][];
        var fine = new int[total];
        var coarse = new int[total];
        int idx = 0;

        // Cross-section offset (radial a, axial b) in the tube's normal plane.
        (double a, double b) CrossSection(EyeStroke s, double minorEff)
        {
            if (s.CrossSection == CrossSectionShape.Ribbon)
            {
                double aR = minorEff * (2.0 * rng.NextDouble() - 1.0);
                double bR = s.RibbonThickness * SyntheticData.SampleStandardNormal(rng);
                return (aR, bR);
            }

            double g = Math.Clamp(s.DensityGradientStrength, 0.0, 1.0);
            double phi;
            while (true)
            {
                phi = 2.0 * Math.PI * rng.NextDouble();
                double weight = (1.0 - g) + g * (1.0 - Math.Cos(phi)) * 0.5;
                if (rng.NextDouble() <= weight) break;
            }

            double rad = minorEff;
            if (s.CrossSection == CrossSectionShape.Solid)
                rad *= Math.Sqrt(rng.NextDouble());

            return (rad * Math.Cos(phi), rad * Math.Sin(phi) * s.ZThickness);
        }

        void RealizeStroke(EyeStroke s)
        {
            double span = s.ArcEnd - s.ArcStart;
            for (int p = 0; p < s.PointCount; p++)
            {
                double theta = s.ArcStart + rng.NextDouble() * span;
                double t = span != 0.0 ? (theta - s.ArcStart) / span : 0.0;
                double taper = (1.0 - s.EndTaper) + s.EndTaper * Math.Sin(Math.PI * t);
                double minorEff = s.MinorRadius * taper;

                var (a, b) = CrossSection(s, minorEff);
                double ringR = s.MajorRadius + a;

                var v = new double[dim];
                v[s.SweepAxisA] += ringR * Math.Cos(theta);
                v[s.SweepAxisB] += ringR * Math.Sin(theta);
                v[s.AxialAxis] += b;
                for (int d = 0; d < dim; d++) v[d] += s.Offset[d];
                for (int d = 0; d < dim; d++) v[d] += structureNoiseSigma * SyntheticData.SampleStandardNormal(rng);

                local[idx] = v;
                fine[idx] = s.Label;
                coarse[idx] = 0;
                idx++;
            }
        }

        foreach (var s in skeleton.Strokes) RealizeStroke(s);

        if (skeleton.Pupil is EyePupil pup)
        {
            for (int i = 0; i < pup.PointCount; i++)
            {
                var dir = new double[dim];
                double n2 = 0.0;
                for (int d = 0; d < dim; d++) { dir[d] = SyntheticData.SampleStandardNormal(rng); n2 += dir[d] * dir[d]; }
                double norm = Math.Sqrt(n2);
                if (norm < 1e-12) { dir[0] = 1.0; norm = 1.0; }
                double rr = pup.Radius * Math.Pow(rng.NextDouble(), 1.0 / dim);

                var v = new double[dim];
                for (int d = 0; d < dim; d++) v[d] = rr * dir[d] / norm;
                local[idx] = v;
                fine[idx] = pup.Label;
                coarse[idx] = 0;
                idx++;
            }
        }

        return (local, fine, coarse);
    }

    public static SyntheticDataset Generate(EyeTorusToyConfig? cfg = null, int? overrideSeed = null)
    {
        cfg ??= new EyeTorusToyConfig();
        int seed = overrideSeed ?? cfg.Seed;
        var rng = new Random(seed);

        var skeleton = BuildSkeleton(cfg);
        int dim = skeleton.Center.Length;
        double[] center = skeleton.Center;
        double rhoMax = skeleton.MaxGeodesicRadius;
        double scale = cfg.GlobalScale;

        var (local, fine, coarse) = SampleLocal(skeleton, cfg.StructureNoiseSigma, rng);

        int totalStructure = local.Length;
        int structureGroups = skeleton.Strokes.Count + (skeleton.Pupil is null ? 0 : 1);
        int bgLabel = structureGroups;
        int backgroundPoints = (int)(totalStructure * cfg.BackgroundDensityRatio);
        int n = totalStructure + backgroundPoints;

        var features = new double[n][];
        var fineLabels = new int[n];
        var coarseLabels = new int[n];

        // Flat realization: exp_center(v) = center + v, with the ρ_max cap as a
        // Euclidean-radius clamp.
        for (int i = 0; i < totalStructure; i++)
        {
            var pt = new double[dim];
            for (int d = 0; d < dim; d++) pt[d] = local[i][d] * scale + center[d];
            CapToShell(pt, center, rhoMax);
            features[i] = pt;
            fineLabels[i] = fine[i];
            coarseLabels[i] = coarse[i];
        }

        int idx = totalStructure;

        // Background blanket: an N-D colored-noise density field (power-law
        // spectral field, rejection-sampled into a log-normal density — pink/brown
        // grows phantom blobs). The noise stream is its own deterministic die so
        // structure sampling stays bit-stable. The box is wide in the first two
        // axes, flatter in the rest.
        double bound = (cfg.CentralMajorR + cfg.LowerMajorR + cfg.LowerMinorR) * scale * 1.4;
        double zBound = bound * 0.6;
        if (backgroundPoints > 0)
        {
            var boxMin = new double[dim];
            var boxMax = new double[dim];
            for (int d = 0; d < dim; d++) { double bb = d < 2 ? bound : zBound; boxMin[d] = -bb; boxMax[d] = bb; }

            var noiseRng = new Xoshiro256PlusPlus(seed ^ unchecked((int)0x9E3779B9));
            var field = SpectralNoiseField.Generate(
                noiseRng, cfg.NoiseGridSize, cfg.SpectralExponent, boxMin, boxMax);

            double contrast = cfg.BackgroundContrast;
            double logMax = contrast * field.Max;
            var cand = new double[dim];

            int placed = 0;
            long guard = 0, guardMax = (long)backgroundPoints * 1000 + 1000;
            while (placed < backgroundPoints && guard++ < guardMax)
            {
                for (int d = 0; d < dim; d++) cand[d] = (rng.NextDouble() * 2 - 1) * (d < 2 ? bound : zBound);
                double accept = contrast == 0.0 ? 1.0 : Math.Exp(contrast * field.Sample(cand) - logMax);
                if (rng.NextDouble() > accept) continue;
                var bgp = new double[dim];
                Array.Copy(cand, bgp, dim);
                features[idx] = bgp;
                fineLabels[idx] = bgLabel;
                coarseLabels[idx] = 1;
                idx++;
                placed++;
            }
            while (placed < backgroundPoints)
            {
                var bgp = new double[dim];
                for (int d = 0; d < dim; d++) bgp[d] = (rng.NextDouble() * 2 - 1) * (d < 2 ? bound : zBound);
                features[idx] = bgp;
                fineLabels[idx] = bgLabel;
                coarseLabels[idx] = 1;
                idx++;
                placed++;
            }
        }

        return new SyntheticDataset
        {
            Features = features,
            Labels = fineLabels,
            ClusterCount = structureGroups,
            LabelsByLevel = new[] { coarseLabels, fineLabels },
            Parameters = new Dictionary<string, object>
            {
                ["generator"] = nameof(EyeTorusToy),
                ["dimension"] = dim,
                ["totalPoints"] = n,
                ["structurePoints"] = totalStructure,
                ["backgroundPoints"] = backgroundPoints,
                ["backgroundLabel"] = bgLabel,
                ["structureGroups"] = structureGroups,
                ["pupilPoints"] = skeleton.Pupil?.PointCount ?? 0,
                ["seed"] = seed,
                ["crossSection"] = cfg.CrossSection.ToString(),
                ["centralCrossSection"] = (cfg.CentralCrossSection ?? cfg.CrossSection).ToString(),
                ["halfArcTaper"] = cfg.HalfArcTaper,
                ["backgroundKind"] = "colored-spectral",
                ["spectralExponent"] = cfg.SpectralExponent,
                ["backgroundContrast"] = cfg.BackgroundContrast,
                ["noiseGridSize"] = cfg.NoiseGridSize,
                ["backgroundDensityRatio"] = cfg.BackgroundDensityRatio,
                ["densityGradientStrength"] = cfg.DensityGradientStrength,
                ["zThickness"] = cfg.ZThickness,
                ["globalScale"] = cfg.GlobalScale,
                ["maxGeodesicRadius"] = double.IsPositiveInfinity(rhoMax) ? "infinity" : rhoMax,
                ["warpStrength"] = cfg.WarpStrength,
                ["reference"] = "Blatt, Wiseman, Domany 1996 PRL 76:3251, Fig. 1 (toroidal lift)"
            },
            Metadata = new SyntheticDatasetMeta(
                GeneratorName: nameof(EyeTorusToy),
                GeometryClass: "Euclidean",
                TopologyTag: "toroidal-hierarchical",
                HierarchyTag: "multi-scale-density",
                GTNumClusters: structureGroups,
                AmbientDimensionality: dim,
                LiteratureReference: "Blatt, Wiseman, Domany 1996 PRL 76:3251, Fig. 1 — canonical SPC \"Toy\" (toroidal lift)")
        };
    }

    /// <summary>
    /// Orthographic projection onto the first two axes (drop the rest) — recovers
    /// the 2-D approximation of the original BWD 1996 Fig. 1 plate, and the naive
    /// geometry-ignoring flatten of a curved or higher-dim eye for the distortion
    /// demo. Labels and level structure carry over; ambient dimensionality → 2.
    /// </summary>
    public static SyntheticDataset FlattenToPlane(SyntheticDataset toy)
    {
        ArgumentNullException.ThrowIfNull(toy);
        var src = toy.Features;
        var flat = new double[src.Length][];
        for (int i = 0; i < src.Length; i++)
            flat[i] = new[] { src[i][0], src[i][1] };

        return new SyntheticDataset
        {
            Features = flat,
            Labels = toy.Labels,
            ClusterCount = toy.ClusterCount,
            LabelsByLevel = toy.LabelsByLevel,
            Parameters = new Dictionary<string, object>(toy.Parameters) { ["projection"] = "orthographic-xy" },
            Metadata = toy.Metadata is null ? null : toy.Metadata with { AmbientDimensionality = 2 }
        };
    }

    private static double[] Resize(double[] src, int dim)
    {
        var dst = new double[dim];
        int copy = Math.Min(src.Length, dim);
        Array.Copy(src, dst, copy);
        return dst;
    }

    private static void CapToShell(double[] pt, double[] center, double rhoMax)
    {
        if (double.IsPositiveInfinity(rhoMax)) return;
        double dist2 = 0.0;
        for (int d = 0; d < pt.Length; d++) { double dd = pt[d] - center[d]; dist2 += dd * dd; }
        double dist = Math.Sqrt(dist2);
        if (dist > rhoMax && dist > 0.0)
        {
            double k = rhoMax / dist;
            for (int d = 0; d < pt.Length; d++) pt[d] = center[d] + (pt[d] - center[d]) * k;
        }
    }
}
