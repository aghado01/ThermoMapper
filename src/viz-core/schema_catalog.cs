namespace Viz
{
    /// <summary>
    /// Registry of generator parameter schemas.
    /// Each generator known to VizApi is listed here so the viewer can build
    /// its regen panel dynamically — no hardcoded HTML per generator.
    ///
    /// The Key strings must match the "generator" value in generator_params.
    /// Add an entry here whenever a new synthetic generator is wired into VizApi.
    /// </summary>
    public static class SchemaCatalog
    {
        /// <summary>
        /// Returns the schema for the named generator, or null if none is registered.
        /// </summary>
        /// <summary>
        /// All generator names known to VizApi, in display order.
        /// The viewer reads this to populate its generator picker dropdown.
        /// Must stay in sync with the switch in BuildPackage (Program.cs).
        /// </summary>
        public static readonly string[] KnownGenerators =
        [
            "CrescentAndEllipsoid",
            "MobiusAndEllipsoid",
            "HyperbolicBlobs",
            "HyperbolicBlattHierarchy",
            "Simplex",
            "GaussianManifold",
            "TwoMoons",
            "BlattHierarchy",
            "BlattThreeCluster",
        ];

        /// <summary>
        /// Returns the schema for the named generator, or null if none is registered.
        /// </summary>
        public static GeneratorParamSchema? ForGenerator(string? generatorName) =>
            generatorName switch
            {
                "CrescentAndEllipsoid" => CrescentAndEllipsoidSchema,
                "MobiusAndEllipsoid" => MobiusAndEllipsoidSchema,
                "HyperbolicBlobs" => HyperbolicBlobsSchema,
                "HyperbolicBlattHierarchy" => HyperbolicBlattHierarchySchema,
                "Simplex" => SimplexSchema,
                "GaussianManifold" => GaussianManifoldSchema,
                "TwoMoons" => TwoMoonsSchema,
                "BlattHierarchy" => BlattHierarchySchema,
                "BlattThreeCluster" => BlattThreeClusterSchema,
                _ => null,
            };

        // ── Shared sections ───────────────────────────────────────────────────

        private static readonly ParamSection SharedSection = new()
        {
            Label = "Shared",
            Params =
            [
                new() { Key = "seed",  Label = "seed",  Type = "int", Min = 0, Max = 99999 },
                new() { Key = "knnK",  Label = "KNN k", Type = "int", Min = 1, Max = 30    },
            ],
        };

        private static readonly ParamSection GraphSection = new()
        {
            Label = "Graph",
            Params =
            [
                new() { Key = "metric",   Label = "metric", Type = "enum",
                    EnumValues = ["euclidean", "manhattan", "minkowski:p=2", "cosine", "poincare", "hamming"] },
                new() { Key = "neighborRule", Label = "rule", Type = "enum",
                        EnumValues = ["Knn", "MutualKnn", "EpsilonBall", "MstAugmented"] },
                new() { Key = "epsilonBallEpsilon", Label = "epsilon", Type = "float",
                        Min = 0.05, Max = 20.0, Step = 0.05 },
                new() { Key = "kernel", Label = "kernel", Type = "enum",
                        EnumValues = ["Gaussian", "Cauchy", "Laplacian", "Linear"] },
                new() { Key = "bandwidth", Label = "bandwidth", Type = "float",
                        Min = 0.0, Max = 10.0, Step = 0.05,
                        Description = "0 = auto-estimate" },
            ],
        };

        // Overlay section: cross-generator controls for post-generation analysis layers.
        // Appears in every generator schema so the same slider works regardless of which
        // dataset is loaded.
        private static readonly ParamSection OverlaySection = new()
        {
            Label = "Overlay",
            Params =
            [
                new() { Key = "gmmComponents", Label = "GMM K", Type = "float",
                        Min = 0, Max = 4, Step = 1.0,
                        DisplayValues = ["1", "2", "4", "8", "16"] },
                new() { Key = "showFlow", Label = "flow", Type = "bool" },
            ],
        };

        // ── CrescentAndEllipsoid ──────────────────────────────────────────────

        private static readonly GeneratorParamSchema CrescentAndEllipsoidSchema = new()
        {
            Sections =
            [
                new()
                {
                    Label  = "Crescent",
                    Params =
                    [
                        new() { Key = "crescentPoints", Label = "N",      Type = "int",   Min = 30,   Max = 1500 },
                        new() { Key = "crescentRadius", Label = "radius", Type = "float", Min = 1.0,  Max = 6.0,  Step = 0.1  },
                        new() { Key = "crescentWidth",   Label = "width",   Type = "float", Min = 0.05, Max = 1.5, Step = 0.05 },
                        new() { Key = "arcHalfAngle",    Label = "arc ½φ",  Type = "float", Min = 0.3,  Max = 3.1, Step = 0.05 },
                    ],
                },
                new()
                {
                    Label  = "Ellipsoid",
                    Params =
                    [
                        new() { Key = "ellipsoidPoints",    Label = "N",          Type = "int",  Min = 30,  Max = 1500 },
                        new() { Key = "ellipsoidAxes",      Label = "axes",       Type = "vec3", Min = 0.1, Max = 5.0, Step = 0.05,
                                VecLabels = ["a", "b", "c"] },
                        new() { Key = "ellipsoidShellMode", Label = "shell mode", Type = "enum",
                                EnumValues = ["Solid", "Gaussian", "Hollow", "Annular"] },
                        new() { Key = "placement",          Label = "placement",  Type = "enum",
                                EnumValues = ["NearOpenFace", "OrthogonalElbowIntersect", "IntersectUpperTip", "IntersectLowerTip"] },
                        new() { Key = "intersectDepth",        Label = "depth",       Type = "float", Min = -3.0, Max = 3.0, Step = 0.1 },
                        new() { Key = "intersectRadialShift",  Label = "radial",      Type = "float", Min = -3.0, Max = 3.0, Step = 0.1 },
                        new() { Key = "gapScale",              Label = "gap×",        Type = "float", Min =  0.2, Max = 3.0, Step = 0.1 },
                    ],
                },
            ],
        };

        // ── HyperbolicBlobs ───────────────────────────────────────────────────
        // Paired metric: Poincaré. Points live strictly inside the open unit ball
        // so the hyperbolic geodesic is exact, not clipped.

        private static readonly GeneratorParamSchema HyperbolicBlobsSchema = new()
        {
            Sections =
            [
                new()
                {
                    Label  = "Hyperbolic Blobs",
                    Params =
                    [
                        new() { Key = "clusterCount",     Label = "K (clusters)",  Type = "int",   Min = 2,  Max = 8                    },
                        new() { Key = "pointsPerCluster", Label = "N per cluster", Type = "int",   Min = 20, Max = 500                  },
                        new() { Key = "dimensions",       Label = "dim",           Type = "enum",
                                EnumValues = ["2", "3"] },
                        new() { Key = "separation",       Label = "separation",    Type = "float", Min = 0.5, Max = 5.0, Step = 0.1     },
                        new() { Key = "spread",           Label = "spread",        Type = "float", Min = 0.05, Max = 1.5, Step = 0.05   },
                    ],
                },
            ],
        };

        // ── HyperbolicBlattHierarchy ─────────────────────────────────────────

        private static readonly GeneratorParamSchema HyperbolicBlattHierarchySchema = new()
        {
            Sections =
            [
                new()
                {
                    Label = "Hierarchy",
                    Params =
                    [
                        new() { Key = "hierarchyPoints", Label = "N total", Type = "int", Min = 50, Max = 5000 },
                        new() { Key = "hierarchyDepth", Label = "depth", Type = "int", Min = 1, Max = 6 },
                        new() { Key = "branchesPerNode", Label = "branches", Type = "int", Min = 2, Max = 6 },
                        new() { Key = "basePointsPerLeaf", Label = "leaf floor", Type = "int", Min = 5, Max = 200 },
                        new() { Key = "dimensions", Label = "dim", Type = "enum", EnumValues = ["2", "3"] },
                    ],
                },
                new()
                {
                    Label = "Scale",
                    Params =
                    [
                        new() { Key = "separation", Label = "separation", Type = "float", Min = 0.5, Max = 6.0, Step = 0.1 },
                        new() { Key = "radiusDecay", Label = "radius decay", Type = "float", Min = 0.2, Max = 0.95, Step = 0.05 },
                        new() { Key = "spread", Label = "noise", Type = "float", Min = 0.01, Max = 1.5, Step = 0.01 },
                    ],
                },
            ],
        };

        // ── Simplex ───────────────────────────────────────────────────────────
        // Paired metric: Fisher–Rao (simplex). Each point is a probability vector;
        // first three categories project as XYZ. With disjointSupports = true the
        // clusters concentrate on different cardinal axes — visually distinctive.

        private static readonly GeneratorParamSchema SimplexSchema = new()
        {
            Sections =
            [
                new()
                {
                    Label  = "Simplex",
                    Params =
                    [
                        new() { Key = "clusterCount",     Label = "K (clusters)",  Type = "int",   Min = 2,   Max = 8                  },
                        new() { Key = "pointsPerCluster", Label = "N per cluster", Type = "int",   Min = 20,  Max = 500                },
                        new() { Key = "categories",       Label = "categories",    Type = "int",   Min = 3,   Max = 30                 },
                        new() { Key = "disjointSupports", Label = "disjoint",      Type = "bool"  },
                        new() { Key = "concentration",    Label = "concentration", Type = "float", Min = 2.0, Max = 200.0, Step = 1.0  },
                    ],
                },
            ],
        };

        // ── GaussianManifold ──────────────────────────────────────────────────
        // Paired metric: Fisher–Rao (half-plane). Each point is a (μ, log σ)
        // pair describing a 1D Gaussian. Centers are placed so Euclidean and
        // Fisher–Rao topologies disagree — the storytelling case for geodesics.

        private static readonly GeneratorParamSchema GaussianManifoldSchema = new()
        {
            Sections =
            [
                new()
                {
                    Label  = "Gaussian Manifold",
                    Params =
                    [
                        new() { Key = "clusterCount",     Label = "K (clusters)",  Type = "int",   Min = 2,   Max = 8                 },
                        new() { Key = "pointsPerCluster", Label = "N per cluster", Type = "int",   Min = 20,  Max = 500               },
                        new() { Key = "clusterRadius",    Label = "spread",        Type = "float", Min = 0.05, Max = 1.0, Step = 0.05 },
                    ],
                },
            ],
        };

        // ── TwoMoons ──────────────────────────────────────────────────────────

        private static readonly GeneratorParamSchema TwoMoonsSchema = new()
        {
            Sections =
            [
                new()
                {
                    Label  = "Two Moons",
                    Params =
                    [
                        new() { Key = "pointsPerMoon", Label = "N per moon", Type = "int", Min = 20, Max = 2000 },
                        new() { Key = "noise", Label = "noise", Type = "float", Min = 0.0, Max = 1.0, Step = 0.01 },
                    ],
                },
            ],
        };

        // ── BlattHierarchy ───────────────────────────────────────────────────

        private static readonly GeneratorParamSchema BlattHierarchySchema = new()
        {
            Sections =
            [
                new()
                {
                    Label  = "Hierarchy",
                    Params =
                    [
                        new() { Key = "coarseClusters", Label = "coarse K", Type = "int", Min = 1, Max = 8 },
                        new() { Key = "mediumPerCoarse", Label = "medium / coarse", Type = "int", Min = 1, Max = 8 },
                        new() { Key = "finePerMedium", Label = "fine / medium", Type = "int", Min = 1, Max = 8 },
                        new() { Key = "pointsPerFine", Label = "N per fine", Type = "int", Min = 5, Max = 500 },
                    ],
                },
                new()
                {
                    Label  = "Scale",
                    Params =
                    [
                        new() { Key = "coarseSeparation", Label = "coarse sep", Type = "float", Min = 1.0, Max = 50.0, Step = 0.1 },
                        new() { Key = "mediumSeparation", Label = "medium sep", Type = "float", Min = 0.1, Max = 20.0, Step = 0.1 },
                        new() { Key = "fineSeparation", Label = "fine sep", Type = "float", Min = 0.05, Max = 5.0, Step = 0.05 },
                        new() { Key = "leafSpread", Label = "leaf spread", Type = "float", Min = 0.01, Max = 2.0, Step = 0.01 },
                    ],
                },
            ],
        };

        // ── BlattThreeCluster ────────────────────────────────────────────────

        private static readonly GeneratorParamSchema BlattThreeClusterSchema = new()
        {
            Sections =
            [
                new()
                {
                    Label  = "Blatt Three Cluster",
                    Params =
                    [
                        new() { Key = "pointsPerCluster", Label = "N per cluster", Type = "int", Min = 10, Max = 2000 },
                        new() { Key = "stdDev", Label = "std dev", Type = "float", Min = 0.05, Max = 5.0, Step = 0.05 },
                    ],
                },
            ],
        };

        // ── MobiusAndEllipsoid ────────────────────────────────────────────────

        private static readonly GeneratorParamSchema MobiusAndEllipsoidSchema = new()
        {
            Sections =
            [
                new()
                {
                    Label  = "Möbius",
                    Params =
                    [
                        new() { Key = "crossSection",   Label = "cross-section", Type = "enum",
                                EnumValues = ["Ribbon", "GaussianIsotropic", "GaussianAnisotropic", "UniformDisk", "Annular"] },
                        new() { Key = "dimensions",     Label = "embedding",     Type = "enum",
                                EnumValues = ["3", "4"] },
                        new() { Key = "mobiusPoints",   Label = "N",             Type = "int",   Min = 50,   Max = 2000             },
                        new() { Key = "spineRadius",    Label = "spine radius",  Type = "float", Min = 1.0,  Max = 8.0,  Step = 0.1  },
                        new() { Key = "halfWidth",      Label = "width",         Type = "float", Min = 0.05, Max = 2.5,  Step = 0.05 },
                        new() { Key = "halfThickness",  Label = "thickness",     Type = "float", Min = 0.01, Max = 1.5,  Step = 0.01 },
                        new() { Key = "noiseSigma",     Label = "noise σ",       Type = "float", Min = 0.0,  Max = 0.5,  Step = 0.01 },
                        new() { Key = "twistCount",     Label = "twists",        Type = "int",   Min = 1,    Max = 5                 },
                        new() { Key = "radialBias",     Label = "radial bias",   Type = "float", Min = 0.0,  Max = 1.0,  Step = 0.05 },
                        new() { Key = "spineShape",     Label = "spine shape",   Type = "enum",
                                EnumValues = ["Circle", "FigureEight"] },
                        new() { Key = "splayFactor",    Label = "splay",         Type = "float", Min = 0.0,  Max = 1.0,  Step = 0.05 },
                    ],
                },
                new()
                {
                    Label  = "Ellipsoid",
                    Params =
                    [
                        new() { Key = "ellipsoidPoints",    Label = "N",          Type = "int",  Min = 30,  Max = 1500 },
                        new() { Key = "ellipsoidAxes",      Label = "axes",       Type = "vec3", Min = 0.1, Max = 3.0, Step = 0.05,
                                VecLabels = ["a", "b", "c"] },
                        new() { Key = "ellipsoidShellMode", Label = "shell mode", Type = "enum",
                                EnumValues = ["Solid", "Gaussian", "Hollow", "Annular"] },
                    ],
                },
                new()
                {
                    Label  = "Placement",
                    Params =
                    [
                        new() { Key = "placement",            Label = "placement",   Type = "enum",
                                EnumValues = ["NearSeam", "CenterCrossOrtho", "PeripheralElbow", "CenterCrossCoPlanar", "Manual"] },
                        new() { Key = "intersectDepth",       Label = "depth",       Type = "float", Min = -3.0, Max = 3.0, Step = 0.1 },
                        new() { Key = "intersectRadialShift", Label = "radial shift",Type = "float", Min = -3.0, Max = 3.0, Step = 0.1 },
                        new() { Key = "gapScale",             Label = "gap×",        Type = "float", Min =  0.2, Max = 3.0, Step = 0.1 },
                    ],
                },
            ],
        };
    }
}
