using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.Json;
using System.Text.Json.Serialization;
using Maths.LinAlg;
using Viz;

namespace Viz.Renderers
{
    /// <summary>
    /// Renders a ScenePackage as a self-contained JSON document.
    ///
    /// Schema contract (stable across versions):
    ///   • "schema_version": integer — bump when breaking changes are made.
    ///   • All layer arrays mirror the VizCore type names exactly (snake_case).
    ///   • Numeric arrays are flat row-major with accompanying shape fields.
    ///   • Null fields are omitted (JsonIgnoreCondition.WhenWritingNull).
    ///   • Enum values are written as strings, not integers.
    /// </summary>
    public sealed class JsonExportRenderTarget : IRenderTarget
    {
        public const int SchemaVersion = 1;

        private static readonly JsonSerializerOptions _indentedOptions = new()
        {
            WriteIndented = true,
            DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
            Converters = { new JsonStringEnumConverter() },
        };

        private static readonly JsonSerializerOptions _compactOptions = new()
        {
            WriteIndented = false,
            DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
            Converters = { new JsonStringEnumConverter() },
        };

        private readonly JsonSerializerOptions _options;

        /// <param name="compact">
        /// When true, emits minified JSON (no indentation or extra whitespace).
        /// Recommended for large scenes (N ≥ 50k, multiple edge layers) where
        /// indentation inflates file size noticeably. Default false = readable.
        /// </param>
        public JsonExportRenderTarget(bool compact = false)
        {
            _options = compact ? _compactOptions : _indentedOptions;
        }

        public void Render(ScenePackage scene, Stream output)
        {
            if (scene is null) throw new ArgumentNullException(nameof(scene));
            if (output is null) throw new ArgumentNullException(nameof(output));

            var doc = BuildDocument(scene);
            string json = JsonSerializer.Serialize(doc, _options);
            byte[] bytes = Encoding.UTF8.GetBytes(json);
            output.Write(bytes, 0, bytes.Length);
        }

        // ---------------------------------------------------------------
        // Top-level document assembly
        // ---------------------------------------------------------------

        private static ScenePackageJson BuildDocument(ScenePackage scene)
        {
            return new ScenePackageJson
            {
                SchemaVersion = SchemaVersion,
                // Prefer the descriptor-supplied title (e.g. "Hyperbolic Blobs");
                // fall back to the active label layer name, then the cloud label.
                Title = !string.IsNullOrEmpty(scene.Title)
                                         ? scene.Title
                                         : scene.ActiveLabelLayers.Count > 0
                                             ? scene.ActiveLabelLayers[0].Name
                                             : scene.Points.Label ?? "Untitled",
                Points = SerializePointCloud(scene.Points),
                LabelLayers = SerializeLabelLayers(scene.ActiveLabelLayers),
                NodeSignalLayers = SerializeNodeSignalLayers(scene.ActiveNodeSignalLayers),
                EdgeLayers = SerializeEdgeLayers(scene.ActiveEdgeLayers),
                TriangleLayers = SerializeTriangleLayers(scene.ActiveTriangleLayers),
                GaussianLayers = SerializeGaussianLayers(scene.ActiveGaussianLayers),
                SpineLayers = SerializeSpineLayers(scene.ActiveSpineLayers),
                VectorFieldLayers = SerializeVectorFieldLayers(scene.ActiveVectorFieldLayers),
                LineFieldLayers = SerializeLineFieldLayers(scene.ActiveLineFieldLayers),
                TemporalSequence = scene.ActiveSequence is null
                                         ? null
                                         : SerializeTemporalSequence(scene.ActiveSequence),
                TemporalFrameIndex = scene.ActiveSequence is null ? null : scene.TemporalFrameIndex,
                Hints = SerializeHints(scene.Hints),
                GeneratorParams = scene.GeneratorParams is not null
                    ? new Dictionary<string, object>(scene.GeneratorParams)
                    : null,
                GeneratorParamSchema = scene.GeneratorParamSchema,
                GeneratorCatalog = SchemaCatalog.KnownGenerators,
            };
        }

        // ---------------------------------------------------------------
        // PointCloud
        // ---------------------------------------------------------------

        private static PointCloudJson SerializePointCloud(PointCloud cloud)
        {
            return new PointCloudJson
            {
                N = cloud.N,
                D = cloud.D,
                Label = cloud.Label,
                Features = cloud.Features.ToArray(), // flat [N×D], stride: i*D+d
            };
        }

        // ---------------------------------------------------------------
        // LabelLayer
        // ---------------------------------------------------------------

        private static List<LabelLayerJson> SerializeLabelLayers(
            IReadOnlyList<LabelLayer> layers)
        {
            var result = new List<LabelLayerJson>(layers.Count);
            foreach (var layer in layers)
                result.Add(new LabelLayerJson
                {
                    Id = layer.Id,
                    Name = layer.Name,
                    Kind = layer.Kind,
                    Labels = layer.Labels.ToArray(),
                });
            return result;
        }

        // ---------------------------------------------------------------
        // NodeSignalLayer  (subsumes the older ScalarLayer)
        // ---------------------------------------------------------------

        private static List<NodeSignalLayerJson> SerializeNodeSignalLayers(
            IReadOnlyList<NodeSignalLayer> layers)
        {
            var result = new List<NodeSignalLayerJson>(layers.Count);
            foreach (var layer in layers)
            {
                var span = layer.Values.Span;
                double min = double.MaxValue, max = double.MinValue;
                foreach (var v in span) { if (v < min) min = v; if (v > max) max = v; }
                result.Add(new NodeSignalLayerJson
                {
                    Id = layer.Id,
                    Name = layer.Name,
                    Source = layer.Source,
                    Values = layer.Values.ToArray(),
                    Min = min,
                    Max = max,
                });
            }
            return result;
        }

        // ---------------------------------------------------------------
        // EdgeLayer
        // ---------------------------------------------------------------

        private static List<EdgeLayerJson> SerializeEdgeLayers(
            IReadOnlyList<EdgeLayer> layers)
        {
            var result = new List<EdgeLayerJson>(layers.Count);
            foreach (var layer in layers)
                result.Add(new EdgeLayerJson
                {
                    Id = layer.Id,
                    Name = layer.Name,
                    Metric = layer.Metric,
                    Proximity = SerializeProximitySpec(layer.Proximity),
                    Src = layer.Src.ToArray(),
                    Dst = layer.Dst.ToArray(),
                    Weight = layer.Weight.ToArray(),
                    // is_false_bridge[i] = EdgeClusterSrc[i] != EdgeClusterDst[i]
                    EdgeClusterSrc = layer.EdgeClusterSrc?.ToArray(),
                    EdgeClusterDst = layer.EdgeClusterDst?.ToArray(),
                });
            return result;
        }

        // ---------------------------------------------------------------
        // TriangleLayer
        // ---------------------------------------------------------------

        private static List<TriangleLayerJson> SerializeTriangleLayers(
            IReadOnlyList<TriangleLayer> layers)
        {
            var result = new List<TriangleLayerJson>(layers.Count);
            foreach (var layer in layers)
                result.Add(new TriangleLayerJson
                {
                    Id = layer.Id,
                    Name = layer.Name,
                    Source = layer.Source,
                    SourceEdgeLayerId = layer.SourceEdgeLayerId,
                    Vertices = layer.Vertices.ToArray(),
                });
            return result;
        }

        private static ProximitySpecJson? SerializeProximitySpec(ProximitySpec? spec) =>
            spec switch
            {
                null => null,
                KnnSpec s => new ProximitySpecJson { Kind = spec.Kind, K = s.K },
                MutualKnnSpec s => new ProximitySpecJson { Kind = spec.Kind, K = s.K },
                EpsilonBallSpec s => new ProximitySpecJson { Kind = spec.Kind, Epsilon = s.Epsilon },
                MstAugmentedSpec s => new ProximitySpecJson { Kind = spec.Kind, K = s.K },
                _ => new ProximitySpecJson { Kind = spec.Kind },
            };

        // ---------------------------------------------------------------
        // GaussianLayer
        // ---------------------------------------------------------------

        private static List<GaussianLayerJson> SerializeGaussianLayers(
            IReadOnlyList<GaussianLayer> layers)
        {
            var result = new List<GaussianLayerJson>(layers.Count);
            foreach (var layer in layers)
            {
                // Pre-compute Cholesky factor L (Σ = LLᵀ) per component so JS
                // can load it directly into a Matrix4 — no math in the browser.
                var cholL = ComputeCholeskiFactors(layer);
                result.Add(new GaussianLayerJson
                {
                    Id = layer.Id,
                    Name = layer.Name,
                    K = layer.K,
                    D = layer.D,
                    Means = layer.Means.ToArray(),            // flat [K×D]
                    Covariances = layer.Covariances.ToArray(), // flat [K×D×D]
                    CholeskyL = cholL,                         // flat [K×3×3], lower-triangular
                    Weights = layer.Weights.ToArray(),
                    ComponentToClusterMap = layer.ComponentToClusterMap?.ToArray(),
                });
            }
            return result;
        }

        // Compute flat [K×3×3] Cholesky factors from a GaussianLayer.
        // The viewer always renders in 3D, so we emit a 3×3 factor regardless of
        // the source covariance dimension:
        //   D == 3: straight copy + Cholesky.
        //   D <  3: embed the D×D covariance into the top-left of a 3×3 block;
        //           fill the remaining diagonal with a tiny ε so the Cholesky is
        //           numerically stable. The resulting ellipsoid is a thin disc /
        //           line in the missing dimensions — exactly the right look for
        //           a low-D point cloud rendered in a 3D scene.
        //   D >  3: take the leading 3×3 block (viewer can only show 3 dims).
        private const double SerializerCovEps = 1e-8;
        private static double[] ComputeCholeskiFactors(GaussianLayer layer)
        {
            const int outDim = 3;
            int srcDim = layer.D;
            var covs = layer.Covariances.Span;
            var result = new double[layer.K * outDim * outDim];
            var tmp = new double[outDim, outDim];
            var chol = new CholeskyDecomposition(outDim);
            int srcStride = srcDim * srcDim;
            int outStride = outDim * outDim;

            for (int k = 0; k < layer.K; k++)
            {
                int covBase = k * srcStride;
                int copyDim = Math.Min(srcDim, outDim);
                for (int r = 0; r < outDim; r++)
                    for (int c = 0; c < outDim; c++)
                        tmp[r, c] = (r < copyDim && c < copyDim)
                            ? covs[covBase + r * srcDim + c]
                            : (r == c ? SerializerCovEps : 0.0);

                chol.Decompose(tmp);
                chol.WriteLTo(result.AsSpan(k * outStride, outStride));
            }
            return result;
        }

        // ---------------------------------------------------------------
        // SpineLayer  (jagged → flat with explicit shape)
        // ---------------------------------------------------------------

        private static List<SpineLayerJson> SerializeSpineLayers(
            IReadOnlyList<SpineLayer> layers)
        {
            var result = new List<SpineLayerJson>(layers.Count);
            foreach (var layer in layers)
            {
                int m = layer.SpineSamples.Length;
                int d = m > 0 ? layer.SpineSamples[0].Length : 0;

                var flatSamples = new double[m * d];
                for (int i = 0; i < m; i++)
                    layer.SpineSamples[i].CopyTo(flatSamples, i * d);

                double[]? flatTangents = null;
                if (layer.TangentBases is not null)
                {
                    flatTangents = new double[m * d * d];
                    for (int i = 0; i < m; i++)
                        for (int r = 0; r < d; r++)
                            layer.TangentBases[i][r].CopyTo(flatTangents, i * d * d + r * d);
                }

                result.Add(new SpineLayerJson
                {
                    Id = layer.Id,
                    Name = layer.Name,
                    ClusterIdx = layer.ClusterIdx,
                    Kind = layer.Kind,
                    SpineSamplesShape = new[] { m, d },
                    SpineSamples = flatSamples,
                    TangentBasesShape = layer.TangentBases is null ? null : new[] { m, d, d },
                    TangentBases = flatTangents,
                });
            }
            return result;
        }

        // ---------------------------------------------------------------
        // TemporalLabelSequence
        // ---------------------------------------------------------------

        private static TemporalSequenceJson SerializeTemporalSequence(
            TemporalLabelSequence seq)
        {
            var frames = new List<LabelLayerJson>(seq.Frames.Count);
            foreach (var frame in seq.Frames)
                frames.Add(new LabelLayerJson
                {
                    Name = frame.Name,
                    Kind = frame.Kind,
                    Labels = frame.Labels.ToArray(),
                });

            return new TemporalSequenceJson
            {
                Name = seq.Name,
                Axis = seq.Axis,
                IndexValues = seq.IndexValues.ToArray(),
                Frames = frames,
            };
        }

        // ---------------------------------------------------------------
        // SceneRenderHints
        // ---------------------------------------------------------------

        private static SceneRenderHintsJson SerializeHints(SceneRenderHints hints) =>
            new SceneRenderHintsJson
            {
                ShowEdgeWeightsAsOpacity = hints.ShowEdgeWeightsAsOpacity,
                ShowGaussianEllipsoids = hints.ShowGaussianEllipsoids,
                ShowSpineOverlays = hints.ShowSpineOverlays,
                ShowTangentBases = hints.ShowTangentBases,
                OverlayComponentAndClusterColoring = hints.OverlayComponentAndClusterColoring,
                AnnotateSpinColorVsEquilibrium = hints.AnnotateSpinColorVsEquilibrium,
                HighlightFalseBridges = hints.HighlightFalseBridges,
                ShowVectorField = hints.ShowVectorField,
            };

        // ---------------------------------------------------------------
        // VectorFieldLayer
        // ---------------------------------------------------------------

        private static List<VectorFieldLayerJson> SerializeVectorFieldLayers(
            IReadOnlyList<VectorFieldLayer> layers)
        {
            var result = new List<VectorFieldLayerJson>(layers.Count);
            foreach (var layer in layers)
                result.Add(new VectorFieldLayerJson
                {
                    Id = layer.Id,
                    Name = layer.Name,
                    N = layer.N,
                    D = layer.D,
                    Vectors = layer.Vectors.ToArray(),
                });
            return result;
        }

        // ---------------------------------------------------------------
        // LineFieldLayer  (unoriented direction per node — d and -d equivalent)
        // ---------------------------------------------------------------

        private static List<LineFieldLayerJson> SerializeLineFieldLayers(
            IReadOnlyList<LineFieldLayer> layers)
        {
            var result = new List<LineFieldLayerJson>(layers.Count);
            foreach (var layer in layers)
                result.Add(new LineFieldLayerJson
                {
                    Id = layer.Id,
                    Name = layer.Name,
                    N = layer.N,
                    D = layer.D,
                    Source = layer.Source,
                    Directions = layer.Directions.ToArray(),
                });
            return result;
        }
    }

    // ===================================================================
    // JSON DTO types  (internal — consumers see only the JSON schema)
    // ===================================================================

    internal sealed class ScenePackageJson
    {
        [JsonPropertyName("schema_version")] public int SchemaVersion { get; set; }
        [JsonPropertyName("title")] public string Title { get; set; } = "";
        [JsonPropertyName("points")] public PointCloudJson Points { get; set; } = null!;
        [JsonPropertyName("label_layers")] public List<LabelLayerJson> LabelLayers { get; set; } = new();
        [JsonPropertyName("node_signal_layers")] public List<NodeSignalLayerJson> NodeSignalLayers { get; set; } = new();
        [JsonPropertyName("edge_layers")] public List<EdgeLayerJson> EdgeLayers { get; set; } = new();
        [JsonPropertyName("triangle_layers")] public List<TriangleLayerJson> TriangleLayers { get; set; } = new();
        [JsonPropertyName("gaussian_layers")] public List<GaussianLayerJson> GaussianLayers { get; set; } = new();
        [JsonPropertyName("spine_layers")] public List<SpineLayerJson> SpineLayers { get; set; } = new();
        [JsonPropertyName("vector_field_layers")] public List<VectorFieldLayerJson> VectorFieldLayers { get; set; } = new();
        [JsonPropertyName("line_field_layers")] public List<LineFieldLayerJson> LineFieldLayers { get; set; } = new();
        [JsonPropertyName("temporal_sequence")] public TemporalSequenceJson? TemporalSequence { get; set; }
        [JsonPropertyName("temporal_frame_index")] public int? TemporalFrameIndex { get; set; }
        [JsonPropertyName("hints")] public SceneRenderHintsJson Hints { get; set; } = null!;
        [JsonPropertyName("generator_params")] public Dictionary<string, object>? GeneratorParams { get; set; }
        [JsonPropertyName("generator_param_schema")] public GeneratorParamSchema? GeneratorParamSchema { get; set; }
        /// <summary>All registered generator names in display order. Populated from SchemaCatalog.KnownGenerators.</summary>
        [JsonPropertyName("generator_catalog")] public string[] GeneratorCatalog { get; set; } = Array.Empty<string>();
    }

    internal sealed class PointCloudJson
    {
        [JsonPropertyName("n")] public int N { get; set; }
        [JsonPropertyName("d")] public int D { get; set; }
        [JsonPropertyName("label")] public string? Label { get; set; }
        /// <summary>Flat [N×D], stride: point i dim d → i*D+d.</summary>
        [JsonPropertyName("features")] public double[] Features { get; set; } = Array.Empty<double>();
    }

    internal sealed class LabelLayerJson
    {
        [JsonPropertyName("id")] public string Id { get; set; } = "";
        [JsonPropertyName("name")] public string Name { get; set; } = "";
        [JsonPropertyName("kind")] public LabelLayerKind Kind { get; set; }
        /// <summary>Length N. -1 = unassigned/noise.</summary>
        [JsonPropertyName("labels")] public int[] Labels { get; set; } = Array.Empty<int>();
    }

    internal sealed class NodeSignalLayerJson
    {
        [JsonPropertyName("id")] public string Id { get; set; } = "";
        [JsonPropertyName("name")] public string Name { get; set; } = "";
        /// <summary>Semantic origin of the signal — filter, eigenfunction, susceptibility, etc.</summary>
        [JsonPropertyName("source")] public ScalarSource Source { get; set; }
        [JsonPropertyName("values")] public double[] Values { get; set; } = Array.Empty<double>();
        /// <summary>Pre-computed min value. C# bakes these so JS never has to scan the array.</summary>
        [JsonPropertyName("min")] public double Min { get; set; }
        /// <summary>Pre-computed max value.</summary>
        [JsonPropertyName("max")] public double Max { get; set; }
    }

    internal sealed class EdgeLayerJson
    {
        [JsonPropertyName("id")] public string Id { get; set; } = "";
        [JsonPropertyName("name")] public string Name { get; set; } = "";
        /// <summary>Canonical metric spec used to build this graph (for example "euclidean" or "minkowski:p=2"). Null = unknown/custom.</summary>
        [JsonPropertyName("metric")] public string? Metric { get; set; }
        /// <summary>Graph construction rule with type-specific parameters. Null = unknown/custom.</summary>
        [JsonPropertyName("proximity")] public ProximitySpecJson? Proximity { get; set; }
        [JsonPropertyName("src")] public int[] Src { get; set; } = Array.Empty<int>();
        [JsonPropertyName("dst")] public int[] Dst { get; set; } = Array.Empty<int>();
        [JsonPropertyName("weight")] public double[] Weight { get; set; } = Array.Empty<double>();
        [JsonPropertyName("edge_cluster_src")] public int[]? EdgeClusterSrc { get; set; }
        [JsonPropertyName("edge_cluster_dst")] public int[]? EdgeClusterDst { get; set; }
    }

    internal sealed class TriangleLayerJson
    {
        [JsonPropertyName("id")] public string Id { get; set; } = "";
        [JsonPropertyName("name")] public string Name { get; set; } = "";
        [JsonPropertyName("source")] public TriangleSource Source { get; set; }
        [JsonPropertyName("source_edge_layer_id")] public string SourceEdgeLayerId { get; set; } = "";
        /// <summary>Flat [T×3] array of sorted vertex triples.</summary>
        [JsonPropertyName("vertices")] public int[] Vertices { get; set; } = Array.Empty<int>();
    }

    /// <summary>
    /// Flat DTO for a ProximitySpec discriminated union.
    /// <c>kind</c> identifies the subtype; only the relevant parameter field is non-null.
    /// </summary>
    internal sealed class ProximitySpecJson
    {
        [JsonPropertyName("kind")] public NeighborRule Kind { get; set; }
        /// <summary>Present for Knn, MutualKnn, MstAugmented. Null for EpsilonBall.</summary>
        [JsonPropertyName("k")] public int? K { get; set; }
        /// <summary>Present for EpsilonBall only. Null for all KNN variants.</summary>
        [JsonPropertyName("epsilon")] public double? Epsilon { get; set; }
    }

    internal sealed class GaussianLayerJson
    {
        [JsonPropertyName("id")] public string Id { get; set; } = "";
        [JsonPropertyName("name")] public string Name { get; set; } = "";
        [JsonPropertyName("k")] public int K { get; set; }
        [JsonPropertyName("d")] public int D { get; set; }
        /// <summary>Flat [K×D]. Mean of component k: k*D .. k*D+D-1.</summary>
        [JsonPropertyName("means")] public double[] Means { get; set; } = Array.Empty<double>();
        /// <summary>Flat [K×D×D]. Cov of component k, row r, col c: k*D*D + r*D + c.</summary>
        [JsonPropertyName("covariances")] public double[] Covariances { get; set; } = Array.Empty<double>();
        /// <summary>Flat [K×3×3] lower-triangular Cholesky factor L (Σ=LLᵀ), geometry-space (always 3×3).
        /// Component k, row r, col c: k*9 + r*3 + c. Upper triangle is 0.
        /// JS reads this directly into Matrix4 to stretch a unit sphere — no math in the browser.</summary>
        [JsonPropertyName("cholesky_l")] public double[] CholeskyL { get; set; } = Array.Empty<double>();
        [JsonPropertyName("weights")] public double[] Weights { get; set; } = Array.Empty<double>();
        /// <summary>null = flat GMM. Non-null = ComponentToClusterMap[k] → cluster id.</summary>
        [JsonPropertyName("component_to_cluster_map")] public int[]? ComponentToClusterMap { get; set; }
    }

    internal sealed class SpineLayerJson
    {
        [JsonPropertyName("id")] public string Id { get; set; } = "";
        [JsonPropertyName("name")] public string Name { get; set; } = "";
        [JsonPropertyName("cluster_idx")] public int ClusterIdx { get; set; }
        [JsonPropertyName("kind")] public SpineLayerKind Kind { get; set; }
        [JsonPropertyName("spine_samples_shape")] public int[] SpineSamplesShape { get; set; } = Array.Empty<int>();
        /// <summary>Flat [M×D]. Sample i, dim d: i*D+d.</summary>
        [JsonPropertyName("spine_samples")] public double[] SpineSamples { get; set; } = Array.Empty<double>();
        [JsonPropertyName("tangent_bases_shape")] public int[]? TangentBasesShape { get; set; }
        /// <summary>Flat [M×D×D]. null when Kind == Arc.</summary>
        [JsonPropertyName("tangent_bases")] public double[]? TangentBases { get; set; }
    }

    internal sealed class TemporalSequenceJson
    {
        [JsonPropertyName("name")] public string Name { get; set; } = "";
        [JsonPropertyName("axis")] public TemporalAxis Axis { get; set; }
        [JsonPropertyName("index_values")] public double[] IndexValues { get; set; } = Array.Empty<double>();
        [JsonPropertyName("frames")] public List<LabelLayerJson> Frames { get; set; } = new();
    }

    internal sealed class SceneRenderHintsJson
    {
        [JsonPropertyName("show_edge_weights_as_opacity")] public bool ShowEdgeWeightsAsOpacity { get; set; }
        [JsonPropertyName("show_gaussian_ellipsoids")] public bool ShowGaussianEllipsoids { get; set; }
        [JsonPropertyName("show_spine_overlays")] public bool ShowSpineOverlays { get; set; }
        [JsonPropertyName("show_tangent_bases")] public bool ShowTangentBases { get; set; }
        [JsonPropertyName("overlay_component_and_cluster_coloring")] public bool OverlayComponentAndClusterColoring { get; set; }
        [JsonPropertyName("annotate_spin_color_vs_equilibrium")] public bool AnnotateSpinColorVsEquilibrium { get; set; }
        [JsonPropertyName("highlight_false_bridges")] public bool HighlightFalseBridges { get; set; }
        [JsonPropertyName("show_vector_field")] public bool ShowVectorField { get; set; }
    }

    internal sealed class VectorFieldLayerJson
    {
        [JsonPropertyName("id")] public string Id { get; set; } = "";
        [JsonPropertyName("name")] public string Name { get; set; } = "";
        [JsonPropertyName("n")] public int N { get; set; }
        [JsonPropertyName("d")] public int D { get; set; }
        /// <summary>Flat [N×D]. Slice [i*D .. i*D+D) is unit tangent at point i. Zero = undefined.</summary>
        [JsonPropertyName("vectors")] public double[] Vectors { get; set; } = Array.Empty<double>();
    }

    /// <summary>
    /// Unoriented direction per node — d and -d are equivalent. Viewer must render as
    /// lines / cylinders / capsules, not arrows. Distinct from VectorFieldLayer (oriented)
    /// and from SpineLayer.TangentBases (curve-anchored, not graph-anchored).
    /// </summary>
    internal sealed class LineFieldLayerJson
    {
        [JsonPropertyName("id")] public string Id { get; set; } = "";
        [JsonPropertyName("name")] public string Name { get; set; } = "";
        [JsonPropertyName("n")] public int N { get; set; }
        [JsonPropertyName("d")] public int D { get; set; }
        /// <summary>Semantic origin — spectral gradient, local PCA, etc.</summary>
        [JsonPropertyName("source")] public LineFieldSource Source { get; set; }
        /// <summary>Flat [N×D]. Slice [i*D .. i*D+D) is the unoriented direction at point i. Zero = undefined.</summary>
        [JsonPropertyName("directions")] public double[] Directions { get; set; } = Array.Empty<double>();
    }
}
