using Graphs;
using Graphs.Distance;
using TDA.Ph;

namespace TDA.DimReduction;

/// <summary>
/// Declarative recipe for the SPRED persistent-homology objective (Yu &amp; You, arXiv:2106.02096).
/// Consumed by <see cref="PersistenceObjective"/>; a fluent shell (if any) belongs at the CLI/REPL
/// boundary, not here.
/// </summary>
public sealed record PersistenceObjectiveConfig
{
    /// <summary><b>Required.</b> The graph-construction recipe applied to the projected cloud (and, by
    /// default, to the ambient reference). Its projection stage should emit a distance graph when
    /// <see cref="Filtration"/> is <see cref="FiltrationWeights.RawDistance"/>.</summary>
    public required GraphCompilerConfig Graph { get; init; }

    /// <summary>Optional distinct recipe for the ambient reference barcode. <c>null</c> reuses
    /// <see cref="Graph"/> — the paper-faithful "same construction on both sides".</summary>
    public GraphCompilerConfig? ReferenceGraph { get; init; }

    /// <summary>Filtration-value source for the Rips complex: raw edge distance, or Laplacian-derived
    /// effective resistance.</summary>
    public FiltrationWeights Filtration { get; init; } = FiltrationWeights.RawDistance;

    /// <summary>Distance metric on the projected space. <c>null</c> = Euclidean.</summary>
    public IDistanceMetric? ProjectedMetric { get; init; }

    /// <summary>Max simplex dimension the Rips complex builds (2 = H0+H1; loops need triangle
    /// fillers).</summary>
    public int MaxDimension { get; init; } = 2;

    /// <summary>Homological dimensions matched, with weights — the paper's <c>(λ, 1−λ)</c> multi-order
    /// combination (§6) generalized to arbitrary orders/weights.</summary>
    public (int Dim, double Weight)[] Dimensions { get; init; } = [(0, 0.5), (1, 0.5)];

    /// <summary>Wasserstein order <c>p</c>. (The ground metric <c>q</c> is fixed at L∞ by
    /// <see cref="DiagramMetrics"/>.)</summary>
    public double WassersteinOrder { get; init; } = 2.0;

    /// <summary>Diagram-distance backend for the per-dimension matching. Exact Hungarian matching
    /// is O((bars)³) and is the block-scale wall (ISOLET brief, "H0 matching-cost gate");
    /// <see cref="DiagramDistanceKind.SlicedWasserstein"/> and
    /// <see cref="DiagramDistanceKind.SinkhornWasserstein"/> are the screening-scale
    /// alternatives.</summary>
    public DiagramDistanceKind DiagramDistance { get; init; } = DiagramDistanceKind.Wasserstein;

    /// <summary>Slice count for <see cref="DiagramDistanceKind.SlicedWasserstein"/> — evenly
    /// spaced, deterministic (no RNG).</summary>
    public int SlicedDirections { get; init; } = 50;

    /// <summary>Entropic regularization for <see cref="DiagramDistanceKind.SinkhornWasserstein"/>,
    /// dimensionless (relative to the largest finite matching cost). Smaller = closer to exact,
    /// slower to converge.</summary>
    public double SinkhornEpsilon { get; init; } = 0.01;

    /// <summary>Iteration cap for <see cref="DiagramDistanceKind.SinkhornWasserstein"/>.</summary>
    public int SinkhornMaxIters { get; init; } = 500;

    /// <summary>Essential-bar policy. <c>null</c> auto-derives <see cref="DiagramMetrics.EssentialPolicy.FinitePenalty"/>
    /// at scale diam(X)/2, so an essential-count mismatch yields a finite (SA-usable) penalty rather
    /// than the <c>+∞</c> of <see cref="DiagramMetrics.EssentialPolicy.InfiniteOnMismatch"/>.</summary>
    public DiagramMetrics.EssentialPolicy? Essential { get; init; }

    /// <summary>Objective value returned when the compiler rejects a projected graph
    /// (<see cref="GraphPathologyException"/>), so the annealer rejects that proposal.</summary>
    public double PathologyPenalty { get; init; } = 1e6;

    /// <summary>PCA-spirit variance regularizer (§6): adds <c>w · tr(P Σ_X Pᵀ)</c>. <b>Negative</b>
    /// rewards variance (PCA maximizes that trace); default 0 (off).</summary>
    public double VarianceRegularizer { get; init; } = 0.0;

    /// <summary>Finite bars with persistence (death − birth) below this are pruned from both barcodes
    /// before Wasserstein matching. 0 = no pruning (exact). Raising it collapses the H1 Hungarian —
    /// the graph-restricted Rips emits O(n) near-diagonal noise-loop bars that dominate matching cost
    /// but barely move the Wasserstein value — so pruning is a near-exact speedup (and a denoiser).
    /// Essential (infinite) bars are always kept.</summary>
    public double MinPersistence { get; init; } = 0.0;
}
