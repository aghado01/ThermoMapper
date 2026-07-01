using Graphs.Primitives;

namespace Graphs.Observables;

/// <summary>
/// A per-node scalar field derived from a per-edge currency over a graph.
/// Mathematically a graph signal <c>f: V → ℝ</c> on the graph's vertex set,
/// computed from an edge field — a currency such as <see cref="Affinities"/> or
/// <see cref="Alignments"/> — together with the <see cref="CsrGraph"/>
/// structure. Each implementation names a distinct measurement; many are
/// routinely composed in a single analysis.
/// </summary>
/// <remarks>
/// <para><b>Tier.</b> Graph-intrinsic by signature: a signal reads only the
/// model-agnostic currency and the graph, never the sampler's accumulator or
/// run result. Any inference strategy that mints the currency (Swendsen–Wang,
/// PKWang, a GMM/HDBSCAN ensemble) feeds the same signal — that is what living
/// in <c>graphs/observables</c> buys.</para>
///
/// <para><b>Pattern role.</b> Unlike a partition strategy (alternatives — pick
/// one per analysis), graph signals are <i>distinct measurements</i> —
/// <see cref="AffinityDegree"/> and <see cref="AffinityBinaryEntropySum"/> do
/// not compete; they answer different questions about the same currency. The
/// "Strategy" suffix is intentionally omitted — these are not interchangeable
/// approaches to one job.</para>
///
/// <para><b>Naming convention for concrete signals.</b>
/// <c>&lt;Subject&gt;&lt;Op&gt;</c> — the family axis is the substrate being
/// measured (Affinity, Alignment, SpinHistogram, Cluster, ...), the
/// secondary axis is the mathematical operation that reduces that substrate to
/// one scalar per node (Degree, BinaryEntropySum, EigenvectorCentrality,
/// CyclomaticComplexity, Compressibility, ...). This inverts a partition
/// strategy's <c>&lt;Op&gt;&lt;Subject&gt;</c>: partitions vary by <i>how</i>,
/// signals vary by <i>what</i>. The naming order in each family puts the
/// varying axis first.</para>
///
/// <para><b>Op precision.</b> Op tokens name the full mathematical operation,
/// not a category: <c>EigenCentrality</c> not <c>Centrality</c>;
/// <c>CyclomaticComplexity</c> not <c>Cyclomatic</c>; <c>Compressibility</c>
/// (the property measured) rather than a specific estimator (LempelZiv etc.) —
/// estimator choice becomes a configuration parameter, not part of the type
/// identity.</para>
///
/// <para><b>Algebraic closure.</b> Graph signals are closed under pointwise
/// combination (sum, scale, normalize, product) of their per-node outputs. This
/// lets compound signals like "z-score of <see cref="AffinityDegree"/>" be
/// expressed as combinators over primitive signals; partitions do not compose
/// this way.</para>
///
/// <para><b>Output contract.</b> Returns a <c>double[]</c> of length
/// <c>graph.NodeCount</c>. Implementations validate the currency's per-slot
/// field against the CSR slot count and throw
/// <see cref="System.ArgumentException"/> on mismatch.</para>
/// </remarks>
/// <typeparam name="TCurrency">
/// The per-edge currency the signal reads — <see cref="Affinities"/> (universal
/// bond-survival) or <see cref="Alignments"/> (SW-native spin-alignment).
/// </typeparam>
public interface IGraphSignal<in TCurrency>
{
    double[] Compute(TCurrency currency, CsrGraph graph);
}
