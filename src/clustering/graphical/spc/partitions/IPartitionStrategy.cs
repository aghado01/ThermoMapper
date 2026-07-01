using Clustering.Primitives;
using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Partitions;

/// <summary>
/// A strategy for partitioning a graph's node set into discrete clusters
/// from the model-agnostic currencies minted at the chosen equilibrium
/// temperature. Implementations are <i>alternatives</i> — pick one per
/// analysis; replacing one with another gives a different answer to the
/// same question.
/// </summary>
/// <remarks>
/// <para><b>Pattern role.</b> Genuinely Strategy-shaped, unlike a graph
/// signal (<c>Graphs.Observables.IGraphSignal</c>, which is plural by nature —
/// many distinct measurements composed in one analysis).</para>
///
/// <para><b>Naming convention for concrete strategies.</b>
/// <c>&lt;Operation&gt;&lt;Subject&gt;</c> — the family axis is the
/// operation (Threshold, Spectral, Consensus, ...), the secondary axis
/// is the substrate the operation acts on (BondFrequency, SpinAgreement,
/// Alignment, ...). This asymmetry with a graph signal's
/// <c>&lt;Subject&gt;&lt;Op&gt;</c>
/// convention is deliberate: partition strategies vary by <i>how</i>;
/// signals vary by <i>what</i>. Naming order follows the conceptual
/// variation axis.</para>
///
/// <para><b>Output contract.</b> Returns a densely-labeled
/// <see cref="Assignment"/>. Implementations that require a currency
/// absent from the sweep throw <see cref="System.InvalidOperationException"/>
/// with a clear message (e.g. <see cref="Strategies.ThresholdSpinAgreement"/>
/// when <paramref name="alignments"/> is null).</para>
/// </remarks>
public interface IPartitionStrategy
{
    Assignment Apply(CsrGraph graph, Affinities affinities, Alignments? alignments, CoMembership? coMembership = null);
}
