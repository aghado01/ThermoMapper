using Graphs.Primitives;

namespace Clustering.Graphical.CumulativeField;

/// <summary>
/// Static-abstract struct interface (mirrors <c>ISwConfig</c>) carrying the
/// only field-specific step: building the per-CSR-slot cumulative-energy ladder
/// <c>Hcum</c>. The survival kernel that consumes <c>Hcum</c> is field-agnostic,
/// so the JIT specializes nothing beyond this method and the symmetrize flag.
/// </summary>
internal interface IField
{
    /// <summary>
    /// Build the temperature-independent <c>Hcum</c> ladder, indexed by CSR
    /// slot. Meaningful entries are written at the <c>j &gt; i</c> slots; the
    /// mirror half stays zero unless the field is directed
    /// (<see cref="DirectedSymmetrize"/>).
    /// </summary>
    static abstract double[] BuildHcum(CsrGraph graph);

    /// <summary>
    /// The same ladder in log space, <c>ln Hcum</c> by CSR slot, built from <see cref="CsrGraph.LogWeights"/> by a
    /// running log-sum-exp, so the rungs below the double range keep their values. A slot with no ladder entry is
    /// <c>−∞</c> (the linear ladder's 0). Requires <see cref="CsrGraph.LogWeights"/>.
    /// </summary>
    static abstract double[] BuildLogHcum(CsrGraph graph);

    /// <summary>
    /// When <see langword="true"/>, both directed slots of each undirected edge
    /// carry a value and the kernel symmetrizes them before clustering
    /// (LocalField). MeanField is symmetric by construction and sets this false.
    /// </summary>
    static abstract bool DirectedSymmetrize { get; }
}
