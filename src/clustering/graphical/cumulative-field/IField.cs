using Graphs.Primitives;

namespace Clustering.Graphical.SPC.Runtime.Core.Solver;

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
    /// When <see langword="true"/>, both directed slots of each undirected edge
    /// carry a value and the kernel symmetrizes them before clustering
    /// (LocalField). MeanField is symmetric by construction and sets this false.
    /// </summary>
    static abstract bool DirectedSymmetrize { get; }
}
