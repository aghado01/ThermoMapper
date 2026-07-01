using System;
using Graphs.Primitives;
using Maths.LinAlg;

namespace Graphs.Observables;

/// <summary>
/// Per-node entry of the top eigenvector of the alignment matrix
/// <c>G_ij = ⟨δ(s_i, s_j)⟩</c> — the <see cref="Alignments"/> currency densified
/// (unit diagonal, the correlation-matrix convention). Reads as eigenvector
/// centrality on the equilibrium alignment graph: nodes whose entries are large
/// participate strongly in the dominant collective mode of aligned spin
/// fluctuations at equilibrium.
/// </summary>
/// <remarks>
/// <para><b>Currency.</b> Consumes <see cref="Alignments"/>, the SW-native
/// spin-alignment channel <c>⟨δ(s_i,s_j)⟩</c> — <i>not</i>
/// <see cref="Affinities"/> (bond-survival). A forward solver (PKWang) draws no
/// spins and mints no alignment, so this signal is meaningful only for
/// samplers that materialize a spin ensemble.</para>
///
/// <para><b>Algorithm.</b>
/// <list type="number">
///   <item>Build the dense N×N symmetric matrix <c>G</c>:
///     <list type="bullet">
///       <item><c>G[i,j] = G[j,i] = Alignments.G[e]</c> for every undirected
///         edge <c>(i,j)</c> in the CSR.</item>
///       <item><c>G[i,i] = 1</c> (a node trivially aligns with itself — the
///         correlation-matrix convention).</item>
///       <item>Off-graph pairs default to 0 (no observation; FK-zero baseline).
///         Pairs whose mutual alignment is significant but whose edge was
///         absent from the proximity graph cannot be recovered from this
///         signal.</item>
///     </list>
///   </item>
///   <item>Call <see cref="DenseEigen.DecomposeSymmetric(double[,], int, double, DenseEigenOptions)"/>;
///     the result's eigenvalues are sorted descending.</item>
///   <item>Take <c>Eigenvectors[0]</c> as the per-node scalar.</item>
///   <item>By Perron-Frobenius (G is non-negative symmetric) the top
///     eigenvector is sign-uniform. If the solver returned the negative branch,
///     flip the sign so all entries are non-negative — gives callers a
///     canonical, interpretable centrality.</item>
/// </list>
/// </para>
///
/// <para><b>Subject / Op decomposition.</b>
/// Subject: <c>Alignment</c> (per-pair <c>G_ij = ⟨δ(s_i,s_j)⟩</c>).
/// Op: <c>EigenCentrality</c> — eigenvector centrality (the top eigenvector of
/// the symmetric non-negative matrix). The Op commits to eigenvector centrality
/// specifically — not betweenness, not Katz, not PageRank — and not merely the
/// abstract category <c>Centrality</c>.</para>
///
/// <para><b>Cost.</b> Dense <c>O(N²)</c> memory and <c>O(N³)</c> compute for the
/// eigendecomposition. Suited to typical patch sizes (N ≲ a few hundred); for
/// global SPC on large graphs a Lanczos-style sparse iteration on the symmetric
/// alignment operator would scale better but is not in the current LinAlg
/// surface — revisit if profiling shows this hot.</para>
/// </remarks>
public sealed class AlignmentEigenCentrality : IGraphSignal<Alignments>
{
    /// <summary>
    /// Maximum Jacobi sweeps for the symmetric eigendecomposition.
    /// Default matches <see cref="DenseEigen"/>'s default.
    /// </summary>
    public int MaxSweeps { get; init; } = 256;

    /// <summary>
    /// Convergence tolerance for the off-diagonal infinity norm.
    /// Default matches <see cref="DenseEigen"/>'s default.
    /// </summary>
    public double Tolerance { get; init; } = 1e-12;

    /// <summary>
    /// Optional eigensolver options (e.g. fast-variant selection). Default
    /// uses the dispatcher's defaults.
    /// </summary>
    public DenseEigenOptions EigenOptions { get; init; } = default;

    public double[] Compute(Alignments alignments, CsrGraph graph)
    {
        ArgumentNullException.ThrowIfNull(alignments);
        ArgumentNullException.ThrowIfNull(graph);

        double[] align = alignments.G;
        if (align.Length != graph.Targets.Length)
            throw new ArgumentException(
                $"Alignments.G length ({align.Length}) does not match CSR slot count " +
                $"({graph.Targets.Length}).", nameof(alignments));

        int n = graph.NodeCount;
        if (n == 0) return Array.Empty<double>();
        if (n == 1) return new double[] { 1.0 };

        // ── Build dense symmetric G ──────────────────────────────────────
        var g = new double[n, n];
        for (int i = 0; i < n; i++)
            g[i, i] = 1.0;  // self-alignment (correlation-matrix convention)

        // Only the upper-triangular CSR slots (j > i) carry currency; mirror to both.
        for (int i = 0; i < n; i++)
        {
            int rowEnd = graph.RowPointers[i + 1];
            for (int e = graph.RowPointers[i]; e < rowEnd; e++)
            {
                int j = graph.Targets[e];
                if (j <= i) continue;
                double gij = align[e];
                g[i, j] = gij;
                g[j, i] = gij;
            }
        }

        // ── Top eigenvector via dispatched symmetric eigendecomposition ──
        EigenResult eig = DenseEigen.DecomposeSymmetric(g, MaxSweeps, Tolerance, EigenOptions);
        double[] top = eig.Eigenvectors[0];
        if (top.Length != n)
            throw new InvalidOperationException(
                $"Eigenvector length ({top.Length}) does not match node count ({n}).");

        // ── Canonicalize sign via Perron-Frobenius (non-negative branch) ─
        // For a non-negative symmetric matrix the leading eigenvector is
        // sign-uniform; if the solver returned the negative branch, flip.
        double signMass = 0.0;
        for (int i = 0; i < n; i++) signMass += top[i];
        if (signMass < 0.0)
        {
            var flipped = new double[n];
            for (int i = 0; i < n; i++) flipped[i] = -top[i];
            return flipped;
        }
        // Defensive copy: top is owned by the EigenResult and may be shared
        // with other consumers of the same decomposition.
        var result = new double[n];
        Array.Copy(top, result, n);
        return result;
    }
}
