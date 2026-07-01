// src/graphs/spectral/MagneticLaplacianOperator.cs
#nullable enable
using System;
using System.Collections.Generic;
using Graphs.Primitives;
using Maths.LinAlg;

namespace Graphs.Spectral;

/// <summary>
/// The magnetic (U(1)-connection) graph Laplacian
/// <c>L^(q) = D − A ∘ exp(i·2πq·Θ)</c> as a matrix-free <see cref="ILinearOperator"/>
/// for <see cref="LOBPCG"/>. <c>D</c> is the symmetrized degree (neighbour count),
/// <c>A</c> the 0/1 symmetrized adjacency, <c>Θ</c> an antisymmetric per-edge
/// orientation (<c>+1</c> for <c>i→j</c>, <c>−1</c> for <c>j→i</c>), and <c>q</c>
/// the charge / field-strength knob. <c>L^(q)</c> is Hermitian PSD (real
/// nonnegative spectrum, Hodge theory intact) and collapses to the ordinary
/// undirected Laplacian at <c>q=0</c>.
/// </summary>
/// <remarks>
/// <para>
/// <b>Key trick — no complex eigensolver.</b> A Hermitian <c>H = A_re + i·B_im</c>
/// (<c>A_re</c> symmetric, <c>B_im</c> antisymmetric) embeds into the real symmetric
/// <c>2n×2n</c> matrix <c>[[A_re, −B_im], [B_im, A_re]]</c>, whose spectrum is exactly
/// that of <c>H</c> with every eigenvalue doubled (eigenvectors <c>[x;y]</c> and
/// <c>[−y;x]</c>). So this operator presents <see cref="Dimension"/> = <c>2n</c> and
/// applies that embedding matrix-free — the existing real <see cref="LOBPCG"/> solves
/// the complex Hermitian eigenproblem unchanged. The block layout per column is
/// <c>[x₀..x₍ₙ₋₁₎, y₀..y₍ₙ₋₁₎]</c> (real part stacked over imaginary part).
/// </para>
/// <para>
/// Seam A of the complex-analytic TDA thread: a U(1)-phase "v1.5" that lets H₁-like
/// cycles feel backbone direction (a directed return encloses magnetic flux) while
/// keeping a clean real spectrum — between the undirected post-check (v1) and full
/// directed path homology (v2). The phase rides a <em>connection</em>, not the chain
/// coefficients (which are topologically inert). Refs: Lieb–Loss 1993 / Sunada 1994
/// (discrete magnetic Laplacian); Fanuel et al., magnetic eigenmaps (ACHA 2018).
/// </para>
/// <para>
/// Scope: scalar U(1), unit edge magnitudes (phase-only), full-CSR streaming.
/// Hermitian half-storage, weighted magnitudes, and a <c>Vector512</c>/PB-CSR kernel
/// are deferred optimizations; vector stalks (U(n) connections) are out of scope.
/// </para>
/// </remarks>
public sealed class MagneticLaplacianOperator : ILinearOperator
{
    private readonly CsrGraph _graph;
    private readonly double[] _degree;   // symmetrized neighbour count, length n
    private readonly int[] _orientation; // per-CSR-slot Θ (+1/−1/0), antisymmetric
    private readonly double[] _cos;      // cos(2πq·Θ) per CSR slot
    private readonly double[] _sin;      // sin(2πq·Θ) per CSR slot
    private readonly double _charge;

    /// <summary>
    /// Wraps a symmetric <paramref name="graph"/> plus a per-CSR-slot orientation
    /// <paramref name="slotOrientation"/> (length = <c>graph.Targets.Length</c>;
    /// <c>+1</c>/<c>−1</c>/<c>0</c>, antisymmetric across reverse slots) at charge
    /// <paramref name="charge"/>. Phases <c>2πq·Θ</c> are precomputed once.
    /// </summary>
    public MagneticLaplacianOperator(CsrGraph graph, int[] slotOrientation, double charge)
    {
        if (slotOrientation is null) throw new ArgumentNullException(nameof(slotOrientation));
        if (slotOrientation.Length != (graph.Targets?.Length ?? 0))
            throw new ArgumentException(
                $"slotOrientation length ({slotOrientation.Length}) must equal CSR slot count ({graph.Targets?.Length ?? 0}).",
                nameof(slotOrientation));

        _graph = graph;
        _charge = charge;
        _orientation = (int[])slotOrientation.Clone();

        int n = graph.NodeCount;
        _degree = new double[n];
        for (int i = 0; i < n; i++)
            _degree[i] = graph.Degree(i); // 0/1 adjacency ⇒ degree = neighbour count

        int slots = slotOrientation.Length;
        _cos = new double[slots];
        _sin = new double[slots];
        double twoPiQ = 2.0 * Math.PI * charge;
        for (int e = 0; e < slots; e++)
        {
            double theta = twoPiQ * slotOrientation[e];
            _cos[e] = Math.Cos(theta);
            _sin[e] = Math.Sin(theta);
        }
    }

    /// <summary>The charge / field-strength knob <c>q</c>.</summary>
    public double Charge => _charge;

    /// <summary>Node count <c>n</c> of the underlying graph (operator dimension is <c>2n</c>).</summary>
    public int NodeCount => _graph.NodeCount;

    /// <summary>Embedding dimension <c>2n</c> — the operator acts on stacked <c>[x;y]</c> reals.</summary>
    public int Dimension => 2 * _graph.NodeCount;

    /// <summary>
    /// Builds an operator for a directed-edge list, the self-contained path (no Rips /
    /// proximity pipeline). Each pair <c>(from,to)</c> becomes a unit undirected edge
    /// carrying <c>Θ=+1</c> on the <c>from→to</c> slot and <c>Θ=−1</c> on its reverse.
    /// Edges must be distinct and must not include both <c>(u,v)</c> and <c>(v,u)</c>
    /// (a symmetric pair is undirected, <c>Θ=0</c>).
    /// </summary>
    public static MagneticLaplacianOperator FromDirectedEdges(
        int nodeCount, ReadOnlySpan<(int from, int to)> directedEdges, double charge)
    {
        var edges = new Edge[directedEdges.Length];
        for (int e = 0; e < directedEdges.Length; e++)
            edges[e] = new Edge(directedEdges[e].from, directedEdges[e].to, 1.0);

        CsrGraph graph = CsrGraph.FromEdges(edges, nodeCount);
        var orientation = new int[graph.Targets.Length]; // 0 = undirected

        foreach (var (from, to) in directedEdges)
        {
            orientation[graph.FindSlot(from, to)] = +1;
            orientation[graph.FindSlot(to, from)] = -1;
        }

        return new MagneticLaplacianOperator(graph, orientation, charge);
    }

    /// <summary>
    /// Builds an operator from a directed <paramref name="backbone"/> (the time arrow;
    /// each <c>(u,v)</c> carries <c>Θ=±1</c>) plus undirected similarity
    /// <paramref name="chords"/> (<c>Θ=0</c>). A return that runs forward along the
    /// backbone and jumps back via a chord encloses flux ∝ the backbone span it
    /// bridges (§4 reach grading). Edges across the union must be distinct.
    /// </summary>
    public static MagneticLaplacianOperator FromBackboneAndChords(
        int nodeCount,
        ReadOnlySpan<(int from, int to)> backbone,
        ReadOnlySpan<(int a, int b)> chords,
        double charge)
    {
        var edges = new Edge[backbone.Length + chords.Length];
        int k = 0;
        foreach (var (from, to) in backbone) edges[k++] = new Edge(from, to, 1.0);
        foreach (var (a, b) in chords) edges[k++] = new Edge(a, b, 1.0);

        CsrGraph graph = CsrGraph.FromEdges(edges, nodeCount);
        var orientation = new int[graph.Targets.Length]; // chords stay 0 (undirected)

        foreach (var (from, to) in backbone)
        {
            orientation[graph.FindSlot(from, to)] = +1;
            orientation[graph.FindSlot(to, from)] = -1;
        }

        return new MagneticLaplacianOperator(graph, orientation, charge);
    }

    public void Apply(ReadOnlySpan<double> block, Span<double> result, int columns)
    {
        int n = _graph.NodeCount;
        int twoN = 2 * n;
        ReadOnlySpan<int> rowPtrs = _graph.RowPointers.AsSpan();
        ReadOnlySpan<int> targets = _graph.Targets.AsSpan();

        for (int c = 0; c < columns; c++)
        {
            int off = c * twoN;          // real part at [off, off+n), imag at [off+n, off+2n)
            for (int i = 0; i < n; i++)
            {
                double xi = block[off + i];
                double yi = block[off + n + i];
                double re = _degree[i] * xi;   // diagonal D is real
                double im = _degree[i] * yi;

                int start = rowPtrs[i];
                int end = rowPtrs[i + 1];
                for (int e = start; e < end; e++)
                {
                    int j = targets[e];
                    double cos = _cos[e];
                    double sin = _sin[e];
                    double xj = block[off + j];
                    double yj = block[off + n + j];

                    // Off-diagonal H_ij = −e^{iθ} = −(cos + i·sin); H_ij·z_j adds to (Hz)_i:
                    //   real:  −cos·xj + sin·yj
                    //   imag:  −sin·xj − cos·yj
                    re += -cos * xj + sin * yj;
                    im += -sin * xj - cos * yj;
                }

                result[off + i] = re;
                result[off + n + i] = im;
            }
        }
    }

    /// <summary>
    /// Materializes the real <c>2n×2n</c> symmetric embedding as a dense
    /// <c>double[,]</c> — the dense oracle peer of <see cref="Apply"/>, for parity
    /// tests and small-problem cross-checks against <see cref="SpectralMath"/>.
    /// </summary>
    public double[,] ToDenseEmbedding()
    {
        int n = _graph.NodeCount;
        int twoN = 2 * n;
        var m = new double[twoN, twoN];

        for (int i = 0; i < n; i++)
        {
            m[i, i] = _degree[i];           // A_re diagonal (top-left block)
            m[n + i, n + i] = _degree[i];   // A_re diagonal (bottom-right block)

            int start = _graph.RowPointers[i];
            int end = _graph.RowPointers[i + 1];
            for (int e = start; e < end; e++)
            {
                int j = _graph.Targets[e];
                double aRe = -_cos[e];      // A_re off-diagonal = −cos θ
                double bIm = -_sin[e];      // B_im off-diagonal = −sin θ
                m[i, j] = aRe;              // [[A_re, −B_im], [B_im, A_re]]
                m[n + i, n + j] = aRe;
                m[i, n + j] = -bIm;
                m[n + i, j] = bIm;
            }
        }

        return m;
    }

    /// <summary>
    /// Signed phase <c>θ = 2πq·Θ</c> (radians) on the directed slot
    /// <c>from→to</c>; throws if no such edge exists.
    /// </summary>
    public double EdgePhase(int from, int to)
        => 2.0 * Math.PI * _charge * _orientation[_graph.FindSlot(from, to)];

    /// <summary>
    /// Magnetic flux enclosed by a closed walk — the Aharonov–Bohm holonomy divided
    /// by <c>2π</c> (holonomy = <c>exp(i·2π·flux)</c>) — summing <see cref="EdgePhase"/>
    /// over consecutive nodes and the wrap from last back to first. A coherently
    /// directed loop accumulates flux ∝ its backbone span (reach); an out-and-back
    /// excursion cancels to zero. This is the discriminator a revisitation count
    /// cannot fake (§3/§4); it vanishes identically at <c>q=0</c>.
    /// </summary>
    public double EnclosedFlux(ReadOnlySpan<int> closedWalk)
    {
        if (closedWalk.Length < 2) return 0.0;

        double phase = 0.0;
        for (int s = 0; s < closedWalk.Length; s++)
            phase += EdgePhase(closedWalk[s], closedWalk[(s + 1) % closedWalk.Length]);

        return phase / (2.0 * Math.PI);
    }
}

/// <summary>
/// Spectral entry points for the magnetic Laplacian, built on the pure
/// <see cref="LOBPCG"/> primitive over the real <c>2n</c> embedding.
/// </summary>
public static class MagneticSpectral
{
    /// <summary>
    /// The <paramref name="k"/> smallest <em>distinct</em> eigenvalues of
    /// <c>L^(q)</c>, ascending. The embedding doubles every eigenvalue, so this
    /// solves for <c>2k</c> pairs and strides by two. The smallest value lifting off
    /// zero at fractional <paramref name="op"/> charge is the directed cycle being
    /// frustrated by the field (Aharonov–Bohm flux).
    /// </summary>
    public static double[] BottomKEigenvalues(
        MagneticLaplacianOperator op, int k, LOBPCG.Options? options = null)
    {
        if (op is null) throw new ArgumentNullException(nameof(op));
        if (k <= 0) return Array.Empty<double>();

        LOBPCG.Result result = LOBPCG.Solve(op, 2 * k, options);
        IReadOnlyList<EigenPair> pairs = result.Eigenpairs; // ascending

        var values = new List<double>(k);
        for (int i = 0; i < pairs.Count && values.Count < k; i += 2)
            values.Add(pairs[i].Lambda);

        return values.ToArray();
    }
}
