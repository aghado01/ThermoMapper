// src/graphs/spectral/GraphSpectralOperators.cs
#nullable enable
using System;
using Graphs.Primitives;
using Maths.LinAlg;

namespace Graphs.Spectral;

/// <summary>
/// Matrix-free graph Laplacian as an <see cref="ILinearOperator"/> for
/// <see cref="LOBPCG"/>. Precomputes the (clamped) weighted degree and its
/// inverse square root once; <see cref="Apply"/> then streams the CSR structure
/// per column with no further allocation.
/// </summary>
/// <remarks>
/// The matrix realized here is bit-for-bit the same as
/// <see cref="GraphLaplacian.BuildDense"/>:
/// <list type="bullet">
/// <item><see cref="LaplacianType.Combinatorial"/> = <c>D − W</c>, with
/// non-positive edge weights clamped to <c>1.0</c> on <em>both</em> the diagonal
/// and the off-diagonal.</item>
/// <item><see cref="LaplacianType.NormalizedSymmetric"/> = <c>I − D^(-1/2) W D^(-1/2)</c>,
/// using raw edge weights off-diagonal (degree clamp only feeds <c>D</c>).</item>
/// </list>
/// Keeping the two construction paths in lockstep is what lets LOBPCG and the
/// dense eigensolver agree on the same graph.
/// </remarks>
public sealed class GraphLaplacianOperator : ILinearOperator
{
    private readonly CsrGraph _graph;
    private readonly LaplacianType _lapType;
    private readonly double[] _degree;
    private readonly double[] _invSqrtDegree;

    public GraphLaplacianOperator(CsrGraph graph, LaplacianType lapType)
    {
        _graph = graph;
        _lapType = lapType;

        int n = graph.NodeCount;
        _degree = new double[n];
        for (int i = 0; i < n; i++)
        {
            int start = graph.RowPointers[i];
            int end = graph.RowPointers[i + 1];
            double sum = 0.0;
            for (int e = start; e < end; e++)
                sum += graph.Weights[e] > 0.0 ? graph.Weights[e] : 1.0;
            _degree[i] = sum;
        }

        _invSqrtDegree = new double[n];
        for (int i = 0; i < n; i++)
            _invSqrtDegree[i] = _degree[i] > 1e-12 ? 1.0 / Math.Sqrt(_degree[i]) : 0.0;
    }

    /// <summary>Clamped weighted degree per node — the Laplacian diagonal <c>D</c>.</summary>
    public double[] Degree => _degree;

    /// <summary>Inverse square-root of <see cref="Degree"/> (0 where degree is ~0).</summary>
    public double[] InverseSqrtDegree => _invSqrtDegree;

    public int Dimension => _graph.NodeCount;

    public void Apply(ReadOnlySpan<double> block, Span<double> result, int columns)
    {
        int n = _graph.NodeCount;
        ReadOnlySpan<int> rowPtrs = _graph.RowPointers.AsSpan();
        ReadOnlySpan<int> targets = _graph.Targets.AsSpan();
        ReadOnlySpan<double> weights = _graph.Weights.AsSpan();

        for (int j = 0; j < columns; j++)
        {
            int colOffset = j * n;

            for (int i = 0; i < n; i++)
            {
                if (_lapType == LaplacianType.NormalizedSymmetric)
                {
                    if (_invSqrtDegree[i] == 0.0)
                    {
                        result[colOffset + i] = block[colOffset + i];
                        continue;
                    }

                    double neighborSum = 0.0;
                    for (int e = rowPtrs[i]; e < rowPtrs[i + 1]; e++)
                    {
                        int target = targets[e];
                        neighborSum += weights[e] * _invSqrtDegree[target] * block[colOffset + target];
                    }

                    result[colOffset + i] = block[colOffset + i] - _invSqrtDegree[i] * neighborSum;
                }
                else
                {
                    double val = _degree[i] * block[colOffset + i];

                    int start = rowPtrs[i];
                    int end = rowPtrs[i + 1];
                    for (int e = start; e < end; e++)
                    {
                        // Clamp to match the diagonal degree and GraphLaplacian.BuildDense.
                        double w = weights[e] > 0.0 ? weights[e] : 1.0;
                        val -= w * block[colOffset + targets[e]];
                    }

                    result[colOffset + i] = val;
                }
            }
        }
    }
}

/// <summary>
/// Jacobi (diagonal) preconditioner for a graph Laplacian:
/// <c>M⁻¹ = D⁻¹</c> for the combinatorial Laplacian, identity for the normalized
/// symmetric Laplacian (whose diagonal is already unit).
/// </summary>
public sealed class JacobiPreconditioner : IPreconditioner
{
    private readonly double[] _degree;
    private readonly LaplacianType _lapType;

    public JacobiPreconditioner(double[] degree, LaplacianType lapType)
    {
        _degree = degree ?? throw new ArgumentNullException(nameof(degree));
        _lapType = lapType;
    }

    public void Apply(ReadOnlySpan<double> residualBlock, Span<double> resultBlock, int columns)
    {
        int n = _degree.Length;

        if (_lapType == LaplacianType.NormalizedSymmetric)
        {
            residualBlock.Slice(0, n * columns).CopyTo(resultBlock);
            return;
        }

        for (int j = 0; j < columns; j++)
        {
            int offset = j * n;
            for (int i = 0; i < n; i++)
            {
                double d = _degree[i];
                if (d < 1e-12) d = 1.0;
                resultBlock[offset + i] = residualBlock[offset + i] / d;
            }
        }
    }
}

/// <summary>
/// Graph spectral entry points built on the pure <see cref="LOBPCG"/> primitive.
/// </summary>
public static class GraphSpectral
{
    /// <summary>
    /// Computes the <paramref name="k"/> smallest eigenpairs of the graph
    /// Laplacian (the low-frequency spectrum: algebraic connectivity, Fiedler
    /// vector, spectral-embedding coordinates).
    /// </summary>
    /// <remarks>
    /// When <paramref name="deflateNullSpace"/> is <see langword="true"/> (default),
    /// the trivial null mode — the constant vector for the combinatorial Laplacian,
    /// <c>D^(1/2)·1</c> for the normalized symmetric one — is deflated as a LOBPCG
    /// constraint, so the returned pairs are the smallest <em>non-trivial</em> modes
    /// (the embedding coordinates) rather than the λ≈0 null eigenvector. Pass
    /// <see langword="false"/> to include the null mode and recover the same set as
    /// a dense bottom-K decomposition.
    /// </remarks>
    public static LOBPCG.Result ComputeBottomK(
        CsrGraph graph,
        int k = 8,
        int maxIter = 300,
        double tol = 1e-9,
        LaplacianType lapType = LaplacianType.Combinatorial,
        int seed = 12345,
        bool deflateNullSpace = true)
    {
        if (graph.NodeCount == 0 || k <= 0)
            return new LOBPCG.Result(Array.Empty<EigenPair>(), 0, 0.0, true);

        var op = new GraphLaplacianOperator(graph, lapType);
        var preconditioner = new JacobiPreconditioner(op.Degree, lapType);
        double[]? nullVector = deflateNullSpace ? BuildTrivialNullVector(op.Degree, lapType) : null;

        var options = new LOBPCG.Options
        {
            MaxIterations = maxIter,
            Tolerance = tol,
            WantLargest = false,
            Constraints = nullVector,
            ConstraintColumns = nullVector is null ? 0 : 1,
            Preconditioner = preconditioner,
            Seed = seed,
        };

        return LOBPCG.Solve(op, k, options);
    }

    /// <summary>
    /// Builds the normalized trivial null vector of the chosen Laplacian, or
    /// <see langword="null"/> for an empty / degenerate graph.
    /// </summary>
    private static double[]? BuildTrivialNullVector(double[] degree, LaplacianType type)
    {
        int n = degree.Length;
        if (n == 0) return null;

        var nullVector = new double[n];
        bool anyNonZero = false;
        for (int i = 0; i < n; i++)
        {
            // Normalized symmetric: null vector is D^(1/2)·1; combinatorial: constant 1.
            double value = type == LaplacianType.NormalizedSymmetric ? Math.Sqrt(degree[i]) : 1.0;
            nullVector[i] = value;
            anyNonZero |= value != 0.0;
        }

        if (!anyNonZero) return null;

        double norm = Math.Sqrt(MatrixOps.Dot(nullVector, nullVector));
        if (norm < 1e-12) return null;

        double invNorm = 1.0 / norm;
        for (int i = 0; i < n; i++)
            nullVector[i] *= invNorm;

        return nullVector;
    }
}
