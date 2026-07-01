// src/graphs/spectral/GraphLaplacian.cs
#nullable enable
using System;
using System.Buffers;
using Graphs.Primitives;

namespace Graphs.Spectral;

public enum LaplacianType { Combinatorial, NormalizedSymmetric }

public enum DenseLaplacianMaterialization
{
    Rectangular = 0,
    FlatColumnMajor = 1,
}

/// <summary>
/// Materialises a <see cref="CsrGraph"/> as a dense graph Laplacian matrix.
/// Two layouts:
///   • <see cref="BuildDense"/> — <c>double[,]</c> row/col-indexed (idiomatic C#).
///   • <see cref="BuildDenseColumnMajor"/> — flat <c>Span&lt;double&gt;</c> column-major
///     (LAPACK / SIMD-eigenkernel layout; the input format <see cref="Maths.LinAlg.DenseEigen"/>
///     prefers when fed a span directly).
///
/// <see cref="LaplacianType.Combinatorial"/> = D − W (weighted degree on diagonal,
/// negative edge weights off-diagonal).<br/>
/// <see cref="LaplacianType.NormalizedSymmetric"/> = I − D^(-1/2) W D^(-1/2).
/// Use Normalized for graphs with highly variable degree (kernel-weighted
/// proximity graphs etc.); Combinatorial for uniformly-weighted or
/// degree-homogeneous graphs.
///
/// Zero or negative edge weights are clamped to 1.0 (treated as unweighted)
/// for the Combinatorial path; the Normalized path uses the raw weight values
/// as given.
/// </summary>
public static class GraphLaplacian
{
    public static double[,] BuildDense(CsrGraph graph, LaplacianType lapType)
    {
        int n = graph.NodeCount;
        var laplacian = new double[n, n];
        var degree = new double[n];

        for (int i = 0; i < n; i++)
        {
            int start = graph.RowPointers[i];
            int end = graph.RowPointers[i + 1];
            double sum = 0.0;
            for (int edge = start; edge < end; edge++)
                sum += graph.Weights[edge] > 0.0 ? graph.Weights[edge] : 1.0;

            degree[i] = sum;
        }

        if (lapType == LaplacianType.NormalizedSymmetric)
        {
            var invSqrtDegree = new double[n];
            for (int i = 0; i < n; i++)
                invSqrtDegree[i] = degree[i] > 1e-12 ? 1.0 / Math.Sqrt(degree[i]) : 0.0;

            for (int i = 0; i < n; i++)
            {
                laplacian[i, i] = 1.0;
                if (invSqrtDegree[i] == 0.0)
                    continue;

                int start = graph.RowPointers[i];
                int end = graph.RowPointers[i + 1];
                for (int edge = start; edge < end; edge++)
                {
                    int target = graph.Targets[edge];
                    double weight = graph.Weights[edge];
                    laplacian[i, target] -= invSqrtDegree[i] * weight * invSqrtDegree[target];
                }
            }

            return laplacian;
        }

        for (int i = 0; i < n; i++)
        {
            laplacian[i, i] = degree[i];
            int start = graph.RowPointers[i];
            int end = graph.RowPointers[i + 1];
            for (int edge = start; edge < end; edge++)
            {
                int target = graph.Targets[edge];
                double weight = graph.Weights[edge] > 0.0 ? graph.Weights[edge] : 1.0;
                laplacian[i, target] -= weight;
            }
        }

        return laplacian;
    }

    public static void BuildDenseColumnMajor(
        CsrGraph graph,
        LaplacianType lapType,
        Span<double> destination)
    {
        int n = graph.NodeCount;
        if (destination.Length != n * n)
            throw new ArgumentException("Destination dimensions must match N x N.", nameof(destination));

        ArrayPool<double> pool = ArrayPool<double>.Shared;
        double[] degreeArray = pool.Rent(n);
        Span<double> degree = degreeArray.AsSpan(0, n);
        degree.Clear();
        destination.Clear();

        try
        {
            for (int i = 0; i < n; i++)
            {
                int start = graph.RowPointers[i];
                int end = graph.RowPointers[i + 1];
                double sum = 0.0;
                for (int edge = start; edge < end; edge++)
                    sum += graph.Weights[edge] > 0.0 ? graph.Weights[edge] : 1.0;

                degree[i] = sum;
            }

            if (lapType == LaplacianType.NormalizedSymmetric)
            {
                double[] invSqrtArray = pool.Rent(n);
                Span<double> invSqrtDegree = invSqrtArray.AsSpan(0, n);

                try
                {
                    for (int i = 0; i < n; i++)
                        invSqrtDegree[i] = degree[i] > 1e-12 ? 1.0 / Math.Sqrt(degree[i]) : 0.0;

                    for (int i = 0; i < n; i++)
                    {
                        destination[i * n + i] = 1.0;
                        if (invSqrtDegree[i] == 0.0)
                            continue;

                        int start = graph.RowPointers[i];
                        int end = graph.RowPointers[i + 1];
                        for (int edge = start; edge < end; edge++)
                        {
                            int target = graph.Targets[edge];
                            double weight = graph.Weights[edge];
                            destination[target * n + i] -= invSqrtDegree[i] * weight * invSqrtDegree[target];
                        }
                    }
                }
                finally
                {
                    pool.Return(invSqrtArray);
                }

                return;
            }

            for (int i = 0; i < n; i++)
            {
                destination[i * n + i] = degree[i];
                int start = graph.RowPointers[i];
                int end = graph.RowPointers[i + 1];
                for (int edge = start; edge < end; edge++)
                {
                    int target = graph.Targets[edge];
                    double weight = graph.Weights[edge] > 0.0 ? graph.Weights[edge] : 1.0;
                    destination[target * n + i] -= weight;
                }
            }
        }
        finally
        {
            pool.Return(degreeArray);
        }
    }
}
