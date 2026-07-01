#nullable enable
using System;
using System.Collections.Generic;

namespace Maths.Topology;

/// <summary>
/// A simple, unoptimized Z/2 linear algebra utility for the naive zigzag oracle.
/// Uses dense bool arrays for obvious correctness over performance.
/// </summary>
public static class Z2LinearAlgebra
{
    /// <summary>
    /// Computes the rank of a dense Z/2 matrix.
    /// Matrix is given as rows x cols, with M[r, c].
    /// </summary>
    public static int Rank(bool[,] matrix)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);
        bool[,] m = (bool[,])matrix.Clone();

        int rank = 0;
        int r = 0;
        for (int c = 0; c < cols && r < rows; c++)
        {
            // Find pivot
            int pivot = -1;
            for (int i = r; i < rows; i++)
            {
                if (m[i, c])
                {
                    pivot = i;
                    break;
                }
            }

            if (pivot == -1) continue;

            // Swap rows
            if (pivot != r)
            {
                for (int j = c; j < cols; j++)
                {
                    (m[r, j], m[pivot, j]) = (m[pivot, j], m[r, j]);
                }
            }

            // Eliminate below
            for (int i = r + 1; i < rows; i++)
            {
                if (m[i, c])
                {
                    for (int j = c; j < cols; j++)
                    {
                        m[i, j] ^= m[r, j];
                    }
                }
            }

            rank++;
            r++;
        }

        return rank;
    }

    /// <summary>
    /// Finds a basis for the nullspace (kernel) of a Z/2 matrix.
    /// Returns a list of basis vectors (each of length = cols).
    /// </summary>
    public static List<bool[]> Nullspace(bool[,] matrix)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);
        bool[,] m = (bool[,])matrix.Clone();

        // Track pivot columns to identify free variables
        int[] pivotCols = new int[rows];
        for (int i = 0; i < rows; i++) pivotCols[i] = -1;

        int r = 0;
        for (int c = 0; c < cols && r < rows; c++)
        {
            int pivot = -1;
            for (int i = r; i < rows; i++)
            {
                if (m[i, c])
                {
                    pivot = i;
                    break;
                }
            }

            if (pivot == -1) continue;

            if (pivot != r)
            {
                for (int j = c; j < cols; j++)
                {
                    (m[r, j], m[pivot, j]) = (m[pivot, j], m[r, j]);
                }
            }

            pivotCols[r] = c;

            // Eliminate BOTH below and above for RREF
            for (int i = 0; i < rows; i++)
            {
                if (i != r && m[i, c])
                {
                    for (int j = c; j < cols; j++)
                    {
                        m[i, j] ^= m[r, j];
                    }
                }
            }

            r++;
        }

        bool[] isFree = new bool[cols];
        for (int j = 0; j < cols; j++) isFree[j] = true;
        for (int i = 0; i < r; i++)
        {
            if (pivotCols[i] != -1) isFree[pivotCols[i]] = false;
        }

        var basis = new List<bool[]>();
        for (int j = 0; j < cols; j++)
        {
            if (isFree[j])
            {
                bool[] vec = new bool[cols];
                vec[j] = true;
                for (int i = 0; i < r; i++)
                {
                    if (pivotCols[i] != -1 && m[i, j])
                    {
                        vec[pivotCols[i]] = true;
                    }
                }
                basis.Add(vec);
            }
        }

        return basis;
    }

    /// <summary>
    /// Reduces a set of vectors modulo the span of another set of vectors.
    /// Vectors are represented as columns. 
    /// Returns a list of the reduced vectors. If a reduced vector becomes zero, it remains in the list as zero.
    /// </summary>
    public static List<bool[]> ReduceModuloSpan(List<bool[]> vectorsToReduce, List<bool[]> spanBasis)
    {
        if (vectorsToReduce.Count == 0) return new List<bool[]>();
        int dim = vectorsToReduce[0].Length;

        // Build a matrix with spanBasis first, then vectorsToReduce
        // We will row reduce to eliminate the spanBasis components from vectorsToReduce.
        // Actually, it's easier to put spanBasis as rows of a matrix, row reduce it,
        // then for each vector, eliminate using the pivot rows.

        bool[,] basisMat = new bool[spanBasis.Count, dim];
        for (int i = 0; i < spanBasis.Count; i++)
        {
            for (int j = 0; j < dim; j++) basisMat[i, j] = spanBasis[i][j];
        }

        int rows = spanBasis.Count;
        int r = 0;
        for (int c = 0; c < dim && r < rows; c++)
        {
            int pivot = -1;
            for (int i = r; i < rows; i++)
            {
                if (basisMat[i, c]) { pivot = i; break; }
            }
            if (pivot == -1) continue;

            if (pivot != r)
            {
                for (int j = c; j < dim; j++)
                {
                    (basisMat[r, j], basisMat[pivot, j]) = (basisMat[pivot, j], basisMat[r, j]);
                }
            }

            for (int i = r + 1; i < rows; i++)
            {
                if (basisMat[i, c])
                {
                    for (int j = c; j < dim; j++) basisMat[i, j] ^= basisMat[r, j];
                }
            }
            r++;
        }

        var result = new List<bool[]>();
        foreach (var vec in vectorsToReduce)
        {
            bool[] v = (bool[])vec.Clone();
            for (int i = 0; i < r; i++)
            {
                // Find leading 1 of basis row
                int leading = -1;
                for (int j = 0; j < dim; j++)
                {
                    if (basisMat[i, j]) { leading = j; break; }
                }

                if (leading != -1 && v[leading])
                {
                    for (int j = leading; j < dim; j++) v[j] ^= basisMat[i, j];
                }
            }
            result.Add(v);
        }

        return result;
    }
}
