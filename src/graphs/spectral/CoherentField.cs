// src/graphs/spectral/CoherentField.cs
#nullable enable
using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using Graphs.Primitives;
using Maths.LinAlg;

namespace Graphs.Spectral;

/// <summary>
/// Builds viz-engine-ready coherent vector fields from low-frequency graph
/// Laplacian eigenmodes.
///
/// Output is a single flat row-major <c>double[]</c> of length <c>N × fieldDim</c>,
/// laid out as <c>[node 0 dim 0, node 0 dim 1, …, node N-1 dim D-1]</c> — the
/// exact format Three.js / WebGL expects for packed per-vertex attribute
/// buffers. Avoids any extra copy on the upload path.
///
/// Algorithm: pull bottom-K eigenpairs via <see cref="Spectral.ComputeBottomK"/>,
/// drop the trivial near-zero (connected-component) modes, then take the next
/// <c>fieldDim</c> non-trivial modes as the field components and L2-normalise
/// per node in place.
/// </summary>
public static class CoherentField
{
    public static double[] Build(
        CsrGraph graph,
        int fieldDim = 2,
        int basisK = 8,
        LaplacianType lapType = LaplacianType.Combinatorial,
        SolverKind solverKind = SolverKind.Auto,
        int seed = 0,
        DenseEigenOptions denseOptions = default,
        DenseLaplacianMaterialization denseMaterialization = DenseLaplacianMaterialization.Rectangular)
    {
        int n = graph.NodeCount;
        if (n == 0) return Array.Empty<double>();
        if (basisK < fieldDim) basisK = fieldDim + 1; // Safeguard spectral search space

        var pairs = Spectral.ComputeBottomK(
            graph,
            seed: seed,
            k: basisK,
            lapType: lapType,
            solverKind: solverKind,
            denseOptions: denseOptions,
            denseMaterialization: denseMaterialization);

        // Drop near-zero constant trivial modes (e.g., the algebraic multiplicity of connected components)
        int validModeStartIndex = 0;
        while (validModeStartIndex < pairs.Count && pairs[validModeStartIndex].Lambda < 1e-7)
        {
            validModeStartIndex++;
        }

        // Allocate a single contiguous flat array for the node vector fields
        // Stored in row-major layout: [Node 0 Dim 0, Node 0 Dim 1, ... Node N Dim M]
        // This mirrors exactly how Three.js/WebGL expects vertex attributes packed in a buffer.
        double[] coherentField = new double[n * fieldDim];

        // Project the modes into the field structure
        for (int d = 0; d < fieldDim; d++)
        {
            int pairIndex = validModeStartIndex + d;

            // If the graph doesn't have enough non-trivial modes, break out early
            if (pairIndex >= pairs.Count) break;

            ReadOnlySpan<double> modeVector = pairs[pairIndex].Vector;

            for (int i = 0; i < n; i++)
            {
                coherentField[i * fieldDim + d] = modeVector[i];
            }
        }

        // In-place vector normalization across the flat block using raw references
        NormalizeFieldInPlace(coherentField, n, fieldDim);

        return coherentField;
    }

    /// <summary>
    /// Highly optimized vector row field normalizer using tiered JIT loop vectorization.
    /// </summary>
    private static void NormalizeFieldInPlace(double[] field, int nodeCount, int fieldDim)
    {
        ref double fieldBase = ref MemoryMarshal.GetReference(field.AsSpan());

        for (int i = 0; i < nodeCount; i++)
        {
            int rowOffset = i * fieldDim;
            double sumSq = 0.0;

            // Loop unrolling for common spatial visualization tracking layouts (2D/3D)
            if (fieldDim == 2)
            {
                double v0 = Unsafe.Add(ref fieldBase, rowOffset);
                double v1 = Unsafe.Add(ref fieldBase, rowOffset + 1);
                sumSq = (v0 * v0) + (v1 * v1);
            }
            else if (fieldDim == 3)
            {
                double v0 = Unsafe.Add(ref fieldBase, rowOffset);
                double v1 = Unsafe.Add(ref fieldBase, rowOffset + 1);
                double v2 = Unsafe.Add(ref fieldBase, rowOffset + 2);
                sumSq = (v0 * v0) + (v1 * v1) + (v2 * v2);
            }
            else
            {
                for (int d = 0; d < fieldDim; d++)
                {
                    double v = Unsafe.Add(ref fieldBase, rowOffset + d);
                    sumSq += v * v;
                }
            }

            if (sumSq > 1e-12)
            {
                double invNorm = 1.0 / Math.Sqrt(sumSq);

                int d = 0;
                if (Vector256.IsHardwareAccelerated && fieldDim >= 4)
                {
                    var vInv = Vector256.Create(invNorm);
                    int simdLen = fieldDim & ~3;
                    for (; d < simdLen; d += 4)
                    {
                        var vVals = Vector256.LoadUnsafe(ref fieldBase, (uint)(rowOffset + d));
                        Vector256.Multiply(vVals, vInv).StoreUnsafe(ref fieldBase, (uint)(rowOffset + d));
                    }
                }
                for (; d < fieldDim; d++)
                {
                    Unsafe.Add(ref fieldBase, rowOffset + d) *= invNorm;
                }
            }
        }
    }
}
