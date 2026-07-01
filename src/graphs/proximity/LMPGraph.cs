using System;
using System.Collections.Generic;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Threading.Tasks;
using Graphs.Primitives;

namespace Graphs.Proximity
{
    public static class LocalMutualProximity
    {
        /// <summary>
        /// Re-weights a CSR graph using the Local Scaling Mutual Proximity algorithm.
        /// Outputs a new CsrGraph sharing the original topology, but with edge weights
        /// replaced by MP probabilities [0.0, 1.0].
        /// </summary>
        /// <param name="graph">The sparse adjacency structure.</param>
        /// <param name="weightsAreCouplings">
        /// If false (raw distances): counts neighbors further away (&gt; weight).
        /// If true (J couplings): counts neighbors with weaker bonds (&lt; weight).
        /// </param>
        /// <param name="protectedEdges">
        /// Undirected vertex pairs that keep their input weight instead of the MP product
        /// (load-bearing H1 cycle edges from PH).
        /// </param>
        public static CsrGraph ApplyLocalScaling(
            CsrGraph graph,
            bool weightsAreCouplings = false,
            IReadOnlySet<(int Lo, int Hi)>? protectedEdges = null)
        {
            double[] mpWeights = new double[graph.Weights.Length];

            Parallel.For(0, graph.NodeCount, u =>
            {
                int rowStartU = graph.RowPointers[u];
                int degreeU = graph.RowPointers[u + 1] - rowStartU;
                if (degreeU == 0) return;

                // The local distribution for node U
                ReadOnlySpan<double> rowU = new ReadOnlySpan<double>(graph.Weights, rowStartU, degreeU);

                for (int e = rowStartU; e < rowStartU + degreeU; e++)
                {
                    int v = graph.Targets[e];
                    double w = graph.Weights[e];

                    // Fetch the local distribution for the target node V
                    int rowStartV = graph.RowPointers[v];
                    int degreeV = graph.RowPointers[v + 1] - rowStartV;
                    ReadOnlySpan<double> rowV = new ReadOnlySpan<double>(graph.Weights, rowStartV, degreeV);

                    int lo = u < v ? u : v;
                    int hi = u < v ? v : u;
                    if (protectedEdges?.Contains((lo, hi)) == true)
                    {
                        mpWeights[e] = w;
                        continue;
                    }

                    // P(U > w) and P(V > w)
                    int countU = weightsAreCouplings ? CountWeakerCouplings(rowU, w) : CountLargerDistances(rowU, w);
                    int countV = weightsAreCouplings ? CountWeakerCouplings(rowV, w) : CountLargerDistances(rowV, w);

                    double pU = (double)countU / degreeU;
                    double pV = (double)countV / degreeV;

                    // Final MP score is the intersection of marginal probabilities
                    mpWeights[e] = pU * pV;
                }
            });

            return new CsrGraph
            {
                Targets = graph.Targets,         // Share the exact same topology array
                RowPointers = graph.RowPointers, // Share the row pointers
                NodeCount = graph.NodeCount,
                Weights = mpWeights              // Inject the new probabilities
            };
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int CountLargerDistances(ReadOnlySpan<double> values, double threshold)
        {
            int count = 0;
            int i = 0;

            if (Vector256.IsHardwareAccelerated && values.Length >= 4)
            {
                Vector256<double> vecTarget = Vector256.Create(threshold);
                ref double ptr = ref MemoryMarshal.GetReference(values);

                for (; i <= values.Length - 4; i += 4)
                {
                    Vector256<double> vecValues = Vector256.LoadUnsafe(ref ptr, (nuint)i);
                    Vector256<double> mask = Vector256.GreaterThan(vecValues, vecTarget);

                    // Extract a 4-bit mask (1 bit per double lane) and pop-count it
                    uint bits = mask.ExtractMostSignificantBits();
                    count += BitOperations.PopCount(bits);
                }
            }

            // Tail loop for non-multiple of 4
            for (; i < values.Length; i++)
            {
                if (values[i] > threshold) count++;
            }

            return count;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static int CountWeakerCouplings(ReadOnlySpan<double> values, double threshold)
        {
            int count = 0;
            int i = 0;

            if (Vector256.IsHardwareAccelerated && values.Length >= 4)
            {
                Vector256<double> vecTarget = Vector256.Create(threshold);
                ref double ptr = ref MemoryMarshal.GetReference(values);

                for (; i <= values.Length - 4; i += 4)
                {
                    Vector256<double> vecValues = Vector256.LoadUnsafe(ref ptr, (nuint)i);
                    Vector256<double> mask = Vector256.LessThan(vecValues, vecTarget);

                    uint bits = mask.ExtractMostSignificantBits();
                    count += BitOperations.PopCount(bits);
                }
            }

            for (; i < values.Length; i++)
            {
                if (values[i] < threshold) count++;
            }

            return count;
        }
    }
}
