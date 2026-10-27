using System;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace Maths.LinAlg
{
    /// <summary>
    /// Sums and dot products by pairwise summation (Higham 1993): a balanced tree over blocks of <see cref="Block"/>
    /// terms, each block summed in four interleaved accumulators. The rounding error grows as
    /// ε·(log₂(N/<see cref="Block"/>) + <see cref="Block"/>/4) instead of left-to-right summation's ε·N, and it does
    /// not pile up when many equal small terms follow one large one, as they do in a near-constant design column
    /// (the BASS port's fixture F-O7). Accumulates the cross-products that become Gram entries. The four accumulators
    /// are the lanes of one 256-bit vector, whose adds and multiplies are exact per lane in IEEE arithmetic, so every
    /// machine, accelerated or not, gets the same bits.
    /// </summary>
    public static class PairwiseSum
    {
        /// <summary>Terms summed directly, in four interleaved accumulators, at each leaf of the tree.</summary>
        public const int Block = 128;

        /// <summary>Σ x_i.</summary>
        public static double Of(ReadOnlySpan<double> x)
        {
            int n = x.Length;
            if (n > Block)
            {
                int half = Split(n);
                return Of(x[..half]) + Of(x[half..]);
            }
            var lanes = Vector256<double>.Zero;
            ref double xr = ref MemoryMarshal.GetReference(x);
            int i = 0;
            for (; i + 3 < n; i += 4) lanes = Vector256.Add(lanes, Vector256.LoadUnsafe(ref xr, (nuint)i));
            double s0 = lanes.GetElement(0);
            for (; i < n; i++) s0 += x[i];
            return (s0 + lanes.GetElement(1)) + (lanes.GetElement(2) + lanes.GetElement(3));
        }

        /// <summary>Σ a_i·b_i.</summary>
        public static double Dot(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        {
            if (a.Length != b.Length) throw new ArgumentException("Vector dimensions must match.");
            int n = a.Length;
            if (n > Block)
            {
                int half = Split(n);
                return Dot(a[..half], b[..half]) + Dot(a[half..], b[half..]);
            }
            var lanes = Vector256<double>.Zero;
            ref double ar = ref MemoryMarshal.GetReference(a);
            ref double br = ref MemoryMarshal.GetReference(b);
            int i = 0;
            for (; i + 3 < n; i += 4)
                lanes = Vector256.Add(lanes, Vector256.Multiply(Vector256.LoadUnsafe(ref ar, (nuint)i), Vector256.LoadUnsafe(ref br, (nuint)i)));
            double s0 = lanes.GetElement(0);
            for (; i < n; i++) s0 += a[i] * b[i];
            return (s0 + lanes.GetElement(1)) + (lanes.GetElement(2) + lanes.GetElement(3));
        }

        // The left part: about half, on a block boundary, and at least one block, so every leaf but the last is full.
        private static int Split(int n) => Math.Max(Block, n / 2 / Block * Block);
    }
}
