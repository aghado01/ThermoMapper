using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// Sums and dot products by pairwise summation (Higham 1993): a balanced tree over blocks of <see cref="Block"/>
    /// terms, each block summed in four interleaved accumulators. The rounding error grows as
    /// ε·(log₂(N/<see cref="Block"/>) + <see cref="Block"/>/4) instead of left-to-right summation's ε·N, and it does
    /// not pile up when many equal small terms follow one large one, as they do in a near-constant design column
    /// (the BASS port's fixture F-O7). Accumulates the cross-products that become Gram entries. Plain double
    /// arithmetic in a fixed order, so every machine gets the same bits.
    /// </summary>
    public static class PairwiseSum
    {
        /// <summary>Terms summed directly, in four interleaved accumulators, at each leaf of the tree.</summary>
        public const int Block = 128;

        /// <summary>Σ x_i.</summary>
        public static double Of(ReadOnlySpan<double> x)
        {
            int n = x.Length;
            if (n <= Block)
            {
                double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
                int i = 0;
                for (; i + 3 < n; i += 4)
                {
                    s0 += x[i];
                    s1 += x[i + 1];
                    s2 += x[i + 2];
                    s3 += x[i + 3];
                }
                for (; i < n; i++) s0 += x[i];
                return (s0 + s1) + (s2 + s3);
            }
            int half = Split(n);
            return Of(x[..half]) + Of(x[half..]);
        }

        /// <summary>Σ a_i·b_i.</summary>
        public static double Dot(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
        {
            if (a.Length != b.Length) throw new ArgumentException("Vector dimensions must match.");
            int n = a.Length;
            if (n <= Block)
            {
                double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
                int i = 0;
                for (; i + 3 < n; i += 4)
                {
                    s0 += a[i] * b[i];
                    s1 += a[i + 1] * b[i + 1];
                    s2 += a[i + 2] * b[i + 2];
                    s3 += a[i + 3] * b[i + 3];
                }
                for (; i < n; i++) s0 += a[i] * b[i];
                return (s0 + s1) + (s2 + s3);
            }
            int half = Split(n);
            return Dot(a[..half], b[..half]) + Dot(a[half..], b[half..]);
        }

        // The left part: about half, on a block boundary, and at least one block, so every leaf but the last is full.
        private static int Split(int n) => Math.Max(Block, n / 2 / Block * Block);
    }
}
