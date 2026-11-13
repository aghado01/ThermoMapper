using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// The principal angles between two subspaces, given orthonormal bases A (p × k_A) and B (p × k_B): min(k_A, k_B)
    /// angles, ascending. The cosines are the singular values of AᵀB and the sines those of the residual (I − LLᵀ)S of
    /// the basis S with fewer columns off the span of the other, L; the i-th smallest angle pairs the i-th largest cosine
    /// with the i-th smallest sine through atan2 (Björck and Golub 1973; Knyazev and Argentati 2002). The cosines hold a
    /// large angle and the sines a small one, so every angle keeps its accuracy, down to about √(p·k)·ε radians, where an
    /// arccos of the cosines alone resolves only about √(2ε) ≈ 2e-8 near zero. Both SVDs are <see cref="ThinSvd"/>'s.
    /// The bases are not checked for orthonormality.
    /// </summary>
    public static class PrincipalAngles
    {
        /// <summary>The angles between the column spans of <paramref name="a"/> and <paramref name="b"/>.</summary>
        public static double[] Between(double[,] a, double[,] b)
        {
            ArgumentNullException.ThrowIfNull(a);
            ArgumentNullException.ThrowIfNull(b);
            int p = a.GetLength(0), ka = a.GetLength(1), kb = b.GetLength(1);
            if (b.GetLength(0) != p) throw new ArgumentException("The bases have the same number of rows.", nameof(b));
            var flatA = new double[p * ka];
            var flatB = new double[p * kb];
            for (int j = 0; j < ka; j++) for (int i = 0; i < p; i++) flatA[i + j * p] = a[i, j];
            for (int j = 0; j < kb; j++) for (int i = 0; i < p; i++) flatB[i + j * p] = b[i, j];
            return Between(flatA, ka, flatB, kb, p);
        }

        /// <summary>
        /// The angles between the column spans of the column-major bases <paramref name="a"/> (p × <paramref name="ka"/>)
        /// and <paramref name="b"/> (p × <paramref name="kb"/>).
        /// </summary>
        public static double[] Between(ReadOnlySpan<double> a, int ka, ReadOnlySpan<double> b, int kb, int p)
        {
            if (p < 1 || ka < 1 || kb < 1) throw new ArgumentOutOfRangeException(nameof(p), "The bases need rows and columns.");
            if (a.Length < p * ka) throw new ArgumentException("The first basis is shorter than p × ka.", nameof(a));
            if (b.Length < p * kb) throw new ArgumentException("The second basis is shorter than p × kb.", nameof(b));

            // S has the fewer columns, L the more; the angles are those of span(S) against span(L).
            bool aSmall = ka <= kb;
            ReadOnlySpan<double> s = aSmall ? a : b, l = aSmall ? b : a;
            int k = Math.Min(ka, kb), kl = Math.Max(ka, kb);

            // C = Lᵀ S (kl × k), whose singular values are the cosines.
            var c = new double[kl, k];
            for (int j = 0; j < k; j++)
                for (int i = 0; i < kl; i++)
                {
                    double dot = 0.0;
                    for (int r = 0; r < p; r++) dot += l[r + i * p] * s[r + j * p];
                    c[i, j] = dot;
                }

            // R = S − L C = (I − LLᵀ) S (p × k), whose singular values are the sines.
            var residual = new double[p, k];
            for (int j = 0; j < k; j++)
                for (int r = 0; r < p; r++)
                {
                    double value = s[r + j * p];
                    for (int i = 0; i < kl; i++) value -= l[r + i * p] * c[i, j];
                    residual[r, j] = value;
                }

            double[] cosines = ThinSvd.Decompose(c, leftVectors: false).SingularValues;      // descending: the smallest angle first
            double[] sines = ThinSvd.Decompose(residual, leftVectors: false).SingularValues; // descending: the largest angle first
            var angles = new double[k];
            for (int i = 0; i < k; i++)
                angles[i] = Math.Atan2(Math.Min(sines[k - 1 - i], 1.0), Math.Min(cosines[i], 1.0));
            Array.Sort(angles);
            return angles;
        }
    }
}
