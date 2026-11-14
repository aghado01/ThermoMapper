using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// A factor's condition-estimate scratch, kept across its factorizations: the three work vectors and E⁻¹'s apply,
    /// v ← D^(½)·A⁻¹·D^(½)·v, over the factor's solve and its √A_jj array, which each factorization refills in place.
    /// </summary>
    internal sealed class ConditionWorkspace
    {
        public ConditionWorkspace(ISpdFactor factor, double[] sqrtDiagonal)
        {
            int n = factor.Dimension;
            Dimension = n;
            X = new double[n];
            Sign = new double[n];
            Z = new double[n];
            ApplyInverse = v =>
            {
                for (int i = 0; i < n; i++) v[i] *= sqrtDiagonal[i];
                factor.SolveInPlace(v);
                for (int i = 0; i < n; i++) v[i] *= sqrtDiagonal[i];
            };
        }

        public int Dimension { get; }

        public double[] X { get; }

        public double[] Sign { get; }

        public double[] Z { get; }

        public Action<double[]> ApplyInverse { get; }
    }

    /// <summary>
    /// Condition estimates through a factor's solves. The equilibrated condition number of a symmetric
    /// positive-definite A is κ₁(E) = ‖E‖₁·‖E⁻¹‖₁, with E = D^(−½)·A·D^(−½) and D = diag(A): the number that governs
    /// Cholesky's accuracy (van der Sluis 1969; Demmel 1989). ‖E‖₁ is read from A when it is factored; ‖E⁻¹‖₁ is
    /// estimated by Higham's (1988) refinement of Hager's (1984) method, as LAPACK's xLACN2 does for xPOCON, in at most
    /// eleven solves. The estimate is ‖E⁻¹v‖₁/‖v‖₁ for vectors it found, so it is a lower bound on the exact norm; in
    /// practice it is rarely below a third of it, but nothing guarantees that.
    /// </summary>
    internal static class ConditionEstimate
    {
        private const int MaxIterations = 5;

        /// <summary>κ₁ of the equilibrated matrix, estimated through <paramref name="factor"/>'s solves, given the
        /// square roots of A's diagonal and ‖E‖₁.</summary>
        public static double Equilibrated(ISpdFactor factor, double[] sqrtDiagonal, double equilibratedNorm1)
        {
            int n = factor.Dimension;
            void ApplyInverse(double[] v)   // v ← E⁻¹v = D^(½)·A⁻¹·D^(½)·v
            {
                for (int i = 0; i < n; i++) v[i] *= sqrtDiagonal[i];
                factor.SolveInPlace(v);
                for (int i = 0; i < n; i++) v[i] *= sqrtDiagonal[i];
            }
            return equilibratedNorm1 * SymmetricNorm1(n, ApplyInverse);
        }

        /// <summary>κ₁ of the equilibrated matrix, as <see cref="Equilibrated(ISpdFactor, double[], double)"/>, through a
        /// factor's own <paramref name="workspace"/>: the same bits, with nothing allocated per call.</summary>
        public static double Equilibrated(ConditionWorkspace workspace, double equilibratedNorm1) =>
            equilibratedNorm1 * SymmetricNorm1(workspace.Dimension, workspace.ApplyInverse, workspace.X, workspace.Sign, workspace.Z);

        /// <summary>A lower-bound estimate of ‖B‖₁ for a symmetric B of order <paramref name="n"/>, given
        /// <paramref name="apply"/>, which overwrites a vector v with Bv.</summary>
        public static double SymmetricNorm1(int n, Action<double[]> apply) =>
            SymmetricNorm1(n, apply, new double[n], new double[n], new double[n]);

        /// <summary>As <see cref="SymmetricNorm1(int, Action{double[]})"/>, in the caller's work vectors of length
        /// <paramref name="n"/>, whose contents on entry are never read.</summary>
        public static double SymmetricNorm1(int n, Action<double[]> apply, double[] x, double[] sign, double[] z)
        {
            if (n == 1)
            {
                x[0] = 1.0;
                apply(x);
                return Math.Abs(x[0]);
            }

            Array.Fill(x, 1.0 / n);
            apply(x);
            double estimate = Norm1(x);
            for (int i = 0; i < n; i++) sign[i] = Sign(x[i]);
            Array.Copy(sign, z, n);
            apply(z);                                  // Bᵀξ, which is Bξ for symmetric B
            int j = ArgMaxAbs(z);

            for (int iteration = 2; iteration <= MaxIterations; iteration++)
            {
                Array.Clear(x);
                x[j] = 1.0;
                apply(x);
                double previous = estimate;
                estimate = Norm1(x);
                bool repeated = true;
                for (int i = 0; i < n && repeated; i++) repeated = Sign(x[i]) == sign[i];
                if (repeated || estimate <= previous)  // converged, or cycling
                {
                    estimate = Math.Max(estimate, previous);
                    break;
                }
                for (int i = 0; i < n; i++) sign[i] = Sign(x[i]);
                Array.Copy(sign, z, n);
                apply(z);
                int last = j;
                j = ArgMaxAbs(z);
                if (Math.Abs(z[last]) == Math.Abs(z[j])) break;
            }

            // Higham's extra vector, x_i = (−1)^i·(1 + i/(n−1)), for the matrices the iteration misjudges.
            for (int i = 0; i < n; i++) x[i] = (i % 2 == 0 ? 1.0 : -1.0) * (1.0 + (double)i / (n - 1));
            apply(x);
            return Math.Max(estimate, 2.0 * Norm1(x) / (3.0 * n));
        }

        /// <summary>‖E‖₁ = max_j Σ_i |A_ij|/√(A_ii·A_jj) over a symmetric matrix's lower band (<paramref name="band"/>
        /// [d, j] = A(j+d, j), d = 0..p), writing √A_jj into <paramref name="sqrtDiagonal"/>; NaN unless every diagonal
        /// entry is finite and positive.</summary>
        public static double EquilibratedNorm1(double[,] band, int n, int p, double[] sqrtDiagonal, double[] columnSums)
        {
            for (int j = 0; j < n; j++)
            {
                double d = band[0, j];
                if (!(d > 0.0 && double.IsFinite(d))) return double.NaN;
                sqrtDiagonal[j] = Math.Sqrt(d);
                columnSums[j] = 1.0;
            }
            for (int j = 0; j < n; j++)
            {
                int dmax = Math.Min(p, n - 1 - j);
                for (int d = 1; d <= dmax; d++)
                {
                    double e = Math.Abs(band[d, j]) / (sqrtDiagonal[j] * sqrtDiagonal[j + d]);
                    columnSums[j] += e;
                    columnSums[j + d] += e;
                }
            }
            double norm = 0.0;
            for (int j = 0; j < n; j++) norm = Math.Max(norm, columnSums[j]);
            return norm;
        }

        private static double Sign(double v) => v >= 0.0 ? 1.0 : -1.0;

        private static double Norm1(double[] v)
        {
            double s = 0.0;
            foreach (double t in v) s += Math.Abs(t);
            return s;
        }

        private static int ArgMaxAbs(double[] v)
        {
            int best = 0;
            for (int i = 1; i < v.Length; i++) if (Math.Abs(v[i]) > Math.Abs(v[best])) best = i;
            return best;
        }
    }
}
