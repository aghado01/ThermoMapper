using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// Σ = D + U·diag(s)·Uᵀ, a positive diagonal plus a low-rank term of rank r ≪ m: BASS calibration's functional
    /// covariance, noise and truncation variance on the diagonal and the emulator's component variances on the basis.
    /// Solves, bᵀΣ⁻¹b, ln|Σ| and draws cost O(m·r) after an O(m·r²) setup, and no m × m matrix is ever formed.
    /// <para>
    /// The factor of I + WᵀW, W = D^(−½)·U·diag(√s), is taken as the thin SVD W = Q·diag(λ)·Pᵀ (<see cref="ThinSvd"/>),
    /// so that with β = D^(−½)b, bᵀΣ⁻¹b = ‖β − QQᵀβ‖² + Σ_k (q_kᵀβ)²/(1 + λ_k²) and ln|Σ| = Σ ln d_i + Σ ln(1 + λ_k²):
    /// sums of nonnegative terms, where the Woodbury form bᵀD⁻¹b − zᵀ(I + WᵀW)⁻¹z cancels when b lies along a heavy
    /// direction.
    /// </para>
    /// Immutable once constructed, so concurrent chains may share it.
    /// </summary>
    public sealed class DiagonalPlusLowRank
    {
        private readonly double[] _sqrtDiagonal;
        private readonly double[][] _q;          // left singular vectors of W, each of length m
        private readonly double[] _lambda2;      // λ_k²

        // The functions of WWᵀ applied through its eigenvectors: g(λ²) − 1 for g = 1/(1 + λ²), √(1 + λ²), 1/√(1 + λ²).
        private enum Function { Inverse, Root, InverseRoot }

        /// <param name="diagonal">D's entries, each finite and positive.</param>
        /// <param name="factor">U, m × r.</param>
        /// <param name="weights">s, r finite nonnegative weights.</param>
        public DiagonalPlusLowRank(ReadOnlySpan<double> diagonal, double[,] factor, ReadOnlySpan<double> weights)
        {
            ArgumentNullException.ThrowIfNull(factor);
            int m = diagonal.Length, r = factor.GetLength(1);
            if (factor.GetLength(0) != m) throw new ArgumentException($"The factor needs {m} rows, one per diagonal entry.", nameof(factor));
            if (weights.Length != r) throw new ArgumentException($"Expected {r} weights, one per factor column.", nameof(weights));

            _sqrtDiagonal = new double[m];
            double logDet = 0.0;
            for (int i = 0; i < m; i++)
            {
                double d = diagonal[i];
                if (!(d > 0.0 && double.IsFinite(d))) throw new ArgumentException($"Diagonal entry {i} is {d}; each must be finite and positive.", nameof(diagonal));
                _sqrtDiagonal[i] = Math.Sqrt(d);
                logDet += Math.Log(d);
            }
            var w = new double[m, r];
            for (int k = 0; k < r; k++)
            {
                double sk = weights[k];
                if (!(sk >= 0.0 && double.IsFinite(sk))) throw new ArgumentException($"Weight {k} is {sk}; each must be finite and nonnegative.", nameof(weights));
                double root = Math.Sqrt(sk);
                for (int i = 0; i < m; i++) w[i, k] = factor[i, k] * root / _sqrtDiagonal[i];
            }

            ThinSvdResult svd = ThinSvd.Decompose(w);
            int t = svd.SingularValues.Length;
            _q = svd.LeftVectors;
            _lambda2 = new double[t];
            for (int k = 0; k < t; k++)
            {
                double l = svd.SingularValues[k];
                _lambda2[k] = l * l;
                logDet += double.LogP1(_lambda2[k]);
            }
            LogDet = logDet;
        }

        /// <summary>m, the order of Σ.</summary>
        public int Dimension => _sqrtDiagonal.Length;

        /// <summary>ln|Σ|.</summary>
        public double LogDet { get; }

        /// <summary>bᵀΣ⁻¹b, as a sum of nonnegative terms.</summary>
        public double InverseQuadraticForm(ReadOnlySpan<double> b)
        {
            int m = RequireLength(b.Length);
            var beta = new double[m];
            for (int i = 0; i < m; i++) beta[i] = b[i] / _sqrtDiagonal[i];
            var c = new double[_q.Length];
            Project(beta, c);                                 // β ← β − QQᵀβ, c = Qᵀβ
            double q = PairwiseSum.Dot(beta, beta);
            for (int k = 0; k < c.Length; k++) q += c[k] * c[k] / (1.0 + _lambda2[k]);
            return q;
        }

        /// <summary>Solves <c>Σ·x = b</c>.</summary>
        public double[] Solve(ReadOnlySpan<double> b)
        {
            double[] x = b.ToArray();
            SolveInPlace(x);
            return x;
        }

        /// <summary>In place: <paramref name="x"/> enters as b and leaves as Σ⁻¹b = D^(−½)(I + WWᵀ)⁻¹D^(−½)b.</summary>
        public void SolveInPlace(Span<double> x)
        {
            int m = RequireLength(x.Length);
            for (int i = 0; i < m; i++) x[i] /= _sqrtDiagonal[i];
            Apply(x, Function.Inverse);
            for (int i = 0; i < m; i++) x[i] /= _sqrtDiagonal[i];
        }

        /// <summary>Maps standard normals <paramref name="z"/> to a draw from N(0, Σ): D^(½)(I + WWᵀ)^(½)z.</summary>
        public void Sample(ReadOnlySpan<double> z, Span<double> result)
        {
            int m = RequireLength(z.Length);
            if (result.Length != m) throw new ArgumentException("The result length must equal the order.", nameof(result));
            z.CopyTo(result);
            Apply(result, Function.Root);
            for (int i = 0; i < m; i++) result[i] *= _sqrtDiagonal[i];
        }

        /// <summary>Maps standard normals <paramref name="z"/> to a draw from N(0, Σ⁻¹): D^(−½)(I + WWᵀ)^(−½)z.</summary>
        public void SampleInnovation(ReadOnlySpan<double> z, Span<double> result)
        {
            int m = RequireLength(z.Length);
            if (result.Length != m) throw new ArgumentException("The result length must equal the order.", nameof(result));
            z.CopyTo(result);
            Apply(result, Function.InverseRoot);
            for (int i = 0; i < m; i++) result[i] /= _sqrtDiagonal[i];
        }

        // x ← g(I + WWᵀ)·x = x + Q·diag(g(1 + λ²) − 1)·Qᵀx. Each g − 1 is written without cancellation for small λ.
        private void Apply(Span<double> x, Function g)
        {
            int t = _q.Length;
            Span<double> c = t <= 64 ? stackalloc double[t] : new double[t];
            for (int k = 0; k < t; k++) c[k] = PairwiseSum.Dot(_q[k], x);
            for (int k = 0; k < t; k++)
            {
                double l2 = _lambda2[k], root = Math.Sqrt(1.0 + l2);
                double f = g switch
                {
                    Function.Inverse => -l2 / (1.0 + l2),
                    Function.Root => l2 / (root + 1.0),
                    _ => -l2 / (root * (root + 1.0)),
                };
                double step = f * c[k];
                if (step != 0.0) MatrixOps.ScaledSubtract(x, _q[k], -step, x.Length);
            }
        }

        // β ← β − QQᵀβ and c = Qᵀβ; a second pass removes what the first left.
        private void Project(double[] beta, double[] c)
        {
            for (int pass = 0; pass < 2; pass++)
                for (int k = 0; k < _q.Length; k++)
                {
                    double ck = PairwiseSum.Dot(_q[k], beta);
                    c[k] += ck;
                    MatrixOps.ScaledSubtract(beta, _q[k], ck, beta.Length);
                }
        }

        private int RequireLength(int length)
        {
            if (length != _sqrtDiagonal.Length) throw new ArgumentException($"Expected length {_sqrtDiagonal.Length}, got {length}.");
            return length;
        }
    }
}
