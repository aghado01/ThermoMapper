using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// In-place lower-triangular Cholesky decomposition Σ = LLᵀ with derived quantities
    /// for Gaussian statistics: ln|Σ|, Σ⁻¹, and Cholesky-factor sampling.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Allocate once per component with <see cref="CholeskyDecomposition(int)"/>, then call
    /// <see cref="Decompose"/> each time the covariance changes. This avoids per-iteration
    /// heap allocations in EM loops.
    /// </para>
    /// <para>
    /// <strong>Diagonal floor.</strong> If a pivot drops to ≤ 1e-12 during decomposition
    /// it is clamped, providing a soft guard against near-singular matrices without
    /// throwing. Callers that require strict positive-definiteness should add a
    /// regularisation term to the diagonal before calling <see cref="Decompose"/>.
    /// </para>
    /// </remarks>
    public sealed class CholeskyDecomposition
    {
        private readonly int _dim;

        // Lower-triangular factor L  (Σ = L·Lᵀ).
        // Upper triangle is always zero; only [i, j] for j ≤ i is written.
        private readonly double[,] _l;

        // L⁻¹ — lower-triangular intermediate produced during Decompose.
        // Not meaningful to callers; used internally to build Σ⁻¹.
        private readonly double[,] _lInv;

        /// <summary>
        /// ln|Σ| = 2 · Σ_i ln(L_{ii}). Valid after <see cref="Decompose"/> has been called.
        /// </summary>
        public double LogDet { get; private set; }

        /// <param name="dim">Dimensionality of the square covariance matrix.</param>
        public CholeskyDecomposition(int dim)
        {
            _dim = dim;
            _l = new double[dim, dim];
            _lInv = new double[dim, dim];
        }

        // ── Decompose ────────────────────────────────────────────────────────────

        /// <summary>
        /// Decomposes <paramref name="sigma"/> into L·Lᵀ, updating <see cref="LogDet"/>
        /// and preparing the internal L⁻¹ required by <see cref="WriteInverseTo"/>
        /// and <see cref="Sample"/>.
        /// </summary>
        /// <param name="sigma">
        /// Square, symmetric, positive-semi-definite covariance matrix of size dim × dim.
        /// Only the lower triangle is read; the upper triangle is ignored.
        /// </param>
        public void Decompose(double[,] sigma)
        {
            double logDetL = 0.0;

            // ── Step 1: L via Cholesky-Banachiewicz ──────────────────────────────
            for (int i = 0; i < _dim; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    double sum = sigma[i, j];
                    for (int k = 0; k < j; k++)
                        sum -= _l[i, k] * _l[j, k];

                    if (i == j)
                    {
                        if (sum <= 1e-12) sum = 1e-12; // diagonal floor
                        _l[i, i] = Math.Sqrt(sum);
                        logDetL += Math.Log(_l[i, i]);
                    }
                    else
                    {
                        _l[i, j] = sum / _l[j, j];
                        _l[j, i] = 0.0; // keep upper triangle zero
                    }
                }
            }

            LogDet = 2.0 * logDetL; // ln|Σ| = 2·ln|L|

            // ── Step 2: L⁻¹ via forward substitution ────────────────────────────
            for (int i = 0; i < _dim; i++)
            {
                _lInv[i, i] = 1.0 / _l[i, i];
                for (int j = 0; j < i; j++)
                {
                    double sum = 0.0;
                    for (int k = j; k < i; k++)
                        sum -= _l[i, k] * _lInv[k, j];
                    _lInv[i, j] = sum / _l[i, i];
                    _lInv[j, i] = 0.0;
                }
            }
        }

        // ── Σ⁻¹ ─────────────────────────────────────────────────────────────────

        /// <summary>
        /// Writes Σ⁻¹ = (L⁻¹)ᵀ · L⁻¹ into <paramref name="target"/>.
        /// Must be called after <see cref="Decompose"/>.
        /// </summary>
        /// <param name="target">Pre-allocated dim × dim matrix to receive Σ⁻¹.</param>
        public void WriteInverseTo(double[,] target)
        {
            // L⁻¹ is lower-triangular: _lInv[k, i] = 0 for k < i.
            // Sum therefore starts at k = j (the larger of i and j).
            for (int i = 0; i < _dim; i++)
            {
                for (int j = i; j < _dim; j++)
                {
                    double sum = 0.0;
                    for (int k = j; k < _dim; k++)
                        sum += _lInv[k, i] * _lInv[k, j];
                    target[i, j] = sum;
                    target[j, i] = sum; // Σ⁻¹ is symmetric
                }
            }
        }

        // ── Solve ─────────────────────────────────────────────────────────────────

        /// <summary>
        /// Solves Σx = b for x = Σ⁻¹b by forward/back substitution on the factor (Σ = L·Lᵀ ⇒ L·y = b, then
        /// Lᵀ·x = y), in O(n²) — the triangular-solve path, which avoids building the full inverse via
        /// <see cref="WriteInverseTo"/> when only a solve is needed. Mirrors <c>BandCholesky.Solve</c>. Must be
        /// called after <see cref="Decompose"/>.
        /// </summary>
        /// <param name="b">Right-hand side of length dim.</param>
        /// <returns>The solution x of length dim.</returns>
        public double[] Solve(double[] b)
        {
            ArgumentNullException.ThrowIfNull(b);
            if (b.Length != _dim) throw new ArgumentException("Right-hand side length must equal the matrix dimension.", nameof(b));
            var x = new double[_dim];

            // Forward substitution: L·y = b (y accumulated in x; L is lower-triangular).
            for (int i = 0; i < _dim; i++)
            {
                double sum = b[i];
                for (int k = 0; k < i; k++) sum -= _l[i, k] * x[k];
                x[i] = sum / _l[i, i];
            }
            // Back substitution: Lᵀ·x = y  (Lᵀ[i,k] = L[k,i]).
            for (int i = _dim - 1; i >= 0; i--)
            {
                double sum = x[i];
                for (int k = i + 1; k < _dim; k++) sum -= _l[k, i] * x[k];
                x[i] = sum / _l[i, i];
            }
            return x;
        }

        // ── L read-out ───────────────────────────────────────────────────────────

        /// <summary>
        /// Writes the lower-triangular Cholesky factor L row-major into <paramref name="dst"/>.
        /// Length of <paramref name="dst"/> must be dim×dim. Upper-triangle entries are written as 0.
        /// Must be called after <see cref="Decompose"/>.
        /// </summary>
        public void WriteLTo(Span<double> dst)
        {
            for (int i = 0; i < _dim; i++)
                for (int j = 0; j < _dim; j++)
                    dst[i * _dim + j] = _l[i, j];
        }

        // ── Sampling ─────────────────────────────────────────────────────────────

        /// <summary>
        /// Returns a sample x = <paramref name="mean"/> + L·z where z ~ N(0, I),
        /// using the Box-Muller transform to generate standard normals.
        /// Must be called after <see cref="Decompose"/>.
        /// </summary>
        /// <param name="rng">Random number generator.</param>
        /// <param name="mean">Mean vector of length dim.</param>
        public double[] Sample(Random rng, double[] mean)
        {
            // Generate z ~ N(0, I) via Box-Muller pairs.
            double[] z = new double[_dim];
            for (int i = 0; i < _dim - 1; i += 2)
            {
                double u1 = 1.0 - rng.NextDouble(); // (0,1] avoids log(0)
                double u2 = 1.0 - rng.NextDouble();
                double mag = Math.Sqrt(-2.0 * Math.Log(u1));
                z[i] = mag * Math.Cos(2.0 * Math.PI * u2);
                z[i + 1] = mag * Math.Sin(2.0 * Math.PI * u2);
            }
            if (_dim % 2 == 1)
            {
                double u1 = 1.0 - rng.NextDouble();
                double u2 = 1.0 - rng.NextDouble();
                z[_dim - 1] = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
            }

            // x = μ + L·z  (L is lower-triangular: _l[i, j] = 0 for j > i)
            double[] x = new double[_dim];
            for (int i = 0; i < _dim; i++)
            {
                x[i] = mean[i];
                for (int j = 0; j <= i; j++)
                    x[i] += _l[i, j] * z[j];
            }
            return x;
        }
    }
}
