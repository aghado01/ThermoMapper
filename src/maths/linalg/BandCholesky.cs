using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// Which square-root convention the band factorization uses. Both factor the same SPD banded matrix; they
    /// differ only in whether a per-pivot <c>sqrt</c> is taken.
    /// </summary>
    public enum BandFactorization
    {
        /// <summary>A = L·Lᵀ, L lower-triangular with positive diagonal — one <c>sqrt</c> per pivot.</summary>
        Cholesky,

        /// <summary>
        /// A = L·D·Lᵀ, L <i>unit</i> lower-triangular and D diagonal — square-root-free. The solve carries no
        /// division through the triangular sweeps (only the diagonal scale by D), so it is the lighter inner-loop
        /// kernel; preferred when many small factorizations dominate (the MCMC pattern).
        /// </summary>
        Ldlt
    }

    /// <summary>
    /// Band factorization of a symmetric positive-definite matrix with lower half-bandwidth p — the Cholesky
    /// family restricted to the band (Golub &amp; Van Loan, <i>Matrix Computations</i> §4.3; the LAPACK
    /// <c>pbtrf</c>/<c>pbtrs</c> pair). The factor inherits the same bandwidth, so factorization is O(n·p²) and a
    /// solve O(n·p) — versus dense O(n³)/O(n²). Two specializations share the storage and the band geometry:
    /// <see cref="BandFactorization.Cholesky"/> (L·Lᵀ) and <see cref="BandFactorization.Ldlt"/> (root-free L·D·Lᵀ).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Built for the free-knot spline normal equations BᵀWB, which are banded at p = spline degree (heptadiagonal
    /// for cubics): factor once, then solve A·x = b by two banded triangular substitutions. The marginal never
    /// needs A⁻¹ — both the quadratic form bᵀA⁻¹b and the posterior coefficients A⁻¹b are recovered from the
    /// <see cref="Solve(double[])"/> alone, so no inverse is ever materialized.
    /// </para>
    /// <para>
    /// Allocate once with <see cref="BandCholesky(int,int,BandFactorization)"/>, then call <see cref="Decompose"/>
    /// each time the matrix changes (the inner-loop MCMC pattern). Mirrors <see cref="CholeskyDecomposition"/>'s
    /// 1e-12 diagonal floor — a soft guard against a near-interpolating (over-knotted) design without throwing;
    /// the floor is applied to the variance so the two conventions clamp to the same D = 1e-12.
    /// </para>
    /// </remarks>
    public sealed class BandCholesky
    {
        private readonly int _n;
        private readonly int _p;                 // lower half-bandwidth (sub-diagonals); L has the same.
        private readonly BandFactorization _mode;

        // Lower-band factor in LAPACK 'L' layout: _l[d, j] = L(j+d, j) for d = 1..p (sub-diagonals).
        // The diagonal slot _l[0, j] holds L(j,j) for Cholesky, or D(j) for LDLᵀ (where L is unit-diagonal).
        private readonly double[,] _l;

        /// <summary>ln|A|. For Cholesky 2·Σ_j ln(L_{jj}); for LDLᵀ Σ_j ln(D_j). Valid after <see cref="Decompose"/>.</summary>
        public double LogDet { get; private set; }

        /// <summary>
        /// True if any pivot was clamped to the diagonal floor during the last factorization — i.e. the matrix is
        /// (numerically) singular / rank-deficient. A consumer should treat the system as unfittable rather than
        /// trust the resulting solve. Set per call by <see cref="Decompose"/> / <see cref="DecomposeBanded"/>.
        /// </summary>
        public bool HitFloor { get; private set; }

        /// <summary>Matrix order n.</summary>
        public int Dimension => _n;

        /// <summary>Lower half-bandwidth p (number of sub-diagonals); for a degree-d spline this is d.</summary>
        public int Bandwidth => _p;

        /// <summary>The factorization convention in force.</summary>
        public BandFactorization Factorization => _mode;

        /// <param name="dim">Matrix order n.</param>
        /// <param name="bandwidth">Lower half-bandwidth p; clamped to n−1 (a dense matrix is the limiting case).</param>
        /// <param name="factorization">Square-root convention; defaults to <see cref="BandFactorization.Cholesky"/>.</param>
        public BandCholesky(int dim, int bandwidth, BandFactorization factorization = BandFactorization.Cholesky)
        {
            if (dim <= 0) throw new ArgumentOutOfRangeException(nameof(dim), "Dimension must be positive.");
            if (bandwidth < 0) throw new ArgumentOutOfRangeException(nameof(bandwidth), "Bandwidth must be non-negative.");
            _n = dim;
            _p = Math.Min(bandwidth, dim - 1);
            _mode = factorization;
            _l = new double[_p + 1, dim];
        }

        // ── Decompose ────────────────────────────────────────────────────────────

        /// <summary>
        /// Decomposes <paramref name="a"/> over the band, updating <see cref="LogDet"/> and preparing the factor
        /// used by <see cref="Solve(double[])"/>. Only the lower band (|i−j| ≤ p) is read; entries outside it are
        /// assumed zero (the band structure is a precondition, not enforced).
        /// </summary>
        /// <param name="a">
        /// Square, symmetric, positive-definite matrix of order n with lower half-bandwidth ≤ p. Only the lower
        /// triangle within the band is read.
        /// </param>
        public void Decompose(double[,] a)
        {
            ArgumentNullException.ThrowIfNull(a);
            for (int j = 0; j < _n; j++)                   // load _l from the dense band, then factor in place
            {
                _l[0, j] = a[j, j];
                int dmax = Math.Min(_p, _n - 1 - j);
                for (int d = 1; d <= dmax; d++)
                    _l[d, j] = a[j + d, j];
            }
            FactorInPlace();
        }

        /// <summary>
        /// Decomposes a matrix already in LAPACK lower-band storage: <paramref name="band"/><c>[d, j] = A(j+d, j)</c>
        /// for d = 0..p. This is the allocation-lean entry point — the caller accumulates the band directly,
        /// never materializing the dense ν×ν matrix. Updates <see cref="LogDet"/>; prepares <see cref="Solve(double[])"/>.
        /// </summary>
        public void DecomposeBanded(double[,] band)
        {
            ArgumentNullException.ThrowIfNull(band);
            for (int j = 0; j < _n; j++)                   // load _l from the supplied band, then factor in place
            {
                _l[0, j] = band[0, j];
                int dmax = Math.Min(_p, _n - 1 - j);
                for (int d = 1; d <= dmax; d++)
                    _l[d, j] = band[d, j];
            }
            FactorInPlace();
        }

        private void FactorInPlace()
        {
            if (_mode == BandFactorization.Ldlt) FactorLdlt();
            else FactorCholesky();
        }

        // A = L·Lᵀ, in place over _l (which holds A on entry). Left-looking column Cholesky restricted to the
        // band: column j reads only already-factored columns k ∈ [max(0, i−p), j) — at most p — so O(p²)/column.
        private void FactorCholesky()
        {
            double logDetL = 0.0;
            HitFloor = false;
            for (int j = 0; j < _n; j++)
            {
                double sum = _l[0, j];                     // A(j,j)
                int kmin = Math.Max(0, j - _p);
                for (int k = kmin; k < j; k++)
                {
                    double ljk = _l[j - k, k];            // L(j,k)
                    sum -= ljk * ljk;
                }
                if (sum <= 1e-12) { sum = 1e-12; HitFloor = true; }   // diagonal floor ⇒ singular
                double ljj = Math.Sqrt(sum);
                _l[0, j] = ljj;
                logDetL += Math.Log(ljj);

                int dmax = Math.Min(_p, _n - 1 - j);
                for (int d = 1; d <= dmax; d++)
                {
                    int i = j + d;
                    double s = _l[d, j];                   // A(i,j)
                    int kk = Math.Max(0, i - _p);
                    for (int k = kk; k < j; k++)
                        s -= _l[i - k, k] * _l[j - k, k];  // L(i,k)·L(j,k)
                    _l[d, j] = s / ljj;
                }
            }
            LogDet = 2.0 * logDetL;                        // ln|A| = 2·ln|L|
        }

        // A = L·D·Lᵀ, L unit lower-triangular, in place over _l. No sqrt: the diagonal slot holds D(j).
        private void FactorLdlt()
        {
            double logDet = 0.0;
            HitFloor = false;
            for (int j = 0; j < _n; j++)
            {
                double dj = _l[0, j];                      // A(j,j)
                int kmin = Math.Max(0, j - _p);
                for (int k = kmin; k < j; k++)
                {
                    double ljk = _l[j - k, k];            // L(j,k)
                    dj -= ljk * ljk * _l[0, k];           // L(j,k)²·D(k)
                }
                if (dj <= 1e-12) { dj = 1e-12; HitFloor = true; }   // diagonal floor ⇒ singular
                _l[0, j] = dj;                             // store D(j)
                logDet += Math.Log(dj);

                int dmax = Math.Min(_p, _n - 1 - j);
                for (int d = 1; d <= dmax; d++)
                {
                    int i = j + d;
                    double s = _l[d, j];                   // A(i,j)
                    int kk = Math.Max(0, i - _p);
                    for (int k = kk; k < j; k++)
                        s -= _l[i - k, k] * _l[j - k, k] * _l[0, k];  // L(i,k)·L(j,k)·D(k)
                    _l[d, j] = s / dj;                     // L(i,j); L(j,j)=1 implicit
                }
            }
            LogDet = logDet;                               // ln|A| = Σ ln D(j)
        }

        // ── Solve ─────────────────────────────────────────────────────────────────

        /// <summary>
        /// Solves A·x = b and returns x, leaving <paramref name="b"/> untouched. Banded forward + back
        /// substitution. Must be called after <see cref="Decompose"/>.
        /// </summary>
        public double[] Solve(double[] b)
        {
            ArgumentNullException.ThrowIfNull(b);
            if (b.Length != _n) throw new ArgumentException("Right-hand side length must equal the matrix order.", nameof(b));
            var x = (double[])b.Clone();
            SolveInPlace(x);
            return x;
        }

        /// <summary>
        /// In-place band solve: <paramref name="x"/> enters as the right-hand side b and exits as A⁻¹b.
        /// Must be called after <see cref="Decompose"/>.
        /// </summary>
        public void SolveInPlace(double[] x)
        {
            ArgumentNullException.ThrowIfNull(x);
            if (x.Length != _n) throw new ArgumentException("Right-hand side length must equal the matrix order.", nameof(x));
            if (_mode == BandFactorization.Ldlt) SolveLdlt(x);
            else SolveCholesky(x);
        }

        // L·Lᵀ solve: divide by L(i,i) in both sweeps.
        private void SolveCholesky(double[] x)
        {
            for (int i = 0; i < _n; i++)               // forward: L·y = b
            {
                double s = x[i];
                int kmin = Math.Max(0, i - _p);
                for (int k = kmin; k < i; k++)
                    s -= _l[i - k, k] * x[k];          // L(i,k)
                x[i] = s / _l[0, i];                   // L(i,i)
            }
            for (int i = _n - 1; i >= 0; i--)          // back: Lᵀ·x = y
            {
                double s = x[i];
                int kmax = Math.Min(_n - 1, i + _p);
                for (int k = i + 1; k <= kmax; k++)
                    s -= _l[k - i, i] * x[k];          // L(k,i) = (Lᵀ)(i,k)
                x[i] = s / _l[0, i];
            }
        }

        // L·D·Lᵀ solve: unit-triangular sweeps carry no division; only the middle scale by D.
        private void SolveLdlt(double[] x)
        {
            for (int i = 0; i < _n; i++)               // forward: L·z = b (unit diagonal)
            {
                double s = x[i];
                int kmin = Math.Max(0, i - _p);
                for (int k = kmin; k < i; k++)
                    s -= _l[i - k, k] * x[k];          // L(i,k)
                x[i] = s;                              // L(i,i) = 1
            }
            for (int i = 0; i < _n; i++)               // diagonal: y = z / D
                x[i] /= _l[0, i];                      // D(i)
            for (int i = _n - 1; i >= 0; i--)          // back: Lᵀ·x = y (unit diagonal)
            {
                double s = x[i];
                int kmax = Math.Min(_n - 1, i + _p);
                for (int k = i + 1; k <= kmax; k++)
                    s -= _l[k - i, i] * x[k];          // L(k,i)
                x[i] = s;
            }
        }

        // ── Gaussian draw ───────────────────────────────────────────────────────

        /// <summary>
        /// Maps i.i.d. standard normals <paramref name="z"/> to a draw <c>v = L⁻ᵀ D^{−½} z ~ N(0, A⁻¹)</c> — the
        /// innovation for posterior sampling <c>β ~ N(μ, A⁻¹)</c> as <c>β = μ + v</c> (so for the conjugate
        /// P-spline draw <c>β ~ N(β̂, σ² A⁻¹)</c>, scale <paramref name="z"/> by σ). LDLᵀ factorization only — the
        /// root-free factor is what makes <c>D^{−½}</c> a single diagonal scale. <paramref name="result"/> and
        /// <paramref name="z"/> may alias. Must be called after <see cref="DecomposeBanded"/>/<see cref="Decompose"/>.
        /// </summary>
        public void SampleInnovation(double[] z, double[] result)
        {
            ArgumentNullException.ThrowIfNull(z);
            ArgumentNullException.ThrowIfNull(result);
            if (_mode != BandFactorization.Ldlt)
                throw new InvalidOperationException("Gaussian sampling requires the Ldlt factorization.");
            if (z.Length != _n || result.Length != _n)
                throw new ArgumentException("Vector length must equal the matrix order.");

            for (int i = 0; i < _n; i++)
                result[i] = z[i] / Math.Sqrt(_l[0, i]);   // D^{−½} z

            for (int i = _n - 1; i >= 0; i--)             // solve Lᵀ result = D^{−½} z (unit diagonal)
            {
                double s = result[i];
                int kmax = Math.Min(_n - 1, i + _p);
                for (int k = i + 1; k <= kmax; k++)
                    s -= _l[k - i, i] * result[k];        // L(k,i)
                result[i] = s;
            }
        }
    }
}
