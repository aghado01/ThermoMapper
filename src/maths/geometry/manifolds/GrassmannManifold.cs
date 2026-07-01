using System;
using Maths.LinAlg;

namespace Maths.Geometry
{
    /// <summary>
    /// Grassmann manifold Gr(r, n): the r-dimensional linear subspaces of Rⁿ. A point is an
    /// n×r orthonormal basis (a Stiefel representative), stored column-major as a flat
    /// length-(n·r) buffer [col0 | col1 | … | col_{r−1}] — matching <see cref="MatrixOps"/>.
    /// Every operation is gauge-invariant: it depends only on the span, so Y and Y·Q represent
    /// the same point for any r×r orthogonal Q.
    ///
    /// <para>Geometry follows Edelman, Arias &amp; Smith (1998) and Bendokat–Zimmermann–Absil
    /// (2020). The geodesic distance is the 2-norm of the principal angles; Log/Exp use a thin
    /// SVD built here from the symmetric eigendecomposition of the r×r Gram matrix
    /// (<see cref="DenseEigen"/>). The canonical metric gives ‖log_Y(Z)‖_F = d(Y, Z), mirroring
    /// the embedded <see cref="SphericalManifold"/>.</para>
    ///
    /// <para><b>Domain.</b> Log is defined inside the injectivity radius — principal angles &lt; π/2
    /// (YᵀZ nonsingular). For the geometric-median / Karcher use case the data subspaces sit near
    /// the running estimate, well inside that ball.</para>
    /// </summary>
    public readonly struct GrassmannManifold : IRiemannianManifold
    {
        public static bool IsFlat => false;

        private readonly int _n;
        private readonly int _r;

        /// <summary>Ambient buffer length n·r (points and tangents are n×r column-major).</summary>
        public int Dimension { get; }

        /// <summary>Intrinsic manifold dimension r(n−r). Feed THIS to curvature/dimension
        /// consumers, never <see cref="Dimension"/> (cf. the SphericalManifold ward).</summary>
        public int IntrinsicDimension => _r * (_n - _r);

        public GrassmannManifold(int ambientN, int subspaceR)
        {
            if (subspaceR < 1 || subspaceR > ambientN)
                throw new ArgumentOutOfRangeException(nameof(subspaceR),
                    "Grassmann subspace dimension must satisfy 1 ≤ r ≤ n.");
            _n = ambientN;
            _r = subspaceR;
            Dimension = ambientN * subspaceR;
        }

        /// <summary>d(Y,Z) = ‖Θ‖₂, principal angles θ_i = arccos(σ_i(YᵀZ)).</summary>
        public double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
        {
            int n = _n, r = _r;
            double[] m = MulAtB(p, n, r, q, r);            // M = YᵀZ  (r×r)
            double[] cos = SingularValues(m, r);           // σ_i(M) = cos θ_i (descending)
            double sumSq = 0.0;
            for (int i = 0; i < r; i++)
            {
                double theta = Math.Acos(Math.Clamp(cos[i], -1.0, 1.0));
                sumSq += theta * theta;
            }
            return Math.Sqrt(sumSq);
        }

        /// <summary>
        /// log_Y(Z) = U₂ Θ U₁ᵀ, where M = YᵀZ = U₁ cosΘ V₁ᵀ and the orthogonal-complement
        /// directions are u_i = ((Z V₁)_i − cosθ_i (Y U₁)_i)/sinθ_i.
        /// </summary>
        public void LogMap(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
        {
            int n = _n, r = _r;
            double[] m = MulAtB(p, n, r, q, r);            // M = YᵀZ (r×r)
            ThinSvd(m, r, r, out double[] u1, out double[] c, out double[] v1);  // M = U₁ diag(c) V₁ᵀ

            double[] zv1 = Mul(q, n, r, v1, r);            // Z V₁ (n×r) → b_i
            double[] yu1 = Mul(p, n, r, u1, r);            // Y U₁ (n×r) → a_i

            double[] u2theta = new double[n * r];
            for (int i = 0; i < r; i++)
            {
                double ci = Math.Clamp(c[i], -1.0, 1.0);
                double th = Math.Acos(ci);
                double sin = Math.Sin(th);
                int col = i * n;
                if (sin < 1e-12) continue;                 // θ_i ≈ 0: no motion, u_i = 0
                double scale = th / sin;                   // fold Θ into U₂ directly
                for (int row = 0; row < n; row++)
                    u2theta[col + row] = (zv1[col + row] - ci * yu1[col + row]) * scale;
            }

            double[] logm = MulABt(u2theta, n, r, u1, r);  // (U₂ diagΘ) U₁ᵀ → n×r
            logm.AsSpan(0, n * r).CopyTo(dst);
        }

        /// <summary>exp_Y(Δ) = (Y V cosΣ + U sinΣ) Vᵀ for Δ = U Σ Vᵀ, re-orthonormalized.</summary>
        public void ExpMap(ReadOnlySpan<double> p, ReadOnlySpan<double> v, Span<double> dst)
        {
            int n = _n, r = _r;
            ThinSvd(v, n, r, out double[] u, out double[] sig, out double[] vmat);  // Δ = U diag(sig) Vᵀ

            double[] yv = Mul(p, n, r, vmat, r);           // Y V (n×r)
            double[] term = new double[n * r];
            for (int i = 0; i < r; i++)
            {
                double cos = Math.Cos(sig[i]);
                double sin = Math.Sin(sig[i]);
                int col = i * n;
                for (int row = 0; row < n; row++)
                    term[col + row] = yv[col + row] * cos + u[col + row] * sin;  // Y V cosΣ + U sinΣ
            }
            double[] res = MulABt(term, n, r, vmat, r);    // · Vᵀ → n×r
            MatrixOps.Orthonormalize(res, n, r);           // guard drift off the manifold
            res.AsSpan(0, n * r).CopyTo(dst);
        }

        public void AddScaled(Span<double> dst, ReadOnlySpan<double> v, double scalar)
        {
            for (int i = 0; i < Dimension; i++) dst[i] += v[i] * scalar;
        }

        /// <summary>Canonical Grassmann metric: ‖Δ‖_F, so ‖log_Y(Z)‖ = d(Y, Z).</summary>
        public double Norm(ReadOnlySpan<double> p, ReadOnlySpan<double> v)
        {
            double s = 0.0;
            for (int i = 0; i < Dimension; i++) s += v[i] * v[i];
            return Math.Sqrt(s);
        }

        // ── small column-major linear algebra ───────────────────────────────────────────

        // Aᵀ B : A is n×ra, B is n×rb (column-major) → ra×rb column-major.
        private static double[] MulAtB(ReadOnlySpan<double> a, int n, int ra, ReadOnlySpan<double> b, int rb)
        {
            var c = new double[ra * rb];
            for (int j = 0; j < rb; j++)
                for (int i = 0; i < ra; i++)
                {
                    double s = 0.0;
                    int ai = i * n, bj = j * n;
                    for (int k = 0; k < n; k++) s += a[ai + k] * b[bj + k];
                    c[i + j * ra] = s;
                }
            return c;
        }

        // A B : A is m×k, B is k×bc (column-major) → m×bc column-major.
        private static double[] Mul(ReadOnlySpan<double> a, int m, int k, ReadOnlySpan<double> b, int bc)
        {
            var c = new double[m * bc];
            for (int j = 0; j < bc; j++)
                for (int l = 0; l < k; l++)
                {
                    double blj = b[l + j * k];
                    if (blj == 0.0) continue;
                    int al = l * m, cj = j * m;
                    for (int i = 0; i < m; i++) c[cj + i] += a[al + i] * blj;
                }
            return c;
        }

        // A Bᵀ : A is m×k, B is p×k (column-major; Bᵀ is k×p) → m×p column-major.
        private static double[] MulABt(ReadOnlySpan<double> a, int m, int k, ReadOnlySpan<double> b, int p)
        {
            var c = new double[m * p];
            for (int l = 0; l < k; l++)
                for (int j = 0; j < p; j++)
                {
                    double bjl = b[j + l * p];
                    if (bjl == 0.0) continue;
                    int al = l * m, cj = j * m;
                    for (int i = 0; i < m; i++) c[cj + i] += a[al + i] * bjl;
                }
            return c;
        }

        // Singular values (descending) of a column-major rows×r matrix, via the r×r Gram.
        private static double[] SingularValues(ReadOnlySpan<double> a, int r)
        {
            var eig = DenseEigen.DecomposeSymmetric(Gram(a, r, r));
            var s = new double[r];
            for (int i = 0; i < r; i++) s[i] = Math.Sqrt(Math.Max(eig.Eigenvalues[i], 0.0));
            return s;
        }

        // Thin SVD A = U diag(S) Vᵀ for column-major A (rows×r, rows ≥ r) via G = AᵀA.
        // S descending; V (r×r) right vectors as columns; U (rows×r) left vectors as columns
        // (U column left zero where S_i ≈ 0).
        private static void ThinSvd(ReadOnlySpan<double> a, int rows, int r,
            out double[] u, out double[] s, out double[] v)
        {
            var eig = DenseEigen.DecomposeSymmetric(Gram(a, rows, r));   // descending
            s = new double[r];
            v = new double[r * r];
            for (int i = 0; i < r; i++)
            {
                s[i] = Math.Sqrt(Math.Max(eig.Eigenvalues[i], 0.0));
                double[] vi = eig.Eigenvectors[i];                       // i-th right singular vector
                for (int row = 0; row < r; row++) v[row + i * r] = vi[row];
            }
            double[] av = Mul(a, rows, r, v, r);                         // A V (rows×r)
            u = new double[rows * r];
            for (int i = 0; i < r; i++)
            {
                if (s[i] < 1e-12) continue;
                int col = i * rows;
                double inv = 1.0 / s[i];
                for (int row = 0; row < rows; row++) u[col + row] = av[col + row] * inv;
            }
        }

        // AᵀA as a dense [r,r] for the symmetric eigensolver. A is column-major rows×r.
        private static double[,] Gram(ReadOnlySpan<double> a, int rows, int r)
        {
            var g = new double[r, r];
            for (int i = 0; i < r; i++)
                for (int j = i; j < r; j++)
                {
                    double dot = 0.0;
                    int ai = i * rows, aj = j * rows;
                    for (int k = 0; k < rows; k++) dot += a[ai + k] * a[aj + k];
                    g[i, j] = g[j, i] = dot;
                }
            return g;
        }
    }
}
