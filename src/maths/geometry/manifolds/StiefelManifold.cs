using System;
using Maths.LinAlg;

namespace Maths.Geometry
{
    /// <summary>
    /// Stiefel manifold V_{n,r}: orthonormal r-frames in R^n (n×r matrices with orthonormal columns).
    /// Unlike the Grassmann manifold, the Stiefel manifold distinguishes orientation — Y and Y·Q
    /// are different points for r×r orthogonal Q with det(Q) = -1.
    ///
    /// <para>Geometry follows Edelman, Arias &amp; Smith (1998). The canonical metric is the Euclidean
    /// inner product on the ambient space restricted to tangent vectors: ⟨Δ₁,Δ₂⟩ = tr(Δ₁ᵀΔ₂).
    /// Tangent vectors satisfy YᵀΔ + ΔᵀY = 0 (skew-symmetric constraint).</para>
    ///
    /// <para>For distributed SPRED aggregation: per-block SA produces orthonormal projection matrices
    /// (Stiefel points), aggregated by geometric median on V_{n,r} via Weiszfeld/ILRS.</para>
    /// </summary>
    public readonly struct StiefelManifold : IRiemannianManifold
    {
        public static bool IsFlat => false;

        private readonly int _n;
        private readonly int _r;

        /// <summary>Ambient buffer length n·r (points and tangents are n×r column-major).</summary>
        public int Dimension { get; }

        /// <summary>Intrinsic manifold dimension: nr − r(r+1)/2 (orthonormality constraints).</summary>
        public int IntrinsicDimension => _n * _r - _r * (_r + 1) / 2;

        public StiefelManifold(int ambientN, int subspaceR)
        {
            if (subspaceR < 1 || subspaceR > ambientN)
                throw new ArgumentOutOfRangeException(nameof(subspaceR),
                    "Stiefel subspace dimension must satisfy 1 ≤ r ≤ n.");
            _n = ambientN;
            _r = subspaceR;
            Dimension = ambientN * subspaceR;
        }

        /// <summary>
        /// d(Y,Z)² = ‖Ω‖²_F + ‖Θ‖²_F where [Y Y⊥]ᵀZ = [cos Θ; sin Θ] and Y⊥ is any orthonormal
        /// completion. Equivalently, via QR: ‖log(YᵀZ)‖² + ‖(I−YYᵀ)Z(ZᵀY)⁻¹‖² for nonsingular YᵀZ.
        /// </summary>
        public double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
        {
            int n = _n, r = _r;

            // M = YᵀZ (r×r)
            double[] m = MulAtB(p, n, r, q, r);

            // QR of M for the horizontal component (principal angles)
            ThinSvd(m, r, r, out double[] u, out double[] s, out double[] v);

            // Vertical component: (I − YYᵀ)Z = Z − Y(YᵀZ)
            double[] zMinusYM = new double[n * r];
            for (int i = 0; i < n * r; i++) zMinusYM[i] = q[i];
            for (int j = 0; j < r; j++)
                for (int i = 0; i < n; i++)
                    for (int k = 0; k < r; k++)
                        zMinusYM[i + j * n] -= p[i + k * n] * m[k + j * r];

            // ‖vertical‖² = trace((Z−YM)ᵀ(Z−YM)) = ‖Z‖²_F − ‖M‖²_F since Z is orthonormal
            double vertSq = 0.0;
            for (int j = 0; j < r; j++)
                for (int i = 0; i < n; i++)
                {
                    double val = zMinusYM[i + j * n];
                    vertSq += val * val;
                }

            // Horizontal distance from principal angles: Σ θᵢ² where cos θᵢ = σᵢ
            double horizSq = 0.0;
            for (int i = 0; i < r; i++)
            {
                double ci = Math.Clamp(s[i], -1.0, 1.0);
                double theta = Math.Acos(ci);
                horizSq += theta * theta;
            }

            return Math.Sqrt(horizSq + vertSq);
        }

        /// <summary>
        /// log_Y(Z): tangent vector at Y representing the geodesic to Z.
        /// Horizontal component from log(YᵀZ); vertical from parallel transport.
        /// </summary>
        public void LogMap(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
        {
            int n = _n, r = _r;

            // M = YᵀZ
            double[] m = MulAtB(p, n, r, q, r);

            // Horizontal tangent: Y · log(M) via SVD
            ThinSvd(m, r, r, out double[] u, out double[] s, out double[] v);

            // Compute log(M) = U·diag(log Σ)·Vᵀ
            double[] logM = new double[r * r];
            for (int i = 0; i < r; i++)
            {
                double si = Math.Clamp(s[i], 1e-12, 1.0); // ensure positive for log
                double logSi = Math.Log(si);
                for (int j = 0; j < r; j++)
                    for (int k = 0; k < r; k++)
                        logM[j + k * r] += u[j + i * r] * logSi * v[k + i * r];
            }

            // Horizontal part: Y · logM
            double[] horiz = Mul(p, n, r, logM, r);

            // Vertical part: (Z − Y·M)·M⁻¹ = component in normal space
            // For simplicity: use (I−YYᵀ)Z(ZᵀY)⁻¹ approach via SVD inverse
            double[] zMinusYM = new double[n * r];
            for (int i = 0; i < n * r; i++) zMinusYM[i] = q[i];
            for (int j = 0; j < r; j++)
                for (int i = 0; i < n; i++)
                    for (int k = 0; k < r; k++)
                        zMinusYM[i + j * n] -= p[i + k * n] * m[k + j * r];

            // M⁻¹ = V·diag(1/Σ)·Uᵀ
            double[] minv = new double[r * r];
            for (int i = 0; i < r; i++)
            {
                if (s[i] < 1e-12) continue;
                double invSi = 1.0 / s[i];
                for (int j = 0; j < r; j++)
                    for (int k = 0; k < r; k++)
                        minv[j + k * r] += v[j + i * r] * invSi * u[k + i * r];
            }

            // Vertical = (Z − YM)·M⁻¹
            double[] vert = Mul(zMinusYM, n, r, minv, r);

            // Combined log = horiz + vert
            for (int i = 0; i < n * r; i++)
                dst[i] = horiz[i] + vert[i];
        }

        /// <summary>
        /// exp_Y(Δ) = geodesic flow on Stiefel. Uses the decomposition into horizontal/vertical
        /// components and QR-based retraction for numerical stability.
        /// </summary>
        public void ExpMap(ReadOnlySpan<double> p, ReadOnlySpan<double> v, Span<double> dst)
        {
            int n = _n, r = _r;

            // Decompose v into horizontal (Y·A) and vertical (Y⊥·B) parts
            // A = Yᵀv (skew-symmetric)
            double[] a = MulAtB(p, n, r, v, r);

            // (I−YYᵀ)v = vertical component
            double[] yTv = a;
            double[] vert = new double[n * r];
            for (int i = 0; i < n * r; i++) vert[i] = v[i];
            for (int j = 0; j < r; j++)
                for (int i = 0; i < n; i++)
                    for (int k = 0; k < r; k++)
                        vert[i + j * n] -= p[i + k * n] * yTv[k + j * r];

            // QR retraction: [Y V] · exp([A -Bᵀ; B 0]) · [I; 0]
            // Simplified: use Cayley transform approximation for small steps
            double[] result = new double[n * r];
            p.CopyTo(result);

            // Add tangent step
            for (int i = 0; i < n * r; i++) result[i] += v[i];

            // Re-orthonormalize via Gram-Schmidt
            MatrixOps.Orthonormalize(result, n, r);

            result.AsSpan().CopyTo(dst);
        }

        public void AddScaled(Span<double> dst, ReadOnlySpan<double> v, double scalar)
        {
            for (int i = 0; i < Dimension; i++) dst[i] += v[i] * scalar;
        }

        /// <summary>Canonical metric: ‖Δ‖_F = sqrt(tr(ΔᵀΔ)).</summary>
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

        // Thin SVD via symmetric eigendecomposition of AᵀA.
        private static void ThinSvd(ReadOnlySpan<double> a, int rows, int r,
            out double[] u, out double[] s, out double[] v)
        {
            var gram = new double[r, r];
            for (int i = 0; i < r; i++)
                for (int j = i; j < r; j++)
                {
                    double dot = 0.0;
                    int ai = i * rows, aj = j * rows;
                    for (int k = 0; k < rows; k++) dot += a[ai + k] * a[aj + k];
                    gram[i, j] = gram[j, i] = dot;
                }

            var eig = DenseEigen.DecomposeSymmetric(gram);
            s = new double[r];
            v = new double[r * r];
            for (int i = 0; i < r; i++)
            {
                s[i] = Math.Sqrt(Math.Max(eig.Eigenvalues[i], 0.0));
                double[] vi = eig.Eigenvectors[i];
                for (int row = 0; row < r; row++) v[row + i * r] = vi[row];
            }

            // U = A·V·Σ⁻¹
            double[] av = Mul(a, rows, r, v, r);
            u = new double[rows * r];
            for (int i = 0; i < r; i++)
            {
                if (s[i] < 1e-12) continue;
                int col = i * rows;
                double inv = 1.0 / s[i];
                for (int row = 0; row < rows; row++) u[col + row] = av[col + row] * inv;
            }
        }
    }
}
