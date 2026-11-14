using System;
using Maths.LinAlg;

namespace Maths.Geometry
{
    /// <summary>
    /// Stiefel manifold V_{n,r}: orthonormal r-frames in Rⁿ (n×r matrices with orthonormal columns), stored column-major as a flat
    /// length-(n·r) buffer, as <see cref="GrassmannManifold"/> stores its representatives. Unlike the Grassmannian it sees the frame:
    /// Y and Y·Q are different points for every r×r orthogonal Q ≠ I.
    ///
    /// <para><b>The canonical metric</b> of Edelman, Arias and Smith (1998), the quotient of the orthogonal group's bi-invariant metric
    /// by O(n − r) (Absil, Mahony and Sepulchre 2008, chapter 5's notes, read; the formula from EAS, cited from memory):
    /// ⟨Δ, Δ⟩ = tr(Δᵀ(I − ½YYᵀ)Δ). A tangent vector Δ = YA + QB, with A skew and Q orthonormal and orthogonal to Y, has squared norm
    /// ½‖A‖² + ‖B‖², so turning the frame within its span weighs half what moving the span does.</para>
    ///
    /// <para><b>Exp</b> is the closed-form geodesic Y(t) = [Y Q] exp(t·[A, −Bᵀ; B, 0]) [I; 0], the first r columns of the orthogonal
    /// group's one-parameter subgroup through the horizontal lift. <b>Log</b> is Zimmermann's iteration (2017, cited from memory):
    /// complete [YᵀZ; QᵀZ] to a rotation V, then turn V's last block by exp(−C) until the lower-right block C of log V vanishes, which
    /// makes log V a horizontal lift and its first block column the tangent of a geodesic from Y to Z. The completion is first aligned
    /// to a symmetric lower block (a Procrustes step), so that V starts near the identity for nearby frames. Q spans only the part of
    /// Z off Y's span, so a pure turn within the span is log(YᵀZ) directly. The exponentials and logarithms of the small skew and
    /// rotation blocks are <see cref="MatrixFunctions.SkewExponential"/> and <see cref="MatrixFunctions.TryRotationLogarithm"/>,
    /// read from symmetric eigensystems, so no angle passes through an arccos near zero and a small distance keeps its digits.</para>
    ///
    /// <para><b>Domain.</b> Log converges for frames that are not far apart: in a sweep over n ≤ 8 and r ≤ 4, every pair at canonical
    /// distance up to 2.7, with the round trip to about 2·10⁻¹⁴, and most at 2.9. It throws <see cref="InvalidOperationException"/> where
    /// the iteration does not converge, where a rotation block reaches a half-turn, or where Z = YQ with det Q = −1: all far-apart
    /// frames, the last on or past the cut locus, where the log is not unique. <see cref="Distance"/> is the norm of the log, so it
    /// throws there too.
    /// Exp takes any tangent vector, and reads a vector off the tangent space through its tangent part.</para>
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

        // A direction of Z off Y's span whose sine is at most this is dropped from the log's normal block: the endpoint moves by at
        // most this much.
        private const double LogRankFloor = 1e-14;

        // Exp keeps a normal direction whose singular value exceeds this fraction of ‖Δ‖_F, above the rounding the projection leaves
        // (about ε‖Δ‖), so a tiny step keeps its shape and a turn within the span (r = n included) adds no spurious direction.
        private const double ExpRankFraction = 1e-13;

        // Zimmermann's iteration converges linearly, so it stops once ‖C‖_F ≤ LogTolerance·(1 + ‖log V‖_F), or once ‖C‖_F stops
        // falling below LogFloorTolerance·(1 + ‖log V‖_F), the rounding floor of the eigensystem-based logarithm.
        private const double LogTolerance = 1e-15;
        private const double LogFloorTolerance = 1e-11;

        private const int LogMaxIterations = 1000;

        /// <summary>d(Y, Z) = ‖log_Y(Z)‖ in the canonical metric, √(½‖A‖² + ‖B‖²) from the log's blocks.</summary>
        public double Distance(ReadOnlySpan<double> p, ReadOnlySpan<double> q)
        {
            Log(p, q, out double[,] a, out _, out double[,] b, out _);
            return Math.Sqrt(0.5 * FrobeniusSquared(a) + FrobeniusSquared(b));
        }

        /// <summary>log_Y(Z) = YA + QB, the tangent of the geodesic from Y to Z found by Zimmermann's iteration.</summary>
        public void LogMap(ReadOnlySpan<double> p, ReadOnlySpan<double> q, Span<double> dst)
        {
            int n = _n, r = _r;
            Log(p, q, out double[,] a, out double[] basis, out double[,] b, out int k);
            for (int j = 0; j < r; j++)
                for (int i = 0; i < n; i++)
                {
                    double s = 0.0;
                    for (int l = 0; l < r; l++) s += p[i + l * n] * a[l, j];
                    for (int l = 0; l < k; l++) s += basis[i + l * n] * b[l, j];
                    dst[i + j * n] = s;
                }
        }

        /// <summary>
        /// exp_Y(Δ) = [Y Q] exp([A, −Rᵀ; R, 0]) [I; 0], with A = skew(YᵀΔ) and QR the normal part (I − YYᵀ)Δ, Q taken from its thin
        /// SVD over the directions it moves.
        /// </summary>
        public void ExpMap(ReadOnlySpan<double> p, ReadOnlySpan<double> v, Span<double> dst)
        {
            int n = _n, r = _r;
            double[,] ytv = AtB(p, v, n, r, r);
            var a = new double[r, r];
            for (int i = 0; i < r; i++)
                for (int j = i + 1; j < r; j++)
                {
                    double s = 0.5 * (ytv[i, j] - ytv[j, i]);
                    a[i, j] = s;
                    a[j, i] = -s;
                }

            double[] w = v.ToArray();
            SubtractSpanTwice(p, w, n, r);
            ThinSvdResult svd = ThinSvd.Decompose(ToMatrix(w, n, r));
            double size = 0.0;
            for (int i = 0; i < Dimension; i++) size += v[i] * v[i];
            double floor = ExpRankFraction * Math.Sqrt(size);
            int k = 0;
            while (k < Math.Min(svd.SingularValues.Length, n - r) && svd.SingularValues[k] > floor) k++;
            double[] basis = NormalBasis(p, svd, k, n, r);
            double[,] rr = AtB(basis, w, n, k, r);

            int m = r + k;
            var omega = new double[m, m];
            for (int i = 0; i < r; i++)
                for (int j = 0; j < r; j++) omega[i, j] = a[i, j];
            for (int i = 0; i < k; i++)
                for (int j = 0; j < r; j++)
                {
                    omega[r + i, j] = rr[i, j];
                    omega[j, r + i] = -rr[i, j];
                }
            double[,] e = MatrixFunctions.SkewExponential(omega);

            for (int j = 0; j < r; j++)
                for (int i = 0; i < n; i++)
                {
                    double s = 0.0;
                    for (int l = 0; l < r; l++) s += p[i + l * n] * e[l, j];
                    for (int l = 0; l < k; l++) s += basis[i + l * n] * e[r + l, j];
                    dst[i + j * n] = s;
                }
        }

        public void AddScaled(Span<double> dst, ReadOnlySpan<double> v, double scalar)
        {
            for (int i = 0; i < Dimension; i++) dst[i] += v[i] * scalar;
        }

        /// <summary>The canonical norm √(tr(Δᵀ(I − ½YYᵀ)Δ)) = √(‖Δ‖²_F − ½‖YᵀΔ‖²_F).</summary>
        public double Norm(ReadOnlySpan<double> p, ReadOnlySpan<double> v)
        {
            double total = 0.0;
            for (int i = 0; i < Dimension; i++) total += v[i] * v[i];
            double[,] ytv = AtB(p, v, _n, _r, _r);
            return Math.Sqrt(Math.Max(total - 0.5 * FrobeniusSquared(ytv), 0.0));
        }

        // ── the log ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────

        // The blocks of log_Y(Z) = YA + QB: A (r×r, skew), the normal basis Q (n×k, column-major) and B (k×r).
        private void Log(ReadOnlySpan<double> p, ReadOnlySpan<double> q, out double[,] a, out double[] basis, out double[,] b,
            out int k)
        {
            int n = _n, r = _r;
            double[,] m = AtB(p, q, n, r, r);
            double[] w = q.ToArray();
            SubtractSpanTwice(p, w, n, r);
            ThinSvdResult svd = ThinSvd.Decompose(ToMatrix(w, n, r));
            k = 0;
            while (k < Math.Min(svd.SingularValues.Length, n - r) && svd.SingularValues[k] > LogRankFloor) k++;
            basis = NormalBasis(p, svd, k, n, r);
            double[,] nBlock = AtB(basis, w, n, k, r);

            if (k == 0)
            {
                // Z = YM with M orthogonal: a turn within the span.
                if (DeterminantSign(m) <= 0)
                    throw new InvalidOperationException(
                        "Stiefel log: Z = YQ with det Q = −1 lies on or past the cut locus of Y; no unique geodesic.");
                if (!MatrixFunctions.TryRotationLogarithm(m, out a))
                    throw new InvalidOperationException(
                        "Stiefel log: Z turns Y's frame by a half-turn, the cut locus; no unique geodesic.");
                b = new double[0, r];
                return;
            }

            int size = r + k;
            var v = new double[size, size];
            for (int j = 0; j < r; j++)
            {
                for (int i = 0; i < r; i++) v[i, j] = m[i, j];
                for (int i = 0; i < k; i++) v[r + i, j] = nBlock[i, j];
            }
            Complete(v, r, k);

            double[,] log = null!;
            bool converged = false;
            double previous = double.PositiveInfinity;
            for (int iteration = 0; iteration < LogMaxIterations; iteration++)
            {
                if (!MatrixFunctions.TryRotationLogarithm(v, out log))
                    throw new InvalidOperationException(
                        "Stiefel log: the completed rotation reached a half-turn; Z is too far from Y.");
                double cSquared = 0.0, total = 0.0;
                for (int i = 0; i < size; i++)
                    for (int j = 0; j < size; j++)
                    {
                        double x = log[i, j] * log[i, j];
                        total += x;
                        if (i >= r && j >= r) cSquared += x;
                    }
                double c = Math.Sqrt(cSquared), scale = 1.0 + Math.Sqrt(total);
                if (c <= LogTolerance * scale || (c >= previous && c <= LogFloorTolerance * scale))
                {
                    converged = true;
                    break;
                }
                previous = c;
                var minusC = new double[k, k];
                for (int i = 0; i < k; i++)
                    for (int j = 0; j < k; j++) minusC[i, j] = -log[r + i, r + j];
                double[,] phi = MatrixFunctions.SkewExponential(minusC);
                TurnLastBlock(v, r, k, phi);
            }
            if (!converged)
                throw new InvalidOperationException(
                    $"Stiefel log: Zimmermann's iteration did not converge in {LogMaxIterations} steps; Z is too far from Y.");

            a = new double[r, r];
            b = new double[k, r];
            for (int i = 0; i < r; i++)
                for (int j = 0; j < r; j++) a[i, j] = log[i, j];
            for (int i = 0; i < k; i++)
                for (int j = 0; j < r; j++) b[i, j] = log[r + i, j];
        }

        // Fills v's last k columns with an orthonormal completion of its first r, aligned so the lower k×k block is symmetric positive
        // semidefinite, and with det v = +1 (flipping the completion's least-weighted direction when it comes out −1).
        private static void Complete(double[,] v, int r, int k)
        {
            int size = r + k;
            var projector = new double[size, size];
            for (int i = 0; i < size; i++)
                for (int j = i; j < size; j++)
                {
                    double s = i == j ? 1.0 : 0.0;
                    for (int l = 0; l < r; l++) s -= v[i, l] * v[j, l];
                    projector[i, j] = projector[j, i] = s;
                }
            EigenResult eigen = Eigen.DecomposeSymmetric(projector, tol: 1e-14);
            int[] order = new int[size];
            for (int i = 0; i < size; i++) order[i] = i;
            Array.Sort(order, (x, y) => eigen.Eigenvalues[y].CompareTo(eigen.Eigenvalues[x]));
            var completion = new double[size, k];
            for (int c = 0; c < k; c++)
            {
                double[] vec = eigen.Eigenvectors[order[c]];
                for (int i = 0; i < size; i++) completion[i, c] = vec[i];
            }

            // Procrustes alignment: with the lower block Y₀ = U S Wᵀ, right-multiplying by W Uᵀ makes it U S Uᵀ.
            var lower = new double[k, k];
            for (int i = 0; i < k; i++)
                for (int j = 0; j < k; j++) lower[i, j] = completion[r + i, j];
            ThinSvdResult svd = ThinSvd.Decompose(lower);
            var turn = new double[k, k];
            for (int i = 0; i < k; i++)
                for (int j = 0; j < k; j++)
                {
                    double s = 0.0;
                    for (int l = 0; l < k; l++) s += svd.RightVectors[l][i] * svd.LeftVectors[l][j];
                    turn[i, j] = s;
                }
            for (int i = 0; i < size; i++)
            {
                var row = new double[k];
                for (int j = 0; j < k; j++)
                {
                    double s = 0.0;
                    for (int l = 0; l < k; l++) s += completion[i, l] * turn[l, j];
                    row[j] = s;
                }
                for (int j = 0; j < k; j++) v[i, r + j] = row[j];
            }

            if (DeterminantSign(v) < 0)
            {
                // Reflect along the lower block's least singular direction u: columns ← columns·(I − 2uuᵀ).
                double[] u = svd.LeftVectors[k - 1];
                for (int i = 0; i < size; i++)
                {
                    double s = 0.0;
                    for (int l = 0; l < k; l++) s += v[i, r + l] * u[l];
                    for (int j = 0; j < k; j++) v[i, r + j] -= 2.0 * s * u[j];
                }
            }
        }

        // v's last k columns ← those columns · phi.
        private static void TurnLastBlock(double[,] v, int r, int k, double[,] phi)
        {
            int size = r + k;
            var row = new double[k];
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < k; j++)
                {
                    double s = 0.0;
                    for (int l = 0; l < k; l++) s += v[i, r + l] * phi[l, j];
                    row[j] = s;
                }
                for (int j = 0; j < k; j++) v[i, r + j] = row[j];
            }
        }

        // ── small column-major linear algebra ──────────────────────────────────────────────────────────────────────────────────────

        // Aᵀ B for column-major A (n×ra) and B (n×rb): an ra×rb matrix.
        private static double[,] AtB(ReadOnlySpan<double> a, ReadOnlySpan<double> b, int n, int ra, int rb)
        {
            var c = new double[ra, rb];
            for (int i = 0; i < ra; i++)
                for (int j = 0; j < rb; j++)
                {
                    double s = 0.0;
                    int ai = i * n, bj = j * n;
                    for (int l = 0; l < n; l++) s += a[ai + l] * b[bj + l];
                    c[i, j] = s;
                }
            return c;
        }

        // w ← w − Y(Yᵀw), twice, so w is orthogonal to Y's span to rounding.
        private static void SubtractSpanTwice(ReadOnlySpan<double> y, double[] w, int n, int r)
        {
            for (int pass = 0; pass < 2; pass++)
            {
                double[,] ytw = AtB(y, w, n, r, r);
                for (int j = 0; j < r; j++)
                    for (int l = 0; l < r; l++)
                    {
                        double c = ytw[l, j];
                        if (c == 0.0) continue;
                        for (int i = 0; i < n; i++) w[i + j * n] -= y[i + l * n] * c;
                    }
            }
        }

        // The first k left singular vectors, made orthonormal and orthogonal to Y by Gram–Schmidt behind Y's columns.
        private static double[] NormalBasis(ReadOnlySpan<double> y, ThinSvdResult svd, int k, int n, int r)
        {
            if (k == 0) return Array.Empty<double>();
            var block = new double[n * (r + k)];
            y.CopyTo(block);
            for (int c = 0; c < k; c++) svd.LeftVectors[c].CopyTo(block, (r + c) * n);
            MatrixOps.Orthonormalize(block, n, r + k);
            return block.AsSpan(r * n, k * n).ToArray();
        }

        private static double[,] ToMatrix(double[] flat, int n, int r)
        {
            var m = new double[n, r];
            for (int j = 0; j < r; j++)
                for (int i = 0; i < n; i++) m[i, j] = flat[i + j * n];
            return m;
        }

        private static double FrobeniusSquared(double[,] m)
        {
            double s = 0.0;
            foreach (double x in m) s += x * x;
            return s;
        }

        // The sign of det(a) by Gaussian elimination with partial pivoting; 0 for a singular a.
        private static int DeterminantSign(double[,] a)
        {
            int size = a.GetLength(0);
            var lu = (double[,])a.Clone();
            int sign = 1;
            for (int col = 0; col < size; col++)
            {
                int pivot = col;
                for (int i = col + 1; i < size; i++)
                    if (Math.Abs(lu[i, col]) > Math.Abs(lu[pivot, col])) pivot = i;
                if (lu[pivot, col] == 0.0) return 0;
                if (pivot != col)
                {
                    sign = -sign;
                    for (int j = 0; j < size; j++) (lu[col, j], lu[pivot, j]) = (lu[pivot, j], lu[col, j]);
                }
                if (lu[col, col] < 0.0) sign = -sign;
                for (int i = col + 1; i < size; i++)
                {
                    double f = lu[i, col] / lu[col, col];
                    if (f == 0.0) continue;
                    for (int j = col; j < size; j++) lu[i, j] -= f * lu[col, j];
                }
            }
            return sign;
        }
    }
}
