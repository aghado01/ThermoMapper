using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// Strict dense Cholesky <c>A = L·Lᵀ</c> of a symmetric positive-definite matrix: the reference factor for dense
    /// Gram matrices (BASS's active-term Gram), the full rebuild that an incremental update is checked against, and the
    /// factor of a Gaussian's covariance (GMM's components, the Bhattacharyya distance, the scene export's ellipsoids).
    /// <list type="bullet">
    /// <item><b>Strict.</b> There is no pivot floor. The first pivot that is not finite and positive stops the
    /// factorization, <see cref="FailedAt"/> says where, and the factor refuses to solve. LinAlg reports; whether a
    /// singular design is admissible under some other prior is the caller's decision.</item>
    /// <item><b>Solve-oriented.</b> Its own reads form no inverse. Solves, <c>bᵀA⁻¹b</c>, the log-determinant and
    /// <c>N(0, A⁻¹)</c> draws all go through the two triangular solves; <c>N(μ, A)</c> draws through L
    /// (<see cref="AddSample"/>). A caller that keeps a precision matrix forms one on demand
    /// (<see cref="TryWriteInverseTo"/>), refused where it overflows.</item>
    /// <item><b>Row at a time.</b> Rows are factored in order (Cholesky–Banachiewicz) into packed row-major storage,
    /// each row the bordered extension of the factor above it: the step a birth reuses.</item>
    /// </list>
    /// Only the lower triangle of the input is read.
    /// </summary>
    public sealed class DenseCholesky : ISpdFactor
    {
        private readonly double[] _l;          // packed lower triangle: L(i, j) at i(i+1)/2 + j
        private readonly double[] _relative;   // d_j / A(j,j) per row
        private readonly double[] _sqrtDiagonal;   // √A(j,j), for the condition estimate
        private readonly double[] _columnSums;     // scratch for ‖D^(−½)·A·D^(−½)‖₁
        private double[]? _inverseFactor;          // L⁻¹, packed as L, formed on demand by TryWriteInverseTo
        private double _equilibratedNorm1 = double.NaN;
        private int _n;
        private bool _valid;

        /// <param name="capacity">Largest matrix order this instance will factor; storage is allocated once.</param>
        public DenseCholesky(int capacity)
        {
            if (capacity < 0 || capacity > 46_340)
                throw new ArgumentOutOfRangeException(nameof(capacity), "Capacity must be in [0, 46340].");
            Capacity = capacity;
            _l = new double[capacity * (capacity + 1) / 2];
            _relative = new double[capacity];
            _sqrtDiagonal = new double[capacity];
            _columnSums = new double[capacity];
        }

        /// <summary>Largest order this instance can factor.</summary>
        public int Capacity { get; }

        /// <summary>Order of the last factored matrix.</summary>
        public int Dimension => _n;

        /// <summary>What the last <see cref="Decompose"/> found, up to and including any failing pivot.</summary>
        public FactorStatus Status { get; private set; }

        /// <summary>Row at which the last factorization stopped; −1 when the matrix is positive definite.</summary>
        public int FailedAt { get; private set; } = -1;

        /// <summary>ln|A| = 2·Σ ln L(i,i); NaN unless the last factorization succeeded.</summary>
        public double LogDet { get; private set; } = double.NaN;

        /// <summary>Factors <paramref name="a"/>, reading its lower triangle. Stops at the first pivot that is not
        /// finite and positive.</summary>
        public void Decompose(double[,] a)
        {
            ArgumentNullException.ThrowIfNull(a);
            int n = a.GetLength(0);
            if (a.GetLength(1) != n) throw new ArgumentException("The matrix must be square.", nameof(a));
            if (n > Capacity) throw new ArgumentException($"Order {n} exceeds the capacity {Capacity}.", nameof(a));

            _n = n;
            _valid = false;
            FailedAt = -1;
            LogDet = double.NaN;
            var tally = PivotTally.Start();
            double logDetL = 0.0;
            for (int i = 0; i < n; i++)
            {
                int ri = Row(i);
                for (int j = 0; j < i; j++)
                {
                    int rj = Row(j);
                    double sum = a[i, j];
                    for (int k = 0; k < j; k++) sum -= _l[ri + k] * _l[rj + k];
                    _l[ri + j] = sum / _l[rj + j];
                }
                double aii = a[i, i];
                double d = aii;
                for (int k = 0; k < i; k++) d -= _l[ri + k] * _l[ri + k];
                tally.Observe(d, aii, floored: false);
                if (!(d > 0.0) || !double.IsFinite(d))
                {
                    FailedAt = i;
                    Status = tally.ToStatus();
                    return;
                }
                double lii = Math.Sqrt(d);
                _l[ri + i] = lii;
                _relative[i] = d / aii;
                logDetL += Math.Log(lii);
            }
            _valid = true;
            Status = tally.ToStatus();
            LogDet = 2.0 * logDetL;
            _equilibratedNorm1 = EquilibratedNorm1(a, n);
        }

        // ‖D^(−½)·A·D^(−½)‖₁ from A's lower triangle, recording √A(j,j); called once A is known to be positive definite.
        private double EquilibratedNorm1(double[,] a, int n)
        {
            for (int j = 0; j < n; j++)
            {
                _sqrtDiagonal[j] = Math.Sqrt(a[j, j]);
                _columnSums[j] = 1.0;
            }
            for (int i = 1; i < n; i++)
                for (int j = 0; j < i; j++)
                {
                    double e = Math.Abs(a[i, j]) / (_sqrtDiagonal[i] * _sqrtDiagonal[j]);
                    _columnSums[i] += e;
                    _columnSums[j] += e;
                }
            double norm = 0.0;
            for (int j = 0; j < n; j++) norm = Math.Max(norm, _columnSums[j]);
            return norm;
        }

        /// <inheritdoc/>
        public double EquilibratedCondition()
        {
            RequireValid();
            return ConditionEstimate.Equilibrated(this, _sqrtDiagonal, _equilibratedNorm1);
        }

        /// <summary>L(i, j); 0 above the diagonal.</summary>
        public double Factor(int i, int j)
        {
            RequireValid();
            if ((uint)i >= (uint)_n || (uint)j >= (uint)_n) throw new ArgumentOutOfRangeException(i < 0 || i >= _n ? nameof(i) : nameof(j));
            return j > i ? 0.0 : _l[Row(i) + j];
        }

        /// <summary>The relative pivot <c>d_j/A(j,j)</c> of row <paramref name="j"/>: for a Gram matrix, the squared sine
        /// of the angle between column j and the span of the earlier columns.</summary>
        public double RelativePivot(int j)
        {
            RequireValid();
            if ((uint)j >= (uint)_n) throw new ArgumentOutOfRangeException(nameof(j));
            return _relative[j];
        }

        /// <summary>Solves <c>A·x = b</c>.</summary>
        public double[] Solve(ReadOnlySpan<double> b)
        {
            RequireValid(b.Length);
            double[] x = b.ToArray();
            SolveInPlace(x);
            return x;
        }

        /// <summary>In place: <paramref name="x"/> enters as b and leaves as <c>A⁻¹b</c>.</summary>
        public void SolveInPlace(Span<double> x)
        {
            ForwardSolveInPlace(x);
            BackSolveInPlace(x);
        }

        /// <summary>In place: <c>L·y = b</c>.</summary>
        public void ForwardSolveInPlace(Span<double> x)
        {
            RequireValid(x.Length);
            for (int i = 0; i < _n; i++)
            {
                int ri = Row(i);
                double s = x[i];
                for (int k = 0; k < i; k++) s -= _l[ri + k] * x[k];
                x[i] = s / _l[ri + i];
            }
        }

        /// <summary>In place: <c>Lᵀ·x = y</c>, row-oriented, so each row of the packed factor is read contiguously.</summary>
        public void BackSolveInPlace(Span<double> x)
        {
            RequireValid(x.Length);
            for (int i = _n - 1; i >= 0; i--)
            {
                int ri = Row(i);
                double xi = x[i] / _l[ri + i];
                x[i] = xi;
                for (int k = 0; k < i; k++) x[k] -= _l[ri + k] * xi;
            }
        }

        /// <summary><c>bᵀA⁻¹b = ‖L⁻¹b‖²</c>: a sum of squares, never the difference a materialized inverse invites.</summary>
        public double InverseQuadraticForm(ReadOnlySpan<double> b)
        {
            RequireValid(b.Length);
            Span<double> y = _n <= 256 ? stackalloc double[_n] : new double[_n];
            b.CopyTo(y);
            ForwardSolveInPlace(y);
            double s = 0.0;
            for (int i = 0; i < _n; i++) s += y[i] * y[i];
            return s;
        }

        /// <summary>
        /// Maps standard normals <paramref name="z"/> to <c>v = L⁻ᵀz ~ N(0, A⁻¹)</c>, the precision-factor draw: for
        /// <c>β ~ N(A⁻¹h, σ²A⁻¹)</c>, add σ·v to the solve of h.
        /// </summary>
        public void SampleInnovation(ReadOnlySpan<double> z, Span<double> result)
        {
            RequireValid(z.Length);
            if (result.Length != _n) throw new ArgumentException("The result length must equal the matrix order.", nameof(result));
            z.CopyTo(result);
            BackSolveInPlace(result);
        }

        /// <summary>
        /// Writes <c>A⁻¹ = (L⁻¹)ᵀ·L⁻¹</c> into <paramref name="target"/>, forming L⁻¹ by forward substitution on demand; true
        /// when every entry is finite. A valid factor can have an inverse past <see cref="double.MaxValue"/>: <c>[1e-320]</c>
        /// factors, but its reciprocal overflows. Then this returns false and <paramref name="target"/> holds no inverse. For a
        /// caller that keeps a precision matrix; a solve or <see cref="InverseQuadraticForm"/> needs none.
        /// </summary>
        public bool TryWriteInverseTo(double[,] target)
        {
            ArgumentNullException.ThrowIfNull(target);
            RequireValid();
            int n = _n;
            if (target.GetLength(0) != n || target.GetLength(1) != n)
                throw new ArgumentException($"The target must be {n} × {n}.", nameof(target));
            _inverseFactor ??= new double[_l.Length];
            double[] inv = _inverseFactor;   // L⁻¹, packed lower triangle as L
            for (int i = 0; i < n; i++)
            {
                int ri = Row(i);
                double lii = _l[ri + i];
                inv[ri + i] = 1.0 / lii;
                if (!double.IsFinite(inv[ri + i])) return false;
                for (int j = 0; j < i; j++)
                {
                    double sum = 0.0;
                    for (int k = j; k < i; k++) sum -= _l[ri + k] * inv[Row(k) + j];
                    inv[ri + j] = sum / lii;
                    if (!double.IsFinite(inv[ri + j])) return false;
                }
            }
            for (int i = 0; i < n; i++)
                for (int j = i; j < n; j++)
                {
                    double sum = 0.0;
                    for (int k = j; k < n; k++) sum += inv[Row(k) + i] * inv[Row(k) + j];
                    if (!double.IsFinite(sum)) return false;
                    target[i, j] = sum;
                    target[j, i] = sum;
                }
            return true;
        }

        /// <summary>
        /// Adds <c>L·z</c> to <paramref name="x"/>, each row accumulated into x in column order: with x the mean and z standard
        /// normals, a draw of <c>N(x, A)</c>, the covariance-factor draw (<see cref="SampleInnovation"/> is the precision-factor
        /// one).
        /// </summary>
        public void AddSample(ReadOnlySpan<double> z, Span<double> x)
        {
            RequireValid(z.Length);
            if (x.Length != _n) throw new ArgumentException($"Expected length {_n}, got {x.Length}.", nameof(x));
            for (int i = 0; i < _n; i++)
            {
                int ri = Row(i);
                for (int j = 0; j <= i; j++) x[i] += _l[ri + j] * z[j];
            }
        }

        private static int Row(int i) => i * (i + 1) / 2;

        private void RequireValid()
        {
            if (_valid) return;
            throw new InvalidOperationException(FailedAt >= 0
                ? $"The matrix is not positive definite: the factorization stopped at pivot {FailedAt}."
                : "No matrix has been factored.");
        }

        private void RequireValid(int length)
        {
            RequireValid();
            if (length != _n) throw new ArgumentException($"Expected length {_n}, got {length}.");
        }
    }
}
