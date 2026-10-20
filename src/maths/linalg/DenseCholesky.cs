using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// Strict dense Cholesky <c>A = L·Lᵀ</c> of a symmetric positive-definite matrix: the reference factor for dense
    /// Gram matrices (BASS's active-term Gram), and the full rebuild that an incremental update is checked against.
    /// <list type="bullet">
    /// <item><b>Strict.</b> There is no pivot floor. The first pivot that is not finite and positive stops the
    /// factorization, <see cref="FailedAt"/> says where, and the factor refuses to solve. LinAlg reports; whether a
    /// singular design is admissible under some other prior is the caller's decision.</item>
    /// <item><b>Solve-oriented.</b> No inverse is formed. Solves, <c>bᵀA⁻¹b</c>, the log-determinant and
    /// <c>N(0, A⁻¹)</c> draws all go through the two triangular solves.</item>
    /// <item><b>Row at a time.</b> Rows are factored in order (Cholesky–Banachiewicz) into packed row-major storage,
    /// each row the bordered extension of the factor above it: the step a birth reuses.</item>
    /// </list>
    /// Only the lower triangle of the input is read. With no pivot near the floor, its factor is bit-identical to
    /// <see cref="CholeskyDecomposition"/>'s, whose floored, inverse-forming path its existing callers keep.
    /// </summary>
    public sealed class DenseCholesky
    {
        private readonly double[] _l;          // packed lower triangle: L(i, j) at i(i+1)/2 + j
        private readonly double[] _relative;   // d_j / A(j,j) per row
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
