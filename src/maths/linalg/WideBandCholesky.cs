using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// Strict Cholesky <c>A = L·Lᵀ</c> of a symmetric positive-definite matrix whose entries vanish beyond a half-bandwidth
    /// <c>p</c> (N-X2): the factor for Grams that are block-banded rather than banded by a spline's degree, such as a tensor
    /// design's in lexicographic order (half-bandwidth about <c>p·ν₂ + p</c>), where <see cref="BandCholesky"/>'s diagonal-major
    /// layout strides across a wide band. Storage is row-packed, <c>L(i, j)</c> for <c>j ∈ [i − p, i]</c> contiguous in row
    /// <c>i</c>, and rows are factored in order (Cholesky–Banachiewicz) so every inner product reads two contiguous runs.
    /// Cost O(n·p²) to factor, O(n·p) per solve, against <see cref="DenseCholesky"/>'s O(n³) and O(n²).
    /// <list type="bullet">
    /// <item><b>The dense factor's conventions.</b> Strict, with no pivot floor: the first pivot that is not finite and
    /// positive stops the factorization, <see cref="FailedAt"/> says where, and the factor refuses to solve. A draw is
    /// <c>L⁻ᵀz</c>. Where the matrix is zero outside the band the two factors do the same operations on the same nonzero
    /// terms in the same order, so the factor, the log-determinant, the solves, the draws and the status are bit-identical.</item>
    /// <item><b>The band is a precondition.</b> Only the lower triangle within the band is read; an entry beyond it is taken
    /// to be zero, and a matrix with one is another matrix.</item>
    /// </list>
    /// <see cref="SpdFactors.ForBand"/> chooses between this and the dense factor by size.
    /// </summary>
    public sealed class WideBandCholesky : ISpdFactor
    {
        private readonly int _p;                  // half-bandwidth: L(i, j) = 0 for j < i − p
        private readonly int _width;              // p + 1 stored entries per row
        private readonly double[] _l;             // row i holds L(i, i − p .. i) at i·(p + 1) + (j − i + p); an early row's first slots are unused
        private readonly double[] _relative;      // d_j / A(j, j) per row
        private readonly double[] _sqrtDiagonal;  // √A(j, j), for the condition estimate
        private readonly double[] _columnSums;    // scratch for ‖D^(−½)·A·D^(−½)‖₁
        private double _equilibratedNorm1 = double.NaN;
        private bool _valid;

        /// <param name="dimension">Matrix order n.</param>
        /// <param name="bandwidth">Half-bandwidth p: the entries <c>A(i, j)</c> with <c>i − j &gt; p</c> are zero. Clamped to n − 1.</param>
        public WideBandCholesky(int dimension, int bandwidth)
        {
            if (dimension <= 0) throw new ArgumentOutOfRangeException(nameof(dimension), "Dimension must be positive.");
            if (bandwidth < 0) throw new ArgumentOutOfRangeException(nameof(bandwidth), "Bandwidth must be non-negative.");
            Dimension = dimension;
            _p = Math.Min(bandwidth, dimension - 1);
            _width = _p + 1;
            if ((long)dimension * _width > Array.MaxLength) throw new ArgumentOutOfRangeException(nameof(bandwidth), "The band is too large to store.");
            _l = new double[dimension * _width];
            _relative = new double[dimension];
            _sqrtDiagonal = new double[dimension];
            _columnSums = new double[dimension];
        }

        /// <summary>Matrix order n.</summary>
        public int Dimension { get; }

        /// <summary>Half-bandwidth p.</summary>
        public int Bandwidth => _p;

        /// <summary>What the last <see cref="Decompose"/> found, up to and including any failing pivot.</summary>
        public FactorStatus Status { get; private set; }

        /// <summary>Row at which the last factorization stopped; −1 when the matrix is positive definite.</summary>
        public int FailedAt { get; private set; } = -1;

        /// <summary>ln|A| = 2·Σ ln L(i,i); NaN unless the last factorization succeeded.</summary>
        public double LogDet { get; private set; } = double.NaN;

        /// <summary>Factors <paramref name="a"/>, reading the lower triangle within the band. Stops at the first pivot that
        /// is not finite and positive.</summary>
        public void Decompose(double[,] a)
        {
            ArgumentNullException.ThrowIfNull(a);
            int n = Dimension;
            if (a.GetLength(0) != n || a.GetLength(1) != n) throw new ArgumentException($"The matrix must be {n} × {n}.", nameof(a));

            _valid = false;
            FailedAt = -1;
            LogDet = double.NaN;
            var tally = PivotTally.Start();
            double logDetL = 0.0;
            for (int i = 0; i < n; i++)
            {
                int ri = Row(i), first = Math.Max(0, i - _p);
                for (int j = first; j < i; j++)
                {
                    int rj = Row(j);
                    double sum = a[i, j];
                    int from = Math.Max(first, Math.Max(0, j - _p));
                    for (int k = from; k < j; k++) sum -= _l[ri + k - i + _p] * _l[rj + k - j + _p];
                    _l[ri + j - i + _p] = sum / _l[rj + _p];
                }
                double aii = a[i, i];
                double d = aii;
                for (int k = first; k < i; k++) d -= _l[ri + k - i + _p] * _l[ri + k - i + _p];
                tally.Observe(d, aii, floored: false);
                if (!(d > 0.0) || !double.IsFinite(d))
                {
                    FailedAt = i;
                    Status = tally.ToStatus();
                    return;
                }
                double lii = Math.Sqrt(d);
                _l[ri + _p] = lii;
                _relative[i] = d / aii;
                logDetL += Math.Log(lii);
            }
            _valid = true;
            Status = tally.ToStatus();
            LogDet = 2.0 * logDetL;
            _equilibratedNorm1 = EquilibratedNorm1(a);
        }

        // ‖D^(−½)·A·D^(−½)‖₁ over the band, recording √A(j,j); called once A is known to be positive definite.
        private double EquilibratedNorm1(double[,] a)
        {
            int n = Dimension;
            for (int j = 0; j < n; j++)
            {
                _sqrtDiagonal[j] = Math.Sqrt(a[j, j]);
                _columnSums[j] = 1.0;
            }
            for (int i = 1; i < n; i++)
                for (int j = Math.Max(0, i - _p); j < i; j++)
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

        /// <summary>L(i, j); 0 outside the band and above the diagonal.</summary>
        public double Factor(int i, int j)
        {
            RequireValid();
            if ((uint)i >= (uint)Dimension || (uint)j >= (uint)Dimension) throw new ArgumentOutOfRangeException(i < 0 || i >= Dimension ? nameof(i) : nameof(j));
            return j > i || i - j > _p ? 0.0 : _l[Row(i) + j - i + _p];
        }

        /// <summary>The relative pivot <c>d_j/A(j,j)</c> of row <paramref name="j"/>.</summary>
        public double RelativePivot(int j)
        {
            RequireValid();
            if ((uint)j >= (uint)Dimension) throw new ArgumentOutOfRangeException(nameof(j));
            return _relative[j];
        }

        /// <inheritdoc/>
        public double[] Solve(ReadOnlySpan<double> b)
        {
            RequireValid(b.Length);
            double[] x = b.ToArray();
            SolveInPlace(x);
            return x;
        }

        /// <inheritdoc/>
        public void SolveInPlace(Span<double> x)
        {
            ForwardSolveInPlace(x);
            BackSolveInPlace(x);
        }

        /// <summary>In place: <c>L·y = b</c>.</summary>
        public void ForwardSolveInPlace(Span<double> x)
        {
            RequireValid(x.Length);
            for (int i = 0; i < Dimension; i++)
            {
                int ri = Row(i), first = Math.Max(0, i - _p);
                double s = x[i];
                for (int k = first; k < i; k++) s -= _l[ri + k - i + _p] * x[k];
                x[i] = s / _l[ri + _p];
            }
        }

        /// <summary>In place: <c>Lᵀ·x = y</c>, row-oriented, so each row of the packed factor is read contiguously.</summary>
        public void BackSolveInPlace(Span<double> x)
        {
            RequireValid(x.Length);
            for (int i = Dimension - 1; i >= 0; i--)
            {
                int ri = Row(i), first = Math.Max(0, i - _p);
                double xi = x[i] / _l[ri + _p];
                x[i] = xi;
                for (int k = first; k < i; k++) x[k] -= _l[ri + k - i + _p] * xi;
            }
        }

        /// <summary><c>bᵀA⁻¹b = ‖L⁻¹b‖²</c>: a sum of squares, never the difference a materialized inverse invites.</summary>
        public double InverseQuadraticForm(ReadOnlySpan<double> b)
        {
            RequireValid(b.Length);
            Span<double> y = Dimension <= 256 ? stackalloc double[Dimension] : new double[Dimension];
            b.CopyTo(y);
            ForwardSolveInPlace(y);
            double s = 0.0;
            for (int i = 0; i < Dimension; i++) s += y[i] * y[i];
            return s;
        }

        /// <summary>Maps standard normals <paramref name="z"/> to <c>v = L⁻ᵀz ~ N(0, A⁻¹)</c>, the precision-factor draw
        /// <see cref="DenseCholesky"/> makes, so one z gives the same v under either factor.</summary>
        public void SampleInnovation(ReadOnlySpan<double> z, Span<double> result)
        {
            RequireValid(z.Length);
            if (result.Length != Dimension) throw new ArgumentException("The result length must equal the matrix order.", nameof(result));
            z.CopyTo(result);
            BackSolveInPlace(result);
        }

        private int Row(int i) => i * _width;

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
            if (length != Dimension) throw new ArgumentException($"Expected length {Dimension}, got {length}.");
        }
    }

    /// <summary>
    /// Chooses the factor of a Gram whose bandwidth its caller declares (N-X2): the wide-band factor when the band saves
    /// both storage and arithmetic, the dense factor otherwise. A factor of order n and half-bandwidth p stores n·(p + 1)
    /// numbers against the dense n(n + 1)/2, and costs about n·(p + 1)² against n³/3, so the band wins on both once
    /// 2(p + 1) ≤ n; below <see cref="MinimumOrder"/> the dense factor's simpler loops win whatever the band. Timed, the band is
    /// already 1.5× faster at 2(p + 1) = n and 16× faster at p + 1 = n/8, so the rule errs toward the dense factor.
    /// </summary>
    public static class SpdFactors
    {
        /// <summary>The smallest order at which a banded Gram takes the wide-band factor.</summary>
        public const int MinimumOrder = 64;

        /// <summary>Whether a Gram of order <paramref name="order"/> and half-bandwidth <paramref name="halfBandwidth"/> takes the
        /// wide-band factor.</summary>
        public static bool UsesWideBand(int order, int halfBandwidth)
            => order >= MinimumOrder && halfBandwidth >= 0 && 2L * (halfBandwidth + 1) <= order;

        /// <summary>
        /// The strict factor of <paramref name="gram"/>, whose entries vanish beyond <paramref name="halfBandwidth"/>: wide-band
        /// by the rule of <see cref="UsesWideBand"/>, else dense. The factor is decomposed; read its
        /// <see cref="ISpdFactor.Status"/> before trusting it.
        /// </summary>
        public static ISpdFactor ForBand(double[,] gram, int halfBandwidth)
        {
            ArgumentNullException.ThrowIfNull(gram);
            int n = gram.GetLength(0);
            if (gram.GetLength(1) != n) throw new ArgumentException("The matrix must be square.", nameof(gram));
            if (UsesWideBand(n, halfBandwidth))
            {
                var band = new WideBandCholesky(n, halfBandwidth);
                band.Decompose(gram);
                return band;
            }
            var dense = new DenseCholesky(n);
            dense.Decompose(gram);
            return dense;
        }
    }
}
