using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// Cholesky with diagonal pivoting of a symmetric positive semidefinite matrix (Harbrecht, Peters and Schneider 2012),
    /// for drawing from N(0, A) when A is singular to working precision, as a smooth kernel's covariance is.
    /// <list type="bullet">
    /// <item><b>The factor.</b> Each step takes the largest remaining diagonal; the factor stops where every remaining
    /// diagonal is at most <c>relativeTolerance</c> times A's largest diagonal. A ≈ L·Lᵀ with L of <see cref="Rank"/>
    /// columns, its rows in A's own order (<see cref="Loading"/>).</item>
    /// <item><b>The remainder.</b> A semidefinite remainder whose diagonals are at most the stop has every entry within it.
    /// One with an entry beyond <c>indefiniteTolerance</c> of the largest diagonal, a diagonal below it or an off-diagonal
    /// past it, is materially indefinite: the factor is refused (<see cref="PositiveSemidefinite"/>), as is one with an entry
    /// that is not finite. <c>[[1, 2], [2, 1]]</c> leaves −3 on the diagonal; <c>[[0, 1], [1, 0]]</c> leaves 1 off it. The
    /// default, √ε, sits above the rounding a covariance formed in double carries: BGP's posterior covariances leave
    /// remainders to 1.2·10⁻¹¹ of their largest variance, past the 10⁻¹² stop, and those are kept.</item>
    /// <item><b>Sampling first.</b> <see cref="Map"/> takes Rank standard normals to a draw of N(0, L·Lᵀ). A solve on a
    /// singular A, which must refuse or project a right-hand side outside L's range, waits for a caller that needs one.</item>
    /// </list>
    /// Only the lower triangle of the input is read.
    /// </summary>
    public sealed class PivotedCholesky
    {
        private readonly double[] _l;           // n × n row-major at stride n, in A's row order; columns ≥ Rank are 0
        private readonly double[] _remaining;
        private readonly int[] _order;
        private int _n;
        private bool _valid, _factored;

        /// <param name="capacity">Largest matrix order this instance will factor; storage is allocated once.</param>
        public PivotedCholesky(int capacity)
        {
            if (capacity < 0 || capacity > 46_340)
                throw new ArgumentOutOfRangeException(nameof(capacity), "Capacity must be in [0, 46340].");
            Capacity = capacity;
            _l = new double[capacity * capacity];
            _remaining = new double[capacity];
            _order = new int[capacity];
        }

        /// <summary>Largest order this instance can factor.</summary>
        public int Capacity { get; }

        /// <summary>Order of the last factored matrix.</summary>
        public int Dimension => _n;

        /// <summary>The columns kept: the factor's rank at the tolerance.</summary>
        public int Rank { get; private set; }

        /// <summary>The absolute tolerance the last factorization stopped at: the relative tolerance times A's largest
        /// diagonal.</summary>
        public double Tolerance { get; private set; }

        /// <summary>The largest |entry| of the remainder A − L·Lᵀ over A's largest diagonal; NaN when an entry is not finite.</summary>
        public double RelativeRemainder { get; private set; } = double.NaN;

        /// <summary>The remainder lies within the tolerance and every entry is finite: the factor holds, and is read.</summary>
        public bool PositiveSemidefinite => _valid;

        /// <summary>
        /// Factors <paramref name="lower"/> to <paramref name="relativeTolerance"/> of its largest diagonal, then checks the
        /// remainder against <paramref name="indefiniteTolerance"/> of it. Never throws on the matrix's values;
        /// <see cref="PositiveSemidefinite"/> says whether the factor holds.
        /// </summary>
        public void Decompose(double[,] lower, double relativeTolerance = 1e-12, double indefiniteTolerance = DefaultIndefiniteTolerance)
        {
            ArgumentNullException.ThrowIfNull(lower);
            int p = lower.GetLength(0);
            if (lower.GetLength(1) != p) throw new ArgumentException("The matrix must be square.", nameof(lower));
            if (p > Capacity) throw new ArgumentException($"Order {p} exceeds the capacity {Capacity}.", nameof(lower));
            if (!(relativeTolerance >= 0.0 && relativeTolerance < 1.0))
                throw new ArgumentOutOfRangeException(nameof(relativeTolerance), "The relative tolerance lies in [0, 1).");
            if (!(indefiniteTolerance >= 0.0 && indefiniteTolerance < 1.0))
                throw new ArgumentOutOfRangeException(nameof(indefiniteTolerance), "The indefinite tolerance lies in [0, 1).");
            _n = p;
            _factored = true;
            _valid = false;
            RelativeRemainder = double.NaN;
            Array.Clear(_l, 0, p * p);
            double[] factor = _l, remaining = _remaining;
            int[] order = _order;
            double largest = 0.0;
            for (int a = 0; a < p; a++)
            {
                remaining[a] = lower[a, a];
                order[a] = a;
                largest = Math.Max(largest, remaining[a]);
            }
            double tolerance = relativeTolerance * largest;
            int k = 0;
            for (; k < p; k++)
            {
                int best = k;
                for (int i = k + 1; i < p; i++) if (remaining[order[i]] > remaining[order[best]]) best = i;
                double pivot = remaining[order[best]];
                if (!(pivot > tolerance)) break;
                (order[k], order[best]) = (order[best], order[k]);
                int at = order[k];
                double root = Math.Sqrt(pivot);
                factor[at * p + k] = root;
                for (int i = k + 1; i < p; i++)
                {
                    int q = order[i];
                    double s = q >= at ? lower[q, at] : lower[at, q];
                    for (int c = 0; c < k; c++) s -= factor[q * p + c] * factor[at * p + c];
                    double v = s / root;
                    factor[q * p + k] = v;
                    remaining[q] -= v * v;
                }
            }
            Rank = k;
            Tolerance = tolerance;
            _valid = Remainder(lower, largest, Math.Max(tolerance, indefiniteTolerance * largest));
        }

        private const double Unit = 2.220446049250313e-16;   // 2^-52

        /// <summary>The default indefinite tolerance, √ε of the largest diagonal.</summary>
        public const double DefaultIndefiniteTolerance = 1.4901161193847656e-8;

        // The remainder over the rows left out, each entry against the indefinite tolerance widened by the factor's rounding,
        // which is within 2(k + 1)·ε·largest by Cauchy–Schwarz on the loadings. Records the largest |entry| relative to the
        // largest diagonal.
        private bool Remainder(double[,] lower, double largest, double tolerance)
        {
            int p = _n, k = Rank;
            double bound = tolerance + 2.0 * (k + 1) * Unit * largest, worst = 0.0;
            for (int a = 0; a < p; a++)
                for (int c = 0; c < k; c++)
                    if (!double.IsFinite(_l[a * p + c])) return NotFinite();
            for (int i = k; i < p; i++)
            {
                int qi = _order[i];
                double d = _remaining[qi];
                if (!double.IsFinite(d)) return NotFinite();
                worst = Math.Max(worst, Math.Abs(d));
                for (int j = k; j < i; j++)
                {
                    int qj = _order[j];
                    double s = qi >= qj ? lower[qi, qj] : lower[qj, qi];
                    for (int c = 0; c < k; c++) s -= _l[qi * p + c] * _l[qj * p + c];
                    if (!double.IsFinite(s)) return NotFinite();
                    worst = Math.Max(worst, Math.Abs(s));
                }
            }
            RelativeRemainder = largest > 0.0 ? worst / largest : worst;
            return worst <= bound;

            bool NotFinite()
            {
                RelativeRemainder = double.NaN;
                return false;
            }
        }

        /// <summary>L(row, column), the row's loading on a kept column, in A's own row order.</summary>
        public double Loading(int row, int column)
        {
            RequireValid();
            if ((uint)row >= (uint)_n) throw new ArgumentOutOfRangeException(nameof(row));
            if ((uint)column >= (uint)Rank) throw new ArgumentOutOfRangeException(nameof(column));
            return _l[row * _n + column];
        }

        /// <summary>
        /// <paramref name="result"/>[a] = Σ_k L(a, k)·<paramref name="z"/>[k] over the <see cref="Rank"/> kept columns, each row
        /// summed from 0 in column order: with z standard normals, a draw of N(0, L·Lᵀ).
        /// </summary>
        public void Map(ReadOnlySpan<double> z, Span<double> result)
        {
            RequireValid();
            if (z.Length != Rank) throw new ArgumentException($"Expected {Rank} standard normals, got {z.Length}.", nameof(z));
            if (result.Length != _n) throw new ArgumentException($"Expected length {_n}, got {result.Length}.", nameof(result));
            for (int a = 0; a < _n; a++)
            {
                double s = 0.0;
                int row = a * _n;
                for (int c = 0; c < z.Length; c++) s += _l[row + c] * z[c];
                result[a] = s;
            }
        }

        private void RequireValid()
        {
            if (_valid) return;
            throw new InvalidOperationException(!_factored ? "No matrix has been factored."
                : double.IsNaN(RelativeRemainder) ? "The factor holds an entry that is not finite."
                : $"The matrix is materially indefinite: its remainder reaches {RelativeRemainder:E3} of its largest diagonal.");
        }
    }
}
