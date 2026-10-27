using System;

namespace Maths.LinAlg
{
    /// <summary>
    /// A strict dense Cholesky A = L·Lᵀ of a symmetric positive-definite matrix whose rows and columns come and go:
    /// BASS's active-term Gram as terms are born, die and change. It keeps A itself, so it can rebuild, and
    /// diag(A⁻¹), so each column's angle to the span of all the others reads in O(1) (the BASS port's order-free
    /// admission gate, D06). Updating follows Golub and Van Loan (2013, §6.5).
    /// <list type="bullet">
    /// <item><b>Births</b> border the factor with a new last row, computed exactly as <see cref="DenseCholesky"/>
    /// computes a row, so a factor grown one row at a time is bit-identical to the strict one.</item>
    /// <item><b>Deaths</b> drop a row of L and restore the triangle with Givens rotations.</item>
    /// <item><b>Changes</b> are a death and then a birth, staged together; the new column goes last.</item>
    /// <item><b>Proposals are private.</b> Each Propose… call computes its candidate beside the accepted factor and
    /// returns a <see cref="FactorProposal"/>. <see cref="Accept"/> commits it; a new proposal, a rebuild or a
    /// decomposition discards it. A rejected proposal leaves every bit of the accepted factor as it was.</item>
    /// <item><b>Rebuilds</b> refactor the kept matrix strictly, every <see cref="RebuildInterval"/> accepted updates
    /// or on demand, and report the drift they removed as a <see cref="FactorDrift"/>.</item>
    /// </list>
    /// LinAlg reports; which proposals to admit is the consumer's decision. Not thread-safe.
    /// </summary>
    public sealed class BorderedCholesky : ISpdFactor
    {
        // A death's downdate of diag(A⁻¹) that keeps less than this share of an entry has cancelled away the digits
        // that matter, so that entry is recomputed from the new factor.
        private const double CancellationGuard = 1e-3;

        private State _accepted;
        private State _staged;                       // a death or change's candidate; a rebuild's fresh factor
        private readonly double[] _row;              // a birth's new factor row, y = L⁻¹u
        private readonly double[] _inverseColumn;    // a birth's A⁻¹u, or a death's A⁻¹e_m
        private readonly double[] _border;           // a birth's u, kept to store into A
        private readonly double[] _scratch;
        private readonly double[] _cos, _sin;        // a death's rotations
        private readonly int[] _cancelled;
        private readonly double[] _sqrtDiagonal, _columnSums;
        private double _borderDiagonal, _pivot;
        private FactorProposal? _pending;
        private bool _valid = true;

        /// <param name="capacity">Largest matrix order this instance will hold; storage is allocated once.</param>
        /// <param name="rebuildInterval">Accepted updates between strict rebuilds; 0 rebuilds only on demand.</param>
        public BorderedCholesky(int capacity, int rebuildInterval = 100)
        {
            if (capacity < 0 || capacity > 46_340)
                throw new ArgumentOutOfRangeException(nameof(capacity), "Capacity must be in [0, 46340].");
            ArgumentOutOfRangeException.ThrowIfNegative(rebuildInterval);
            Capacity = capacity;
            RebuildInterval = rebuildInterval;
            _accepted = new State(capacity);
            _staged = new State(capacity);
            _row = new double[capacity];
            _inverseColumn = new double[capacity];
            _border = new double[capacity];
            _scratch = new double[capacity + 1];
            _cos = new double[capacity];
            _sin = new double[capacity];
            _cancelled = new int[capacity];
            _sqrtDiagonal = new double[capacity];
            _columnSums = new double[capacity];
            Refresh();
        }

        /// <summary>Largest order this instance can hold.</summary>
        public int Capacity { get; }

        /// <summary>Accepted updates between strict rebuilds; 0 rebuilds only on demand.</summary>
        public int RebuildInterval { get; }

        /// <summary>Order of the accepted matrix.</summary>
        public int Dimension => _accepted.N;

        /// <summary>The accepted factor's pivots, d_j = L(j,j)², against the diagonal of A.</summary>
        public FactorStatus Status { get; private set; }

        /// <summary>ln|A| = 2·Σ ln L(i,i); 0 for the empty matrix, NaN after a failed <see cref="Decompose"/>.</summary>
        public double LogDet { get; private set; }

        /// <summary>Row at which the last <see cref="Decompose"/> stopped; −1 when the factor is valid.</summary>
        public int FailedAt { get; private set; } = -1;

        /// <summary>Updates accepted since the last strict factorization.</summary>
        public int UpdatesSinceRebuild { get; private set; }

        /// <summary>What the most recent rebuild found; null before the first.</summary>
        public FactorDrift? LastDrift { get; private set; }

        /// <summary>The proposal <see cref="Accept"/> would commit, if any.</summary>
        public FactorProposal? Pending => _pending;

        /// <summary>Empties the factor: order 0, valid.</summary>
        public void Clear()
        {
            _pending = null;
            _accepted.N = 0;
            _valid = true;
            FailedAt = -1;
            UpdatesSinceRebuild = 0;
            LastDrift = null;
            Refresh();
        }

        /// <summary>Factors <paramref name="a"/> strictly, reading its lower triangle, and keeps it. Stops at the first
        /// pivot that is not finite and positive, as <see cref="DenseCholesky"/> does.</summary>
        public void Decompose(double[,] a)
        {
            ArgumentNullException.ThrowIfNull(a);
            int n = a.GetLength(0);
            if (a.GetLength(1) != n) throw new ArgumentException("The matrix must be square.", nameof(a));
            if (n > Capacity) throw new ArgumentException($"Order {n} exceeds the capacity {Capacity}.", nameof(a));

            _pending = null;
            UpdatesSinceRebuild = 0;
            LastDrift = null;
            State s = _accepted;
            s.N = n;
            for (int i = 0; i < n; i++)
                for (int j = 0; j <= i; j++) s.A[Row(i) + j] = a[i, j];
            if (!FactorStrict(s, out FactorStatus status, out int failedAt))
            {
                _valid = false;
                FailedAt = failedAt;
                Status = status;
                LogDet = double.NaN;
                return;
            }
            FreshInverseDiagonal(s);
            _valid = true;
            FailedAt = -1;
            Refresh();
        }

        /// <summary>
        /// Proposes a new last row and column: <paramref name="border"/> is its entries against the current columns,
        /// in order, and <paramref name="diagonal"/> its diagonal entry. O(n²).
        /// </summary>
        public FactorProposal ProposeBirth(ReadOnlySpan<double> border, double diagonal)
        {
            RequireValid();
            int n = _accepted.N;
            if (n == Capacity) throw new InvalidOperationException($"The factor is at its capacity {Capacity}.");
            if (border.Length != n) throw new ArgumentException($"Expected a border of length {n}, got {border.Length}.", nameof(border));

            _pending = null;
            border.CopyTo(_border);
            _borderDiagonal = diagonal;
            _pivot = Border(_accepted, border, diagonal);
            return Propose(FactorChange.Birth, n, _pivot, diagonal, Math.Log(_pivot), BirthMinLeaveOneOut(_accepted, _pivot, diagonal));
        }

        /// <summary>Proposes removing row and column <paramref name="index"/>; the later ones move up. Deaths never
        /// lose positive definiteness. O(n²).</summary>
        public FactorProposal ProposeDeath(int index)
        {
            RequireValid();
            RequireIndex(index);
            _pending = null;
            StageDeath(index);
            return Propose(FactorChange.Death, index, double.NaN, double.NaN, StagedLogDet() - LogDet, MinLeaveOneOut(_staged));
        }

        /// <summary>
        /// Proposes replacing row and column <paramref name="index"/>: it is removed and the new one added last, with
        /// <paramref name="border"/> its entries against the remaining columns in their new order. O(n²).
        /// </summary>
        public FactorProposal ProposeChange(int index, ReadOnlySpan<double> border, double diagonal)
        {
            RequireValid();
            RequireIndex(index);
            int n = _accepted.N;
            if (border.Length != n - 1) throw new ArgumentException($"Expected a border of length {n - 1}, got {border.Length}.", nameof(border));

            _pending = null;
            StageDeath(index);
            border.CopyTo(_border);
            _borderDiagonal = diagonal;
            _pivot = Border(_staged, border, diagonal);
            double logDetChange = StagedLogDet() + Math.Log(_pivot) - LogDet;
            return Propose(FactorChange.Change, index, _pivot, diagonal, logDetChange, BirthMinLeaveOneOut(_staged, _pivot, diagonal));
        }

        /// <summary>
        /// bᵀA′⁻¹b for the proposed matrix A′, with <paramref name="b"/> in its order: one entry per column after
        /// the change. A birth and a change add the new row's contribution to a solve with the factor beneath it.
        /// </summary>
        public double ProposedInverseQuadraticForm(ReadOnlySpan<double> b)
        {
            if (_pending is not FactorProposal p) throw new InvalidOperationException("Nothing is proposed.");
            if (!p.PositiveDefinite) throw new InvalidOperationException("The proposal is not positive definite.");
            State beneath = p.Change == FactorChange.Birth ? _accepted : _staged;
            int n = beneath.N;
            int expected = p.Change == FactorChange.Death ? n : n + 1;
            if (b.Length != expected) throw new ArgumentException($"Expected length {expected}, got {b.Length}.", nameof(b));

            Span<double> y = _scratch.AsSpan(0, n);
            b[..n].CopyTo(y);
            ForwardSolve(beneath, y);
            double q = 0.0;
            for (int i = 0; i < n; i++) q += y[i] * y[i];
            if (p.Change == FactorChange.Death) return q;
            double t = b[n];
            for (int k = 0; k < n; k++) t -= _row[k] * y[k];
            return q + t * t / _pivot;
        }

        /// <summary>Commits the pending proposal. Every <see cref="RebuildInterval"/> accepted updates, a strict
        /// rebuild follows.</summary>
        public void Accept()
        {
            if (_pending is not FactorProposal p) throw new InvalidOperationException("Nothing is proposed.");
            if (!p.PositiveDefinite) throw new InvalidOperationException("The proposal is not positive definite.");
            switch (p.Change)
            {
                case FactorChange.Birth:
                    Append(_accepted);
                    break;
                case FactorChange.Death:
                    Swap();
                    break;
                case FactorChange.Change:
                    Append(_staged);
                    Swap();
                    break;
            }
            _pending = null;
            UpdatesSinceRebuild++;
            if (RebuildInterval > 0 && UpdatesSinceRebuild >= RebuildInterval) Rebuild();
            else Refresh();
        }

        /// <summary>
        /// Refactors the kept matrix strictly and recomputes diag(A⁻¹), replacing the incremental factor when the
        /// matrix is still numerically positive definite and keeping it otherwise. Returns, and records as
        /// <see cref="LastDrift"/>, how far the incremental factor had drifted. O(n³).
        /// </summary>
        public FactorDrift Rebuild()
        {
            RequireValid();
            _pending = null;
            UpdatesSinceRebuild = 0;
            _staged.CopyFrom(_accepted);
            if (!FactorStrict(_staged, out _, out _))
            {
                LastDrift = new FactorDrift(false, double.NaN, double.NaN);
                return LastDrift.Value;
            }
            FreshInverseDiagonal(_staged);

            int n = _accepted.N;
            double factorDrift = 0.0, inverseDrift = 0.0;
            for (int i = 0; i < n; i++)
            {
                int ri = Row(i);
                double scale = Math.Sqrt(_accepted.A[ri + i]);
                for (int j = 0; j <= i; j++)
                    factorDrift = Math.Max(factorDrift, Math.Abs(_accepted.L[ri + j] - _staged.L[ri + j]) / scale);
                double fresh = _staged.InverseDiagonal[i];
                inverseDrift = Math.Max(inverseDrift, Math.Abs(_accepted.InverseDiagonal[i] - fresh) / fresh);
            }
            Swap();
            Refresh();
            LastDrift = new FactorDrift(true, factorDrift, inverseDrift);
            return LastDrift.Value;
        }

        /// <summary>A(i, j) of the accepted matrix.</summary>
        public double Entry(int i, int j)
        {
            RequireValid();
            RequireIndex(i);
            RequireIndex(j);
            return i >= j ? _accepted.A[Row(i) + j] : _accepted.A[Row(j) + i];
        }

        /// <summary>L(i, j); 0 above the diagonal.</summary>
        public double Factor(int i, int j)
        {
            RequireValid();
            RequireIndex(i);
            RequireIndex(j);
            return j > i ? 0.0 : _accepted.L[Row(i) + j];
        }

        /// <summary>The relative pivot <c>d_j/A(j,j)</c> of row <paramref name="j"/>: for a Gram matrix, the squared sine
        /// of the angle between column j and the span of the earlier columns.</summary>
        public double RelativePivot(int j)
        {
            RequireValid();
            RequireIndex(j);
            int rj = Row(j);
            double l = _accepted.L[rj + j];
            return l * l / _accepted.A[rj + j];
        }

        /// <summary>(A⁻¹)(j, j), kept through every update.</summary>
        public double InverseDiagonal(int j)
        {
            RequireValid();
            RequireIndex(j);
            return _accepted.InverseDiagonal[j];
        }

        /// <summary>
        /// 1/(A(j,j)·(A⁻¹)(j,j)): the relative pivot column j would have if it were factored last. For a Gram matrix
        /// it is the squared sine of the angle between column j and the span of all the others, the reciprocal of
        /// column j's variance inflation factor, and it does not depend on column order.
        /// </summary>
        public double LeaveOneOutPivot(int j)
        {
            RequireValid();
            RequireIndex(j);
            return 1.0 / (_accepted.A[Row(j) + j] * _accepted.InverseDiagonal[j]);
        }

        /// <summary>The smallest <see cref="LeaveOneOutPivot"/>; +∞ for the empty matrix.</summary>
        public double MinLeaveOneOutPivot()
        {
            RequireValid();
            return MinLeaveOneOut(_accepted);
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
            RequireValid(x.Length);
            ForwardSolve(_accepted, x);
            BackSolve(_accepted, x);
        }

        /// <summary>In place: <c>L·y = b</c>.</summary>
        public void ForwardSolveInPlace(Span<double> x)
        {
            RequireValid(x.Length);
            ForwardSolve(_accepted, x);
        }

        /// <summary>In place: <c>Lᵀ·x = y</c>.</summary>
        public void BackSolveInPlace(Span<double> x)
        {
            RequireValid(x.Length);
            BackSolve(_accepted, x);
        }

        /// <summary><c>bᵀA⁻¹b = ‖L⁻¹b‖²</c>.</summary>
        public double InverseQuadraticForm(ReadOnlySpan<double> b)
        {
            RequireValid(b.Length);
            int n = _accepted.N;
            Span<double> y = _scratch.AsSpan(0, n);
            b.CopyTo(y);
            ForwardSolve(_accepted, y);
            double s = 0.0;
            for (int i = 0; i < n; i++) s += y[i] * y[i];
            return s;
        }

        /// <summary>Maps standard normals <paramref name="z"/> to <c>v = L⁻ᵀz ~ N(0, A⁻¹)</c>.</summary>
        public void SampleInnovation(ReadOnlySpan<double> z, Span<double> result)
        {
            RequireValid(z.Length);
            if (result.Length != _accepted.N) throw new ArgumentException("The result length must equal the matrix order.", nameof(result));
            z.CopyTo(result);
            BackSolve(_accepted, result);
        }

        /// <inheritdoc/>
        public double EquilibratedCondition()
        {
            RequireValid();
            int n = _accepted.N;
            if (n == 0) return 1.0;
            double[] a = _accepted.A;
            for (int j = 0; j < n; j++)
            {
                _sqrtDiagonal[j] = Math.Sqrt(a[Row(j) + j]);
                _columnSums[j] = 1.0;
            }
            for (int i = 1; i < n; i++)
                for (int j = 0; j < i; j++)
                {
                    double e = Math.Abs(a[Row(i) + j]) / (_sqrtDiagonal[i] * _sqrtDiagonal[j]);
                    _columnSums[i] += e;
                    _columnSums[j] += e;
                }
            double norm = 0.0;
            for (int j = 0; j < n; j++) norm = Math.Max(norm, _columnSums[j]);
            return ConditionEstimate.Equilibrated(this, _sqrtDiagonal, norm);
        }

        // y = L⁻¹u row by row, in DenseCholesky's order, so a birth's row is the row the strict factor would compute;
        // returns the pivot c − ‖y‖² and leaves A⁻¹u = L⁻ᵀy in _inverseColumn.
        private double Border(State s, ReadOnlySpan<double> u, double c)
        {
            int n = s.N;
            double[] l = s.L;
            for (int j = 0; j < n; j++)
            {
                int rj = Row(j);
                double sum = u[j];
                for (int k = 0; k < j; k++) sum -= _row[k] * l[rj + k];
                _row[j] = sum / l[rj + j];
            }
            double d = c;
            for (int k = 0; k < n; k++) d -= _row[k] * _row[k];
            Span<double> w = _inverseColumn.AsSpan(0, n);
            _row.AsSpan(0, n).CopyTo(w);
            BackSolve(s, w);
            return d;
        }

        // Appends the bordered row computed by Border on s, and adds its rank-one term to diag(A⁻¹): the bordered
        // inverse gains w·wᵀ/d in its old block and 1/d in its new corner.
        private void Append(State s)
        {
            int n = s.N;
            int rn = Row(n);
            _row.AsSpan(0, n).CopyTo(s.L.AsSpan(rn, n));
            s.L[rn + n] = Math.Sqrt(_pivot);
            _border.AsSpan(0, n).CopyTo(s.A.AsSpan(rn, n));
            s.A[rn + n] = _borderDiagonal;
            for (int j = 0; j < n; j++) s.InverseDiagonal[j] += _inverseColumn[j] * _inverseColumn[j] / _pivot;
            s.InverseDiagonal[n] = 1.0 / _pivot;
            s.N = n + 1;
        }

        private double BirthMinLeaveOneOut(State s, double d, double c)
        {
            if (!IsPivot(d)) return d / c;
            double min = d / c;
            for (int j = 0; j < s.N; j++)
            {
                double inverse = s.InverseDiagonal[j] + _inverseColumn[j] * _inverseColumn[j] / d;
                min = Math.Min(min, 1.0 / (s.A[Row(j) + j] * inverse));
            }
            return min;
        }

        // The accepted factor with row and column m removed, into _staged.
        private void StageDeath(int m)
        {
            int n = _accepted.N;
            Span<double> v = _inverseColumn.AsSpan(0, n);
            v.Clear();
            v[m] = 1.0;
            ForwardSolveFrom(_accepted, v, m);
            BackSolve(_accepted, v);                 // v = A⁻¹e_m

            _staged.CopyFrom(_accepted);
            DeleteFactorRow(_staged, m);
            DeleteMatrixRow(_staged, m);

            // The inverse of A without row and column m is (A⁻¹) without them, less v·vᵀ/v_m.
            double vm = v[m];
            int count = 0, w = 0;
            double[] inverse = _staged.InverseDiagonal;
            for (int j = 0; j < n; j++)
            {
                if (j == m) continue;
                double old = inverse[j];
                double updated = old - v[j] * v[j] / vm;
                inverse[w] = updated;
                if (!(updated > CancellationGuard * old)) _cancelled[count++] = w;
                w++;
            }
            _staged.N = n - 1;
            for (int t = 0; t < count; t++) inverse[_cancelled[t]] = ExactInverseDiagonal(_staged, _cancelled[t]);
        }

        // Deletes row m of L, leaving rows below it one column past their diagonal, and rotates columns (j, j+1),
        // j = m, m+1, …, to clear that entry. Each later row takes every earlier rotation in turn, and the rotation
        // that clears it, before it moves up into its new packed place.
        private void DeleteFactorRow(State s, int m)
        {
            int n = s.N;
            double[] l = s.L;
            for (int r = m; r < n - 1; r++)
            {
                Span<double> x = _scratch.AsSpan(0, r + 2);
                l.AsSpan(Row(r + 1), r + 2).CopyTo(x);
                for (int j = m; j < r; j++)
                {
                    double a = x[j], b = x[j + 1];
                    x[j] = _cos[j] * a + _sin[j] * b;
                    x[j + 1] = _cos[j] * b - _sin[j] * a;
                }
                double rho = double.Hypot(x[r], x[r + 1]);   // x[r + 1] is the old diagonal, so rho > 0
                _cos[r] = x[r] / rho;
                _sin[r] = x[r + 1] / rho;
                x[r] = rho;
                x[..(r + 1)].CopyTo(l.AsSpan(Row(r), r + 1));
            }
        }

        // Removes row and column m from the packed lower triangle of A.
        private static void DeleteMatrixRow(State s, int m)
        {
            double[] a = s.A;
            for (int i = m + 1; i < s.N; i++)
            {
                int source = Row(i), target = Row(i - 1);
                Array.Copy(a, source, a, target, m);
                Array.Copy(a, source + m + 1, a, target + m, i - m);
            }
        }

        // (A⁻¹)(j, j) = ‖L⁻¹e_j‖², by a forward solve that starts at row j.
        private double ExactInverseDiagonal(State s, int j)
        {
            int n = s.N;
            double[] l = s.L;
            double sum = 0.0;
            for (int i = j; i < n; i++)
            {
                int ri = Row(i);
                double t = i == j ? 1.0 : 0.0;
                for (int k = j; k < i; k++) t -= l[ri + k] * _scratch[k];
                _scratch[i] = t / l[ri + i];
                sum += _scratch[i] * _scratch[i];
            }
            return sum;
        }

        private void FreshInverseDiagonal(State s)
        {
            for (int j = 0; j < s.N; j++) s.InverseDiagonal[j] = ExactInverseDiagonal(s, j);
        }

        // DenseCholesky's strict factorization, reading the kept packed matrix.
        private static bool FactorStrict(State s, out FactorStatus status, out int failedAt)
        {
            int n = s.N;
            double[] l = s.L, a = s.A;
            var tally = PivotTally.Start();
            for (int i = 0; i < n; i++)
            {
                int ri = Row(i);
                for (int j = 0; j < i; j++)
                {
                    int rj = Row(j);
                    double sum = a[ri + j];
                    for (int k = 0; k < j; k++) sum -= l[ri + k] * l[rj + k];
                    l[ri + j] = sum / l[rj + j];
                }
                double aii = a[ri + i];
                double d = aii;
                for (int k = 0; k < i; k++) d -= l[ri + k] * l[ri + k];
                tally.Observe(d, aii, floored: false);
                if (!IsPivot(d))
                {
                    status = tally.ToStatus();
                    failedAt = i;
                    return false;
                }
                l[ri + i] = Math.Sqrt(d);
            }
            status = tally.ToStatus();
            failedAt = -1;
            return true;
        }

        private static void ForwardSolve(State s, Span<double> x) => ForwardSolveFrom(s, x, 0);

        // L·y = b, for a b that is zero above row 'from'.
        private static void ForwardSolveFrom(State s, Span<double> x, int from)
        {
            double[] l = s.L;
            for (int i = from; i < s.N; i++)
            {
                int ri = Row(i);
                double t = x[i];
                for (int k = from; k < i; k++) t -= l[ri + k] * x[k];
                x[i] = t / l[ri + i];
            }
        }

        private static void BackSolve(State s, Span<double> x)
        {
            double[] l = s.L;
            for (int i = s.N - 1; i >= 0; i--)
            {
                int ri = Row(i);
                double xi = x[i] / l[ri + i];
                x[i] = xi;
                for (int k = 0; k < i; k++) x[k] -= l[ri + k] * xi;
            }
        }

        private static double MinLeaveOneOut(State s)
        {
            double min = double.PositiveInfinity;
            for (int j = 0; j < s.N; j++) min = Math.Min(min, 1.0 / (s.A[Row(j) + j] * s.InverseDiagonal[j]));
            return min;
        }

        private double StagedLogDet()
        {
            double sum = 0.0;
            for (int j = 0; j < _staged.N; j++) sum += Math.Log(_staged.L[Row(j) + j]);
            return 2.0 * sum;
        }

        private FactorProposal Propose(FactorChange change, int index, double pivot, double diagonal, double logDetChange, double minLeaveOneOut)
        {
            bool positiveDefinite = change == FactorChange.Death || IsPivot(pivot);
            double relative = change == FactorChange.Death ? double.NaN : diagonal > 0.0 ? pivot / diagonal : 0.0;
            var p = new FactorProposal(change, index, positiveDefinite, relative, logDetChange, minLeaveOneOut);
            _pending = p;
            return p;
        }

        // The accepted factor's pivots and log-determinant, after any change.
        private void Refresh()
        {
            var tally = PivotTally.Start();
            double logDetL = 0.0;
            for (int j = 0; j < _accepted.N; j++)
            {
                int rj = Row(j);
                double l = _accepted.L[rj + j];
                tally.Observe(l * l, _accepted.A[rj + j], floored: false);
                logDetL += Math.Log(l);
            }
            Status = tally.ToStatus();
            LogDet = 2.0 * logDetL;
        }

        private void Swap() => (_accepted, _staged) = (_staged, _accepted);

        private static bool IsPivot(double d) => d > 0.0 && double.IsFinite(d);

        private static int Row(int i) => i * (i + 1) / 2;

        private void RequireValid()
        {
            if (_valid) return;
            throw new InvalidOperationException($"The matrix is not positive definite: the factorization stopped at pivot {FailedAt}.");
        }

        private void RequireValid(int length)
        {
            RequireValid();
            if (length != _accepted.N) throw new ArgumentException($"Expected length {_accepted.N}, got {length}.");
        }

        private void RequireIndex(int i)
        {
            if ((uint)i >= (uint)_accepted.N) throw new ArgumentOutOfRangeException(nameof(i), $"Index {i} is outside [0, {_accepted.N}).");
        }

        // One complete factor: L and A in packed lower-triangular rows, and diag(A⁻¹).
        private sealed class State
        {
            public readonly double[] L;
            public readonly double[] A;
            public readonly double[] InverseDiagonal;
            public int N;

            public State(int capacity)
            {
                L = new double[capacity * (capacity + 1) / 2];
                A = new double[capacity * (capacity + 1) / 2];
                InverseDiagonal = new double[capacity];
            }

            public void CopyFrom(State other)
            {
                int length = Row(other.N);
                Array.Copy(other.L, L, length);
                Array.Copy(other.A, A, length);
                Array.Copy(other.InverseDiagonal, InverseDiagonal, other.N);
                N = other.N;
            }
        }
    }
}
