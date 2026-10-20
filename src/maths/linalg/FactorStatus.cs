namespace Maths.LinAlg
{
    /// <summary>
    /// What a symmetric factorization found about its input, reported beside the factor so the consumer decides
    /// admission. The pivots are the Schur-complement diagonals d_j before any flooring: D(j) for L·D·Lᵀ, L(j,j)²
    /// for L·Lᵀ. For a symmetric positive-definite matrix every d_j lies in [λ_min, λ_max], so
    /// <see cref="PivotRatio"/> ≥ 1/κ(A): a small ratio proves ill-conditioning, but a ratio near 1 does not prove
    /// good conditioning. <see cref="MinRelativePivot"/> is the scale-free counterpart.
    /// </summary>
    /// <param name="Finite">Every pivot was finite (a NaN or infinite input entry makes a later pivot non-finite).</param>
    /// <param name="FlooredPivots">Pivots at or below the factorization's absolute floor, which were clamped to it.</param>
    /// <param name="MinPivot">The smallest pivot before flooring; ≤ 0 means the matrix is not numerically positive definite.</param>
    /// <param name="MaxPivot">The largest pivot before flooring.</param>
    /// <param name="MinRelativePivot">
    /// The smallest d_j/A_jj: the pivots of the equilibrated matrix D^(−½)·A·D^(−½), D = diag(A). Each is at most 1,
    /// and for a Gram matrix it is the squared sine of the angle between column j and the span of the earlier
    /// columns. It is at least λ_min of the equilibrated matrix, whose λ_max is at least 1, so 1/MinRelativePivot is a
    /// lower bound on that matrix's condition number. Unlike <see cref="PivotRatio"/> it does not change when columns
    /// are rescaled; it does depend on column order. ≤ 0 means not positive definite.
    /// </param>
    public readonly record struct FactorStatus(bool Finite, int FlooredPivots, double MinPivot, double MaxPivot, double MinRelativePivot)
    {
        /// <summary>Smallest over largest pivot before flooring. For a positive-definite matrix it is at least 1/κ(A),
        /// so 1/PivotRatio is a lower bound on the condition number.</summary>
        public double PivotRatio => MinPivot / MaxPivot;

        /// <summary>Every pivot was finite and strictly positive: the input is numerically positive definite. The
        /// factor may still have clamped pivots below the absolute floor (<see cref="FlooredPivots"/>).</summary>
        public bool PositiveDefinite => Finite && MinPivot > 0.0;
    }

    // Accumulates a factorization's pivots, before flooring, into a FactorStatus.
    internal struct PivotTally
    {
        private bool _nonFinite;
        private int _floored;
        private double _min, _max, _minRelative;

        public static PivotTally Start() => new()
        {
            _min = double.PositiveInfinity, _max = double.NegativeInfinity, _minRelative = double.PositiveInfinity,
        };

        // pivot: d_j before flooring; diagonal: A(j,j). A non-positive diagonal already makes d_j ≤ 0.
        public void Observe(double pivot, double diagonal, bool floored)
        {
            if (floored) _floored++;
            if (!double.IsFinite(pivot)) { _nonFinite = true; return; }
            if (pivot < _min) _min = pivot;
            if (pivot > _max) _max = pivot;
            double relative = diagonal > 0.0 ? pivot / diagonal : 0.0;
            if (relative < _minRelative) _minRelative = relative;
        }

        public readonly FactorStatus ToStatus() => new(!_nonFinite, _floored, _min, _max, _minRelative);
    }
}
