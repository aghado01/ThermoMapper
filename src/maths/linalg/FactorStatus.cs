namespace Maths.LinAlg
{
    /// <summary>
    /// What a symmetric factorization found about its input, reported beside the factor so the consumer decides
    /// admission. The pivots are the Schur-complement diagonals d_j before any flooring: D(j) for L·D·Lᵀ, L(j,j)²
    /// for L·Lᵀ. For a symmetric positive-definite matrix every d_j lies in [λ_min, λ_max], so
    /// <see cref="PivotRatio"/> ≥ 1/κ(A): a small ratio proves ill-conditioning, but a ratio near 1 does not prove
    /// good conditioning.
    /// </summary>
    /// <param name="Finite">Every pivot was finite (a NaN or infinite input entry makes a later pivot non-finite).</param>
    /// <param name="FlooredPivots">Pivots at or below the factorization's absolute floor, which were clamped to it.</param>
    /// <param name="MinPivot">The smallest pivot before flooring; ≤ 0 means the matrix is not numerically positive definite.</param>
    /// <param name="MaxPivot">The largest pivot before flooring.</param>
    public readonly record struct FactorStatus(bool Finite, int FlooredPivots, double MinPivot, double MaxPivot)
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
        private double _min, _max;

        public static PivotTally Start() => new() { _min = double.PositiveInfinity, _max = double.NegativeInfinity };

        public void Observe(double pivot, bool floored)
        {
            if (floored) _floored++;
            if (!double.IsFinite(pivot)) { _nonFinite = true; return; }
            if (pivot < _min) _min = pivot;
            if (pivot > _max) _max = pivot;
        }

        public readonly FactorStatus ToStatus() => new(!_nonFinite, _floored, _min, _max);
    }
}
