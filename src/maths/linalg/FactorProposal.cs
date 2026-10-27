namespace Maths.LinAlg
{
    /// <summary>How a proposal changes the factored matrix's rows and columns.</summary>
    public enum FactorChange
    {
        /// <summary>A new last row and column.</summary>
        Birth,

        /// <summary>One row and column removed; the later ones move up.</summary>
        Death,

        /// <summary>One row and column removed and a new one added last: a staged death, then a birth.</summary>
        Change,
    }

    /// <summary>
    /// What a proposed change to a <see cref="BorderedCholesky"/> would do, computed beside the accepted factor
    /// without touching it, so the consumer decides admission before anything is committed.
    /// </summary>
    /// <param name="Change">The kind of change.</param>
    /// <param name="Index">The row added (a birth: the current order) or removed (a death or a change).</param>
    /// <param name="PositiveDefinite">The changed matrix is numerically positive definite: the new pivot is finite
    /// and positive. Only such a proposal can be accepted.</param>
    /// <param name="RelativePivot">The new row's pivot over its diagonal entry, d/c: the squared sine of the angle
    /// between the new column and the span of the others. NaN for a death.</param>
    /// <param name="LogDetChange">ln|A′| − ln|A|.</param>
    /// <param name="MinLeaveOneOutPivot">The smallest leave-one-out pivot of the changed matrix
    /// (<see cref="BorderedCholesky.LeaveOneOutPivot"/>): what an order-free admission gate reads.</param>
    public readonly record struct FactorProposal(
        FactorChange Change, int Index, bool PositiveDefinite, double RelativePivot, double LogDetChange,
        double MinLeaveOneOutPivot);

    /// <summary>
    /// How far an incrementally updated factor had moved from a fresh strict factorization of the same matrix.
    /// </summary>
    /// <param name="PositiveDefinite">The strict factorization succeeded; when it did not, the incremental factor is
    /// kept and the other fields are NaN.</param>
    /// <param name="Factor">max |ΔL(i,j)|/√A(i,i): each row of L has norm √A(i,i), so this is scale-free.</param>
    /// <param name="InverseDiagonal">max |Δ(A⁻¹)(j,j)|/(A⁻¹)(j,j), against the fresh values.</param>
    public readonly record struct FactorDrift(bool PositiveDefinite, double Factor, double InverseDiagonal);
}
