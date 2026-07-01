namespace Graphs.Primitives;

/// <summary>
/// How the two directed values of an undirected edge — one from each endpoint's
/// perspective — collapse into the single value a consumer reads. The shared
/// reconciliation policy for per-edge field symmetrization
/// (<see cref="EdgeFieldSymmetrization"/>); the value-level kin of the set-level
/// kNN symmetrization in <c>Graphs.Neighbors</c>. Inert for fields that are
/// symmetric by construction.
/// </summary>
public enum SymmetrizationRule
{
    /// <summary>min — kept only where BOTH directions rank it strong (mutual-kNN
    /// flavored; density-robust, splits more at density boundaries).</summary>
    Mutual,

    /// <summary>max — kept where EITHER direction ranks it strong (chains across
    /// density gradients, merges more, closer to plain single-linkage).</summary>
    Inclusive,

    /// <summary>arithmetic mean of the two directed values — smooths the
    /// asymmetry without committing to AND/OR semantics.</summary>
    Mean,
}
