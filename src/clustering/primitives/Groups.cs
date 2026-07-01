using System;

namespace Clustering.Primitives;

/// <summary>
/// The soft, pre-resolution clustering structure: an <c>N×K</c> membership field where
/// <c>M[i,k] ≥ 0</c> is the strength of point <c>i</c>'s membership in group <c>k</c>. The lingua franca
/// that GMM (responsibilities), HDBSCAN (condensed-tree memberships), and SPC (basins over the alignment
/// landscape) all reduce to before a resolution strategy collapses it to an <see cref="Assignment"/>.
/// A single "group" is a column; the collection is this matrix.
/// </summary>
/// <remarks>
/// Stored flat row-major (<c>M[i*K + k]</c>) for HPC — contiguous and SIMD/BLAS-friendly; the resolution
/// ops (argmax, threshold, row-normalize, row-entropy) are row/column reductions. Soft membership and
/// overlap (a point strong in several groups) are native; the crisp <see cref="Assignment"/> is a
/// per-row reduction. The buffer is <b>aliased</b>, not copied — the caller must not mutate it afterward.
/// </remarks>
public sealed class Groups
{
    private readonly double[] _m;   // row-major, length PointCount * GroupCount

    /// <summary>Number of points (rows), <c>N</c>.</summary>
    public int PointCount { get; }

    /// <summary>Number of groups (columns), <c>K</c>.</summary>
    public int GroupCount { get; }

    /// <summary>
    /// Wraps a flat row-major membership buffer of length <paramref name="pointCount"/> ·
    /// <paramref name="groupCount"/> (<c>M[i*K + k]</c>). The buffer is aliased, not copied.
    /// </summary>
    public Groups(double[] membership, int pointCount, int groupCount)
    {
        ArgumentNullException.ThrowIfNull(membership);
        if (pointCount < 0) throw new ArgumentOutOfRangeException(nameof(pointCount));
        if (groupCount < 0) throw new ArgumentOutOfRangeException(nameof(groupCount));
        if (membership.Length != (long)pointCount * groupCount)
            throw new ArgumentException(
                $"membership length ({membership.Length}) must equal pointCount*groupCount " +
                $"({(long)pointCount * groupCount}).", nameof(membership));

        _m = membership;
        PointCount = pointCount;
        GroupCount = groupCount;
    }

    /// <summary>Membership of point <paramref name="i"/> in group <paramref name="k"/>.</summary>
    public double this[int i, int k] => _m[i * GroupCount + k];

    /// <summary>Point <paramref name="i"/>'s membership row (length <see cref="GroupCount"/>).</summary>
    public ReadOnlySpan<double> Row(int i) => _m.AsSpan(i * GroupCount, GroupCount);

    /// <summary>The raw flat row-major membership buffer (length <c>N·K</c>).</summary>
    public ReadOnlySpan<double> Membership => _m;

    /// <summary>
    /// MAP / argmax resolution: assign each point to its highest-membership group — the degenerate baseline
    /// resolution strategy (MATLAB's <c>soft2hard</c>). Never abstains: every point takes its row's argmax
    /// (ties resolve to the first/lowest group id). When <see cref="GroupCount"/> is 0, every point is left
    /// <see cref="Assignment.Unassigned"/>. Richer strategies (threshold-to-abstain, modal ascent, EOM)
    /// live in the resolution layer.
    /// </summary>
    public Assignment Argmax()
    {
        var labels = new int[PointCount];
        if (GroupCount == 0)
        {
            Array.Fill(labels, Assignment.Unassigned);
            return new Assignment { Labels = labels, Count = 0 };
        }

        for (int i = 0; i < PointCount; i++)
        {
            int baseIdx = i * GroupCount;
            int best = 0;
            double bestVal = _m[baseIdx];
            for (int k = 1; k < GroupCount; k++)
            {
                double v = _m[baseIdx + k];
                if (v > bestVal) { bestVal = v; best = k; }
            }
            labels[i] = best;
        }
        return new Assignment { Labels = labels, Count = GroupCount };
    }
}
