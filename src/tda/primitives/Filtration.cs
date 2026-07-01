#nullable enable
using System;
using System.Collections.Generic;
using Graphs.Primitives;

namespace TDA.Primitives;

// ── NerveFiltrationFrame ──────────────────────────────────────────────────────

/// <summary>
/// One frame in an ordered nerve filtration.
/// <para>
/// <see cref="Nerve"/> is the nerve (or proximity-graph skeleton) at this parameter value.
/// <see cref="NodeMemberIndices"/>: <c>NodeMemberIndices[i]</c> is the array of original
/// data-point indices belonging to nerve node <c>i</c>. Stored alongside the nerve
/// so cross-frame component matching can use point identity rather than just graph topology.
/// </para>
/// <para>
/// For Mapper filtrations: nerve = Mapper nerve graph; node members = Mapper cluster members.
/// For future PH filtrations: nerve = simplicial complex skeleton; node members = simplex vertices.
/// </para>
/// </summary>
public sealed record NerveFiltrationFrame(
    double ParameterValue,
    CsrGraph Nerve,
    int[][] NodeMemberIndices,
    int FrameIndex);

// ── NerveFiltration ───────────────────────────────────────────────────────────

/// <summary>
/// An ordered sequence of nerve frames over a scalar parameter axis (e.g., SPC temperature T,
/// epsilon distance threshold, or Mapper cover scale).
/// <para>
/// Frames must be in non-decreasing <see cref="NerveFiltrationFrame.ParameterValue"/> order.
/// </para>
/// <para>
/// First concrete consumer: Persistent Mapper over the SPC T-sweep. The second consumer
/// (classical PH on distance-matrix filtrations) is out of current scope but the type
/// is designed to accommodate it without modification.
/// </para>
/// </summary>
public sealed class NerveFiltration
{
    public IReadOnlyList<NerveFiltrationFrame> Frames { get; }
    public string ParameterLabel { get; }

    public NerveFiltration(IReadOnlyList<NerveFiltrationFrame> frames, string parameterLabel = "T")
    {
        ArgumentNullException.ThrowIfNull(frames);
        for (int i = 1; i < frames.Count; i++)
            if (frames[i].ParameterValue < frames[i - 1].ParameterValue)
                throw new ArgumentException(
                    $"NerveFiltration frames must be in non-decreasing parameter order: " +
                    $"frame {i - 1} has parameter {frames[i - 1].ParameterValue} " +
                    $"but frame {i} has {frames[i].ParameterValue}.");
        Frames = frames;
        ParameterLabel = parameterLabel;
    }

    /// <summary>
    /// Compute the <see cref="NerveDiff"/> for each consecutive frame pair.
    /// Returns an empty list when the filtration has fewer than two frames.
    /// </summary>
    public IReadOnlyList<NerveDiff> ComputeDiffs()
    {
        if (Frames.Count < 2) return Array.Empty<NerveDiff>();
        var diffs = new NerveDiff[Frames.Count - 1];
        for (int i = 0; i < diffs.Length; i++)
            diffs[i] = NerveDiff.Compute(Frames[i], Frames[i + 1]);
        return diffs;
    }
}
