using System;
using System.Collections;
using System.Collections.Generic;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// One chain's retained draws, in draw order: the record "record, then reduce" keeps beside the shell's traces, so an
/// engine's summaries are views over it and a new readout needs no rerun. Every draw is kept by default; an
/// <see cref="Every"/> above 1 keeps one in k, a memory option that never stands in for the record. <see cref="Copy"/>
/// gives an independent store for a checkpoint, with the thinning's phase: the list is copied and the draws are shared,
/// since a retained draw is never written again.
/// </summary>
public sealed class DrawStore<TDraw> : IReadOnlyList<TDraw>
{
    private readonly List<TDraw> _draws = new();

    /// <param name="every">Keep one draw in this many, the first included; 1, the default, keeps all.</param>
    public DrawStore(int every = 1)
    {
        if (every < 1) throw new ArgumentOutOfRangeException(nameof(every), "A store keeps at least one draw in every.");
        Every = every;
    }

    /// <summary>The thinning: one draw kept in this many.</summary>
    public int Every { get; }

    /// <summary>Draws offered so far, kept or not: the thinning's phase, which a copy carries.</summary>
    public long Offered { get; private set; }

    /// <inheritdoc/>
    public int Count => _draws.Count;

    /// <inheritdoc/>
    public TDraw this[int index] => _draws[index];

    /// <summary>Offer the next draw; it is kept when the thinning's phase allows.</summary>
    public void Add(TDraw draw)
    {
        if (Offered % Every == 0) _draws.Add(draw);
        Offered++;
    }

    /// <summary>Append another store's kept draws after this one's, as when chains' records are pooled.</summary>
    public void AppendFrom(DrawStore<TDraw> other)
    {
        ArgumentNullException.ThrowIfNull(other);
        _draws.AddRange(other._draws);
        Offered += other.Offered;
    }

    /// <summary>An independent store with the same draws, thinning and phase.</summary>
    public DrawStore<TDraw> Copy()
    {
        var copy = new DrawStore<TDraw>(Every) { Offered = Offered };
        copy._draws.AddRange(_draws);
        return copy;
    }

    /// <summary>The kept draws, in order.</summary>
    public TDraw[] ToArray() => _draws.ToArray();

    /// <inheritdoc/>
    public IEnumerator<TDraw> GetEnumerator() => _draws.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}
