#nullable enable
using System;
using System.Runtime.CompilerServices;

namespace Hashish;

/// <summary>Count-Min Sketch for streaming approximate frequency estimation.</summary>
public sealed class CountMin
{
    private readonly long[,] _table;
    private readonly int _width;
    private readonly int _depth;

    public CountMin(int width, int depth)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(width, 0, nameof(width));
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(depth, 0, nameof(depth));

        _width = width;
        _depth = depth;
        _table = new long[depth, width];
    }

    public int Width => _width;
    public int Depth => _depth;
    public long TotalCount { get; private set; }

    public static CountMin Create(double epsilon, double delta)
    {
        if (epsilon <= 0.0 || epsilon >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(epsilon), "Must be in (0, 1).");
        if (delta <= 0.0 || delta >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(delta), "Must be in (0, 1).");

        int width = (int)Math.Ceiling(Math.E / epsilon);
        int depth = (int)Math.Ceiling(Math.Log(1.0 / delta));
        return new CountMin(width, depth);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Add(string value, long count = 1)
    {
        ArgumentNullException.ThrowIfNull(value);
        AddHash(SeededHash.Fnv1a(value.AsSpan(), SeededHash.Seed(0)), count);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Add(ReadOnlySpan<byte> value, long count = 1)
        => AddHash(SeededHash.Fnv1a(value, SeededHash.Seed(0)), count);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long Estimate(string value)
    {
        ArgumentNullException.ThrowIfNull(value);
        return EstimateHash(SeededHash.Fnv1a(value.AsSpan(), SeededHash.Seed(0)));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long Estimate(ReadOnlySpan<byte> value)
        => EstimateHash(SeededHash.Fnv1a(value, SeededHash.Seed(0)));

    public void Clear()
    {
        Array.Clear(_table);
        TotalCount = 0;
    }

    private void AddHash(ulong hash, long count)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(count, 0, nameof(count));

        for (int row = 0; row < _depth; row++)
        {
            int column = Column(hash, row);
            _table[row, column] += count;
        }

        TotalCount += count;
    }

    private long EstimateHash(ulong hash)
    {
        long estimate = long.MaxValue;
        for (int row = 0; row < _depth; row++)
        {
            int column = Column(hash, row);
            estimate = Math.Min(estimate, _table[row, column]);
        }

        return estimate == long.MaxValue ? 0 : estimate;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private int Column(ulong hash, int row)
        => (int)(SeededHash.Mix64(hash ^ SeededHash.Seed(row + 1)) % (uint)_width);
}
