#nullable enable
using System;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace Hashish;

/// <summary>HyperLogLog cardinality estimator for approximate distinct counts.</summary>
public sealed class HyperLogLog
{
    private readonly byte[] _registers;
    private readonly int _precision;
    private readonly int _registerCount;

    public HyperLogLog(int precision = 14)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(precision, 4, nameof(precision));
        ArgumentOutOfRangeException.ThrowIfGreaterThan(precision, 18, nameof(precision));

        _precision = precision;
        _registerCount = 1 << precision;
        _registers = new byte[_registerCount];
    }

    public int Precision => _precision;
    public int RegisterCount => _registerCount;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Add(string value)
    {
        ArgumentNullException.ThrowIfNull(value);
        AddHash(SeededHash.Fnv1a(value.AsSpan(), SeededHash.Seed(0)));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Add(ReadOnlySpan<byte> value)
        => AddHash(SeededHash.Fnv1a(value, SeededHash.Seed(0)));

    public double Estimate()
    {
        double sum = 0.0;
        int zeroRegisters = 0;

        for (int i = 0; i < _registers.Length; i++)
        {
            int rank = _registers[i];
            sum += Math.Pow(2.0, -rank);
            if (rank == 0)
                zeroRegisters++;
        }

        double alpha = Alpha(_registerCount);
        double estimate = alpha * _registerCount * _registerCount / sum;

        if (estimate <= 2.5 * _registerCount && zeroRegisters > 0)
            return _registerCount * Math.Log(_registerCount / (double)zeroRegisters);

        return estimate;
    }

    public void Merge(HyperLogLog other)
    {
        ArgumentNullException.ThrowIfNull(other);
        if (other._precision != _precision)
            throw new ArgumentException("HyperLogLog precision must match.", nameof(other));

        for (int i = 0; i < _registers.Length; i++)
            _registers[i] = Math.Max(_registers[i], other._registers[i]);
    }

    public void Clear() => Array.Clear(_registers);

    private void AddHash(ulong hash)
    {
        int index = (int)(hash & (ulong)(_registerCount - 1));
        ulong remainder = hash >> _precision;
        int rank = Rank(remainder);
        if (rank > _registers[index])
            _registers[index] = (byte)rank;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private int Rank(ulong remainder)
    {
        int width = 64 - _precision;
        if (remainder == 0UL)
            return width + 1;

        return BitOperations.LeadingZeroCount(remainder) - _precision + 1;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Alpha(int m) => m switch
    {
        16 => 0.673,
        32 => 0.697,
        64 => 0.709,
        _ => 0.7213 / (1.0 + 1.079 / m)
    };
}
