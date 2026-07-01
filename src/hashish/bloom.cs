#nullable enable
using System;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace Hashish;

/// <summary>Classic Bloom filter for approximate membership queries.</summary>
public sealed class BloomFilter
{
    private readonly ulong[] _bits;
    private readonly int _bitCount;
    private readonly int _hashCount;

    public BloomFilter(int bitCount, int hashCount)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(bitCount, 0, nameof(bitCount));
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(hashCount, 0, nameof(hashCount));

        _bitCount = bitCount;
        _hashCount = hashCount;
        _bits = new ulong[(bitCount + 63) / 64];
    }

    public int BitCount => _bitCount;
    public int HashCount => _hashCount;
    public long Insertions { get; private set; }

    public static BloomFilter Create(int expectedItems, double falsePositiveProbability)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(expectedItems, 0, nameof(expectedItems));
        if (falsePositiveProbability <= 0.0 || falsePositiveProbability >= 1.0)
            throw new ArgumentOutOfRangeException(nameof(falsePositiveProbability), "Must be in (0, 1).");

        double ln2 = Math.Log(2.0);
        int bits = (int)Math.Ceiling(-(expectedItems * Math.Log(falsePositiveProbability)) / (ln2 * ln2));
        int hashes = Math.Max(1, (int)Math.Round((bits / (double)expectedItems) * ln2));
        return new BloomFilter(bits, hashes);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Add(string value)
    {
        ArgumentNullException.ThrowIfNull(value);
        AddHash(SeededHash.Fnv1a(value.AsSpan(), SeededHash.Seed(0)));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Add(ReadOnlySpan<byte> value)
        => AddHash(SeededHash.Fnv1a(value, SeededHash.Seed(0)));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Contains(string value)
    {
        ArgumentNullException.ThrowIfNull(value);
        return ContainsHash(SeededHash.Fnv1a(value.AsSpan(), SeededHash.Seed(0)));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public bool Contains(ReadOnlySpan<byte> value)
        => ContainsHash(SeededHash.Fnv1a(value, SeededHash.Seed(0)));

    public double FillRatio()
    {
        long setBits = 0;
        for (int i = 0; i < _bits.Length; i++)
            setBits += BitOperations.PopCount(_bits[i]);

        return (double)setBits / _bitCount;
    }

    public void Clear()
    {
        Array.Clear(_bits);
        Insertions = 0;
    }

    private void AddHash(ulong hash)
    {
        ulong h1 = hash;
        ulong h2 = SeededHash.Mix64(hash ^ 0x9e3779b97f4a7c15UL) | 1UL;

        for (int i = 0; i < _hashCount; i++)
        {
            int bit = (int)((h1 + ((ulong)i * h2)) % (uint)_bitCount);
            int word = bit >> 6;
            _bits[word] |= 1UL << (bit & 63);
        }

        Insertions++;
    }

    private bool ContainsHash(ulong hash)
    {
        ulong h1 = hash;
        ulong h2 = SeededHash.Mix64(hash ^ 0x9e3779b97f4a7c15UL) | 1UL;

        for (int i = 0; i < _hashCount; i++)
        {
            int bit = (int)((h1 + ((ulong)i * h2)) % (uint)_bitCount);
            int word = bit >> 6;
            if ((_bits[word] & (1UL << (bit & 63))) == 0UL)
                return false;
        }

        return true;
    }
}
