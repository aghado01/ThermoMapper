#nullable enable
using System;
using System.Runtime.CompilerServices;

namespace Hashish;

internal static class SeededHash
{
    private const ulong FnvOffsetBasis = 14695981039346656037UL;
    private const ulong FnvPrime = 1099511628211UL;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static ulong Mix64(ulong value)
    {
        value ^= value >> 30;
        value *= 0xbf58476d1ce4e5b9UL;
        value ^= value >> 27;
        value *= 0x94d049bb133111ebUL;
        value ^= value >> 31;
        return value;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static ulong Seed(int seed)
        => Mix64(FnvOffsetBasis ^ (uint)seed);

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static ulong Fnv1a(ReadOnlySpan<char> value, ulong seed)
    {
        ulong hash = seed;
        for (int i = 0; i < value.Length; i++)
        {
            hash ^= value[i];
            hash *= FnvPrime;
        }

        return Mix64(hash);
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static ulong Fnv1a(ReadOnlySpan<byte> value, ulong seed)
    {
        ulong hash = seed;
        for (int i = 0; i < value.Length; i++)
        {
            hash ^= value[i];
            hash *= FnvPrime;
        }

        return Mix64(hash);
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static ulong Fnv1a(ReadOnlySpan<uint> value, ulong seed)
    {
        ulong hash = seed;
        for (int i = 0; i < value.Length; i++)
        {
            uint item = value[i];
            hash ^= item & 0xFFU;
            hash *= FnvPrime;
            hash ^= (item >> 8) & 0xFFU;
            hash *= FnvPrime;
            hash ^= (item >> 16) & 0xFFU;
            hash *= FnvPrime;
            hash ^= (item >> 24) & 0xFFU;
            hash *= FnvPrime;
        }

        return Mix64(hash);
    }
}
