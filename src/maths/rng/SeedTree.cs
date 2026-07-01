using System;

namespace Maths.Rng;

/// <summary>
/// Deterministic seed expansion: grows one master seed into N decorrelated child seeds via SplitMix64 (the xoshiro
/// initializer of Blackman &amp; Vigna). Any fan-out — an ensemble of chains, a tempering ladder, a pool of parallel
/// workers — becomes reproducible from a single integer with every branch on an independent stream, the provenance
/// a run records as <c>{requested master, resolved children}</c>. Lives beside <see cref="Xoshiro256PlusPlus"/> in
/// <see cref="Maths.Rng"/> so every layer (samplers, regression) shares one mixer rather than re-inlining it.
/// </summary>
public static class SeedTree
{
    /// <summary>Derive <paramref name="count"/> decorrelated seeds from <paramref name="master"/>.</summary>
    public static int[] Derive(int master, int count)
    {
        if (count < 0) throw new ArgumentOutOfRangeException(nameof(count));
        var seeds = new int[count];
        ulong state = unchecked((ulong)master);
        for (int i = 0; i < count; i++)
            seeds[i] = unchecked((int)SplitMix64(ref state));
        return seeds;
    }

    private static ulong SplitMix64(ref ulong state)
    {
        unchecked
        {
            state += 0x9e3779b97f4a7c15UL;
            ulong z = state;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9UL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebUL;
            return z ^ (z >> 31);
        }
    }
}
