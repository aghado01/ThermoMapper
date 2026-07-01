using System;
using System.Numerics;
using System.Security.Cryptography;

namespace Maths.Rng;

/// <summary>
/// xoshiro256++ pseudo-random generator (Blackman &amp; Vigna 2018) with
/// raw 4-word state save/restore. .NET 6+'s unseeded <see cref="Random"/>
/// and <see cref="Random.Shared"/> use a xoshiro256 family generator
/// internally, but do not expose its state — this class exists so SPC
/// checkpoint artifacts can roundtrip RNG state bit-identically across
/// save/restore boundaries.
/// </summary>
/// <remarks>
/// <para><b>Variant.</b> This is xoshiro256++ — the scrambler is
/// <c>rotl(s0 + s3, 23) + s0</c>. Equivalent statistical quality to
/// xoshiro256** but a different output function; the two are not
/// interchangeable. The state-advance is identical across both variants.</para>
///
/// <para><b>Seeding.</b> An <see langword="int"/> seed is expanded to four
/// state words via SplitMix64 (Vigna's recommended xoshiro initializer).
/// A <see langword="null"/> seed draws 32 bytes from
/// <see cref="RandomNumberGenerator"/> for OS-entropy bootstrapping.
/// The all-zero state is a fixed point of the recurrence and is guarded
/// against in both seeding paths.</para>
///
/// <para><b>Provisional placement.</b> Lives under
/// <c>Clustering.Graphical.SPC.Potts</c> for now to keep the SPC spine compiling,
/// but xoshiro is a generic numerics primitive — should migrate to
/// <c>src/numerics</c> (or wherever shared numerical primitives land)
/// when that namespace is established.</para>
///
/// <para><b>Reference.</b>
/// https://prng.di.unimi.it/ — public-domain reference implementation.</para>
/// </remarks>
public sealed class Xoshiro256PlusPlus
{
    private ulong _s0, _s1, _s2, _s3;

    /// <summary>
    /// Construct from an optional 32-bit seed. <see langword="null"/>
    /// draws from <see cref="RandomNumberGenerator"/> for non-reproducible
    /// runs; pass a concrete seed for reproducibility across runs (and
    /// for replica generation in the SPC scheduler).
    /// </summary>
    public Xoshiro256PlusPlus(int? seed = null)
    {
        if (seed is int s)
        {
            ulong x = unchecked((ulong)s);
            _s0 = SplitMix64(ref x);
            _s1 = SplitMix64(ref x);
            _s2 = SplitMix64(ref x);
            _s3 = SplitMix64(ref x);
        }
        else
        {
            Span<byte> bytes = stackalloc byte[32];
            RandomNumberGenerator.Fill(bytes);
            _s0 = BitConverter.ToUInt64(bytes.Slice(0, 8));
            _s1 = BitConverter.ToUInt64(bytes.Slice(8, 8));
            _s2 = BitConverter.ToUInt64(bytes.Slice(16, 8));
            _s3 = BitConverter.ToUInt64(bytes.Slice(24, 8));
        }

        // All-zero state is a fixed point of the recurrence; vanishingly
        // unlikely from RandomNumberGenerator.Fill, theoretically reachable
        // from a pathological int seed. Reseed s0 with the golden-ratio
        // constant if it happens.
        if ((_s0 | _s1 | _s2 | _s3) == 0)
            _s0 = 0x9e3779b97f4a7c15UL;
    }

    private Xoshiro256PlusPlus(ulong s0, ulong s1, ulong s2, ulong s3)
    {
        _s0 = s0; _s1 = s1; _s2 = s2; _s3 = s3;
    }

    /// <summary>Snapshot the current 4-word state for checkpointing.</summary>
    public (ulong s0, ulong s1, ulong s2, ulong s3) SaveState()
        => (_s0, _s1, _s2, _s3);

    /// <summary>
    /// Overwrite this generator's state in place. The next call to
    /// <see cref="NextUInt64"/> produces output bit-identical to the
    /// generator at the moment the state was saved. Used by checkpoint
    /// restore where the consumer wants to keep the existing instance
    /// (and its readonly field membership) rather than allocate a fresh
    /// one via <see cref="Restore"/>.
    /// </summary>
    public void LoadState(ulong s0, ulong s1, ulong s2, ulong s3)
    {
        if ((s0 | s1 | s2 | s3) == 0)
            throw new ArgumentException(
                "All-zero state is a fixed point of xoshiro256++ and would lock the generator.");
        _s0 = s0; _s1 = s1; _s2 = s2; _s3 = s3;
    }

    /// <summary>
    /// Reconstruct a generator from a saved state. The restored instance
    /// produces output bit-identical to the generator at the moment
    /// <see cref="SaveState"/> was called.
    /// </summary>
    public static Xoshiro256PlusPlus Restore(ulong s0, ulong s1, ulong s2, ulong s3)
    {
        if ((s0 | s1 | s2 | s3) == 0)
            throw new ArgumentException(
                "All-zero state is a fixed point of xoshiro256++ and would lock the generator.");
        return new Xoshiro256PlusPlus(s0, s1, s2, s3);
    }

    /// <summary>Advance the state and return the next 64-bit output.</summary>
    public ulong NextUInt64()
    {
        ulong result = BitOperations.RotateLeft(_s0 + _s3, 23) + _s0;

        ulong t = _s1 << 17;
        _s2 ^= _s0;
        _s3 ^= _s1;
        _s1 ^= _s2;
        _s0 ^= _s3;
        _s2 ^= t;
        _s3 = BitOperations.RotateLeft(_s3, 45);

        return result;
    }

    /// <summary>
    /// Uniform <see langword="double"/> in <c>[0, 1)</c>. Uses the top 53
    /// bits of one generator output — the canonical conversion for IEEE-754
    /// double precision.
    /// </summary>
    public double NextDouble()
        => (NextUInt64() >> 11) * (1.0 / (1UL << 53));

    /// <summary>
    /// Uniform <see langword="int"/> in <c>[0, q)</c>. Multiplication-shift
    /// on the upper 32 bits of a state output — bias is <c>q / 2^32</c>,
    /// negligible for Potts-scale <c>q</c> (≤ ~30).
    /// </summary>
    public int NextInt(int q)
    {
        if (q <= 0) throw new ArgumentOutOfRangeException(nameof(q));
        return (int)(((NextUInt64() >> 32) * (ulong)(uint)q) >> 32);
    }

    /// <summary>
    /// SplitMix64 step — Vigna's recommended xoshiro seed initializer.
    /// Each call advances <paramref name="state"/> by the golden-ratio
    /// constant and returns a well-distributed 64-bit word.
    /// </summary>
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
