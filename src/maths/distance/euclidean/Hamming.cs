using System.Numerics;

namespace Maths.Distance.Euclidean;

/// <summary>
/// 64-bit packed Hamming distance — XOR + POPCNT on a pair of ulongs.
/// </summary>
public static class Hamming
{
    public static double Distance(ulong a, ulong b)
        => BitOperations.PopCount(a ^ b);
}
