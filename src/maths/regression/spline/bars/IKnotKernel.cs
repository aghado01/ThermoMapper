using Maths.Rng;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// An elementary knot proposal kernel <c>h(x | center)</c> — the swappable locality policy of the free-knot
/// moves (DMGK's declared proposal component). <see cref="UniformKernel"/> ignores the center (the Denison
/// baseline); <see cref="LocalBetaKernel"/> concentrates proposals near the center. Birth/death consume it as
/// a mixture over existing knots (<see cref="ProposalMath.LogMixtureDensity"/>); relocate uses it directly and
/// asymmetrically.
/// </summary>
public interface IKnotKernel
{
    /// <summary>Draw a candidate knot near <paramref name="center"/>.</summary>
    double Sample(double center, Xoshiro256PlusPlus rng);

    /// <summary>Log density <c>log h(x | center)</c> of proposing <paramref name="x"/> from that center.</summary>
    double LogDensity(double x, double center);
}
