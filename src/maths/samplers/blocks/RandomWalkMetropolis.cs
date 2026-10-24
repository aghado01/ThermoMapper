using System;
using Maths.Distributions;
using Maths.Rng;

namespace Maths.Samplers.Blocks;

/// <summary>
/// One random-walk Metropolis step on a scalar: propose <c>x′ = x + s·z</c>, z ~ N(0, 1), and accept with probability
/// min(1, exp(f(x′) − f(x))) against the cached log target f(x). The uniform is drawn on every step, as
/// <c>log(1 − U) &lt; f(x′) − f(x)</c>, the form BAPS's λ walk uses. BGP's bandwidth walk skips the uniform when the ratio is
/// at least 1, which consumes the stream differently, so it keeps its own step.
/// </summary>
public static class RandomWalkMetropolis
{
    /// <summary>Take one step from <paramref name="x"/>, whose log target is <paramref name="logTarget"/>; on acceptance
    /// both move to the proposal's. Returns whether the proposal was accepted.</summary>
    public static bool Step(ref double x, ref double logTarget, double step, Func<double, double> logDensity, Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(logDensity);
        ArgumentNullException.ThrowIfNull(rng);
        double proposal = x + step * Normal.SampleStandard(rng);
        double proposalTarget = logDensity(proposal);
        if (!(Math.Log(1.0 - rng.NextDouble()) < proposalTarget - logTarget)) return false;
        x = proposal;
        logTarget = proposalTarget;
        return true;
    }
}
