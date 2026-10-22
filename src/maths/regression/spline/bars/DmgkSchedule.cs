using System;
using Maths.Samplers.Rjmcmc;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Green's (1995, §4.1) prior-aware birth/death/relocate schedule, which Denison et al. (1998) and
/// DiMatteo–Genovese–Kass (2001) use with <c>c = 0.4</c>: <c>b_k = c·min{1, p(k+1)/p(k)}</c>,
/// <c>d_k = c·min{1, p(k−1)/p(k)}</c>, relocate = <c>1 − b_k − d_k</c>. It makes <c>b_k·p(k) = d_{k+1}·p(k+1)</c>, so
/// with the engine's selection ratio <c>j_{m'}(x')/j_m(x)</c> the prior cancels from the acceptance.
/// The prior moves from the acceptance into the proposal rate, and because acceptance is capped at 1 this does not
/// speed dimension mixing: the birth/death flow becomes <c>c·min(p_k, p_{k+1})·min(Ã_f, Ã_r)</c> against the
/// constant schedule's <c>c·min(p_k·Ã_f, p_{k+1}·Ã_r)</c>, where Ã is likelihood times proposal density. The two agree
/// under the prior alone, Green's motivating case; with a likelihood, at equal <c>c</c> the flow is never higher, and
/// the unused mass goes to relocation. Birth weight vanishes where the prior has no mass (<c>p(k+1)=0</c>), so a hard
/// cap (e.g. <see cref="UniformPrior"/>) costs no rejected proposals. Wraps existing knot moves via
/// <see cref="ScheduledMove{TState}"/>; the proposal mechanics are unchanged.
/// </summary>
public static class DmgkSchedule
{
    /// <param name="c">Schedule constant in (0, 1], which scales the birth+death proposal mass; DMGK use 0.4. Green
/// instead takes the largest c with b_k + d_k ≤ 0.9 for all k (<see cref="GreenConstant"/>), which can exceed 0.5. The
/// weights sum to one, as the prior's cancellation needs, only while b_k + d_k ≤ 1 at every reachable k.</param>
    public static IRjMove<KnotConfig>[] Wrap(
        IRjMove<KnotConfig> birth, IRjMove<KnotConfig> death, IRjMove<KnotConfig> relocate,
        IComplexityPrior prior, double c = 0.4)
    {
        ArgumentNullException.ThrowIfNull(birth);
        ArgumentNullException.ThrowIfNull(death);
        ArgumentNullException.ThrowIfNull(relocate);
        ArgumentNullException.ThrowIfNull(prior);
        if (!(c > 0.0 && c <= 1.0)) throw new ArgumentOutOfRangeException(nameof(c), "DMGK constant must be in (0, 1].");

        double Birth(KnotConfig s) => c * Ratio(prior, s.Count, s.Count + 1);
        double Death(KnotConfig s) => s.Count == 0 ? 0.0 : c * Ratio(prior, s.Count, s.Count - 1);
        double Relocate(KnotConfig s) => Math.Max(0.0, 1.0 - Birth(s) - Death(s));

        return new IRjMove<KnotConfig>[]
        {
            new ScheduledMove<KnotConfig>(birth, Birth),
            new ScheduledMove<KnotConfig>(death, Death),
            new ScheduledMove<KnotConfig>(relocate, Relocate),
        };
    }

    /// <summary>
    /// Green's (1995, §4.1) choice: the largest c with <c>b_k + d_k ≤ budget</c> for every k in [0, <paramref name="maxKnots"/>],
    /// i.e. <c>budget / max_k(min{1, p(k+1)/p(k)} + min{1, p(k−1)/p(k)})</c>. Green used a budget of 0.9, which leaves
    /// relocation at least 0.1.
    /// </summary>
    public static double GreenConstant(IComplexityPrior prior, int maxKnots, double budget = 0.9)
    {
        ArgumentNullException.ThrowIfNull(prior);
        if (maxKnots < 0) throw new ArgumentOutOfRangeException(nameof(maxKnots));
        if (!(budget > 0.0 && budget <= 1.0)) throw new ArgumentOutOfRangeException(nameof(budget), "The budget must be in (0, 1].");
        double worst = 0.0;
        for (int k = 0; k <= maxKnots; k++)
        {
            if (double.IsNegativeInfinity(prior.LogPrior(k))) continue;   // outside the prior's support
            double sum = Ratio(prior, k, k + 1) + (k == 0 ? 0.0 : Ratio(prior, k, k - 1));
            worst = Math.Max(worst, sum);
        }
        return Math.Min(1.0, budget / worst);
    }

    // min{1, p(to)/p(from)} from the log-prior difference; 0 when the target carries no prior mass.
    private static double Ratio(IComplexityPrior prior, int from, int to)
    {
        double delta = prior.LogPrior(to) - prior.LogPrior(from);
        if (double.IsNegativeInfinity(delta)) return 0.0;
        return delta >= 0.0 ? 1.0 : Math.Exp(delta);
    }
}
