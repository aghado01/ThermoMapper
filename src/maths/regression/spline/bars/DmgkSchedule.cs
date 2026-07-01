using System;
using Maths.Samplers.Rjmcmc;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// The DiMatteo–Genovese–Kass (2001) birth/death/relocate schedule: state-dependent move-selection weights that
/// fold the knot-count prior into the proposal — <c>b_k = c·min{1, p(k+1)/p(k)}</c>, <c>d_k = c·min{1,
/// p(k−1)/p(k)}</c>, relocate = <c>1 − b_k − d_k</c>. Because the Green-general engine applies the selection
/// ratio <c>j_{m'}(x')/j_m(x)</c>, the prior then cancels between the posterior and selection terms, so a
/// dimension change is accepted on the <i>likelihood</i> evidence alone — the DMGK mixing improvement. Birth
/// weight vanishes where the prior has no mass (<c>p(k+1)=0</c>), so a hard cap (e.g. <see cref="UniformPrior"/>)
/// is honored automatically. Wraps existing knot moves via <see cref="ScheduledMove{TState}"/>; the proposal
/// mechanics are unchanged.
/// </summary>
public static class DmgkSchedule
{
    /// <param name="c">DMGK constant in (0, 0.5] (paper default 0.4) — caps the birth+death proposal mass.</param>
    public static IRjMove<KnotConfig>[] Wrap(
        IRjMove<KnotConfig> birth, IRjMove<KnotConfig> death, IRjMove<KnotConfig> relocate,
        IComplexityPrior prior, double c = 0.4)
    {
        ArgumentNullException.ThrowIfNull(birth);
        ArgumentNullException.ThrowIfNull(death);
        ArgumentNullException.ThrowIfNull(relocate);
        ArgumentNullException.ThrowIfNull(prior);
        if (!(c > 0.0 && c <= 0.5)) throw new ArgumentOutOfRangeException(nameof(c), "DMGK constant must be in (0, 0.5].");

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

    // min{1, p(to)/p(from)} from the log-prior difference; 0 when the target carries no prior mass.
    private static double Ratio(IComplexityPrior prior, int from, int to)
    {
        double delta = prior.LogPrior(to) - prior.LogPrior(from);
        if (double.IsNegativeInfinity(delta)) return 0.0;
        return delta >= 0.0 ? 1.0 : Math.Exp(delta);
    }
}
