namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// One Swendsen–Wang equilibration's MC budget — burn-in cycles discarded, then measurement
/// cycles accumulated. A single bundled primitive rather than two loose ints: the sweep and the
/// chosen-T equilibrium each carry an instance (<c>SweepBudget</c> / <c>EquilibriumBudget</c>), and
/// the CLI binds one <c>burnin=…,cycles=…</c> value into it at the shell boundary.
/// </summary>
public readonly record struct RunBudget(int BurnIn, int Cycles);
