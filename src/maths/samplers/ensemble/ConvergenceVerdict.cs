using System;
using System.Collections.Generic;

namespace Maths.Samplers.Ensemble;

/// <summary>A run's convergence outcome (D10), read by <see cref="ConvergenceVerdict"/>.</summary>
public enum ConvergenceOutcome
{
    /// <summary>Every gating group passes, and at least one of them holds an estimate.</summary>
    Converged,

    /// <summary>A gating group fails: a functional's R̂ misses the target or its ESS the gate, or its reading is frozen
    /// or unavailable.</summary>
    NotConverged,

    /// <summary>Nothing could be estimated: one chain, or every gating functional constant.</summary>
    Undetermined,
}

/// <summary>
/// A named group of a model's functionals, by index, whose agreement the verdict reads together (D10's agreement groups:
/// BASS's predictive, level and structural). A model declares its groups (<see cref="IEnsembleModel{TDraw}.Groups"/>),
/// and every functional belongs to exactly one.
/// </summary>
public sealed record FunctionalGroup(string Name, IReadOnlyList<int> Functionals)
{
    /// <summary>The name of the one group a model that declares none gets.</summary>
    public const string AllName = "all";

    /// <summary>One group, <see cref="AllName"/>, holding every one of <paramref name="functionals"/> functionals.</summary>
    public static FunctionalGroup All(int functionals)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(functionals);
        var indices = new int[functionals];
        for (int f = 0; f < functionals; f++) indices[f] = f;
        return new FunctionalGroup(AllName, indices);
    }
}

/// <summary>
/// One functional's reading: its index and name, its R̂ (rank-normalized when it is monitored, classic otherwise) and its
/// ESS, the less favourable of bulk and tail; null for an unmonitored functional, which keeps no trace.
/// </summary>
public sealed record FunctionalReading(int Index, string Name, RHatDiagnostic RHat, EssDiagnostic? Ess);

/// <summary>
/// One group's verdict: its outcome under the stop's rules, whether it gates an R̂ stop (it is in the run's stopping
/// scope), and the functionals that failed, in index order.
/// </summary>
public sealed record GroupVerdict(string Name, ConvergenceOutcome Outcome, bool GatesStop, IReadOnlyList<FunctionalReading> Failing);

/// <summary>
/// A run's convergence verdict (D10), read at the end of every <see cref="ChainEnsemble"/> run, and between rounds when an
/// R̂ target is set, by B01's rules applied to each of the model's groups:
/// <list type="bullet">
/// <item>a monitored functional passes when its rank-normalized R̂ meets <see cref="Target"/> and its bulk and tail ESS
/// reach <see cref="MinEssPerChain"/> per chain; an unmonitored one passes on its classic R̂;</item>
/// <item>a frozen or unavailable reading fails; a constant functional passes but cannot certify a group on its own, so a
/// group with no estimate is <see cref="ConvergenceOutcome.Undetermined"/>, and so is every group of one chain;</item>
/// <item>any failing gating group makes the run <see cref="ConvergenceOutcome.NotConverged"/>; otherwise every gating group
/// undetermined makes it <see cref="ConvergenceOutcome.Undetermined"/>; otherwise it is
/// <see cref="ConvergenceOutcome.Converged"/>.</item>
/// </list>
/// <see cref="Groups"/> reports every group, gating or not, in the model's order. <see cref="Target"/> is the run's R̂
/// target when it set one, and <see cref="ChainEnsemble.ReferenceRHatTarget"/> otherwise.
/// </summary>
public sealed record ConvergenceVerdict(
    ConvergenceOutcome Outcome, double Target, double MinEssPerChain, IReadOnlyList<GroupVerdict> Groups)
{
    /// <summary>The verdict of the group named <paramref name="name"/>.</summary>
    public GroupVerdict Group(string name)
    {
        ArgumentNullException.ThrowIfNull(name);
        foreach (GroupVerdict g in Groups)
            if (g.Name == name) return g;
        throw new KeyNotFoundException($"The verdict has no group named \"{name}\".");
    }

    // The rule, which the stop between rounds and the verdict at the end both read: each group judged on the functionals'
    // classic R̂ (every functional) and rank-normalized diagnostics (the leading monitored ones), over `chains` chains.
    internal static ConvergenceVerdict Of(
        ReadOnlySpan<RHatDiagnostic> classic, ReadOnlySpan<ConvergenceDiagnostic> monitored, int chains, IReadOnlyList<string> names,
        IReadOnlyList<FunctionalGroup> groups, bool[] gates, double target, double minEssPerChain)
    {
        var verdicts = new GroupVerdict[groups.Count];
        bool anyFails = false, anyCertifies = false;
        for (int g = 0; g < verdicts.Length; g++)
        {
            var failing = new List<FunctionalReading>();
            bool anyEstimate = false;
            if (chains >= 2)
                foreach (int f in groups[g].Functionals)
                {
                    bool passes;
                    RHatDiagnostic rhat;
                    EssDiagnostic? ess = null;
                    if (f < monitored.Length)
                    {
                        ConvergenceDiagnostic d = monitored[f];
                        rhat = d.RHat;
                        ess = d.Ess;
                        passes = rhat.MeetsTarget(target) && d.Ess.MeetsMinimum(minEssPerChain * chains);
                    }
                    else
                    {
                        rhat = classic[f];
                        passes = rhat.MeetsTarget(target);
                    }
                    if (!passes) failing.Add(new FunctionalReading(f, names[f], rhat, ess));
                    anyEstimate |= rhat.State == DiagnosticState.Estimated;
                }
            failing.Sort((a, b) => a.Index.CompareTo(b.Index));
            ConvergenceOutcome outcome = failing.Count > 0 ? ConvergenceOutcome.NotConverged
                : anyEstimate ? ConvergenceOutcome.Converged
                : ConvergenceOutcome.Undetermined;
            verdicts[g] = new GroupVerdict(groups[g].Name, outcome, gates[g], failing);
            if (!gates[g]) continue;
            anyFails |= outcome == ConvergenceOutcome.NotConverged;
            anyCertifies |= outcome == ConvergenceOutcome.Converged;
        }
        ConvergenceOutcome run = anyFails ? ConvergenceOutcome.NotConverged
            : anyCertifies ? ConvergenceOutcome.Converged
            : ConvergenceOutcome.Undetermined;
        return new ConvergenceVerdict(run, target, minEssPerChain, verdicts);
    }
}
