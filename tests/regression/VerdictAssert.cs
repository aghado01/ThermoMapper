using System.Collections.Generic;
using System.Linq;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>Checks on a <see cref="ConvergenceVerdict"/> shared by the ensemble clients' tests.</summary>
internal static class VerdictAssert
{
    /// <summary>
    /// The verdict is the shared rule read on <paramref name="readings"/>, each functional's R̂ and ESS by name (a null ESS
    /// for an unmonitored functional): every functional named once across the groups, each failing exactly when its R̂
    /// misses the target or its ESS the gate over <paramref name="chains"/> chains, with the reading the result reports.
    /// </summary>
    public static void Consistent(ConvergenceVerdict v, Dictionary<string, (RHatDiagnostic RHat, EssDiagnostic? Ess)> readings, int chains)
    {
        var failing = v.Groups.SelectMany(g => g.Failing).ToDictionary(f => f.Name);
        foreach ((string name, (RHatDiagnostic rhat, EssDiagnostic? ess)) in readings)
        {
            bool fails = !(rhat.MeetsTarget(v.Target) && (ess is not { } e || e.MeetsMinimum(v.MinEssPerChain * chains)));
            Assert.True(fails == failing.ContainsKey(name), $"{name}: R̂ {rhat}, ESS {ess}; the verdict {(fails ? "passes" : "fails")} it");
            if (!fails) continue;
            Assert.Equal(rhat, failing[name].RHat);
            Assert.Equal(ess, failing[name].Ess);
        }
        Assert.All(failing.Keys, name => Assert.True(readings.ContainsKey(name), $"the verdict fails {name}, which the result does not report"));
        foreach (GroupVerdict g in v.Groups)
            if (g.Failing.Count > 0) Assert.Equal(ConvergenceOutcome.NotConverged, g.Outcome);
        Assert.Equal(v.Groups.Any(g => g.GatesStop && g.Outcome == ConvergenceOutcome.NotConverged), v.Outcome == ConvergenceOutcome.NotConverged);
    }

    /// <summary>Two verdicts alike: outcome, target, gate, and each group's outcome, gating and failing readings.</summary>
    public static void Same(ConvergenceVerdict a, ConvergenceVerdict b)
    {
        Assert.Equal((a.Outcome, a.Target, a.MinEssPerChain), (b.Outcome, b.Target, b.MinEssPerChain));
        Assert.Equal(a.Groups.Count, b.Groups.Count);
        for (int g = 0; g < a.Groups.Count; g++)
        {
            Assert.Equal((a.Groups[g].Name, a.Groups[g].Outcome, a.Groups[g].GatesStop), (b.Groups[g].Name, b.Groups[g].Outcome, b.Groups[g].GatesStop));
            Assert.Equal(a.Groups[g].Failing, b.Groups[g].Failing);
        }
    }
}
