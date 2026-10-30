using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// D10's verdict on BASS (N14): each retained draw's log posterior among the monitored quantities, the three groups
/// (predictive, level, structural) and their names, the stopping scope, and the verdict and mode-separation readings on
/// <see cref="BassDiagnostics"/>; on S10's fits, the verdict each group reads where F-M2's pass read its R̂.
/// </summary>
public sealed class BassVerdictTests
{
    private readonly ITestOutputHelper _out;
    public BassVerdictTests(ITestOutputHelper output) => _out = output;

    /// <summary>
    /// The log posterior the run monitors is the chain's own log density at its state: each chain started and burned by the
    /// run's own client, and beside it the same chain replayed as the run takes it, untempered or its ladder burned and
    /// scanned; at three draws per chain the states match and the monitored value equals the replayed chain's b = 1 density.
    /// </summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void LogPosterior_IsTheChainsOwnDensity(bool tempered)
    {
        (BassConfig config, BassData data) = BassBitIdentityTests.Fixture(BassMixing.Friedman, "equal", tempered);
        BassSampling s = config.Sampling;
        BassModel model = BassModel.Create(config, data);
        var ensemble = new BassEnsemble(model, Bass.Fingerprint(config, data), new BassRunOptions(), 0, s.Samples);
        int lp = ensemble.FunctionalNames.ToList().IndexOf("log_posterior");
        int[] seeds = SeedTree.Derive(s.MasterSeed, 2 * s.Chains);
        for (int c = 0; c < s.Chains; c++)
        {
            IChainRun<BassState> run = ensemble.StartChain(seeds[c], seeds[s.Chains + c]);
            run.Burn(s.BurnIn);
            Func<BassState> step;
            Func<double> own;
            if (tempered)
            {
                BassLadder ladder = BassMixing.Ladder(model, c, s.Tempering!);
                BassMixing.Burn(ladder, s.BurnIn);
                step = ladder.Step;
                own = () => ladder.Moves.CurrentLogDensity;
            }
            else
            {
                BassChain chain = BassMixing.Chain(model, c);
                for (int k = 0; k < s.BurnIn; k++) chain.Step();
                step = chain.Step;
                own = () => chain.Moves.CurrentLogDensity;
            }
            var functionals = new double[ensemble.FunctionalDim];
            for (int k = 0; k < s.Samples; k++)
            {
                BassState state = run.Chain.Step();
                run.Accumulate(state, functionals);
                BassState replayed = step();
                if (k is not (0 or 25 or 49)) continue;
                Assert.Equal(BassMixing.Key(replayed.Terms), BassMixing.Key(state.Terms));
                Assert.Equal(BitConverter.DoubleToInt64Bits(replayed.LogTau), BitConverter.DoubleToInt64Bits(state.LogTau));
                double density = own();
                _out.WriteLine($"{(tempered ? "tempered" : "untempered")}, chain {c}, draw {k}: {functionals[lp]:R} against the chain's {density:R}");
                Assert.Equal(BitConverter.DoubleToInt64Bits(density), BitConverter.DoubleToInt64Bits(functionals[lp]));
            }
        }
    }

    [Fact]
    public void Functionals_AreNamed_AndGroupedByD10()
    {
        (BassConfig config, BassData data) = BassBitIdentityTests.Fixture(BassMixing.Friedman, "equal");
        BassEnsemble sampled = Ensemble(config, data);
        IReadOnlyList<BassMonitoredPoint> points = BassModel.Create(config, data).Resolved.MonitoredPoints;
        Assert.Equal(["M", "tau", "log_posterior", "sigma2", .. points.Select(Name)], sampled.FunctionalNames);
        Assert.Equal(sampled.MonitoredDim, sampled.FunctionalDim);
        Assert.Equal(4 + points.Count, sampled.FunctionalDim);
        AssertGroups(sampled, structural: ["M", "tau"]);

        BassEnsemble held = Ensemble(config with { Coefficients = BassCoefficientPrior.G(TauPrior.Held()) }, data);
        Assert.Equal(["M", "log_posterior", "sigma2", .. points.Select(Name)], held.FunctionalNames);
        AssertGroups(held, structural: ["M"]);

        BassData functional = BassRunTests.FunctionalData();
        BassEnsemble direct = Ensemble(BassRunTests.Config(), functional);
        IReadOnlyList<BassMonitoredPoint> at = BassModel.Create(BassRunTests.Config(), functional).Resolved.MonitoredPoints;
        Assert.Equal(at.Select(p => $"fitted[{p.Case},{p.Location}]"), direct.FunctionalNames.Skip(4));
        AssertGroups(direct, structural: ["M", "tau"]);

        static void AssertGroups(BassEnsemble ensemble, string[] structural)
        {
            IReadOnlyList<string> names = ensemble.FunctionalNames;
            IReadOnlyList<FunctionalGroup> groups = ensemble.Groups;
            Assert.Equal([BassDiagnostics.PredictiveGroup, BassDiagnostics.LevelGroup, BassDiagnostics.StructuralGroup], groups.Select(g => g.Name));
            Assert.Equal(["sigma2", .. names.Where(n => n.StartsWith("fitted[", StringComparison.Ordinal))], groups[0].Functionals.Select(f => names[f]));
            Assert.Equal(["log_posterior"], groups[1].Functionals.Select(f => names[f]));
            Assert.Equal(structural, groups[2].Functionals.Select(f => names[f]));
            Assert.Equal(Enumerable.Range(0, ensemble.FunctionalDim), groups.SelectMany(g => g.Functionals).Order());
        }
    }

    /// <summary>
    /// A fit carries its verdict, read by the shared rule on its diagnostics: the three groups in order, each gating the
    /// stop, against 1.01 without a target; each failing reading the group's diagnostic, and a quantity failing exactly
    /// when its R̂ misses the target or its ESS the gate. The log posterior's diagnostic and the mode-separation readings
    /// match those recomputed from the retained draws.
    /// </summary>
    [Fact]
    public void Fit_CarriesTheVerdict_AndTheModeSeparation()
    {
        (BassConfig config, BassData data) = BassBitIdentityTests.Fixture(BassMixing.Mixed, "equal");
        BassFit fit = Bass.Fit(config, data);
        BassDiagnostics d = fit.Diagnostics;
        ConvergenceVerdict v = d.Verdict;
        Assert.Equal(ChainEnsemble.ReferenceRHatTarget, v.Target);
        Assert.Equal(config.Sampling.MinEssPerChain, v.MinEssPerChain);
        Assert.Equal([BassDiagnostics.PredictiveGroup, BassDiagnostics.LevelGroup, BassDiagnostics.StructuralGroup], v.Groups.Select(g => g.Name));
        Assert.All(v.Groups, g => Assert.True(g.GatesStop));
        Assert.Equal(fit.Resolved.StoppingGroups, v.Groups.Select(g => g.Name));

        int chains = config.Sampling.Chains;
        bool Fails(ConvergenceDiagnostic q) => !(q.RHat.MeetsTarget(v.Target) && q.Ess.MeetsMinimum(v.MinEssPerChain * chains));
        (string Name, ConvergenceDiagnostic Diagnostic)[] predictive =
            [("sigma2", d.Sigma2), .. d.Fitted.Select(f => (Name(f.Point), f.Diagnostic))];
        (string, ConvergenceDiagnostic)[] level = [("log_posterior", d.LogPosterior)];
        (string, ConvergenceDiagnostic)[] structural = [("M", d.Count), ("tau", d.Tau!.Value)];
        foreach (((string Name, ConvergenceDiagnostic Diagnostic)[] quantities, string group) in
                 new[] { (predictive, BassDiagnostics.PredictiveGroup), (level, BassDiagnostics.LevelGroup), (structural, BassDiagnostics.StructuralGroup) })
        {
            GroupVerdict g = v.Group(group);
            var expected = quantities.Where(q => Fails(q.Diagnostic)).ToArray();
            Assert.Equal(expected.Select(q => q.Name), g.Failing.Select(r => r.Name));
            Assert.Equal(expected.Select(q => q.Diagnostic.RHat), g.Failing.Select(r => r.RHat));
            Assert.Equal(expected.Select(q => (EssDiagnostic?)q.Diagnostic.Ess), g.Failing.Select(r => r.Ess));
            Assert.Equal(expected.Length > 0 ? ConvergenceOutcome.NotConverged : ConvergenceOutcome.Converged, g.Outcome);
            _out.WriteLine($"{group}: {g.Outcome}; " + string.Join(", ", quantities.Select(q => $"{q.Name} R̂ {q.Diagnostic.RHat:F3} ESS {q.Diagnostic.Ess:F0}")));
        }
        Assert.Equal(v.Groups.Any(g => g.Outcome == ConvergenceOutcome.NotConverged) ? ConvergenceOutcome.NotConverged : ConvergenceOutcome.Converged, v.Outcome);

        // The log posterior and the mode separation, from the retained draws.
        BassModel model = BassModel.Create(config, data);
        double[][] logPosterior = BassMovePassTests.LogPosteriors(model, fit);
        ConvergenceDiagnostic reference = ChainDiagnostics.Diagnose(logPosterior);
        Assert.Equal(reference.RHat.Value, d.LogPosterior.RHat.Value, 9);
        Assert.Equal(reference.BulkEss.Value, d.LogPosterior.BulkEss.Value, 6);
        Assert.Equal(reference.TailEss.Value, d.LogPosterior.TailEss.Value, 6);

        double[] chainLevel = [.. logPosterior.Select(c => c.Average())];
        double[][] chainFitted = [.. Enumerable.Range(0, chains).Select(c => d.Fitted.Select(f => FittedTrace(model, fit, c, f.Point).Average()).ToArray())];
        double gap = chainLevel.Max() - chainLevel.Min();
        double widest = 0.0;
        for (int i = 0; i < chains; i++)
            for (int j = i + 1; j < chains; j++)
                widest = Math.Max(widest, Math.Sqrt(chainFitted[i].Zip(chainFitted[j], (a, b) => (a - b) * (a - b)).Sum()));
        _out.WriteLine($"mode separation: {d.ModeSeparation.LogPosteriorGap:F3} nats, fitted distance {d.ModeSeparation.FittedDistance:F4}");
        Assert.Equal(gap, d.ModeSeparation.LogPosteriorGap, 1e-8 * (1 + Math.Abs(chainLevel.Max())));
        Assert.Equal(widest, d.ModeSeparation.FittedDistance, 1e-10 * (1 + chainFitted.SelectMany(f => f).Max(Math.Abs)));

        BassFit one = Bass.Fit(config with { Sampling = config.Sampling with { Chains = 1 } }, data);
        Assert.Equal(ConvergenceOutcome.Undetermined, one.Diagnostics.Verdict.Outcome);
        Assert.True(double.IsNaN(one.Diagnostics.ModeSeparation.LogPosteriorGap) && double.IsNaN(one.Diagnostics.ModeSeparation.FittedDistance));
    }

    [Fact]
    public void StoppingScope_IsChecked_AndEchoed()
    {
        Assert.Equal(BassStoppingScope.All, new BassSampling().StoppingScope);
        Assert.Contains("\"StoppingScope\":\"All\"", new BassConfig().ToJson());
        var predictive = new BassConfig { Sampling = new BassSampling { StoppingScope = BassStoppingScope.Predictive } };
        string json = predictive.ToJson();
        Assert.Contains("\"StoppingScope\":\"Predictive\"", json);
        Assert.Equal(predictive, BassConfig.FromJson(json));
        Assert.NotEqual(Bass.Fingerprint(new BassConfig(), BassSweepTests.Data()), Bass.Fingerprint(predictive, BassSweepTests.Data()));
        Assert.StartsWith("Sampling.StoppingScope",
            Assert.Throws<ArgumentException>(new BassConfig { Sampling = new BassSampling { StoppingScope = (BassStoppingScope)7 } }.Validate).Message);

        BassConfig config = BassRunTests.Config(chains: 2, burn: 0, samples: 1);
        Assert.Equal([BassDiagnostics.PredictiveGroup, BassDiagnostics.LevelGroup, BassDiagnostics.StructuralGroup],
                     BassModel.Create(config, BassSweepTests.Data()).Resolved.StoppingGroups);
        BassConfig scoped = config with { Sampling = config.Sampling with { StoppingScope = BassStoppingScope.Predictive } };
        BassFit fit = Bass.Fit(scoped, BassSweepTests.Data());
        Assert.Equal([BassDiagnostics.PredictiveGroup], fit.Resolved.StoppingGroups);
        Assert.Equal([true, false, false], fit.Diagnostics.Verdict.Groups.Select(g => g.GatesStop));
    }

    /// <summary>
    /// The predictive scope against every group on a fit with an R̂ target: the draws are the same chains', so the predictive
    /// stop, reading a subset of the groups on the same rounds, ends no later; it ends Converged, and its verdict still
    /// reports the level and structural groups, which do not gate it.
    /// </summary>
    [Fact]
    public void PredictiveScope_StopsNoLaterThanAll_AndReportsEveryGroup()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig config = BassRunTests.Config(chains: 3, samples: 4000);
        config = config with { Sampling = config.Sampling with { RHatTarget = 1.05, BatchSize = 100, MinEssPerChain = 50 } };
        BassFit all = Bass.Fit(config, data);
        BassFit predictive = Bass.Fit(config with { Sampling = config.Sampling with { StoppingScope = BassStoppingScope.Predictive } }, data);
        _out.WriteLine($"every group: {all.Diagnostics.SamplesUsed} draws per chain, {all.Diagnostics.Verdict.Outcome}; " +
                       $"predictive: {predictive.Diagnostics.SamplesUsed}, {predictive.Diagnostics.Verdict.Outcome}");
        foreach (GroupVerdict g in predictive.Diagnostics.Verdict.Groups)
            _out.WriteLine($"  {g.Name}: {g.Outcome}, gates {g.GatesStop}; failing {string.Join(", ", g.Failing.Select(r => $"{r.Name} R̂ {r.RHat:F3} ESS {r.Ess:F0}"))}");

        Assert.True(predictive.Diagnostics.SamplesUsed <= all.Diagnostics.SamplesUsed);
        Assert.True(predictive.Diagnostics.StoppedEarly);
        Assert.Equal(ConvergenceOutcome.Converged, predictive.Diagnostics.Verdict.Outcome);
        Assert.Equal(1.05, predictive.Diagnostics.Verdict.Target);
        Assert.Equal(ConvergenceOutcome.Converged, predictive.Diagnostics.Verdict.Group(BassDiagnostics.PredictiveGroup).Outcome);
        Assert.False(predictive.Diagnostics.Verdict.Group(BassDiagnostics.LevelGroup).GatesStop);
        Assert.False(predictive.Diagnostics.Verdict.Group(BassDiagnostics.StructuralGroup).GatesStop);
        if (all.Diagnostics.StoppedEarly) Assert.Equal(ConvergenceOutcome.Converged, all.Diagnostics.Verdict.Outcome);
        else Assert.Equal(ConvergenceOutcome.NotConverged, all.Diagnostics.Verdict.Outcome);
        for (int c = 0; c < 3; c++)
            for (int k = 0; k < predictive.Diagnostics.SamplesUsed; k++) Assert.Equal(all.Draws[c][k].Sigma2, predictive.Draws[c][k].Sigma2);
    }

    /// <summary>
    /// S10's untempered Friedman, mixed and functional fits at their seeds (F-M2's A0 at k = 0), and Friedman tempered (A6):
    /// each quantity's R̂ is the one the pass recorded, and the verdict fails every group, as the pass read them (D10). The
    /// mixed fit's σ² passes inside its failing predictive group.
    /// </summary>
    [Fact]
    public void S10_Fits_FailEveryGroup()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassMovePassTests.Fixture[] fixtures = BassMovePassTests.Fixtures();
        // The pass's R̂ for M, σ², the log posterior and the three fitted values (local-move-pass-results.csv, k = 0).
        var cases = new (string Fixture, string Arm, double[] RHat, string[] PredictiveFailing)[]
        {
            ("friedman", "A0", [2.3477, 1.0831, 2.0500, 1.7465, 1.3173, 1.9615], ["sigma2", "f1", "f2", "f3"]),
            ("mixed", "A0", [1.5373, 1.0073, 1.6460, 1.2993, 1.0121, 1.0874], ["f1", "f2", "f3"]),
            ("functional", "A0", [2.5478, 1.0226, 2.3927, 1.1340, 1.4973, 1.1403], ["sigma2", "f1", "f2", "f3"]),
            ("friedman", "A6", [1.4901, 1.0768, 1.2637, 1.2348, 1.1886, 1.0487], ["sigma2", "f1", "f2", "f3"]),
        };
        foreach ((string name, string armName, double[] rhat, string[] predictiveFailing) in cases)
        {
            BassMovePassTests.Fixture f = fixtures.Single(x => x.Name == name);
            BassMovePassTests.Arm arm = BassMovePassTests.Arms.Single(a => a.Name == armName);
            BassFit fit = Bass.Fit(BassMovePassTests.Config(f, 0, arm), f.Data);
            BassDiagnostics d = fit.Diagnostics;
            ConvergenceVerdict v = d.Verdict;
            _out.WriteLine($"## {f.Label}, {arm.Name} ({arm.Label}): {v.Outcome}; mode separation {d.ModeSeparation.LogPosteriorGap:F1} nats, " +
                           $"fitted distance {d.ModeSeparation.FittedDistance:F3}");
            foreach (GroupVerdict g in v.Groups)
                _out.WriteLine($"- {g.Name}: {g.Outcome}; " + string.Join(", ", g.Failing.Select(r => $"{r.Name} R̂ {r.RHat:F3} ESS {r.Ess:F0}")));

            ConvergenceDiagnostic[] read = [d.Count, d.Sigma2, d.LogPosterior, .. d.Fitted.Select(x => x.Diagnostic)];
            for (int q = 0; q < rhat.Length; q++) Assert.Equal(rhat[q], read[q].RHat.Value, 1e-4);   // the CSV's, cut to four places

            Assert.Equal(ConvergenceOutcome.NotConverged, v.Outcome);
            Assert.All(v.Groups, g => Assert.Equal(ConvergenceOutcome.NotConverged, g.Outcome));
            string[] fittedNames = [.. d.Fitted.Select(x => Name(x.Point))];
            string Named(string q) => q switch { "f1" => fittedNames[0], "f2" => fittedNames[1], "f3" => fittedNames[2], _ => q };
            Assert.Equal(predictiveFailing.Select(Named), v.Group(BassDiagnostics.PredictiveGroup).Failing.Select(r => r.Name));
            Assert.Equal(["log_posterior"], v.Group(BassDiagnostics.LevelGroup).Failing.Select(r => r.Name));
            Assert.Contains("M", v.Group(BassDiagnostics.StructuralGroup).Failing.Select(r => r.Name));
            Assert.True(d.ModeSeparation.LogPosteriorGap > 0.0);
        }
    }

    // A monitored point's functional name: "fitted[i]" at case i, "fitted[i,r]" at location r too.
    private static string Name(BassMonitoredPoint p) => p.Location is int r ? $"fitted[{p.Case},{r}]" : $"fitted[{p.Case}]";

    private static BassEnsemble Ensemble(BassConfig config, BassData data)
        => new(BassModel.Create(config, data), Bass.Fingerprint(config, data), new BassRunOptions(), 0, config.Sampling.Samples);

    // Chain c's fitted value at `point` in each retained draw, summed from a₀ in the ids' order, as the run sums.
    private static double[] FittedTrace(BassModel model, BassFit fit, int c, BassMonitoredPoint point)
    {
        BassInputs inputs = model.Inputs;
        var column = new Dictionary<int, double>();
        double At(int id)
        {
            if (column.TryGetValue(id, out double value)) return value;
            var a = new double[inputs.Cases.Cases];
            model.Basis.CaseColumn(fit.Terms[id], inputs.Cases, a);
            value = a[point.Case];
            if (point.Location is int r)
            {
                var h = new double[inputs.Outputs!.Locations];
                model.Basis.OutputColumn(fit.Terms[id], inputs.Outputs, h);
                value *= h[r];
            }
            return column[id] = value;
        }
        return [.. fit.Draws[c].Select(d => d.Structure.TermIds.Select((id, k) => d.Coefficients[k + 1] * At(id)).Aggregate(d.Coefficients[0], (f, t) => f + t))];
    }
}
