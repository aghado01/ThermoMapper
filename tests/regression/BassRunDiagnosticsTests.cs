using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's run diagnostics (stage 4's S7): the totals are the chains' sums, over the sampling phase, where every sweep is
/// a move attempt or an empty birth; the support tally's attempts are the births and changes, its refusals the support
/// refusals, and its rate rises with the order; the mass at the cap is 0 when the cap cannot bind and the share of draws
/// there when it does; the fitted values' diagnostics are those of each draw's function recomputed from the basis; the
/// R̂ stop ends a run early; fixed g has no τ diagnostic; and every chain starts in F.
/// </summary>
public sealed class BassRunDiagnosticsTests
{
    private readonly ITestOutputHelper _out;
    public BassRunDiagnosticsTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData("continuous")]
    [InlineData("categorical")]
    [InlineData("functional")]
    public void Diagnostics_AreTheChainsSums(string fixture)
    {
        BassConfig config = BassRunTests.Config(chains: 3, burn: 300, samples: 400, thin: 2);
        BassData data = fixture switch
        {
            "categorical" => BassPriorSamplerTests.MixedData(),
            "functional" => BassRunTests.FunctionalData(),
            _ => BassSweepTests.Data(large: true),
        };
        BassFit fit = Bass.Fit(config, data);
        BassDiagnostics d = fit.Diagnostics;
        Assert.Equal(3, d.Chains.Count);
        foreach (BassRefusal reason in Enum.GetValues<BassRefusal>())
            Assert.Equal(d.Chains.Sum(c => c.Refusals.Count(reason)), d.Refusals.Count(reason));
        Assert.Equal(d.Chains.Sum(c => c.EmptyBirths), d.EmptyBirths);
        foreach (BassOrderSupport e in d.Support.Entries)
        {
            Assert.Equal(d.Chains.Sum(c => c.Support.Attempts(e.Part, e.Order)), e.Attempts);
            Assert.Equal(d.Chains.Sum(c => c.Support.Refused(e.Part, e.Order)), e.Refused);
        }
        for (int i = 0; i < d.Moves.Count; i++)
        {
            Assert.Equal(d.Chains.Sum(c => c.Moves[i].Attempts), d.Moves[i].Attempts);
            Assert.Equal(d.Chains.Sum(c => c.Moves[i].Accepted), d.Moves[i].Accepted);
        }
        Assert.Equal((double)d.Moves.Sum(m => m.Accepted) / d.Moves.Sum(m => m.Attempts), d.AcceptanceRate, 12);
        Assert.Equal(d.Chains.Average(c => c.MassAtCap), d.MassAtCap, 12);   // equal draws per chain

        BassStructurePrior prior = BassModel.Create(config, data).Prior;
        foreach (BassChainAccount c in d.Chains)
        {
            // The sampling phase: under the equal schedule each sweep attempts a move or draws an empty birth.
            long sweeps = (long)config.Sampling.Samples * config.Sampling.Thin;
            Assert.Equal(sweeps, c.Moves.Sum(m => m.Attempts) + c.EmptyBirths);
            // Each birth and change counts its candidate once at each part's order; the refusals are the support ones.
            long candidates = c.Moves.Where(m => m.Key is "birth" or "change").Sum(m => m.Attempts);
            foreach (BassPart part in Enum.GetValues<BassPart>())
            {
                if (!prior.Has(part)) continue;
                IEnumerable<BassOrderSupport> orders = c.Support.Entries.Where(e => e.Part == part);
                Assert.Equal(candidates, orders.Sum(e => e.Attempts));
                Assert.Equal(c.Refusals.Count(BassRefusal.Support), orders.Sum(e => e.Refused));
            }
            Assert.True(c.Start.Attempts >= 1);
        }
        _out.WriteLine($"{fixture}: moves {string.Join(", ", d.Moves.Select(m => $"{m.Key} {m.Accepted}/{m.Attempts}"))}; " +
                       $"refusals {string.Join(", ", d.Refusals.Entries.Where(e => e.Count > 0).Select(e => $"{e.Category} {e.Count}"))}; " +
                       $"empty births {d.EmptyBirths}; support {string.Join(", ", d.Support.Entries.Where(e => e.Attempts > 0).Select(e => $"{e.Part} {e.Order}: {e.RefusalRate:F3}"))}");
        if (fixture == "continuous") Assert.Equal(0, d.EmptyBirths);
        else Assert.True(d.EmptyBirths > 0, "two parts: a birth's draw can come out the intercept");
    }

    [Fact]
    public void SupportRefusals_RiseWithTheOrder()
    {
        BassFit fit = Bass.Fit(BassRunTests.Config(chains: 2, burn: 300, samples: 1500), BassSweepTests.Data(large: true));
        BassSupportTally s = fit.Diagnostics.Support;
        Assert.Equal([1, 2], s.Entries.Where(e => e.Part == BassPart.Continuous && e.Attempts > 0).Select(e => e.Order));
        double first = s.Entries.Single(e => e.Order == 1).RefusalRate, second = s.Entries.Single(e => e.Order == 2).RefusalRate;
        _out.WriteLine($"support refusal rate: order 1 {first:F3}, order 2 {second:F3}");
        // A product of hinges is nonzero only where each factor is, so an order-2 candidate fails at least as often.
        Assert.True(second > first, $"order 2 {second:F3} against order 1 {first:F3}");
    }

    [Fact]
    public void MassAtCap_ReadsTheDraws()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassFit loose = Bass.Fit(BassRunTests.Config(maxBasis: 60), data);   // the cap cannot bind
        Assert.Equal(0.0, loose.Diagnostics.MassAtCap);
        Assert.All(loose.Diagnostics.Chains, c => Assert.Equal(0.0, c.MassAtCap));
        Assert.Equal(0, loose.Diagnostics.Refusals.Count(BassRefusal.Cap));

        BassFit tight = Bass.Fit(BassRunTests.Config(maxBasis: 1), data);   // the data want more than one term
        double atCap = tight.Draws.Sum(c => c.Count(d => d.Structure.TermIds.Length == 1)) / (double)tight.Draws.Sum(c => c.Count);
        _out.WriteLine($"mass at a cap of one term: {tight.Diagnostics.MassAtCap:F3}");
        Assert.Equal(atCap, tight.Diagnostics.MassAtCap, 12);
        Assert.True(atCap > 0.9);
    }

    [Theory]
    [InlineData("continuous")]
    [InlineData("functional")]
    public void Fitted_AreTheDrawsFunctionsAtTheDefaultPoints(string fixture)
    {
        BassConfig config = BassRunTests.Config(chains: 3, samples: 300);
        BassData data = fixture == "functional" ? BassRunTests.FunctionalData() : BassSweepTests.Data(large: true);
        BassFit fit = Bass.Fit(config, data);
        BassMonitoredPoint[] expected = fixture == "functional"
            ? [new(6, 3), new(12, 3), new(18, 3)]     // 25 cases: positions 6, 12 and 18; seven locations: the middle, 3
            : [new(19, null), new(37, null), new(56, null)];   // 75 cases: 18.5, 37 and 55.5, rounded away from zero
        Assert.Equal(expected, fit.Diagnostics.Fitted.Select(f => f.Point));

        // Each draw's f = a₀ + Σⱼ aⱼBⱼ at the point, from the basis's full training columns.
        var model = BassModel.Create(config, data);
        BassBasis basis = model.Basis;
        BassInputs inputs = model.Inputs;
        var columns = new Dictionary<int, double[]>();
        double ColumnAt(int id, BassMonitoredPoint point)
        {
            if (!columns.TryGetValue(id, out double[]? a))
            {
                a = new double[inputs.Cases.Cases];
                basis.CaseColumn(fit.Terms[id], inputs.Cases, a);
                columns[id] = a;
            }
            double value = a[point.Case];
            if (point.Location is int r)
            {
                var h = new double[inputs.Outputs!.Locations];
                basis.OutputColumn(fit.Terms[id], inputs.Outputs, h);
                value *= h[r];
            }
            return value;
        }
        for (int p = 0; p < expected.Length; p++)
        {
            // Summed from a₀ in the ids' order, as the run sums: the tail R̂ folds the draws about their median, where an
            // order of rounding can swap two near-ties.
            double[][] traces = [.. fit.Draws.Select(chain => chain.Select(d => d.Structure.TermIds
                .Select((id, k) => d.Coefficients[k + 1] * ColumnAt(id, expected[p])).Aggregate(d.Coefficients[0], (f, t) => f + t)).ToArray())];
            ConvergenceDiagnostic reference = ChainDiagnostics.Diagnose(traces), actual = fit.Diagnostics.Fitted[p].Diagnostic;
            _out.WriteLine($"{expected[p]}: R̂ {actual.RHat}, bulk ESS {actual.BulkEss}, tail ESS {actual.TailEss}");
            Assert.Equal(reference.RHat.State, actual.RHat.State);
            Assert.Equal(reference.RHat.Value, actual.RHat.Value, 9);
            Assert.Equal(reference.BulkEss.Value, actual.BulkEss.Value, 6);
            Assert.Equal(reference.TailEss.Value, actual.TailEss.Value, 6);
        }
        // M, τ and σ² likewise, from the draws.
        Assert.Equal(ChainDiagnostics.Diagnose([.. fit.Draws.Select(c => c.Select(d => (double)d.Structure.TermIds.Length).ToArray())]), fit.Diagnostics.Count);
        Assert.Equal(ChainDiagnostics.Diagnose([.. fit.Draws.Select(c => c.Select(d => d.Tau).ToArray())]), fit.Diagnostics.Tau);
        Assert.Equal(ChainDiagnostics.Diagnose([.. fit.Draws.Select(c => c.Select(d => d.Sigma2).ToArray())]), fit.Diagnostics.Sigma2);
    }

    [Fact]
    public void RHatStop_EndsTheRunEarly()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig config = BassRunTests.Config(chains: 3, samples: 4000);
        BassFit full = Bass.Fit(config, data);
        Assert.False(full.Diagnostics.StoppedEarly);
        Assert.Equal(4000, full.Diagnostics.SamplesUsed);

        BassFit stopped = Bass.Fit(config with { Sampling = config.Sampling with { RHatTarget = 1.2, BatchSize = 100, MinEssPerChain = 20 } }, data);
        _out.WriteLine($"stopped after {stopped.Diagnostics.SamplesUsed} draws per chain");
        Assert.True(stopped.Diagnostics.StoppedEarly);
        Assert.True(stopped.Diagnostics.SamplesUsed < 4000);
        Assert.Equal(0, stopped.Diagnostics.SamplesUsed % 100);
        Assert.All(stopped.Draws, chain => Assert.Equal(stopped.Diagnostics.SamplesUsed, chain.Count));
        Assert.True(stopped.Diagnostics.Count.RHat.MeetsTarget(1.2));
        // The stopped run's draws are the full run's first ones: the stop only ends sampling.
        for (int c = 0; c < 3; c++)
            for (int d = 0; d < stopped.Diagnostics.SamplesUsed; d++) Assert.Equal(full.Draws[c][d].Sigma2, stopped.Draws[c][d].Sigma2);
    }

    [Fact]
    public void DefaultPoints_SpreadThroughTheCases()
    {
        BassResolved resolved = Bass.Fit(BassRunTests.Config(chains: 1, burn: 0, samples: 1), BassSweepTests.Data()).Resolved;
        static int[] Cases(BassResolved r) => [.. BassDiagnostics.DefaultPoints(r, false).Select(p => p.Case)];
        Assert.Equal([0], Cases(resolved with { Cases = 1 }));
        Assert.Equal([0, 1], Cases(resolved with { Cases = 2 }));
        Assert.Equal([1, 2], Cases(resolved with { Cases = 3 }));   // ½ and 1½ round away from zero; 1 is the middle
        Assert.Equal([4, 9, 13], Cases(resolved));                  // 18 cases
        Assert.All(BassDiagnostics.DefaultPoints(resolved, false), p => Assert.Null(p.Location));
        Assert.All(BassDiagnostics.DefaultPoints(resolved with { Locations = 6 }, true), p => Assert.Equal(2, p.Location));
        Assert.All(BassDiagnostics.DefaultPoints(resolved with { Locations = 1 }, true), p => Assert.Equal(0, p.Location));
    }

    [Fact]
    public void SupportTally_MergesOnlyItsOwnOrders()
    {
        var one = BassModel.Create(BassRunTests.Config(), BassSweepTests.Data());
        var mixed = BassModel.Create(BassRunTests.Config(), BassPriorSamplerTests.MixedData());
        var tally = new BassSupportTally(one.Prior);
        tally.Add(new BassTerm([new HingeFactor(0, 1, 0), new HingeFactor(1, -1, 2)]), true);
        tally.MergeFrom(tally.Copy());
        Assert.Equal(2, tally.Attempts(BassPart.Continuous, 2));
        Assert.Equal(2, tally.Refused(BassPart.Continuous, 2));
        Assert.Equal([1, 2], tally.Entries.Select(e => e.Order));   // one part: no order 0
        Assert.Throws<ArgumentException>(() => tally.MergeFrom(new BassSupportTally(mixed.Prior)));
        tally.Clear();
        Assert.All(tally.Entries, e => Assert.True(e.Attempts == 0 && double.IsNaN(e.RefusalRate)));
    }

    [Fact]
    public void FixedG_HasNoTauDiagnostic()
    {
        BassConfig config = BassRunTests.Config() with { Coefficients = new FixedG() };
        BassFit fit = Bass.Fit(config, BassSweepTests.Data(large: true));
        Assert.Null(fit.Diagnostics.Tau);
        Assert.Equal(3, fit.Diagnostics.Fitted.Count);
        Assert.All(fit.Diagnostics.Chains, c => Assert.Equal(1.0 / 75, c.Start.Tau));
    }

    [Fact]
    public void Starts_LieInF()
    {
        // Under Jeffreys on data a single hinge fits exactly, F excludes the exact fits; no chain starts at one.
        BassConfig config = BassSweepTests.Config() with
        {
            Sigma2 = Sigma2Prior.Jeffreys,
            Sampling = new BassSampling { Chains = 12, MasterSeed = 5, BurnIn = 0, Samples = 1 },
        };
        BassData data = BassSweepTests.Data(exact: true);
        BassFit fit = Bass.Fit(config, data);
        var model = BassModel.Create(config, data);
        var f = BassTermSpace.F(model).Select(s => BassPriorSamplerTests.Key(s.Terms)).ToHashSet();
        foreach (BassChainAccount c in fit.Diagnostics.Chains)
        {
            Assert.Contains(BassPriorSamplerTests.Key(c.Start.Terms), f);
            BassDesign design = model.NewDesign();
            Assert.Null(design.Load(c.Start.Terms));
            Assert.Null(model.Target.Refusal(BassStructureFit.Of(design, 0.0)));
        }
        _out.WriteLine($"starts: {string.Join("; ", fit.Diagnostics.Chains.Select(c => $"M = {c.Start.Terms.Count} after {c.Start.Attempts}"))}");
        Assert.True(fit.Diagnostics.Chains.Select(c => BassPriorSamplerTests.Key(c.Start.Terms)).Distinct().Count() > 1, "the starts are dispersed");
    }
}
