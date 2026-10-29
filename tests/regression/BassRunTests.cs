using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's run (stage 4's S7): <see cref="Bass.Fit"/> on the chain ensemble. Each chain, rebuilt by hand from the seed
/// tree's streams, gives the fit's draws exactly, through the term tables' merge; the result does not depend on the
/// parallelism bound; thinning keeps the trajectory; every option the config offers runs end to end; the run's options are
/// checked, cancellation stops it and progress counts every sweep. The posterior itself is checked against the oracle in
/// <c>BassFitParityTests</c>.
/// </summary>
public sealed class BassRunTests
{
    private readonly ITestOutputHelper _out;
    public BassRunTests(ITestOutputHelper output) => _out = output;

    internal static BassConfig Config(int chains = 3, int burn = 200, int samples = 150, int thin = 1, int seed = 11, int maxBasis = 4)
        => new()
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = maxBasis },
            Sampling = new BassSampling { Chains = chains, MasterSeed = seed, BurnIn = burn, Samples = samples, Thin = thin },
        };

    // A direct functional response: x₀ and x₁ on a five-value grid (n = 25), t on seven locations (N = 175).
    internal static BassData FunctionalData(int seed = 13)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] t = [.. Enumerable.Range(0, 7).Select(r => r / 6.0)];
        var x0 = new double[25];
        var x1 = new double[25];
        var y = new double[25, 7];
        for (int i = 0; i < 25; i++)
        {
            x0[i] = (i / 5) / 4.0;
            x1[i] = (i % 5) / 4.0;
            for (int r = 0; r < 7; r++)
                y[i, r] = 1 + 2 * x0[i] + 1.5 * t[r] * x1[i] + 2 * Math.Max(x0[i] - 0.5, 0) * Math.Max(t[r] - 0.5, 0) + 0.2 * (rng.NextDouble() - 0.5);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)],
                            new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void Fit_RebuildsFromTheSeedTreeByHand(int thin)
    {
        BassConfig config = Config(thin: thin);
        BassData data = BassSweepTests.Data(large: true);
        BassFit fit = Bass.Fit(config, data);
        BassSampling s = config.Sampling;
        Assert.Equal(s.Chains, fit.Draws.Count);
        Assert.Null(fit.DrawWeights);
        Assert.Null(fit.Checkpoint);
        Assert.Equal(new RunIdentity(s.MasterSeed, s.Chains, Bass.Fingerprint(config, data)), fit.Identity);

        // The shell's seed tree: chain streams first, then readout streams; each chain seed fans out to the start, the
        // moves and τ's block, and its readout seed drives composition.
        var model = BassModel.Create(config, data);
        int[] all = SeedTree.Derive(s.MasterSeed, 2 * s.Chains);
        for (int c = 0; c < s.Chains; c++)
        {
            int[] seeds = SeedTree.Derive(all[c], 3);
            BassState start = new BassPriorSampler(model).Draw(new Xoshiro256PlusPlus(seeds[0])).State;
            var chain = new BassChain(model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(seeds[1]), new Xoshiro256PlusPlus(seeds[2]));
            var composition = new Xoshiro256PlusPlus(all[s.Chains + c]);
            var table = new BassTermTable();
            for (int i = 0; i < s.BurnIn; i++) chain.Step();
            Assert.Equal(s.Samples, fit.Draws[c].Count);
            for (int d = 0; d < s.Samples; d++)
            {
                for (int t = 0; t < thin; t++) chain.Step();
                BassDraw expected = chain.Compose(composition, table), actual = fit.Draws[c][d];
                Assert.Equal(expected.Tau, actual.Tau);
                Assert.Equal(expected.Sigma2, actual.Sigma2);
                Assert.Equal(expected.Lambda, actual.Lambda);
                // The same terms with the same coefficients, whatever their ids in the merged table.
                Assert.Equal(Paired(expected, table), Paired(actual, fit.Terms));
                Assert.Equal(expected.Coefficients[0], actual.Coefficients[0]);
                Assert.True(actual.Structure.TermIds.SequenceEqual(actual.Structure.TermIds.Order()), "ids are sorted");
            }
        }
        Assert.Equal(fit.Terms.Count, Enumerable.Range(0, fit.Terms.Count).Select(id => fit.Terms[id]).Distinct().Count());
    }

    [Fact]
    public void Fit_DoesNotDependOnTheParallelismBound()
    {
        BassConfig config = Config(chains: 4, samples: 120);
        BassData data = BassSweepTests.Data(large: true);
        BassFit reference = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        foreach (int degree in new[] { 2, 3, -1 })
            AssertSameFit(reference, Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = degree }));
    }

    [Fact]
    public void Fit_DoesNotDependOnTheRunsOtherOptions()
    {
        BassConfig config = Config(chains: 3, samples: 100);
        BassData data = BassSweepTests.Data(large: true);
        BassFit plain = Bass.Fit(config, data);
        Assert.Null(plain.ResidualDependence);
        var options = new BassRunOptions { Checkpoint = true, Progress = new Inline(_ => { }), ProgressInterval = 7 };
        BassFit other = Bass.Fit(config, data, options);
        Assert.NotNull(other.Checkpoint);
        AssertSameFit(plain, other);
    }

    [Fact]
    public void Thinning_KeepsTheTrajectory()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassFit thinned = Bass.Fit(Config(chains: 2, samples: 60, thin: 3), data);
        BassFit every = Bass.Fit(Config(chains: 2, samples: 180, thin: 1), data);
        for (int c = 0; c < 2; c++)
            for (int d = 0; d < 60; d++)
            {
                BassDraw a = thinned.Draws[c][d], b = every.Draws[c][3 * d + 2];
                Assert.Equal(a.Tau, b.Tau);
                Assert.Equal(Structure(a, thinned.Terms), Structure(b, every.Terms));
            }
    }

    [Theory]
    [InlineData("defaults")]
    [InlineData("flat")]
    [InlineData("zero")]
    [InlineData("jeffreys")]
    [InlineData("inverse gamma")]
    [InlineData("fixed g")]
    [InlineData("fixed g, declared")]
    [InlineData("cancelling")]
    [InlineData("categorical")]
    [InlineData("categorical, uniform subsets")]
    [InlineData("functional")]
    [InlineData("column mass")]
    [InlineData("quadratic")]
    public void Fit_RunsEveryOption(string option)
    {
        BassConfig config = Config(chains: 2, burn: 150, samples: 100);
        BassData data = BassSweepTests.Data(large: true);
        switch (option)
        {
            case "flat": config = config with { Intercept = GPriorCentre.Flat }; break;
            case "zero": config = config with { Intercept = GPriorCentre.Zero }; break;
            case "jeffreys": config = config with { Sigma2 = Sigma2Prior.Jeffreys }; break;
            case "inverse gamma": config = config with { Sigma2 = Sigma2Prior.InverseGamma(2, 0.5) }; break;
            case "fixed g": config = config with { Coefficients = new FixedG() }; break;
            case "fixed g, declared": config = config with { Coefficients = new FixedG(0.05) }; break;
            case "cancelling": config = config with { Sampling = config.Sampling with { Moves = new PriorCancelling(0.4) } }; break;
            case "categorical":
            case "categorical, uniform subsets":
                data = BassPriorSamplerTests.MixedData();
                if (option.EndsWith("subsets")) config = config with { Structure = config.Structure with { CategoricalSubsets = CategoricalSubsetMeasure.Uniform } };
                break;
            case "functional": data = FunctionalData(); break;
            case "column mass": config = config with { Structure = config.Structure with { Support = new SupportRule(SupportKind.ColumnMass) } }; break;
            case "quadratic": config = config with { Structure = config.Structure with { Degree = 2 } }; break;
        }
        BassFit fit = Bass.Fit(config, data);
        var model = BassModel.Create(config, data);
        double? fixedTau = model.Target.FixedTau;
        int moved = 0;
        foreach (DrawStore<BassDraw> chain in fit.Draws)
        {
            Assert.Equal(100, chain.Count);
            foreach (BassDraw d in chain)
            {
                int m = d.Structure.TermIds.Length;
                Assert.InRange(m, 0, fit.Resolved.MaxBasis);
                Assert.Equal(m + 1, d.Coefficients.Length);
                Assert.True(d.Sigma2 > 0 && d.Lambda > 0 && d.Tau > 0 && d.Coefficients.All(double.IsFinite));
                if (fixedTau is { } tau) Assert.Equal(tau, d.Tau);
                // Every retained structure lies in F.
                BassTerm[] terms = [.. d.Structure.TermIds.Select(id => fit.Terms[id])];
                BassDesign design = model.NewDesign();
                Assert.Null(design.Load(terms));
                Assert.Null(model.Target.Refusal(BassStructureFit.Of(design, 0.0)));
                if (m > 0) moved++;
            }
        }
        _out.WriteLine($"{option}: {fit.Terms.Count} terms; acceptance {fit.Diagnostics.AcceptanceRate:F3}; draws with terms {moved}");
        Assert.True(moved > 0, "the chains visit structures with terms");
    }

    [Fact]
    public void Fit_ChecksItsOptions()
    {
        BassData data = BassSweepTests.Data();
        BassConfig config = Config(chains: 2, burn: 10, samples: 10);
        Assert.Throws<ArgumentOutOfRangeException>(() => Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = -2 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => Bass.Fit(config, data, new BassRunOptions { ProgressInterval = 0 }));
        var e = Assert.Throws<ArgumentException>(() => Bass.Fit(config with { Sampling = config.Sampling with { Chains = 1, RHatTarget = 1.1 } }, data));
        Assert.Contains("Sampling.Chains", e.Message);
        Assert.Throws<ArgumentNullException>(() => Bass.Fit(null!, data));
        Assert.Throws<ArgumentNullException>(() => Bass.Fit(config, null!));
    }

    [Fact]
    public void Fit_ThrowsTheChainsSharedFailure()
    {
        // F is the intercept alone, about 6·10⁻⁹ of the proposal (BassPriorSamplerTests' attempt-limit fixture), so every
        // chain's start refuses; the caller sees the sampler's refusal, not the fan-out's wrapper.
        BassConfig config = BassSweepTests.Config() with
        {
            Count = new BassCountPrior(1000, 1e-3),
            Structure = BassSweepTests.Config().Structure with { Support = new SupportRule(Minimum: 18) },
            Sampling = new BassSampling { Chains = 2, BurnIn = 1, Samples = 1 },
        };
        var e = Assert.Throws<InvalidOperationException>(() => Bass.Fit(config, BassSweepTests.Data()));
        Assert.Contains("No collection in F", e.Message);
    }

    [Fact]
    public void Fit_StopsOnCancellation()
    {
        BassData data = BassSweepTests.Data();
        using var before = new CancellationTokenSource();
        before.Cancel();
        Assert.Throws<OperationCanceledException>(() => Bass.Fit(Config(), data, new BassRunOptions { Cancellation = before.Token }));

        // Cancelled mid-run, at the first progress report, from a chain's thread.
        using var during = new CancellationTokenSource();
        var options = new BassRunOptions { Cancellation = during.Token, Progress = new Inline(_ => during.Cancel()), ProgressInterval = 50 };
        Assert.Throws<OperationCanceledException>(() => Bass.Fit(Config(burn: 5_000), data, options));
    }

    [Fact]
    public void Progress_CountsEverySweep()
    {
        var reports = new List<BassProgress>();
        var options = new BassRunOptions { Progress = new Inline(p => { lock (reports) reports.Add(p); }), ProgressInterval = 100 };
        BassConfig config = Config(chains: 3, burn: 130, samples: 70, thin: 2);
        Bass.Fit(config, BassSweepTests.Data(), options);
        long planned = 3 * (130 + 70 * 2);
        Assert.All(reports, r => Assert.Equal(planned, r.PlannedSweeps));
        Assert.Equal(planned, reports[^1].Sweeps);   // the last, from the caller once the chains are done
        Assert.Equal(Enumerable.Range(1, (int)(planned / 100)).Select(k => 100L * k), reports.Take(reports.Count - 1).Select(r => r.Sweeps).Order());
    }

    // The draw's terms by value with their coefficients.
    private static List<(string Term, double Coefficient)> Paired(BassDraw draw, BassTermTable terms)
        => [.. draw.Structure.TermIds.Select((id, k) => (terms[id].ToString(), draw.Coefficients[k + 1])).OrderBy(p => p.Item1, StringComparer.Ordinal)];

    private static string Structure(BassDraw draw, BassTermTable terms) => BassPriorSamplerTests.Key(draw.Structure.TermIds.Select(id => terms[id]));

    internal static void AssertSameFit(BassFit a, BassFit b)
    {
        Assert.Equal(a.Terms.Count, b.Terms.Count);
        for (int id = 0; id < a.Terms.Count; id++) Assert.Equal(a.Terms[id], b.Terms[id]);
        Assert.Equal(a.Draws.Count, b.Draws.Count);
        for (int c = 0; c < a.Draws.Count; c++)
        {
            Assert.Equal(a.Draws[c].Count, b.Draws[c].Count);
            for (int d = 0; d < a.Draws[c].Count; d++)
            {
                BassDraw x = a.Draws[c][d], y = b.Draws[c][d];
                Assert.Equal(x.Structure.TermIds, y.Structure.TermIds);
                Assert.Equal(x.Tau, y.Tau);
                Assert.Equal(x.Sigma2, y.Sigma2);
                Assert.Equal(x.Coefficients, y.Coefficients);
                Assert.Equal(x.Lambda, y.Lambda);
            }
        }
        AssertSameDiagnostics(a.Diagnostics, b.Diagnostics);
    }

    internal static void AssertSameDiagnostics(BassDiagnostics a, BassDiagnostics b)
    {
        Assert.Equal(a.Count, b.Count);
        Assert.Equal(a.Tau, b.Tau);
        Assert.Equal(a.LogPosterior, b.LogPosterior);
        Assert.Equal(a.Sigma2, b.Sigma2);
        Assert.Equal(a.Fitted, b.Fitted);
        Assert.Equal(a.ModeSeparation, b.ModeSeparation);
        VerdictAssert.Same(a.Verdict, b.Verdict);
        Assert.Equal(a.Moves, b.Moves);
        Assert.Equal(a.AcceptanceRate, b.AcceptanceRate);
        Assert.Equal(a.Refusals.Entries, b.Refusals.Entries);
        Assert.Equal(a.EmptyBirths, b.EmptyBirths);
        Assert.Equal(a.Support.Entries, b.Support.Entries);
        Assert.Equal(a.MassAtCap, b.MassAtCap);
        Assert.Equal(a.SamplesUsed, b.SamplesUsed);
        Assert.Equal(a.StoppedEarly, b.StoppedEarly);
        Assert.Equal(a.Chains.Count, b.Chains.Count);
        for (int c = 0; c < a.Chains.Count; c++)
        {
            BassChainAccount x = a.Chains[c], y = b.Chains[c];
            Assert.Equal(x.Start.Terms, y.Start.Terms);
            Assert.Equal(x.Start.Tau, y.Start.Tau);
            Assert.Equal(x.Start.Attempts, y.Start.Attempts);
            Assert.Equal(x.Moves, y.Moves);
            Assert.Equal(x.Refusals.Entries, y.Refusals.Entries);
            Assert.Equal(x.EmptyBirths, y.EmptyBirths);
            Assert.Equal(x.Support.Entries, y.Support.Entries);
            Assert.Equal(x.MassAtCap, y.MassAtCap);
        }
    }

    // Reports on the calling thread, as the run's callers see them.
    private sealed class Inline(Action<BassProgress> report) : IProgress<BassProgress>
    {
        public void Report(BassProgress value) => report(value);
    }
}
