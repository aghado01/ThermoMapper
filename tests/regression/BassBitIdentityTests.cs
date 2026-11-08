using System;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The bit-identity baseline of F-M2's Phase 1 (the local-move options brief, Step 0): short fits hashed at the commit
/// before any option was built, so that every option left off keeps today's draws bit for bit.
/// <list type="bullet">
/// <item>S10's Friedman, mixed and functional fixtures, untempered, on both move schedules: each retained draw's M, τ, σ²,
/// λ, coefficients and terms, chain by chain.</item>
/// <item>Friedman tempered with a self-sizing ladder: the same, and each chain's frozen rungs.</item>
/// <item>A tempered BARS fit, which calls the same ladder: its fit and sd on the grid, and the ladder's readout.</item>
/// </list>
/// A hash moves only when a change moves the draws; regenerating one is a deliberate, called-out act.
/// <para><b>Re-baselined 2026-11-04 (the efficiency pass's lever 1).</b> In a plain chain of a scalar fit of at most 16 continuous predictors and
/// at least 500 cases, which the untempered Friedman fixture is (the mixed, functional and tempered ones are not: a ladder's replicas, handed
/// states over, keep the dots), a one-way hinge's cross products with the held columns come from prefix tables. They differ from the dots by
/// rounding, at most 4·K·ε of the columns' norms, so the draws differ in their last digits and two hashes moved: untempered equal
/// ee7d60a2738b56d8 → 07d9e2ca9822490c, untempered cancelling 414e30ea867f9058 → 9e1e7fb5e46ca40b. With <c>BASS_TABLES=off</c> in the
/// environment all eight former values hold, which is the check that the pass's other levers left every bit; the structure sequences with
/// tables on and off are identical (<c>BassPrefixTableTests</c>).</para>
/// </summary>
[Trait("Tier", "Pinned")]
public sealed class BassBitIdentityTests
{
    private readonly ITestOutputHelper _out;
    public BassBitIdentityTests(ITestOutputHelper output) => _out = output;

    private const int BurnIn = 300, Samples = 100, Chains = 2;

    // S10's fixture, short, on the schedule; tempered by a self-sizing ladder when asked.
    internal static (BassConfig Config, BassData Data) Fixture(string name, string schedule, bool tempered = false)
    {
        BassMixing.Fixture f = BassMixing.Get(name, tempered ? new TemperingSpec() : null);
        MoveSchedule moves = schedule == "equal" ? new EqualOverAvailable() : new PriorCancelling();
        BassConfig config = f.Config with
        {
            Sampling = f.Config.Sampling with
            {
                Chains = Chains, BurnIn = tempered ? 255 : BurnIn, Samples = tempered ? 50 : Samples, Moves = moves,
            },
        };
        return (config, f.Data);
    }

    // FNV-1a over every retained draw, chain by chain: M, τ, σ², λ, the coefficients, and each term by value; each
    // chain's move counts, refusals and empty births; and a tempered fit's frozen rungs.
    internal static string Hash(BassFit fit)
    {
        var hash = new RunFingerprint();
        foreach (DrawStore<BassDraw> chain in fit.Draws)
        {
            hash.Add(chain.Count);
            foreach (BassDraw d in chain)
            {
                hash.Add(d.Structure.TermIds.Length).Add(d.Tau).Add(d.Sigma2).Add(d.Lambda).Add(d.Coefficients);
                foreach (int id in d.Structure.TermIds) hash.Add(fit.Terms[id].ToString());
            }
        }
        foreach (BassChainAccount account in fit.Diagnostics.Chains)
        {
            foreach (MoveStat m in account.Moves) hash.Add(m.Key).Add((ulong)m.Attempts).Add((ulong)m.Accepted);
            foreach (BassRefusal reason in Enum.GetValues<BassRefusal>()) hash.Add((ulong)account.Refusals.Count(reason));
            hash.Add((ulong)account.EmptyBirths);
        }
        if (fit.Diagnostics.Tempering is { } t)
            foreach (double[] rungs in t.Ladders) hash.Add(rungs);
        return hash.ToString();
    }

    [Theory]
    [InlineData(BassMixing.Friedman, "equal", "07d9e2ca9822490c")]
    [InlineData(BassMixing.Friedman, "cancelling", "9e1e7fb5e46ca40b")]
    [InlineData(BassMixing.Mixed, "equal", "137dad3cbeb33acb")]
    [InlineData(BassMixing.Mixed, "cancelling", "72fda364ee08f1f6")]
    [InlineData(BassMixing.Functional, "equal", "3e5f43ed894f60bc")]
    [InlineData(BassMixing.Functional, "cancelling", "e57e8cbe0c39d339")]
    public void Untempered_KeepsTheBaseline(string name, string schedule, string expected)
    {
        (BassConfig config, BassData data) = Fixture(name, schedule);
        string hash = Hash(Bass.Fit(config, data));
        _out.WriteLine($"{name}, {schedule}: {hash}");
        Assert.Equal(expected, hash);
    }

    [Fact]
    public void Tempered_SelfSizing_KeepsTheBaseline()
    {
        (BassConfig config, BassData data) = Fixture(BassMixing.Friedman, "equal", tempered: true);
        BassFit fit = Bass.Fit(config, data);
        string hash = Hash(fit);
        _out.WriteLine($"Friedman tempered: {hash}; levels {string.Join(", ", fit.Diagnostics.Tempering!.Levels)}");
        Assert.Equal("0d067c924599ed93", hash);
    }

    [Fact]
    public void Bars_Tempered_KeepsTheBaseline()
    {
        var x = new double[120];
        var y = new double[120];
        for (int i = 0; i < x.Length; i++)
        {
            x[i] = i / 119.0;
            y[i] = Math.Sin(6 * x[i]) + 0.3 * Math.Cos(17 * x[i]) + (i % 13 == 0 ? 1.5 : 0.0);
        }
        var ensemble = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(5.0), new LocalBetaKernel(50.0),
                                        temperLevels: 4);
        BarsResult result = ensemble.Run(x, y, grid: x, chains: 2, masterSeed: 6, burn: 300, samples: 200);
        var hash = new RunFingerprint();
        hash.Add(result.Fit).Add(result.FitSd).Add(result.MeanKnots).Add(result.SamplesUsed);
        TemperingReadout ladder = result.Tempering!;
        hash.Add(ladder.Barrier).Add(ladder.RoundTripsPerScan);
        foreach (double[] rungs in ladder.Ladders) hash.Add(rungs);
        _out.WriteLine($"BARS tempered: {hash}");
        Assert.Equal("7edb59818f141da9", hash.ToString());
    }

    // SMC fits, hashed on fc39acc before SMC's self-flags were built (the KHAOS trap's decision, 2026-11-08), so the flags,
    // which read the population, are seen to change no draw, weight or normalizer: the draws as Hash reads them, each draw's
    // weight, and each replicate's ln Ẑ, its genealogy error and its stages' b and sweeps. A hinge fit on BassSmcTests' fixture
    // and a KHAOS fit on slice 0's exact-law fixture under the modified g-prior.
    [Theory]
    [InlineData("hinge", "3440c6398ce3c0d5")]
    [InlineData("khaos", "a51939ff36591784")]
    public void Smc_KeepsTheBaseline(string fixture, string expected)
    {
        var smc = new SmcSpec { Particles = 300, Replicates = 2, MutationSteps = 3 };
        (BassConfig config, BassData data) = fixture == "hinge"
            ? (BassSweepTests.Config(maxBasis: 4), BassSweepTests.Data(large: true))
            : (KhaosReachTests.Config("modified g"), KhaosChainTests.Data());
        config = config with { Sampling = config.Sampling with { MasterSeed = 13, Smc = smc } };
        BassFit fit = Bass.Fit(config, data);
        var hash = new RunFingerprint();
        hash.Add(Hash(fit)).Add(fit.DrawWeights!);
        foreach (BassSmcReplicate r in fit.Diagnostics.Smc!.Replicates)
        {
            hash.Add(r.LogNormalizer).Add(r.LogNormalizerError);
            foreach (var stage in r.Stages) hash.Add(stage.Temperature.InverseTemperature).Add(stage.Steps);
        }
        _out.WriteLine($"SMC, {fixture}: {hash}");
        Assert.Equal(expected, hash.ToString());
    }
}
