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
/// </summary>
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
    [InlineData(BassMixing.Friedman, "equal", "ee7d60a2738b56d8")]
    [InlineData(BassMixing.Friedman, "cancelling", "414e30ea867f9058")]
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
}
