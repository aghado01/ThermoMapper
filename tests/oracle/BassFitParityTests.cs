using System;
using System.Linq;
using System.Text.Json;
using Maths.Regression.Spline.Bass;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// BASS's run (stage 4's S7) end to end against the stage-4 oracle's exact posterior over an enumerable structure space
/// with data (<c>bass_enumeration_oracle.R</c>, <c>posterior</c>): <see cref="Bass.Fit"/> with four chains from exact
/// prior draws on F, against the exact count law, τ's posterior mean and every term's inclusion probability, read from
/// the fit's retained draws through its merged term table. Fixture A at M_max = 3, and fixture C with its categorical
/// predictor at M_max = 2. Batch means across the chains give each comparison's error. Opt-in: it runs where the R
/// toolchain and the BASS checkout are present.
/// </summary>
public sealed class BassFitParityTests
{
    private readonly ITestOutputHelper _out;
    public BassFitParityTests(ITestOutputHelper output) => _out = output;

    private const int Chains = 4, Burn = 3_000, Samples = 40_000, Batches = 20;

    [Theory]
    [InlineData("A")]
    [InlineData("C")]
    public void Fit_MatchesTheOraclesPosterior(string fixture)
    {
        if (!BassTargetParityTests.Available) return;
        bool categorical = fixture == "C";
        double[][] rows = categorical ? BassTargetParityTests.FixtureC() : BassTargetParityTests.FixtureA();
        int maxBasis = categorical ? 2 : 3;
        JsonElement r = categorical
            ? BassTargetParityTests.Posterior(rows, "M_max=2", "K_max=1", "categorical=1")
            : BassTargetParityTests.Posterior(rows, "M_max=3");
        BassConfig config = BassTargetParityTests.Config(maxInteraction: categorical ? 1 : 2);
        config = config with
        {
            Structure = config.Structure with { MaxBasis = maxBasis },
            Sampling = config.Sampling with { Chains = Chains, MasterSeed = 2026, BurnIn = Burn, Samples = Samples },
        };
        BassFit fit = Bass.Fit(config, BassTargetParityTests.Data(rows, categorical));

        JsonElement p = r.GetProperty("posterior");
        double[] countLaw = [.. p.GetProperty("count").GetProperty("posterior").EnumerateArray().Select(v => v.GetDouble())];
        double[] inclusion = [.. p.GetProperty("inclusion").EnumerateArray().Select(v => v.GetDouble())];
        double tauMean = p.GetProperty("tau").GetProperty("mean").GetDouble();
        BassTerm[] terms = BassTargetParityTests.Terms(r);
        var index = terms.Select((t, i) => (t, i)).ToDictionary(e => e.t, e => e.i);
        Assert.Equal(maxBasis + 1, countLaw.Length);
        Assert.Equal(Chains, fit.Draws.Count);

        var countBatches = new double[countLaw.Length, Chains * Batches];
        var termBatches = new double[terms.Length, Chains * Batches];
        var tauBatches = new double[Chains * Batches];
        const int perBatch = Samples / Batches;
        for (int c = 0; c < Chains; c++)
        {
            DrawStore<BassDraw> draws = fit.Draws[c];
            Assert.Equal(Samples, draws.Count);
            for (int i = 0; i < Samples; i++)
            {
                BassDraw d = draws[i];
                int b = c * Batches + i / perBatch;
                countBatches[d.Structure.TermIds.Length, b] += 1.0 / perBatch;
                foreach (int id in d.Structure.TermIds)
                {
                    Assert.True(index.TryGetValue(fit.Terms[id], out int j), $"the fit holds {fit.Terms[id]}, which the oracle's space lacks");
                    termBatches[j, b] += 1.0 / perBatch;
                }
                tauBatches[b] += d.Tau / perBatch;
            }
        }

        for (int m = 0; m < countLaw.Length; m++)
        {
            (double mean, double error) = Summary(Row(countBatches, m));
            _out.WriteLine($"{fixture}, M = {m}: {mean:F4} against {countLaw[m]:F4} (± {error:F4})");
            Assert.True(Math.Abs(mean - countLaw[m]) <= 5 * error + 0.003, $"M = {m}: {mean:F4} against {countLaw[m]:F4} (± {error:F4})");
        }
        (double tau, double tauError) = Summary(tauBatches);
        _out.WriteLine($"{fixture}, E[τ]: {tau:G5} against {tauMean:G5} (± {tauError:G2})");
        Assert.True(Math.Abs(tau - tauMean) <= 5 * tauError, $"E[τ] {tau:G5} against {tauMean:G5} (± {tauError:G2})");
        double worst = 0;
        for (int j = 0; j < terms.Length; j++)
        {
            (double mean, double error) = Summary(Row(termBatches, j));
            double gap = Math.Abs(mean - inclusion[j]);
            worst = Math.Max(worst, gap / (5 * error + 0.005));
            Assert.True(gap <= 5 * error + 0.005, $"{terms[j]}: inclusion {mean:F4} against {inclusion[j]:F4} (± {error:F4})");
        }
        _out.WriteLine($"{fixture}: {terms.Length} terms' inclusions, worst gap {worst:F2} of its allowance; acceptance {fit.Diagnostics.AcceptanceRate:F3}");
    }

    private static double[] Row(double[,] batches, int row) => [.. Enumerable.Range(0, batches.GetLength(1)).Select(b => batches[row, b])];

    private static (double Mean, double Error) Summary(double[] batchMeans)
    {
        double mean = batchMeans.Average();
        double variance = batchMeans.Sum(v => (v - mean) * (v - mean)) / (batchMeans.Length - 1);
        return (mean, Math.Sqrt(variance / batchMeans.Length));
    }
}
