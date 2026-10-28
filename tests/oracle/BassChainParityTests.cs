using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// BASS's sweep (stage 4's S6) against the stage-4 oracle's exact posterior over an enumerable structure space with data
/// (<c>bass_enumeration_oracle.R</c>, <c>posterior</c>): four chains on D01's schedule with τ's slice step, from the
/// intercept alone, against the exact count law, τ's posterior mean and every term's inclusion probability. Fixture A
/// at M_max = 3 (its mass spread over M = 1–3, with exactly tied states: x and 1 − x span one space with the
/// intercept, so tied terms' inclusions are compared each, the chain moving between them by changes), and fixture C
/// with its categorical predictor at M_max = 2. Batch means across the chains give each comparison's error. Opt-in: it
/// runs where the R toolchain and the BASS checkout are present.
/// </summary>
public sealed class BassChainParityTests
{
    private readonly ITestOutputHelper _out;
    public BassChainParityTests(ITestOutputHelper output) => _out = output;

    private const int Chains = 4, Burn = 3_000, Sweeps = 40_000, Batches = 20;

    [Theory]
    [InlineData("A")]
    [InlineData("C")]
    public void Chain_MatchesTheOraclesPosterior(string fixture)
    {
        if (!BassTargetParityTests.Available) return;
        bool categorical = fixture == "C";
        double[][] rows = categorical ? BassTargetParityTests.FixtureC() : BassTargetParityTests.FixtureA();
        JsonElement r = categorical
            ? BassTargetParityTests.Posterior(rows, "M_max=2", "K_max=1", "categorical=1")
            : BassTargetParityTests.Posterior(rows, "M_max=3");
        BassConfig config = BassTargetParityTests.Config(maxInteraction: categorical ? 1 : 2);
        config = config with { Structure = config.Structure with { MaxBasis = categorical ? 2 : 3 } };
        var model = BassModel.Create(config, BassTargetParityTests.Data(rows, categorical));

        JsonElement p = r.GetProperty("posterior");
        double[] countLaw = [.. p.GetProperty("count").GetProperty("posterior").EnumerateArray().Select(v => v.GetDouble())];
        double[] inclusion = [.. p.GetProperty("inclusion").EnumerateArray().Select(v => v.GetDouble())];
        double tauMean = p.GetProperty("tau").GetProperty("mean").GetDouble();
        BassTerm[] terms = BassTargetParityTests.Terms(r);
        var index = terms.Select((t, i) => (t, i)).ToDictionary(e => e.t, e => e.i);
        Assert.Equal(model.MaxBasis + 1, countLaw.Length);

        // Batch means, chain by chain: of each M's indicator, each term's, and τ.
        var countBatches = new double[countLaw.Length, Chains * Batches];
        var termBatches = new double[terms.Length, Chains * Batches];
        var tauBatches = new double[Chains * Batches];
        int perBatch = Sweeps / Batches;
        for (int c = 0; c < Chains; c++)
        {
            var chain = new BassChain(model, [], Math.Log(1.0 / rows.Length), new Xoshiro256PlusPlus(500 + c), new Xoshiro256PlusPlus(600 + c));
            for (int i = 0; i < Burn; i++) chain.Step();
            for (int i = 0; i < Sweeps; i++)
            {
                BassState s = chain.Step();
                int b = c * Batches + i / perBatch;
                countBatches[s.Count, b] += 1.0 / perBatch;
                foreach (BassTerm t in s.Terms)
                {
                    Assert.True(index.TryGetValue(t, out int j), $"the chain holds {t}, which the oracle's space lacks");
                    termBatches[j, b] += 1.0 / perBatch;
                }
                tauBatches[b] += s.Tau / perBatch;
            }
            _out.WriteLine($"chain {c}: moves {string.Join(", ", chain.Moves.MoveStats().Select(m => $"{m.Key} {m.Accepted}/{m.Attempts}"))}; " +
                           $"refusals {string.Join(", ", chain.Structure.Refusals.Entries.Where(e => e.Count > 0).Select(e => $"{e.Category} {e.Count}"))}");
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
        _out.WriteLine($"{fixture}: {terms.Length} terms' inclusions, worst gap {worst:F2} of its allowance");
    }

    private static double[] Row(double[,] batches, int row) => [.. Enumerable.Range(0, batches.GetLength(1)).Select(b => batches[row, b])];

    private static (double Mean, double Error) Summary(double[] batchMeans)
    {
        double mean = batchMeans.Average();
        double variance = batchMeans.Sum(v => (v - mean) * (v - mean)) / (batchMeans.Length - 1);
        return (mean, Math.Sqrt(variance / batchMeans.Length));
    }
}
