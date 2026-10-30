using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Stage 8's qualification (C4's step 4; D10's stage 8 question): does an SMC population find the better representation
/// of S10's fixtures, which untempered MCMC chains miss by 12–42 nats (F-M2's Phase 0, M4), and keep it? Exploration: it
/// reads, and sets no default.
/// <list type="bullet">
/// <item><b>Runs.</b> SMC's provisional defaults (<see cref="SmcSpec"/>: 1,000 particles, four replicates, conditional ESS
/// 0.99, five sweeps per stage) on S10's Friedman, mixed and functional fixtures with S10's master seeds, at N = 250,
/// 1,000 and 4,000 and at S10's own size, where M4's chains are the reference. Each fixture's cases come in its stream's
/// order, so a larger N extends S10's data; the functional fixture keeps its 20 locations, so its N counts n_x · 20.</item>
/// <item><b>Read</b> per run: each replicate's stage count, and the run's time; each replicate's weighted mean and median
/// of the log posterior M4 reads (the log marginal plus the log structure prior, τ's density left out), against M4's
/// chains at S10's size; the replicates' spread of the posterior means of M, σ², the log posterior and the fitted
/// values, in posterior-sd units; each replicate's families and Kish ESS; the evidence with its pooled error; σ's 95%
/// interval and whether it covers the truth; and each replicate's top structures' weights.</item>
/// </list>
/// A qualification run (<see cref="BassQualification.Enabled"/>); each is heavy.
/// </summary>
public sealed class BassSmcQualificationTests
{
    private readonly ITestOutputHelper _out;
    public BassSmcQualificationTests(ITestOutputHelper output) => _out = output;

    // M4's chains' median log posteriors at S10's sizes (local-move-measurements-brief.md, M4), best first.
    private static readonly Dictionary<string, double[]> ChainMedians = new()
    {
        [BassMixing.Friedman] = [-892.4, -916.0, -926.2, -933.9],
        [BassMixing.Mixed] = [-309.9, -309.9, -312.5, -331.9],
        [BassMixing.Functional] = [181.1, 168.6, 152.7, 151.7],
    };

    [Theory]
    [InlineData(BassMixing.Friedman, 250)]
    [InlineData(BassMixing.Friedman, 500)]
    [InlineData(BassMixing.Friedman, 1000)]
    [InlineData(BassMixing.Friedman, 4000)]
    [InlineData(BassMixing.Mixed, 250)]
    [InlineData(BassMixing.Mixed, 300)]
    [InlineData(BassMixing.Mixed, 1000)]
    [InlineData(BassMixing.Mixed, 4000)]
    [InlineData(BassMixing.Functional, 260)]
    [InlineData(BassMixing.Functional, 1000)]
    [InlineData(BassMixing.Functional, 1200)]
    [InlineData(BassMixing.Functional, 4000)]
    public void Smc_OnS10sFixtures(string name, int observations)
    {
        if (!BassQualification.Enabled) return;   // heavy
        (BassData data, double sigma, int seed, bool s10) = Fixture(name, observations);
        BassConfig config = BassQualificationRecoveryTests.Config(seed);
        config = config with { Sampling = config.Sampling with { Smc = new SmcSpec() } };
        var clock = Stopwatch.StartNew();
        BassFit fit = Bass.Fit(config, data);
        double seconds = clock.Elapsed.TotalSeconds;
        var model = BassModel.Create(config, data);
        BassSmcReadout smc = fit.Diagnostics.Smc!;
        Assert.True(smc.ReachedPosterior);
        int r = smc.Replicates.Count;
        string Line(IEnumerable<double> values, string format) => string.Join(", ", values.Select(v => v.ToString(format, CultureInfo.InvariantCulture)));

        // Each draw's M4 log posterior, by structure: its fit once per distinct structure, then at the draw's τ.
        var fits = new Dictionary<string, BassStructureFit?>();
        double[] w = (double[])fit.DrawWeights!.Clone();
        var level = new double[w.Length];
        var keys = new string[w.Length];
        int s = 0, unloaded = 0;
        foreach (BassDraw d in fit.Draws.SelectMany(c => c))
        {
            BassTerm[] terms = [.. d.Structure.TermIds.Select(id => fit.Terms[id])];
            string key = BassMixing.Key(terms);
            if (!fits.TryGetValue(key, out BassStructureFit? loaded))
                fits.Add(key, loaded = BassMixing.Loaded(model, terms, Math.Log(d.Tau), out _)?.Fit);
            keys[s] = key;
            if (loaded is not { } structure)
            {
                level[s++] = double.NaN;   // refused when loaded whole, near D06's floor: left out of the readings
                unloaded++;
                continue;
            }
            BassEvaluation e = model.Target.SamplesTau ? model.Target.Evaluate(structure, Math.Log(d.Tau)) : model.Target.Evaluate(structure);
            level[s++] = e.LogMarginal + model.Target.LogStructurePrior(structure);
        }
        if (unloaded > 0) _out.WriteLine($"  {unloaded} draws' structures refused a whole load; their log posteriors are left out");
        for (int i = 0; i < w.Length; i++)
            if (double.IsNaN(level[i])) w[i] = 0.0;
        int n = fit.Draws[0].Count;
        var means = new double[r];
        var medians = new double[r];
        for (int c = 0; c < r; c++)
        {
            var sample = new DrawSample();
            double mean = 0.0;
            for (int i = c * n; i < (c + 1) * n; i++)
            {
                sample.Add(level[i], w[i]);
                mean += w[i] * r * level[i];
            }
            means[c] = mean;
            medians[c] = sample.Quantile(0.5);
        }

        _out.WriteLine($"{name}, N = {data.Response.Observations}: {seconds:F0} s; stages [{string.Join(", ", smc.Replicates.Select(x => x.StageCount))}]; " +
                       $"prior attempts per particle [{Line(smc.Replicates.Select(x => (double)x.PriorAttempts / n), "F2")}]");
        _out.WriteLine($"  log posterior (M4's) per replicate: mean [{Line(means, "F1")}], median [{Line(medians, "F1")}]" +
                       (s10 ? $"; M4's chain medians [{Line(ChainMedians[name], "F1")}]" : ""));
        BassSmcSpread[] fitted = [.. smc.Spread.Where(q => q.Name.StartsWith("fitted", StringComparison.Ordinal))];
        _out.WriteLine($"  spread (sd units): M {smc.Quantity("M").Spread:F3}, σ² {smc.Quantity("sigma2").Spread:F3}, log posterior " +
                       $"{smc.Quantity("log_posterior").Spread:F3}, fitted [{Line(fitted.Select(q => q.Spread), "F3")}]" +
                       (model.Target.SamplesTau ? $", τ {smc.Quantity("tau").Spread:F3}" : ""));
        _out.WriteLine($"  M per replicate [{Line(smc.Quantity("M").ReplicateMeans, "F2")}]; families [{string.Join(", ", smc.Replicates.Select(x => x.Families))}]; " +
                       $"Kish ESS [{Line(smc.Replicates.Select(x => x.EffectiveSampleSize), "F0")}]; acceptance [{Line(smc.Replicates.Select(x => x.Account.AcceptanceRate), "F3")}]");
        BassEvidence evidence = fit.Evidence!;
        _out.WriteLine($"  evidence {evidence.LogEvidence:F2} ± {evidence.MonteCarloError:F2} (within {evidence.SteppingStone.WithinChainError:F2}, between " +
                       $"{evidence.SteppingStone.BetweenChainError:F2}); per replicate [{Line(smc.Replicates.Select(x => x.LogNormalizer), "F2")}] ± " +
                       $"[{Line(smc.Replicates.Select(x => x.LogNormalizerError), "F2")}]");
        var sd = new DrawSample();
        for (int i = 0; i < w.Length; i++) sd.Add(Math.Sqrt(fit.Draws[i / n][i % n].Sigma2), w[i]);
        double lo = sd.Quantile(0.025), hi = sd.Quantile(0.975);
        _out.WriteLine($"  σ's 95% interval [{lo:F4}, {hi:F4}] against {sigma}: {(lo <= sigma && sigma <= hi ? "covered" : "missed")}");

        // Each replicate's top structures, by weight, with their variables and M4's log posterior at the replicate's mean.
        BassMixing.Fixture f = BassMixing.Get(name);
        for (int c = 0; c < r; c++)
        {
            var share = new Dictionary<string, (double Weight, double Level, BassDraw Draw)>();
            for (int i = c * n; i < (c + 1) * n; i++)
            {
                share.TryGetValue(keys[i], out var acc);
                share[keys[i]] = (acc.Weight + w[i] * r, acc.Level + w[i] * r * level[i], fit.Draws[c][i % n]);
            }
            var top = share.OrderByDescending(x => x.Value.Weight).Take(3).Select(x =>
            {
                string variables = string.Join(" + ", x.Value.Draw.Structure.TermIds.Select(id => f.Variables(fit.Terms[id]))
                                                        .GroupBy(v => v).OrderBy(g => g.Key, StringComparer.Ordinal)
                                                        .Select(g => g.Count() > 1 ? $"{g.Count()} × {g.Key}" : g.Key));
                return $"{x.Value.Weight:P1} at {x.Value.Level / x.Value.Weight:F1}: {variables}";
            });
            _out.WriteLine($"  replicate {c}: {share.Count} structures; {string.Join("; ", top)}");
        }
    }

    // The fixture at N observations: its data, its noise sd, S10's master seed, and whether N is S10's own size.
    private static (BassData Data, double Sigma, int Seed, bool S10) Fixture(string name, int observations) => name switch
    {
        BassMixing.Friedman => (BassQualificationRecoveryTests.Friedman(observations).Data, 1.0, 53, observations == 500),
        BassMixing.Mixed => (Mixed(observations), 0.5, 57, observations == 300),
        BassMixing.Functional => (Functional(observations / 20), 0.2, 61, observations == 1200),
        _ => throw new ArgumentOutOfRangeException(nameof(name)),
    };

    // S10's mixed fixture at n cases (BassQualificationRecoveryTests.Mixed is n = 300), the cases in its stream's order.
    private static BassData Mixed(int n)
    {
        var rng = new Xoshiro256PlusPlus(seed: 55);
        var x = new double[3][];
        for (int v = 0; v < 3; v++) x[v] = new double[n];
        var z = new int[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            for (int v = 0; v < 3; v++) x[v][i] = rng.NextDouble();
            z[i] = (int)(3 * rng.NextDouble());
            y[i] = 2 * x[0][i] + (z[i] != 2 ? 3 * x[1][i] : 2 * Math.Sin(2 * Math.PI * x[2][i])) + 0.5 * BassQualification.Gaussian(rng);
        }
        return new BassData(
            [new ContinuousPredictor("x1", x[0]), new ContinuousPredictor("x2", x[1]), new ContinuousPredictor("x3", x[2]),
             new CategoricalPredictor("z", z, ["a", "b", "c"])],
            new ScalarResponse(y));
    }

    // S10's direct functional fixture at n cases on its 20 locations (BassQualificationRecoveryTests.Functional is n = 60).
    private static BassData Functional(int n)
    {
        var rng = new Xoshiro256PlusPlus(seed: 59);
        const int m = 20;
        double[] t = [.. Enumerable.Range(0, m).Select(r => r / (m - 1.0))];
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            for (int r = 0; r < m; r++)
                y[i, r] = 1 + 2 * x0[i] + 1.5 * t[r] * x1[i] + 3 * Math.Max(x0[i] - 0.5, 0) * Math.Sin(Math.PI * t[r]) + 0.2 * BassQualification.Gaussian(rng);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)],
                            new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
    }

    [Fact]
    public void Fixtures_AtS10sSizes_AreS10s()
    {
        // The generators above extend S10's: at its sizes they give its data exactly.
        Assert.Equal(Bass.Fingerprint(new BassConfig(), BassQualificationRecoveryTests.Mixed().Data), Bass.Fingerprint(new BassConfig(), Mixed(300)));
        Assert.Equal(Bass.Fingerprint(new BassConfig(), BassQualificationRecoveryTests.Functional().Data), Bass.Fingerprint(new BassConfig(), Functional(60)));
    }
}
