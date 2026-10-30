using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Friedman's band and modes (the Friedman band brief): four readings on S10's Friedman fixture (n = 500, p = 10, data
/// seed 51, master seed 53), for the user's call among more particles or a lineage-keeping schedule, a compound
/// replacement move, or nothing.
/// <list type="bullet">
/// <item><b>The basins.</b> A structure's basin is its structural key (<see cref="BassMixing.Fixture.Basin"/>), the
/// multiset of its terms' variable sets, signs and knots ignored. Every reading labels the basins M4's untempered chains
/// hold, A, B, … (the best chain's first), so its tables share one legend.</item>
/// <item><b>R1, ℓ at fixed b.</b> M5's dense fixed ladder (its Part 2, built from Part 1's self-sized rungs as M5 builds
/// it): at every rung and sampled scan, ℓ₁ (the untempered log marginal, whatever the rung's b) and M, read after the
/// scan's swaps. Per rung: the histograms pooled and per chain, the moments, Sarle's bimodality coefficient, and two
/// Gaussians against one by BIC (M dequantized); the pairs' Σr (Λ̂) over the band and the ladder; the refusals by kind
/// and cause. Part 1's cold rungs give S10's tempered chains' basins.</item>
/// <item><b>R2, R3, SMC.</b> Fitting SMC at its defaults but the particle count (1,000, 4,000, 16,000) or the mutation
/// steps (20 against 5), four replicates: each replicate's ln Ẑ, families over b, stages, kernel acceptance by kind, log
/// posterior (M4's), σ's interval, and its weighted basin fractions with each basin's log mass ln Ẑ + ln f.</item>
/// <item><b>R4, the refusals.</b> M4's four untempered chains, their refusals by kind and cause over the sampling phase,
/// the numerical refusals split into repeated terms, the angle floor and failed factors; their basins; and, where a cause
/// refuses more than 1% of a kind's proposals, the chains rerun with the floor at 10⁻¹⁴ and 10⁻¹⁸, or the support
/// minimum halved.</item>
/// </list>
/// It reads, it does not judge (<see cref="BassQualification.Enabled"/>). It asserts its own soundness: the harness's
/// chains are the runs' (M, draw for draw), the refusal recorder leaves a ladder's path unchanged, every refusal is read
/// once with its cause, the weights are normalized and the basins' fractions sum to one.
/// </summary>
public sealed class BassFriedmanBandTests
{
    private const double Top = 0.40, Bottom = 0.10, Step = 0.01, BandLow = 0.16, BandHigh = 0.32;
    private const int Bins = 24;
    private readonly ITestOutputHelper _out;
    public BassFriedmanBandTests(ITestOutputHelper output) => _out = output;

    // ---- R1 ----

    [Fact]
    public void Band_LogLikelihoodAtFixedB()
    {
        if (!BassQualification.Enabled) return;   // heavy
        Legend legend = Legend.Of(BassMixing.Get(BassMixing.Friedman));
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman, new TemperingSpec());
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        int chains = fixture.Config.Sampling.Chains;
        var classify = new Classifier(model);

        // Part 1, S10's self-sized ladders (M5's Part 1): the frozen rungs, the pairs' rejections, the cold rung per scan.
        var clock = Stopwatch.StartNew();
        var coarse = new Ladder[chains];
        Parallel.For(0, chains, c =>
        {
            BassLadder ladder = BassMixing.Ladder(model, c, fixture.Config.Sampling.Tempering!);
            BassMixing.Burn(ladder);
            coarse[c] = Sample(fixture, ladder, null);
        });
        double coarseSeconds = clock.Elapsed.TotalSeconds;
        // They are S10's tempered chains: its run's retained draws hold the cold rung's M, draw for draw.
        BassFit tempered = Bass.Fit(fixture.Config, fixture.Data);
        for (int c = 0; c < chains; c++) Assert.Equal(coarse[c].M[0], tempered.Draws[c].Select(d => d.Structure.TermIds.Length));

        // Part 2, M5's dense ladders, each replica's refusals recorded after burn-in; beside them chain 0's twin, unrecorded.
        clock.Restart();
        var dense = new Ladder[chains];
        Ladder? unrecorded = null;
        Parallel.For(0, chains + 1, k =>
        {
            if (k < chains) dense[k] = Dense(fixture, model, k, coarse[k].Rungs, classify);
            else unrecorded = Dense(fixture, model, 0, coarse[0].Rungs, null);
        });
        double denseSeconds = clock.Elapsed.TotalSeconds;
        // The recorders leave the path unchanged: chain 0's dense ladder without them holds the same M and ℓ₁ throughout.
        Ladder twin = unrecorded!;
        for (int i = 0; i < twin.Rungs.Length; i++)
        {
            Assert.Equal(twin.M[i], dense[0].M[i]);
            Assert.Equal(twin.LogLikelihood[i].Select(BassMixing.Bits), dense[0].LogLikelihood[i].Select(BassMixing.Bits));
        }

        _out.WriteLine($"## R1, Friedman tempered: Part 1 self-sized ({coarseSeconds:F1} s), Part 2 dense through [{Bottom}, {Top}] at Δb = {Step} " +
                       $"({denseSeconds:F1} s); {BassMixing.BurnIn:N0} scans of burn-in, {BassMixing.Samples:N0} sampled; ℓ₁ is the untempered log marginal");
        legend.Print(_out);
        Barriers(coarse, dense);
        double[] window = [.. Enumerable.Range(0, (int)Math.Round((Top - Bottom) / Step) + 1).Select(k => Math.Round(Top - k * Step, 2))];
        LogLikelihoods(dense, window);
        Terms(dense, window);
        ByTerms(dense, window);
        Power(dense, window);
        LadderRefusals(dense, window);

        // S10's tempered chains (Part 1's cold rungs) and the dense ladders' cold rungs, by basin.
        Basins("S10's tempered chains (Part 1, b = 1)", legend,
               [.. coarse.Select((l, c) => ($"chain {c}", l.ColdBasins, Uniform(l.ColdBasins.Length)))]);
        Basins("the dense ladders' cold rungs (Part 2, b = 1)", legend,
               [.. dense.Select((l, c) => ($"chain {c}", l.ColdBasins, Uniform(l.ColdBasins.Length)))]);
    }

    // A ladder's sampling phase: its rungs, cold first; per rung and scan, ℓ₁ and M after the scan's swaps; the refusals
    // per rung (null unless recorded); the cold rung's basin per scan; the pairs' rejections and round trips.
    private sealed record Ladder(double[] Rungs, double[][] LogLikelihood, int[][] M, Refusals[]? Refusals, string[] ColdBasins,
                                 double[] Rejection, double RoundTrips, double Ideal);

    // M5's Part 2: the chain's own rungs above the window, the window at Δb = 0.01, then its own rungs below it; burned as
    // M5 burns it, then sampled.
    private static Ladder Dense(BassMixing.Fixture fixture, BassModel model, int chain, double[] frozen, Classifier? classify)
    {
        const double slack = 1e-9;
        int count = (int)Math.Round((Top - Bottom) / Step) + 1;
        double[] rungs = [.. frozen.Where(b => b > Top + slack), .. Enumerable.Range(0, count).Select(k => Math.Round(Top - k * Step, 2)),
                          .. frozen.Where(b => b < Bottom - slack)];
        var spec = new TemperingSpec { Levels = rungs.Length, Adaptive = false, MaxLevels = Math.Max(64, rungs.Length) };
        BassLadder ladder = BassMixing.Ladder(model, chain, spec, rungs);
        BassMixing.Burn(ladder);
        return Sample(fixture, ladder, classify);
    }

    // The sampling phase on a ladder past its burn-in; with `classify`, each replica's target behind a recorder, and each
    // scan's proposals read per rung.
    private static Ladder Sample(BassMixing.Fixture fixture, BassLadder ladder, Classifier? classify)
    {
        IReadOnlyList<ReversibleJumpChain<BassState>> replicas = ladder.Tempering.Replicas;
        int levels = replicas.Count;
        Assert.Equal(levels, ladder.Structures.Count);
        RefusalRecorder[]? recorders = null;
        Refusals[]? refusals = null;
        if (classify is not null)
        {
            recorders = new RefusalRecorder[levels];
            refusals = new Refusals[levels];
            for (int i = 0; i < levels; i++)
            {
                recorders[i] = new RefusalRecorder(ladder.Model.Target);
                replicas[i].ReplaceTarget(recorders[i]);
                refusals[i] = new Refusals();
            }
        }
        var ell = new double[levels][];
        var m = new int[levels][];
        for (int i = 0; i < levels; i++)
        {
            ell[i] = new double[BassMixing.Samples];
            m[i] = new int[BassMixing.Samples];
        }
        var cold = new string[BassMixing.Samples];
        var before = new Counters[levels];
        for (int s = 0; s < BassMixing.Samples; s++)
        {
            if (recorders is not null)
                for (int i = 0; i < levels; i++) before[i] = Counters.Of(ladder.Structures[i]);
            ladder.Step();
            for (int i = 0; i < levels; i++)
            {
                if (recorders is not null) Read(classify!, ladder.Structures[i], before[i], recorders[i], refusals![i]);
                var e = (BassEvaluation)replicas[i].CurrentEvaluation;
                ell[i][s] = e.LogMarginal;
                m[i][s] = replicas[i].Current.Count;
            }
            cold[s] = fixture.Basin(ladder.Current.Terms);
        }
        // ℓ₁ is the untempered log likelihood whatever the rung: the integrated form's at b = 1.
        for (int i = 0; i < levels; i++)
            Assert.Equal(BassMixing.Bits(((BassEvaluation)replicas[i].CurrentEvaluation).LogLikelihood(1.0)), BassMixing.Bits(ell[i][^1]));
        BassTemperingReadout readout = BassTemperingReadout.Of([ladder]);
        double[] rejection = readout.PairRejection[0];
        Assert.All(rejection, r => Assert.InRange(r, 0.0, 1.0 + 1e-12));
        return new Ladder(ladder.Tempering.InverseTemperatures, ell, m, refusals, cold, rejection, readout.RoundTripsPerScan, readout.IdealRoundTripRate);
    }

    // Λ̂ = Σr per chain on both ladders: over the band, the window and the whole ladder.
    private void Barriers(Ladder[] coarse, Ladder[] dense)
    {
        _out.WriteLine("");
        _out.WriteLine($"### Λ̂ = Σr (a pair straddling an edge counted by its share of Δb inside)");
        _out.WriteLine($"chain | Part 1 rungs | Σr band [{BandLow}, {BandHigh}] | Σr window [{Bottom}, {Top}] | Σr ladder | round trips per scan (ideal) | " +
                       $"Part 2 rungs | Σr band | Σr window | Σr ladder | round trips per scan (ideal) | Part 2's largest pair in the window");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < coarse.Length; c++)
        {
            Ladder a = coarse[c], d = dense[c];
            int top = Enumerable.Range(0, d.Rejection.Length).Where(i => d.Rungs[i] <= Top + 1e-9 && d.Rungs[i + 1] >= Bottom - 1e-9).MaxBy(i => d.Rejection[i]);
            _out.WriteLine($"{c} | {a.Rungs.Length} | {Sum(a, BandLow, BandHigh):F2} | {Sum(a, Bottom, Top):F2} | {a.Rejection.Sum():F2} | {a.RoundTrips:G3} ({a.Ideal:G3}) | " +
                           $"{d.Rungs.Length} | {Sum(d, BandLow, BandHigh):F2} | {Sum(d, Bottom, Top):F2} | {d.Rejection.Sum():F2} | {d.RoundTrips:G3} ({d.Ideal:G3}) | " +
                           $"[{d.Rungs[top + 1]:G3}, {d.Rungs[top]:G3}] {d.Rejection[top]:F3}");
        }
    }

    // Σr over the ladder's pairs within [lo, hi], each pair weighted by the share of its Δb inside (M5's).
    private static double Sum(Ladder ladder, double lo, double hi)
    {
        double sum = 0.0;
        for (int i = 0; i < ladder.Rejection.Length; i++)
        {
            double upper = ladder.Rungs[i], lower = ladder.Rungs[i + 1];
            double inside = Math.Min(upper, hi) - Math.Max(lower, lo);
            if (inside > 0.0) sum += inside / (upper - lower) * ladder.Rejection[i];
        }
        return sum;
    }

    // ℓ₁ per rung: the window's rungs pooled over the chains, then every rung of each chain.
    private void LogLikelihoods(Ladder[] dense, double[] window)
    {
        _out.WriteLine("");
        _out.WriteLine($"### ℓ₁ at the window's rungs, pooled over the four chains ({4 * BassMixing.Samples:N0} states; histogram over the rung's pooled range, {Bins} bins)");
        _out.WriteLine("b | mean | sd | skew | ex. kurt | BC | ΔBIC (2 − 1) | minority w | D | means | sds | ESS | range | histogram | chains' means");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        var pooled = new Reading[window.Length];
        Parallel.For(0, window.Length, k =>
        {
            double[][] chains = [.. dense.Select(d => d.LogLikelihood[Array.IndexOf(d.Rungs, window[k])])];
            pooled[k] = Reading.Of(chains);
        });
        for (int k = 0; k < window.Length; k++)
        {
            double[][] chains = [.. dense.Select(d => d.LogLikelihood[Array.IndexOf(d.Rungs, window[k])])];
            Reading r = pooled[k];
            _out.WriteLine($"{window[k]:F2} | {r.Row("F1")} | [{r.Min:F0}, {r.Max:F0}] | `{Histogram(chains.SelectMany(x => x), r.Min, r.Max)}` | " +
                           string.Join(", ", chains.Select(x => x.Average().ToString("F0"))));
        }

        _out.WriteLine("");
        _out.WriteLine("### ℓ₁ at the window's rungs, per chain (each histogram over the rung's pooled range)");
        _out.WriteLine("b | " + string.Join(" | ", Enumerable.Range(0, dense.Length).Select(c => $"chain {c}")));
        _out.WriteLine("---" + string.Concat(Enumerable.Repeat(" | ---", dense.Length)));
        foreach (double b in window)
        {
            double[][] chains = [.. dense.Select(d => d.LogLikelihood[Array.IndexOf(d.Rungs, b)])];
            double lo = chains.Min(x => x.Min()), hi = chains.Max(x => x.Max());
            _out.WriteLine($"{b:F2} | " + string.Join(" | ", chains.Select(x => $"`{Histogram(x, lo, hi)}`")));
        }

        for (int c = 0; c < dense.Length; c++)
        {
            Ladder d = dense[c];
            var readings = new Reading[d.Rungs.Length];
            Parallel.For(0, d.Rungs.Length, i => readings[i] = Reading.Of([d.LogLikelihood[i]]));
            _out.WriteLine("");
            _out.WriteLine($"### ℓ₁ per rung, chain {c} ({d.Rungs.Length} rungs, {BassMixing.Samples:N0} scans; histogram over the chain's own range at the rung)");
            _out.WriteLine("rung | b | mean | sd | skew | ex. kurt | BC | ΔBIC (2 − 1) | minority w | D | means | sds | ESS | switches | range | histogram");
            _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
            for (int i = 0; i < d.Rungs.Length; i++)
            {
                Reading r = readings[i];
                _out.WriteLine($"{i} | {d.Rungs[i]:G4} | {r.Row("F1")} | {r.Switches} | [{r.Min:F0}, {r.Max:F0}] | `{Histogram(d.LogLikelihood[i], r.Min, r.Max)}`");
            }
        }
    }

    // M per rung: its values' counts, and the same readings on M dequantized by a uniform on (−½, ½).
    private void Terms(Ladder[] dense, double[] window)
    {
        double[][][] jittered = Dequantized(dense);
        _out.WriteLine("");
        _out.WriteLine("### M at the window's rungs, pooled over the four chains (the readings on M + U(−½, ½); BC on M itself)");
        _out.WriteLine("b | counts | mean | sd | skew | ex. kurt | BC | ΔBIC (2 − 1) | minority w | D | means | sds | ESS | chains' means");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        var pooled = new Reading[window.Length];
        var raw = new Moments[window.Length];
        Parallel.For(0, window.Length, k =>
        {
            int[] at = [.. dense.Select(d => Array.IndexOf(d.Rungs, window[k]))];
            pooled[k] = Reading.Of([.. at.Select((i, c) => jittered[c][i])]);
            raw[k] = Moments.Of([.. at.SelectMany((i, c) => dense[c].M[i].Select(v => (double)v))]);
        });
        for (int k = 0; k < window.Length; k++)
        {
            int[] at = [.. dense.Select(d => Array.IndexOf(d.Rungs, window[k]))];
            IEnumerable<int> values = at.SelectMany((i, c) => dense[c].M[i]);
            Reading r = pooled[k];
            _out.WriteLine($"{window[k]:F2} | {Counts(values)} | {raw[k].Mean:F2} | {raw[k].Sd:F2} | {raw[k].Skewness:F2} | {raw[k].Kurtosis:F2} | {raw[k].Bimodality:F3} | " +
                           $"{r.Mixture.DeltaBic:F1} | {r.Mixture.Minority:F3} | {r.Mixture.Separation:F2} | {r.Mixture.Mean1:F2}, {r.Mixture.Mean2:F2} | " +
                           $"{r.Mixture.Sd1:F2}, {r.Mixture.Sd2:F2} | {r.Ess:F0} | {string.Join(", ", at.Select((i, c) => dense[c].M[i].Average().ToString("F2")))}");
        }

        for (int c = 0; c < dense.Length; c++)
        {
            Ladder d = dense[c];
            var readings = new Reading[d.Rungs.Length];
            Parallel.For(0, d.Rungs.Length, i => readings[i] = Reading.Of([jittered[c][i]]));
            _out.WriteLine("");
            _out.WriteLine($"### M per rung, chain {c}");
            _out.WriteLine("rung | b | counts | mean | sd | BC | ΔBIC (2 − 1) | minority w | D | means | ESS | switches");
            _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
            for (int i = 0; i < d.Rungs.Length; i++)
            {
                Moments mm = Moments.Of([.. d.M[i].Select(v => (double)v)]);
                Reading r = readings[i];
                _out.WriteLine($"{i} | {d.Rungs[i]:G4} | {Counts(d.M[i])} | {mm.Mean:F2} | {mm.Sd:F2} | {mm.Bimodality:F3} | {r.Mixture.DeltaBic:F1} | " +
                               $"{r.Mixture.Minority:F3} | {r.Mixture.Separation:F2} | {r.Mixture.Mean1:F2}, {r.Mixture.Mean2:F2} | {r.Ess:F0} | {r.Switches}");
            }
        }
    }

    // ℓ₁ given M at the window's rungs, the chains pooled: each M holding 2% or more of the rung's states, its share, its
    // ℓ₁'s mean and sd, its share of the chains' states (chain by chain), and its histogram over the rung's pooled range.
    private void ByTerms(Ladder[] dense, double[] window)
    {
        _out.WriteLine("");
        _out.WriteLine("### ℓ₁ given M at the window's rungs, pooled over the four chains (each M at 2% or more of the rung's states; histogram over the rung's pooled range)");
        _out.WriteLine("b | M | share | per chain | ℓ₁ mean | ℓ₁ sd | ℓ₁ range | histogram");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | ---");
        foreach (double b in window)
        {
            int[] at = [.. dense.Select(d => Array.IndexOf(d.Rungs, b))];
            (int M, double Ell, int Chain)[] states = [.. at.SelectMany((i, c) => dense[c].M[i].Select((m, s) => (m, dense[c].LogLikelihood[i][s], c)))];
            double lo = states.Min(x => x.Ell), hi = states.Max(x => x.Ell);
            foreach (var g in states.GroupBy(x => x.M).OrderBy(g => g.Key))
            {
                if (g.Count() < 0.02 * states.Length) continue;
                double[] ell = [.. g.Select(x => x.Ell)];
                Moments m = Moments.Of(ell);
                string chains = string.Join(", ", Enumerable.Range(0, dense.Length).Select(c => $"{(double)g.Count(x => x.Chain == c) / BassMixing.Samples:P0}"));
                _out.WriteLine($"{b:F2} | {g.Key} | {(double)g.Count() / states.Length:P1} | {chains} | {m.Mean:F1} | {m.Sd:F1} | [{ell.Min():F0}, {ell.Max():F0}] | `{Histogram(ell, lo, hi)}`");
            }
        }
    }

    // Each chain's M at each rung plus a uniform on (−½, ½), from streams derived per chain and rung.
    private static double[][][] Dequantized(Ladder[] dense)
    {
        int[] chainSeeds = SeedTree.Derive(71, dense.Length);
        var jittered = new double[dense.Length][][];
        for (int c = 0; c < dense.Length; c++)
        {
            int[] seeds = SeedTree.Derive(chainSeeds[c], dense[c].Rungs.Length);
            jittered[c] = new double[dense[c].Rungs.Length][];
            for (int i = 0; i < dense[c].Rungs.Length; i++)
            {
                var rng = new Xoshiro256PlusPlus(seeds[i]);
                jittered[c][i] = [.. dense[c].M[i].Select(v => v + rng.NextDouble() - 0.5)];
            }
        }
        return jittered;
    }

    // The tests' sensitivity: iid samples of two Gaussians, the pooled ℓ₁ at the band's edges (b = 0.32 and 0.16), mixed
    // at each minority weight, against one Gaussian; the share of 100 samples each test flags, at a chain's length, the
    // four chains' and 50 (the band's ESS per chain is 3–50).
    private void Power(Ladder[] dense, double[] window)
    {
        double[] High() => [.. dense.SelectMany(d => d.LogLikelihood[Array.IndexOf(d.Rungs, BandHigh)])];
        double[] Low() => [.. dense.SelectMany(d => d.LogLikelihood[Array.IndexOf(d.Rungs, BandLow)])];
        Moments hi = Moments.Of(High()), lo = Moments.Of(Low());
        double[] weights = [0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50];
        int[] sizes = [50, BassMixing.Samples, 4 * BassMixing.Samples];
        const int reps = 100;
        _out.WriteLine("");
        _out.WriteLine($"### The tests' sensitivity: iid mixtures of the b = {BandHigh} phase N({hi.Mean:F1}, {hi.Sd:F1}²) and the b = {BandLow} phase " +
                       $"N({lo.Mean:F1}, {lo.Sd:F1}²), {reps} samples per cell; the share flagged by BC > 5/9, by ΔBIC > 10, and by ΔBIC > 10 with D > 2");
        _out.WriteLine("n | minority | its weight | mixture sd | BC > 5/9 | ΔBIC > 10 | ΔBIC > 10 and D > 2 | median minority w fitted");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | ---");
        int[] seeds = SeedTree.Derive(73, sizes.Length * 2 * weights.Length);
        int cell = 0;
        foreach (int n in sizes)
            foreach (bool minorityLow in new[] { true, false })
                foreach (double w in weights)
                {
                    int seed = seeds[cell++];
                    var flags = new (bool Bc, bool Bic, bool Both, double Minority)[reps];
                    int[] repSeeds = SeedTree.Derive(seed, reps);
                    Parallel.For(0, reps, k =>
                    {
                        var rng = new Xoshiro256PlusPlus(repSeeds[k]);
                        var x = new double[n];
                        for (int i = 0; i < n; i++)
                        {
                            bool low = rng.NextDouble() < (minorityLow ? w : 1.0 - w);
                            x[i] = low ? lo.Mean + lo.Sd * BassQualification.Gaussian(rng) : hi.Mean + hi.Sd * BassQualification.Gaussian(rng);
                        }
                        Moments mm = Moments.Of(x);
                        TwoGaussians g = TwoGaussians.Fit(x);
                        flags[k] = (mm.Bimodality > 5.0 / 9.0, g.DeltaBic > 10, g.DeltaBic > 10 && g.Separation > 2, g.Minority);
                    });
                    double wl = minorityLow ? w : 1.0 - w;
                    double mean = wl * lo.Mean + (1 - wl) * hi.Mean;
                    double sd = Math.Sqrt(wl * (lo.Sd * lo.Sd + (lo.Mean - mean) * (lo.Mean - mean)) + (1 - wl) * (hi.Sd * hi.Sd + (hi.Mean - mean) * (hi.Mean - mean)));
                    _out.WriteLine($"{n:N0} | {(minorityLow ? $"b = {BandLow}" : $"b = {BandHigh}")} | {w:F2} | {sd:F1} | {flags.Count(f => f.Bc) / (double)reps:P0} | " +
                                   $"{flags.Count(f => f.Bic) / (double)reps:P0} | {flags.Count(f => f.Both) / (double)reps:P0} | {BassMixing.Quantile(flags.Select(f => f.Minority), 0.5):F3}");
                }
    }

    // The dense ladders' refusals by rung: the window's rungs pooled over the chains, the cold rung, and the rungs above
    // and below the window pooled.
    private void LadderRefusals(Ladder[] dense, double[] window)
    {
        _out.WriteLine("");
        _out.WriteLine("### Refusals along the dense ladders, as shares of each kind's proposals, the four chains pooled");
        _out.WriteLine(Refusals.Header("b"));
        _out.WriteLine(Refusals.Rule);
        Refusals Pool(Func<Ladder, int, bool> keep)
        {
            var pooled = new Refusals();
            foreach (Ladder d in dense)
                for (int i = 0; i < d.Rungs.Length; i++)
                    if (keep(d, i)) pooled.MergeFrom(d.Refusals![i]);
            return pooled;
        }
        foreach (string row in Pool((d, i) => i == 0).Rows("1 (cold)")) _out.WriteLine(row);
        foreach (string row in Pool((d, i) => i > 0 && d.Rungs[i] > Top + 1e-9).Rows($"({Top}, 1) pooled")) _out.WriteLine(row);
        foreach (double b in window)
            foreach (string row in Pool((d, i) => d.Rungs[i] == b).Rows($"{b:F2}")) _out.WriteLine(row);
        foreach (string row in Pool((d, i) => d.Rungs[i] < Bottom - 1e-9).Rows($"below {Bottom} pooled")) _out.WriteLine(row);
        Refusals all = Pool((_, _) => true);
        foreach (string row in all.Rows("every rung")) _out.WriteLine(row);
        _out.WriteLine("");
        _out.WriteLine(all.Details());
        for (int c = 0; c < dense.Length; c++)
        {
            var chain = new Refusals();
            foreach (Refusals r in dense[c].Refusals!) chain.MergeFrom(r);
            _out.WriteLine($"chain {c}, every rung: " + string.Join("; ", chain.Rows("").Select(r => r.TrimStart(' ', '|'))));
        }
    }

    // ---- R2, R3 ----

    [Theory]
    [InlineData(1000, 5)]
    [InlineData(4000, 5)]
    [InlineData(16000, 5)]
    [InlineData(1000, 20)]
    public void Smc_AgainstParticlesAndMutation(int particles, int mutationSteps)
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman);
        Legend legend = Legend.Of(fixture);
        BassConfig config = fixture.Config with
        {
            Sampling = fixture.Config.Sampling with { Smc = new SmcSpec { Particles = particles, MutationSteps = mutationSteps } },
        };
        var clock = Stopwatch.StartNew();
        BassFit fit = Bass.Fit(config, fixture.Data);
        double seconds = clock.Elapsed.TotalSeconds;
        BassModel model = BassModel.Create(config, fixture.Data);
        BassSmcReadout smc = fit.Diagnostics.Smc!;
        int r = smc.Replicates.Count, n = fit.Draws[0].Count;
        Assert.Equal(particles, n);

        // The weights: each replicate's sum 1/R; read within a replicate as R·w.
        double[] w = fit.DrawWeights!;
        Assert.Equal(r * n, w.Length);
        for (int c = 0; c < r; c++) Assert.True(Math.Abs(w.Skip(c * n).Take(n).Sum() * r - 1.0) <= 1e-12, $"replicate {c}'s weights sum to {w.Skip(c * n).Take(n).Sum() * r:R} of 1");

        // Each draw's basin and log posterior (M4's: the log marginal at its τ plus the log structure prior), its structure
        // loaded whole once; a structure refused whole is left out of the log posterior's readings (C4's rule).
        var basins = new string[w.Length];
        var level = new double[w.Length];
        var fits = new Dictionary<string, BassStructureFit?>();
        int s = 0, unloaded = 0;
        foreach (BassDraw d in fit.Draws.SelectMany(x => x))
        {
            BassTerm[] terms = [.. d.Structure.TermIds.Select(id => fit.Terms[id])];
            basins[s] = fixture.Basin(terms);
            string key = BassMixing.Key(terms);
            if (!fits.TryGetValue(key, out BassStructureFit? loaded))
                fits.Add(key, loaded = BassMixing.Loaded(model, terms, Math.Log(d.Tau), out _)?.Fit);
            if (loaded is { } structure)
            {
                BassEvaluation e = model.Target.Evaluate(structure, Math.Log(d.Tau));
                level[s] = e.LogMarginal + model.Target.LogStructurePrior(structure);
            }
            else
            {
                level[s] = double.NaN;
                unloaded++;
            }
            s++;
        }

        _out.WriteLine($"## R{(mutationSteps == 5 ? 2 : 3)}, Friedman SMC: {particles:N0} particles, {mutationSteps} mutation steps, {r} replicates; {seconds:F0} s; " +
                       $"reached b = 1: {smc.ReachedPosterior}; stages [{string.Join(", ", smc.Replicates.Select(x => x.StageCount))}] of a budget of {config.Sampling.Smc!.MaxStages:N0}");
        if (unloaded > 0) _out.WriteLine($"  {unloaded} draws' structures refused a whole load; their log posteriors are left out");
        legend.Print(_out);

        _out.WriteLine("");
        _out.WriteLine("replicate | seed | stages | ln Ẑ ± genealogy error | final families | Kish ESS | log post mean | log post median | σ's 95% interval | prior attempts per particle");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | --- | ---");
        for (int c = 0; c < r; c++)
        {
            BassSmcReplicate rep = smc.Replicates[c];
            var lp = new DrawSample();
            var sigma = new DrawSample();
            double mean = 0.0, mass = 0.0;
            for (int i = c * n; i < (c + 1) * n; i++)
            {
                sigma.Add(Math.Sqrt(fit.Draws[c][i - c * n].Sigma2), w[i]);
                if (double.IsNaN(level[i])) continue;
                lp.Add(level[i], w[i]);
                mean += w[i] * level[i];
                mass += w[i];
            }
            _out.WriteLine($"{c} | {rep.Seed} | {rep.StageCount} | {rep.LogNormalizer:F2} ± {rep.LogNormalizerError:F2} | {rep.Families} | {rep.EffectiveSampleSize:F0} | " +
                           $"{mean / mass:F1} | {lp.Quantile(0.5):F1} | [{sigma.Quantile(0.025):F4}, {sigma.Quantile(0.975):F4}] | {(double)rep.PriorAttempts / n:F2}");
        }
        double[] logZ = [.. smc.Replicates.Select(x => x.LogNormalizer)];
        double zMean = logZ.Average(), zSd = Math.Sqrt(logZ.Sum(z => (z - zMean) * (z - zMean)) / (r - 1));
        var pooledSigma = new DrawSample();
        for (int i = 0; i < w.Length; i++) pooledSigma.Add(Math.Sqrt(fit.Draws[i / n][i % n].Sigma2), w[i]);
        _out.WriteLine("");
        _out.WriteLine($"ln Ẑ: replicates' mean {zMean:F2}, sd {zSd:F2}, range {logZ.Max() - logZ.Min():F2}" +
                       (fit.Evidence is { } ev
                           ? $"; pooled evidence {ev.LogEvidence:F2} ± {ev.MonteCarloError:F2} (within {ev.SteppingStone.WithinChainError:F2}, between {ev.SteppingStone.BetweenChainError:F2})"
                           : "; no evidence (a replicate short of b = 1)"));
        _out.WriteLine($"σ's 95% interval pooled [{pooledSigma.Quantile(0.025):F4}, {pooledSigma.Quantile(0.975):F4}] against 1; M4's chains' median log posteriors " +
                       $"[{string.Join(", ", legend.Medians.Select(x => x.ToString("F1")))}]");

        // Families over b: at the first stage reaching each mark, and where they first fall to two and to one.
        double[] marks = [0.05, 0.10, BandLow, 0.24, BandHigh, Top, 0.60, 1.0];
        _out.WriteLine("");
        _out.WriteLine("replicate | " + string.Join(" | ", marks.Select(b => $"families at b ≥ {b:G2}")) + " | b at ≤ 2 families | b at 1 family | stages resampled");
        _out.WriteLine("---" + string.Concat(Enumerable.Repeat(" | ---", marks.Length + 3)));
        for (int c = 0; c < r; c++)
        {
            IReadOnlyList<SmcStageReport> stages = smc.Replicates[c].Stages;
            string At(double b)
            {
                foreach (SmcStageReport st in stages)
                    if (st.Temperature.InverseTemperature >= b - 1e-12) return $"{st.Families}";
                return "—";
            }
            string First(int families)
            {
                foreach (SmcStageReport st in stages)
                    if (st.Families <= families) return $"{st.Temperature.InverseTemperature:F3}";
                return "—";
            }
            _out.WriteLine($"{c} | {string.Join(" | ", marks.Select(At))} | {First(2)} | {First(1)} | {stages.Count(st => st.Resampled)}");
        }

        // The kernel's acceptance by kind.
        _out.WriteLine("");
        _out.WriteLine("replicate | sweeps | " + string.Join(" | ", smc.Replicates[0].Account.Kinds.Entries.Where(k => k.Attempts > 0).Select(k => $"{k.Kind}: proposals, accepted, refused")));
        _out.WriteLine("---" + string.Concat(Enumerable.Repeat(" | ---", 1 + smc.Replicates[0].Account.Kinds.Entries.Count(k => k.Attempts > 0))));
        for (int c = 0; c < r; c++)
        {
            BassSmcAccount account = smc.Replicates[c].Account;
            _out.WriteLine($"{c} | {account.Sweeps:N0} | " + string.Join(" | ", account.Kinds.Entries.Where(k => smc.Replicates[0].Account.Kinds[k.Kind].Attempts > 0)
                .Select(k => $"{k.Attempts:N0}, {k.AcceptanceRate:P2}, {(k.Attempts > 0 ? (double)k.Refused / k.Attempts : double.NaN):P2}")));
        }

        // The basins: each replicate's weighted fractions and each basin's log mass ln Ẑ + ln f.
        var sets = new (string Name, string[] Basins, double[] Weights)[r];
        for (int c = 0; c < r; c++) sets[c] = ($"replicate {c}", basins[(c * n)..((c + 1) * n)], [.. w[(c * n)..((c + 1) * n)].Select(x => x * r)]);
        Basins($"the replicates' final populations (b = 1), weighted", legend, sets, logZ, level);
    }

    // ---- R4 ----

    [Fact]
    public void Refusals_OnM4sChains()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassMixing.Fixture fixture = BassMixing.Get(BassMixing.Friedman);
        var clock = Stopwatch.StartNew();
        BassModel model = BassModel.Create(fixture.Config, fixture.Data);
        Run run = Run.Of(fixture, model);
        double seconds = clock.Elapsed.TotalSeconds;
        // The harness's chains are M4's, S10's: the run's retained draws hold the same M, draw for draw.
        BassFit fit = Bass.Fit(fixture.Config, fixture.Data);
        for (int c = 0; c < run.States.Length; c++) Assert.Equal(run.States[c].Select(h => h.M), fit.Draws[c].Select(d => d.Structure.TermIds.Length));
        var legend = new Legend(run);

        // The structural key ignores signs and knots, and reads the variables: every hinge of a retained structure flipped
        // and moved keeps its basin; one hinge moved to another variable leaves it.
        BassTerm[] held = [.. fit.Draws[1][0].Structure.TermIds.Select(id => fit.Terms[id])];
        BassTerm Moved(BassTerm t) => new([.. t.Continuous.Select(h => new HingeFactor(h.Variable, -h.Sign, (h.Knot + 7) % BassMixing.Knots(model, BassPart.Continuous, h.Variable)))],
                                          t.Categorical, t.Functional);
        Assert.Equal(fixture.Basin(held), fixture.Basin(held.Select(Moved)));
        HingeFactor first = held[0].Continuous[0];
        int other = Enumerable.Range(0, fixture.Continuous.Length).First(v => held[0].Continuous.All(h => h.Variable != v));
        BassTerm[] elsewhere = [BassMixing.Replaced(held[0], new BassMixing.FactorAt(BassPart.Continuous, 0), new HingeFactor(other, first.Sign, 0), null), .. held[1..]];
        Assert.NotEqual(fixture.Basin(held), fixture.Basin(elsewhere));

        _out.WriteLine($"## R4, Friedman untempered: M4's four chains, {BassMixing.Samples:N0} sampled sweeps each ({seconds:F1} s); the support minimum " +
                       $"{model.Resolved.SupportMinimum} nonzero points, the angle floor {fixture.Config.Structure.AngleFloor:G2}");
        legend.Print(_out);
        _out.WriteLine("");
        _out.WriteLine("### Refusals over the sampling phase, as shares of each kind's proposals");
        _out.WriteLine(Refusals.Header("chain"));
        _out.WriteLine(Refusals.Rule);
        var pooled = new Refusals();
        for (int c = 0; c < run.Refusals.Length; c++)
        {
            foreach (string row in run.Refusals[c].Rows($"{c}")) _out.WriteLine(row);
            pooled.MergeFrom(run.Refusals[c]);
        }
        foreach (string row in pooled.Rows("pooled")) _out.WriteLine(row);
        _out.WriteLine("");
        _out.WriteLine(pooled.Details());
        for (int c = 0; c < run.Refusals.Length; c++) _out.WriteLine($"chain {c}: {run.Refusals[c].Details()}");

        Basins("M4's chains (b = 1)", legend, [.. run.States.Select((h, c) => ($"chain {c}", h.Select(x => x.Basin).ToArray(), Uniform(h.Length)))],
               level: null, levels: [.. run.States.Select(h => h.Select(x => x.LogPosterior).ToArray())]);

        // The reruns the brief's rule asks for: a cause refusing more than 1% of a kind's proposals, pooled over the chains.
        double angle = pooled.WorstShare(Cause.AngleFloor), support = pooled.WorstShare(Cause.Support);
        double repeatAndAngle = pooled.WorstShare(Cause.AngleFloor, Cause.Repeat, Cause.Factor);
        _out.WriteLine("");
        _out.WriteLine($"The rule: the angle floor refuses at most {angle:P3} of a kind's proposals ({repeatAndAngle:P3} with the repeated terms and failed factors, " +
                       $"which no floor admits); the support rule at most {support:P3}.");
        var reruns = new List<(string Label, BassConfig Config)>();
        BassStructureSpec structure = fixture.Config.Structure;
        if (angle > 0.01)
        {
            reruns.Add(("angle floor 1e-14", fixture.Config with { Structure = structure with { AngleFloor = 1e-14 } }));
            reruns.Add(("angle floor 1e-18", fixture.Config with { Structure = structure with { AngleFloor = 1e-18 } }));
        }
        if (support > 0.01)
            reruns.Add(($"support minimum {model.Resolved.SupportMinimum / 2}",
                        fixture.Config with { Structure = structure with { Support = structure.Support with { Minimum = model.Resolved.SupportMinimum / 2 } } }));
        if (reruns.Count == 0)
        {
            _out.WriteLine("Neither cause refuses more than 1%: no rerun.");
            return;
        }

        _out.WriteLine("");
        _out.WriteLine("### The reruns: M4's chains, the same streams, under a changed support F (each a different posterior: they measure the rule's share of the stall)");
        _out.WriteLine("setting | chain | median log post (5–95%) | below the best chain's median | M mean (range) | top basin, its share | the run's R̂ of M | R̂ of the log posterior | s");
        _out.WriteLine("--- | --- | --- | --- | --- | --- | --- | --- | ---");
        Row(" S10's", run, fit, legend, seconds);
        var runs = new List<(string Label, Run Run)>();
        foreach ((string label, BassConfig config) in reruns)
        {
            clock.Restart();
            BassModel changed = BassModel.Create(config, fixture.Data);
            Run rerun = Run.Of(fixture, changed);
            double s = clock.Elapsed.TotalSeconds;
            BassFit refit = Bass.Fit(config, fixture.Data);
            for (int c = 0; c < rerun.States.Length; c++) Assert.Equal(rerun.States[c].Select(h => h.M), refit.Draws[c].Select(d => d.Structure.TermIds.Length));
            Row(label, rerun, refit, legend, s);
            runs.Add((label, rerun));
        }
        foreach ((string label, Run rerun) in runs)
        {
            _out.WriteLine("");
            _out.WriteLine($"{label}: refusals pooled over the chains");
            _out.WriteLine(Refusals.Header("setting"));
            _out.WriteLine(Refusals.Rule);
            var p = new Refusals();
            foreach (Refusals x in rerun.Refusals) p.MergeFrom(x);
            foreach (string row in p.Rows(label)) _out.WriteLine(row);
            Basins($"the chains under {label}", legend, [.. rerun.States.Select((h, c) => ($"chain {c}", h.Select(x => x.Basin).ToArray(), Uniform(h.Length)))],
                   level: null, levels: [.. rerun.States.Select(h => h.Select(x => x.LogPosterior).ToArray())]);
        }
    }

    // A run's row per chain: its log posterior's median and band, M, its most held basin, and the run's R̂s.
    private void Row(string label, Run run, BassFit fit, Legend legend, double seconds)
    {
        double[] medians = [.. run.States.Select(h => BassMixing.Quantile(h.Select(x => x.LogPosterior), 0.5))];
        for (int c = 0; c < run.States.Length; c++)
        {
            Sampled[] h = run.States[c];
            var top = h.GroupBy(x => x.Basin).MaxBy(g => g.Count())!;
            _out.WriteLine($"{label} | {c} | {medians[c]:F1} ({BassMixing.Quantile(h.Select(x => x.LogPosterior), 0.05):F1}, {BassMixing.Quantile(h.Select(x => x.LogPosterior), 0.95):F1}) | " +
                           $"{medians.Max() - medians[c]:F1} | {h.Average(x => x.M):F2} ({h.Min(x => x.M)}–{h.Max(x => x.M)}) | {legend.Name(top.Key)}, {(double)top.Count() / h.Length:P1} | " +
                           $"{fit.Diagnostics.Count.RHat.Value:F2} | {fit.Diagnostics.LogPosterior.RHat.Value:F2} | {seconds:F1}");
        }
    }

    // ---- The basins ----

    // One sampled state: its M, log posterior (M4's) and basin.
    private sealed record Sampled(int M, double LogPosterior, string Basin);

    // Untempered chains of `model` started, seeded and burned as the run does it, each with a refusal recorder in its
    // target's place: per sampled state its M, log posterior and basin, and the refusals over the sampling phase.
    private sealed record Run(Sampled[][] States, Refusals[] Refusals)
    {
        public static Run Of(BassMixing.Fixture fixture, BassModel model)
        {
            int chains = model.Config.Sampling.Chains;
            var states = new Sampled[chains][];
            var refusals = new Refusals[chains];
            var classify = new Classifier(model);
            Parallel.For(0, chains, c =>
            {
                (BassState start, int[] seeds) = BassMixing.Start(model, c);
                var recorder = new RefusalRecorder(model.Target);
                var chain = new BassChain(model, start.Terms, start.LogTau, new Xoshiro256PlusPlus(seeds[1]), new Xoshiro256PlusPlus(seeds[2]), recorder);
                for (int s = 0; s < BassMixing.BurnIn; s++)
                {
                    chain.Step();
                    recorder.Take();
                }
                chain.Structure.ResetTallies();
                var tally = new Refusals();
                states[c] = new Sampled[BassMixing.Samples];
                for (int s = 0; s < BassMixing.Samples; s++)
                {
                    Counters before = Counters.Of(chain.Structure);
                    BassState state = chain.Step();
                    Read(classify, chain.Structure, before, recorder, tally);
                    var evaluation = (BassEvaluation)chain.Moves.CurrentEvaluation;
                    states[c][s] = new Sampled(state.Count, evaluation.LogMarginal + model.Target.LogStructurePrior(state.Fit), fixture.Basin(state.Terms));
                }
                refusals[c] = tally;
            });
            return new Run(states, refusals);
        }
    }

    // M4's basins: the keys its untempered chains hold in at least 1% of some chain's sampling phase, labelled A, B, …,
    // the best chain's (highest median log posterior) by share first, then each other chain's in chain order.
    private sealed class Legend
    {
        private readonly Dictionary<string, string> _labels = [];
        private readonly Sampled[][] _states;

        public Legend(Run run)
        {
            _states = run.States;
            Medians = [.. _states.Select(h => BassMixing.Quantile(h.Select(x => x.LogPosterior), 0.5))];
            Best = Array.IndexOf(Medians, Medians.Max());
            IEnumerable<int> order = [Best, .. Enumerable.Range(0, _states.Length).Where(c => c != Best)];
            foreach (int c in order)
                foreach (var g in _states[c].GroupBy(x => x.Basin).OrderByDescending(g => g.Count()).ThenBy(g => g.Key, StringComparer.Ordinal))
                    if (g.Count() >= 0.01 * _states[c].Length && !_labels.ContainsKey(g.Key))
                        _labels.Add(g.Key, Letter(_labels.Count));
        }

        public static Legend Of(BassMixing.Fixture fixture) => new(Run.Of(fixture, BassModel.Create(fixture.Config, fixture.Data)));

        /// <summary>Each chain's median log posterior.</summary>
        public double[] Medians { get; }

        public int Best { get; }

        public string Name(string basin) => _labels.TryGetValue(basin, out string? label) ? label : basin;

        /// <summary>The labelled basins, in label order.</summary>
        public IEnumerable<string> Keys => _labels.Keys;

        public string Order(string basin) => _labels.TryGetValue(basin, out string? label) ? $"0{label.Length:D2}{label}" : $"1{basin}";

        private static string Letter(int k) => k < 26 ? $"{(char)('A' + k)}" : $"{(char)('A' + k / 26 - 1)}{(char)('A' + k % 26)}";

        public void Print(ITestOutputHelper output)
        {
            output.WriteLine("");
            output.WriteLine($"### The basins M4's chains hold (chain {Best} the best, median log posterior {Medians[Best]:F1}; the others' " +
                             $"{string.Join(", ", Enumerable.Range(0, Medians.Length).Where(c => c != Best).Select(c => $"chain {c} {Medians[c]:F1}"))})");
            output.WriteLine($"label | basin | M | {string.Join(" | ", Enumerable.Range(0, _states.Length).Select(c => $"share in chain {c}"))} | median log post where held");
            output.WriteLine("--- | --- | ---" + string.Concat(Enumerable.Repeat(" | ---", _states.Length + 1)));
            foreach ((string basin, string label) in _labels)
            {
                Sampled[] held = [.. _states.SelectMany(h => h).Where(x => x.Basin == basin)];
                output.WriteLine($"{label} | {basin} | {held[0].M} | " +
                                 string.Join(" | ", _states.Select(h => $"{(double)h.Count(x => x.Basin == basin) / h.Length:P1}")) +
                                 $" | {BassMixing.Quantile(held.Select(x => x.LogPosterior), 0.5):F1}");
            }
        }
    }

    // A basin's M, its terms counted from its key.
    private static int TermCount(string basin)
        => basin == "∅" ? 0 : basin.Split(" + ").Sum(t => t.Contains(" × ") ? int.Parse(t[..t.IndexOf(' ')]) : 1);

    // Equal weights, each 1/n.
    private static double[] Uniform(int n) => [.. Enumerable.Repeat(1.0 / n, n)];

    // The basins' weighted fractions in each named sample, every basin at 0.1% or more in some sample, the rest as
    // "other"; with each sample's ln Ẑ, each basin's log mass ln Ẑ + ln f; with per-state log posteriors (pooled, or per
    // sample), each basin's median among them.
    private void Basins(string title, Legend legend, (string Name, string[] Basins, double[] Weights)[] sets, double[]? logZ = null,
                        double[]? level = null, double[][]? levels = null)
    {
        var fractions = new Dictionary<string, double>[sets.Length];
        for (int k = 0; k < sets.Length; k++)
        {
            (_, string[] basins, double[] weights) = sets[k];
            Assert.Equal(basins.Length, weights.Length);
            Assert.True(Math.Abs(weights.Sum() - 1.0) <= 1e-9, $"{sets[k].Name}'s weights sum to {weights.Sum():R}");
            fractions[k] = [];
            for (int i = 0; i < basins.Length; i++) fractions[k][basins[i]] = fractions[k].GetValueOrDefault(basins[i]) + weights[i];
            Assert.True(Math.Abs(fractions[k].Values.Sum() - 1.0) <= 1e-9);
        }
        string[] shown = [.. fractions.SelectMany(f => f).Where(x => x.Value >= 0.001).Select(x => x.Key).Distinct().OrderBy(legend.Order, StringComparer.Ordinal)];
        _out.WriteLine("");
        _out.WriteLine($"### Basins: {title}{(logZ is null ? "" : "; ln Ẑ_basin = ln Ẑ + ln f")}");
        string massHeader = logZ is null ? "" : " | " + string.Join(" | ", sets.Select(x => $"ln Ẑ_basin, {x.Name}"));
        string levelHeader = level is not null || levels is not null ? " | median log post where held" : "";
        _out.WriteLine($"basin | M | {string.Join(" | ", sets.Select(x => x.Name))} | pooled{massHeader}{levelHeader}");
        int columns = 3 + sets.Length + (logZ is null ? 0 : sets.Length) + (levelHeader.Length > 0 ? 1 : 0);
        _out.WriteLine("---" + string.Concat(Enumerable.Repeat(" | ---", columns - 1)));
        foreach (string basin in shown)
        {
            double pooled = fractions.Average(f => f.GetValueOrDefault(basin));
            string masses = logZ is null ? "" : " | " + string.Join(" | ", fractions.Select((f, k) => f.GetValueOrDefault(basin) > 0 ? $"{logZ[k] + Math.Log(f[basin]):F2}" : "—"));
            string median = "";
            if (level is not null)
            {
                IEnumerable<double> held = sets.SelectMany((x, k) => x.Basins.Select((b, i) => (b, i: k * x.Basins.Length + i))).Where(t => t.b == basin)
                                               .Select(t => level[t.i]).Where(double.IsFinite);
                median = $" | {BassMixing.Quantile(held, 0.5):F1}";
            }
            else if (levels is not null)
                median = " | " + string.Join(", ", sets.Select((x, k) => x.Basins.Select((b, i) => (b, i)).Where(t => t.b == basin).Select(t => levels[k][t.i]).ToArray())
                                                       .Select(v => v.Length == 0 ? "—" : BassMixing.Quantile(v, 0.5).ToString("F1")));
            _out.WriteLine($"{legend.Name(basin)} | {TermCount(basin)} | " +
                           $"{string.Join(" | ", fractions.Select(f => $"{f.GetValueOrDefault(basin):P1}"))} | {pooled:P1}{masses}{median}");
        }
        double[] other = [.. fractions.Select(f => f.Where(x => !shown.Contains(x.Key)).Sum(x => x.Value))];
        _out.WriteLine($"other | — | {string.Join(" | ", other.Select(x => $"{x:P1}"))} | {other.Average():P1}" +
                       (logZ is null ? "" : " | " + string.Join(" | ", other.Select((x, k) => x > 0 ? $"{logZ[k] + Math.Log(x):F2}" : "—"))) +
                       (levelHeader.Length > 0 ? " | —" : ""));
        if (logZ is not null)
            for (int k = 0; k < sets.Length; k++)
            {
                // The basins' masses add up to the replicate's: ln Σ exp(ln Ẑ + ln f) = ln Ẑ.
                double total = LogSumExp(fractions[k].Values.Where(f => f > 0).Select(f => logZ[k] + Math.Log(f)));
                Assert.True(Math.Abs(total - logZ[k]) <= 1e-9 * Math.Max(1.0, Math.Abs(logZ[k])), $"{sets[k].Name}'s basins' log masses add to {total:R}, its ln Ẑ is {logZ[k]:R}");
            }
        string[] absent = [.. legend.Keys.Where(b => !shown.Contains(b)).Select(legend.Name)];
        if (absent.Length > 0) _out.WriteLine($"M4's basins under 0.1% in every sample here: {string.Join(", ", absent)}");

        // A coarser view beside the key: the share by the count of x1·x2 terms, the count that sets M4's modes apart.
        _out.WriteLine("");
        _out.WriteLine("x1·x2 terms | " + string.Join(" | ", sets.Select(x => x.Name)));
        _out.WriteLine("---" + string.Concat(Enumerable.Repeat(" | ---", sets.Length)));
        int[] counts = [.. fractions.SelectMany(f => f.Keys).Select(Pairs).Distinct().Order()];
        foreach (int k in counts)
            _out.WriteLine($"{k} | " + string.Join(" | ", fractions.Select(f => $"{f.Where(x => Pairs(x.Key) == k).Sum(x => x.Value):P1}")));
    }

    // The basin's count of x1·x2 terms.
    private static int Pairs(string basin)
        => basin.Split(" + ").Sum(t => t == "x1·x2" ? 1 : t.EndsWith(" × x1·x2", StringComparison.Ordinal) ? int.Parse(t[..t.IndexOf(' ')]) : 0);

    private static double LogSumExp(IEnumerable<double> values)
    {
        double[] v = [.. values];
        double max = v.Max();
        return max + Math.Log(v.Sum(x => Math.Exp(x - max)));
    }

    // ---- The refusals ----

    private static readonly BassMoveKind[] KindsOf = Enum.GetValues<BassMoveKind>();
    private static readonly BassRefusal[] ReasonsOf = Enum.GetValues<BassRefusal>();
    private static readonly Cause[] CausesOf = Enum.GetValues<Cause>();

    // A refusal's cause: the design's or the target's reason, the numerical one split into a repeated term, the angle
    // floor and a failed factor; a prune that would leave the intercept.
    private enum Cause { Cap, Support, Repeat, AngleFloor, Factor, VarianceInflation, ExactFit, Intercept }

    // A moves' tallies as counts: attempts and refusals by kind, refusals by reason.
    private readonly record struct Counters(long[] Attempts, long[] Refused, long[] Reasons)
    {
        public static Counters Of(BassStructureMoves moves) => new(
            [.. KindsOf.Select(k => moves.Kinds[k].Attempts)], [.. KindsOf.Select(k => moves.Kinds[k].Refused)], [.. ReasonsOf.Select(moves.Refusals.Count)]);
    }

    // A target that passes every evaluation through and keeps the candidate outside F it was last shown: a sweep proposes
    // one structure, and τ's slice evaluates only states in F, so a sweep shows it at most one.
    private sealed class RefusalRecorder(BassTarget target) : ISamplingTarget<BassState>
    {
        private BassState? _refused;
        private int _count;

        public ITargetEvaluation Evaluate(BassState state)
        {
            BassEvaluation evaluation = target.Evaluate(state);
            if (evaluation.Refusal is not null)
            {
                _refused = state;
                _count++;
            }
            return evaluation;
        }

        // The candidate refused since the last call, or null.
        public BassState? Take()
        {
            Assert.InRange(_count, 0, 1);
            BassState? refused = _refused;
            _refused = null;
            _count = 0;
            return refused;
        }
    }

    // One sweep's proposal, read from the change in its moves' tallies and the recorder: its kind, and its refusal's cause.
    private static void Read(Classifier classify, BassStructureMoves moves, Counters before, RefusalRecorder recorder, Refusals into)
    {
        Counters after = Counters.Of(moves);
        int kind = Changed(before.Attempts, after.Attempts), reason = Changed(before.Reasons, after.Reasons);
        BassState? refused = recorder.Take();
        if (kind < 0)
        {
            Assert.Null(refused);
            Assert.Equal(-1, reason);
            return;
        }
        into.Attempts[kind]++;
        if (after.Refused[kind] == before.Refused[kind])
        {
            Assert.Null(refused);
            Assert.Equal(-1, reason);
            return;
        }
        Assert.Equal(before.Refused[kind] + 1, after.Refused[kind]);
        if (reason < 0)
        {
            Assert.Null(refused);   // a prune of a term's only factor: no candidate
            into.Refused[kind, (int)Cause.Intercept]++;
            return;
        }
        Assert.NotNull(refused);
        Assert.Equal((BassRefusal?)ReasonsOf[reason], classify.Model.Target.Refusal(refused.Fit));
        into.Refused[kind, (int)classify.Classify(refused, ReasonsOf[reason], into)]++;
    }

    // The one index whose count rose by one, or −1 when none did.
    private static int Changed(long[] before, long[] after)
    {
        int changed = -1;
        for (int k = 0; k < before.Length; k++)
        {
            long d = after[k] - before[k];
            Assert.InRange(d, 0, 1);
            if (d == 0) continue;
            Assert.Equal(-1, changed);
            changed = k;
        }
        return changed;
    }

    // Each refusal's cause. A numerical refusal is a repeated term when the new term (the candidate's last) repeats one
    // held; otherwise the candidate is loaded whole at an angle floor of 1e-300, and a failed factor there is a failed
    // factor, an admitted candidate the angle floor with its smallest leave-one-out angle. A support refusal records the
    // new term's nonzero count.
    private sealed class Classifier(BassModel model)
    {
        private readonly BassModel _permissive = BassModel.Create(model.Config with { Structure = model.Config.Structure with { AngleFloor = 1e-300 } }, model.Data);

        public BassModel Model => model;

        public Cause Classify(BassState candidate, BassRefusal reason, Refusals into)
        {
            BassTerm added = candidate.Terms[^1];
            switch (reason)
            {
                case BassRefusal.Cap: return Cause.Cap;
                case BassRefusal.VarianceInflation: return Cause.VarianceInflation;
                case BassRefusal.ExactFit: return Cause.ExactFit;
                case BassRefusal.Support:
                    var column = new double[model.Basis.Inputs.Cases.Cases];
                    model.Basis.CaseColumn(added, model.Basis.Inputs.Cases, column);
                    int nonzero = column.Count(v => v != 0.0);
                    Assert.True(nonzero < model.Resolved.SupportMinimum, $"A support refusal's term has {nonzero} nonzero points.");
                    into.Supports.Add(nonzero);
                    return Cause.Support;
            }
            Assert.Equal(BassRefusal.Numerical, reason);
            for (int j = 0; j < candidate.Count - 1; j++)
                if (candidate.Terms[j].Equals(added)) return Cause.Repeat;
            BassDesign whole = _permissive.NewDesign();
            if (whole.Load(candidate.Terms) is { } refusal)
            {
                Assert.Equal(BassRefusal.Numerical, refusal);
                return Cause.Factor;
            }
            double angle = whole.MinAngle;
            if (angle >= model.Config.Structure.AngleFloor) into.AdmittedWhole++;
            into.Angles.Add(angle);
            return Cause.AngleFloor;
        }
    }

    // Proposals and refusals by kind and cause, with the angle refusals' smallest angles and the support refusals' nonzero
    // counts.
    private sealed class Refusals
    {
        public long[] Attempts { get; } = new long[KindsOf.Length];
        public long[,] Refused { get; } = new long[KindsOf.Length, CausesOf.Length];
        public List<double> Angles { get; } = [];
        public List<int> Supports { get; } = [];
        public int AdmittedWhole { get; set; }

        public void MergeFrom(Refusals other)
        {
            for (int k = 0; k < KindsOf.Length; k++)
            {
                Attempts[k] += other.Attempts[k];
                for (int c = 0; c < CausesOf.Length; c++) Refused[k, c] += other.Refused[k, c];
            }
            Angles.AddRange(other.Angles);
            Supports.AddRange(other.Supports);
            AdmittedWhole += other.AdmittedWhole;
        }

        public static string Header(string first) => $"{first} | kind | proposals | refused | " + string.Join(" | ", CausesOf.Select(c => $"{c}"));

        public static string Rule => "---" + string.Concat(Enumerable.Repeat(" | ---", 3 + CausesOf.Length));

        public IEnumerable<string> Rows(string label)
        {
            for (int k = 0; k < KindsOf.Length; k++)
            {
                if (Attempts[k] == 0) continue;
                long refused = Enumerable.Range(0, CausesOf.Length).Sum(c => Refused[k, c]);
                yield return $"{label} | {KindsOf[k]} | {Attempts[k]:N0} | {Share(refused, Attempts[k])} | " +
                             string.Join(" | ", Enumerable.Range(0, CausesOf.Length).Select(c => Share(Refused[k, c], Attempts[k])));
            }
        }

        // The largest share of a kind's proposals the causes refuse together.
        public double WorstShare(params Cause[] causes)
        {
            double worst = 0.0;
            for (int k = 0; k < KindsOf.Length; k++)
                if (Attempts[k] > 0) worst = Math.Max(worst, causes.Sum(c => Refused[k, (int)c]) / (double)Attempts[k]);
            return worst;
        }

        public string Details()
        {
            string angles = Angles.Count == 0
                ? "no angle-floor refusal"
                : $"{Angles.Count} angle-floor refusals, their smallest angles: {Angles.Count(a => a < 1e-18)} below 1e-18, {Angles.Count(a => a >= 1e-18 && a < 1e-14)} in [1e-18, 1e-14), " +
                  $"{Angles.Count(a => a >= 1e-14 && a < 1e-10)} in [1e-14, 1e-10), {AdmittedWhole} at or above the floor when loaded whole (min {Angles.Min():G3}, max {Angles.Max():G3})";
            string supports = Supports.Count == 0
                ? "no support refusal"
                : $"{Supports.Count} support refusals, the new term's nonzero points: {Supports.Count(x => x == 0)} at 0, {Supports.Count(x => x > 0 && x < 10)} in [1, 10), " +
                  $"{Supports.Count(x => x >= 10)} in [10, 20) (admitted at half the minimum)";
            return $"{angles}; {supports}";
        }

        private static string Share(long count, long of) => count == 0 ? "0" : $"{100.0 * count / of:F2}%";
    }

    // ---- The readings' statistics ----

    // One sample's readings, pooled over its sequences: the moments, two Gaussians against one, the ESS of the mean (the
    // sequences as chains), and the switches between the two components along each sequence.
    private sealed record Reading(Moments Moments, TwoGaussians Mixture, double Ess, int Switches, double Min, double Max)
    {
        public static Reading Of(double[][] sequences)
        {
            double[] all = [.. sequences.SelectMany(x => x)];
            Moments m = Moments.Of(all);
            TwoGaussians g = TwoGaussians.Fit(all);
            int switches = 0;
            if (double.IsFinite(g.DeltaBic))
                foreach (double[] x in sequences)
                    for (int i = 1; i < x.Length; i++)
                        if (g.Component(x[i]) != g.Component(x[i - 1])) switches++;
            EssDiagnostic ess = ChainDiagnostics.MeanEss(sequences);
            return new Reading(m, g, ess.Value, switches, all.Min(), all.Max());
        }

        public string Row(string format)
            => $"{Moments.Mean.ToString(format)} | {Moments.Sd.ToString(format)} | {Moments.Skewness:F2} | {Moments.Kurtosis:F2} | {Moments.Bimodality:F3} | " +
               $"{Mixture.DeltaBic:F1} | {Mixture.Minority:F3} | {Mixture.Separation:F2} | {Mixture.Mean1.ToString(format)}, {Mixture.Mean2.ToString(format)} | " +
               $"{Mixture.Sd1.ToString(format)}, {Mixture.Sd2.ToString(format)} | {Ess:F0}";
    }

    // A sample's mean, sd, sample-size-adjusted skewness G₁ and excess kurtosis G₂, and Sarle's bimodality coefficient
    // (G₁² + 1)/(G₂ + 3(n − 1)²/((n − 2)(n − 3))), which exceeds 5/9, the uniform's, for bimodal or heavily skewed samples.
    private readonly record struct Moments(int N, double Mean, double Sd, double Skewness, double Kurtosis)
    {
        public double Bimodality => (Skewness * Skewness + 1.0) / (Kurtosis + 3.0 * (N - 1.0) * (N - 1.0) / ((N - 2.0) * (N - 3.0)));

        public static Moments Of(double[] x)
        {
            int n = x.Length;
            double mean = x.Average(), m2 = 0.0, m3 = 0.0, m4 = 0.0;
            foreach (double v in x)
            {
                double d = v - mean, d2 = d * d;
                m2 += d2;
                m3 += d2 * d;
                m4 += d2 * d2;
            }
            m2 /= n;
            m3 /= n;
            m4 /= n;
            if (!(m2 > 0.0)) return new Moments(n, mean, 0.0, double.NaN, double.NaN);
            double g1 = m3 / Math.Pow(m2, 1.5), g2 = m4 / (m2 * m2) - 3.0;
            double skew = g1 * Math.Sqrt(n * (n - 1.0)) / (n - 2.0);
            double kurt = ((n + 1.0) * g2 + 6.0) * (n - 1.0) / ((n - 2.0) * (n - 3.0));
            return new Moments(n, mean, Math.Sqrt(m2 * n / (n - 1.0)), skew, kurt);
        }
    }

    // Two Gaussians against one by maximum likelihood: EM from five splits of the sorted sample (at its 10, 25, 50, 75
    // and 90% points), the best kept, each variance floored at 10⁻⁶ of the sample's. ΔBIC = BIC₁ − BIC₂ (2 and 5
    // parameters), positive favouring two; the components ordered by mean; Ashman's D, the means' separation over the sds'
    // root mean square, above 2 for a clean separation.
    private readonly record struct TwoGaussians(double DeltaBic, double Weight, double Mean1, double Sd1, double Mean2, double Sd2)
    {
        public double Minority => Math.Min(Weight, 1.0 - Weight);

        public double Separation => Math.Abs(Mean2 - Mean1) / Math.Sqrt((Sd1 * Sd1 + Sd2 * Sd2) / 2.0);

        // The component more responsible for x.
        public int Component(double x)
            => Math.Log(Weight) + LogNormal(x, Mean1, Sd1 * Sd1) >= Math.Log(1.0 - Weight) + LogNormal(x, Mean2, Sd2 * Sd2) ? 0 : 1;

        public static TwoGaussians Fit(double[] x)
        {
            int n = x.Length;
            double mean = x.Average(), variance = x.Sum(v => (v - mean) * (v - mean)) / n;
            if (!(variance > 0.0)) return new TwoGaussians(double.NaN, double.NaN, mean, 0.0, mean, 0.0);
            double one = -0.5 * n * (Math.Log(2.0 * Math.PI * variance) + 1.0), floor = 1e-6 * variance;
            double[] sorted = [.. x.Order()];
            double best = double.NegativeInfinity;
            TwoGaussians fitted = default;
            foreach (double q in new[] { 0.10, 0.25, 0.50, 0.75, 0.90 })
            {
                int k = Math.Clamp((int)Math.Round(q * n), 2, n - 2);
                double w = (double)k / n, m1 = sorted[..k].Average(), m2 = sorted[k..].Average();
                double v1 = Math.Max(floor, sorted[..k].Sum(v => (v - m1) * (v - m1)) / k), v2 = Math.Max(floor, sorted[k..].Sum(v => (v - m2) * (v - m2)) / (n - k));
                double logLikelihood = double.NegativeInfinity;
                for (int iteration = 0; iteration < 1000; iteration++)
                {
                    double s1 = 0.0, x1 = 0.0, xx1 = 0.0, s2 = 0.0, x2 = 0.0, xx2 = 0.0, total = 0.0;
                    double lw1 = Math.Log(w), lw2 = Math.Log(1.0 - w);
                    foreach (double v in x)
                    {
                        double a = lw1 + LogNormal(v, m1, v1), b = lw2 + LogNormal(v, m2, v2), hi = Math.Max(a, b);
                        double lse = hi + Math.Log(Math.Exp(a - hi) + Math.Exp(b - hi)), r = Math.Exp(a - lse);
                        total += lse;
                        s1 += r; x1 += r * v; xx1 += r * v * v;
                        s2 += 1.0 - r; x2 += (1.0 - r) * v; xx2 += (1.0 - r) * v * v;
                    }
                    bool done = total - logLikelihood < 1e-7;
                    logLikelihood = total;
                    if (done || s1 < 1e-9 || s2 < 1e-9) break;
                    w = s1 / n;
                    m1 = x1 / s1;
                    m2 = x2 / s2;
                    v1 = Math.Max(floor, xx1 / s1 - m1 * m1);
                    v2 = Math.Max(floor, xx2 / s2 - m2 * m2);
                }
                if (logLikelihood > best && w > 0.0 && w < 1.0)
                {
                    best = logLikelihood;
                    fitted = m1 <= m2
                        ? new TwoGaussians(0.0, w, m1, Math.Sqrt(v1), m2, Math.Sqrt(v2))
                        : new TwoGaussians(0.0, 1.0 - w, m2, Math.Sqrt(v2), m1, Math.Sqrt(v1));
                }
            }
            double bic1 = -2.0 * one + 2.0 * Math.Log(n), bic2 = -2.0 * best + 5.0 * Math.Log(n);
            return fitted with { DeltaBic = bic1 - bic2 };
        }

        private static double LogNormal(double x, double mean, double variance)
            => -0.5 * (Math.Log(2.0 * Math.PI * variance) + (x - mean) * (x - mean) / variance);
    }

    private const string Blocks = "▁▂▃▄▅▆▇█";

    // A text histogram over [lo, hi]: "·" for an empty bin, then eighths of the fullest.
    private static string Histogram(IEnumerable<double> values, double lo, double hi, int bins = Bins)
    {
        var counts = new int[bins];
        foreach (double v in values)
            counts[hi > lo ? Math.Clamp((int)((v - lo) / (hi - lo) * bins), 0, bins - 1) : 0]++;
        int max = counts.Max();
        return string.Concat(counts.Select(k => k == 0 ? '·' : Blocks[Math.Min(Blocks.Length - 1, (int)((k - 1) * (double)Blocks.Length / max))]));
    }

    // M's values with their counts: "6:12 7:803 8:1,100".
    private static string Counts(IEnumerable<int> values)
        => string.Join(" ", values.GroupBy(v => v).OrderBy(g => g.Key).Select(g => $"{g.Key}:{g.Count()}"));
}
