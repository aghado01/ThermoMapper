using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The nested conditional reference sampler (stage 7, C3; D09) on N14's ensemble: the pooled cut posterior against the
/// exact equal mixture of its conditionals, each by quadrature; every inner run's verdict kept, a run that does not
/// converge reported and left out of the pool unless asked; the outer loop's draws; and the error finite inner runs
/// leave, measured against the exact mixture as their length varies.
/// </summary>
public sealed class BassCalibrationNestedTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationNestedTests(ITestOutputHelper output) => _out = output;

    private static CalibrationModel Model(CalibrationSampling sampling, ReferenceSampling reference, int categoricals = 0)
    {
        TestEmulator emulator = CalibrationFixture.Emulator(categoricals);
        CalibrationConfig config = CalibrationFixture.Config(emulator, CalibrationFixture.Groups(one: true)) with { Sampling = sampling, Reference = reference };
        return new CalibrationModel(config, CalibrationFixture.Observations(emulator), emulator);
    }

    // The cut posterior's E[θ] in unit coordinates: the equal mixture of the conditionals', each by quadrature.
    private static double ExactMean(CalibrationModel model)
        => Enumerable.Range(0, model.Emulator.DrawCount).Average(s => BassCalibrationSamplerTests.Quadrature(model, model.Target(s), 1.0).MeanTheta);

    private static (double Mean, double Sd) Weighted(IReadOnlyList<CalibrationDraw> draws, double[] weights)
    {
        double[] theta = [.. draws.Select(d => CalibrationFixture.Theta.ToUnit(d.Theta[0]))];
        double total = weights.Sum(), mean = theta.Zip(weights, (t, w) => t * w).Sum() / total;
        double variance = theta.Zip(weights, (t, w) => w * (t - mean) * (t - mean)).Sum() / total;
        return (mean, Math.Sqrt(variance));
    }

    [Fact]
    public void Pooled_IsTheEqualMixtureOfTheConverged_Conditionals()
    {
        // θ and the level are coupled here, and untempered one draw's θ group missed the gate within 6,000 draws (R̂ 1.033):
        // the inner runs take a likelihood-path ladder.
        var sampling = new CalibrationSampling { MasterSeed = 7, BurnIn = 1500, DecorrelationEvery = 25, Tempering = [1.0, 0.6, 0.3, 0.1] };
        var reference = new ReferenceSampling { Chains = 3, MaxSamples = 6000, BatchSize = 500 };
        CalibrationModel model = Model(sampling, reference, categoricals: 1);
        NestedCalibration result = NestedCalibrationSampler.Run(model);

        Assert.Equal(3, result.Runs.Count);
        Assert.Equal([0, 1, 2], result.Runs.Select(r => r.EmulatorDraw).Order());
        foreach (InnerRun run in result.Runs)
        {
            _out.WriteLine($"e_{run.EmulatorDraw}: {run.Verdict.Outcome} after {run.SamplesPerChain} draws per chain " +
                           $"({string.Join(", ", run.Verdict.Groups.Select(g => $"{g.Name} {g.Outcome}"))}), acceptance {run.AcceptanceRate:F2}; failing: " +
                           string.Join(", ", run.Verdict.Groups.SelectMany(g => g.Failing).Select(f => $"{f.Name} R̂ {f.RHat.Value:F3} ESS {f.Ess?.ToString() ?? "-"}")));
            Assert.Equal(ConvergenceOutcome.Converged, run.Verdict.Outcome);
            Assert.Equal(["theta", "variances", "discrepancy", "fit"], run.Verdict.Groups.Select(g => g.Name));
            Assert.Equal(3, run.Chains.Count);
            Assert.All(run.Chains, c => Assert.Equal(run.SamplesPerChain, c.Count));
            Assert.All(run.Chains.SelectMany(c => c), d => Assert.Equal(run.EmulatorDraw, d.EmulatorDraw));
            Assert.Equal((1500, 1), (run.BurnIn, run.Thin));
        }
        Assert.Empty(result.NotConverged);
        (IReadOnlyList<CalibrationDraw> draws, double[] weights) = result.Pooled();
        Assert.Equal(1.0, weights.Sum(), 1e-12);
        (double mean, double sd) = Weighted(draws, weights);
        // Each run's weight is a third, whatever its length.
        foreach (InnerRun run in result.Runs)
            Assert.Equal(1.0 / 3, Enumerable.Range(0, draws.Count).Where(i => draws[i].EmulatorDraw == run.EmulatorDraw).Sum(i => weights[i]), 1e-12);
        double exact = ExactMean(model);
        // The mixture's Monte Carlo error: each run's own, a third of it each, with an integrated autocorrelation of at
        // most 10 scans (the θ walk's, measured in the sampler's tests).
        double mcse = Math.Sqrt(result.Runs.Sum(r => sd * sd / (r.Chains.Sum(c => c.Count) / 10.0)) / 9.0);
        _out.WriteLine($"pooled E[θ] {mean:F4} vs the exact mixture {exact:F4} (sd {sd:F3}, MCSE ≈ {mcse:F4})");
        Assert.InRange(mean - exact, -5 * mcse - 2.0 / 600, 5 * mcse + 2.0 / 600);
        Assert.InRange(draws.Min(d => d.Theta[0]), CalibrationFixture.Theta.Lower, CalibrationFixture.Theta.Upper);
    }

    [Fact]
    public void NotConverged_IsReported_AndNotPooledUnlessAsked()
    {
        // A budget too short for the ESS gate: every run ends NotConverged.
        var sampling = new CalibrationSampling { MasterSeed = 9, BurnIn = 100 };
        CalibrationModel model = Model(sampling, new ReferenceSampling { Chains = 2, MaxSamples = 30, BatchSize = 30, OuterDraws = 2 });
        NestedCalibration result = NestedCalibrationSampler.Run(model);
        Assert.Equal(2, result.Runs.Count);
        Assert.All(result.Runs, r => Assert.Equal(ConvergenceOutcome.NotConverged, r.Verdict.Outcome));
        Assert.Equal(2, result.NotConverged.Count);
        Assert.Empty(result.Pooled().Draws);
        (IReadOnlyList<CalibrationDraw> all, double[] weights) = result.Pooled(includeNotConverged: true);
        Assert.Equal(2 * 2 * 30, all.Count);
        Assert.Equal(1.0, weights.Sum(), 1e-12);
        foreach (InnerRun run in result.Runs)
            foreach (GroupVerdict g in run.Verdict.Groups.Where(g => g.Outcome == ConvergenceOutcome.NotConverged))
                _out.WriteLine($"e_{run.EmulatorDraw}, {g.Name}: {string.Join(", ", g.Failing.Select(f => $"{f.Name} R̂ {f.RHat.Value:F3}"))}");
    }

    [Fact]
    public void OuterDraws_AreASampleOfTheRetainedDraws()
    {
        TestEmulator emulator = CalibrationFixture.Emulator();
        var rng = new Xoshiro256PlusPlus(3);
        Assert.Equal([0, 1, 2], NestedCalibrationSampler.OuterDraws(emulator, null, rng));
        Assert.Equal([0, 1, 2], NestedCalibrationSampler.OuterDraws(emulator, 5, rng));
        var counts = new int[3];
        for (int i = 0; i < 6000; i++)
        {
            int[] two = NestedCalibrationSampler.OuterDraws(emulator, 2, rng);
            Assert.Equal(2, two.Distinct().Count());
            foreach (int d in two) counts[d]++;
        }
        foreach (int c in counts) Assert.InRange(c / 12000.0, 1 / 3.0 - 0.02, 1 / 3.0 + 0.02);
        Assert.Contains("OuterDraws", Assert.Throws<ArgumentException>(() => new ReferenceSampling { OuterDraws = 0 }.Validate()).Message);
        Assert.Contains("Chains", Assert.Throws<ArgumentException>(() => new ReferenceSampling { Chains = 1 }.Validate()).Message);
    }

    /// <summary>An emulator whose draws carry weights (SMC's), over a <see cref="TestEmulator"/>.</summary>
    private sealed class WeightedEmulator(TestEmulator inner, double[] weights) : ICalibrationEmulator
    {
        public IReadOnlyList<CalibrationParameter> Parameters => inner.Parameters;
        public int DrawCount => inner.DrawCount;
        public IReadOnlyList<double>? DrawWeights => weights;
        public int Outputs => inner.Outputs;
        public EmulatorNoise Noise(int draw) => inner.Noise(draw);
        public double[,]? TruncationResiduals() => inner.TruncationResiduals();
        public void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output) => inner.Predict(draw, continuous, levels, output);
    }

    [Fact]
    public void OuterDraws_FollowTheEmulatorsWeights()
    {
        var emulator = new WeightedEmulator(CalibrationFixture.Emulator(), [1.0, 3.0, 0.0]);
        var rng = new Xoshiro256PlusPlus(5);
        var counts = new int[3];
        foreach (int d in NestedCalibrationSampler.OuterDraws(emulator, 8000, rng)) counts[d]++;
        Assert.Equal(0, counts[2]);
        Assert.InRange(counts[1] / 8000.0, 0.75 - 0.02, 0.75 + 0.02);
        Assert.Equal(3, NestedCalibrationSampler.OuterDraws(emulator, null, rng).Length);
        Assert.Throws<ArgumentException>(() => NestedCalibrationSampler.OuterDraws(new WeightedEmulator(CalibrationFixture.Emulator(), [0.0, 0.0, 0.0]), 2, rng));
    }

    [Fact]
    public void Config_WithTheSamplersSettings_RoundTripsThroughJson()
    {
        TestEmulator emulator = CalibrationFixture.Emulator();
        CalibrationConfig config = CalibrationFixture.Config(emulator) with
        {
            Sampling = new CalibrationSampling { MasterSeed = 4, BurnIn = 10, Samples = 20, Thin = 2, DecorrelationEvery = 0, Tempering = [1.0, 0.3, 0.0] },
            Reference = new ReferenceSampling { OuterDraws = 2, Chains = 3, MaxSamples = 50, BatchSize = 25, RHatTarget = 1.05, MinEssPerChain = 10 },
        };
        string json = config.ToJson();
        CalibrationConfig back = CalibrationConfig.FromJson(json);
        Assert.Equal(json, back.ToJson());
        Assert.Equal([1.0, 0.3, 0.0], back.Sampling.Tempering);
        Assert.Equal(2, back.Reference.OuterDraws);
        Assert.Contains("CalibrationConfig.Reference.RHatTarget", Assert.Throws<ArgumentException>(() =>
            CalibrationConfig.FromJson(json.Replace("\"RHatTarget\":1.05", "\"RHatTarget\":1"))).Message);
        Assert.Contains("CalibrationConfig.Sampling.Tempering", Assert.Throws<ArgumentException>(() =>
            CalibrationConfig.FromJson(json.Replace("[1,0.3,0]", "[1,0.3,0.3]"))).Message);
    }

    [Fact]
    public void InnerLength_TheErrorItLeaves()
    {
        // Fixed inner lengths (an unreachable ESS gate disables the early stop, so every run is pooled explicitly), three
        // replicate seeds each, against the exact mixture: the pooled E[θ] and sd as the inner length grows.
        CalibrationModel probe = Model(new CalibrationSampling(), new ReferenceSampling());
        double exact = ExactMean(probe);
        int[] lengths = [25, 100, 400, 1600];
        var rows = new List<(int Length, double Rms, double MeanSd, double Spread)>();
        foreach (int length in lengths)
        {
            var errors = new List<double>();
            var sds = new List<double>();
            for (int rep = 0; rep < 3; rep++)
            {
                var sampling = new CalibrationSampling { MasterSeed = 100 * length + rep, BurnIn = 500, DecorrelationEvery = 25 };
                var reference = new ReferenceSampling { Chains = 2, MaxSamples = length, BatchSize = length, MinEssPerChain = 1e9 };
                NestedCalibration result = NestedCalibrationSampler.Run(Model(sampling, reference));
                (IReadOnlyList<CalibrationDraw> draws, double[] weights) = result.Pooled(includeNotConverged: true);
                (double mean, double sd) = Weighted(draws, weights);
                errors.Add(mean - exact);
                sds.Add(sd);
            }
            double rms = Math.Sqrt(errors.Average(e => e * e));
            rows.Add((length, rms, sds.Average(), errors.Max() - errors.Min()));
        }
        _out.WriteLine($"exact mixture E[θ] = {exact:F4} (unit coordinates); inner length → RMS error of the pooled E[θ] over 3 seeds, mean pooled sd, spread");
        foreach ((int length, double rms, double sd, double spread) in rows)
            _out.WriteLine($"  {length,5} draws per chain × 2 chains × 3 draws: RMS {rms:F4}, sd {sd:F4}, spread {spread:F4}");
        // Each length's error stays within what its draws can resolve: 5 MCSE at an integrated autocorrelation of 10.
        foreach ((int length, double rms, double sd, _) in rows)
            Assert.InRange(rms, 0, 5 * sd / Math.Sqrt(3 * 2 * length / 10.0) + 2.0 / 600);
    }
}
