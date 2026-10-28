using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers.Evidence;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Samplers.Tests;

/// <summary>
/// Evidence over a tempering ladder (N11), on targets whose tempered normalizer is known at every b:
/// <list type="bullet">
/// <item>a conjugate Normal prior and likelihood, sampled exactly at each rung, so the estimators' bias and their reported
/// errors are read apart from any sampler's error;</item>
/// <item>an enumerable target whose likelihood is not linear in b, sampled by traced tempering ladders whose local moves
/// leave each rung's trace autocorrelated, so the effective sample size enters the error.</item>
/// </list>
/// </summary>
public sealed class LadderEvidenceTests
{
    private readonly ITestOutputHelper _out;
    public LadderEvidenceTests(ITestOutputHelper output) => _out = output;

    // x ~ N(0, 3²) and one observation y = 1.5 ~ N(x, 0.25²): a likelihood far narrower than the prior, so the mean slope
    // climbs steeply near b = 0.
    private static readonly Conjugate Normal1 = new(3.0, 0.25, 1.5);

    // Friel and Pettitt's (2008) schedule, b_k = (k/K)⁵, ascending.
    private static double[] PowerLadder(int intervals)
    {
        var b = new double[intervals + 1];
        for (int k = 0; k <= intervals; k++) b[k] = Math.Pow((double)k / intervals, 5.0);
        return b;
    }

    [Fact]
    public void PathSampling_DiscretizationError_FallsWithRungs_AndItsDiagnosticTracksIt()
    {
        // Each rung's exact law as a three-point Gauss–Hermite population, exact for the slope's mean and variance (degree
        // ≤ 4), so the trapezoid's error is the discretization's alone.
        double truth = Normal1.LogZ(1.0) - Normal1.LogZ(0.0);
        double previous = double.PositiveInfinity;
        foreach (int intervals in new[] { 4, 8, 16, 32, 64 })
        {
            double[] b = PowerLadder(intervals);
            var rungs = new RungDraws[b.Length];
            for (int k = 0; k < b.Length; k++) rungs[k] = Normal1.Quadrature(b[k]);
            LadderEstimate e = LadderEvidence.PathSampling(b, rungs);
            for (int k = 0; k < b.Length; k++) Assert.Equal(Normal1.MeanSlope(b[k]), e.Rungs[k].Value, 1e-9);

            double error = Math.Abs(e.LogRatio - truth), corrected = Math.Abs(e.LogRatio + e.DiscretizationError - truth);
            _out.WriteLine($"K = {intervals,2}: trapezoid error {e.LogRatio - truth,10:G4}, diagnostic {e.DiscretizationError,10:G4}, corrected error {e.LogRatio + e.DiscretizationError - truth,10:G4}");
            Assert.True(error < previous, $"K = {intervals}: the error {error:G4} should fall below {previous:G4}");
            if (intervals >= 8)   // at four intervals the correction overshoots: 2.53 against an error of 1.44
                Assert.True(corrected < 0.2 * error, $"K = {intervals}: the corrected error {corrected:G4} vs the trapezoid's {error:G4}");
            previous = error;
        }
    }

    [Fact]
    public void ExactDraws_EstimatorsAgreeWithTheTruth_AndTheirErrorsCoverAtTheNominalRate()
    {
        // 400 replicates of 16 intervals, 1000 exact draws per rung. Stepping-stone's estimand is the truth itself; path
        // sampling's is the trapezoid of the exact mean slopes, its discretization bias reported apart.
        const int replicates = 400, draws = 1000;
        double[] b = PowerLadder(16);
        double truth = Normal1.LogZ(1.0) - Normal1.LogZ(0.0);
        double trapezoid = 0.0;
        for (int k = 0; k + 1 < b.Length; k++) trapezoid += 0.5 * (b[k + 1] - b[k]) * (Normal1.MeanSlope(b[k]) + Normal1.MeanSlope(b[k + 1]));

        var ss = new Calibration(truth);
        var ps = new Calibration(trapezoid);
        int[] seeds = SeedTree.Derive(20261027, replicates);
        for (int r = 0; r < replicates; r++)
        {
            var rng = new Xoshiro256PlusPlus(seeds[r]);
            var slopes = new RungDraws[b.Length];
            var ratios = new RungDraws?[b.Length];
            for (int k = 0; k < b.Length; k++)
            {
                var slope = new double[draws];
                var ratio = new double[draws];
                for (int i = 0; i < draws; i++)
                {
                    double l = Normal1.LogLikelihood(Normal1.Draw(b[k], rng));
                    slope[i] = l;
                    ratio[i] = k + 1 < b.Length ? (b[k + 1] - b[k]) * l : 0.0;
                }
                slopes[k] = RungDraws.FromChains(slope);
                ratios[k] = k + 1 < b.Length ? RungDraws.FromChains(ratio) : null;
            }
            ss.Add(LadderEvidence.SteppingStone(b, ratios));
            ps.Add(LadderEvidence.PathSampling(b, slopes));
        }

        _out.WriteLine($"truth {truth:F5}; path sampling's estimand (exact trapezoid) {trapezoid:F5}, its discretization bias {trapezoid - truth:F5}");
        _out.WriteLine($"stepping-stone: {ss}");
        _out.WriteLine($"path sampling:  {ps}");
        ss.AssertCalibrated("stepping-stone");
        ps.AssertCalibrated("path sampling");
        Assert.True(Math.Abs(ss.Bias) < 0.25 * ss.MeanError, $"stepping-stone's bias {ss.Bias:G4} against its error {ss.MeanError:G4}");
    }

    [Fact]
    public void TracedLadders_OnAnEnumerableTarget_GiveCalibratedEvidence()
    {
        // 200 replicates of two traced ladders (a reference rung at b = 0 and seven geometric rungs), 3000 scans each.
        // Unit-step moves leave the rung traces autocorrelated; the errors divide by effective sample sizes. The likelihood
        // is not linear in b, so the trapezoid's bias is far larger than its diagnostic reads (the curvature term).
        const int replicates = 200, scans = 3000;
        double[] betas = ParallelTempering<int>.ReferenceLadder(8, 0.01);
        double truth = Bumps.LogZ(1.0) - Bumps.LogZ(0.0), trapezoid = 0.0;
        for (int k = 0; k + 1 < betas.Length; k++)
            trapezoid += 0.5 * (betas[k] - betas[k + 1]) * (Bumps.MeanSlope(betas[k]) + Bumps.MeanSlope(betas[k + 1]));

        var ss = new Calibration(truth);
        var ps = new Calibration(trapezoid);
        double smallestEssFraction = 1.0, diagnostic = 0.0;
        int[] seeds = SeedTree.Derive(1027, replicates);
        for (int r = 0; r < replicates; r++)
        {
            int[] ladderSeeds = SeedTree.Derive(seeds[r], 2);
            var ladders = new ParallelTempering<int>[2];
            for (int l = 0; l < 2; l++)
            {
                ladders[l] = Bumps.Ladder(betas, ladderSeeds[l]);
                for (int i = 0; i < 500; i++) ladders[l].Step();
                ladders[l].TraceRungs();
                for (int i = 0; i < scans; i++) ladders[l].Step();
            }
            LadderEstimate s = LadderEvidence.FromLadders(ladders);
            LadderEstimate p = LadderEvidence.FromLadders(ladders, LadderEstimator.PathSampling);
            ss.Add(s);
            ps.Add(p);
            diagnostic += p.DiscretizationError / replicates;
            foreach (RungEstimate rung in p.Rungs) smallestEssFraction = Math.Min(smallestEssFraction, rung.EffectiveSampleSize / (2.0 * scans));
        }

        _out.WriteLine($"truth {truth:F5}; exact trapezoid {trapezoid:F5}; mean discretization diagnostic {diagnostic:F5}; smallest slope ESS / draws {smallestEssFraction:F3}");
        _out.WriteLine($"stepping-stone: {ss}");
        _out.WriteLine($"path sampling:  {ps}");
        Assert.True(smallestEssFraction < 0.5, $"the rung traces should be autocorrelated (ESS/draws {smallestEssFraction:F3})");
        ss.AssertCalibrated("stepping-stone", coverageFloor: 0.88);
        ps.AssertCalibrated("path sampling", coverageFloor: 0.88);
    }

    [Fact]
    public void SmcStages_ReproduceTheRunsNormalizer_AndEstimateTheEvidence()
    {
        // The conjugate Normal tempered by SMC from its prior: stepping-stone over the traced stages is the run's own
        // normalizer, term for term, and path sampling over the same stages estimates the same evidence.
        double truth = Normal1.LogZ(1.0) - Normal1.LogZ(0.0);
        SmcResult<double[]> result = SmcSampler.RunTempered(ConjugateTarget, TemperingPath.Likelihood, 0.9, new PriorParticles(),
            new RandomWalkKernel(), new SmcOptions { Particles = 2000, MutationSteps = 5, Resampling = ResamplingScheme.Multinomial, TraceStages = true }, seed: 31);
        IReadOnlyList<SmcStageTrace> stages = result.StageTraces!;
        Assert.Equal(result.Stages.Count, stages.Count);
        for (int k = 0; k < stages.Count; k++)
        {
            Assert.Equal(result.Stages[k].Temperature, stages[k].Temperature);
            Assert.Equal(k + 1 < stages.Count, stages[k].LogRatios is not null);
        }

        LadderEstimate ss = LadderEvidence.FromSmc(result), ps = LadderEvidence.FromSmc(result, LadderEstimator.PathSampling);
        _out.WriteLine($"{stages.Count} stages; truth {truth:F5}; SMC normalizer {result.LogNormalizerRatio:F5}; " +
                       $"stepping-stone {ss.LogRatio:F5} ± {ss.MonteCarloError:F5}; path sampling {ps.LogRatio:F5} ± {ps.MonteCarloError:F5}, diagnostic {ps.DiscretizationError:F5}");
        Assert.Equal(result.LogNormalizerRatio, ss.LogRatio, 1e-10);
        Assert.True(Math.Abs(ss.LogRatio - truth) < 4.0 * ss.MonteCarloError, $"stepping-stone {ss.LogRatio:F5} against {truth:F5}");
        Assert.True(Math.Abs(ps.LogRatio + ps.DiscretizationError - truth) < 4.0 * ps.MonteCarloError + 0.25 * Math.Abs(ps.DiscretizationError),
                    $"path sampling {ps.LogRatio:F5} (diagnostic {ps.DiscretizationError:F5}) against {truth:F5}");
        Assert.Throws<InvalidOperationException>(() => LadderEvidence.FromSmc(SmcSampler.RunTempered(ConjugateTarget, TemperingPath.Likelihood, 0.9,
            new PriorParticles(), new RandomWalkKernel(), new SmcOptions { Particles = 200 }, seed: 31)));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(5)]
    public void SmcNormalizerError_SeesSharedAncestry_AndCovers(int moves)
    {
        // The probe that found the stages' weights blind to shared ancestry (Kish: coverage 0.795 at one move per stage),
        // rerun with the genealogy error: the conjugate Normal tempered from its prior at a conditional ESS of 0.9, 1000
        // particles, multinomial resampling below half the population, 200 seeds. The Kish reading is logged beside it.
        double truth = Normal1.LogZ(1.0) - Normal1.LogZ(0.0);
        var genealogy = new Calibration(truth);
        var kish = new Calibration(truth);
        foreach (int seed in SeedTree.Derive(2027, 200))
        {
            SmcResult<double[]> result = SmcSampler.RunTempered(ConjugateTarget, TemperingPath.Likelihood, 0.9, new PriorParticles(),
                new RandomWalkKernel(), new SmcOptions { Particles = 1000, MutationSteps = moves, Resampling = ResamplingScheme.Multinomial, TraceStages = true }, seed);
            LadderEstimate ss = LadderEvidence.FromSmc(result);
            Assert.Equal(result.LogNormalizerError, ss.MonteCarloError);
            genealogy.Add(ss);
            double weightsOnly = 0.0;   // the stages' own errors combined as independent: the Kish reading
            foreach (RungEstimate rung in ss.Rungs) weightsOnly += rung.MonteCarloError * rung.MonteCarloError;
            kish.Add(ss.LogRatio, Math.Sqrt(weightsOnly));
        }
        _out.WriteLine($"{moves} moves per stage. Genealogy: {genealogy}. Kish: {kish}.");
        genealogy.AssertCalibrated("the genealogy error");
    }

    [Fact]
    public void SmcNormalizerError_IsNaN_AfterAResamplingItDoesNotDescribe()
    {
        // Systematic resampling's offspring are not multinomial's, so the estimate stops at the first such resampling,
        // and a run resumed from just after it keeps it stopped whatever its options.
        var options = new SmcOptions { Particles = 300, MutationSteps = 1, Resampling = ResamplingScheme.Systematic, TraceStages = true };
        SmcRun<double[]> run = SmcRun<double[]>.Start(Stage(0.0), new PriorParticles(), new RandomWalkKernel(), options, seed: 9);
        int first = -1;
        SmcCheckpoint<double[]>? atFirst = null;
        for (int s = 1; s <= 10; s++)
        {
            run.Advance(Stage(s / 10.0));
            SmcStageReport report = run.Reports[^1];
            if (first < 0)
            {
                Assert.True(double.IsFinite(report.NormalizerRelativeVariance), $"stage {s}, before any resampling");
                if (report.Resampled)
                {
                    first = s;
                    atFirst = run.Checkpoint();
                }
            }
            else Assert.True(double.IsNaN(report.NormalizerRelativeVariance), $"stage {s}, after the first resampling at {first}");
        }
        Assert.True(first is > 0 and < 10, $"the case should resample before its last stage, not at {first}");
        Assert.True(double.IsNaN(run.Result().LogNormalizerError));
        Assert.True(double.IsNaN(LadderEvidence.FromSmc(run.Result()).MonteCarloError));

        Assert.False(atFirst!.GenealogyHolds);
        Assert.True(double.IsFinite(atFirst.Reports[^1].NormalizerRelativeVariance));   // read before that stage's resampling
        SmcRun<double[]> resumed = SmcRun<double[]>.Resume(atFirst, new PriorParticles(), new RandomWalkKernel(),
            options with { Resampling = ResamplingScheme.Multinomial });
        resumed.Advance(Stage(1.0));
        Assert.True(double.IsNaN(resumed.Reports[^1].NormalizerRelativeVariance));

        static SmcStage<double[]> Stage(double b) => new(ConjugateTarget, new Temperature(b, TemperingPath.Likelihood));
    }

    [Fact]
    public void SmcStageTraces_ResumeExactly()
    {
        var options = new SmcOptions { Particles = 300, MutationSteps = 2, TraceStages = true };
        SmcRun<double[]> straight = SmcRun<double[]>.Start(Stage(0.0), new PriorParticles(), new RandomWalkKernel(), options, seed: 8);
        SmcRun<double[]> first = SmcRun<double[]>.Start(Stage(0.0), new PriorParticles(), new RandomWalkKernel(), options, seed: 8);
        foreach (double b in new[] { 0.01, 0.05 }) { straight.Advance(Stage(b)); first.Advance(Stage(b)); }
        SmcCheckpoint<double[]> checkpoint = first.Checkpoint();
        Assert.Throws<ArgumentException>(() => SmcRun<double[]>.Resume(checkpoint, new PriorParticles(), new RandomWalkKernel(),
                                                                        options with { TraceStages = false }));
        SmcRun<double[]> resumed = SmcRun<double[]>.Resume(checkpoint, new PriorParticles(), new RandomWalkKernel(), options);
        foreach (double b in new[] { 0.2, 1.0 }) { straight.Advance(Stage(b)); resumed.Advance(Stage(b)); }
        Assert.Equal(straight.StageTraces!.Count, resumed.StageTraces!.Count);
        for (int k = 0; k < straight.StageTraces.Count; k++)
        {
            SmcStageTrace a = straight.StageTraces[k], c = resumed.StageTraces[k];
            Assert.Equal(a.Temperature, c.Temperature);
            Assert.Equal(a.Weights, c.Weights);
            Assert.Equal(a.Slopes, c.Slopes);
            Assert.Equal(a.LogRatios, c.LogRatios);
        }
        Assert.Equal(LadderEvidence.FromSmc(straight.Result()).LogRatio, LadderEvidence.FromSmc(resumed.Result()).LogRatio);

        static SmcStage<double[]> Stage(double b) => new(ConjugateTarget, new Temperature(b, TemperingPath.Likelihood));
    }

    // Normal1 on real-vector particles, for SMC: x[0] ~ N(0, 3²), one observation 1.5 ~ N(x, 0.25²), separable.
    private static readonly ISamplingTarget<double[]> ConjugateTarget = SamplingTarget.Separable<double[]>(
        x => -0.5 * Math.Log(2.0 * Math.PI * 9.0) - x[0] * x[0] / 18.0, x => Normal1.LogLikelihood(x[0]));

    private sealed class PriorParticles : ISmcParticles<double[]>
    {
        public double[] Draw(Xoshiro256PlusPlus rng) => [3.0 * Normal.SampleStandard(rng)];
        public double[] Fork(double[] state) => (double[])state.Clone();
    }

    [Fact]
    public void TraceRungs_RecordsEachRungsSlopeAndLogRatio_AfterTheScansSwaps()
    {
        double[] betas = ParallelTempering<int>.ReferenceLadder(5, 0.05);
        ParallelTempering<int> ladder = Bumps.Ladder(betas, 3);
        for (int i = 0; i < 50; i++) ladder.Step();
        Assert.Null(ladder.RungTraces);
        ladder.TraceRungs();
        Assert.True(ladder.TracesSlopes);
        Assert.Throws<InvalidOperationException>(() => ladder.AdaptLadder());
        for (int i = 0; i < 200; i++)
        {
            ladder.Step();
            for (int r = 0; r < betas.Length; r++)
            {
                ReversibleJumpChain<int> rung = ladder.Replicas[r];
                double[] slope = ladder.RungTraces![r].Draws(RungFunctional.Slope);
                double[] ratio = ladder.RungTraces[r].Draws(RungFunctional.LogRatio);
                Assert.Equal(i + 1, slope.Length);
                Assert.Equal(Bumps.Slope(rung.Current, betas[r]), slope[^1], 1e-12);
                if (r == 0) Assert.True(double.IsNaN(ratio[^1]));
                else Assert.Equal(Bumps.LogLikelihood(rung.Current, betas[r - 1]) - Bumps.LogLikelihood(rung.Current, betas[r]), ratio[^1], 1e-12);
            }
        }
    }

    [Fact]
    public void TracedLadder_ResumesWithItsTraces()
    {
        double[] betas = ParallelTempering<int>.ReferenceLadder(5, 0.05);
        ParallelTempering<int> straight = Bumps.Ladder(betas, 11), first = Bumps.Ladder(betas, 11);
        straight.TraceRungs();
        first.TraceRungs();
        for (int i = 0; i < 300; i++) { straight.Step(); first.Step(); }
        ParallelTempering<int> resumed = ParallelTempering<int>.Restore(Bumps.Moves, new Bumps(), first.Snapshot(), new Bumps.Prior());
        for (int i = 0; i < 300; i++) { straight.Step(); resumed.Step(); }
        Assert.True(resumed.TracesSlopes);
        for (int r = 0; r < betas.Length; r++)
            for (int f = 0; f < 2; f++)
                Assert.Equal(straight.RungTraces![r].Draws(f), resumed.RungTraces![r].Draws(f));
        Assert.Equal(LadderEvidence.FromLadders([straight]).LogRatio, LadderEvidence.FromLadders([resumed]).LogRatio);
    }

    [Fact]
    public void ATargetWithoutASlope_LeavesSteppingStone()
    {
        // A likelihood-path target whose evaluations refuse the slope: the traces hold NaN there and path sampling refuses;
        // stepping-stone reads ℓ_b alone. On the powered path the slope is ln prior + ℓ₁, which every target has.
        double[] betas = [1.0, 0.5, 0.1, 0.0];
        foreach (TemperingPath path in new[] { TemperingPath.Likelihood, TemperingPath.Powered })
        {
            double[] ladderBetas = path == TemperingPath.Powered ? [1.0, 0.5, 0.1] : betas;
            var ladder = new ParallelTempering<int>(Bumps.Moves, new SlopeFree(), 5, ladderBetas, path, masterSeed: 2);
            ladder.TraceRungs();
            for (int i = 0; i < 100; i++) ladder.Step();
            Assert.Equal(path == TemperingPath.Powered, ladder.TracesSlopes);
            Assert.True(double.IsFinite(LadderEvidence.FromLadders([ladder]).LogRatio));
            if (path == TemperingPath.Likelihood)
            {
                Assert.True(double.IsNaN(ladder.RungTraces![1].Draws(RungFunctional.Slope)[0]));
                Assert.Throws<InvalidOperationException>(() => LadderEvidence.FromLadders([ladder], LadderEstimator.PathSampling));
            }
            else Assert.True(double.IsFinite(LadderEvidence.FromLadders([ladder], LadderEstimator.PathSampling).LogRatio));
        }
        Assert.Throws<InvalidOperationException>(() => LadderEvidence.FromLadders([Bumps.Ladder(betas, 1)]));   // untraced
    }

    [Fact]
    public void Inputs_AreValidated_AndTheLadderMayRunEitherWay()
    {
        var rng = new Xoshiro256PlusPlus(5);
        double[] b = [0.0, 0.2, 0.6, 1.0];
        var slopes = new RungDraws[b.Length];
        var ratios = new RungDraws?[b.Length];
        for (int k = 0; k < b.Length; k++)
        {
            var draws = new double[200];
            for (int i = 0; i < draws.Length; i++) draws[i] = Normal1.LogLikelihood(Normal1.Draw(b[k], rng));
            slopes[k] = RungDraws.FromChains(draws);
            ratios[k] = k + 1 < b.Length ? RungDraws.FromChains(Array.ConvertAll(draws, l => (b[k + 1] - b[k]) * l)) : null;
        }
        double[] down = [1.0, 0.6, 0.2, 0.0];
        RungDraws[] slopesDown = [slopes[3], slopes[2], slopes[1], slopes[0]];
        RungDraws?[] ratiosDown = [null, ratios[2], ratios[1], ratios[0]];
        Assert.Equal(LadderEvidence.PathSampling(b, slopes), LadderEvidence.PathSampling(down, slopesDown), Equivalent);
        Assert.Equal(LadderEvidence.SteppingStone(b, ratios), LadderEvidence.SteppingStone(down, ratiosDown), Equivalent);

        Assert.Throws<ArgumentException>(() => LadderEvidence.PathSampling([0.0, 0.6, 0.2, 1.0], slopes));
        Assert.Throws<ArgumentException>(() => LadderEvidence.PathSampling([0.0, 1.0], slopes));
        Assert.Throws<ArgumentException>(() => LadderEvidence.SteppingStone(b, [null, ratios[1], ratios[2], null]));
        Assert.Throws<ArgumentException>(() => LadderEvidence.PathSampling(b, [slopes[0], slopes[1], slopes[2], RungDraws.FromChains([1.0, double.NaN])]));
        Assert.Throws<ArgumentException>(() => RungDraws.FromChains([1.0, 2.0], [1.0]));
        Assert.Throws<ArgumentException>(() => RungDraws.FromPopulation([1.0, 2.0], [0.5]));
        Assert.Throws<ArgumentException>(() => RungDraws.FromPopulation([1.0, 2.0], [0.0, 0.0]));
    }

    private static readonly EstimateComparer Equivalent = new();

    // The same estimate up to rounding: the ladder's direction only reverses the order of the sums.
    private sealed class EstimateComparer : IEqualityComparer<LadderEstimate>
    {
        public bool Equals(LadderEstimate? x, LadderEstimate? y)
            => x is not null && y is not null && Math.Abs(x.LogRatio - y.LogRatio) < 1e-12 * Math.Max(1.0, Math.Abs(x.LogRatio))
               && Math.Abs(x.MonteCarloError - y.MonteCarloError) < 1e-12 * Math.Max(1.0, x.MonteCarloError)
               && x.Rungs.Count == y.Rungs.Count;

        public int GetHashCode(LadderEstimate obj) => 0;
    }

    // z = (estimate − estimand)/error across replicates: the coverage of ±1.96 errors and the mean of z².
    private sealed class Calibration(double estimand)
    {
        private readonly List<double> _z = new(), _estimates = new(), _errors = new();

        public void Add(LadderEstimate e) => Add(e.LogRatio, e.MonteCarloError);

        public void Add(double estimate, double error)
        {
            _estimates.Add(estimate);
            _errors.Add(error);
            _z.Add((estimate - estimand) / error);
        }

        public double Coverage => _z.FindAll(z => Math.Abs(z) <= 1.96).Count / (double)_z.Count;
        public double MeanSquare { get { double s = 0.0; foreach (double z in _z) s += z * z; return s / _z.Count; } }
        public double Bias { get { double s = 0.0; foreach (double v in _estimates) s += v; return s / _estimates.Count - estimand; } }
        public double MeanError { get { double s = 0.0; foreach (double v in _errors) s += v; return s / _errors.Count; } }

        public void AssertCalibrated(string name, double coverageFloor = 0.91)
        {
            Assert.True(Coverage >= coverageFloor && Coverage <= 0.99, $"{name}: coverage {Coverage:F3} of the nominal 0.95");
            Assert.True(MeanSquare is > 0.7 and < 1.3, $"{name}: mean z² {MeanSquare:F3}, where a calibrated error gives 1");
        }

        public override string ToString()
            => $"coverage {Coverage:F3}, mean z² {MeanSquare:F3}, bias {Bias:G4}, mean error {MeanError:G4} over {_z.Count}";
    }

    // The conjugate Normal on the likelihood path: q_b(x) = N(x; 0, s0²)·N(y; x, σ²)^b, complete.
    private sealed record Conjugate(double PriorSd, double Sigma, double Y)
    {
        private double Lambda(double b) => b / (Sigma * Sigma);
        private double Precision(double b) => 1.0 / (PriorSd * PriorSd) + Lambda(b);
        private double Mean(double b) => Lambda(b) * Y / Precision(b);
        private double LogNormalizer => -0.5 * Math.Log(2.0 * Math.PI * Sigma * Sigma);

        public double LogLikelihood(double x) => LogNormalizer - (x - Y) * (x - Y) / (2.0 * Sigma * Sigma);

        public double LogZ(double b)
        {
            double t0 = 1.0 / (PriorSd * PriorSd), lambda = Lambda(b);
            return b * LogNormalizer + 0.5 * Math.Log(t0 / (t0 + lambda)) - 0.5 * Y * Y * t0 * lambda / (t0 + lambda);
        }

        public double MeanSlope(double b)
        {
            double offset = Mean(b) - Y;
            return LogNormalizer - (1.0 / Precision(b) + offset * offset) / (2.0 * Sigma * Sigma);
        }

        public double Draw(double b, Xoshiro256PlusPlus rng) => Mean(b) + Normal.SampleStandard(rng) / Math.Sqrt(Precision(b));

        public RungDraws Quadrature(double b)
        {
            double m = Mean(b), h = Math.Sqrt(3.0 / Precision(b));
            return RungDraws.FromPopulation([LogLikelihood(m - h), LogLikelihood(m), LogLikelihood(m + h)], [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0]);
        }
    }

    // States 0..23 under a uniform prior, with ℓ_b(k) = 3·ln(1 + b·A_k) − b·C_k: two bumps, not linear in b, whose slope
    // the evaluation supplies in closed form.
    private sealed class Bumps : ISamplingTarget<int>
    {
        private const int States = 24;
        private static readonly double[] A = new double[States], C = new double[States];

        static Bumps()
        {
            for (int k = 0; k < States; k++)
            {
                A[k] = 40.0 * Math.Exp(-(k - 5) * (k - 5) / 6.0) + 90.0 * Math.Exp(-(k - 18) * (k - 18) / 4.0);
                C[k] = 0.2 * k;
            }
        }

        public static double LogLikelihood(int k, double b) => 3.0 * Math.Log(1.0 + b * A[k]) - b * C[k];
        public static double Slope(int k, double b) => 3.0 * A[k] / (1.0 + b * A[k]) - C[k];

        public static double LogZ(double b)
        {
            double max = double.NegativeInfinity;
            for (int k = 0; k < States; k++) max = Math.Max(max, LogLikelihood(k, b));
            double sum = 0.0;
            for (int k = 0; k < States; k++) sum += Math.Exp(LogLikelihood(k, b) - max);
            return max + Math.Log(sum / States);
        }

        public static double MeanSlope(double b)
        {
            double logZ = LogZ(b), mean = 0.0;
            for (int k = 0; k < States; k++) mean += Math.Exp(LogLikelihood(k, b) - logZ) / States * Slope(k, b);
            return mean;
        }

        public static IReadOnlyList<IRjMove<int>> Moves(int replica) => [new UnitStep("up", "down", +1), new UnitStep("down", "up", -1)];

        public static ParallelTempering<int> Ladder(double[] betas, int seed)
            => new(Moves, new Bumps(), 5, betas, TemperingPath.Likelihood, seed, reference: betas[^1] == 0.0 ? new Prior() : null);

        public ITargetEvaluation Evaluate(int k) => new Evaluation(k);

        private sealed class Evaluation(int k) : ITargetEvaluation
        {
            private bool Inside => k is >= 0 and < States;
            public double LogPrior => Inside ? -Math.Log(States) : double.NegativeInfinity;
            public double LogLikelihood(double b) => Inside ? Bumps.LogLikelihood(k, b) : double.NegativeInfinity;
            public double LogLikelihoodSlope(double b) => Inside ? Slope(k, b) : double.NaN;
        }

        // Exact draws from the uniform prior, the likelihood path's b = 0.
        public sealed class Prior : IReferenceSampler<int>
        {
            public TemperingPath Path => TemperingPath.Likelihood;
            public long Attempts { get; private set; }

            public (int State, ITargetEvaluation Evaluation) Draw(Xoshiro256PlusPlus rng)
            {
                Attempts++;
                int k = Math.Min(States - 1, (int)(rng.NextDouble() * States));
                return (k, new Evaluation(k));
            }
        }
    }

    // Bumps without the slope: its evaluations keep the interface's refusal.
    private sealed class SlopeFree : ISamplingTarget<int>
    {
        private readonly Bumps _inner = new();
        public ITargetEvaluation Evaluate(int k) => new Plain(_inner.Evaluate(k));

        private sealed class Plain(ITargetEvaluation inner) : ITargetEvaluation
        {
            public double LogPrior => inner.LogPrior;
            public double LogLikelihood(double b) => inner.LogLikelihood(b);
        }
    }

    private sealed class UnitStep(string key, string reverseKey, int delta) : IRjMove<int>
    {
        public string Key => key;
        public string ReverseKey => reverseKey;
        public double Weight(int state) => 1.0;
        public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng) => new(current + delta, 0.0, 0.0);
    }
}
