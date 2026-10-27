using System;
using System.Collections.Generic;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Smc;
using Xunit;

namespace Maths.Samplers.Tests;

/// <summary>
/// The SMC core against references that do not come from it. Weight arithmetic at extreme scales; each resampling
/// scheme's offspring counts average N·W_i, and systematic and residual keep their deterministic parts. On a conjugate
/// Gaussian the adaptive tempered run recovers the closed-form posterior and evidence; on a bimodal likelihood it
/// recovers the modes' masses. Without resampling or moves the evidence telescopes to the plain importance-sampling
/// estimate exactly. Results do not depend on parallelism, and a checkpoint resumes the uninterrupted run bit for bit.
/// </summary>
public sealed class SmcTests
{
    [Fact]
    public void Weights_NormalizeAtExtremeScales()
    {
        var w = new double[3];
        double logSum = SmcWeights.Normalize([1000.0, 0.0, -1000.0], w);
        Assert.Equal(1000.0, logSum, 12);
        Assert.Equal(1.0, w[0], 15);
        Assert.Equal(1.0, SmcWeights.EffectiveSampleSize(w), 12);

        var uniform = new double[8];
        Array.Fill(uniform, 1.0 / 8);
        Assert.Equal(8.0, SmcWeights.EffectiveSampleSize(uniform), 12);
        double[] increments = [0.3, -1.2, 5.0, 0.0, 2.2, -7.0, 1.0, 0.5];
        double cess = SmcWeights.ConditionalEss(uniform, increments);
        var shifted = new double[8];
        for (int i = 0; i < 8; i++) shifted[i] = increments[i] + 700.0;
        Assert.Equal(cess, SmcWeights.ConditionalEss(uniform, shifted), 10);
        Assert.Equal(8.0, SmcWeights.ConditionalEss(uniform, new double[8]), 12);
        Assert.Throws<InvalidOperationException>(() => SmcWeights.Normalize([double.NegativeInfinity, double.NegativeInfinity], new double[2]));
    }

    [Theory]
    [InlineData(ResamplingScheme.Multinomial)]
    [InlineData(ResamplingScheme.Stratified)]
    [InlineData(ResamplingScheme.Systematic)]
    [InlineData(ResamplingScheme.Residual)]
    public void Resampling_OffspringAverageNTimesTheWeights(ResamplingScheme scheme)
    {
        const int n = 40, reps = 20_000;
        var rng = new Xoshiro256PlusPlus(seed: 1201);
        var weights = new double[n];
        double total = 0.0;
        for (int i = 0; i < n; i++)
        {
            weights[i] = i % 7 == 0 ? 0.0 : rng.NextDouble() * rng.NextDouble();
            total += weights[i];
        }
        for (int i = 0; i < n; i++) weights[i] /= total;

        var sums = new double[n];
        var ancestors = new int[n];
        var counts = new int[n];
        for (int r = 0; r < reps; r++)
        {
            Resampling.Ancestors(scheme, weights, rng, ancestors);
            Array.Clear(counts);
            for (int k = 0; k < n; k++)
            {
                if (k > 0) Assert.True(ancestors[k] >= ancestors[k - 1]);
                Assert.True(weights[ancestors[k]] > 0.0);
                counts[ancestors[k]]++;
            }
            for (int i = 0; i < n; i++)
            {
                double expected = n * weights[i];
                if (scheme == ResamplingScheme.Systematic)
                    Assert.True(counts[i] >= Math.Floor(expected) && counts[i] <= Math.Ceiling(expected), $"particle {i}: {counts[i]} for {expected:F3}");
                if (scheme == ResamplingScheme.Residual) Assert.True(counts[i] >= Math.Floor(expected));
                sums[i] += counts[i];
            }
        }
        for (int i = 0; i < n; i++)
        {
            double expected = n * weights[i];
            double sd = Math.Sqrt(expected * (1.0 - weights[i]) / reps);
            Assert.True(Math.Abs(sums[i] / reps - expected) <= 5.0 * sd + 1e-12, $"particle {i}: mean {sums[i] / reps:F4} for {expected:F4}");
        }
    }

    [Fact]
    public void ConjugateGaussian_RecoversPosteriorAndEvidence()
    {
        // x ~ N(0, I₃), y_j | x ~ N(x_j, 0.5²): the posterior is N(0.8·y, 0.2·I) and the evidence Π N(y_j | 0, 1.25).
        (ISamplingTarget<double[]> target, double[] y) = Conjugate();
        var options = new SmcOptions { Particles = 2000, MutationSteps = 5 };
        SmcResult<double[]> result = SmcSampler.RunTempered(target, TemperingPath.Likelihood, 0.95, new GaussianParticles(3, 1.0), new RandomWalkKernel(), options, seed: 1202);

        double logZ = 0.0;
        foreach (double yj in y) logZ += -0.5 * Math.Log(2.0 * Math.PI * 1.25) - yj * yj / (2.0 * 1.25);
        Assert.True(Math.Abs(result.LogNormalizerRatio - logZ) <= 0.15, $"log Z {result.LogNormalizerRatio:F4} vs {logZ:F4}");
        for (int j = 0; j < 3; j++)
        {
            (double mean, double variance) = Moments(result, j);
            Assert.True(Math.Abs(mean - 0.8 * y[j]) <= 0.06, $"mean[{j}] {mean:F4} vs {0.8 * y[j]:F4}");
            Assert.True(Math.Abs(variance - 0.2) <= 0.04, $"var[{j}] {variance:F4}");
        }
        IReadOnlyList<SmcStageReport> stages = result.Stages;
        Assert.Equal(1.0, stages[^1].Temperature.InverseTemperature);
        for (int s = 1; s < stages.Count - 1; s++)
            Assert.True(Math.Abs(stages[s].ConditionalEss - 0.95 * 2000) <= 1.0, $"stage {s}: CESS {stages[s].ConditionalEss:F2}");
        Assert.True(stages[^1].ConditionalEss >= 0.95 * 2000 - 1.0);
        Assert.Equal(2000, stages[0].Families);
        Assert.True(stages[^1].Families >= 1 && stages[^1].Families <= 2000);
    }

    [Fact]
    public void BimodalLikelihood_RecoversTheModesMasses()
    {
        // x ~ N(0, 3²), likelihood 0.3·N(x | −2, 0.3²) + 0.7·N(x | 2, 0.3²): the modes carry 0.3 and 0.7 of the posterior.
        static double Log(double x, double m) => -0.5 * Math.Log(2.0 * Math.PI * 0.09) - (x - m) * (x - m) / 0.18;
        ISamplingTarget<double[]> target = SamplingTarget.Separable<double[]>(
            x => -0.5 * Math.Log(2.0 * Math.PI * 9.0) - x[0] * x[0] / 18.0,
            x =>
            {
                double a = Math.Log(0.3) + Log(x[0], -2.0), b = Math.Log(0.7) + Log(x[0], 2.0);
                double max = Math.Max(a, b);
                return max + Math.Log(Math.Exp(a - max) + Math.Exp(b - max));
            });
        var options = new SmcOptions { Particles = 4000, MutationSteps = 5 };
        SmcResult<double[]> result = SmcSampler.RunTempered(target, TemperingPath.Likelihood, 0.9, new GaussianParticles(1, 3.0), new RandomWalkKernel(), options, seed: 1203);
        double right = 0.0;
        for (int i = 0; i < result.Particles.Count; i++) if (result.Particles[i][0] > 0.0) right += result.Weights[i];
        Assert.True(Math.Abs(right - 0.7) <= 0.04, $"mass of the right mode {right:F4}");
    }

    [Fact]
    public void WithoutResamplingOrMoves_EvidenceTelescopesToImportanceSampling()
    {
        (ISamplingTarget<double[]> target, _) = Conjugate();
        var stages = new List<SmcStage<double[]>>();
        foreach (double b in new[] { 0.0, 0.1, 0.3, 0.6, 1.0 }) stages.Add(new SmcStage<double[]>(target, new Temperature(b, TemperingPath.Likelihood)));
        var options = new SmcOptions { Particles = 500, MutationSteps = 0, ResampleBelow = 0.0 };
        SmcResult<double[]> result = SmcSampler.Run(stages, new GaussianParticles(3, 1.0), new RandomWalkKernel(), options, seed: 1204);

        // The same prior draws: log mean_i exp(ℓ(x_i)).
        var logLik = new double[result.Particles.Count];
        double max = double.NegativeInfinity;
        for (int i = 0; i < logLik.Length; i++)
        {
            logLik[i] = target.Evaluate(result.Particles[i]).LogLikelihood(1.0);
            max = Math.Max(max, logLik[i]);
        }
        double sum = 0.0;
        foreach (double l in logLik) sum += Math.Exp(l - max);
        double importance = max + Math.Log(sum / logLik.Length);
        Assert.True(Math.Abs(result.LogNormalizerRatio - importance) <= 1e-12 * Math.Abs(importance), $"{result.LogNormalizerRatio:R} vs {importance:R}");
    }

    [Fact]
    public void Results_DoNotDependOnParallelism()
    {
        (ISamplingTarget<double[]> target, _) = Conjugate();
        SmcResult<double[]> Run(int threads) => SmcSampler.RunTempered(target, TemperingPath.Likelihood, 0.9, new GaussianParticles(3, 1.0),
            new RandomWalkKernel(), new SmcOptions { Particles = 300, MutationSteps = 3, MaxDegreeOfParallelism = threads }, seed: 1205);
        AssertSame(Run(1), Run(4));
    }

    [Fact]
    public void Checkpoint_ResumesTheUninterruptedRun()
    {
        (ISamplingTarget<double[]> target, _) = Conjugate();
        var particles = new GaussianParticles(3, 1.0);
        var options = new SmcOptions { Particles = 300, MutationSteps = 2 };
        var schedule = new List<SmcStage<double[]>>();
        for (int s = 0; s <= 10; s++) schedule.Add(new SmcStage<double[]>(target, new Temperature(s / 10.0, TemperingPath.Likelihood)));

        SmcRun<double[]> run = SmcRun<double[]>.Start(schedule[0], particles, new RandomWalkKernel(), options, seed: 1206);
        SmcCheckpoint<double[]>? checkpoint = null;
        for (int s = 1; s <= 10; s++)
        {
            run.Advance(schedule[s]);
            if (s == 4) checkpoint = run.Checkpoint();
        }
        SmcRun<double[]> resumed = SmcRun<double[]>.Resume(checkpoint!, particles, new RandomWalkKernel(), options);
        for (int s = 5; s <= 10; s++) resumed.Advance(schedule[s]);
        AssertSame(run.Result(), resumed.Result());
    }

    private static (ISamplingTarget<double[]> Target, double[] Y) Conjugate()
    {
        double[] y = [1.0, -0.5, 2.0];
        const double s2 = 0.25;
        ISamplingTarget<double[]> target = SamplingTarget.Separable<double[]>(
            x =>
            {
                double v = -1.5 * Math.Log(2.0 * Math.PI);
                foreach (double xj in x) v -= 0.5 * xj * xj;
                return v;
            },
            x =>
            {
                double v = -1.5 * Math.Log(2.0 * Math.PI * s2);
                for (int j = 0; j < 3; j++) v -= (y[j] - x[j]) * (y[j] - x[j]) / (2.0 * s2);
                return v;
            });
        return (target, y);
    }

    private static (double Mean, double Variance) Moments(SmcResult<double[]> result, int j)
    {
        double mean = 0.0, second = 0.0;
        for (int i = 0; i < result.Particles.Count; i++)
        {
            mean += result.Weights[i] * result.Particles[i][j];
            second += result.Weights[i] * result.Particles[i][j] * result.Particles[i][j];
        }
        return (mean, second - mean * mean);
    }

    private static void AssertSame(SmcResult<double[]> a, SmcResult<double[]> b)
    {
        Assert.Equal(a.LogNormalizerRatio, b.LogNormalizerRatio);
        Assert.Equal(a.Weights, b.Weights);
        Assert.Equal(a.Particles.Count, b.Particles.Count);
        for (int i = 0; i < a.Particles.Count; i++) Assert.Equal(a.Particles[i], b.Particles[i]);
        Assert.Equal(a.Stages.Count, b.Stages.Count);
        for (int s = 0; s < a.Stages.Count; s++) Assert.Equal(a.Stages[s], b.Stages[s]);
    }

    // Particles drawn from N(0, sd²·I_d).
    private sealed class GaussianParticles(int dimension, double sd) : ISmcParticles<double[]>
    {
        public double[] Draw(Xoshiro256PlusPlus rng)
        {
            var x = new double[dimension];
            for (int j = 0; j < dimension; j++) x[j] = sd * Normal.SampleStandard(rng);
            return x;
        }

        public double[] Fork(double[] state) => (double[])state.Clone();
    }
}
