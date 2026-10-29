using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The conditional sampler for θ given one emulator draw (stage 7; D09): the whole chain's laws (θ's walk, its
/// adaptation, R's decorrelation step, the levels, the variances on ln σ² and the discrepancy's draw) against quadrature
/// on a small fixture, untempered and at every rung of a likelihood-path ladder ending at b = 0, with v integrated out
/// in closed form; the swap's ratio, (b₁ − b₂)(ℓ(x₂) − ℓ(x₁)); the propriety rule; and the outputs in physical
/// coordinates, with the burn-in, thinning and emulator draw explicit.
/// </summary>
public sealed class BassCalibrationSamplerTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationSamplerTests(ITestOutputHelper output) => _out = output;

    private const int M = CalibrationFixture.M;

    // One continuous θ, one categorical (3 levels) or none, one inverse-gamma group, the functional discrepancy.
    internal static CalibrationModel Model(int categoricals, CalibrationSampling sampling)
    {
        TestEmulator emulator = CalibrationFixture.Emulator(categoricals);
        CalibrationConfig config = CalibrationFixture.Config(emulator, CalibrationFixture.Groups(one: true)) with { Sampling = sampling };
        return new CalibrationModel(config, CalibrationFixture.Observations(emulator), emulator);
    }

    /// <summary>
    /// Quadrature of the rung-b target over (θ, ln σ², level) with v integrated out: p_b ∝ π(σ²)·|Σ|^((1−b)/2)·
    /// N(y; prediction + μ_δ, Σ/b + s_v²DDᵀ), its marginal moments of θ's unit coordinate, ln σ² and each level; at b = 0
    /// the prior's.
    /// </summary>
    internal static (double MeanTheta, double MedianTheta, double MeanLogVariance, double[] Levels) Quadrature(CalibrationModel model, CalibrationTarget target, double b)
    {
        int levels = model.Layout.Categorical.Count == 0 ? 1 : model.Layout.Categorical[0].Levels.Length;
        const int thetaPoints = 600, uPoints = 260;
        const double uLower = -14, uUpper = 0;
        double hu = (uUpper - uLower) / uPoints;
        double[,] d = model.Discrepancy!;
        double s2 = model.DiscrepancyScale * model.DiscrepancyScale;
        var noise = (ComponentNoise)target.Noise;
        double[,] residuals = model.Emulator.TruncationResiduals()!;
        int n = residuals.GetLength(0);
        var fixedPart = new double[M, M];   // the components and the residual set, which σ² does not enter
        for (int j = 0; j < M; j++)
            for (int l = 0; l < M; l++)
            {
                double s = 0;
                for (int k = 0; k < noise.Variances.Length; k++) s += noise.Factor[j, k] * noise.Variances[k] * noise.Factor[l, k];
                for (int i = 0; i < n; i++) s += residuals[i, j] * residuals[i, l] / n;
                fixedPart[j, l] = s;
            }
        // Per σ²: (1−b)/2·ln|Σ| and the factor of Σ/b + s²DDᵀ.
        var weightOfU = new double[uPoints];
        var factors = new DenseCholesky[uPoints];
        for (int iu = 0; iu < uPoints; iu++)
        {
            double u = uLower + (iu + 0.5) * hu, sigma2 = Math.Exp(u);
            weightOfU[iu] = model.Measurement.LogPrior(0, sigma2) + u;
            if (b == 0) continue;
            var sigma = new double[M, M];
            var marginal = new double[M, M];
            for (int j = 0; j < M; j++)
                for (int l = 0; l < M; l++)
                {
                    sigma[j, l] = fixedPart[j, l] + (j == l ? sigma2 : 0);
                    marginal[j, l] = sigma[j, l] / b + s2 * (d[j, 0] * d[l, 0] + d[j, 1] * d[l, 1]);
                }
            var fs = new DenseCholesky(M);
            fs.Decompose(sigma);
            weightOfU[iu] += 0.5 * (1 - b) * fs.LogDet;
            factors[iu] = new DenseCholesky(M);
            factors[iu].Decompose(marginal);
        }
        var logW = new double[levels, thetaPoints, uPoints];
        var pred = new double[M];
        double max = double.NegativeInfinity;
        for (int c = 0; c < levels; c++)
            for (int it = 0; it < thetaPoints; it++)
            {
                double unit = (it + 0.5) / thetaPoints;
                model.Emulator.Predict(target.Draw, [CalibrationFixture.Theta.ToPhysical(unit)], levels > 1 ? [c] : [], pred);
                double[] r = [.. Enumerable.Range(0, M).Select(j => model.Observations[j] - pred[j] - 0.01 * j)];
                for (int iu = 0; iu < uPoints; iu++)
                {
                    double lw = weightOfU[iu];
                    if (b > 0) lw += -0.5 * (factors[iu].LogDet + factors[iu].InverseQuadraticForm(r));
                    logW[c, it, iu] = lw;
                    max = Math.Max(max, lw);
                }
            }
        double total = 0, meanTheta = 0, meanU = 0;
        var thetaMass = new double[thetaPoints];
        var levelMass = new double[levels];
        for (int c = 0; c < levels; c++)
            for (int it = 0; it < thetaPoints; it++)
                for (int iu = 0; iu < uPoints; iu++)
                {
                    double w = Math.Exp(logW[c, it, iu] - max);
                    total += w;
                    meanTheta += w * (it + 0.5) / thetaPoints;
                    meanU += w * (uLower + (iu + 0.5) * hu);
                    thetaMass[it] += w;
                    levelMass[c] += w;
                }
        double cumulative = 0, median = 1;
        for (int it = 0; it < thetaPoints; it++)
            if ((cumulative += thetaMass[it] / total) >= 0.5) { median = (it + 0.5) / thetaPoints; break; }
        return (meanTheta / total, median, meanU / total, [.. levelMass.Select(v => v / total)]);
    }

    private void Compare(string label, (double MeanTheta, double MedianTheta, double MeanLogVariance, double[] Levels) q, IReadOnlyList<CalibrationState> draws)
    {
        double[] theta = [.. draws.Select(s => s.Unit[0])];
        double[] below = [.. theta.Select(t => t <= q.MedianTheta ? 1.0 : 0.0)];
        double[] u = [.. draws.Select(s => Math.Log(s.Variances[0]))];
        _out.WriteLine($"{label}: E[θ] {theta.Average():F4} vs {q.MeanTheta:F4} (MCSE {CalibrationFixture.Mcse(theta):F4}); P(θ ≤ median) {below.Average():F3} " +
                       $"(MCSE {CalibrationFixture.Mcse(below):F3}); E[ln σ²] {u.Average():F3} vs {q.MeanLogVariance:F3} (MCSE {CalibrationFixture.Mcse(u):F3})");
        const double grid = 2.0 / 600;   // the θ grid's resolution
        Assert.InRange(theta.Average() - q.MeanTheta, -5 * CalibrationFixture.Mcse(theta) - grid, 5 * CalibrationFixture.Mcse(theta) + grid);
        Assert.InRange(below.Average() - 0.5, -5 * CalibrationFixture.Mcse(below) - 0.01, 5 * CalibrationFixture.Mcse(below) + 0.01);
        Assert.InRange(u.Average() - q.MeanLogVariance, -5 * CalibrationFixture.Mcse(u) - 0.02, 5 * CalibrationFixture.Mcse(u) + 0.02);
        for (int c = 0; c < q.Levels.Length && draws[0].Levels.Length > 0; c++)
        {
            double[] at = [.. draws.Select(s => s.Levels[0] == c ? 1.0 : 0.0)];
            Assert.InRange(at.Average() - q.Levels[c], -5 * CalibrationFixture.Mcse(at) - 0.005, 5 * CalibrationFixture.Mcse(at) + 0.005);
        }
    }

    [Fact]
    public void Untempered_TheChainsLawIsTheConditional()
    {
        var sampling = new CalibrationSampling { MasterSeed = 3, BurnIn = 3000, Samples = 40000, DecorrelationEvery = 25 };
        CalibrationModel model = Model(categoricals: 1, sampling);
        ConditionalRun run = CalibrationSampler.Run(model, 1, sampling);
        Assert.Equal(1, run.EmulatorDraw);
        Assert.Equal(40000, run.Draws.Count);
        Assert.True(double.IsNaN(run.SwapRate));
        CalibrationState[] states = [.. run.Draws.Select(d => new CalibrationState([CalibrationFixture.Theta.ToUnit(d.Theta[0])], d.Levels, d.Variances, d.Coefficients))];
        Compare("untempered", Quadrature(model, model.Target(1), 1.0), states);
        foreach (MoveStat stat in run.Statistics) _out.WriteLine($"  {stat.Key}: {stat.Accepted}/{stat.Attempts}");
        Assert.Equal(["theta", "theta.adapt", "theta.decorrelate", "levels (enumerated)", "variances", "discrepancy"], run.Statistics.Select(s => s.Key));
    }

    [Fact]
    public void Tempered_EveryRungsLawIsItsLikelihoodPathTarget()
    {
        double[] ladder = [1.0, 0.6, 0.25, 0.0];
        var sampling = new CalibrationSampling { MasterSeed = 5, BurnIn = 3000, Samples = 30000, DecorrelationEvery = 25, Tempering = ladder };
        CalibrationModel model = Model(categoricals: 0, sampling);
        CalibrationTarget target = model.Target(0);
        var chain = new CalibrationChain(target, sampling, model.Start(new Xoshiro256PlusPlus(9)), 11);
        chain.Burn(sampling.BurnIn);
        var draws = new List<CalibrationState>[ladder.Length];
        for (int i = 0; i < ladder.Length; i++) draws[i] = new List<CalibrationState>(sampling.Samples);
        for (int s = 0; s < sampling.Samples; s++)
        {
            chain.Step();
            for (int i = 0; i < ladder.Length; i++) draws[i].Add(chain.Ladder!.Replicas[i].Current);
        }
        _out.WriteLine($"swap acceptance {chain.Ladder!.SwapAccepts / (double)chain.Ladder.SwapAttempts:F2}, round trips {chain.Ladder.RoundTrips}");
        for (int i = 0; i < ladder.Length; i++) Compare($"b = {ladder[i]}", Quadrature(model, target, ladder[i]), draws[i]);
        Assert.True(chain.Ladder.RoundTrips > 10);
    }

    [Fact]
    public void Multivariate_TheChainsLawIsTheConditional()
    {
        // The multivariate formulation: no discrepancy to draw; quadrature over (θ, ln σ²) of π(σ²)·N(y; prediction + μ_δ, Σ).
        double[][] correlation = [.. Enumerable.Range(0, M).Select(i => Enumerable.Range(0, M).Select(j => Math.Pow(0.4, Math.Abs(i - j))).ToArray())];
        double[][] delta = [.. Enumerable.Range(0, M).Select(i => Enumerable.Range(0, M).Select(j => (i == j ? 1e-4 : 0.0) + 0.002 * Math.Exp(-0.5 * (i - j) * (i - j))).ToArray())];
        double[] mean = [.. Enumerable.Range(0, M).Select(j => 0.01 * j)];
        var sampling = new CalibrationSampling { MasterSeed = 21, BurnIn = 3000, Samples = 40000, DecorrelationEvery = 25 };
        TestEmulator emulator = CalibrationFixture.Emulator();
        CalibrationConfig config = CalibrationFixture.Config(emulator, CalibrationFixture.Groups(one: true), new MultivariateFormulation(correlation, delta, mean))
                                   with { Sampling = sampling };
        var model = new CalibrationModel(config, CalibrationFixture.Observations(emulator), emulator);
        CalibrationTarget target = model.Target(1);
        ConditionalRun run = CalibrationSampler.Run(model, 1, sampling);
        Assert.DoesNotContain(run.Statistics, s => s.Key == "discrepancy");

        var noise = (ComponentNoise)target.Noise;
        double[,] residuals = emulator.TruncationResiduals()!;
        int n = residuals.GetLength(0);
        const int thetaPoints = 600, uPoints = 260;
        const double uLower = -14, uUpper = 0, hu = (uUpper - uLower) / uPoints;
        var logW = new double[thetaPoints, uPoints];
        var factors = new DenseCholesky[uPoints];
        for (int iu = 0; iu < uPoints; iu++)
        {
            double s2 = Math.Exp(uLower + (iu + 0.5) * hu);
            var sigma = new double[M, M];
            for (int j = 0; j < M; j++)
                for (int l = 0; l < M; l++)
                {
                    double s = s2 * correlation[j][l] + delta[j][l];
                    for (int k = 0; k < 2; k++) s += noise.Factor[j, k] * noise.Variances[k] * noise.Factor[l, k];
                    for (int i = 0; i < n; i++) s += residuals[i, j] * residuals[i, l] / n;
                    sigma[j, l] = s;
                }
            factors[iu] = new DenseCholesky(M);
            factors[iu].Decompose(sigma);
        }
        var pred = new double[M];
        double max = double.NegativeInfinity;
        for (int it = 0; it < thetaPoints; it++)
        {
            emulator.Predict(1, [CalibrationFixture.Theta.ToPhysical((it + 0.5) / thetaPoints)], [], pred);
            double[] r = [.. Enumerable.Range(0, M).Select(j => model.Observations[j] - pred[j] - mean[j])];
            for (int iu = 0; iu < uPoints; iu++)
            {
                double u = uLower + (iu + 0.5) * hu;
                logW[it, iu] = model.Measurement.LogPrior(0, Math.Exp(u)) + u - 0.5 * (factors[iu].LogDet + factors[iu].InverseQuadraticForm(r));
                max = Math.Max(max, logW[it, iu]);
            }
        }
        double total = 0, meanTheta = 0, meanU = 0;
        for (int it = 0; it < thetaPoints; it++)
            for (int iu = 0; iu < uPoints; iu++)
            {
                double w = Math.Exp(logW[it, iu] - max);
                total += w;
                meanTheta += w * (it + 0.5) / thetaPoints;
                meanU += w * (uLower + (iu + 0.5) * hu);
            }
        double[] theta = [.. run.Draws.Select(d => CalibrationFixture.Theta.ToUnit(d.Theta[0]))];
        double[] logVariance = [.. run.Draws.Select(d => Math.Log(d.Variances[0]))];
        _out.WriteLine($"multivariate: E[θ] {theta.Average():F4} vs {meanTheta / total:F4} (MCSE {CalibrationFixture.Mcse(theta):F4}); " +
                       $"E[ln σ²] {logVariance.Average():F3} vs {meanU / total:F3} (MCSE {CalibrationFixture.Mcse(logVariance):F3})");
        Assert.InRange(theta.Average() - meanTheta / total, -5 * CalibrationFixture.Mcse(theta) - 2.0 / thetaPoints, 5 * CalibrationFixture.Mcse(theta) + 2.0 / thetaPoints);
        Assert.InRange(logVariance.Average() - meanU / total, -5 * CalibrationFixture.Mcse(logVariance) - 0.02, 5 * CalibrationFixture.Mcse(logVariance) + 0.02);
    }

    [Fact]
    public void Swap_IsTheLikelihoodPathsRatio()
    {
        CalibrationModel model = Model(0, new CalibrationSampling());
        CalibrationTarget target = model.Target(2);
        var rng = new Xoshiro256PlusPlus(31);
        for (int trial = 0; trial < 20; trial++)
        {
            CalibrationState x1 = model.Start(rng).WithCoefficients([rng.NextDouble() - 0.5, rng.NextDouble() - 0.5]);
            CalibrationState x2 = model.Start(rng);
            CalibrationEvaluation e1 = target.EvaluateState(x1), e2 = target.EvaluateState(x2);
            double b1 = 0.3 + 0.7 * rng.NextDouble(), b2 = b1 * rng.NextDouble();
            Temperature t1 = new(b1, TemperingPath.Likelihood), t2 = new(b2, TemperingPath.Likelihood);
            // ParallelTempering's ratio, d₁(x₂) + d₂(x₁) − d₁(x₁) − d₂(x₂): the priors and v's prior cancel, so it is the
            // likelihood's alone, with no Jacobian (R's swap omits them from a target that tempers them).
            double ratio = t1.LogDensity(e2) + t2.LogDensity(e1) - t1.LogDensity(e1) - t2.LogDensity(e2);
            Assert.Equal((b1 - b2) * (e2.LogLikelihoodValue - e1.LogLikelihoodValue), ratio, 1e-9);
        }
    }

    [Fact]
    public void ProprietyRule_ALadderToZeroNeedsProperPriors()
    {
        TestEmulator emulator = CalibrationFixture.Emulator();
        var sampling = new CalibrationSampling { Tempering = [1.0, 0.5, 0.0] };
        CalibrationConfig config = CalibrationFixture.Config(emulator) with { Sampling = sampling };
        double[] y = CalibrationFixture.Observations(emulator);
        var start = new CalibrationState([0.5], [], [0.01, 0.01], [0.0, 0.0]);
        _ = new CalibrationChain(new CalibrationModel(config, y, emulator).Target(0), sampling, start, 1);

        var improper = new CalibrationModel(config, y, emulator,
            new MeasurementModel([0, 0, 0, 0, 1, 1, 1, 1], [new RsKernelPrior(), new InverseGammaVariancePrior(8, 0.05)]));
        Assert.Contains("improper", Assert.Throws<ArgumentException>(() => new CalibrationChain(improper.Target(0), sampling, start, 1)).Message);
        _ = new CalibrationChain(improper.Target(0), sampling with { Tempering = [1.0, 0.5, 0.1] }, start, 1);
        Assert.Contains("descend", Assert.Throws<ArgumentException>(() => new CalibrationSampling { Tempering = [1.0, 0.5, 0.5] }.Validate()).Message);
        Assert.Contains("cold", Assert.Throws<ArgumentException>(() => new CalibrationSampling { Tempering = [0.9, 0.5] }.Validate()).Message);
    }

    private sealed class RsKernelPrior : IVariancePrior
    {
        public double LogDensity(double sigma2) => sigma2 > 0 ? -Math.Log(1 + sigma2) : double.NegativeInfinity;
        public double Sample(Xoshiro256PlusPlus rng) => throw new NotSupportedException();
        public bool IsProper => false;
    }

    [Fact]
    public void Outputs_ArePhysical_Reproducible_AndThinned()
    {
        var sampling = new CalibrationSampling { MasterSeed = 13, BurnIn = 200, Samples = 50, Thin = 3 };
        CalibrationModel model = Model(1, sampling);
        ConditionalRun a = CalibrationSampler.Run(model, 2, sampling), b = CalibrationSampler.Run(model, 2, sampling);
        Assert.Equal(3, a.Thin);
        Assert.Equal(200, a.BurnIn);
        Assert.Equal(50, a.Draws.Count);
        for (int i = 0; i < a.Draws.Count; i++)
        {
            CalibrationDraw d = a.Draws[i];
            Assert.Equal(2, d.EmulatorDraw);
            Assert.InRange(d.Theta[0], CalibrationFixture.Theta.Lower, CalibrationFixture.Theta.Upper);
            Assert.Equal(d.Theta, b.Draws[i].Theta);
            Assert.Equal(d.Variances, b.Draws[i].Variances);
            var state = new CalibrationState([CalibrationFixture.Theta.ToUnit(d.Theta[0])], d.Levels, d.Variances, d.Coefficients);
            Assert.Equal(model.Target(2).EvaluateState(state).LogLikelihoodValue, d.LogLikelihood, 1e-9);
        }
        Assert.NotEqual(a.Draws[0].Theta, CalibrationSampler.Run(model, 2, sampling with { MasterSeed = 14 }).Draws[0].Theta);
        Assert.Contains("Thin", Assert.Throws<ArgumentException>(() => CalibrationSampler.Run(model, 0, sampling with { Thin = 0 })).Message);
    }
}
