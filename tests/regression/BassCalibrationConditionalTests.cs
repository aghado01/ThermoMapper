using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// A small calibration problem for the conditional and sampler checks: m = 8 outputs, θ on [2, 5], optionally one or two
/// categorical parameters (3 and 4 levels), two variance groups (an inverse gamma and a half-Cauchy) or one, the
/// functional discrepancy on two columns with s_v = 0.3 and a mean, and five truncation residual vectors.
/// </summary>
internal static class CalibrationFixture
{
    public const int M = 8;

    public static readonly ContinuousParameter Theta = new("theta", 2, 5);
    public static readonly CategoricalParameter Mode = new("mode", ["a", "b", "c"]);
    public static readonly CategoricalParameter Phase = new("phase", ["p", "q", "r", "s"]);

    public static TestEmulator Emulator(int categoricals = 0, bool residuals = true, bool independent = false)
    {
        var offset = new double[M];
        var scale = new double[M];
        var basis = new double[M, 2];
        var rng = new Xoshiro256PlusPlus(41);
        var res = new double[5, M];
        for (int l = 0; l < M; l++)
        {
            double t = l / (M - 1.0);
            offset[l] = 1 + 0.2 * l;
            scale[l] = 0.8 + 0.1 * l;
            basis[l, 0] = Math.Sin(Math.PI * t) + 0.2;
            basis[l, 1] = Math.Cos(Math.PI * t);
            for (int i = 0; i < 5; i++) res[i, l] = scale[l] * 0.04 * Normal.SampleStandard(rng);
        }
        List<CalibrationParameter> parameters = [Theta];
        if (categoricals >= 1) parameters.Add(Mode);
        if (categoricals >= 2) parameters.Add(Phase);
        double[][] variances = independent ? [[0.01], [0.02], [0.015]] : [[0.02, 0.05], [0.03, 0.04], [0.01, 0.06]];
        return new TestEmulator(parameters, offset, scale, basis, (draw, theta, levels) =>
        {
            double l0 = levels.Length > 0 ? levels[0] : 0, l1 = levels.Length > 1 ? levels[1] : 0;
            return [Math.Sin(theta[0]) + 0.4 * l0 - 0.3 * l1 + 0.05 * draw, Math.Cos(0.7 * theta[0]) - 0.25 * l0 + 0.2 * l1];
        }, variances, residuals ? res : null, independent);
    }

    public static GaussianMeasurement Groups(bool one = false, VariancePrior? only = null)
        => one
            ? new GaussianMeasurement([new VarianceGroup([.. Enumerable.Range(0, M)], only ?? new InverseGammaVariance(6, 0.08))])
            : new GaussianMeasurement([new VarianceGroup([0, 1, 2, 3], new InverseGammaVariance(8, 0.05)), new VarianceGroup([4, 5, 6, 7], new HalfCauchyDeviation(0.05))]);

    public static FunctionalDiscrepancy Discrepancy()
        => new([.. Enumerable.Range(0, M).Select(j => new[] { 1.0, j / (M - 1.0) - 0.5 })], 0.3, [.. Enumerable.Range(0, M).Select(j => 0.01 * j)]);

    public static CalibrationConfig Config(TestEmulator emulator, GaussianMeasurement? groups = null, CovarianceFormulation? covariance = null, int cap = 1024)
        => new()
        {
            Parameters = emulator.Parameters,
            Measurement = groups ?? Groups(),
            Covariance = covariance ?? new FunctionalFormulation { Discrepancy = Discrepancy() },
            EnumerationCap = cap,
        };

    // Observations at θ = 3.4 and levels (1, 2), with the outputs' noise.
    public static double[] Observations(TestEmulator emulator, int seed = 7)
    {
        var y = new double[M];
        int[] levels = [.. new[] { 1, 2 }.Take(emulator.Parameters.Count - 1)];
        emulator.Predict(0, [3.4], levels, y);
        var rng = new Xoshiro256PlusPlus(seed);
        for (int j = 0; j < M; j++) y[j] += 0.05 * Normal.SampleStandard(rng) + 0.02 * Math.Sin(j);
        return y;
    }

    public static CalibrationModel Model(int categoricals = 0, GaussianMeasurement? groups = null, CovarianceFormulation? covariance = null, int cap = 1024)
    {
        TestEmulator emulator = Emulator(categoricals);
        return new CalibrationModel(Config(emulator, groups, covariance, cap), Observations(emulator), emulator);
    }

    /// <summary>A move that is never available, so a chain's steps are its blocks'.</summary>
    public sealed class Stay : IRjMove<CalibrationState>
    {
        public string Key => "stay";
        public string ReverseKey => "stay";
        public double Weight(CalibrationState state) => 1.0;
        public Proposal<CalibrationState>? Propose(CalibrationState current, Xoshiro256PlusPlus rng) => null;
    }

    public static ReversibleJumpChain<CalibrationState> Chain(CalibrationTarget target, CalibrationState start, double b, int seed = 1)
        => new([new Stay()], target, start, new Xoshiro256PlusPlus(seed), new Temperature(b, TemperingPath.Likelihood));

    // Batch-means Monte Carlo standard error of the mean of `x`.
    public static double Mcse(IReadOnlyList<double> x, int batches = 40)
    {
        int size = x.Count / batches;
        double mean = x.Take(size * batches).Average();
        double[] means = [.. Enumerable.Range(0, batches).Select(b => x.Skip(b * size).Take(size).Average())];
        return Math.Sqrt(means.Sum(v => (v - mean) * (v - mean)) / (batches - 1) / batches);
    }
}

/// <summary>
/// Calibration's conditional model and its conditionals (stage 7, C2; D09), each against an independent computation on
/// small fixtures: the target's log density against a dense Gaussian with every prior written out; the caches kept and
/// invalidated by what changed; the discrepancy's Gaussian conditional against the joint-Gaussian (Kalman) form at each
/// rung; the variance groups' moves on ln σ² against quadrature, untempered, tempered and at b = 0, under both priors;
/// the enumerated categorical conditional against brute force; the single-site updates above the cap against
/// enumeration, with their mixing measured; the propriety rule; and R's (1 + σ²)⁻¹ kernel as a defect fixture, improper
/// and unit-dependent where the half-Cauchy is neither.
/// </summary>
public sealed class BassCalibrationConditionalTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationConditionalTests(ITestOutputHelper output) => _out = output;

    private const int M = CalibrationFixture.M;

    private static double[,] DenseSigma(CalibrationModel model, CalibrationTarget target, CalibrationState state)
    {
        var sigma = new double[M, M];
        var perOutput = new double[M];
        model.Measurement.OutputVariances(state.Variances, perOutput);
        var noise = (ComponentNoise)target.Noise;
        double[,] residuals = model.Emulator.TruncationResiduals()!;
        int n = residuals.GetLength(0);
        for (int j = 0; j < M; j++)
            for (int l = 0; l < M; l++)
            {
                double s = j == l ? perOutput[j] : 0.0;
                for (int k = 0; k < noise.Variances.Length; k++) s += noise.Factor[j, k] * noise.Variances[k] * noise.Factor[l, k];
                for (int i = 0; i < n; i++) s += residuals[i, j] * residuals[i, l] / n;
                sigma[j, l] = s;
            }
        return sigma;
    }

    private static double DenseLogNormal(double[] r, double[,] sigma)
    {
        var f = new DenseCholesky(M);
        f.Decompose(sigma);
        return -0.5 * (M * Math.Log(2 * Math.PI) + f.LogDet + f.InverseQuadraticForm(r));
    }

    private static double[] Inverse(double[,] a, double[] b)
    {
        var f = new DenseCholesky(a.GetLength(0));
        f.Decompose(a);
        return f.Solve(b);
    }

    [Fact]
    public void Target_IsTheDenseGaussianWithEveryPriorWrittenOut()
    {
        CalibrationModel model = CalibrationFixture.Model(categoricals: 2);
        CalibrationTarget target = model.Target(1);
        var state = new CalibrationState([0.47], [2, 1], [0.004, 0.009], [0.2, -0.4]);
        CalibrationEvaluation e = target.EvaluateState(state);

        double theta = CalibrationFixture.Theta.ToPhysical(0.47);
        var pred = new double[M];
        model.Emulator.Predict(1, [theta], [2, 1], pred);
        double[,] d = model.Discrepancy!;
        double[] r = [.. Enumerable.Range(0, M).Select(j => model.Observations[j] - pred[j] - 0.01 * j - d[j, 0] * 0.2 - d[j, 1] * -0.4)];
        double likelihood = DenseLogNormal(r, DenseSigma(model, target, state));
        double a = 4, b = a * 0.05 * 0.05;   // IG(ν/2, (ν/2)ŝ²), ν = 8
        double ig = a * Math.Log(b) - SpecialFunctions.LogGamma(a) - (a + 1) * Math.Log(0.004) - b / 0.004;
        double hc = Math.Log(2 / (Math.PI * 0.05)) - Math.Log(1 + 0.009 / (0.05 * 0.05)) - Math.Log(2) - 0.5 * Math.Log(0.009);
        double v = -Math.Log(2 * Math.PI * 0.09) - (0.04 + 0.16) / (2 * 0.09);
        double prior = -Math.Log(3) - Math.Log(4) + ig + hc + v;
        Assert.Equal(likelihood, e.LogLikelihoodValue, 1e-10);
        Assert.Equal(prior, e.LogPrior, 1e-10);
        Assert.Equal(0.3 * likelihood, e.LogLikelihood(0.3), 1e-10);
        Assert.Equal(e.LogLikelihoodValue, e.LogLikelihoodSlope(0.3));
        Assert.Equal([theta], model.Physical(state));

        // Off the support: no likelihood work, so no prediction.
        var emulator = (TestEmulator)model.Emulator;
        int calls = emulator.Predictions;
        foreach (CalibrationState off in new[] { state.WithUnit([1.2]), state.WithVariances([0.0, 0.01]), state.WithLevels([3, 0]) })
        {
            CalibrationEvaluation o = target.EvaluateState(off);
            Assert.Equal(double.NegativeInfinity, o.LogPrior);
            Assert.Equal(double.NegativeInfinity, Temperature.Posterior.LogDensity(o));
        }
        Assert.Equal(calls, emulator.Predictions);
    }

    [Fact]
    public void Target_UnderTheMultivariateFormulation_AndWithoutABasis_IsTheDenseGaussian()
    {
        // Multivariate: AR(1) correlation, a fixed Σ_δ and mean, the residual set formed whole.
        double[][] correlation = [.. Enumerable.Range(0, M).Select(i => Enumerable.Range(0, M).Select(j => Math.Pow(0.5, Math.Abs(i - j))).ToArray())];
        double[][] delta = [.. Enumerable.Range(0, M).Select(i => Enumerable.Range(0, M).Select(j => (i == j ? 0.002 : 0.0) + 0.001 * Math.Exp(-Math.Abs(i - j))).ToArray())];
        double[] mean = [.. Enumerable.Range(0, M).Select(j => -0.02 * j)];
        CalibrationModel model = CalibrationFixture.Model(covariance: new MultivariateFormulation(correlation, delta, mean));
        Assert.Equal(0, model.Coefficients);
        CalibrationTarget target = model.Target(1);
        var state = new CalibrationState([0.3], [], [0.004, 0.009], []);
        double[,] sigma = DenseSigma(model, target, state);
        var perOutput = new double[M];
        model.Measurement.OutputVariances(state.Variances, perOutput);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < M; j++)
            {
                sigma[i, j] += Math.Sqrt(perOutput[i]) * correlation[i][j] * Math.Sqrt(perOutput[j]) + delta[i][j];
                if (i == j) sigma[i, j] -= perOutput[i];   // DenseSigma put the independent variances on the diagonal
            }
        var pred = new double[M];
        model.Emulator.Predict(1, model.Physical(state), [], pred);
        double expected = DenseLogNormal([.. Enumerable.Range(0, M).Select(j => model.Observations[j] - pred[j] - mean[j])], sigma);
        Assert.Equal(expected, target.EvaluateState(state).LogLikelihoodValue, 1e-10);

        // An emulator without a basis or residuals: σ²_e·I on the diagonal, with the functional discrepancy.
        TestEmulator independent = CalibrationFixture.Emulator(residuals: false, independent: true);
        var direct = new CalibrationModel(CalibrationFixture.Config(independent), CalibrationFixture.Observations(independent), independent);
        var s = new CalibrationState([0.3], [], [0.004, 0.009], [0.2, 0.1]);
        independent.Predict(2, direct.Physical(s), [], pred);
        double[] d = direct.Delta(s);
        var diagonal = new double[M, M];
        direct.Measurement.OutputVariances(s.Variances, perOutput);
        for (int j = 0; j < M; j++) diagonal[j, j] = perOutput[j] + 0.015;
        Assert.Equal(DenseLogNormal([.. Enumerable.Range(0, M).Select(j => direct.Observations[j] - pred[j] - d[j])], diagonal),
                     direct.Target(2).EvaluateState(s).LogLikelihoodValue, 1e-10);
    }

    [Fact]
    public void Config_WithEveryDeclaration_RoundTripsThroughJson()
    {
        TestEmulator emulator = CalibrationFixture.Emulator(categoricals: 2);
        double[][] identity = [.. Enumerable.Range(0, M).Select(i => Enumerable.Range(0, M).Select(j => i == j ? 1.0 : 0.0).ToArray())];
        foreach (CovarianceFormulation formulation in new CovarianceFormulation[]
                 {
                     new FunctionalFormulation { Discrepancy = CalibrationFixture.Discrepancy(), Truncation = TruncationTerm.Diagonal },
                     new FunctionalFormulation(),
                     new MultivariateFormulation(identity, identity, new double[M]),
                 })
        {
            CalibrationConfig config = CalibrationFixture.Config(emulator, covariance: formulation, cap: 64);
            string json = config.ToJson();
            CalibrationConfig back = CalibrationConfig.FromJson(json);
            Assert.Equal(json, back.ToJson());
            Assert.Equal(64, back.EnumerationCap);
            Assert.Equal(formulation.GetType(), back.Covariance.GetType());
            _ = new CalibrationModel(back, CalibrationFixture.Observations(emulator), emulator);
        }
        string missing = CalibrationFixture.Config(emulator).ToJson().Replace("\"Scale\":0.3", "\"Scale\":null");
        Assert.Contains("CalibrationConfig.Covariance.Discrepancy.Scale", Assert.Throws<ArgumentException>(() => CalibrationConfig.FromJson(missing)).Message);
        Assert.Contains("EnumerationCap", Assert.Throws<ArgumentException>(() => (CalibrationFixture.Config(emulator) with { EnumerationCap = 0 }).Validate()).Message);
        Assert.Contains("other parameters", Assert.Throws<ArgumentException>(() =>
            new CalibrationModel(CalibrationFixture.Config(CalibrationFixture.Emulator()), CalibrationFixture.Observations(emulator), emulator)).Message);
        Assert.Contains("observations", Assert.Throws<ArgumentException>(() =>
            new CalibrationModel(CalibrationFixture.Config(emulator), new double[M - 1], emulator)).Message);
    }

    [Fact]
    public void Caches_FollowWhatChanged()
    {
        CalibrationModel model = CalibrationFixture.Model(categoricals: 1);
        var emulator = (TestEmulator)model.Emulator;
        CalibrationTarget target = model.Target(0), other = model.Target(2);
        var state = new CalibrationState([0.3], [1], [0.004, 0.009], [0.1, 0.1]);
        CalibrationEvaluation e = target.EvaluateState(state);
        int calls = emulator.Predictions;

        // v: prediction and factor kept, likelihood recomputed.
        CalibrationState v = state.WithCoefficients([0.5, -0.2]);
        CalibrationEvaluation ev = target.EvaluateState(v);
        Assert.Equal(calls, emulator.Predictions);
        Assert.Same(e.Factor, ev.Factor);
        Assert.Same(e.Prediction, ev.Prediction);
        Assert.NotEqual(e.LogLikelihoodValue, ev.LogLikelihoodValue);
        Assert.Equal(target.EvaluateState(new CalibrationState([0.3], [1], [0.004, 0.009], [0.5, -0.2])).LogLikelihoodValue, ev.LogLikelihoodValue, 1e-12);
        calls = emulator.Predictions;   // the fresh state's evaluation predicted once
        // σ²: prediction kept, factor rebuilt.
        CalibrationEvaluation es = target.EvaluateState(state.WithVariances([0.02, 0.009]));
        Assert.Equal(calls, emulator.Predictions);
        Assert.NotSame(e.Factor, es.Factor);
        Assert.Same(e.Prediction, es.Prediction);
        // θ and the levels: factor kept, prediction recomputed.
        CalibrationEvaluation et = target.EvaluateState(state.WithUnit([0.8]));
        Assert.Equal(calls + 1, emulator.Predictions);
        Assert.Same(e.Factor, et.Factor);
        CalibrationEvaluation el = target.EvaluateState(state.WithLevels([2]));
        Assert.Equal(calls + 2, emulator.Predictions);
        Assert.Same(e.Factor, el.Factor);
        // Another emulator draw: both recomputed, though the state held them for the first.
        CalibrationEvaluation eo = other.EvaluateState(state);
        Assert.Equal(calls + 3, emulator.Predictions);
        Assert.NotSame(e.Factor, eo.Factor);
        CalibrationEvaluation fresh = other.EvaluateState(new CalibrationState([0.3], [1], [0.004, 0.009], [0.1, 0.1]));
        Assert.Equal(fresh.LogLikelihoodValue, eo.LogLikelihoodValue);
        Assert.NotEqual(e.LogLikelihoodValue, eo.LogLikelihoodValue);
    }

    [Theory]
    [InlineData(1.0)]
    [InlineData(0.35)]
    [InlineData(0.0)]
    public void Discrepancy_DrawsFromItsGaussianConditional(double b)
    {
        CalibrationModel model = CalibrationFixture.Model();
        CalibrationTarget target = model.Target(0);
        var state = new CalibrationState([0.6], [], [0.003, 0.006], [0.0, 0.0]);
        var block = new DiscrepancyBlock(target);
        (DenseCholesky precision, double[] rhs) = block.Conditional(state, b);
        double[] mean = precision.Solve(rhs);

        // The joint-Gaussian form: v | y ~ N(s²Dᵀ(Σ/b + s²DDᵀ)⁻¹r₀, s²I − s⁴Dᵀ(Σ/b + s²DDᵀ)⁻¹D); the prior at b = 0.
        double s2 = 0.09;
        double[,] d = model.Discrepancy!;
        double[] expectedMean = [0.0, 0.0];
        double[,] expectedCov = { { s2, 0 }, { 0, s2 } };
        if (b > 0)
        {
            double[,] sigma = DenseSigma(model, target, state);
            var c = new double[M, M];
            for (int i = 0; i < M; i++)
                for (int j = 0; j < M; j++) c[i, j] = sigma[i, j] / b + s2 * (d[i, 0] * d[j, 0] + d[i, 1] * d[j, 1]);
            var pred = new double[M];
            model.Emulator.Predict(0, model.Physical(state), [], pred);
            double[] r0 = [.. Enumerable.Range(0, M).Select(j => model.Observations[j] - pred[j] - 0.01 * j)];
            double[] cr = Inverse(c, r0);
            double[][] cd = [Inverse(c, [.. Enumerable.Range(0, M).Select(j => d[j, 0])]), Inverse(c, [.. Enumerable.Range(0, M).Select(j => d[j, 1])])];
            for (int k = 0; k < 2; k++)
            {
                expectedMean[k] = s2 * Enumerable.Range(0, M).Sum(j => d[j, k] * cr[j]);
                for (int l = 0; l < 2; l++) expectedCov[k, l] -= s2 * s2 * Enumerable.Range(0, M).Sum(j => d[j, k] * cd[l][j]);
            }
        }
        for (int k = 0; k < 2; k++)
        {
            Assert.Equal(expectedMean[k], mean[k], 1e-10);
            double[] unit = [k == 0 ? 1 : 0, k == 1 ? 1 : 0];
            double[] column = precision.Solve(unit);
            for (int l = 0; l < 2; l++) Assert.Equal(expectedCov[l, k], column[l], 1e-11);
        }

        // Its draws: each update draws afresh from the conditional, which v does not enter.
        ReversibleJumpChain<CalibrationState> chain = CalibrationFixture.Chain(target, state, b);
        var rng = new Xoshiro256PlusPlus(5);
        const int n = 6000;
        double[][] draws = new double[2][];
        draws[0] = new double[n];
        draws[1] = new double[n];
        for (int i = 0; i < n; i++)
        {
            Assert.True(block.Update(chain, rng));
            draws[0][i] = chain.Current.Coefficients[0];
            draws[1][i] = chain.Current.Coefficients[1];
            Assert.Equal(target.EvaluateState(new CalibrationState([0.6], [], [0.003, 0.006], [draws[0][i], draws[1][i]])).LogLikelihoodValue,
                         ((CalibrationEvaluation)chain.CurrentEvaluation).LogLikelihoodValue, 1e-9);
        }
        for (int k = 0; k < 2; k++)
            Assert.InRange(draws[k].Average() - expectedMean[k], -5 * Math.Sqrt(expectedCov[k, k] / n), 5 * Math.Sqrt(expectedCov[k, k] / n));
    }

    // The block's chain over ln σ² against quadrature of p(ln σ² | ·) ∝ L(σ²)^b·p(σ²)·σ² on a fine grid.
    [Theory]
    [InlineData(1.0, false)]
    [InlineData(0.4, true)]
    [InlineData(0.0, true)]
    [InlineData(1.0, true)]
    public void Variances_MoveOnLogSigma2_ToTheirConditional(double b, bool halfCauchy)
    {
        VariancePrior prior = halfCauchy ? new HalfCauchyDeviation(0.08) : new InverseGammaVariance(6, 0.08);
        CalibrationModel model = CalibrationFixture.Model(groups: CalibrationFixture.Groups(one: true, only: prior));
        CalibrationTarget target = model.Target(0);
        var start = new CalibrationState([0.45], [], [0.01], [0.1, -0.05]);

        // Wide enough for the half-Cauchy's tails at b = 0, which decay as e^(−|u|/2) in u = ln σ².
        const double lower = -45, upper = 35, h = 0.004;
        int grid = (int)((upper - lower) / h) + 1;
        var logDensity = new double[grid];
        for (int i = 0; i < grid; i++)
        {
            double u = lower + i * h;
            CalibrationEvaluation e = target.EvaluateState(start.WithVariances([Math.Exp(u)]));
            logDensity[i] = new Temperature(b, TemperingPath.Likelihood).LogDensity(e) + u;
        }
        double max = logDensity.Max();
        double[] w = [.. logDensity.Select(l => Math.Exp(l - max))];
        double total = w.Sum();
        double mean = Enumerable.Range(0, grid).Sum(i => w[i] * (lower + i * h)) / total;
        double cumulative = 0, median = upper;
        for (int i = 0; i < grid; i++) { cumulative += w[i] / total; if (cumulative >= 0.5) { median = lower + i * h; break; } }

        ReversibleJumpChain<CalibrationState> chain = CalibrationFixture.Chain(target, start, b);
        var block = new VarianceBlock(target) { Adapting = true };
        var scan = new BlockScan<ReversibleJumpChain<CalibrationState>>([block], chain, new Xoshiro256PlusPlus(17), c => c);
        for (int i = 0; i < 2000; i++) scan.Step();
        block.Adapting = false;
        const int n = 40000;
        var u2 = new double[n];
        var below = new double[n];
        for (int i = 0; i < n; i++)
        {
            scan.Step();
            u2[i] = Math.Log(chain.Current.Variances[0]);
            below[i] = u2[i] <= median ? 1 : 0;
        }
        double gapMean = u2.Average() - mean, gapMedian = below.Average() - 0.5;
        _out.WriteLine($"b = {b}, {(halfCauchy ? "half-Cauchy" : "inverse gamma")}: E[ln σ²] {u2.Average():F3} vs {mean:F3} (MCSE {CalibrationFixture.Mcse(u2):F4}); " +
                       $"P(below the median) {below.Average():F3} (MCSE {CalibrationFixture.Mcse(below):F4}); step {block.Steps[0]:F2}, acceptance {block.Accepts / (double)block.Attempts:F2}");
        Assert.InRange(gapMean, -5 * CalibrationFixture.Mcse(u2) - 2 * h, 5 * CalibrationFixture.Mcse(u2) + 2 * h);
        Assert.InRange(gapMedian, -5 * CalibrationFixture.Mcse(below) - 0.003, 5 * CalibrationFixture.Mcse(below) + 0.003);
    }

    [Theory]
    [InlineData(1.0)]
    [InlineData(0.5)]
    [InlineData(0.0)]
    public void Categorical_EnumeratedConditional_IsBruteForce(double b)
    {
        CalibrationModel model = CalibrationFixture.Model(categoricals: 2);
        CalibrationTarget target = model.Target(2);
        var state = new CalibrationState([0.55], [0, 0], [0.003, 0.006], [0.1, 0.05]);
        var block = new CategoricalBlock(target, model.Config.EnumerationCap);
        Assert.True(block.Enumerated);
        double[] probabilities = block.Conditional(state, b);

        double[,] sigma = DenseSigma(model, target, state);
        double[] delta = model.Delta(state);
        var log = new double[12];
        var pred = new double[M];
        for (int c = 0; c < 12; c++)
        {
            model.Emulator.Predict(2, model.Physical(state), [c / 4, c % 4], pred);
            log[c] = b * DenseLogNormal([.. Enumerable.Range(0, M).Select(j => model.Observations[j] - pred[j] - delta[j])], sigma);
        }
        double max = log.Max();
        double total = log.Sum(l => Math.Exp(l - max));
        for (int c = 0; c < 12; c++) Assert.Equal(Math.Exp(log[c] - max) / total, probabilities[c], 1e-12);

        ReversibleJumpChain<CalibrationState> chain = CalibrationFixture.Chain(target, state, b);
        var rng = new Xoshiro256PlusPlus(23);
        const int n = 20000;
        var counts = new int[12];
        for (int i = 0; i < n; i++)
        {
            block.Update(chain, rng);
            counts[chain.Current.Levels[0] * 4 + chain.Current.Levels[1]]++;
            if (i % 997 == 0)
                Assert.Equal(target.EvaluateState(new CalibrationState([0.55], chain.Current.Levels.ToArray(), [0.003, 0.006], [0.1, 0.05])).LogLikelihoodValue,
                             ((CalibrationEvaluation)chain.CurrentEvaluation).LogLikelihoodValue, 1e-10);
        }
        for (int c = 0; c < 12; c++)
        {
            double p = probabilities[c], se = Math.Sqrt(p * (1 - p) / n);
            Assert.InRange(counts[c] / (double)n - p, -5 * se - 1e-4, 5 * se + 1e-4);
        }
    }

    [Fact]
    public void Categorical_AboveTheCap_SingleSiteUpdatesReachTheEnumeratedLaw()
    {
        CalibrationModel capped = CalibrationFixture.Model(categoricals: 2, cap: 8);
        CalibrationTarget target = capped.Target(1);
        var block = new CategoricalBlock(target, capped.Config.EnumerationCap);
        Assert.False(block.Enumerated);
        Assert.Equal("levels (single site)", block.Name);
        var state = new CalibrationState([0.55], [0, 0], [0.01, 0.02], [0.1, 0.05]);
        double[] exact = new CategoricalBlock(target, 1024).Conditional(state, 1.0);

        ReversibleJumpChain<CalibrationState> chain = CalibrationFixture.Chain(target, state, 1.0);
        var scan = new BlockScan<ReversibleJumpChain<CalibrationState>>([block], chain, new Xoshiro256PlusPlus(29), c => c);
        const int n = 60000;
        var indicators = new double[12][];
        for (int c = 0; c < 12; c++) indicators[c] = new double[n];
        for (int i = 0; i < n; i++)
        {
            scan.Step();
            indicators[chain.Current.Levels[0] * 4 + chain.Current.Levels[1]][i] = 1;
        }
        int modal = Array.IndexOf(exact, exact.Max());
        double tau = IntegratedAutocorrelation(indicators[modal]);
        _out.WriteLine($"single site: acceptance {block.SiteAccepts / (double)block.SiteAttempts:F2}; modal combination {modal} (p {exact[modal]:F3}), " +
                       $"integrated autocorrelation time {tau:F1} scans, ESS {n / tau:F0} of {n}");
        for (int c = 0; c < 12; c++)
        {
            double mcse = CalibrationFixture.Mcse(indicators[c]);
            Assert.InRange(indicators[c].Average() - exact[c], -5 * mcse - 0.002, 5 * mcse + 0.002);
        }
    }

    // Sokal's windowed estimate: 1 + 2Σ ρ_k up to the first window with k ≥ 5τ.
    private static double IntegratedAutocorrelation(double[] x)
    {
        double mean = x.Average(), c0 = x.Sum(v => (v - mean) * (v - mean)) / x.Length;
        if (c0 == 0) return 1;
        double tau = 1;
        for (int k = 1; k < x.Length / 4; k++)
        {
            double ck = 0;
            for (int i = 0; i + k < x.Length; i++) ck += (x[i] - mean) * (x[i + k] - mean);
            tau += 2 * ck / x.Length / c0;
            if (k >= 5 * tau) break;
        }
        return tau;
    }

    /// <summary>R's (1 + σ²)⁻¹ kernel on σ² (C-O1's <c>s2.df</c> switch): the defect fixture, improper.</summary>
    private sealed class RsKernel : IVariancePrior
    {
        public double LogDensity(double sigma2) => sigma2 > 0 ? -Math.Log(1 + sigma2) : double.NegativeInfinity;
        public double Sample(Xoshiro256PlusPlus rng) => throw new NotSupportedException("An improper prior has no draws.");
        public bool IsProper => false;
    }

    [Fact]
    public void ProprietyRule_ABZeroRungNeedsEveryPriorProper()
    {
        TestEmulator emulator = CalibrationFixture.Emulator();
        CalibrationConfig config = CalibrationFixture.Config(emulator);
        double[] y = CalibrationFixture.Observations(emulator);
        var proper = new CalibrationModel(config, y, emulator);
        var state = new CalibrationState([0.5], [], [0.01, 0.01], [0.0, 0.0]);
        Assert.True(proper.PriorsProper);
        Assert.Equal(0.0, proper.Target(0).EvaluateState(state).LogLikelihood(0.0));

        var improper = new CalibrationModel(config, y, emulator, new MeasurementModel([0, 0, 0, 0, 1, 1, 1, 1], [new RsKernel(), new InverseGammaVariancePrior(8, 0.05)]));
        Assert.False(improper.PriorsProper);
        CalibrationEvaluation e = improper.Target(0).EvaluateState(state);
        Assert.Throws<ArgumentOutOfRangeException>(() => e.LogLikelihood(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Temperature(0.0, TemperingPath.Likelihood).LogDensity(e));
        Assert.True(double.IsFinite(e.LogLikelihood(0.2)));
    }

    [Fact]
    public void RsKernel_IsImproperAndUnitDependent_WhereTheHalfCauchyIsNeither()
    {
        // Improper: ∫₀^X (1 + σ²)⁻¹ dσ² = ln(1 + X), unbounded.
        foreach (double x in new[] { 1e2, 1e6, 1e12 })
        {
            const double h = 0.001;
            double upper = Math.Log(x), mass = 0;
            for (double u = -40; u < upper; u += h) mass += Math.Exp(new RsKernel().LogDensity(Math.Exp(u)) + u) * h;
            Assert.Equal(Math.Log(1 + x), mass, 2e-3);
        }
        // Unit-dependent: σ's posterior from ten draws of N(0, σ²), the data in units c times larger, reported in the
        // original units. Under R's kernel the median moves (nearly flat on σ² where σ² ≪ 1, nearly 1/σ² where σ² ≫ 1);
        // under the half-Cauchy with its scale in the data's units it does not.
        double[] data = [0.3, -0.12, 0.25, 0.05, -0.4, 0.18, -0.07, 0.33, -0.21, 0.11];
        double Median(Func<double, double, double> logPrior, double c)
        {
            double ss = data.Sum(v => v * v) * c * c, centre = Math.Log(ss / data.Length);
            const double h = 0.0005;
            var w = new List<(double U, double W)>();
            for (double u = centre - 25; u < centre + 25; u += h)
                w.Add((u, Math.Exp(logPrior(Math.Exp(u), c) + u - 0.5 * data.Length * u - 0.5 * ss / Math.Exp(u))));
            double total = w.Sum(t => t.W), cumulative = 0;
            foreach ((double u, double weight) in w)
                if ((cumulative += weight / total) >= 0.5) return Math.Exp(0.5 * u) / c;
            return double.NaN;
        }
        var rs = new RsKernel();
        double rsSmall = Median((s2, _) => rs.LogDensity(s2), 1e-3), rsLarge = Median((s2, _) => rs.LogDensity(s2), 1e3);
        double hcSmall = Median((s2, c) => new HalfCauchyVariancePrior(0.2 * c).LogDensity(s2), 1e-3);
        double hcLarge = Median((s2, c) => new HalfCauchyVariancePrior(0.2 * c).LogDensity(s2), 1e3);
        _out.WriteLine($"posterior median of σ in the original units: R's kernel {rsSmall:F4} (c = 10⁻³) vs {rsLarge:F4} (c = 10³); " +
                       $"half-Cauchy {hcSmall:F4} vs {hcLarge:F4}");
        Assert.True(Math.Abs(rsSmall - rsLarge) / rsSmall > 0.05);
        Assert.Equal(hcSmall, hcLarge, 1e-4);
    }
}
