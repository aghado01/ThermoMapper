using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The cut SMC's bridged target and its kernel (stage 8, C5; D09): the bridge's density against the geometric mixture
/// of its two conditionals; its states' caches holding both draws; the prior draw; the kernel's scan at t = 0 and 1
/// against the single-draw conditional by quadrature; and at t = 0.3 and 0.7 the whole scan against the bridged law by
/// quadrature with v integrated out, then each block alone against its own bridged conditional: the walk and the
/// variance step by quadrature, the levels' enumeration against brute force, the single-site steps above the cap
/// against enumeration, and v's bridged Gaussian against its closed form.
/// </summary>
public sealed class BassCalibrationBridgeTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationBridgeTests(ITestOutputHelper output) => _out = output;

    private const int M = CalibrationFixture.M;

    private static SmcStage<CalibrationState> Stage(ISamplingTarget<CalibrationState> target, double t) => new(target, new Temperature(t, TemperingPath.Likelihood));

    [Fact]
    public void Bridge_IsTheGeometricMixtureOfItsTwoConditionals()
    {
        CalibrationModel model = CalibrationFixture.Model(categoricals: 1);
        CalibrationTarget a = model.Target(0), b = model.Target(2);
        var bridge = new CalibrationBridge(a, b);
        var rng = new Xoshiro256PlusPlus(3);
        double worst = 0;
        for (int i = 0; i < 20; i++)
        {
            CalibrationState x = model.DrawPrior(rng);
            double pa = Temperature.LogPosterior(a.EvaluateState(new CalibrationState(x.Unit.ToArray(), x.Levels.ToArray(), x.Variances.ToArray(), x.Coefficients.ToArray())));
            double pb = Temperature.LogPosterior(b.EvaluateState(new CalibrationState(x.Unit.ToArray(), x.Levels.ToArray(), x.Variances.ToArray(), x.Coefficients.ToArray())));
            ITargetEvaluation e = bridge.Evaluate(x);
            foreach (double t in new[] { 0.0, 0.3, 0.7, 1.0 })
            {
                double expected = (1 - t) * pa + t * pb, actual = Stage(bridge, t).LogDensity(e);
                worst = Math.Max(worst, Math.Abs(actual - expected) / (1 + Math.Abs(expected)));
                Assert.Equal(expected, actual, 1e-12 * (1 + Math.Abs(expected)));
            }
            // A step from t to t′ on the bridge reweights by (t′ − t)·Δ.
            double delta = ((BridgeEvaluation)e).Difference;
            Assert.Equal(pb - pa, delta, 1e-12 * (1 + Math.Abs(pa)));
            Assert.Equal(0.4 * delta, Stage(bridge, 0.7).LogDensity(e) - Stage(bridge, 0.3).LogDensity(e), 1e-12 * (1 + Math.Abs(pa)));
            Assert.Equal(delta, Stage(bridge, 0.5).Temperature.LogDensitySlope(e));
        }
        _out.WriteLine($"bridge against (1 − t)·log p̃_a + t·log p̃_b: worst relative gap {worst:E1}");
        ITargetEvaluation off = bridge.Evaluate(new CalibrationState([1.2], [0], [0.01, 0.01], [0, 0]));
        Assert.True(double.IsNegativeInfinity(Stage(bridge, 0.0).LogDensity(off)));
        Assert.True(double.IsNegativeInfinity(Stage(bridge, 0.6).LogDensity(off)));
        CalibrationModel other = CalibrationFixture.Model(categoricals: 1);
        Assert.Throws<ArgumentException>(() => new CalibrationBridge(a, other.Target(1)));
    }

    [Fact]
    public void Caches_ABridgeReadsBothDrawsFromTheState()
    {
        CalibrationModel model = CalibrationFixture.Model(categoricals: 1);
        var emulator = (TestEmulator)model.Emulator;
        CalibrationTarget a = model.Target(0), b = model.Target(2), third = model.Target(1);
        var bridge = new CalibrationBridge(a, b);
        var state = new CalibrationState([0.4], [1], [0.004, 0.009], [0.1, 0.1]);
        int calls = emulator.Predictions;
        BridgeEvaluation e = bridge.EvaluateState(state);
        Assert.Equal(calls + 2, emulator.Predictions);
        BridgeEvaluation again = bridge.EvaluateState(state);
        Assert.Equal(calls + 2, emulator.Predictions);
        Assert.Same(e.From.Factor, again.From.Factor);
        Assert.Same(e.To.Factor, again.To.Factor);
        // v: both kept; σ²: both predictions kept, both factors rebuilt; θ: both factors kept, both predictions made.
        BridgeEvaluation ev = bridge.EvaluateState(state.WithCoefficients([0.3, -0.1]));
        Assert.Equal(calls + 2, emulator.Predictions);
        Assert.Same(e.To.Prediction, ev.To.Prediction);
        BridgeEvaluation es = bridge.EvaluateState(state.WithVariances([0.02, 0.009]));
        Assert.Equal(calls + 2, emulator.Predictions);
        Assert.NotSame(e.From.Factor, es.From.Factor);
        Assert.NotSame(e.To.Factor, es.To.Factor);
        BridgeEvaluation et = bridge.EvaluateState(state.WithUnit([0.7]));
        Assert.Equal(calls + 4, emulator.Predictions);
        Assert.Same(e.From.Factor, et.From.Factor);
        Assert.Same(e.To.Factor, et.To.Factor);
        // A fork carries the values and caches; a third draw displaces the older of the two.
        CalibrationState fork = new CutSmcParticles(model).Fork(state);
        Assert.NotSame(state, fork);
        Assert.Same(e.From.Prediction, bridge.EvaluateState(fork).From.Prediction);
        Assert.Equal(calls + 4, emulator.Predictions);
        third.EvaluateState(fork);
        Assert.Equal(calls + 5, emulator.Predictions);
        BridgeEvaluation displaced = bridge.EvaluateState(fork);   // e_a displaced e_b, which then displaced the third
        Assert.Equal(calls + 7, emulator.Predictions);
        Assert.Equal(e.From.LogLikelihoodValue, displaced.From.LogLikelihoodValue);
        Assert.Equal(e.To.LogLikelihoodValue, displaced.To.LogLikelihoodValue);
        Assert.Same(e.From.Prediction, bridge.EvaluateState(state).From.Prediction);   // the original's caches are its own
    }

    [Fact]
    public void DrawPrior_IsThePrior_AndRefusedWhenImproper()
    {
        CalibrationModel model = CalibrationFixture.Model(categoricals: 1);
        var rng = new Xoshiro256PlusPlus(5);
        const int n = 40000;
        CalibrationState[] draws = [.. Enumerable.Range(0, n).Select(_ => model.DrawPrior(rng))];
        double[] theta = [.. draws.Select(d => d.Unit[0])];
        Assert.InRange(theta.Average() - 0.5, -5 * Math.Sqrt(1 / 12.0 / n), 5 * Math.Sqrt(1 / 12.0 / n));
        Assert.InRange(theta.Select(t => (t - 0.5) * (t - 0.5)).Average() - 1 / 12.0, -0.002, 0.002);
        for (int l = 0; l < 3; l++)
            Assert.InRange(draws.Count(d => d.Levels[0] == l) / (double)n - 1 / 3.0, -0.012, 0.012);
        for (int k = 0; k < 2; k++)
        {
            double[] v = [.. draws.Select(d => d.Coefficients[k])];
            Assert.InRange(v.Average(), -5 * 0.3 / Math.Sqrt(n), 5 * 0.3 / Math.Sqrt(n));
            Assert.InRange(v.Select(x => x * x).Average() / 0.09, 1 - 5 * Math.Sqrt(2.0 / n), 1 + 5 * Math.Sqrt(2.0 / n));
        }
        // σ²_g from each group's prior: the IG's ln σ² mean, ln(νs²/2) − ψ(ν/2), at ν = 8 and s = 0.05.
        double expected = Math.Log(4 * 0.0025) - (1 + 1 / 2.0 + 1 / 3.0 - 0.5772156649015329);
        double[] u = [.. draws.Select(d => Math.Log(d.Variances[0]))];
        double trigamma = Math.PI * Math.PI / 6 - 1 - 1 / 4.0 - 1 / 9.0;
        Assert.InRange(u.Average() - expected, -5 * Math.Sqrt(trigamma / n), 5 * Math.Sqrt(trigamma / n));

        TestEmulator emulator = CalibrationFixture.Emulator();
        var improper = new CalibrationModel(CalibrationFixture.Config(emulator), CalibrationFixture.Observations(emulator), emulator,
            new MeasurementModel([0, 0, 0, 0, 1, 1, 1, 1], [new ImproperPrior(), new InverseGammaVariancePrior(8, 0.05)]));
        Assert.Throws<InvalidOperationException>(() => improper.DrawPrior(rng));
        Assert.Throws<ArgumentException>(() => new CutSmcParticles(improper));
    }

    private sealed class ImproperPrior : IVariancePrior
    {
        public double LogDensity(double sigma2) => sigma2 > 0 ? -Math.Log(1 + sigma2) : double.NegativeInfinity;
        public double Sample(Xoshiro256PlusPlus rng) => throw new NotSupportedException();
        public bool IsProper => false;
    }

    // One continuous θ, one categorical of 3 levels, one inverse-gamma group, the functional discrepancy (the sampler
    // tests' fixture), bridged from draw 0 to draw 2.
    private static CalibrationModel Model() => BassCalibrationSamplerTests.Model(categoricals: 1, new CalibrationSampling());

    // A population near e's conditional, from a short run of the stage-7 sampler, for the kernel's scales.
    private static CalibrationState[] Population(CalibrationModel model, int draw)
    {
        var sampling = new CalibrationSampling { MasterSeed = 13, BurnIn = 1000, Samples = 400, DecorrelationEvery = 25 };
        return [.. CalibrationSampler.Run(model, draw, sampling).Draws.Select(d =>
            new CalibrationState([CalibrationFixture.Theta.ToUnit(d.Theta[0])], d.Levels, d.Variances, d.Coefficients))];
    }

    private static CutSmcKernel Prepared(CalibrationModel model, SmcStage<CalibrationState> stage, CalibrationState[] population)
    {
        var kernel = new CutSmcKernel(model);
        kernel.Prepare(stage, population, [.. Enumerable.Repeat(1.0 / population.Length, population.Length)]);
        return kernel;
    }

    // The kernel run as one chain from `start`, keeping every sweep's state after `burn`.
    private static List<CalibrationState> Chain(CutSmcKernel kernel, SmcStage<CalibrationState> stage, CalibrationState start, int burn, int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        CalibrationState x = start;
        ITargetEvaluation e = stage.Target.Evaluate(x);
        var kept = new List<CalibrationState>(n);
        for (int i = 0; i < burn + n; i++)
        {
            SmcMove<CalibrationState> move = kernel.Move(x, e, stage, rng);
            (x, e) = (move.State, move.Evaluation);
            if (i >= burn) kept.Add(x);
        }
        return kept;
    }

    private void Compare(string label, (double MeanTheta, double MedianTheta, double MeanLogVariance, double[] Levels) q, IReadOnlyList<CalibrationState> draws)
    {
        double[] theta = [.. draws.Select(s => s.Unit[0])];
        double[] below = [.. theta.Select(t => t <= q.MedianTheta ? 1.0 : 0.0)];
        double[] u = [.. draws.Select(s => Math.Log(s.Variances[0]))];
        _out.WriteLine($"{label}: E[θ] {theta.Average():F4} vs {q.MeanTheta:F4} (MCSE {CalibrationFixture.Mcse(theta):F4}); P(θ ≤ median) {below.Average():F3} " +
                       $"(MCSE {CalibrationFixture.Mcse(below):F3}); E[ln σ²] {u.Average():F3} vs {q.MeanLogVariance:F3} (MCSE {CalibrationFixture.Mcse(u):F3}); " +
                       $"levels {string.Join(" ", Enumerable.Range(0, q.Levels.Length).Select(c => $"{draws.Count(s => s.Levels[0] == c) / (double)draws.Count:F3}/{q.Levels[c]:F3}"))}");
        const double grid = 2.0 / 600;
        Assert.InRange(theta.Average() - q.MeanTheta, -5 * CalibrationFixture.Mcse(theta) - grid, 5 * CalibrationFixture.Mcse(theta) + grid);
        Assert.InRange(below.Average() - 0.5, -5 * CalibrationFixture.Mcse(below) - 0.01, 5 * CalibrationFixture.Mcse(below) + 0.01);
        Assert.InRange(u.Average() - q.MeanLogVariance, -5 * CalibrationFixture.Mcse(u) - 0.02, 5 * CalibrationFixture.Mcse(u) + 0.02);
        for (int c = 0; c < q.Levels.Length; c++)
        {
            double[] at = [.. draws.Select(s => s.Levels[0] == c ? 1.0 : 0.0)];
            Assert.InRange(at.Average() - q.Levels[c], -5 * CalibrationFixture.Mcse(at) - 0.005, 5 * CalibrationFixture.Mcse(at) + 0.005);
        }
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(1.0)]
    public void Kernel_AtTheBridgesEnds_IsTheSingleDrawConditional(double t)
    {
        CalibrationModel model = Model();
        CalibrationTarget a = model.Target(0), b = model.Target(2);
        CalibrationTarget end = t == 0 ? a : b;
        SmcStage<CalibrationState> stage = Stage(new CalibrationBridge(a, b), t);
        CalibrationState[] population = Population(model, end.Draw);
        CutSmcKernel kernel = Prepared(model, stage, population);
        List<CalibrationState> draws = Chain(kernel, stage, population[0], 500, 30000, 31);
        Compare($"bridge at t = {t}", BassCalibrationSamplerTests.Quadrature(model, end, 1.0), draws);
    }

    [Fact]
    public void Kernel_OnTheFirstLeg_IsTheRungsConditional()
    {
        // The first leg's stages: one draw on D05's likelihood path at a rung, the same scan aimed at it.
        CalibrationModel model = Model();
        CalibrationTarget target = model.Target(2);
        SmcStage<CalibrationState> rung = Stage(target, 0.5);
        CalibrationState[] population = Population(model, 2);
        Compare("e_2 at b = 0.5", BassCalibrationSamplerTests.Quadrature(model, target, 0.5), Chain(Prepared(model, rung, population), rung, population[0], 500, 30000, 37));
    }

    /// <summary>
    /// Quadrature of the bridged law over (θ, ln σ², level), v integrated out: the stage's density
    /// log π + w_a·ℓ_a + w_b·ℓ_b is Gaussian in v, so with P = I/s_v² + Σ_k w_k·DᵀΣ_k⁻¹D and h = Σ_k w_k·DᵀΣ_k⁻¹r_k,
    /// r_k = y − prediction_k − μ_δ, its v-marginal is π(σ²)·exp(−½Σ_k w_k(r_kᵀΣ_k⁻¹r_k + ln|Σ_k|) + ½hᵀP⁻¹h − ½ln|P|),
    /// times e^u for u = ln σ². Dense Σ_k throughout.
    /// </summary>
    internal static (double MeanTheta, double MedianTheta, double MeanLogVariance, double[] Levels) BridgeQuadrature(
        CalibrationModel model, CalibrationTarget a, CalibrationTarget b, double wa, double wb)
    {
        int levels = model.Layout.Categorical.Count == 0 ? 1 : model.Layout.Categorical[0].Levels.Length;
        const int thetaPoints = 600, uPoints = 260;
        const double uLower = -14, uUpper = 0;
        double hu = (uUpper - uLower) / uPoints;
        double[,] d = model.Discrepancy!;
        double s2 = model.DiscrepancyScale * model.DiscrepancyScale;
        CalibrationTarget[] targets = [a, b];
        double[] weights = [wa, wb];
        double[][,] fixedParts = [.. targets.Select(t => FixedSigma(model, t))];
        var factors = new DenseCholesky[2, uPoints];
        var solvedD = new double[2, uPoints][][];
        var pInverse = new double[uPoints][,];
        var baseWeight = new double[uPoints];
        for (int iu = 0; iu < uPoints; iu++)
        {
            double u = uLower + (iu + 0.5) * hu;
            baseWeight[iu] = model.Measurement.LogPrior(0, Math.Exp(u)) + u;
            double[,] p = { { 1 / s2, 0 }, { 0, 1 / s2 } };
            for (int k = 0; k < 2; k++)
            {
                var sigma = (double[,])fixedParts[k].Clone();
                for (int j = 0; j < M; j++) sigma[j, j] += Math.Exp(u);
                factors[k, iu] = new DenseCholesky(M);
                factors[k, iu].Decompose(sigma);
                baseWeight[iu] -= 0.5 * weights[k] * factors[k, iu].LogDet;
                solvedD[k, iu] = [factors[k, iu].Solve([.. Enumerable.Range(0, M).Select(j => d[j, 0])]), factors[k, iu].Solve([.. Enumerable.Range(0, M).Select(j => d[j, 1])])];
                for (int r = 0; r < 2; r++)
                    for (int c = 0; c < 2; c++) p[r, c] += weights[k] * Enumerable.Range(0, M).Sum(j => d[j, r] * solvedD[k, iu][c][j]);
            }
            double det = p[0, 0] * p[1, 1] - p[0, 1] * p[1, 0];
            baseWeight[iu] -= 0.5 * Math.Log(det);
            pInverse[iu] = new[,] { { p[1, 1] / det, -p[0, 1] / det }, { -p[1, 0] / det, p[0, 0] / det } };
        }
        var logW = new double[levels, thetaPoints, uPoints];
        var pred = new double[M];
        double max = double.NegativeInfinity;
        for (int c = 0; c < levels; c++)
            for (int it = 0; it < thetaPoints; it++)
            {
                double unit = (it + 0.5) / thetaPoints;
                var r = new double[2][];
                for (int k = 0; k < 2; k++)
                {
                    model.Emulator.Predict(targets[k].Draw, [CalibrationFixture.Theta.ToPhysical(unit)], levels > 1 ? [c] : [], pred);
                    r[k] = [.. Enumerable.Range(0, M).Select(j => model.Observations[j] - pred[j] - 0.01 * j)];
                }
                for (int iu = 0; iu < uPoints; iu++)
                {
                    double lw = baseWeight[iu];
                    var h = new double[2];
                    for (int k = 0; k < 2; k++)
                    {
                        if (weights[k] == 0) continue;
                        lw -= 0.5 * weights[k] * factors[k, iu].InverseQuadraticForm(r[k]);
                        for (int q = 0; q < 2; q++) h[q] += weights[k] * Enumerable.Range(0, M).Sum(j => solvedD[k, iu][q][j] * r[k][j]);
                    }
                    double[,] pi = pInverse[iu];
                    lw += 0.5 * (h[0] * (pi[0, 0] * h[0] + pi[0, 1] * h[1]) + h[1] * (pi[1, 0] * h[0] + pi[1, 1] * h[1]));
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

    // The part of draw e's Σ that σ² does not enter: its components and the residual set's second moment, dense.
    internal static double[,] FixedSigma(CalibrationModel model, CalibrationTarget target)
    {
        var noise = (ComponentNoise)target.Noise;
        double[,] residuals = model.Emulator.TruncationResiduals()!;
        int n = residuals.GetLength(0);
        var sigma = new double[M, M];
        for (int j = 0; j < M; j++)
            for (int l = 0; l < M; l++)
            {
                double s = 0;
                for (int k = 0; k < noise.Variances.Length; k++) s += noise.Factor[j, k] * noise.Variances[k] * noise.Factor[l, k];
                for (int i = 0; i < n; i++) s += residuals[i, j] * residuals[i, l] / n;
                sigma[j, l] = s;
            }
        return sigma;
    }

    // Σ_e at the state's group variances, dense.
    private static double[,] DenseSigma(CalibrationModel model, CalibrationTarget target, CalibrationState state)
    {
        double[,] sigma = FixedSigma(model, target);
        var perOutput = new double[M];
        model.Measurement.OutputVariances(state.Variances, perOutput);
        for (int j = 0; j < M; j++) sigma[j, j] += perOutput[j];
        return sigma;
    }

    [Fact]
    public void BridgeQuadrature_AtAnEnd_IsTheSamplersQuadrature()
    {
        CalibrationModel model = Model();
        var mine = BridgeQuadrature(model, model.Target(0), model.Target(2), 0.0, 1.0);
        var theirs = BassCalibrationSamplerTests.Quadrature(model, model.Target(2), 1.0);
        Assert.Equal(theirs.MeanTheta, mine.MeanTheta, 1e-9);
        Assert.Equal(theirs.MeanLogVariance, mine.MeanLogVariance, 1e-9);
        for (int c = 0; c < 3; c++) Assert.Equal(theirs.Levels[c], mine.Levels[c], 1e-9);
        var tempered = BridgeQuadrature(model, model.Target(0), model.Target(2), 0.0, 0.4);
        Assert.Equal(BassCalibrationSamplerTests.Quadrature(model, model.Target(2), 0.4).MeanTheta, tempered.MeanTheta, 1e-9);
    }

    [Theory]
    [InlineData(0.3)]
    [InlineData(0.7)]
    public void Kernel_InsideTheBridge_TheScansLawIsTheBridgedLaw(double t)
    {
        CalibrationModel model = Model();
        CalibrationTarget a = model.Target(0), b = model.Target(2);
        SmcStage<CalibrationState> stage = Stage(new CalibrationBridge(a, b), t);
        CalibrationState[] population = Population(model, 0);
        List<CalibrationState> draws = Chain(Prepared(model, stage, population), stage, population[^1], 500, 30000, 41);
        Compare($"bridge at t = {t}", BridgeQuadrature(model, a, b, 1 - t, t), draws);
    }

    [Theory]
    [InlineData(0.3)]
    [InlineData(0.7)]
    public void EachBlock_LeavesTheBridgeInvariant(double t)
    {
        CalibrationModel model = Model();
        CalibrationTarget a = model.Target(0), b = model.Target(2);
        var bridge = new CalibrationBridge(a, b);
        SmcStage<CalibrationState> stage = Stage(bridge, t);
        CutSmcKernel kernel = Prepared(model, stage, Population(model, 0));
        LikelihoodTerm[] terms = kernel.Terms(stage);
        Assert.Equal([(0, 1 - t), (2, t)], terms.Select(x => (x.Target.Draw, x.Weight)));
        var start = new CalibrationState([0.45], [1], [0.006], [0.05, -0.02]);
        WalkAndVariance($"t = {t}", kernel, stage, start);

        // The levels alone: the bridged enumeration against brute force, and the draws' frequencies.
        double[] probabilities = CategoricalBlock.Enumerate(model, start, terms, [[0], [1], [2]]).Probabilities;
        double[,] sa = DenseSigma(model, a, start), sb = DenseSigma(model, b, start);
        double[] delta = model.Delta(start);
        double[] log = new double[3];
        var pred = new double[M];
        for (int c = 0; c < 3; c++)
        {
            model.Emulator.Predict(0, model.Physical(start), [c], pred);
            log[c] = (1 - t) * DenseLogNormal([.. Enumerable.Range(0, M).Select(j => model.Observations[j] - pred[j] - delta[j])], sa);
            model.Emulator.Predict(2, model.Physical(start), [c], pred);
            log[c] += t * DenseLogNormal([.. Enumerable.Range(0, M).Select(j => model.Observations[j] - pred[j] - delta[j])], sb);
        }
        double top = log.Max(), sum = log.Sum(l => Math.Exp(l - top));
        for (int c = 0; c < 3; c++) Assert.Equal(Math.Exp(log[c] - top) / sum, probabilities[c], 1e-12);
        double[] levelDraws = Run(kernel, stage, start, 20000, 53, (ref CutSmcKernel.Sweep s, Xoshiro256PlusPlus r) => { kernel.Levels(ref s, terms, r); return true; }, s => s.Levels[0]);
        for (int c = 0; c < 3; c++)
        {
            double f = levelDraws.Count(l => l == c) / 20000.0, se = Math.Sqrt(probabilities[c] * (1 - probabilities[c]) / 20000);
            Assert.InRange(f - probabilities[c], -5 * se - 1e-4, 5 * se + 1e-4);
        }

        // v alone: the bridged Gaussian against its closed form from the dense Σ's, and the draws' means.
        (DenseCholesky precision, double[] rhs) = DiscrepancyBlock.Conditional(model, start, terms);
        double[,] dd = model.Discrepancy!;
        double s2 = model.DiscrepancyScale * model.DiscrepancyScale;
        double[,] p = { { 1 / s2, 0 }, { 0, 1 / s2 } };
        var h = new double[2];
        foreach ((CalibrationTarget target, double w, double[,] sigma) in new[] { (a, 1 - t, sa), (b, t, sb) })
        {
            model.Emulator.Predict(target.Draw, model.Physical(start), [1], pred);
            double[] r0 = [.. Enumerable.Range(0, M).Select(j => model.Observations[j] - pred[j] - 0.01 * j)];
            double[][] sd = [Inverse(sigma, [.. Enumerable.Range(0, M).Select(j => dd[j, 0])]), Inverse(sigma, [.. Enumerable.Range(0, M).Select(j => dd[j, 1])])];
            for (int k = 0; k < 2; k++)
            {
                h[k] += w * Enumerable.Range(0, M).Sum(j => sd[k][j] * r0[j]);
                for (int l = 0; l < 2; l++) p[k, l] += w * Enumerable.Range(0, M).Sum(j => dd[j, l] * sd[k][j]);
            }
        }
        double[] expectedMean = Inverse(p, h);
        double[] gotMean = precision.Solve(rhs);
        for (int k = 0; k < 2; k++)
        {
            Assert.Equal(expectedMean[k], gotMean[k], 1e-10);
            double[] column = precision.Solve([k == 0 ? 1 : 0, k == 1 ? 1 : 0]);
            double[] expectedColumn = Inverse(p, [k == 0 ? 1 : 0, k == 1 ? 1 : 0]);
            for (int l = 0; l < 2; l++) Assert.Equal(expectedColumn[l], column[l], 1e-11);
        }
        for (int k = 0; k < 2; k++)
        {
            int index = k;
            double[] v = Run(kernel, stage, start, 20000, 59 + k, (ref CutSmcKernel.Sweep s, Xoshiro256PlusPlus r) => { kernel.Discrepancy(ref s, terms, r); return true; }, s => s.Coefficients[index]);
            double se = Math.Sqrt(Inverse(p, [k == 0 ? 1 : 0, k == 1 ? 1 : 0])[k] / v.Length);
            Assert.InRange(v.Average() - expectedMean[k], -5 * se, 5 * se);
        }
    }

    // θ's walk alone against quadrature of the bridged density over θ's unit coordinate, then the (first) variance step
    // alone against quadrature over u = ln σ² with its Jacobian, everything else held at `start`.
    private void WalkAndVariance(string label, CutSmcKernel kernel, SmcStage<CalibrationState> stage, CalibrationState start)
    {
        const int n = 40000, points = 4000;
        double[] grid = [.. Enumerable.Range(0, points).Select(i => (i + 0.5) / points)];
        double[] logDensity = [.. grid.Select(u => stage.LogDensity(stage.Target.Evaluate(start.WithUnit([u]))))];
        (double mean, double median) = Moments(grid, logDensity);
        double[] walked = Run(kernel, stage, start, n, 43, (ref CutSmcKernel.Sweep s, Xoshiro256PlusPlus r) => kernel.Walk(ref s, r), s => s.Unit[0]);
        Check($"{label}, walk", walked, mean, median, 1.0 / points);

        double[] us = [.. Enumerable.Range(0, 20000).Select(i => -16 + (i + 0.5) * 0.001)];
        double[] logU = [.. us.Select(u =>
        {
            double[] variances = start.Variances.ToArray();
            variances[0] = Math.Exp(u);
            return stage.LogDensity(stage.Target.Evaluate(start.WithVariances(variances))) + u;
        })];
        (mean, median) = Moments(us, logU);
        double[] varied = Run(kernel, stage, start, n, 47, (ref CutSmcKernel.Sweep s, Xoshiro256PlusPlus r) => kernel.Variances(ref s, r), s => Math.Log(s.Variances[0]));
        Check($"{label}, variance", varied, mean, median, 0.002);
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void Multivariate_AndWithoutADiscrepancy_TheBlocksKeepTheBridge(bool multivariate)
    {
        // The multivariate formulation (no v to draw, a fixed Σ_δ and mean), or the functional one without a discrepancy
        // on an emulator without a basis (σ²_e·I): the bridge's density, then the walk and the variance step at t = 0.5.
        CovarianceFormulation covariance = multivariate
            ? new MultivariateFormulation(
                [.. Enumerable.Range(0, M).Select(i => Enumerable.Range(0, M).Select(j => Math.Pow(0.4, Math.Abs(i - j))).ToArray())],
                [.. Enumerable.Range(0, M).Select(i => Enumerable.Range(0, M).Select(j => (i == j ? 1e-4 : 0.0) + 0.002 * Math.Exp(-0.5 * (i - j) * (i - j))).ToArray())],
                [.. Enumerable.Range(0, M).Select(j => 0.01 * j)])
            : new FunctionalFormulation();
        TestEmulator emulator = CalibrationFixture.Emulator(residuals: multivariate, independent: !multivariate);
        CalibrationConfig config = CalibrationFixture.Config(emulator, CalibrationFixture.Groups(one: true), covariance);
        var model = new CalibrationModel(config, CalibrationFixture.Observations(emulator), emulator);
        Assert.Equal(0, model.Coefficients);
        CalibrationTarget a = model.Target(0), b = model.Target(2);
        var bridge = new CalibrationBridge(a, b);
        var rng = new Xoshiro256PlusPlus(7);
        for (int i = 0; i < 10; i++)
        {
            CalibrationState x = model.DrawPrior(rng);
            double pa = Temperature.LogPosterior(a.EvaluateState(new CalibrationState(x.Unit.ToArray(), [], x.Variances.ToArray(), [])));
            double pb = Temperature.LogPosterior(b.EvaluateState(new CalibrationState(x.Unit.ToArray(), [], x.Variances.ToArray(), [])));
            foreach (double t in new[] { 0.3, 0.7 })
                Assert.Equal((1 - t) * pa + t * pb, Stage(bridge, t).LogDensity(bridge.Evaluate(x)), 1e-12 * (1 + Math.Abs(pa)));
        }
        SmcStage<CalibrationState> stage = Stage(bridge, 0.5);
        var sampling = new CalibrationSampling { MasterSeed = 17, BurnIn = 1000, Samples = 400, DecorrelationEvery = 25 };
        CalibrationState[] population = [.. CalibrationSampler.Run(model, 0, sampling).Draws.Select(d =>
            new CalibrationState([CalibrationFixture.Theta.ToUnit(d.Theta[0])], d.Levels, d.Variances, d.Coefficients))];
        CutSmcKernel kernel = Prepared(model, stage, population);
        WalkAndVariance(multivariate ? "multivariate" : "no discrepancy, σ²_e·I", kernel, stage, population[^1]);
        CalibrationState x0 = population[^1];
        SmcMove<CalibrationState> move = kernel.Move(x0, stage.Target.Evaluate(x0), stage, new Xoshiro256PlusPlus(3));
        Assert.Empty(move.State.Coefficients.ToArray());
    }

    [Theory]
    [InlineData(0.3)]
    [InlineData(0.7)]
    public void Levels_AboveTheCap_SingleSiteStepsReachTheBridgedEnumeration(double t)
    {
        CalibrationModel capped = CalibrationFixture.Model(categoricals: 2, cap: 8);
        CalibrationTarget a = capped.Target(0), b = capped.Target(2);
        var bridge = new CalibrationBridge(a, b);
        SmcStage<CalibrationState> stage = Stage(bridge, t);
        var population = new[] { new CalibrationState([0.4], [0, 0], [0.01, 0.02], [0.1, 0.05]), new CalibrationState([0.6], [1, 1], [0.012, 0.018], [0.1, 0.05]) };
        CutSmcKernel kernel = Prepared(capped, stage, population);
        Assert.False(kernel.Enumerated);
        var state = new CalibrationState([0.55], [0, 0], [0.01, 0.02], [0.1, 0.05]);
        int[][] combinations = [.. Enumerable.Range(0, 12).Select(c => new[] { c / 4, c % 4 })];
        double[] exact = CategoricalBlock.Enumerate(capped, state, kernel.Terms(stage), combinations).Probabilities;
        const int n = 60000;
        double[] chain = Run(kernel, stage, state, n, 61, (ref CutSmcKernel.Sweep s, Xoshiro256PlusPlus r) => { kernel.Levels(ref s, [], r); return true; }, s => s.Levels[0] * 4 + s.Levels[1]);
        CutSmcBlockRates rates = kernel.StageRates[^1];
        _out.WriteLine($"t = {t}, single site: acceptance {rates.Sites:F2} over {rates.SiteAttempts} steps");
        for (int c = 0; c < 12; c++)
        {
            double[] at = [.. chain.Select(x => x == c ? 1.0 : 0.0)];
            Assert.InRange(at.Average() - exact[c], -5 * CalibrationFixture.Mcse(at) - 0.002, 5 * CalibrationFixture.Mcse(at) + 0.002);
        }
    }

    [Fact]
    public void Kernel_RefusesStagesItDoesNotMove()
    {
        CalibrationModel model = Model(), other = Model();
        var kernel = new CutSmcKernel(model);
        CalibrationState[] population = [model.DrawPrior(new Xoshiro256PlusPlus(1))];
        Assert.Throws<InvalidOperationException>(() => kernel.Prepare(new SmcStage<CalibrationState>(model.Target(0), new Temperature(0.5, TemperingPath.Powered)), population, [1.0]));
        Assert.Throws<InvalidOperationException>(() => kernel.Prepare(Stage(other.Target(0), 0.5), population, [1.0]));
        Assert.Throws<InvalidOperationException>(() => kernel.Prepare(Stage(SamplingTarget.FromLogDensity<CalibrationState>(_ => 0.0), 0.5), population, [1.0]));
        // A bridge's ends and a bridge from a draw to itself are single conditionals.
        Assert.Equal([(0, 1.0)], kernel.Terms(Stage(new CalibrationBridge(model.Target(0), model.Target(2)), 0.0)).Select(x => (x.Target.Draw, x.Weight)));
        Assert.Equal([(2, 1.0)], kernel.Terms(Stage(new CalibrationBridge(model.Target(0), model.Target(2)), 1.0)).Select(x => (x.Target.Draw, x.Weight)));
        CalibrationTarget same = model.Target(1);
        Assert.Equal([(1, 1.0)], kernel.Terms(Stage(new CalibrationBridge(same, same), 0.4)).Select(x => (x.Target.Draw, x.Weight)));
        Assert.Empty(kernel.Terms(Stage(model.Target(1), 0.0)));
    }

    private delegate bool BlockStep(ref CutSmcKernel.Sweep sweep, Xoshiro256PlusPlus rng);

    // One block repeated from `start`, the sweep's state handed on, one reading per step.
    private static double[] Run(CutSmcKernel kernel, SmcStage<CalibrationState> stage, CalibrationState start, int n, int seed,
                                BlockStep step, Func<CalibrationState, double> read)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var sweep = new CutSmcKernel.Sweep(start, stage.Target.Evaluate(start), stage);
        var values = new double[n];
        for (int i = 0; i < n; i++)
        {
            step(ref sweep, rng);
            if (i % 97 == 0)
                Assert.Equal(stage.LogDensity(stage.Target.Evaluate(new CalibrationState(sweep.State.Unit.ToArray(), sweep.State.Levels.ToArray(),
                sweep.State.Variances.ToArray(), sweep.State.Coefficients.ToArray()))), sweep.Density, 1e-9);
            values[i] = read(sweep.State);
        }
        return values;
    }

    private static (double Mean, double Median) Moments(double[] grid, double[] logDensity)
    {
        double max = logDensity.Max();
        double[] w = [.. logDensity.Select(l => Math.Exp(l - max))];
        double total = w.Sum(), mean = grid.Zip(w, (x, v) => x * v).Sum() / total, cumulative = 0;
        for (int i = 0; i < grid.Length; i++)
            if ((cumulative += w[i] / total) >= 0.5) return (mean, grid[i]);
        return (mean, grid[^1]);
    }

    private void Check(string label, double[] chain, double mean, double median, double resolution)
    {
        double[] below = [.. chain.Select(x => x <= median ? 1.0 : 0.0)];
        _out.WriteLine($"{label}: mean {chain.Average():F4} vs {mean:F4} (MCSE {CalibrationFixture.Mcse(chain):F4}); P(below the median) {below.Average():F3} (MCSE {CalibrationFixture.Mcse(below):F3})");
        Assert.InRange(chain.Average() - mean, -5 * CalibrationFixture.Mcse(chain) - resolution, 5 * CalibrationFixture.Mcse(chain) + resolution);
        Assert.InRange(below.Average() - 0.5, -5 * CalibrationFixture.Mcse(below) - 0.003, 5 * CalibrationFixture.Mcse(below) + 0.003);
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
}
