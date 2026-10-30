using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass.Calibration;
using Maths.Rng;
using Maths.Samplers.Smc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The cut SMC (stage 8, C5; D09) against quadrature and stage 7's nested reference (C3): the pooled cut posterior
/// against the exact equal mixture of its conditionals, and each population against its own; the pooled means and
/// variances of θ, the levels, σ² and v against C3's at a fixed inner length, over paired replicates (the same outer
/// draws), with a discrepancy, two groups and a level, on BASS's basis and direct fits, under the multivariate
/// formulation and above the enumeration cap, and each μ_s against C3's inner run on the same e_s; one draw (the first
/// leg alone, resampled at every stage); the cut kept against a draw of far lower evidence; a bridge whose conditional
/// ESS collapses; weighted draws; determinism in the seed whatever the parallelism; the config's round trip and
/// refusals; and a bridge's first step.
/// </summary>
public sealed class BassCalibrationCutSmcTests
{
    private readonly ITestOutputHelper _out;
    public BassCalibrationCutSmcTests(ITestOutputHelper output) => _out = output;

    // Stage 7's pooled-posterior fixture: one θ, a categorical of 3 levels, one inverse-gamma group, the discrepancy, three draws.
    private static CalibrationModel Model(int seed, CutSmcSpec smc, int categoricals = 1)
    {
        TestEmulator emulator = CalibrationFixture.Emulator(categoricals);
        CalibrationConfig config = CalibrationFixture.Config(emulator, CalibrationFixture.Groups(one: true)) with
        {
            Sampling = new CalibrationSampling { MasterSeed = seed },
            Smc = smc,
        };
        return new CalibrationModel(config, CalibrationFixture.Observations(emulator), emulator);
    }

    /// <summary>A pooled posterior's means and variances: each θ in physical coordinates, an indicator per level of each
    /// categorical parameter but the last, each group's σ², each discrepancy coefficient v.</summary>
    internal sealed record Summary(string[] Names, double[] Mean, double[] Variance);

    internal static Summary Summarize(CalibrationModel model, IReadOnlyList<CalibrationDraw> draws, IReadOnlyList<double> weights)
    {
        var names = new List<string>();
        var columns = new List<Func<CalibrationDraw, double>>();
        ParameterLayout layout = model.Layout;
        for (int j = 0; j < layout.Continuous.Count; j++) { int k = j; names.Add(layout.Continuous[j].Name); columns.Add(d => d.Theta[k]); }
        for (int j = 0; j < layout.Categorical.Count; j++)
            for (int l = 0; l + 1 < layout.Categorical[j].Levels.Length; l++) { int k = j, c = l; names.Add($"{layout.Categorical[j].Name}={c}"); columns.Add(d => d.Levels[k] == c ? 1 : 0); }
        for (int g = 0; g < model.Measurement.Groups; g++) { int k = g; names.Add($"σ²[{g}]"); columns.Add(d => d.Variances[k]); }
        for (int q = 0; q < model.Coefficients; q++) { int k = q; names.Add($"v[{q}]"); columns.Add(d => d.Coefficients[k]); }
        double total = weights.Sum();
        var mean = new double[columns.Count];
        var variance = new double[columns.Count];
        for (int c = 0; c < columns.Count; c++)
        {
            double m = 0;
            for (int i = 0; i < draws.Count; i++) m += weights[i] * columns[c](draws[i]);
            m /= total;
            double v = 0;
            for (int i = 0; i < draws.Count; i++) v += weights[i] * Math.Pow(columns[c](draws[i]) - m, 2);
            (mean[c], variance[c]) = (m, v / total);
        }
        return new Summary([.. names], mean, variance);
    }

    /// <summary>
    /// Paired replicates of two estimators (the same seed, so the same outer draws): for each mean and variance, the
    /// replicates' differences must average within <paramref name="sigmas"/> standard errors of 0, the error read from the
    /// differences' own spread over the replicates.
    /// </summary>
    private void Agree(string label, IReadOnlyList<Summary> a, IReadOnlyList<Summary> b, double sigmas = 6)
    {
        int r = a.Count;
        string[] names = a[0].Names;
        var lines = new List<string>();
        foreach ((string moment, Func<Summary, double[]> read) in new (string, Func<Summary, double[]>)[] { ("mean", x => x.Mean), ("var", x => x.Variance) })
            for (int q = 0; q < names.Length; q++)
            {
                double[] diff = [.. Enumerable.Range(0, r).Select(k => read(a[k])[q] - read(b[k])[q])];
                double mean = diff.Average(), se = Math.Sqrt(diff.Sum(x => (x - mean) * (x - mean)) / (r - 1) / r);
                double level = Enumerable.Range(0, r).Average(k => read(b[k])[q]);
                lines.Add($"  {moment} {names[q]}: {Enumerable.Range(0, r).Average(k => read(a[k])[q]):G5} vs {level:G5}, difference {mean:G3} (SE {se:G2}, {(se > 0 ? mean / se : 0):F1} SE)");
                Assert.True(Math.Abs(mean) <= sigmas * se + 1e-9 * (1 + Math.Abs(level)), $"{label}: {moment} of {names[q]} differs by {mean:G3}, {mean / se:F1} SE");
            }
        _out.WriteLine($"{label}, {r} paired replicates:");
        foreach (string line in lines) _out.WriteLine(line);
    }

    // C3 at a fixed inner length (an unreachable ESS gate, every run pooled), on the model's own seed and outer draws.
    private static (IReadOnlyList<CalibrationDraw> Draws, double[] Weights, NestedCalibration Result) Reference(CalibrationModel model, int samples)
    {
        CalibrationConfig config = model.Config with
        {
            Reference = model.Config.Reference with { Chains = 2, MaxSamples = samples, BatchSize = samples, MinEssPerChain = 1e9, OuterDraws = model.Config.Smc?.OuterDraws },
        };
        NestedCalibration result = NestedCalibrationSampler.Run(new CalibrationModel(config, model.Observations.ToArray(), model.Emulator));
        (IReadOnlyList<CalibrationDraw> draws, double[] weights) = result.Pooled(includeNotConverged: true);
        return (draws, weights, result);
    }

    // What every run must account for: weights summing to 1, every visit's mass 1/S, each population on its draw, the
    // order a permutation from the first leg's draw, the families reported, the stages adding up.
    private static void Account(CutSmcCalibration result)
    {
        (IReadOnlyList<CalibrationDraw> draws, double[] weights) = result.Pooled();
        int s = result.Populations.Count, offset = 0;
        Assert.Equal(result.OuterDraws.Count, s);
        Assert.Equal(1.0, weights.Sum(), 1e-12);
        Assert.Equal(0, result.Order[0]);
        Assert.Equal(Enumerable.Range(0, s), result.Order.Order());
        Assert.Equal(result.Order.Select(k => result.OuterDraws[k]), result.Populations.Select(p => p.EmulatorDraw));
        foreach (CutSmcPopulation p in result.Populations)
        {
            Assert.Equal(1.0, p.Weights.Sum(), 1e-12);
            Assert.Equal(1.0 / s, weights.Skip(offset).Take(p.Draws.Count).Sum(), 1e-12);
            Assert.All(p.Draws, d => Assert.Equal(p.EmulatorDraw, d.EmulatorDraw));
            Assert.InRange(p.Leg.Families, 1, p.Draws.Count);
            Assert.True(p.Leg.StageCount >= 1);
            Assert.Equal(p.EmulatorDraw, p.Leg.To);
            offset += p.Draws.Count;
        }
        Assert.Null(result.Populations[0].Leg.From);
        Assert.Equal(result.TotalStages, result.Populations.Sum(p => p.Leg.StageCount));
        Assert.Equal(result.Families, result.Populations[^1].Leg.Families);
    }

    private static string Describe(CutSmcCalibration result)
        => $"order {string.Join(",", result.Populations.Select(p => p.EmulatorDraw))}, tour {result.TourLength:F2}, {result.TotalStages} stages ("
           + string.Join(" + ", result.Populations.Select(p => p.Leg.StageCount)) + $"), least ESS {result.Populations.Min(p => p.Leg.MinEss):F0}, families {result.Families}";

    private static double UnitMean(IReadOnlyList<CalibrationDraw> draws, IReadOnlyList<double> weights)
        => Enumerable.Range(0, draws.Count).Sum(i => weights[i] * CalibrationFixture.Theta.ToUnit(draws[i].Theta[0])) / weights.Sum();

    [Fact]
    public void Pooled_IsTheExactEqualMixture_ByQuadrature()
    {
        // Stage 7's pooled-posterior fixture, each conditional by quadrature with v integrated out: the pooled E[θ],
        // E[ln σ²] and level probabilities against their equal mixture, and each population against its own conditional.
        CalibrationModel probe = Model(0, new CutSmcSpec());
        var exact = Enumerable.Range(0, 3).Select(s => BassCalibrationSamplerTests.Quadrature(probe, probe.Target(s), 1.0)).ToArray();
        const int r = 6;
        var pooledTheta = new double[r];
        var pooledU = new double[r];
        var pooledLevels = new double[r][];
        var byDraw = new double[3][];
        var logZ = new double[3][];
        for (int s = 0; s < 3; s++) (byDraw[s], logZ[s]) = (new double[r], new double[r]);
        for (int rep = 0; rep < r; rep++)
        {
            CalibrationModel model = Model(100 + rep, new CutSmcSpec { Particles = 500 });
            CutSmcCalibration result = CutSmcSampler.Run(model);
            Account(result);
            _out.WriteLine($"replicate {rep}: {Describe(result)}");
            (IReadOnlyList<CalibrationDraw> draws, double[] weights) = result.Pooled();
            pooledTheta[rep] = UnitMean(draws, weights);
            pooledU[rep] = Enumerable.Range(0, draws.Count).Sum(i => weights[i] * Math.Log(draws[i].Variances[0]));
            pooledLevels[rep] = [.. Enumerable.Range(0, 3).Select(c => Enumerable.Range(0, draws.Count).Where(i => draws[i].Levels[0] == c).Sum(i => weights[i]))];
            foreach (CutSmcPopulation p in result.Populations) (byDraw[p.EmulatorDraw][rep], logZ[p.EmulatorDraw][rep]) = (UnitMean(p.Draws, p.Weights), p.LogEvidence);
        }
        const double grid = 2.0 / 600;
        void Near(string label, double[] values, double target, double resolution)
        {
            double mean = values.Average(), se = Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (r - 1) / r);
            _out.WriteLine($"{label}: {mean:F4} vs {target:F4} (SE {se:F4} over {r} replicates, {(mean - target) / se:F1} SE)");
            Assert.InRange(mean - target, -6 * se - resolution, 6 * se + resolution);
        }
        Near("pooled E[θ]", pooledTheta, exact.Average(q => q.MeanTheta), grid);
        Near("pooled E[ln σ²]", pooledU, exact.Average(q => q.MeanLogVariance), 0.02);
        for (int c = 0; c < 3; c++) Near($"pooled P(level {c})", [.. pooledLevels.Select(l => l[c])], exact.Average(q => q.Levels[c]), 0.005);
        for (int s = 0; s < 3; s++) Near($"μ_{s} E[θ]", byDraw[s], exact[s].MeanTheta, grid);
        // The run's normalizer after each leg is log p(y | e_s): the first leg's evidence, then each bridge's ratio, so it
        // checks every incremental weight on the way.
        for (int s = 0; s < 3; s++) Near($"log p(y | e_{s})", logZ[s], LogEvidence(probe, probe.Target(s)), 0.01);
    }

    /// <summary>log p(y | e) by quadrature over (θ, ln σ², level), v integrated out in closed form, N(y; prediction + μ_δ,
    /// Σ + s_v²DDᵀ), with every prior's normalization: θ uniform on the unit box, the level 1/L, σ² its group's.</summary>
    internal static double LogEvidence(CalibrationModel model, CalibrationTarget target)
    {
        int m = CalibrationFixture.M, levels = model.Layout.Categorical.Count == 0 ? 1 : model.Layout.Categorical[0].Levels.Length;
        const int thetaPoints = 600, uPoints = 400;
        const double uLower = -14, uUpper = 0;
        double hu = (uUpper - uLower) / uPoints, ht = 1.0 / thetaPoints;
        double[,] d = model.Discrepancy!;
        double s2 = model.DiscrepancyScale * model.DiscrepancyScale;
        double[,] fixedPart = BassCalibrationBridgeTests.FixedSigma(model, target);
        var factors = new DenseCholesky[uPoints];
        var baseWeight = new double[uPoints];
        for (int iu = 0; iu < uPoints; iu++)
        {
            double u = uLower + (iu + 0.5) * hu;
            var sigma = new double[m, m];
            for (int j = 0; j < m; j++)
                for (int l = 0; l < m; l++)
                    sigma[j, l] = fixedPart[j, l] + (j == l ? Math.Exp(u) : 0) + s2 * (d[j, 0] * d[l, 0] + d[j, 1] * d[l, 1]);
            factors[iu] = new DenseCholesky(m);
            factors[iu].Decompose(sigma);
            baseWeight[iu] = model.Measurement.LogPrior(0, Math.Exp(u)) + u + Math.Log(hu) + Math.Log(ht) - Math.Log(levels)
                             - 0.5 * (m * Math.Log(2 * Math.PI) + factors[iu].LogDet);
        }
        var terms = new List<double>(levels * thetaPoints * uPoints);
        var pred = new double[m];
        double[] mean = model.Covariance.Mean ?? new double[m];
        for (int c = 0; c < levels; c++)
            for (int it = 0; it < thetaPoints; it++)
            {
                model.Emulator.Predict(target.Draw, [CalibrationFixture.Theta.ToPhysical((it + 0.5) * ht)], levels > 1 ? [c] : [], pred);
                double[] r = [.. Enumerable.Range(0, m).Select(j => model.Observations[j] - pred[j] - mean[j])];
                for (int iu = 0; iu < uPoints; iu++) terms.Add(baseWeight[iu] - 0.5 * factors[iu].InverseQuadraticForm(r));
            }
        double max = terms.Max();
        return max + Math.Log(terms.Sum(t => Math.Exp(t - max)));
    }

    // Runs `r` paired replicates on models from `model(k)`: the cut SMC in turn on every core, the nested reference at a
    // fixed length concurrently.
    private (CutSmcCalibration[] Smc, (IReadOnlyList<CalibrationDraw> Draws, double[] Weights, NestedCalibration Result)[] Reference) Replicates(
        int r, Func<int, CalibrationModel> model, int samples)
    {
        CalibrationModel[] models = [.. Enumerable.Range(0, r).Select(model)];
        var reference = new (IReadOnlyList<CalibrationDraw>, double[], NestedCalibration)[r];
        Task references = Task.Run(() => Parallel.For(0, r, k => reference[k] = Reference(models[k], samples)));
        var smc = new CutSmcCalibration[r];
        for (int k = 0; k < r; k++)
        {
            smc[k] = CutSmcSampler.Run(models[k]);
            Account(smc[k]);
            _out.WriteLine($"replicate {k}: {Describe(smc[k])}");
        }
        references.Wait();
        return (smc, reference);
    }

    [Fact]
    public void AgreesWithTheNestedReference_WithADiscrepancy_TwoGroups_AndALevel()
    {
        // Stage 7's conditional fixture whole: an inverse-gamma and a half-Cauchy group, the discrepancy on two columns,
        // a categorical of 3 levels, three draws. Pooled, then each μ_s against C3's inner run on the same e_s.
        TestEmulator emulator = CalibrationFixture.Emulator(1);
        double[] y = CalibrationFixture.Observations(emulator);
        const int r = 6;
        var (smc, reference) = Replicates(r, k => new CalibrationModel(CalibrationFixture.Config(emulator) with
        {
            Sampling = new CalibrationSampling { MasterSeed = 200 + k, BurnIn = 500, DecorrelationEvery = 25 },
            Smc = new CutSmcSpec { Particles = 400 },
        }, y, emulator), 1000);
        CalibrationModel model = new(CalibrationFixture.Config(emulator), y, emulator);
        Agree("pooled", [.. smc.Select(x => { var (d, w) = x.Pooled(); return Summarize(model, d, w); })],
              [.. reference.Select(x => Summarize(model, x.Draws, x.Weights))]);
        for (int s = 0; s < 3; s++)
        {
            int draw = s;
            Summary[] a = [.. smc.Select(x => x.Populations.Single(p => p.EmulatorDraw == draw)).Select(p => Summarize(model, p.Draws, p.Weights))];
            Summary[] b = [.. reference.Select(x => x.Result.Runs.Single(run => run.EmulatorDraw == draw)).Select(run =>
            {
                CalibrationDraw[] d = [.. run.Chains.SelectMany(c => c)];
                return Summarize(model, d, [.. Enumerable.Repeat(1.0, d.Length)]);
            })];
            // θ's mean and spread asserted, the rest read.
            Agree($"μ_{s}: θ", [.. a.Select(Theta)], [.. b.Select(Theta)]);
        }
    }

    private static Summary Theta(Summary x) => new([x.Names[0]], [x.Mean[0]], [x.Variance[0]]);

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void AgreesWithTheNestedReference_OnBassFits(bool basis)
    {
        // Stage 7's recovery fixture: BASS fitted on a PCA basis, or directly; field data at θ* = (1.3, 0.4), six outer
        // draws of 200.
        ContinuousParameter[] parameters = [BassCalibrationRecoveryTests.Theta1, BassCalibrationRecoveryTests.Theta2];
        ICalibrationEmulator emulator = basis
            ? new BassBasisEmulator(BassCalibrationRecoveryTests.BasisFit.Value, parameters)
            : new BassFitEmulator(BassCalibrationRecoveryTests.DirectFit.Value.Fit, BassCalibrationRecoveryTests.DirectFit.Value.Data, parameters);
        double[] y = BassCalibrationRecoveryTests.Field(1, [.. Enumerable.Repeat(0.05, BassCalibrationRecoveryTests.M)]);
        CalibrationConfig Config(int seed) => new()
        {
            Parameters = parameters,
            Measurement = BassCalibrationRecoveryTests.OneGroup(),
            Covariance = new FunctionalFormulation(),
            Sampling = new CalibrationSampling { MasterSeed = seed, BurnIn = 2000, DecorrelationEvery = 50 },
            Smc = new CutSmcSpec { OuterDraws = 6 },   // the default population: at 400 particles the basis fit's var θ₁ ran 3% low
        };
        // C3 at 4,000 draws a chain: at 1,000 its direct-fit E[θ₁] sat 2.8 SE off its own 8,000-draw value.
        const int r = 6;
        var (smc, reference) = Replicates(r, k => new CalibrationModel(Config(300 + k), y, emulator), 4000);
        for (int k = 0; k < r; k++)
            Assert.Equal(reference[k].Result.Runs.Select(run => run.EmulatorDraw).Order(), smc[k].OuterDraws.Order());   // the same outer draws
        var model = new CalibrationModel(Config(0), y, emulator);
        Agree(basis ? "basis" : "direct", [.. smc.Select(x => { var (d, w) = x.Pooled(); return Summarize(model, d, w); })],
              [.. reference.Select(x => Summarize(model, x.Draws, x.Weights))]);
    }

    [Theory]
    [InlineData(true)]
    [InlineData(false)]
    public void AgreesWithTheNestedReference_MultivariateOrAboveTheCap(bool multivariate)
    {
        // The multivariate formulation (no v to draw), or two categorical parameters above an enumeration cap of 8
        // (single-site level steps, in both samplers).
        TestEmulator emulator = CalibrationFixture.Emulator(multivariate ? 0 : 2);
        double[] y = CalibrationFixture.Observations(emulator);
        int m = CalibrationFixture.M;
        CovarianceFormulation covariance = multivariate
            ? new MultivariateFormulation(
                [.. Enumerable.Range(0, m).Select(i => Enumerable.Range(0, m).Select(j => Math.Pow(0.4, Math.Abs(i - j))).ToArray())],
                [.. Enumerable.Range(0, m).Select(i => Enumerable.Range(0, m).Select(j => (i == j ? 1e-4 : 0.0) + 0.002 * Math.Exp(-0.5 * (i - j) * (i - j))).ToArray())],
                [.. Enumerable.Range(0, m).Select(j => 0.01 * j)])
            : new FunctionalFormulation { Discrepancy = CalibrationFixture.Discrepancy() };
        CalibrationConfig Config(int seed) => CalibrationFixture.Config(emulator, CalibrationFixture.Groups(one: true), covariance, cap: 8) with
        {
            Sampling = new CalibrationSampling { MasterSeed = seed, BurnIn = 500, DecorrelationEvery = 25 },
            Smc = multivariate ? new CutSmcSpec() : new CutSmcSpec { Particles = 300 },
        };
        const int r = 6;
        var (smc, reference) = Replicates(r, k => new CalibrationModel(Config(400 + k), y, emulator), 1500);
        var model = new CalibrationModel(Config(0), y, emulator);
        Assert.Equal(!multivariate, smc.All(x => x.Populations.All(p => p.Leg.Blocks.SiteAttempts > 0)));
        Agree(multivariate ? "multivariate" : "above the cap", [.. smc.Select(x => { var (d, w) = x.Pooled(); return Summarize(model, d, w); })],
              [.. reference.Select(x => Summarize(model, x.Draws, x.Weights))]);
    }

    [Fact]
    public void OneDraw_IsTheFirstLegAlone()
    {
        // S = 1: the first leg only, resampling at every stage (every offspring after the first a fork).
        CalibrationModel model = Model(5, new CutSmcSpec { Particles = 600, OuterDraws = 1, ResampleBelow = 1.0 });
        CutSmcCalibration result = CutSmcSampler.Run(model);
        Account(result);
        Assert.Single(result.Populations);
        Assert.Null(result.Distances);
        Assert.True(double.IsNaN(result.TourLength));
        CutSmcPopulation p = result.Populations[0];
        Assert.True(p.Leg.Stages.All(stage => stage.Resampled));
        Assert.True(p.Leg.Families < 600);
        var q = BassCalibrationSamplerTests.Quadrature(model, model.Target(p.EmulatorDraw), 1.0);
        double mean = UnitMean(p.Draws, p.Weights);
        double sd = Math.Sqrt(Enumerable.Range(0, p.Draws.Count).Sum(i => p.Weights[i] * Math.Pow(CalibrationFixture.Theta.ToUnit(p.Draws[i].Theta[0]) - mean, 2)));
        // Every lineage but the families' is a copy: the error is at most the families' independent draws'.
        double se = sd / Math.Sqrt(p.Leg.Families);
        _out.WriteLine($"S = 1: {Describe(result)}; E[θ] {mean:F4} vs {q.MeanTheta:F4} (sd {sd:F3}, SE ≤ {se:F4}); log Z {p.LogEvidence:F3} ± {p.LogEvidenceError:F3}");
        Assert.InRange(mean - q.MeanTheta, -6 * se - 2.0 / 600, 6 * se + 2.0 / 600);
        Assert.True(double.IsFinite(p.LogEvidenceError));
    }

    /// <summary>An emulator over a <see cref="TestEmulator"/> whose draw 1 predicts as draw 0 at θ moved by
    /// <paramref name="shift"/> (physical), and whose draw 2 adds <paramref name="misfit"/>·cos(πj) to output j; its draws
    /// weighted when <paramref name="weights"/> is given.</summary>
    private sealed class AlteredEmulator(TestEmulator inner, double shift = 0.0, double misfit = 0.0, double[]? weights = null) : ICalibrationEmulator
    {
        public IReadOnlyList<CalibrationParameter> Parameters => inner.Parameters;
        public int DrawCount => inner.DrawCount;
        public IReadOnlyList<double>? DrawWeights => weights;
        public int Outputs => inner.Outputs;
        public EmulatorNoise Noise(int draw) => inner.Noise(draw);
        public double[,]? TruncationResiduals() => inner.TruncationResiduals();

        public void Predict(int draw, ReadOnlySpan<double> continuous, ReadOnlySpan<int> levels, Span<double> output)
        {
            if (draw == 1 && shift != 0.0) inner.Predict(0, [continuous[0] + shift], levels, output);
            else inner.Predict(draw, continuous, levels, output);
            if (draw == 2) for (int j = 0; j < output.Length; j++) output[j] += misfit * Math.Cos(Math.PI * j);
        }
    }

    private static CalibrationModel Altered(AlteredEmulator emulator, TestEmulator inner, int seed, CutSmcSpec smc)
        => new(CalibrationFixture.Config(inner, CalibrationFixture.Groups(one: true)) with { Sampling = new CalibrationSampling { MasterSeed = seed }, Smc = smc },
               CalibrationFixture.Observations(inner), emulator);

    [Fact]
    public void TheCutIsKept_NoStageReweightsTheEmulator()
    {
        // Draw 2 misfits the data by ±0.15 on alternate outputs, out of the discrepancy's span: its evidence falls well
        // below the others', so a joint posterior would all but drop it, and the cut pools it at 1/S all the same.
        TestEmulator inner = CalibrationFixture.Emulator(1);
        var emulator = new AlteredEmulator(inner, misfit: 0.15);
        CutSmcCalibration result = CutSmcSampler.Run(Altered(emulator, inner, 7, new CutSmcSpec { Particles = 400 }));
        Account(result);
        double[] logZ = [.. Enumerable.Range(0, 3).Select(d => result.Populations.Single(p => p.EmulatorDraw == d).LogEvidence)];
        _out.WriteLine($"{Describe(result)}; log p(y | e_s): {string.Join(", ", logZ.Select(z => z.ToString("F2")))}");
        double joint = 1 / (1 + Math.Exp(logZ[0] - logZ[2]) + Math.Exp(logZ[1] - logZ[2]));
        _out.WriteLine($"draw 2's weight: {joint:P2} in a joint posterior, 1/3 in the cut");
        Assert.True(joint < 0.02);
        (IReadOnlyList<CalibrationDraw> draws, double[] weights) = result.Pooled();
        for (int d = 0; d < 3; d++) Assert.Equal(1.0 / 3, Enumerable.Range(0, draws.Count).Where(i => draws[i].EmulatorDraw == d).Sum(i => weights[i]), 1e-12);
    }

    [Fact]
    public void ABridgeWhoseConditionalEssCollapses_TakesSmallStepsAndArrives()
    {
        // Draw 1's conditional sits 0.8 lower in θ (about five posterior sds): its bridge's increments spread widely, the
        // conditional ESS of a whole step collapses, and the bridge crosses in many small steps, each at the target.
        TestEmulator inner = CalibrationFixture.Emulator(1);
        var emulator = new AlteredEmulator(inner, shift: 0.8);
        var spec = new CutSmcSpec { Particles = 800 };
        CalibrationModel model = Altered(emulator, inner, 9, spec);
        CutSmcCalibration result = CutSmcSampler.Run(model);
        Account(result);
        _out.WriteLine(Describe(result));
        foreach (CutSmcPopulation p in result.Populations)
            _out.WriteLine($"  to e_{p.EmulatorDraw}: {p.Leg.StageCount} stages, first t {p.Leg.FirstStep:G3}, d {p.Leg.Distance:F2}, least ESS {p.Leg.MinEss:F0}, " +
                           $"least conditional ESS {p.Leg.MinConditionalEss:F1}, families {p.Leg.Families}, walk {p.Leg.Blocks.Theta:F2}, resampled {p.Leg.Stages.Count(x => x.Resampled)}");
        CutSmcPopulation far = result.Populations.Single(p => p.EmulatorDraw == 1);
        CutSmcPopulation near = result.Populations.Skip(1).Single(p => p.EmulatorDraw != 1);
        Assert.True(far.Leg.StageCount > 3 * near.Leg.StageCount);
        Assert.True(far.Leg.FirstStep < 0.1);
        foreach (CutSmcPopulation p in result.Populations.Skip(1))
            Assert.All(p.Leg.Stages, x => Assert.True(x.ConditionalEss >= spec.TargetConditionalEss * spec.Particles * (1 - 1e-6)));
        // The far draw's population is its conditional.
        var q = BassCalibrationSamplerTests.Quadrature(model, model.Target(1), 1.0);
        double mean = UnitMean(far.Draws, far.Weights);
        double sd = Math.Sqrt(Enumerable.Range(0, far.Draws.Count).Sum(i => far.Weights[i] * Math.Pow(CalibrationFixture.Theta.ToUnit(far.Draws[i].Theta[0]) - mean, 2)));
        double se = sd / Math.Sqrt(Math.Min(far.Leg.Families, 1 / far.Weights.Sum(w => w * w)));
        _out.WriteLine($"  e_1: E[θ] {mean:F4} vs {q.MeanTheta:F4} (sd {sd:F3}, SE ≤ {se:F4})");
        Assert.InRange(mean - q.MeanTheta, -6 * se - 2.0 / 600, 6 * se + 2.0 / 600);
        // Too few stages for the crossing: the run stops, naming the leg.
        CalibrationModel tight = Altered(emulator, inner, 9, spec with { MaxStagesPerLeg = 8 });
        Assert.Contains("MaxStagesPerLeg", Assert.Throws<InvalidOperationException>(() => CutSmcSampler.Run(tight)).Message);
    }

    [Fact]
    public void WeightedDraws_TheResampledSetIsToured_AndPooledAtOneOverSPerVisit()
    {
        // Weights 0.2, 0.5, 0.3 on the three draws, eight outer draws resampled by weight: the pooled E[θ] against the
        // weighted mixture of the conditionals by quadrature, each draw's pooled mass its share of the visits.
        TestEmulator inner = CalibrationFixture.Emulator(1);
        double[] w = [0.2, 0.5, 0.3];
        var emulator = new AlteredEmulator(inner, weights: w);
        CalibrationModel probe = Altered(emulator, inner, 0, new CutSmcSpec());
        double exact = Enumerable.Range(0, 3).Sum(s => w[s] * BassCalibrationSamplerTests.Quadrature(probe, probe.Target(s), 1.0).MeanTheta);
        const int r = 6;
        var means = new double[r];
        for (int rep = 0; rep < r; rep++)
        {
            CutSmcCalibration result = CutSmcSampler.Run(Altered(emulator, inner, 500 + rep, new CutSmcSpec { Particles = 300, OuterDraws = 8 }));
            Account(result);
            Assert.Equal(NestedCalibrationSampler.OuterDraws(emulator, 8, new Xoshiro256PlusPlus(SeedTree.Derive(500 + rep, 3)[0])), result.OuterDraws);
            int[] visited = [.. result.Populations.Select(p => p.EmulatorDraw)];
            for (int d = 0; d < 3; d++)
            {
                int[] at = [.. Enumerable.Range(0, 8).Where(k => visited[k] == d)];
                if (at.Length > 0) Assert.Equal(at.Length - 1, at[^1] - at[0]);   // a draw's copies visited in one run
            }
            (IReadOnlyList<CalibrationDraw> draws, double[] weights) = result.Pooled();
            for (int d = 0; d < 3; d++)
                Assert.Equal(visited.Count(v => v == d) / 8.0, Enumerable.Range(0, draws.Count).Where(i => draws[i].EmulatorDraw == d).Sum(i => weights[i]), 1e-12);
            means[rep] = UnitMean(draws, weights);
            _out.WriteLine($"replicate {rep}: visits {string.Join(",", visited)}; {Describe(result)}; E[θ] {means[rep]:F4}");
        }
        double mean = means.Average(), se = Math.Sqrt(means.Sum(v => (v - mean) * (v - mean)) / (r - 1) / r);
        _out.WriteLine($"pooled E[θ] {mean:F4} vs the weighted mixture {exact:F4} (SE {se:F4} over {r} replicates, outer resampling included)");
        Assert.InRange(mean - exact, -6 * se - 2.0 / 600, 6 * se + 2.0 / 600);
    }

    [Fact]
    public void Run_IsDeterministicInTheSeed_WhateverTheParallelism()
    {
        CalibrationModel model = Model(13, new CutSmcSpec { Particles = 200 });
        CutSmcCalibration a = CutSmcSampler.Run(model, 1), b = CutSmcSampler.Run(model, 6);
        CutSmcCalibration c = CutSmcSampler.Run(Model(14, new CutSmcSpec { Particles = 200 }), 6);
        (IReadOnlyList<CalibrationDraw> da, double[] wa) = a.Pooled();
        (IReadOnlyList<CalibrationDraw> db, double[] wb) = b.Pooled();
        Assert.Equal(a.Order, b.Order);
        Assert.Equal(wa, wb);
        for (int i = 0; i < da.Count; i++)
        {
            Assert.Equal(da[i].Theta, db[i].Theta);
            Assert.Equal(da[i].Variances, db[i].Variances);
            Assert.Equal(da[i].Coefficients, db[i].Coefficients);
            Assert.Equal(da[i].Levels, db[i].Levels);
            Assert.Equal(da[i].LogLikelihood, db[i].LogLikelihood);
        }
        Assert.Equal(a.Populations.Select(p => p.LogEvidence), b.Populations.Select(p => p.LogEvidence));
        Assert.NotEqual(da[0].Theta, c.Pooled().Draws[0].Theta);
        // A random order: the same first leg, a permutation after it, no distances.
        CutSmcCalibration random = CutSmcSampler.Run(Model(13, new CutSmcSpec { Particles = 200, Order = CutSmcOrder.Random }), 6);
        Account(random);
        Assert.Null(random.Distances);
        Assert.Equal(a.Populations[0].Draws.Select(x => x.Theta[0]), random.Populations[0].Draws.Select(x => x.Theta[0]));
        Assert.Equal(random.Order, CutSmcSampler.Run(Model(13, new CutSmcSpec { Particles = 200, Order = CutSmcOrder.Random }), 2).Order);
    }

    [Fact]
    public void Config_RoundTripsThroughJson_AndRefusesBadSettings()
    {
        TestEmulator emulator = CalibrationFixture.Emulator();
        CalibrationConfig config = CalibrationFixture.Config(emulator) with
        {
            Smc = new CutSmcSpec
            {
                Particles = 700, TargetConditionalEss = 0.95, MutationSteps = 3, ResampleBelow = 0.4, Resampling = ResamplingScheme.Systematic,
                OuterDraws = 2, Order = CutSmcOrder.Random, MaxStagesPerLeg = 250,
            },
        };
        string json = config.ToJson();
        CalibrationConfig back = CalibrationConfig.FromJson(json);
        Assert.Equal(json, back.ToJson());
        Assert.Equal(config.Smc, back.Smc);
        Assert.Contains("\"Order\":\"Random\"", json);
        // Off by default: null, and a config written before stage 8 reads as null.
        CalibrationConfig plain = CalibrationFixture.Config(emulator);
        Assert.Null(CalibrationConfig.FromJson(plain.ToJson()).Smc);
        Assert.Contains(",\"Smc\":null", plain.ToJson());
        Assert.Null(CalibrationConfig.FromJson(plain.ToJson().Replace(",\"Smc\":null", "")).Smc);
        Assert.Equal(new CutSmcSpec { Particles = 1000, TargetConditionalEss = 0.99, MutationSteps = 5, ResampleBelow = 0.5, Resampling = ResamplingScheme.Multinomial }, new CutSmcSpec());
        foreach ((string from, string to, string field) in new[]
        {
            ("\"Particles\":700", "\"Particles\":1", "Particles"),
            ("\"TargetConditionalEss\":0.95", "\"TargetConditionalEss\":1", "TargetConditionalEss"),
            ("\"MutationSteps\":3", "\"MutationSteps\":0", "MutationSteps"),
            ("\"ResampleBelow\":0.4", "\"ResampleBelow\":1.5", "ResampleBelow"),
            ("\"OuterDraws\":2,\"Order\"", "\"OuterDraws\":0,\"Order\"", "OuterDraws"),
            ("\"MaxStagesPerLeg\":250", "\"MaxStagesPerLeg\":0", "MaxStagesPerLeg"),
        })
        {
            Assert.Contains(from, json);
            Assert.Contains($"CalibrationConfig.Smc.{field}", Assert.Throws<ArgumentException>(() => CalibrationConfig.FromJson(json.Replace(from, to))).Message);
        }
        Assert.ThrowsAny<Exception>(() => CalibrationConfig.FromJson(json.Replace("\"Random\"", "\"Sideways\"")));
    }

    [Fact]
    public void Run_RefusesWhatItCannotRun()
    {
        TestEmulator emulator = CalibrationFixture.Emulator();
        double[] y = CalibrationFixture.Observations(emulator);
        Assert.Contains("NestedCalibrationSampler", Assert.Throws<ArgumentException>(() =>
            CutSmcSampler.Run(new CalibrationModel(CalibrationFixture.Config(emulator), y, emulator))).Message);
        var improper = new CalibrationModel(CalibrationFixture.Config(emulator) with { Smc = new CutSmcSpec() }, y, emulator,
            new MeasurementModel([0, 0, 0, 0, 1, 1, 1, 1], [new ImproperPrior(), new InverseGammaVariancePrior(8, 0.05)]));
        Assert.Contains("prior", Assert.Throws<ArgumentException>(() => CutSmcSampler.Run(improper)).Message);
        Assert.Throws<ArgumentOutOfRangeException>(() => CutSmcSampler.Run(Model(1, new CutSmcSpec()), 0));
    }

    private sealed class ImproperPrior : IVariancePrior
    {
        public double LogDensity(double sigma2) => sigma2 > 0 ? -Math.Log(1 + sigma2) : double.NegativeInfinity;
        public double Sample(Xoshiro256PlusPlus rng) => throw new NotSupportedException();
        public bool IsProper => false;
    }

    [Fact]
    public void FirstStep_IsTheBisectionOnTheCachedIncrements()
    {
        // The largest t₁ keeping the conditional ESS at the target, as SmcSampler's bisection finds on a run.
        var rng = new Xoshiro256PlusPlus(3);
        double[] weights = [.. Enumerable.Repeat(1.0 / 500, 500)];
        double[] difference = [.. Enumerable.Range(0, 500).Select(_ => 20 * Distributions.Normal.SampleStandard(rng))];
        double t = CutSmcSampler.FirstStep(weights, difference, 0.99);
        double Cess(double s) => SmcWeights.ConditionalEss(weights, [.. difference.Select(d => s * d)]);
        Assert.InRange(t, 0, 1);
        Assert.True(Cess(t) >= 0.99 * 500 - 1e-9);
        Assert.True(Cess(t + 1e-9) < 0.99 * 500 + 1e-6);
        Assert.Equal(1.0, CutSmcSampler.FirstStep(weights, new double[500], 0.99));   // a bridge between equal conditionals crosses at once
        Assert.Equal(1.0, CutSmcSampler.FirstStep(weights, [.. difference.Select(d => 1e-4 * d)], 0.99));
    }
}
