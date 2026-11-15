using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The engines' evidences for model averaging (E16, CBMA's evidence from the engines' fits; <see cref="BassEvidence.ForAveraging"/>,
/// <see cref="BarsEvidence.ForAveraging"/>). Each engine's value is pinned against the g-prior's marginal written from its definition, a
/// dense multivariate t with every constant, −(m/2) ln 2π and ½ Σ ln wᵢ included: BARS's evidence with no knots, BASS's exact enumeration
/// unweighted and under precision weights, and a fitted BASS evidence within its error. On shared data the two engines' evidences read one
/// response, agree far inside the 2π constant either would drop, and average through CBMA.
/// </summary>
public sealed class EngineEvidenceTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    private static double LogSumExp(IEnumerable<double> values)
    {
        double[] v = [.. values];
        double top = v.Max();
        return top + Math.Log(v.Sum(t => Math.Exp(t - top)));
    }

    /// <summary>
    /// The g-prior's marginal from its definition: y | σ² ~ N(Zμ₀, σ²(W⁻¹ + Z G⁻¹ Zᵀ/τ)), G = ZᵀWZ, σ² ~ IG(a, b), so y is multivariate t:
    /// ln p = −(m/2) ln 2π − ½ ln|Σ| + a ln b + ln Γ(a + m/2) − ln Γ(a) − (a + m/2) ln(b + dᵀΣ⁻¹d/2), d = y − Zμ₀, by a dense Cholesky.
    /// </summary>
    internal static double DenseMarginal(double[][] columns, double[] y, double[]? w, double tau, double[] mean, double a, double b)
    {
        int m = y.Length, k = columns.Length;
        double W(int i) => w?[i] ?? 1.0;
        var g = new double[k, k];
        for (int p = 0; p < k; p++)
            for (int q = 0; q < k; q++)
                for (int i = 0; i < m; i++) g[p, q] += columns[p][i] * W(i) * columns[q][i];
        double[,] gi = Inverse(g);
        var sigma = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
            {
                double s = i == j ? 1.0 / W(i) : 0.0;
                for (int p = 0; p < k; p++)
                    for (int q = 0; q < k; q++) s += columns[p][i] * gi[p, q] * columns[q][j] / tau;
                sigma[i, j] = s;
            }
        // Cholesky Σ = LLᵀ: ln|Σ| = 2 Σ ln Lᵢᵢ and dᵀΣ⁻¹d = ‖L⁻¹d‖².
        var l = new double[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = sigma[i, j];
                for (int r = 0; r < j; r++) s -= l[i, r] * l[j, r];
                l[i, j] = i == j ? Math.Sqrt(s) : s / l[j, j];
            }
        double logDet = 0.0, quad = 0.0;
        var z = new double[m];
        for (int i = 0; i < m; i++)
        {
            double s = y[i] - mean[i];
            for (int r = 0; r < i; r++) s -= l[i, r] * z[r];
            z[i] = s / l[i, i];
            quad += z[i] * z[i];
            logDet += 2.0 * Math.Log(l[i, i]);
        }
        return -0.5 * m * Math.Log(2.0 * Math.PI) - 0.5 * logDet + a * Math.Log(b) + SpecialFunctions.LogGamma(a + 0.5 * m)
               - SpecialFunctions.LogGamma(a) - (a + 0.5 * m) * Math.Log(b + 0.5 * quad);
    }

    private static double[,] Inverse(double[,] a)
    {
        int k = a.GetLength(0);
        var m = new double[k, 2 * k];
        for (int i = 0; i < k; i++) { for (int j = 0; j < k; j++) m[i, j] = a[i, j]; m[i, k + i] = 1.0; }
        for (int c = 0; c < k; c++)
        {
            int pivot = c;
            for (int r = c + 1; r < k; r++) if (Math.Abs(m[r, c]) > Math.Abs(m[pivot, c])) pivot = r;
            for (int j = 0; j < 2 * k; j++) (m[c, j], m[pivot, j]) = (m[pivot, j], m[c, j]);
            double d = m[c, c];
            for (int j = 0; j < 2 * k; j++) m[c, j] /= d;
            for (int r = 0; r < k; r++)
            {
                if (r == c) continue;
                double f = m[r, c];
                for (int j = 0; j < 2 * k; j++) m[r, j] -= f * m[c, j];
            }
        }
        var inverse = new double[k, k];
        for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) inverse[i, j] = m[i, k + j];
        return inverse;
    }

    [Fact]
    public void Bars_NoKnotEvidence_IsTheDenseClosedForm()
    {
        // With no knots BARS's evidence is the empty cubic design's marginal (BarsEvidenceTests); here against the definition: τ = 1/m, the
        // centre the median projected on the design, σ² ~ IG(a, b) as the run resolved it from the response.
        var rng = new Xoshiro256PlusPlus(2026_1154);
        const int n = 60;
        double[] x = [.. Enumerable.Range(0, n).Select(i => (i + 0.5) / n)], y = [.. x.Select(v => Math.Sin(2.0 * Math.PI * v) + 0.2 * Gaussian(rng))];
        var model = new WeightedNormalModel();
        var ensemble = new BarsEnsemble(new SplineBasis(3), model, new UniformPrior(0), new LocalBetaKernel(50.0), temperLevels: 5,
                                        temperReferenceRung: true, evidence: true);
        BarsResult result = ensemble.Run(x, y, grid: x, chains: 2, masterSeed: 5, burn: 200, samples: 300);
        double[,] design = new SplineBasis(3).Design(new KnotConfig([]), x);
        double[][] columns = [.. Enumerable.Range(0, design.GetLength(1)).Select(j => Enumerable.Range(0, n).Select(i => design[i, j]).ToArray())];
        (InverseGamma? sigma2, double level) = model.ResolvedPriors(y);
        double[] centre = Projected(columns, level);
        double dense = DenseMarginal(columns, y, null, 1.0 / n, centre, sigma2!.Value.Shape, sigma2.Value.Scale);
        ModelEvidence evidence = result.Evidence!.ForAveraging(y);
        output.WriteLine($"BARS with no knots: ln p(y) {evidence.LogEvidence:F10}, the definition {dense:F10}");
        Assert.True(Math.Abs(evidence.LogEvidence - dense) <= 1e-9 * (1.0 + Math.Abs(dense)), $"{evidence.LogEvidence:R} against {dense:R}");
        Assert.Equal(0.0, evidence.MonteCarloError);
        Assert.Equal(ModelEvidence.KeyOf(y), evidence.ResponseKey);
        Assert.StartsWith("BARS: ", evidence.Priors);
        Assert.Empty(evidence.Flags);
    }

    // Zμ₀ = level·Z G⁻¹ Zᵀ1, the level projected on the design's span, G = ZᵀZ.
    private static double[] Projected(double[][] columns, double level)
    {
        int m = columns[0].Length, k = columns.Length;
        var g = new double[k, k];
        var c = new double[k];
        for (int p = 0; p < k; p++)
        {
            c[p] = columns[p].Sum();
            for (int q = 0; q < k; q++) for (int i = 0; i < m; i++) g[p, q] += columns[p][i] * columns[q][i];
        }
        double[,] gi = Inverse(g);
        var coefficients = new double[k];
        for (int p = 0; p < k; p++) for (int q = 0; q < k; q++) coefficients[p] += gi[p, q] * c[q];
        return [.. Enumerable.Range(0, m).Select(i => level * Enumerable.Range(0, k).Sum(p => columns[p][i] * coefficients[p]))];
    }

    // One predictor, at most one basis function, the g-prior at τ held at 0.7, σ² ~ IG(3, 2), centred at zero: proper priors, so ln p(y) is
    // defined, and F small enough to enumerate.
    private static (BassConfig Config, BassData Data, double[] Y) BassFixture(double[]? weights, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        const int n = 20;
        double[] x = [.. Enumerable.Range(0, n).Select(i => (i + 0.5) / n)], y = [.. x.Select(v => 0.8 + 1.4 * Math.Max(v - 0.4, 0.0) + 0.3 * Gaussian(rng))];
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxBasis = 1, MaxInteraction = 1 },
            Coefficients = BassCoefficientPrior.G(TauPrior.Held(0.7)),
            Intercept = GPriorCentre.Zero,
            Sigma2 = Sigma2Prior.InverseGamma(3, 2),
            Sampling = new BassSampling { Chains = 4, MasterSeed = seed, BurnIn = 1000, Samples = 4000, ReportEvidence = true,
                                          Tempering = new TemperingSpec { ReferenceRung = true } },
        };
        var data = new BassData([new ContinuousPredictor("x", x, new InputDomain(0, 1))], new ScalarResponse(y), null,
                                weights is null ? null : new CaseWeights(weights, WeightMeaning.Precision, WeightedPriorGram.Weighted));
        return (config, data, y);
    }

    // ln p(y) over F from the definition: each structure's design, the intercept and its term's evaluated column, under the dense marginal,
    // weighed by the structure's prior normalized on F. Under weights the design holds its columns whitened, W^½z, so each is read back as
    // z = (W^½z)/√wᵢ; its centring does not matter, the intercept being in the span.
    private static double DenseEvidence(BassModel model, double[] y, double[]? weights)
    {
        List<(BassTerm[] Terms, double LogPrior)> states = BassTermSpace.F(model);
        BassDesign design = model.NewDesign();
        var terms = new List<double>();
        foreach ((BassTerm[] structure, double logPrior) in states)
        {
            Assert.Null(design.Load(structure));
            double[][] cases = design.Snapshot(shareColumns: true).Columns!.Cases;
            double[][] columns = [[.. y.Select(_ => 1.0)], .. cases.Take(structure.Length).Select(c => c.Select((v, i) => weights is null ? v : v / Math.Sqrt(weights[i])).ToArray())];
            terms.Add(logPrior + DenseMarginal(columns, y, weights, 0.7, new double[y.Length], 3.0, 2.0));
        }
        return LogSumExp(terms) - LogSumExp(states.Select(s => s.LogPrior));
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Bass_TheEnumeratedEvidence_IsTheDenseClosedForm(bool weighted)
    {
        // BASS's exact law over F reads its own collapsed marginal (BassTemperingTests.Exact); here the same sum from the definition, so the
        // constants −(m/2) ln 2π and, under precision weights, ½ Σ ln wᵢ are pinned without Monte Carlo.
        var rng = new Xoshiro256PlusPlus(2026_1155);
        double[]? weights = weighted ? [.. Enumerable.Range(0, 20).Select(_ => 0.5 + 1.5 * rng.NextDouble())] : null;
        (BassConfig config, BassData data, double[] y) = BassFixture(weights, 2026_1156);
        BassModel model = BassModel.Create(config, data);
        double engine = BassTemperingTests.Exact(model, 1.0)[0].LogZ, dense = DenseEvidence(model, y, weights);
        output.WriteLine($"{(weighted ? "weighted" : "unweighted")}: BASS's enumeration {engine:F10}, the definition {dense:F10}, {BassTermSpace.F(model).Count} structures");
        Assert.True(Math.Abs(engine - dense) <= 1e-9 * (1.0 + Math.Abs(dense)), $"{engine:R} against {dense:R}");
    }

    [Fact]
    public void Bass_TheFittedEvidence_MeetsTheDenseClosedForm_AndCarriesItsKey()
    {
        (BassConfig config, BassData data, double[] y) = BassFixture(null, 2026_1157);
        BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        double dense = DenseEvidence(BassModel.Create(config, data), y, null);
        ModelEvidence evidence = fit.Evidence!.ForAveraging(data);
        output.WriteLine($"BASS's tempered evidence {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4}, the definition {dense:F4}; flags [{string.Join("; ", evidence.Flags)}]");
        Assert.True(Math.Abs(evidence.LogEvidence - dense) <= 4 * evidence.MonteCarloError + 0.01, $"{evidence.LogEvidence:F4} against {dense:F4}");
        Assert.Equal(ModelEvidence.KeyOf(y), evidence.ResponseKey);
        Assert.StartsWith("BASS: ", evidence.Priors);

        // The weights and their meaning are part of the key.
        double[] w = [.. y.Select((_, i) => 1.0 + 0.1 * i)];
        var weighed = new BassData(data.Predictors, data.Response, null, new CaseWeights(w, WeightMeaning.Precision));
        var counted = new BassData(data.Predictors, data.Response, null, new CaseWeights(w, WeightMeaning.Replication));
        Assert.NotEqual(evidence.ResponseKey, fit.Evidence.ForAveraging(weighed).ResponseKey);
        Assert.NotEqual(fit.Evidence.ForAveraging(weighed).ResponseKey, fit.Evidence.ForAveraging(counted).ResponseKey);
    }

    [Fact]
    public void OnSharedData_TheEnginesEvidencesReadOneResponse_AndAverage()
    {
        // BASS and BARS on one response, n 60: their evidences share a key, differ far less than the (n/2) ln 2π = 55 nats either would drop,
        // and CBMA averages their full-step seams with them.
        var rng = new Xoshiro256PlusPlus(2026_1158);
        const int n = 60;
        double[] x = [.. Enumerable.Range(0, n).Select(i => (i + 0.5) / n)], y = [.. x.Select(v => Math.Sin(2.0 * Math.PI * v) + 0.25 * Gaussian(rng))];
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxBasis = 6, MaxInteraction = 1 },
            Intercept = GPriorCentre.Mean,
            Sigma2 = Sigma2Prior.InverseGamma(3, 0.2),
            Sampling = new BassSampling { Chains = 2, MasterSeed = 7, BurnIn = 800, Samples = 1500, ReportEvidence = true,
                                          Tempering = new TemperingSpec { ReferenceRung = true } },
        };
        var data = new BassData([new ContinuousPredictor("x", x, new InputDomain(0, 1))], new ScalarResponse(y));
        BassFit bass = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 2 });
        var bars = new BarsEnsemble(new SplineBasis(3), new WeightedNormalModel(), new PoissonPrior(3.0), new LocalBetaKernel(50.0), temperLevels: 6,
                                    temperReferenceRung: true, evidence: true);
        BarsResult barsResult = bars.Run(x, y, grid: x, chains: 2, masterSeed: 9, burn: 800, samples: 1500);
        ModelEvidence fromBass = bass.Evidence!.ForAveraging(data), fromBars = barsResult.Evidence!.ForAveraging(y);
        output.WriteLine($"BASS ln p(y) {fromBass.LogEvidence:F3} ± {fromBass.MonteCarloError:F3}; BARS {fromBars.LogEvidence:F3} ± {fromBars.MonteCarloError:F3}");
        Assert.Equal(fromBass.ResponseKey, fromBars.ResponseKey);
        Assert.True(Math.Abs(fromBass.LogEvidence - fromBars.LogEvidence) < 0.35 * 0.5 * n * Math.Log(2.0 * Math.PI),
                    "the engines' evidences should differ far less than the Gaussian constant");

        double[] cases = [0.15, 0.6];
        double[] grid = [.. Enumerable.Range(0, 81).Select(g => -2.0 + 0.05 * g)];
        var bassModel = new AveragedModel(bass.PointwiseLogDensity(data), bass.CandidateLogDensity(data, new PredictionCases([new ContinuousColumn("x", cases)]), [grid, grid]), fromBass);
        var barsModel = new AveragedModel(bars.PointwiseLogDensity(barsResult, x, y), bars.CandidateLogDensity(barsResult, cases, [grid, grid]), fromBars);
        AveragedPrediction averaged = ModelAveraging.FullStep([bassModel, barsModel], 0.1);
        EvidenceReading reading = averaged.Evidence!;
        output.WriteLine($"P(BASS | Z) {reading.Posterior[0]:F3}, over ±2 errors [{reading.PosteriorLow[0]:F3}, {reading.PosteriorHigh[0]:F3}]; fragile {string.Join(", ", reading.Fragile)}");
        Assert.Equal(1.0, reading.Posterior.Sum(), 12);
        Assert.Equal(2, averaged.Prediction.Sets.Count);
        Assert.All(averaged.Prediction.Sets, s => Assert.False(s.WholeSpace));
        Assert.StartsWith("BASS: ", reading.Priors[0]);
        Assert.StartsWith("BARS: ", reading.Priors[1]);

        // The other engine's response refused: BARS read on a shifted copy.
        ModelEvidence elsewhere = barsResult.Evidence.ForAveraging([.. y.Select(v => v + 1.0)]);
        Assert.Throws<ArgumentException>(() => ModelAveraging.FullStep([bassModel, barsModel with { Evidence = elsewhere }]));
    }
}
