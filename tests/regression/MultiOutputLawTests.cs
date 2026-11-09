using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E34's chain and draws against their laws (the multi-output brief, slice 1).
/// <list type="bullet">
/// <item><b>The law by enumeration:</b> on a small F (M ≤ 2, one-way hinges) at K = 2 and 3, a chain's count law, E[ln τ] and its six
/// likeliest structures against each structure's posterior with τ integrated by the trapezoid.</item>
/// <item><b>Geweke's joint test</b> (qualification): priors that read no response, the zero centre and a declared IW(ν₀, Ψ₀), τ sampled;
/// the prior's draws of (structure, τ, Σ, A) and Y against the successive chain of one sweep and one composition per step.</item>
/// <item><b>Recovery:</b> a fit of three outputs from surfaces outside the model with a known, correlated Σ.</item>
/// <item><b>The cost row</b> (qualification, indicative): a sweep at K = 1, 2, 5 and 10 against scalar BASS.</item>
/// </list>
/// </summary>
public sealed class MultiOutputLawTests(ITestOutputHelper output)
{
    // n = 24 on a 6 × 4 grid; K outputs from surfaces outside the model and correlated noise.
    private static BassData SmallData(int k, int seed)
    {
        const int n = 24;
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x0 = [.. Enumerable.Range(0, n).Select(i => i % 6 / 5.0)];
        double[] x1 = [.. Enumerable.Range(0, n).Select(i => i / 6 / 3.0)];
        var y = new double[n, k];
        for (int i = 0; i < n; i++)
        {
            double shared = 0.15 * Normal.SampleStandard(rng);
            for (int l = 0; l < k; l++)
                y[i, l] = (l + 1) * (l switch { 0 => Math.Sin(2.5 * x0[i]), 1 => Math.Abs(x1[i] - 0.4), _ => x0[i] * x1[i] } + shared * (l == 1 ? -1.0 : 1.0)
                                     + 0.1 * Normal.SampleStandard(rng));
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)],
                            new MultiOutputResponse(y, [.. Enumerable.Range(0, k).Select(l => $"y{l}")]));
    }

    private static readonly BassConfig SmallConfig = new() { Structure = new BassStructureSpec { MaxBasis = 2, MaxInteraction = 1 } };

    private static string Key(IEnumerable<BassTerm> terms) => string.Join(" + ", terms.Select(t => t.ToString()).OrderBy(s => s, StringComparer.Ordinal));

    // ── The law by enumeration ────────────────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [InlineData(2, 20261120)]
    [InlineData(3, 20261121)]
    public void Chain_MeetsTheExactLaw_ByEnumeration(int k, int seed)
    {
        BassModel model = BassModel.Create(SmallConfig, SmallData(k, seed));
        BassTarget target = model.Target;
        List<(BassTerm[] Terms, double LogPrior)> states = BassTermSpace.F(model);
        BassDesign design = model.NewDesign();
        const double lowest = -60.0, step = 0.05;
        double[] nodes = [.. Enumerable.Range(0, (int)Math.Round((8.0 - lowest) / step) + 1).Select(i => lowest + step * i)];
        var logWeight = new double[states.Count];
        var logTauGiven = new double[states.Count];
        var index = new Dictionary<string, int>();
        var values = new double[nodes.Length];
        for (int s = 0; s < states.Count; s++)
        {
            BassTerm[] terms = states[s].Terms;
            index.Add(Key(terms), s);
            Assert.Null(design.Load(terms));
            BassStructureFit fit = BassStructureFit.Of(design, terms.Sum(model.Prior.LogPrior));
            for (int j = 0; j < nodes.Length; j++)
            {
                BassEvaluation e = target.Evaluate(fit, nodes[j]);
                values[j] = e.OffSupport ? double.NegativeInfinity : e.LogPrior + e.LogLikelihood(1.0);
            }
            double m = values.Max(), sum = 0.0, u = 0.0;
            for (int j = 0; j < nodes.Length; j++)
            {
                double w = Math.Exp(values[j] - m);
                sum += w;
                u += w * nodes[j];
            }
            logWeight[s] = m + Math.Log(sum * step);
            logTauGiven[s] = u / sum;
        }
        double top = logWeight.Max();
        double[] pi = [.. logWeight.Select(w => Math.Exp(w - top))];
        double total = pi.Sum();
        for (int s = 0; s < pi.Length; s++) pi[s] /= total;

        var chain = new BassChain(model, [], Math.Log(1.0 / 24), new Xoshiro256PlusPlus(seed), new Xoshiro256PlusPlus(seed + 1));
        for (int t = 0; t < 5_000; t++) chain.Step();
        const int steps = 200_000, batches = 40;
        var visits = new int[steps];
        var logTau = new double[steps];
        for (int t = 0; t < steps; t++)
        {
            BassState s = chain.Step();
            Assert.True(index.TryGetValue(Key(s.Terms), out visits[t]), $"the chain holds a structure outside the enumeration: {Key(s.Terms)}");
            logTau[t] = s.LogTau;
        }

        var readings = new List<(string Name, double Exact, double[] Series)>();
        for (int m = 0; m <= model.MaxBasis; m++)
        {
            int mm = m;
            readings.Add(($"P(M = {m})", Enumerable.Range(0, pi.Length).Where(s => states[s].Terms.Length == mm).Sum(s => pi[s]),
                          [.. visits.Select(v => states[v].Terms.Length == mm ? 1.0 : 0.0)]));
        }
        readings.Add(("E[ln τ]", Enumerable.Range(0, pi.Length).Sum(s => pi[s] * logTauGiven[s]), logTau));
        foreach (int s in Enumerable.Range(0, pi.Length).OrderByDescending(s => pi[s]).Take(6))
        {
            int ss = s;
            readings.Add(($"P({(states[s].Terms.Length == 0 ? "∅" : Key(states[s].Terms))})", pi[s], [.. visits.Select(v => v == ss ? 1.0 : 0.0)]));
        }
        double worst = 0.0;
        foreach ((string name, double exact, double[] series) in readings)
        {
            double mean = series.Average(), se = KhaosChainTests.BatchError(series, batches);
            double z = se > 0.0 ? (mean - exact) / se : (Math.Abs(mean - exact) < 1e-3 ? 0.0 : double.PositiveInfinity);
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{name,-48} chain {mean:F4}  exact {exact:G6}  z {z,6:F2}");
        }
        output.WriteLine($"K = {k}: |F| = {states.Count}, {steps} steps, worst |z| {worst:F2}");
        Assert.True(worst < 4.5, $"the chain's worst |z| {worst:F2} against the exact law");
    }

    // ── Geweke (qualification) ────────────────────────────────────────────────────────────────────────────────────────────

    private const int GewekeN = 12;
    private static readonly double[,] GewekePsi = { { 1.2, -0.4 }, { -0.4, 0.8 } };

    private static BassConfig GewekeConfig => new()
    {
        Structure = new BassStructureSpec { MaxBasis = 3, MaxInteraction = 2 },
        Intercept = GPriorCentre.Zero,
        OutputCovariance = OutputCovariancePrior.InverseWishart(6.0, GewekePsi),
    };

    private static BassData GewekeData(double[,] y)
    {
        double[] x0 = [.. Enumerable.Range(0, GewekeN).Select(i => (i * 7 % GewekeN) / (GewekeN - 1.0))];
        double[] x1 = [.. Enumerable.Range(0, GewekeN).Select(i => i / (GewekeN - 1.0))];
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new MultiOutputResponse(y, ["a", "b"]));
    }

    private static readonly string[] GewekeNames = ["M", "ln tau", "ln S11", "ln S22", "atanh rho", "atan f0(0)", "atan f1(6)", "x0 in"];

    private static double[] Functions(BassState state, double[,] sigma, double[,] f)
    {
        double rho = sigma[1, 0] / Math.Sqrt(sigma[0, 0] * sigma[1, 1]);
        return [state.Count, state.LogTau, Math.Log(sigma[0, 0]), Math.Log(sigma[1, 1]), Math.Atanh(rho), Math.Atan(f[0, 0]), Math.Atan(f[GewekeN / 2, 1]),
                state.Terms.Any(t => t.Continuous.Any(h => h.Variable == 0)) ? 1.0 : 0.0];
    }

    // f = 1αᵀ + X_cβ_c at the design's structure, from the centred coefficients (α first, then the terms in the design's order).
    private static double[,] Fitted(BassModel model, BassDesign design, double[,] centred)
    {
        int n = model.Data.Response.Cases, k = centred.GetLength(1);
        double[,] xc = MultiOutputCoreTests.Columns(model, design);
        var f = new double[n, k];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < k; l++)
            {
                double v = centred[0, l];
                for (int j = 0; j < design.Count; j++) v += xc[i, j] * centred[j + 1, l];
                f[i, l] = v;
            }
        return f;
    }

    private static double[,] Observe(double[,] f, double[,] sigma, Xoshiro256PlusPlus rng)
    {
        int n = f.GetLength(0), k = f.GetLength(1);
        double l11 = Math.Sqrt(sigma[0, 0]), l21 = sigma[1, 0] / l11, l22 = Math.Sqrt(sigma[1, 1] - l21 * l21);
        var y = new double[n, k];
        for (int i = 0; i < n; i++)
        {
            double z0 = Normal.SampleStandard(rng), z1 = Normal.SampleStandard(rng);
            y[i, 0] = f[i, 0] + l11 * z0;
            y[i, 1] = f[i, 1] + l21 * z0 + l22 * z1;
        }
        return y;
    }

    [Theory]
    [InlineData(20261122)]
    [InlineData(20261123)]
    public void Geweke_TheSweepAndTheDrawLeaveTheJointInvariant(int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int draws = 100_000, steps = 150_000, batches = 50;
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(GewekeConfig, GewekeData(new double[,] { { 0.3, -1.0 }, { 1.1, 0.4 }, { -0.2, 0.9 }, { 0.5, 0.1 }, { 1.4, -0.6 }, { 0.0, 0.2 },
                                                                                    { -0.9, 1.3 }, { 0.6, -0.3 }, { 0.2, 0.7 }, { -1.2, 0.0 }, { 0.8, 0.5 }, { -0.4, -0.8 } }));
        var sampler = new BassPriorSampler(prior);
        InverseWishart sigmaPrior = prior.Resolved.Outputs!.Prior;

        // The prior's draw: (structure, τ) from the reference sampler, Σ from its prior, A | Σ ~ MN(0, (τ·blockdiag(N, G_c))⁻¹, Σ).
        (BassState S, double[,] Sigma, double[,] F) Draw()
        {
            BassState state = sampler.DrawState(rng).State;
            double[,] sigma = sigmaPrior.Sample(rng);
            BassDesign d = prior.NewDesign();
            Assert.Null(d.Load(state.Terms));
            int m = d.Count;
            double tau = Math.Exp(state.LogTau), l11 = Math.Sqrt(sigma[0, 0]), l21 = sigma[1, 0] / l11, l22 = Math.Sqrt(sigma[1, 1] - l21 * l21);
            var centred = new double[m + 1, 2];
            var v = new double[2][];
            for (int l = 0; l < 2; l++)
            {
                v[l] = new double[m];
                if (m == 0) continue;
                double[] z = [.. Enumerable.Range(0, m).Select(_ => Normal.SampleStandard(rng))];
                d.Factor.SampleInnovation(z, v[l]);
            }
            for (int j = 0; j < m; j++)
            {
                centred[j + 1, 0] = l11 * v[0][j] / Math.Sqrt(tau);
                centred[j + 1, 1] = (l21 * v[0][j] + l22 * v[1][j]) / Math.Sqrt(tau);
            }
            double u0 = Normal.SampleStandard(rng), u1 = Normal.SampleStandard(rng), s0 = 1.0 / Math.Sqrt(tau * GewekeN);
            centred[0, 0] = s0 * l11 * u0;
            centred[0, 1] = s0 * (l21 * u0 + l22 * u1);
            return (state, sigma, Fitted(prior, d, centred));
        }
        var marginal = new double[draws][];
        for (int t = 0; t < draws; t++)
        {
            (BassState s, double[,] sigma, double[,] f) = Draw();
            marginal[t] = Functions(s, sigma, f);
        }

        (BassState S, double[,] Sigma, double[,] F) theta = Draw();
        var successive = new double[steps][];
        for (int t = 0; t < steps; t++)
        {
            double[,] y = Observe(theta.F, theta.Sigma, rng);
            BassModel model = BassModel.Create(GewekeConfig, GewekeData(y));
            var chain = new BassChain(model, theta.S.Terms, theta.S.LogTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
            BassState state = chain.Step();
            BassDesign d = chain.Structure.DesignAt(state);
            (double[,] sigma, double[,] centred) = BassOutputDraws.Sample(model.Target, d, state.Fit, Math.Exp(state.LogTau), rng);
            theta = (state, sigma, Fitted(model, d, centred));
            successive[t] = Functions(state, sigma, theta.F);
        }

        double worst = 0.0;
        for (int g = 0; g < GewekeNames.Length; g++)
        {
            double[] a = [.. marginal.Select(r => r[g])], b = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1);
            double meanB = b.Average(), se = Math.Sqrt(varA / a.Length + Math.Pow(KhaosChainTests.BatchError(b, batches), 2));
            double z = (meanA - meanB) / se;
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{GewekeNames[g],-11} marginal-conditional {meanA,10:F5}  successive {meanB,10:F5}  z {z,6:F2}");
        }
        output.WriteLine($"seed {seed}: worst |z| {worst:F2}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ── Recovery ──────────────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Fit_RecoversAKnownSigma_AndSurfacesOutsideTheModel()
    {
        // Three surfaces no hinge sum holds exactly, and Σ with sds 0.2, 0.3, 0.15 and correlations 0.3, −0.2, 0.
        const int n = 300, k = 3;
        double[,] truth = { { 0.04, 0.018, -0.006 }, { 0.018, 0.09, 0.0 }, { -0.006, 0.0, 0.0225 } };
        var rng = new Xoshiro256PlusPlus(20261124);
        double[] x0 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        double[] x1 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        double Surface(int l, int i) => l switch
        {
            0 => 0.8 * Math.Sin(2.0 * Math.PI * x0[i]) + x1[i],
            1 => Math.Exp(-3.0 * (x0[i] - 0.4) * (x0[i] - 0.4)) - 0.5 * x1[i] * x1[i],
            _ => Math.Abs(x0[i] - x1[i]),
        };
        double[,] root = ReferenceAlgebra.Cholesky(truth)!;
        var y = new double[n, k];
        for (int i = 0; i < n; i++)
        {
            double[] z = [Normal.SampleStandard(rng), Normal.SampleStandard(rng), Normal.SampleStandard(rng)];
            for (int a = 0; a < k; a++)
            {
                double e = 0.0;
                for (int d = 0; d <= a; d++) e += root[a, d] * z[d];
                y[i, a] = Surface(a, i) + e;
            }
        }
        var data = new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new MultiOutputResponse(y, ["a", "b", "c"]));
        var config = new BassConfig { Sampling = new BassSampling { Chains = 2, BurnIn = 2000, Samples = 600, MasterSeed = 20261124 } };
        BassFit fit = Bass.Fit(config, data);
        BassModel model = BassModel.Create(config, data);

        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        var mean = new double[k, k];
        var squares = new double[k, k];
        var fitted = new double[n, k];
        var column = new double[n];
        foreach (BassDraw draw in draws)
        {
            for (int a = 0; a < k; a++)
                for (int d = 0; d < k; d++) { mean[a, d] += draw.Sigma![a, d] / draws.Length; squares[a, d] += draw.Sigma[a, d] * draw.Sigma[a, d] / draws.Length; }
            int[] ids = draw.Structure.TermIds;
            for (int l = 0; l < k; l++)
                for (int i = 0; i < n; i++) fitted[i, l] += draw.OutputCoefficients![0, l] / draws.Length;
            for (int j = 0; j < ids.Length; j++)
            {
                model.Basis.CaseColumn(fit.Terms[ids[j]], model.Inputs.Cases, column);
                for (int l = 0; l < k; l++)
                    for (int i = 0; i < n; i++) fitted[i, l] += draw.OutputCoefficients![j + 1, l] * column[i] / draws.Length;
            }
        }
        for (int a = 0; a < k; a++)
            for (int d = 0; d <= a; d++)
            {
                double sd = Math.Sqrt(Math.Max(squares[a, d] - mean[a, d] * mean[a, d], 0.0));
                double scale = Math.Sqrt(truth[a, a] * truth[d, d]);
                output.WriteLine($"Σ[{a},{d}]: posterior {mean[a, d]:F5} ± {sd:F5}, truth {truth[a, d]:F5}");
                Assert.True(Math.Abs(mean[a, d] - truth[a, d]) <= 4.0 * sd + 0.02 * scale, $"Σ[{a},{d}] {mean[a, d]:F5} ± {sd:F5} against {truth[a, d]:F5}");
            }
        for (int l = 0; l < k; l++)
        {
            double rmse = Math.Sqrt(Enumerable.Range(0, n).Average(i => (fitted[i, l] - Surface(l, i)) * (fitted[i, l] - Surface(l, i))));
            output.WriteLine($"output {l}: the posterior mean surface's RMSE against the truth {rmse:F4}, noise sd {Math.Sqrt(truth[l, l]):F4}");
            Assert.True(rmse < 0.6 * Math.Sqrt(truth[l, l]), $"output {l}: RMSE {rmse:F4}");
        }
    }

    // ── The cost row (qualification, indicative) ───────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Cost_TheSweepAtSeveralOutputs_AgainstScalarBass()
    {
        if (!BassQualification.Enabled) return;   // a measurement
        const int n = 1000, sweeps = 3000;
        var rng = new Xoshiro256PlusPlus(20261125);
        double[][] x = [.. Enumerable.Range(0, 4).Select(_ => Enumerable.Range(0, n).Select(_ => rng.NextDouble()).ToArray())];
        double[] y0 = [.. Enumerable.Range(0, n).Select(i => Math.Sin(3.0 * x[0][i]) + x[1][i] * x[2][i] + 0.1 * Normal.SampleStandard(rng))];
        List<BassPredictor> predictors = [.. x.Select((v, j) => (BassPredictor)new ContinuousPredictor($"x{j}", v))];
        var rows = new List<string>();
        double Measure(BassResponse response)
        {
            BassModel model = BassModel.Create(new BassConfig { Structure = new BassStructureSpec { MaxBasis = 20 } }, new BassData(predictors, response));
            var chain = new BassChain(model, [], Math.Log(1.0 / n), new Xoshiro256PlusPlus(7), new Xoshiro256PlusPlus(8));
            for (int t = 0; t < 1000; t++) chain.Step();
            long start = BassCost.Tick();
            for (int t = 0; t < sweeps; t++) chain.Step();
            return (BassCost.Tick() - start) / BassCost.TicksPerSecond / sweeps * 1e6;
        }
        double scalar = Measure(new ScalarResponse(y0));
        rows.Add($"scalar {scalar:F1} µs");
        foreach (int k in new[] { 1, 2, 5, 10 })
        {
            var y = new double[n, k];
            for (int i = 0; i < n; i++)
                for (int l = 0; l < k; l++) y[i, l] = y0[i] * (1.0 + 0.1 * l) + 0.1 * Normal.SampleStandard(rng);
            double cost = Measure(new MultiOutputResponse(y, [.. Enumerable.Range(0, k).Select(l => $"y{l}")]));
            rows.Add($"K = {k} {cost:F1} µs ({cost / scalar:F2}×)");
        }
        output.WriteLine($"clock: {(BassCost.Cycles ? "thread cycles" : "wall")} (indicative); a sweep on n = {n}, p = 4, M_max = 20: {string.Join("; ", rows)}");
    }
}
