using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's sweep (stage 4, S6; D01, D02, D04). Under a constant likelihood the chain samples the prior on F, which this
/// test enumerates: the count law, the order law and τ's prior, on both move schedules (F-O3's and F-O4's laws from the
/// chain, the boundary ratios corrected against R's fixed thirds), with λ's retained draws against the Gamma mixture.
/// For a fixed structure, the collapsed sweep's τ against quadrature of its marginal, and its composition draws
/// (σ², the coefficients) against a long joint-Gibbs reference that samples the coefficients, σ² and τ in turn. The
/// design follows the chain through accepts, rejects and swapped-in states; the schedules' boundary weights; and a
/// retained draw's ids and uncentred coefficients. The posterior with data is checked against the stage-4 oracle in
/// <c>BassChainParityTests</c>.
/// </summary>
public sealed class BassSweepTests
{
    private readonly ITestOutputHelper _out;
    public BassSweepTests(ITestOutputHelper output) => _out = output;

    // x₀ and x₁ on {0, ½, 1}, the grid twice (N = 18), or on a five-value grid three times (N = 75).
    private static BassData Data(bool large = false, int seed = 5, bool exact = false)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        int levels = large ? 5 : 3, reps = large ? 3 : 2;
        var x0 = new List<double>();
        var x1 = new List<double>();
        var y = new List<double>();
        for (int rep = 0; rep < reps; rep++)
            for (int i = 0; i < levels; i++)
                for (int j = 0; j < levels; j++)
                {
                    double a = i / (levels - 1.0), b = j / (levels - 1.0);
                    x0.Add(a);
                    x1.Add(b);
                    double noise = Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());
                    y.Add(exact ? 1 + 2 * Math.Max(a - 0.5, 0) : 1 + 2 * a + 3 * Math.Max(a - 0.5, 0) * Math.Max(0.5 - b, 0) + 0.3 * noise);
                }
        return new BassData([new ContinuousPredictor("x0", [.. x0]), new ContinuousPredictor("x1", [.. x1])], new ScalarResponse([.. y]));
    }

    // A direct functional response: x₀ and x₁ on a five-value grid (n = 25), one output variable t on seven locations
    // (N = 175), Y = 1 + 2x₀ + 1.5·t·x₁ + 2(x₀ − ½)₊(t − ½)₊ + noise.
    private static BassData FunctionalData(int seed = 13)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        const int n = 25, m = 7;
        double[] t = [.. Enumerable.Range(0, m).Select(r => r / (m - 1.0))];
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
        {
            x0[i] = (i / 5) / 4.0;
            x1[i] = (i % 5) / 4.0;
            for (int r = 0; r < m; r++)
            {
                double noise = Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());
                y[i, r] = 1 + 2 * x0[i] + 1.5 * t[r] * x1[i] + 2 * Math.Max(x0[i] - 0.5, 0) * Math.Max(t[r] - 0.5, 0) + 0.2 * noise;
            }
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)],
                            new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
    }

    private static BassConfig Config(MoveSchedule? moves = null, int maxBasis = 3)
        => new()
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = maxBasis },
            Sampling = new BassSampling { Moves = moves ?? new EqualOverAvailable() },
        };

    // The prior alone: the structure prior and τ's, with a constant likelihood.
    private sealed class PriorOnly(BassTarget target) : ISamplingTarget<BassState>
    {
        public ITargetEvaluation Evaluate(BassState state) => new Flat(target.Evaluate(state).LogPrior);

        private sealed record Flat(double LogPrior) : ITargetEvaluation
        {
            public double LogLikelihood(double inverseTemperature) => double.IsNegativeInfinity(LogPrior) ? double.NegativeInfinity : 0.0;
        }
    }

    // Every structure of at most M_max distinct terms in F, by the design's own admission, with its prior weight.
    private static List<(BassTerm[] Terms, double LogPrior)> EnumerateF(BassModel model)
    {
        BassStructurePrior prior = model.Prior;
        var terms = new List<BassTerm>();
        int[] knots = [.. Enumerable.Range(0, 2).Select(v => model.Inputs.ContinuousKnots(v).Length)];
        foreach (int v in new[] { 0, 1 })
            foreach (int s in new[] { -1, 1 })
                for (int k = 0; k < knots[v]; k++) terms.Add(new BassTerm([new HingeFactor(v, s, k)]));
        foreach (int s0 in new[] { -1, 1 })
            foreach (int s1 in new[] { -1, 1 })
                for (int k0 = 0; k0 < knots[0]; k0++)
                    for (int k1 = 0; k1 < knots[1]; k1++)
                        terms.Add(new BassTerm([new HingeFactor(0, s0, k0), new HingeFactor(1, s1, k1)]));
        BassDesign design = model.NewDesign();
        var states = new List<(BassTerm[], double)>();
        void Visit(int from, List<BassTerm> chosen)
        {
            BassTerm[] structure = [.. chosen];
            if (design.Load(structure) is null)
                states.Add((structure, prior.LogCollectionWeight(structure.Length) + structure.Sum(prior.LogPrior)));
            else return;   // F is downward closed (D06), so no superset of a refused structure is in it
            if (chosen.Count == model.MaxBasis) return;
            for (int i = from; i < terms.Count; i++)
            {
                chosen.Add(terms[i]);
                Visit(i + 1, chosen);
                chosen.RemoveAt(chosen.Count - 1);
            }
        }
        Visit(0, []);
        return states;
    }

    private static string OrderKey(IReadOnlyList<BassTerm> terms)
        => $"{terms.Count(t => t.Continuous.Count == 1)},{terms.Count(t => t.Continuous.Count == 2)}";

    // A batch-means standard error over every chain's trace of an indicator or value.
    private static double BatchError(IReadOnlyList<double[]> traces, int batches = 25)
    {
        var means = new List<double>();
        foreach (double[] t in traces)
        {
            int size = t.Length / batches;
            for (int b = 0; b < batches; b++) means.Add(t.Skip(b * size).Take(size).Average());
        }
        double mean = means.Average();
        return Math.Sqrt(means.Sum(m => (m - mean) * (m - mean)) / (means.Count - 1) / means.Count);
    }

    [Theory]
    [InlineData("equal")]
    [InlineData("cancelling")]
    public void ConstantLikelihood_SamplesThePriorOnF(string schedule)
    {
        BassData data = Data();
        var model = BassModel.Create(Config(schedule == "equal" ? new EqualOverAvailable() : new PriorCancelling(0.4)), data);
        List<(BassTerm[] Terms, double LogPrior)> f = EnumerateF(model);
        double z = f.Max(s => s.LogPrior);
        double total = f.Sum(s => Math.Exp(s.LogPrior - z));
        var countLaw = new double[model.MaxBasis + 1];
        var orderLaw = new Dictionary<string, double>();
        foreach ((BassTerm[] terms, double lp) in f)
        {
            double p = Math.Exp(lp - z) / total;
            countLaw[terms.Length] += p;
            orderLaw[OrderKey(terms)] = orderLaw.GetValueOrDefault(OrderKey(terms)) + p;
        }
        double h1 = model.Prior.Count.Shape, h2 = model.Prior.Count.Rate;
        double lambdaMean = Enumerable.Range(0, countLaw.Length).Sum(m => countLaw[m] * (h1 + m) / (h2 + 1));

        const int chains = 4, burn = 2_000, sweeps = 60_000;
        var countTraces = new List<double[]>[countLaw.Length];
        for (int m = 0; m < countLaw.Length; m++) countTraces[m] = [];
        var orderHits = orderLaw.Keys.ToDictionary(k => k, _ => 0.0);
        var logTau = new List<double[]>();
        double lambdaSum = 0;
        int lambdaDraws = 0;
        for (int c = 0; c < chains; c++)
        {
            var chain = new BassChain(model, [], Math.Log(1.0 / 18), new Xoshiro256PlusPlus(100 + c), new Xoshiro256PlusPlus(200 + c),
                                      new PriorOnly(model.Target));
            var table = new BassTermTable();
            var composer = new Xoshiro256PlusPlus(300 + c);
            for (int i = 0; i < burn; i++) chain.Step();
            var counts = new double[countLaw.Length][];
            for (int m = 0; m < countLaw.Length; m++) counts[m] = new double[sweeps];
            var u = new double[sweeps];
            for (int i = 0; i < sweeps; i++)
            {
                BassState s = chain.Step();
                counts[s.Count][i] = 1;
                orderHits[OrderKey(s.Terms)] += 1;
                u[i] = s.LogTau;
                if (i % 20 == 0)
                {
                    lambdaSum += chain.Compose(composer, table).Lambda;
                    lambdaDraws++;
                }
            }
            for (int m = 0; m < countLaw.Length; m++) countTraces[m].Add(counts[m]);
            logTau.Add(u);
            Assert.True(chain.Structure.Refusals.Total > 0);
        }

        _out.WriteLine($"{schedule}: {f.Count} states in F");
        for (int m = 0; m < countLaw.Length; m++)
        {
            double observed = countTraces[m].Sum(t => t.Sum()) / (chains * sweeps), error = BatchError(countTraces[m]);
            _out.WriteLine($"  M = {m}: {observed:F4} against {countLaw[m]:F4} (± {error:F4})");
            Assert.True(Math.Abs(observed - countLaw[m]) <= 5 * error + 0.002, $"M = {m}: {observed:F4} against {countLaw[m]:F4} (± {error:F4})");
        }
        // The boundary ratio π(1)/π(0), which R's fixed thirds would make three times the target's (F-O4).
        double ratio = countTraces[1].Sum(t => t.Sum()) / countTraces[0].Sum(t => t.Sum());
        Assert.True(Math.Abs(ratio / (countLaw[1] / countLaw[0]) - 1) <= 0.05, $"π(1)/π(0) = {ratio:F3} against {countLaw[1] / countLaw[0]:F3}");
        foreach ((string key, double p) in orderLaw)
        {
            double observed = orderHits[key] / (chains * sweeps);
            Assert.True(Math.Abs(observed - p) <= 0.01, $"orders ({key}): {observed:F4} against {p:F4}");
        }
        // τ's law is its prior: E[ln τ] = ψ(a) − ln(aN).
        double expectedLogTau = SpecialFunctions.Digamma(0.5) - Math.Log(0.5 * 18);
        double meanLogTau = logTau.Average(t => t.Average());
        Assert.True(Math.Abs(meanLogTau - expectedLogTau) <= 5 * BatchError(logTau) + 0.01, $"E[ln τ] = {meanLogTau:F3} against {expectedLogTau:F3}");
        // λ at retention: the Gamma(h₁ + M, h₂ + 1) mixture over the count law.
        double lambda = lambdaSum / lambdaDraws;
        Assert.True(Math.Abs(lambda - lambdaMean) <= 0.03 * lambdaMean, $"E[λ] = {lambda:F4} against {lambdaMean:F4}");
    }

    // A palette that never moves the structure, so a scan runs τ alone.
    private sealed class Stay : IRjMove<BassState>
    {
        public string Key => "stay";
        public string ReverseKey => "stay";
        public double Weight(BassState state) => 1.0;
        public Proposal<BassState>? Propose(BassState current, Xoshiro256PlusPlus rng) => null;
    }

    [Fact]
    public void FixedStructure_TauAndComposition_MatchQuadratureAndAJointGibbsReference()
    {
        BassData data = Data(large: true);
        var model = BassModel.Create(Config(maxBasis: 10), data);
        BassTerm[] structure = [new BassTerm([new HingeFactor(0, 1, 2)]), new BassTerm([new HingeFactor(0, 1, 2), new HingeFactor(1, -1, 2)])];
        var moves = new BassStructureMoves(model, structure, Math.Log(0.02));
        var chain = new ReversibleJumpChain<BassState>([new Stay()], model.Target, moves.Start, new Xoshiro256PlusPlus(7));
        var scan = new BlockScan<ReversibleJumpChain<BassState>>([new BassTauBlock(model.Target)], chain, new Xoshiro256PlusPlus(8), c => c);
        var composer = new Xoshiro256PlusPlus(9);
        var table = new BassTermTable();
        const int sweeps = 100_000;
        int n = model.Resolved.Observations;
        var tau = new double[sweeps];
        var sigma = new double[sweeps];
        var a = new double[3][];
        for (int j = 0; j < 3; j++) a[j] = new double[sweeps];
        for (int i = 0; i < 1_000; i++) scan.Step();
        for (int i = 0; i < sweeps; i++)
        {
            BassState s = scan.Step().Current;
            BassDraw draw = moves.Compose(s, composer, table);
            tau[i] = draw.Tau;
            sigma[i] = draw.Sigma2;
            for (int j = 0; j < 3; j++) a[j][i] = draw.Coefficients[j];
        }

        // τ | structure, y by quadrature of the collapsed marginal on u.
        BassStructureFit fit = moves.Start.Fit;
        var nodes = BassTargetTests.Nodes(Math.Log(1.0 / n) - 80, Math.Log(1.0 / n) + 8).ToList();
        double[] logw = [.. nodes.Select(p => Math.Log(p.W) + model.Target.LogTauDensity(p.U) + model.Target.Evaluate(fit, p.U).LogMarginal)];
        double top = logw.Max();
        double norm = logw.Sum(l => Math.Exp(l - top));
        double exactTau = nodes.Select((p, i) => Math.Exp(p.U) * Math.Exp(logw[i] - top)).Sum() / norm;

        // The joint Gibbs reference on the explicit design [1, Z_c], ids in the draws' order.
        (double[] gTau, double[] gSigma, double[][] gA) = JointGibbs(model, structure, table, 400_000, seed: 31);

        void Check(string name, double[] collapsed, double[] gibbs)
        {
            double c = collapsed.Average(), g = gibbs.Average();
            double error = Math.Sqrt(Math.Pow(BatchError([collapsed]), 2) + Math.Pow(BatchError([gibbs]), 2));
            _out.WriteLine($"{name}: collapsed {c:G6}, joint Gibbs {g:G6} (± {error:G2})");
            Assert.True(Math.Abs(c - g) <= 5 * error, $"{name}: collapsed {c:G6} against joint Gibbs {g:G6} (± {error:G2})");
        }
        _out.WriteLine($"E[τ]: collapsed {tau.Average():G6}, quadrature {exactTau:G6}");
        Assert.True(Math.Abs(tau.Average() - exactTau) <= 5 * BatchError([tau]), $"E[τ] {tau.Average():G6} against {exactTau:G6}");
        Check("τ", tau, gTau);
        Check("σ²", sigma, gSigma);
        Check("σ⁴", [.. sigma.Select(v => v * v)], [.. gSigma.Select(v => v * v)]);
        for (int j = 0; j < 3; j++) Check($"a{j}", a[j], gA[j]);
        Check("a₁²", [.. a[1].Select(v => v * v)], [.. gA[1].Select(v => v * v)]);
    }

    // Gibbs on (a, σ², τ) for a fixed structure under D03's defaults (level centring, σ² scaled, Zellner–Siow τ), on
    // the explicit centred design X = [1, Z_c] with prior mean μ₀ = ȳe₀ and G = XᵀX:
    //   a | σ², τ ~ N((â + τμ₀)/(1 + τ), σ²G⁻¹/(1 + τ)),
    //   σ² | a, τ ~ IG(c + (N + k)/2, cs² + ‖y − Xa‖²/2 + τ(a − μ₀)ᵀG(a − μ₀)/2),
    //   τ | a, σ² ~ Gamma(a_τ + k/2, rate a_τN + (a − μ₀)ᵀG(a − μ₀)/(2σ²)).
    // The coefficients are returned uncentred, in the composed draws' order.
    private static (double[] Tau, double[] Sigma2, double[][] A) JointGibbs(BassModel model, BassTerm[] structure, BassTermTable table,
                                                                            int sweeps, int seed)
    {
        int n = model.Resolved.Observations, k = structure.Length + 1;
        double[] y = [.. model.Data.Response.Values.ToArray()];
        double ybar = model.Resolved.ResponseMean;
        var columns = new double[structure.Length][];
        var means = new double[structure.Length];
        for (int j = 0; j < structure.Length; j++)
        {
            columns[j] = new double[n];
            model.Basis.CaseColumn(structure[j], model.Inputs.Cases, columns[j]);
            means[j] = columns[j].Average();
        }
        var x = new double[n, k];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = 1;
            for (int j = 0; j < structure.Length; j++) x[i, j + 1] = columns[j][i] - means[j];
        }
        var g = new double[k, k];
        var xty = new double[k];
        for (int r = 0; r < k; r++)
        {
            for (int i = 0; i < n; i++) xty[r] += x[i, r] * y[i];
            for (int c = 0; c < k; c++)
                for (int i = 0; i < n; i++) g[r, c] += x[i, r] * x[i, c];
        }
        double[,] l = Cholesky(g);
        double[] aHat = Solve(l, xty);
        var mu = new double[k];
        mu[0] = ybar;
        InverseGamma sigmaPrior = model.Resolved.Sigma2!.Value;
        Gamma tauPrior = model.Resolved.TauPrior!.Value;
        var rng = new Xoshiro256PlusPlus(seed);
        double tau = 1.0 / n, sigma2 = 0.1;
        var a = new double[k];
        var outTau = new double[sweeps];
        var outSigma = new double[sweeps];
        var outA = new double[k][];
        for (int j = 0; j < k; j++) outA[j] = new double[sweeps];
        int[] ids = [.. structure.Select(t => table.Intern(t))];
        int[] order = [.. Enumerable.Range(0, structure.Length).OrderBy(j => ids[j])];
        for (int it = -2_000; it < sweeps; it++)
        {
            // a | σ², τ: mean (â + τμ₀)/(1 + τ), covariance σ²/(1 + τ)·G⁻¹ = (L Lᵀ)⁻¹·scale, drawn by back-substitution.
            double scale = Math.Sqrt(sigma2 / (1 + tau));
            var z = new double[k];
            for (int j = 0; j < k; j++) z[j] = Normal.SampleStandard(rng);
            double[] w = BackSolve(l, z);
            for (int j = 0; j < k; j++) a[j] = (aHat[j] + tau * mu[j]) / (1 + tau) + scale * w[j];
            double rss = 0, quad = 0;
            for (int i = 0; i < n; i++)
            {
                double fitted = 0;
                for (int j = 0; j < k; j++) fitted += x[i, j] * a[j];
                rss += (y[i] - fitted) * (y[i] - fitted);
            }
            for (int r = 0; r < k; r++)
                for (int c = 0; c < k; c++) quad += (a[r] - mu[r]) * g[r, c] * (a[c] - mu[c]);
            sigma2 = InverseGamma.ShapeScale(sigmaPrior.Shape + 0.5 * (n + k), sigmaPrior.Scale + 0.5 * rss + 0.5 * tau * quad).Sample(rng);
            tau = Gamma.ShapeRate(tauPrior.Shape + 0.5 * k, tauPrior.Rate + 0.5 * quad / sigma2).Sample(rng);
            if (it < 0) continue;
            outTau[it] = tau;
            outSigma[it] = sigma2;
            double intercept = a[0];
            for (int q = 0; q < structure.Length; q++)
            {
                int j = order[q];
                outA[q + 1][it] = a[j + 1];
                intercept -= a[j + 1] * means[j];
            }
            outA[0][it] = intercept;
        }
        return (outTau, outSigma, outA);
    }

    private static double[,] Cholesky(double[,] a)
    {
        int k = a.GetLength(0);
        var l = new double[k, k];
        for (int i = 0; i < k; i++)
            for (int j = 0; j <= i; j++)
            {
                double s = a[i, j];
                for (int q = 0; q < j; q++) s -= l[i, q] * l[j, q];
                l[i, j] = i == j ? Math.Sqrt(s) : s / l[j, j];
            }
        return l;
    }

    // G⁻¹b with G = LLᵀ.
    private static double[] Solve(double[,] l, double[] b)
    {
        int k = b.Length;
        var v = new double[k];
        for (int i = 0; i < k; i++)
        {
            double s = b[i];
            for (int q = 0; q < i; q++) s -= l[i, q] * v[q];
            v[i] = s / l[i, i];
        }
        return BackSolve(l, v);
    }

    // L⁻ᵀz: a draw with covariance G⁻¹ when z is standard normal.
    private static double[] BackSolve(double[,] l, double[] z)
    {
        int k = z.Length;
        var w = new double[k];
        for (int i = k - 1; i >= 0; i--)
        {
            double s = z[i];
            for (int q = i + 1; q < k; q++) s -= l[q, i] * w[q];
            w[i] = s / l[i, i];
        }
        return w;
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Design_FollowsTheChain(bool functional)
    {
        BassData data = functional ? FunctionalData() : Data(large: true);
        var model = BassModel.Create(Config(maxBasis: 12), data);
        var chain = new BassChain(model, [], Math.Log(0.02), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        var other = new BassChain(model, [], Math.Log(0.02), new Xoshiro256PlusPlus(3), new Xoshiro256PlusPlus(4));
        int checks = 0;
        for (int i = 1; i <= 4_000; i++)
        {
            chain.Step();
            other.Step();
            if (i == 2_000)
            {
                // A swapped-in state, as tempering hands one over: the design loads it whole.
                BassState foreign = other.Current;
                chain.Moves.SetState(foreign, model.Target.Evaluate(foreign));
            }
            if (i % 100 != 0) continue;
            BassState s = chain.Current;
            BassDesign live = chain.Structure.DesignAt(s);
            BassDesign fresh = model.NewDesign();
            Assert.Null(fresh.Load(s.Terms));
            Assert.Equal(s.Terms, live.Terms);
            double scale = live.ResponseSquares;
            Assert.True(Math.Abs(live.ResidualSquares - fresh.ResidualSquares) <= 1e-10 * scale);
            Assert.True(Math.Abs(live.FittedSquares - fresh.FittedSquares) <= 1e-10 * scale);
            Assert.True(Math.Abs(s.Fit.ResidualSquares - fresh.ResidualSquares) <= 1e-10 * scale);
            Assert.Equal(s.Count, s.Fit.Terms);
            Assert.Equal(s.Terms.Sum(model.Prior.LogPrior), s.Fit.LogTermPrior, 10);
            checks++;
        }
        Assert.Equal(40, checks);
        Assert.True(chain.Current.Count > 0);
        _out.WriteLine($"M = {chain.Current.Count}; moves {string.Join(", ", chain.Moves.MoveStats().Select(m => $"{m.Key} {m.Accepted}/{m.Attempts}"))}; " +
                       $"refusals {string.Join(", ", chain.Structure.Refusals.Entries.Select(e => $"{e.Category} {e.Count}"))}; empty births {chain.Structure.EmptyBirths}");
    }

    [Fact]
    public void Schedules_AtTheBoundaries()
    {
        BassData data = Data();
        foreach (MoveSchedule schedule in new MoveSchedule[] { new EqualOverAvailable(), new PriorCancelling(0.4) })
        {
            var model = BassModel.Create(Config(schedule, maxBasis: 2), data);
            var moves = new BassStructureMoves(model, [], Math.Log(0.05));
            BassState empty = moves.Start;
            BassTerm one = new([new HingeFactor(0, 1, 1)]), two = new([new HingeFactor(1, 1, 1)]);
            BassState single = new BassStructureMoves(model, [one], 0).Start;
            BassState full = new BassStructureMoves(model, [one, two], 0).Start;
            double[] W(BassState s) => [.. moves.Moves.Select(m => m.Weight(s))];
            double[] w0 = W(empty), w1 = W(single), w2 = W(full);
            if (schedule is EqualOverAvailable)
            {
                Assert.Equal([1.0, 0.0, 0.0], w0);
                Assert.Equal([1.0, 1.0, 1.0], w1);
                Assert.Equal([0.0, 1.0, 1.0], w2);
            }
            else
            {
                // b_M p(M) = d_{M+1} p(M+1), and each M's weights sum to one.
                double r0 = Math.Exp(model.Prior.LogCountRatio(0)), r1 = Math.Exp(model.Prior.LogCountRatio(1));
                Assert.Equal(0.4 * Math.Min(1, r0), w0[0], 15);
                Assert.Equal(0.4 * Math.Min(1, 1 / r0), w1[1], 15);
                Assert.Equal(r0, w0[0] / w1[1], 12);
                Assert.Equal(0.0, w2[0]);
                Assert.Equal(0.4 * Math.Min(1, 1 / r1), w2[1], 15);
                foreach (double[] w in new[] { w0, w1, w2 }) Assert.Equal(1.0, w.Sum(), 15);
            }
        }
    }

    [Fact]
    public void Compose_AlignsOnDirectFunctionalData()
    {
        // A case-only term, an output-only term (its case side is 1) and a mixed one: the draw's a₀ + Σ aⱼ·aⱼ(xᵢ)hⱼ(t_r)
        // against the same draw in the centred basis, α + Σ βⱼ(aⱼ(xᵢ)hⱼ(t_r) − āⱼh̄ⱼ), at all N = n·m observations.
        BassData data = FunctionalData();
        var model = BassModel.Create(Config(maxBasis: 10), data);
        BassTerm[] structure = [new BassTerm([new HingeFactor(0, 1, 2)]), new BassTerm(functional: [new HingeFactor(0, 1, 3)]),
                                new BassTerm([new HingeFactor(1, -1, 3)], functional: [new HingeFactor(0, -1, 4)])];
        var moves = new BassStructureMoves(model, structure, Math.Log(0.01));
        var table = new BassTermTable();
        table.Intern(structure[2]);
        BassDraw draw = moves.Compose(moves.Start, new Xoshiro256PlusPlus(5), table);

        BassDesign design = moves.DesignAt(moves.Start);
        CollapsedGPrior block = CollapsedGPrior.FromCentredLeastSquares(design.Factor, design.Coefficients.ToArray(), design.ResidualSquares,
                                                                         design.FittedSquares, design.Observations, design.ResponseMean, GPriorCentre.Mean);
        var centred = new double[4];
        Assert.Equal(block.At(moves.Start.Tau).Sample(model.Resolved.Sigma2, new Xoshiro256PlusPlus(5), centred), draw.Sigma2);
        int n = model.Inputs.Cases.Cases, m = model.Inputs.Outputs!.Locations;
        double[][] a = [.. design.Terms.Select(t => { var c = new double[n]; model.Basis.CaseColumn(t, model.Inputs.Cases, c); return c; })];
        double[][] h = [.. design.Terms.Select(t => { var c = new double[m]; model.Basis.OutputColumn(t, model.Inputs.Outputs!, c); return c; })];
        Assert.All(a[1], v => Assert.Equal(1.0, v));
        for (int i = 0; i < n; i++)
            for (int r = 0; r < m; r++)
            {
                double viaCentred = centred[0], viaDraw = draw.Coefficients[0];
                for (int j = 0; j < 3; j++) viaCentred += centred[j + 1] * (a[j][i] * h[j][r] - a[j].Average() * h[j].Average());
                for (int q = 0; q < 3; q++)
                {
                    int j = design.Terms.ToList().IndexOf(table[draw.Structure.TermIds[q]]);
                    viaDraw += draw.Coefficients[q + 1] * a[j][i] * h[j][r];
                }
                Assert.True(Math.Abs(viaCentred - viaDraw) <= 1e-12 * Math.Max(1, Math.Abs(viaCentred)), $"({i}, {r}): {viaDraw:R} against {viaCentred:R}");
            }
    }

    [Fact]
    public void Jeffreys_RefusesExactFitsInTheChain()
    {
        // y = 1 + 2(x₀ − ½)₊ exactly: under Jeffreys every structure that spans it lies outside F, so the chain refuses
        // those candidates by reason and never holds one; the scaled default keeps them, with finite marginals.
        BassData data = Data(large: true, exact: true);
        BassConfig jeffreys = Config(maxBasis: 4) with { Sigma2 = Sigma2Prior.Jeffreys };
        var model = BassModel.Create(jeffreys, data);
        var chain = new BassChain(model, [], Math.Log(0.02), new Xoshiro256PlusPlus(21), new Xoshiro256PlusPlus(22));
        for (int i = 0; i < 4_000; i++)
            Assert.Null(model.Target.Refusal(chain.Step().Fit));
        long refused = chain.Structure.Refusals.Count(BassRefusal.ExactFit);
        _out.WriteLine($"Jeffreys: {refused} exact fits refused; final M = {chain.Current.Count}");
        Assert.True(refused > 0);

        var scaled = BassModel.Create(Config(maxBasis: 4), data);
        var kept = new BassChain(scaled, [], Math.Log(0.02), new Xoshiro256PlusPlus(21), new Xoshiro256PlusPlus(22));
        double smallest = double.PositiveInfinity;
        for (int i = 0; i < 4_000; i++) smallest = Math.Min(smallest, kept.Step().Fit.ResidualSquares);
        Assert.Equal(0, kept.Structure.Refusals.Count(BassRefusal.ExactFit));
        Assert.True(smallest < 1e-20, $"the scaled chain's smallest ‖r‖² is {smallest:E2}");
    }

    [Fact]
    public void Compose_AlignsIdsAndUncentredCoefficients()
    {
        BassData data = Data(large: true);
        var model = BassModel.Create(Config(maxBasis: 10), data);
        // Terms given out of id order, so the draw must permute them.
        BassTerm[] structure = [new BassTerm([new HingeFactor(1, -1, 3)]), new BassTerm([new HingeFactor(0, 1, 2)]),
                                new BassTerm([new HingeFactor(0, -1, 1), new HingeFactor(1, 1, 1)])];
        var moves = new BassStructureMoves(model, structure, Math.Log(0.03));
        var table = new BassTermTable();
        table.Intern(structure[2]);
        table.Intern(structure[0]);
        BassDraw draw = moves.Compose(moves.Start, new Xoshiro256PlusPlus(77), table);
        Assert.Equal([0, 1, 2], draw.Structure.TermIds);
        Assert.Equal(4, draw.Coefficients.Length);
        Assert.Equal(0.03, draw.Tau, 15);

        // The same draw in the centred basis, by the block on the same stream: α + Z_c β_c must equal a₀ + Σ aⱼBⱼ.
        BassDesign design = moves.DesignAt(moves.Start);
        CollapsedGPrior block = CollapsedGPrior.FromCentredLeastSquares(design.Factor, design.Coefficients.ToArray(), design.ResidualSquares,
                                                                         design.FittedSquares, design.Observations, design.ResponseMean, GPriorCentre.Mean);
        var centred = new double[4];
        double sigma2 = block.At(moves.Start.Tau).Sample(model.Resolved.Sigma2, new Xoshiro256PlusPlus(77), centred);
        Assert.Equal(sigma2, draw.Sigma2);
        int n = model.Resolved.Observations;
        var columns = design.Terms.Select(t => { var c = new double[n]; model.Basis.CaseColumn(t, model.Inputs.Cases, c); return c; }).ToArray();
        for (int i = 0; i < n; i++)
        {
            double viaCentred = centred[0], viaDraw = draw.Coefficients[0];
            for (int j = 0; j < 3; j++) viaCentred += centred[j + 1] * (columns[j][i] - columns[j].Average());
            for (int q = 0; q < 3; q++)
            {
                BassTerm t = table[draw.Structure.TermIds[q]];
                int j = design.Terms.ToList().IndexOf(t);
                viaDraw += draw.Coefficients[q + 1] * columns[j][i];
            }
            Assert.True(Math.Abs(viaCentred - viaDraw) <= 1e-12 * Math.Max(1, Math.Abs(viaCentred)), $"case {i}: {viaDraw:R} against {viaCentred:R}");
        }
        Assert.True(draw.Lambda > 0);
    }
}
