using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Maths.Samplers.FlatHistogram;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The structure chain's true spectral gap on the enumerated fixtures (the spectral-gap brief; gap G1).
/// <list type="bullet">
/// <item>Always run: the solver against closed forms and the dense solver, and SC-2 and the brief's checks on a small
/// enumerated kernel.</item>
/// <item>Qualification runs (<see cref="BassQualification.Enabled"/>): each move kernel and mix on the base fixtures; the
/// ladders; R̂ and the flat histogram against the exact chain. They read, they do not judge, beyond the brief's checks.
/// M ≤ 4 also needs <c>BASS_GAP_WIDE</c> = 1. Rows go to <c>BASS_GAP_OUT</c>'s folder, or <c>bass-gap</c> under the
/// temporary folder.</item>
/// </list>
/// </summary>
public sealed class BassKernelSpectrumTests
{
    private readonly ITestOutputHelper _out;
    public BassKernelSpectrumTests(ITestOutputHelper output) => _out = output;

    // ---- The solver against closed forms ----

    [Theory]
    [InlineData(0.3, 0.1)]
    [InlineData(1e-40, 1e-60)]
    public void TwoStates_ReadTheirRatesSum(double a, double b)
    {
        double logTotal = LogMath.Add(Math.Log(a), Math.Log(b));
        var flows = new Dictionary<(int, int), double> { [(0, 1)] = Math.Log(a), [(1, 0)] = Math.Log(b) };
        var chain = ReversibleChain.Of([Math.Log(b) - logTotal, Math.Log(a) - logTotal], [flows]);
        BassKernelSpectrum.Reading reading = BassKernelSpectrum.Read(chain);
        Assert.Equal(a + b, reading.Gap!.Value, (a + b) * 1e-10);
    }

    [Fact]
    public void BirthDeath_ReadsItsClosedForm()
    {
        const int n = 40;
        var flows = new Dictionary<(int, int), double>();
        for (int i = 0; i + 1 < n; i++) flows[(i, i + 1)] = flows[(i + 1, i)] = Math.Log(0.25);
        var chain = ReversibleChain.Of(new double[n], [flows]);
        BassKernelSpectrum.Reading reading = BassKernelSpectrum.Read(chain);
        // The path's lazy walk: λ_k = ½ + ½cos(πk/n).
        Assert.Equal(0.5 * (1 - Math.Cos(Math.PI / n)), reading.Gap!.Value, 1e-10);
        Assert.Equal(0.5 * (1 - Math.Cos(2 * Math.PI / n)), reading.NextGap!.Value, 1e-8);
    }

    // Two wells of two states joined by a rate ε near 10⁻⁴⁰, far under the dense solver's reach: the antisymmetric sector
    // gives γ = ε/(½ + ε + √(¼ + ε²)), with no subtraction. The set {0, 1} bounds it at L = 1/ε (SC-2).
    [Fact]
    public void TwoWells_ReadATinyGap()
    {
        const double epsilon = 1e-40;
        var flows = new Dictionary<(int, int), double>
        {
            [(0, 1)] = Math.Log(0.5), [(1, 0)] = Math.Log(0.5), [(2, 3)] = Math.Log(0.5), [(3, 2)] = Math.Log(0.5),
            [(1, 2)] = Math.Log(epsilon), [(2, 1)] = Math.Log(epsilon),
        };
        var chain = ReversibleChain.Of(new double[4], [flows]);
        BassKernelSpectrum.Reading reading = BassKernelSpectrum.Read(chain);
        double exact = epsilon / (0.5 + epsilon + Math.Sqrt(0.25 + epsilon * epsilon));
        Assert.True(Math.Abs(reading.Gap!.Value / exact - 1) <= 1e-10, $"γ reads {reading.Gap} against {exact}");

        var survey = new SetSurvey(chain);
        survey.Read("well", "{0, 1}", [0, 1]);
        survey.ReadMergeTree();
        survey.ReadEigenvectorCuts(reading.Vector!);
        Assert.Equal(-Math.Log(epsilon), survey.Bottleneck!.LogBound, 1e-9);
        Assert.True(survey.Bottleneck.LogBound <= -Math.Log(reading.Gap!.Value) + 1e-9);
    }

    // A random reversible chain whose π spans 30 nats: the elimination against the dense solver, τ_int against the
    // dense spectrum, the law against matrix powers, and SC-2 on random sets.
    [Fact]
    public void Elimination_AgreesWithTheDenseSolver()
    {
        const int n = 60;
        var rng = new Xoshiro256PlusPlus(71);
        double[] logPi = [.. Enumerable.Range(0, n).Select(_ => 30 * rng.NextDouble())];
        var flows = new Dictionary<(int, int), double>();
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
            {
                if (j > i + 1 && rng.NextDouble() > 0.15) continue;
                double q = Math.Log((0.1 + 0.9 * rng.NextDouble()) / n);
                flows[(i, j)] = q + Math.Min(0.0, logPi[j] - logPi[i]);
                flows[(j, i)] = q + Math.Min(0.0, logPi[i] - logPi[j]);
            }
        var chain = ReversibleChain.Of(logPi, [flows]);
        BassKernelSpectrum.Reading reading = BassKernelSpectrum.Read(chain);
        (double lambda2, double lambdaMin, double topOff, double vectorOff) = BassKernelSpectrum.Dense(chain);
        Assert.True(topOff <= 1e-10 && vectorOff <= 1e-8, $"the top eigenpair sits {topOff:E2}, {vectorOff:E2} off");
        Assert.True(lambdaMin > -1.0);
        Assert.True(Math.Abs(reading.Gap!.Value / (1 - lambda2) - 1) <= 1e-6, $"{reading.Gap} against {1 - lambda2}");

        double[] f = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        EigenResult eigen = DenseEigen.DecomposeSymmetric(chain.Symmetrized(), 256, 1e-14);
        double[] pi = chain.Pi();
        double mean = Enumerable.Range(0, n).Sum(i => pi[i] * f[i]), variance = 0.0, weighted = 0.0;
        for (int k = 1; k < n; k++)
        {
            double c = Enumerable.Range(0, n).Sum(i => eigen.Eigenvectors[k][i] * Math.Sqrt(pi[i]) * (f[i] - mean));
            variance += c * c;
            weighted += c * c * (1 + eigen.Eigenvalues[k]) / (1 - eigen.Eigenvalues[k]);
        }
        Assert.Equal(weighted / variance, reading.Solver!.IntegratedTime(f), weighted / variance * 1e-6);

        double[] start = new double[n];
        start[0] = 1;
        (double[] tv, _) = BassKernelSpectrum.Law(chain, start, [40], null);
        double[,] s = chain.Symmetrized();
        var power = new double[n];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++)
                power[i] += Math.Pow(eigen.Eigenvalues[k], 40) * eigen.Eigenvectors[k][0] * eigen.Eigenvectors[k][i] * Math.Sqrt(pi[i] / pi[0]);
        Assert.Equal(0.5 * Enumerable.Range(0, n).Sum(i => Math.Abs(power[i] - pi[i])), tv[0], 1e-8);

        var survey = new SetSurvey(chain);
        for (int trial = 0; trial < 2_000; trial++)
            survey.Read("random", "", [.. Enumerable.Range(0, n).Where(_ => rng.NextDouble() < 0.5)]);
        survey.ReadMergeTree();
        survey.ReadEigenvectorCuts(reading.Vector!);
        Assert.True(survey.Bottleneck!.LogBound <= -Math.Log(reading.Gap!.Value) + 1e-6, $"L(A) reads {survey.Bottleneck.LogBound} over ln t_rel {-Math.Log(reading.Gap!.Value)}");
    }

    // Two modes, 0 and the likelier 1, with a saddle 2 between: 2 joins 0 at conductance e⁻³², so the saddle's depth is 2;
    // then {0, 2} joins 1 at e⁻⁴⁰, so mode 0's depth is 30 (masses and conductances shift together under normalization).
    [Fact]
    public void ModeDepths_ReadTheLikelierTopsSurvival()
    {
        double[] logPi = [-10.0, 0.0, -30.0];
        var flows = new Dictionary<(int, int), double>
        {
            [(0, 2)] = -22.0, [(2, 0)] = -22.0 + logPi[0] - logPi[2],
            [(1, 2)] = -40.0, [(2, 1)] = -40.0 + logPi[1] - logPi[2],
        };
        var chain = ReversibleChain.Of(logPi, [flows]);
        (Dictionary<int, (double Depth, int Mode)> deepest, int modes, double[] depthOf) = BassKernelSpectrum.ModeDepths(chain, [1, 0, 2]);
        Assert.Equal(2, modes);
        Assert.Equal(0, deepest[1].Mode);
        Assert.Equal(30.0, deepest[1].Depth, 6);
        Assert.Equal(2.0, deepest[2].Depth, 6);
        Assert.False(deepest.ContainsKey(0));   // the top never joins a likelier one
        Assert.True(double.IsNaN(depthOf[1]));
        Assert.Equal(30.0, depthOf[0], 6);
        Assert.Equal(2.0, depthOf[2], 6);
    }

    [Fact]
    public void Scaling_ReadsItsRule()
    {
        double[] n = [18, 36, 72, 144, 288, 576];
        Assert.Equal("linear", BassKernelSpectrum.Scaling(n, [.. n.Select(v => 3 + 0.2 * v)]).Reading);
        Assert.Equal("logarithmic", BassKernelSpectrum.Scaling(n, [.. n.Select(v => 3 + 0.5 * Math.Log(v))]).Reading);
        Assert.Equal("flat", BassKernelSpectrum.Scaling(n, [.. n.Select(v => 3 + 0.001 * Math.Log(v))]).Reading);
        Assert.Equal("neither", BassKernelSpectrum.Scaling(n, [3, 5, 4, 8, 6, 7]).Reading);
    }

    // ---- The fixtures ----

    internal const string InsideF = "inside F", OutsideF = "outside F";

    /// <summary>The L × L grid repeated <paramref name="reps"/> times, the noise one stream, so a rung extends the one
    /// below and at L = 3, r = 2, seed 5 inside F is <see cref="BassSweepTests.Data"/>.</summary>
    internal static BassData Grid(string surface, int levels, int reps, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
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
                    double mean = surface == InsideF
                        ? 1 + 2 * a + 3 * Math.Max(a - 0.5, 0) * Math.Max(0.5 - b, 0)
                        : 2 * Math.Sin(Math.PI * a * b) + 4 * (b - 0.5) * (b - 0.5) + a;
                    y.Add(mean + 0.3 * BassQualification.Gaussian(rng));
                }
        return new BassData([new ContinuousPredictor("x0", [.. x0]), new ContinuousPredictor("x1", [.. x1])], new ScalarResponse([.. y]));
    }

    // BassKernelEnumeration.CategoricalData, its 24 cases repeated r times.
    private static BassData CategoricalGrid(int reps, int seed = 17)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        int n = 24 * reps;
        var x = new double[n];
        var z1 = new int[n];
        var z2 = new int[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i % 3) / 2.0;
            z1[i] = (i / 3) % 2;
            z2[i] = (i / 6) % 4;
            y[i] = 1 + x[i] + (z2[i] <= 1 ? 1.0 : 0.0) + 0.5 * (z1[i] == 0 ? x[i] : 0.0) + 0.3 * BassQualification.Gaussian(rng);
        }
        return new BassData([new ContinuousPredictor("x0", x), new CategoricalPredictor("z1", z1, ["a", "b"]),
                             new CategoricalPredictor("z2", z2, ["p", "q", "r", "s"])], new ScalarResponse(y));
    }

    // BassKernelEnumeration.FunctionalData, its six cases repeated r times.
    private static BassData FunctionalGrid(int reps, int seed = 19)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] t1 = [0, 0.5, 1, 0.5], t2 = [0, 0, 1, 1];
        int n = 6 * reps;
        var x = new double[n];
        var y = new double[n, t1.Length];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i % 3) / 2.0;
            for (int r = 0; r < t1.Length; r++) y[i, r] = 1 + x[i] + t1[r] * x[i] + 0.5 * t2[r] + 0.2 * BassQualification.Gaussian(rng);
        }
        return new BassData([new ContinuousPredictor("x0", x)],
                            new FunctionalResponse(y, [new ContinuousPredictor("t1", t1), new ContinuousPredictor("t2", t2)]));
    }

    private static BassCoefficientPrior HeldUnitInformation => BassCoefficientPrior.G(TauPrior.Held());

    // Ladder R holds the support threshold at 0.1·n, past R's cap of 20, so F stays put as n grows.
    private static BassConfig GridConfig(int cap, double? support) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = cap, Support = new SupportRule(Minimum: support) },
        Coefficients = HeldUnitInformation,
    };

    private static BassConfig CategoricalConfig(int cases, CategoricalSubsetMeasure measure) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxCategoricalInteraction = 1, MaxBasis = 2, CategoricalSubsets = measure, Support = new SupportRule(Minimum: 0.1 * cases) },
        Coefficients = HeldUnitInformation,
    };

    private static BassConfig FunctionalConfig(int cases) => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 1, MaxFunctionalInteraction = 1, MaxBasis = 2, Support = new SupportRule(Minimum: 0.1 * cases) },
        Coefficients = HeldUnitInformation,
    };

    /// <summary>One enumerated kernel: its model, F, the flows by group, the whole mix, and the key-basin names.</summary>
    private sealed record Built(string Name, BassModel Model, BassConfig Config, BassData Data, BassKernelEnumeration.Space Space,
                                Dictionary<string, Dictionary<(int From, int To), double>> Flows, ReversibleChain Mix, BassMixing.Fixture Named)
    {
        public double[] LogPi => [.. Space.States.Select(s => s.LogPi)];
        public string Basin(int state) => Named.Basin(Space.States[state].Terms);
        public string Terms(int state) => BassMixing.Key(Space.States[state].Terms);

        // ln of the structure's prior weight on F: the count weight and the terms' measure.
        public double LogPrior(int state)
        {
            BassTerm[] held = Space.States[state].Terms;
            return Model.Prior.LogCollectionWeight(held.Length) + held.Sum(Model.Prior.LogPrior);
        }

        // The prior on F, the run's start.
        public double[] Prior()
        {
            double[] log = [.. Enumerable.Range(0, Space.States.Count).Select(LogPrior)];
            double total = LogMath.Sum(log);
            return [.. log.Select(v => Math.Exp(v - total))];
        }

        // SC-4's slope: the smallest over the state's terms of [ln m(S) − ln m(S without it)] / N; NaN for the intercept.
        public double LoadBearingSlope(int state, double n)
        {
            BassTerm[] held = Space.States[state].Terms;
            if (held.Length == 0) return double.NaN;
            double logM = Space.States[state].LogPi - LogPrior(state), least = double.PositiveInfinity;
            foreach (BassTerm t in held)
            {
                int without = Space.Find(held.Where(u => !u.Equals(t)));
                Assert.True(without >= 0, "F is downward closed");
                least = Math.Min(least, logM - (Space.States[without].LogPi - LogPrior(without)));
            }
            return least / n;
        }
    }

    private static Built Build(string name, BassConfig config, BassData data, LocalSteps? steps = null)
    {
        var model = BassModel.Create(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        Dictionary<string, Dictionary<(int From, int To), double>> flows =
            BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.All(space, key, s, steps));
        BassKernelEnumeration.Check(space, flows);
        var mix = ReversibleChain.Of([.. space.States.Select(s => s.LogPi)], flows.Values);
        return new Built(name, model, config, data, space, flows, mix, new BassMixing.Fixture(name, data, config, null));
    }

    private static SetSurvey Survey(Built built, double[]? vector)
    {
        var survey = new SetSurvey(built.Mix);
        foreach (IGrouping<string, int> basin in Enumerable.Range(0, built.Mix.Size).GroupBy(built.Basin)) survey.Read("key basin", basin.Key, [.. basin]);
        foreach (IGrouping<int, int> level in Enumerable.Range(0, built.Mix.Size).GroupBy(i => built.Space.States[i].Terms.Length))
            survey.Read("M level", $"M = {level.Key}", [.. level]);
        survey.ReadMergeTree();
        if (vector is not null) survey.ReadEigenvectorCuts(vector);
        return survey;
    }

    // ---- A small enumerated kernel, always ----

    // The base grid at M ≤ 2: both solvers, SC-2 over every family, Dirichlet forms adding, and F held along the replicates.
    [Fact]
    public void SmallKernel_HoldsTheBriefsChecks()
    {
        Built built = Build("continuous, M ≤ 2", GridConfig(2, 1.8), Grid(InsideF, 3, 2, 5));
        Assert.Single(built.Mix.Classes());
        BassKernelSpectrum.Reading reading = BassKernelSpectrum.Read(built.Mix);
        (double lambda2, double lambdaMin, double topOff, double vectorOff) = BassKernelSpectrum.Dense(built.Mix);
        Assert.True(topOff <= 1e-10 && vectorOff <= 1e-8, $"the top eigenpair sits {topOff:E2}, {vectorOff:E2} off");
        Assert.True(Math.Abs(reading.Gap!.Value / (1 - lambda2) - 1) <= 1e-6, $"{reading.Gap} against {1 - lambda2}");
        SetSurvey survey = Survey(built, reading.Vector);
        Assert.True(survey.Bottleneck!.LogBound <= -Math.Log(reading.Gap!.Value) + 1e-6);
        double sum = 0.0;
        foreach ((string group, Dictionary<(int From, int To), double> flows) in built.Flows)
        {
            var kernel = ReversibleChain.Of(built.LogPi, [flows]);
            if (kernel.Classes().Length == 1) sum += BassKernelSpectrum.Read(kernel).Gap!.Value;
        }
        Assert.True(reading.Gap!.Value >= sum * (1 - 1e-8));
        _out.WriteLine($"|F| {built.Mix.Size}; γ {reading.Gap:E4}; λ_min {lambdaMin:F4}; bottleneck {survey.Bottleneck.Family}, L {Math.Exp(survey.Bottleneck.LogBound):E3}");
        Assert.Equal(built.Mix.Size, Build("continuous, M ≤ 2, r = 4", GridConfig(2, 3.6), Grid(InsideF, 3, 4, 5)).Mix.Size);
        Assert.Equal(BassSweepTests.Data().Response.Observations, Grid(InsideF, 3, 2, 5).Response.Observations);
    }

    // ---- Each move kernel and mix on the base fixtures ----

    private static readonly (string Mix, string Schedule, bool Steps, bool Family)[] Mixes =
    [
        ("default", "equal", false, false), ("cancelling", "cancelling", false, false),
        ("steps", "equal", true, false), ("family", "equal", false, true), ("all", "equal", true, true),
    ];

    // The dense solver takes about 19 minutes at 2,263 states, so past MixDenseStates only the default mix pays for it.
    private const int DenseStates = 2_500, MixDenseStates = 1_000, ClassDenseStates = 400;

    [Theory]
    [InlineData("continuous")]
    [InlineData("categorical, R's measure")]
    [InlineData("categorical, uniform")]
    [InlineData("functional")]
    public void Kernels_OnTheBaseFixtures(string fixture)
    {
        if (!BassQualification.Enabled) return;
        foreach ((string mix, string schedule, bool withSteps, bool family) in Mixes)
        {
            LocalSteps? steps = withSteps ? new LocalSteps() : null;
            (BassConfig config, BassData data) = BassLocalStepTests.Fixture(fixture, schedule, s => s with { LocalSteps = steps, FactorMoves = family ? new FactorMoves() : null });
            Built built = Build(fixture, config, data, steps);
            Assert.Single(built.Mix.Classes());
            BassKernelSpectrum.Reading reading = BassKernelSpectrum.Read(built.Mix, Threads);
            double gap = reading.Gap!.Value;
            Assert.True(reading.Residual <= 1e-6, $"residual {reading.Residual:E2}");
            string denseGap = "", lambdaMin = "";
            if (built.Mix.Size <= MixDenseStates || built.Mix.Size <= DenseStates && mix == "default")
            {
                (double lambda2, double min, double topOff, double vectorOff) = BassKernelSpectrum.Dense(built.Mix);
                Assert.True(topOff <= 1e-10 && vectorOff <= 1e-8, $"{fixture}, {mix}: the top eigenpair sits {topOff:E2}, {vectorOff:E2} off");
                if (gap >= BassKernelSpectrum.DenseFloor) Assert.True(Math.Abs(gap / (1 - lambda2) - 1) <= 1e-6, $"{fixture}, {mix}: {gap} against {1 - lambda2}");
                denseGap = F(1 - lambda2);
                lambdaMin = F(min);
            }
            SetSurvey survey = Survey(built, reading.Vector);
            Assert.True(survey.Bottleneck!.LogBound <= -Math.Log(gap) + 1e-6, $"{fixture}, {mix}: SC-2 fails");
            Row("kernels", KernelHeader, [fixture, mix, "mix", I(built.Mix.Size), I(built.Mix.Transitions), "1", F(gap), F(1 / gap), reading.NextGap is { } next ? F(next) : "", denseGap, lambdaMin, "", "", "",
                                          .. SetCells(built, survey.Bottleneck, gap), .. SetCells(built, survey.RelevantBottleneck, gap), I(survey.Tried)]);
            _out.WriteLine($"{fixture}, {mix}: |F| {built.Mix.Size}; γ {gap:E4}, t_rel {1 / gap:E3}, next γ {reading.NextGap:E3}; dense {denseGap}, λ_min {lambdaMin}");
            WriteSet("bottleneck", built, survey.Bottleneck, gap);
            WriteSet("relevant", built, survey.RelevantBottleneck, gap);

            double sum = 0.0;
            foreach ((string group, Dictionary<(int From, int To), double> flows) in built.Flows)
            {
                var kernel = ReversibleChain.Of(built.LogPi, [flows]);
                int[][] classes = kernel.Classes();
                int[][] moving = [.. classes.Where(c => c.Length > 1)];
                double[] gaps = [.. moving.Select(c => BassKernelSpectrum.Read(kernel.On(c)).Gap!.Value)];
                foreach (int[] c in moving.Where(c => c.Length <= ClassDenseStates).Take(20))
                {
                    ReversibleChain on = kernel.On(c);
                    (double lambda2, _, _, _) = BassKernelSpectrum.Dense(on);
                    double elimination = BassKernelSpectrum.Read(on).Gap!.Value;
                    if (elimination >= BassKernelSpectrum.DenseFloor) Assert.True(Math.Abs(elimination / (1 - lambda2) - 1) <= 1e-6, $"{group}: {elimination} against {1 - lambda2}");
                }
                if (classes.Length == 1) sum += gaps[0];
                double[] classMass = [.. moving.Select(c => LogMath.Sum(c.Select(i => kernel.LogPi[i]).ToArray()))];
                int heaviest = moving.Length == 0 ? -1 : Array.IndexOf(classMass, classMass.Max());
                int[] heavy = heaviest < 0 ? [] : moving[heaviest];
                string partition = "";
                if (group == BassKernelEnumeration.Redraw)
                {
                    int keys = Enumerable.Range(0, built.Mix.Size).GroupBy(built.Basin).Count();
                    bool same = classes.All(c => c.Select(built.Basin).Distinct().Count() == 1) && classes.Length == keys;
                    partition = $"{classes.Length} classes, {keys} key basins, {(same ? "the same partition" : "not the same partition")}";
                    double worst = gaps.Length > 0 ? gaps.Min() : double.NaN;
                    string slowMode = gap < 0.1 * worst ? "between basins" : gap <= 2 * worst && worst <= 2 * gap ? "within" : "mixed";
                    _out.WriteLine($"  {partition}; worst class γ {worst:E3}: {slowMode}");
                    partition += $"; {slowMode}";
                }
                Row("kernels", KernelHeader, [fixture, mix, group, I(kernel.Size), I(kernel.Transitions), I(classes.Length), classes.Length == 1 ? F(gaps[0]) : "0", "", "", "", "",
                                              gaps.Length > 0 ? F(gaps.Min()) : "", heaviest >= 0 ? F(gaps[heaviest]) : "", partition, .. new string[12 + 12], I(heavy.Length)]);
                _out.WriteLine($"  {group}: {kernel.Transitions} transitions, {classes.Length} classes; " +
                               (classes.Length == 1 ? $"γ {gaps[0]:E4}" : $"worst class γ {(gaps.Length > 0 ? gaps.Min() : double.NaN):E3}, the heaviest class's {(heaviest >= 0 ? gaps[heaviest] : double.NaN):E3} on {heavy.Length} states"));
            }
            Assert.True(gap >= sum * (1 - 1e-8), $"{fixture}, {mix}: γ_mix {gap} under the groups' sum {sum}");
        }
    }

    private void WriteSet(string what, Built built, SetReading? set, double gap)
    {
        if (set is null)
        {
            _out.WriteLine($"  {what}: no set holds 1% each side");
            return;
        }
        int mode = set.Members.MaxBy(i => built.Mix.LogPi[i]);
        _out.WriteLine($"  {what}: {set.Family} {set.Label}, {set.Size} states; π(A) {Math.Exp(set.LogMass):E3}, ln flow {set.LogFlow:F2}, ln L {set.LogBound:F3}, L·γ {Math.Exp(set.LogBound) * gap:F4}, " +
                       $"D {set.Barrier:F2}; mode {built.Terms(mode)} in {built.Basin(mode)}; {Composition(built, set)}");
    }

    // The set's largest key basin, by share of its mass.
    private static string Composition(Built built, SetReading set)
    {
        (string Key, double Log) top = set.Members.GroupBy(built.Basin)
            .Select(g => (g.Key, LogMath.Sum(g.Select(i => built.Mix.LogPi[i]).ToArray()))).MaxBy(g => g.Item2);
        return $"{Math.Exp(top.Log - set.LogMass):P1} in {top.Key}";
    }

    private static string[] SetCells(Built built, SetReading? set, double gap)
    {
        if (set is null) return new string[12];
        int mode = set.Members.MaxBy(i => built.Mix.LogPi[i]);
        return [set.Family, I(set.Size), F(Math.Exp(set.LogMass)), F(set.LogMass), F(set.LogFlow), F(set.LogBound), double.IsNaN(gap) ? "" : F(Math.Exp(set.LogBound) * gap),
                F(set.Barrier), F(set.ComplementBarrier), built.Terms(mode), built.Basin(mode), Composition(built, set)];
    }

    private static readonly string[] KernelHeader =
    [
        "fixture", "mix", "group", "states", "transitions", "classes", "gap", "t_rel", "next_gap", "dense_gap", "lambda_min", "worst_class_gap", "heaviest_class_gap", "note",
        "bn_family", "bn_size", "bn_mass", "bn_ln_mass", "bn_ln_flow", "bn_ln_L", "bn_L_gap", "bn_D", "bn_D_other", "bn_mode", "bn_basin", "bn_composition",
        "rel_family", "rel_size", "rel_mass", "rel_ln_mass", "rel_ln_flow", "rel_ln_L", "rel_L_gap", "rel_D", "rel_D_other", "rel_mode", "rel_basin", "rel_composition", "sets_or_heaviest_size",
    ];

    // ---- One sweep at held τ is one step of K ----

    // The production chain's M autocovariance at lags 1, 2, 5 and 20 against K^k's: T·γ means nothing otherwise.
    [Fact]
    public void Chain_IsOneStepOfK()
    {
        if (!BassQualification.Enabled) return;
        (BassConfig config, BassData data) = BassLocalStepTests.Fixture("continuous", "equal", s => s);
        Built built = Build("continuous", config, data);
        double[] pi = built.Mix.Pi();
        double[] m = [.. built.Space.States.Select(s => (double)s.Terms.Length)];
        double mean = pi.Zip(m).Sum(p => p.First * p.Second);
        double[] centred = [.. m.Select(v => v - mean)];
        int[] lags = [1, 2, 5, 20];
        var exact = new double[lags.Length];
        double[] power = [.. centred], next = new double[power.Length];
        for (int k = 1, l = 0; l < lags.Length; k++)
        {
            built.Mix.Apply(power, next);
            (power, next) = (next, power);
            if (k != lags[l]) continue;
            for (int i = 0; i < pi.Length; i++) exact[l] += pi[i] * centred[i] * power[i];
            l++;
        }
        const int chains = 8, burn = 2_000, sweeps = 200_000;
        var estimates = new double[chains][];
        Parallel.For(0, chains, new ParallelOptions { MaxDegreeOfParallelism = Threads }, c =>
        {
            int[] seeds = SeedTree.Derive(6100 + c, 2);
            var chain = new BassChain(built.Model, [], built.Space.States[0].Loaded.LogTau, new Xoshiro256PlusPlus(seeds[0]), new Xoshiro256PlusPlus(seeds[1]));
            for (int i = 0; i < burn; i++) chain.Step();
            var trace = new double[sweeps];
            for (int i = 0; i < sweeps; i++) trace[i] = chain.Step().Count - mean;
            estimates[c] = [.. lags.Select(lag => Enumerable.Range(0, sweeps - lag).Sum(i => trace[i] * trace[i + lag]) / (sweeps - lag))];
        });
        for (int l = 0; l < lags.Length; l++)
        {
            double[] values = [.. estimates.Select(e => e[l])];
            double average = values.Average(), error = Math.Sqrt(values.Sum(v => (v - average) * (v - average)) / (chains - 1) / chains);
            _out.WriteLine($"lag {lags[l]}: chain {average:F5} ± {error:F5} against K^k's {exact[l]:F5}");
            Assert.True(Math.Abs(average - exact[l]) <= 5 * error + 1e-3, $"lag {lags[l]}: {average} ± {error} against {exact[l]}");
        }
    }

    // ---- The ladders ----

    private static readonly int[] GridReps = [2, 4, 8, 16, 32, 64, 128, 256];
    private static readonly int[] CaseReps = [1, 2, 4, 8, 16, 32, 64];
    private static readonly int[] Seeds = [5, 11, 23];

    // The exact law's clocks: burn-in's end and the last sweep at each length the brief runs.
    private static readonly int[] Clocks = [500, 700, 5_000, 7_000, 50_000, 70_000];

    /// <summary>One ladder rung's exact reading.</summary>
    private sealed record Rung(string Ladder, string Fixture, string Surface, int Cap, int Seed, int R, int N, Dictionary<string, double> ExactBasins,
                               int States, double? Gap, double? RelevantBound, double Barrier, double[] Tv, double[] BottleneckMass, double TauM,
                               Dictionary<int, double> DepthBySlack);

    private static readonly string[] LadderHeader =
    [
        "ladder", "fixture", "surface", "cap", "seed", "r", "N", "states", "transitions", "classes", "gap", "ln_gap", "t_rel", "next_gap", "residual", "top_state", "pi_top", "top_slope",
        "bn_family", "bn_size", "bn_mass", "bn_ln_mass", "bn_ln_flow", "bn_ln_L", "bn_L_gap", "bn_D", "bn_D_other", "bn_mode", "bn_basin", "bn_composition", "bn_slope",
        "rel_family", "rel_size", "rel_mass", "rel_ln_mass", "rel_ln_flow", "rel_ln_L", "rel_L_gap", "rel_D", "rel_D_other", "rel_mode", "rel_basin", "rel_composition",
        "tv_500", "tv_700", "tv_5000", "tv_7000", "tv_50000", "tv_70000", "bn_mu_500", "bn_mu_5000", "bn_mu_50000", "tau_M", "sets",
        "local_modes", "depth_s0", "depth_s0_mode", "depth_s1", "depth_s1_mode", "depth_s2", "depth_s2_mode",
    ];

    private static readonly ConcurrentDictionary<string, Rung> Rungs = new();

    private static Rung ReadRung(string ladder, string fixture, string surface, int cap, int seed, int r, int n, BassConfig config, BassData data, int threads, bool farClocks)
    {
        return Rungs.GetOrAdd($"{ladder}|{fixture}|{surface}|{cap}|{seed}|{r}", _ =>
        {
            Built built = Build(fixture, config, data);
            int[][] classes = built.Mix.Classes();
            Assert.Single(classes);
            BassKernelSpectrum.Reading reading = BassKernelSpectrum.Read(built.Mix, threads);
            double gap = reading.Gap ?? double.NaN;
            if (reading.Gap is not null) Assert.True(reading.Residual <= 1e-6, $"{fixture}, r = {r}: residual {reading.Residual:E2}");
            SetSurvey survey = Survey(built, reading.Vector);
            if (reading.Gap is not null) Assert.True(survey.Bottleneck!.LogBound <= -Math.Log(gap) + 1e-6, $"{fixture}, r = {r}: SC-2 fails, ln L {survey.Bottleneck.LogBound} over ln t_rel {-Math.Log(gap)}");
            int[] clocks = farClocks ? Clocks : [.. Clocks.Where(c => c <= 7_000)];
            SetReading bottleneck = survey.Bottleneck!;
            (double[] tv, double[] mass) = BassKernelSpectrum.Law(built.Mix, built.Prior(), clocks, bottleneck.Members);
            double[] m = [.. built.Space.States.Select(s => (double)s.Terms.Length)];
            double tauM = reading.Solver?.IntegratedTime(m) ?? double.NaN;
            int top = Enumerable.Range(0, built.Mix.Size).MaxBy(i => built.Mix.LogPi[i]);
            int bnMode = bottleneck.Members.MaxBy(i => built.Mix.LogPi[i]);
            string Clock(double[] values, int at) => Array.IndexOf(clocks, at) is >= 0 and var k ? F(values[k]) : "";
            int[] slack = [.. built.Space.States.Select(s => built.Model.MaxBasis - s.Terms.Length)];
            var depths = BassKernelSpectrum.ModeDepths(built.Mix, slack);
            Dictionary<int, (double Depth, int Mode)> deepest = depths.Deepest;
            int localModes = depths.LocalModes;
            IEnumerable<string> DepthCells(int s) => deepest.TryGetValue(s, out (double Depth, int Mode) d) ? [F(d.Depth), built.Terms(d.Mode)] : ["", ""];
            string[] cells =
            [
                ladder, fixture, surface, I(cap), I(seed), I(r), I(n), I(built.Mix.Size), I(built.Mix.Transitions), I(classes.Length),
                reading.Gap is null ? "" : F(gap), reading.Gap is null ? "" : F(Math.Log(gap)), reading.Gap is null ? "" : F(1 / gap), reading.NextGap is { } next ? F(next) : "", F(reading.Residual),
                built.Terms(top), F(Math.Exp(built.Mix.LogPi[top])), F(built.LoadBearingSlope(top, n)),
                .. SetCells(built, bottleneck, gap), F(built.LoadBearingSlope(bnMode, n)),
                .. SetCells(built, survey.RelevantBottleneck, gap),
                .. Clocks.Select(c => Clock(tv, c)), Clock(mass, 500), Clock(mass, 5_000), Clock(mass, 50_000), F(tauM), I(survey.Tried),
                I(localModes), .. DepthCells(0), .. DepthCells(1), .. DepthCells(2),
            ];
            Row("ladder", LadderHeader, cells);
            var exactBasins = new Dictionary<string, double>();
            for (int i = 0; i < built.Mix.Size; i++)
            {
                string key = built.Basin(i);
                exactBasins[key] = exactBasins.GetValueOrDefault(key) + Math.Exp(built.Mix.LogPi[i]);
            }
            double[] clockedTv = [.. Clocks.Select(c => Array.IndexOf(clocks, c) is >= 0 and var k ? tv[k] : double.NaN)];
            double[] clockedMass = [.. Clocks.Select(c => Array.IndexOf(clocks, c) is >= 0 and var k ? mass[k] : double.NaN)];
            return new Rung(ladder, fixture, surface, cap, seed, r, n, exactBasins, built.Mix.Size, reading.Gap,
                            survey.RelevantBottleneck is { } rel ? Math.Exp(rel.LogBound) : null, bottleneck.Barrier, clockedTv, clockedMass, tauM,
                            deepest.ToDictionary(d => d.Key, d => d.Value.Depth));
        });
    }

    private static (string Fixture, string Surface, int Cap, int Seed, int R, int N, BassConfig Config, BassData Data) GridRung(string surface, int cap, int seed, int r)
    {
        int n = 9 * r;
        return ($"grid M ≤ {cap}", surface, cap, seed, r, n, GridConfig(cap, 0.1 * n), Grid(surface, 3, r, seed));
    }

    public static IEnumerable<object[]> GridCells()
    {
        foreach (string surface in new[] { InsideF, OutsideF })
            foreach (int cap in new[] { 2, 3 })
                foreach (int seed in Seeds) yield return [surface, cap, seed];
    }

    [Theory]
    [MemberData(nameof(GridCells))]
    public void LadderR_OnTheGrid(string surface, int cap, int seed)
    {
        if (!BassQualification.Enabled) return;
        Rung[] rungs = ReadLadder("R", GridReps, r => GridRung(surface, cap, seed, r), seed == 5);
        Assert.All(rungs, rung => Assert.Equal(rungs[0].States, rung.States));
        ReadCell("R", rungs);
    }

    [Theory]
    [InlineData("categorical, R's measure")]
    [InlineData("categorical, uniform")]
    [InlineData("functional")]
    public void LadderR_OnTheCases(string fixture)
    {
        if (!BassQualification.Enabled) return;
        Rung[] rungs = ReadLadder("R", CaseReps, r => CaseRung(fixture, r), true);
        Assert.All(rungs, rung => Assert.Equal(rungs[0].States, rung.States));
        ReadCell("R", rungs);
    }

    private static (string Fixture, string Surface, int Cap, int Seed, int R, int N, BassConfig Config, BassData Data) CaseRung(string fixture, int r) => fixture switch
    {
        "functional" => (fixture, "as built", 2, 19, r, 24 * r, FunctionalConfig(6 * r), FunctionalGrid(r)),
        "categorical, uniform" => (fixture, "as built", 2, 17, r, 24 * r, CategoricalConfig(24 * r, CategoricalSubsetMeasure.Uniform), CategoricalGrid(r)),
        _ => (fixture, "as built", 2, 17, r, 24 * r, CategoricalConfig(24 * r, CategoricalSubsetMeasure.SizeUniform), CategoricalGrid(r)),
    };

    // M ≤ 4 on seed 5 at r = 2, 8, 32, 128, heavy: its own switch inside the qualification.
    [Theory]
    [InlineData(InsideF)]
    [InlineData(OutsideF)]
    public void LadderR_AtTheWiderCap(string surface)
    {
        if (!BassQualification.Enabled || Environment.GetEnvironmentVariable("BASS_GAP_WIDE") is not "1") return;
        var rungs = new List<Rung>();
        foreach (int r in new[] { 2, 8, 32, 128 })
        {
            (string fixture, string s, int cap, int seed, int rr, int n, BassConfig config, BassData data) = GridRung(surface, 4, 5, r);
            rungs.Add(ReadRung("R", fixture, s, cap, seed, rr, n, config, data, Threads, false));
            _out.WriteLine(Line(rungs[^1]));
        }
        ReadCell("R", [.. rungs]);
    }

    private static readonly string[] CapDepthHeader = ["surface", "seed", "r", "N", "cap", "states", "free", "deepest", "mode", "tracked_from", "tracked_M", "depth"];

    // Each mode's depth at caps 3, 4 and 5 on the same data, from the kernel's edges alone, so M ≤ 5 is within reach. The
    // modes the narrower caps trap are followed up: M ≤ 4's mode with one free slot, which M ≤ 5 gives two.
    [Theory]
    [InlineData(InsideF)]
    [InlineData(OutsideF)]
    public void Depths_AcrossCaps(string surface)
    {
        if (!BassQualification.Enabled || Environment.GetEnvironmentVariable("BASS_GAP_WIDE") is not "1") return;
        const int seed = 5;
        Parallel.ForEach(new[] { 2, 8, 32, 128, 512, 2048 }, new ParallelOptions { MaxDegreeOfParallelism = Math.Min(Threads, 4) }, r =>
        {
            int n = 9 * r;
            BassData data = Grid(surface, 3, r, seed);
            var followed = new List<(string From, BassTerm[] Terms)>();
            foreach (int cap in new[] { 3, 4, 5 })
            {
                Built built = Build($"grid M ≤ {cap}", GridConfig(cap, 0.1 * n), data);
                int[] slack = [.. built.Space.States.Select(s => cap - s.Terms.Length)];
                (Dictionary<int, (double Depth, int Mode)> deepest, _, double[] depthOf) = BassKernelSpectrum.ModeDepths(built.Mix, slack);
                string[] Cells(string free, double value, string mode, string from, string m) => [surface, I(seed), I(r), I(n), I(cap), I(built.Mix.Size), free, F(value), mode, from, m, ""];
                foreach ((int free, (double depth, int mode)) in deepest.OrderBy(d => d.Key)) Row("capdepths", CapDepthHeader, Cells(I(free), depth, built.Terms(mode), "", ""));
                foreach ((string from, BassTerm[] terms) in followed)
                {
                    int at = built.Space.Find(terms);
                    Row("capdepths", CapDepthHeader, [surface, I(seed), I(r), I(n), I(cap), I(built.Mix.Size), I(cap - terms.Length), "", BassMixing.Key(terms), from, I(terms.Length), at >= 0 ? F(depthOf[at]) : "refused"]);
                }
                foreach (int free in new[] { 0, 1 }.Where(deepest.ContainsKey))
                    if (cap < 5) followed.Add(($"cap {cap}, {free} free", built.Space.States[deepest[free].Mode].Terms));
            }
        });
    }

    public static IEnumerable<object[]> LevelCells()
    {
        foreach (string surface in new[] { InsideF, OutsideF })
            foreach (int seed in Seeds) yield return [surface, seed];
    }

    // Ladder L: the grid's side grows at two replicates, M ≤ 2, the default support rule; F grows.
    [Theory]
    [MemberData(nameof(LevelCells))]
    public void LadderL_OnTheGrid(string surface, int seed)
    {
        if (!BassQualification.Enabled) return;
        var rungs = new List<Rung>();
        foreach (int levels in new[] { 3, 4, 5, 6 })
        {
            int n = 2 * levels * levels;
            rungs.Add(ReadRung("L", "grid M ≤ 2, L", surface, 2, seed, levels, n, GridConfig(2, null), Grid(surface, levels, 2, seed), Threads, false));
            _out.WriteLine(Line(rungs[^1]));
        }
        ReadCell("L", [.. rungs]);
    }

    private Rung[] ReadLadder(string ladder, int[] reps, Func<int, (string Fixture, string Surface, int Cap, int Seed, int R, int N, BassConfig Config, BassData Data)> rung, bool farClocks)
    {
        var read = new Rung[reps.Length];
        Parallel.For(0, reps.Length, new ParallelOptions { MaxDegreeOfParallelism = Threads }, k =>
        {
            (string fixture, string surface, int cap, int seed, int r, int n, BassConfig config, BassData data) = rung(reps[k]);
            read[k] = ReadRung(ladder, fixture, surface, cap, seed, r, n, config, data, 1, farClocks);
        });
        foreach (Rung r in read) _out.WriteLine(Line(r));
        return read;
    }

    private static string Line(Rung r)
        => $"{r.Ladder} {r.Fixture}, {r.Surface}, seed {r.Seed}: r {r.R}, N {r.N}, |F| {r.States}; γ {(r.Gap is { } g ? g.ToString("E3", CultureInfo.InvariantCulture) : "below range")}; D {r.Barrier:F2}; relevant L {(r.RelevantBound is { } l ? l.ToString("E3", CultureInfo.InvariantCulture) : "none")}; TV at 5,000 {r.Tv[2]:F4}; τ_int(M) {r.TauM:E3}";

    // The cell's scaling: −ln γ, D of the bottleneck, and ln L of the relevant bottleneck, each by the brief's rule.
    private void ReadCell(string ladder, Rung[] rungs)
    {
        void One(string what, Func<Rung, double?> value)
        {
            Rung[] held = [.. rungs.Where(r => value(r) is not null)];
            string cell = $"{rungs[0].Fixture}, {rungs[0].Surface}, seed {rungs[0].Seed}";
            if (held.Length < 2)
            {
                _out.WriteLine($"{cell}, {what}: too few rungs ({held.Length})");
                Row("scaling", ScalingHeader, [ladder, rungs[0].Fixture, rungs[0].Surface, I(rungs[0].Cap), I(rungs[0].Seed), what, I(held.Length), "", "", "", "", "", "too few"]);
                return;
            }
            double[] n = [.. held.Select(r => (double)r.N)], y = [.. held.Select(r => value(r)!.Value)];
            (double b, double rssLinear, double c, double rssLog, double[] increments, string reading) = BassKernelSpectrum.Scaling(n, y);
            if (ladder == "L") reading = rssLinear < rssLog ? "linear fit lower" : "logarithmic fit lower";
            _out.WriteLine($"{cell}, {what}: {reading}; b {b:F4}, RSS {rssLinear:F3}; c {c:F3}, RSS {rssLog:F3}; Δ {string.Join(", ", increments.Select(d => d.ToString("F2", CultureInfo.InvariantCulture)))}");
            Row("scaling", ScalingHeader, [ladder, rungs[0].Fixture, rungs[0].Surface, I(rungs[0].Cap), I(rungs[0].Seed), what, I(held.Length), F(b), F(rssLinear), F(c), F(rssLog),
                                           string.Join(" ", increments.Select(F)), reading]);
        }
        One("-ln gap", r => r.Gap is { } g ? -Math.Log(g) : null);
        One("barrier", r => r.Barrier);
        One("ln relevant L", r => r.RelevantBound is { } l ? Math.Log(l) : null);
        for (int s = 0; s <= 2; s++)
        {
            int slots = s;
            One($"depth, {slots} free", r => r.DepthBySlack.TryGetValue(slots, out double d) ? d : null);
        }
    }

    private static readonly string[] ScalingHeader = ["ladder", "fixture", "surface", "cap", "seed", "quantity", "rungs", "b", "rss_linear", "c", "rss_log", "increments", "reading"];

    // ---- R̂ against the exact chain ----

    private static readonly int[] RHatSeeds = [7101, 7102, 7103, 7104];

    private static readonly string[] RHatHeader =
    [
        "fixture", "surface", "cap", "seed", "r", "N", "burn", "kept", "master", "predictive_rhat", "predictive_by", "predictive_ess", "level_rhat", "structural_rhat",
        "count_ess", "exact_count_ess", "gap", "T_gap", "relevant_L", "T_over_L", "tv_burn", "tv_end", "bn_mu_burn", "truth", "gate", "verdict_outcome",
    ];

    private void ReadRHat((string Fixture, string Surface, int Cap, int Seed, int R, int N, BassConfig Config, BassData Data) spec, bool farClocks, (int Burn, int Kept)[] lengths)
    {
        Rung rung = ReadRung("R", spec.Fixture, spec.Surface, spec.Cap, spec.Seed, spec.R, spec.N, spec.Config, spec.Data, 1, farClocks);
        (int Burn, int Kept, int Master)[] runs = [.. lengths.SelectMany(l => RHatSeeds.Select(m => (l.Burn, l.Kept, m)))];
        Parallel.ForEach(runs, new ParallelOptions { MaxDegreeOfParallelism = Threads }, run =>
        {
            BassConfig config = spec.Config with { Sampling = spec.Config.Sampling with { Chains = 4, BurnIn = run.Burn, Samples = run.Kept, MasterSeed = run.Master, ReportEvidence = false } };
            BassDiagnostics d = Bass.Fit(config, spec.Data).Diagnostics;
            (double predictive, string by, double ess) = BassDeterministicProfileTests.Predictive(d);
            double level = d.LogPosterior.RHat.State == DiagnosticState.Estimated ? d.LogPosterior.RHat.Value : d.LogPosterior.RHat.State == DiagnosticState.Constant ? 1.0 : double.PositiveInfinity;
            double structural = d.Count.RHat.State == DiagnosticState.Estimated ? d.Count.RHat.Value : d.Count.RHat.State == DiagnosticState.Constant ? 1.0 : double.PositiveInfinity;
            double tvBurn = rung.Tv[Array.IndexOf(Clocks, run.Burn)], tvEnd = rung.Tv[Array.IndexOf(Clocks, run.Burn + run.Kept)];
            double tOverL = rung.RelevantBound is { } l ? run.Kept / l : double.PositiveInfinity;
            string truth = tvBurn < 0.05 && tOverL >= 10 ? "mixed" : tvBurn > 0.25 || tOverL < 1 ? "not mixed" : "between";
            string gate = predictive <= 1.05 ? "pass" : "fail";
            Row("rhat", RHatHeader,
            [
                spec.Fixture, spec.Surface, I(spec.Cap), I(spec.Seed), I(spec.R), I(spec.N), I(run.Burn), I(run.Kept), I(run.Master), F(predictive), by, F(ess), F(level), F(structural),
                d.Count.BulkEss.State == DiagnosticState.Estimated ? F(d.Count.BulkEss.Value) : "", F(4.0 * run.Kept / rung.TauM), rung.Gap is { } g ? F(g) : "", rung.Gap is { } g2 ? F(run.Kept * g2) : "",
                rung.RelevantBound is { } rel ? F(rel) : "", F(tOverL), F(tvBurn), F(tvEnd), F(rung.BottleneckMass[Array.IndexOf(Clocks, run.Burn)]), truth, gate, d.Verdict.Outcome.ToString(),
            ]);
        });
    }

    [Theory]
    [MemberData(nameof(GridCells))]
    public void RHat_OnTheGrid(string surface, int cap, int seed)
    {
        if (!BassQualification.Enabled) return;
        (int, int)[] lengths = seed == 5 ? [(500, 200), (5_000, 2_000), (50_000, 20_000)] : [(500, 200), (5_000, 2_000)];
        foreach (int r in GridReps) ReadRHat(GridRung(surface, cap, seed, r), seed == 5, lengths);
    }

    [Theory]
    [InlineData("categorical, R's measure")]
    [InlineData("categorical, uniform")]
    [InlineData("functional")]
    public void RHat_OnTheCases(string fixture)
    {
        if (!BassQualification.Enabled) return;
        foreach (int r in CaseReps) ReadRHat(CaseRung(fixture, r), true, [(500, 200), (5_000, 2_000)]);
    }

    // ---- The flat histogram against the exact chain ----

    private static readonly int[] FlatSeeds = [7201, 7202];

    private static readonly string[] FlatHeader =
    [
        "fixture", "surface", "cap", "seed", "r", "N", "route", "master", "exact_ln_Z", "ln_Z", "error", "z", "evidence_agrees", "exact_top_basin", "exact_top_mass", "top_basin", "top_mass",
        "basin_tv", "basins_agree", "round_trips", "gap", "sweeps_gap", "seconds", "refusal",
    ];

    private void ReadFlat((string Fixture, string Surface, int Cap, int Seed, int R, int N, BassConfig Config, BassData Data) spec)
    {
        Rung rung = ReadRung("R", spec.Fixture, spec.Surface, spec.Cap, spec.Seed, spec.R, spec.N, spec.Config, spec.Data, 1, spec.Seed == 5);
        BassFlatHistogramSpace space = BassFlatHistogramSpace.Of(BassModel.Create(spec.Config, spec.Data));
        double exact = space.LogEvidence(1.0);
        const long learning = 30_000;
        const int production = 30_000, thin = 10, chains = 4;
        var pilot = new BassFlatHistogramSpec
        {
            Run = new FlatHistogramSpec { Chains = chains, LearningIterations = learning, ProductionDraws = production, Start = FlatHistogramStart.Wham },
            Pilot = new BassFlatHistogramPilot { Chains = 4, BurnIn = 5_000, TracedScans = 2_000, Tempering = new TemperingSpec { Form = TemperedForm.Linear, ReferenceRung = true } },
        };
        // The window declared from F's exact density of states, for the rungs whose pilot route refuses.
        var declared = new BassFlatHistogramSpec { Run = BassFlatHistogramTests.DeclaredRun(BassFlatHistogramTests.ExactShape(space, 12), learning, production, chains) };
        foreach ((string route, int master) in new[] { "pilot", "declared" }.SelectMany(route => FlatSeeds.Select(m => (route, m))))
        {
            BassConfig config = spec.Config with { Sampling = spec.Config.Sampling with { MasterSeed = master, Thin = thin, FlatHistogram = route == "pilot" ? pilot : declared } };
            var clock = System.Diagnostics.Stopwatch.StartNew();
            BassFit fit;
            try
            {
                fit = Bass.Fit(config, spec.Data);
            }
            catch (Exception refused) when (refused is InvalidOperationException or OverflowException)
            {
                Row("flat", FlatHeader,
                [
                    spec.Fixture, spec.Surface, I(spec.Cap), I(spec.Seed), I(spec.R), I(spec.N), route, I(master), F(exact), "", "", "", "refused", "", "", "", "",
                    "", "refused", "", rung.Gap is { } g0 ? F(g0) : "", rung.Gap is { } g1 ? F((learning + production) * thin * g1) : "", F(clock.Elapsed.TotalSeconds), $"{refused.GetType().Name}: {refused.Message}",
                ]);
                _out.WriteLine($"{spec.Fixture}, {spec.Surface}, r {spec.R}, {route}, master {master}: refused: {refused.Message}");
                continue;
            }
            clock.Stop();
            var named = new BassMixing.Fixture(spec.Fixture, spec.Data, spec.Config, null);
            var masses = new Dictionary<string, double>();
            double[] w = fit.DrawWeights!;
            int at = 0;
            double total = 0.0;
            foreach (IReadOnlyList<BassDraw> chain in fit.Draws)
                foreach (BassDraw draw in chain)
                {
                    string key = named.Basin(draw.Structure.TermIds.Select(id => fit.Terms[id]));
                    masses[key] = masses.GetValueOrDefault(key) + w[at];
                    total += w[at++];
                }
            double tv = 0.5 * rung.ExactBasins.Keys.Union(masses.Keys).Sum(k => Math.Abs(rung.ExactBasins.GetValueOrDefault(k) - masses.GetValueOrDefault(k) / total));
            KeyValuePair<string, double> exactTop = rung.ExactBasins.MaxBy(b => b.Value), top = masses.MaxBy(b => b.Value);
            BassEvidence? evidence = fit.Evidence;
            double error = evidence?.MonteCarloError ?? double.NaN, value = evidence?.LogEvidence ?? double.NaN;
            bool evidenceAgrees = Math.Abs(value - exact) <= Math.Max(3 * error, 0.1);
            bool basinsAgree = tv <= 0.05 && top.Key == exactTop.Key;
            long sweeps = (learning + production) * thin;
            Row("flat", FlatHeader,
            [
                spec.Fixture, spec.Surface, I(spec.Cap), I(spec.Seed), I(spec.R), I(spec.N), route, I(master), F(exact), F(value), F(error), F((value - exact) / error), evidenceAgrees ? "yes" : "no",
                exactTop.Key, F(exactTop.Value), top.Key, F(top.Value / total), F(tv), basinsAgree ? "yes" : "no", I((int)(fit.Diagnostics.FlatHistogram!.Production?.TotalRoundTrips ?? 0)),
                rung.Gap is { } g ? F(g) : "", rung.Gap is { } g2 ? F(sweeps * g2) : "", F(clock.Elapsed.TotalSeconds), "",
            ]);
            _out.WriteLine($"{spec.Fixture}, {spec.Surface}, r {spec.R}, {route}, master {master}: ln Ẑ {value:F3} ± {error:F3} against {exact:F3}; basin TV {tv:F3}; top {top.Key} against {exactTop.Key}; {clock.Elapsed.TotalSeconds:F0} s");
        }
    }

    [Theory]
    [InlineData(InsideF)]
    [InlineData(OutsideF)]
    public void FlatHistogram_OnTheGrid(string surface)
    {
        if (!BassQualification.Enabled) return;
        foreach (int r in new[] { 2, 8, 32, 128 }) ReadFlat(GridRung(surface, 3, 5, r));
    }

    [Theory]
    [InlineData("categorical, R's measure")]
    [InlineData("functional")]
    public void FlatHistogram_OnTheCases(string fixture)
    {
        if (!BassQualification.Enabled) return;
        ReadFlat(CaseRung(fixture, 1));
    }

    // ---- The rows ----

    private static int Threads => int.TryParse(Environment.GetEnvironmentVariable("BASS_GAP_THREADS"), NumberStyles.Integer, CultureInfo.InvariantCulture, out int t) && t > 0 ? t : 6;

    private static string F(double v) => double.IsNaN(v) ? "" : v.ToString("R", CultureInfo.InvariantCulture);

    private static string I(int v) => v.ToString(CultureInfo.InvariantCulture);

    private static readonly object RowLock = new();

    private static void Row(string file, string[] header, IEnumerable<string?> cells)
    {
        string folder = Environment.GetEnvironmentVariable("BASS_GAP_OUT") is { Length: > 0 } d ? d : Path.Combine(Path.GetTempPath(), "bass-gap");
        static string Quote(string? v) => v is null ? "" : v.Contains(',') || v.Contains('"') ? $"\"{v.Replace("\"", "\"\"")}\"" : v;
        lock (RowLock)
        {
            Directory.CreateDirectory(folder);
            string path = Path.Combine(folder, $"spectral-gap-{file}.csv");
            if (!File.Exists(path)) File.WriteAllText(path, string.Join(",", header) + "\n");
            string[] row = [.. cells.Select(Quote)];
            Assert.Equal(header.Length, row.Length);
            File.AppendAllText(path, string.Join(",", row) + "\n");
        }
    }
}
