using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime;
using System.Text;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Repo.TestHarness;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-R1, the deterministic-simulator profile's pre-registered comparison (the deterministic-profile brief, decided
/// 2026-10-30): eight arms, the coefficient prior's form and scale crossed with two count priors, on four noise-free
/// benchmarks of three designs each, read by the holdout error of the posterior mean surface and the latent band's
/// coverage, and judged by the rule in F-M2's form.
/// <list type="bullet">
/// <item><b>Arms</b> (<see cref="Arms"/>): S10's configuration (four chains, 5,000 burn-in, 2,000 draws, untempered) with
/// the brief's coefficient prior and count prior. A0 is today's default.</item>
/// <item><b>Benchmarks</b> (<see cref="NoiseFree"/>): Friedman's f₁ (p = 10, five inert), Ishigami (a = 7, b = 0.1, on
/// [−π, π]³), the borehole function (p = 8, Morris, Mitchell and Ylvisaker 1993's ranges), each at N = 500 from a random
/// Latin hypercube; and S10's direct functional mean on 60 uniform cases × 20 locations. The noisy controls
/// (<see cref="NoisyControls"/>) are Friedman with N(0, 1) noise and Ishigami with N(0, 0.37²), on the same designs.</item>
/// <item><b>Holdouts:</b> 10,000 uniform points per scalar surface, 200 uniform cases at the 20 training locations for the
/// functional one, one per surface from <see cref="HoldoutSeed"/>, shared by every arm, design and noisy twin.</item>
/// <item><b>Seeds,</b> declared before any run: a (surface, design k) has data seed 31,000 + 100·surface + k, noise seed
/// 32,000 + 100·surface + k (the noisy twins only) and master seed 33,000 + 100·surface + k, so a noisy twin shares its
/// noise-free twin's design and master seed. The arms of a (benchmark, design) run in an order drawn from
/// <see cref="OrderSeed"/>.</item>
/// <item><b>Readings</b> (<see cref="Reading"/>): NRMSE, the holdout RMSE of the posterior mean (the latent draws' mean)
/// over the truth's population sd on the holdout; the share of holdout points whose truth lies in the latent draws'
/// pointwise 95% band (type-7 quantiles), never the predictive one; the predictive group's rank-normalized R̂ (the largest
/// over σ² and the default monitored fitted values) and its least ESS; D10's verdict; M's median and 95% interval and
/// whether any draw reached the cap; σ's median; the whole <see cref="Bass.Fit"/>'s seconds; and for the functional
/// benchmark the lag-one residual-dependence check on the training residuals.</item>
/// <item><b>The gate</b> (<see cref="Rule"/>): a (benchmark, design) where any arm's predictive R̂ exceeds 1.05 is left out
/// for every arm. A reading with no estimate (frozen, unavailable, or a group all constant) cannot certify agreement and
/// counts as above the gate. A benchmark with two or more designs left out, or none kept, is unscored.</item>
/// <item><b>The rule</b> (<see cref="Rule"/>): r_f = exp(mean over kept designs of ln NRMSE(arm) − the same for A0); the
/// score the geometric mean of r_f over the scored benchmarks; an arm clears at a score of at most 0.9 with no r_f above
/// 1.05 and its design-averaged coverage no more than 0.05 below A0's on any scored benchmark; the lowest clearing score
/// wins. Each threshold is met by a value equal to it up to 10⁻¹² (<see cref="Slack"/>), so rounding cannot decide.</item>
/// </list>
/// Every fit runs alone on one thread (MaxDegreeOfParallelism 1) after an untimed warm-up. BASS_FR1_SCALE divides the
/// lengths and BASS_FR1_DESIGNS keeps the first designs, for the 1/10 trial. The per-fit CSV and the report land under
/// <c>artifacts/test-runs/BassDeterministicProfile/</c>. The comparison and the checks that fit run as qualification runs
/// (<see cref="BassQualification.Enabled"/>); the benchmarks', the arms', the metrics' and the rule's checks run always.
/// </summary>
public sealed class BassDeterministicProfileTests
{
    /// <summary>The designs per benchmark, the scalar design size, the holdout sizes and the functional fixture's shape.</summary>
    internal const int Designs = 3, DesignSize = 500, HoldoutSize = 10_000, FunctionalCases = 60, FunctionalHoldoutCases = 200, Locations = 20;

    /// <summary>The seeds' bases (see the class summary), the holdouts' base, and the seed the arms' shuffle derives from.</summary>
    internal const int DataSeed = 31_000, NoiseSeed = 32_000, MasterSeed = 33_000, HoldoutSeed = 34_000, OrderSeed = 20261030;

    /// <summary>The gate on the predictive R̂; the rule's score, its ceiling on any r_f, and its coverage margin.</summary>
    internal const double Gate = 1.05, Clears = 0.9, Ceiling = 1.05, CoverageMargin = 0.05;

    /// <summary>The rounding a threshold allows a value equal to it.</summary>
    internal const double Slack = 1e-12;

    /// <summary>Ishigami's a and b.</summary>
    internal const double IshigamiA = 7.0, IshigamiB = 0.1;

    // Holdout points per prediction call, so the latent draws stay near 64 MB.
    private const int ChunkPoints = 1000;

#if DEBUG
    private const string Configuration = "Debug";
#else
    private const string Configuration = "Release";
#endif

    private readonly ITestOutputHelper _out;
    public BassDeterministicProfileTests(ITestOutputHelper output) => _out = output;

    /// <summary>An arm: the coefficient prior and the count prior on S10's sampling.</summary>
    internal sealed record Arm(string Name, string Label, BassCoefficientPrior Coefficients, BassCountPrior Count)
    {
        public bool Ridge => Coefficients.Form is RidgeForm;

        public bool Held => Coefficients.Scale.IsHeld;
    }

    /// <summary>The brief's table, A0 the default.</summary>
    internal static readonly Arm[] Arms =
    [
        new("A0", "g, Zellner–Siow, λ (10, 10)", BassCoefficientPrior.G(), new BassCountPrior(10, 10)),
        new("A1", "g, τ held at 1/N, λ (10, 10)", BassCoefficientPrior.G(TauPrior.Held()), new BassCountPrior(10, 10)),
        new("A2", "ridge, τ sampled, λ (10, 10)", BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()), new BassCountPrior(10, 10)),
        new("A3", "ridge, τ held at 1, λ (10, 10)", BassCoefficientPrior.Ridge(TauPrior.Held(1.0)), new BassCountPrior(10, 10)),
        new("A4", "g, Zellner–Siow, λ (1, 1)", BassCoefficientPrior.G(), new BassCountPrior(1, 1)),
        new("A5", "g, τ held at 1/N, λ (1, 1)", BassCoefficientPrior.G(TauPrior.Held()), new BassCountPrior(1, 1)),
        new("A6", "ridge, τ sampled, λ (1, 1)", BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()), new BassCountPrior(1, 1)),
        new("A7", "ridge, τ held at 1, λ (1, 1)", BassCoefficientPrior.Ridge(TauPrior.Held(1.0)), new BassCountPrior(1, 1)),
    ];

    /// <summary>The surfaces, by index: the seeds read it.</summary>
    internal const int FriedmanSurface = 0, IshigamiSurface = 1, BoreholeSurface = 2, FunctionalSurface = 3;

    /// <summary>A benchmark: its surface and its noise's sd, 0 for noise-free.</summary>
    internal sealed record Benchmark(string Name, string Label, int Surface, double NoiseSd)
    {
        public bool Functional => Surface == FunctionalSurface;
    }

    /// <summary>The four noise-free benchmarks the rule reads.</summary>
    internal static readonly Benchmark[] NoiseFree =
    [
        new("friedman", "Friedman", FriedmanSurface, 0.0),
        new("ishigami", "Ishigami", IshigamiSurface, 0.0),
        new("borehole", "borehole", BoreholeSurface, 0.0),
        new("functional", "functional", FunctionalSurface, 0.0),
    ];

    /// <summary>The noisy controls, read and not ruled.</summary>
    internal static readonly Benchmark[] NoisyControls =
    [
        new("friedman-noisy", "Friedman, σ = 1", FriedmanSurface, 1.0),
        new("ishigami-noisy", "Ishigami, σ = 0.37", IshigamiSurface, 0.37),
    ];

    // The scalar surfaces' inputs: names and ranges, each uniform.
    private static readonly (string Name, double Low, double High)[] FriedmanInputs =
        [.. Enumerable.Range(1, 10).Select(v => ($"x{v}", 0.0, 1.0))];

    private static readonly (string Name, double Low, double High)[] IshigamiInputs =
        [("x1", -Math.PI, Math.PI), ("x2", -Math.PI, Math.PI), ("x3", -Math.PI, Math.PI)];

    // r_w, r, T_u, H_u, T_l, H_l, L, K_w on Morris, Mitchell and Ylvisaker's ranges.
    private static readonly (string Name, double Low, double High)[] BoreholeInputs =
        [("rw", 0.05, 0.15), ("r", 100, 50_000), ("Tu", 63_070, 115_600), ("Hu", 990, 1110), ("Tl", 63.1, 116), ("Hl", 700, 820),
         ("L", 1120, 1680), ("Kw", 9855, 12_045)];

    internal static (string Name, double Low, double High)[] Inputs(int surface) => surface switch
    {
        FriedmanSurface => FriedmanInputs,
        IshigamiSurface => IshigamiInputs,
        BoreholeSurface => BoreholeInputs,
        _ => throw new ArgumentOutOfRangeException(nameof(surface)),
    };

    /// <summary>Friedman's f₁, S10's expression: 10 sin(πx₁x₂) + 20(x₃ − ½)² + 10x₄ + 5x₅.</summary>
    internal static double Friedman(ReadOnlySpan<double> x)
        => 10 * Math.Sin(Math.PI * x[0] * x[1]) + 20 * Math.Pow(x[2] - 0.5, 2) + 10 * x[3] + 5 * x[4];

    /// <summary>Ishigami: sin x₁ + a sin² x₂ + b x₃⁴ sin x₁.</summary>
    internal static double Ishigami(ReadOnlySpan<double> x)
        => Math.Sin(x[0]) + IshigamiA * Math.Pow(Math.Sin(x[1]), 2) + IshigamiB * Math.Pow(x[2], 4) * Math.Sin(x[0]);

    /// <summary>The borehole's flow, 2πT_u(H_u − H_l) / (ln(r/r_w)·(1 + 2LT_u/(ln(r/r_w)·r_w²·K_w) + T_u/T_l)).</summary>
    internal static double Borehole(ReadOnlySpan<double> x)
    {
        double rw = x[0], r = x[1], tu = x[2], hu = x[3], tl = x[4], hl = x[5], l = x[6], kw = x[7];
        double lnr = Math.Log(r / rw);
        return 2 * Math.PI * tu * (hu - hl) / (lnr * (1 + 2 * l * tu / (lnr * rw * rw * kw) + tu / tl));
    }

    /// <summary>S10's direct functional mean, 1 + 2x₀ + 1.5·t·x₁ + 3(x₀ − ½)₊·sin(πt).</summary>
    internal static double Functional(double x0, double x1, double t)
        => 1 + 2 * x0 + 1.5 * t * x1 + 3 * Math.Max(x0 - 0.5, 0) * Math.Sin(Math.PI * t);

    internal static double Evaluate(int surface, ReadOnlySpan<double> x) => surface switch
    {
        FriedmanSurface => Friedman(x),
        IshigamiSurface => Ishigami(x),
        BoreholeSurface => Borehole(x),
        _ => throw new ArgumentOutOfRangeException(nameof(surface)),
    };

    /// <summary>The training locations, S10's t on 20 points of [0, 1].</summary>
    internal static double[] TrainingLocations() => [.. Enumerable.Range(0, Locations).Select(r => r / (Locations - 1.0))];

    /// <summary>A design's data and its truth at the training points, case-major with locations inner.</summary>
    internal sealed record Problem(BassData Data, double[] Truth);

    /// <summary>A holdout: its input columns by name, its truth case-major with locations inner, and its locations.</summary>
    internal sealed record Holdout(string[] Names, double[][] Columns, double[] Truth, int Locations)
    {
        public int Cases => Columns[0].Length;

        /// <summary>Cases [start, start + count) as prediction cases.</summary>
        public PredictionCases Slice(int start, int count)
            => new([.. Names.Select((name, v) => (PredictionColumn)new ContinuousColumn(name, Columns[v][start..(start + count)]))]);

        public string Hash()
        {
            var h = new RunFingerprint().Add(Locations);
            for (int v = 0; v < Names.Length; v++) h.Add(Names[v]).Add(Columns[v]);
            return h.Add(Truth).ToString();
        }
    }

    /// <summary>
    /// A random Latin hypercube of <paramref name="n"/> points over <paramref name="ranges"/>: for each input in turn, a
    /// Fisher–Yates permutation of the n cells, then a uniform position within each point's cell.
    /// </summary>
    internal static double[][] LatinHypercube(int n, (string Name, double Low, double High)[] ranges, Xoshiro256PlusPlus rng)
    {
        var x = new double[ranges.Length][];
        var cells = new int[n];
        for (int v = 0; v < ranges.Length; v++)
        {
            for (int i = 0; i < n; i++) cells[i] = i;
            for (int i = n - 1; i > 0; i--)
            {
                int j = (int)(rng.NextDouble() * (i + 1));
                (cells[i], cells[j]) = (cells[j], cells[i]);
            }
            (_, double low, double high) = ranges[v];
            x[v] = new double[n];
            for (int i = 0; i < n; i++) x[v][i] = low + (high - low) * ((cells[i] + rng.NextDouble()) / n);
        }
        return x;
    }

    /// <summary>Benchmark <paramref name="b"/>'s design <paramref name="k"/>: the scalar surfaces from a Latin hypercube of
    /// 500, the functional one from 60 uniform cases at the training locations; a noisy twin adds its noise.</summary>
    internal static Problem Design(Benchmark b, int k)
    {
        var rng = new Xoshiro256PlusPlus(DataSeed + 100 * b.Surface + k);
        var noise = new Xoshiro256PlusPlus(NoiseSeed + 100 * b.Surface + k);
        if (b.Functional)
        {
            double[] t = TrainingLocations();
            var x0 = new double[FunctionalCases];
            var x1 = new double[FunctionalCases];
            var y = new double[FunctionalCases, Locations];
            var truth = new double[FunctionalCases * Locations];
            for (int i = 0; i < FunctionalCases; i++)
            {
                x0[i] = rng.NextDouble();
                x1[i] = rng.NextDouble();
                for (int r = 0; r < Locations; r++)
                {
                    truth[i * Locations + r] = Functional(x0[i], x1[i], t[r]);
                    y[i, r] = truth[i * Locations + r] + (b.NoiseSd > 0 ? b.NoiseSd * BassQualification.Gaussian(noise) : 0.0);
                }
            }
            return new Problem(new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)],
                                            new FunctionalResponse(y, [new ContinuousPredictor("t", t)])), truth);
        }
        (string Name, double Low, double High)[] inputs = Inputs(b.Surface);
        double[][] x = LatinHypercube(DesignSize, inputs, rng);
        var f = new double[DesignSize];
        var yScalar = new double[DesignSize];
        var point = new double[inputs.Length];
        for (int i = 0; i < DesignSize; i++)
        {
            for (int v = 0; v < inputs.Length; v++) point[v] = x[v][i];
            f[i] = Evaluate(b.Surface, point);
            yScalar[i] = f[i] + (b.NoiseSd > 0 ? b.NoiseSd * BassQualification.Gaussian(noise) : 0.0);
        }
        return new Problem(new BassData([.. inputs.Select((input, v) => (BassPredictor)new ContinuousPredictor(input.Name, x[v]))], new ScalarResponse(yScalar)), f);
    }

    /// <summary>Surface <paramref name="surface"/>'s holdout: 10,000 uniform points, or 200 uniform functional cases at the
    /// training locations, drawn point by point from its seed.</summary>
    internal static Holdout HoldoutOf(int surface)
    {
        var rng = new Xoshiro256PlusPlus(HoldoutSeed + surface);
        if (surface == FunctionalSurface)
        {
            double[] t = TrainingLocations();
            var x0 = new double[FunctionalHoldoutCases];
            var x1 = new double[FunctionalHoldoutCases];
            var truth = new double[FunctionalHoldoutCases * Locations];
            for (int i = 0; i < FunctionalHoldoutCases; i++)
            {
                x0[i] = rng.NextDouble();
                x1[i] = rng.NextDouble();
                for (int r = 0; r < Locations; r++) truth[i * Locations + r] = Functional(x0[i], x1[i], t[r]);
            }
            return new Holdout(["x0", "x1"], [x0, x1], truth, Locations);
        }
        (string Name, double Low, double High)[] inputs = Inputs(surface);
        double[][] columns = [.. inputs.Select(_ => new double[HoldoutSize])];
        var f = new double[HoldoutSize];
        var point = new double[inputs.Length];
        for (int i = 0; i < HoldoutSize; i++)
        {
            for (int v = 0; v < inputs.Length; v++)
            {
                (_, double low, double high) = inputs[v];
                point[v] = columns[v][i] = low + (high - low) * rng.NextDouble();
            }
            f[i] = Evaluate(surface, point);
        }
        return new Holdout([.. inputs.Select(input => input.Name)], columns, f, 1);
    }

    /// <summary>Arm <paramref name="arm"/> on benchmark <paramref name="b"/>'s design <paramref name="k"/>: S10's
    /// configuration at its master seed with the arm's priors, its lengths divided by <paramref name="scale"/>.</summary>
    internal static BassConfig Config(Benchmark b, int k, Arm arm, int scale = 1)
    {
        BassConfig s10 = BassQualificationRecoveryTests.Config(MasterSeed + 100 * b.Surface + k);
        return s10 with
        {
            Coefficients = arm.Coefficients,
            Count = arm.Count,
            Sampling = s10.Sampling with { BurnIn = s10.Sampling.BurnIn / scale, Samples = s10.Sampling.Samples / scale },
        };
    }

    // ---------------------------------------------------------------------------------------------------------------
    // The checks before the runs.

    /// <summary>
    /// The benchmarks are what they say: Ishigami's mean a/2 and variance a²/8 + bπ⁴/5 + b²π⁸/18 + ½ by Monte Carlo over 10⁶
    /// points within 5 SE; Friedman's f₁ and the functional mean against S10's generators at their own inputs, to rounding;
    /// the borehole increasing in r_w and H_u and decreasing in H_l and L at 1,000 random points. The designs are Latin
    /// (one point per cell of every input), within their ranges and reproducible, and the holdouts are fixed.
    /// </summary>
    [Fact]
    public void Benchmarks_AreWhatTheySay()
    {
        // Ishigami's moments.
        var rng = new Xoshiro256PlusPlus(seed: 20261030);
        const int draws = 1_000_000;
        var values = new double[draws];
        var x = new double[3];
        for (int i = 0; i < draws; i++)
        {
            for (int v = 0; v < 3; v++) x[v] = -Math.PI + 2 * Math.PI * rng.NextDouble();
            values[i] = Ishigami(x);
        }
        double mean = values.Average();
        double m2 = values.Average(v => (v - mean) * (v - mean)), m4 = values.Average(v => Math.Pow(v - mean, 4));
        double expectedMean = IshigamiA / 2;
        double expectedVariance = IshigamiA * IshigamiA / 8 + IshigamiB * Math.Pow(Math.PI, 4) / 5 + IshigamiB * IshigamiB * Math.Pow(Math.PI, 8) / 18 + 0.5;
        double seMean = Math.Sqrt(m2 / draws), seVariance = Math.Sqrt((m4 - m2 * m2) / draws);
        Assert.InRange(Math.Abs(mean - expectedMean) / seMean, 0.0, 5.0);
        Assert.InRange(Math.Abs(m2 - expectedVariance) / seVariance, 0.0, 5.0);
        Assert.InRange(expectedVariance, 13.84, 13.85);
        _out.WriteLine($"Ishigami: mean {mean:F5} against {expectedMean} ({Math.Abs(mean - expectedMean) / seMean:F2} SE), variance {m2:F4} " +
                       $"against {expectedVariance:F4} ({Math.Abs(m2 - expectedVariance) / seVariance:F2} SE), 10⁶ points");

        // Friedman's f₁ at S10's inputs.
        (BassData s10, double[] s10Truth) = BassQualificationRecoveryTests.Friedman();
        double worst = 0.0;
        var p10 = new double[10];
        for (int i = 0; i < s10Truth.Length; i++)
        {
            for (int v = 0; v < 10; v++) p10[v] = ((ContinuousPredictor)s10.Predictors[v]).Values.Span[i];
            worst = Math.Max(worst, Math.Abs(Friedman(p10) - s10Truth[i]));
        }
        Assert.True(worst <= 1e-12, $"Friedman's f₁ differs from S10's by {worst:G3}");
        _out.WriteLine($"Friedman: within {worst:G3} of S10's truth at its 500 inputs");

        // The functional mean at S10's inputs.
        (BassData functional, double[] functionalTruth) = BassQualificationRecoveryTests.Functional();
        ReadOnlySpan<double> fx0 = ((ContinuousPredictor)functional.Predictors[0]).Values.Span, fx1 = ((ContinuousPredictor)functional.Predictors[1]).Values.Span;
        double[] t = TrainingLocations();
        Assert.Equal(t, ((FunctionalResponse)functional.Response).OutputVariables[0].Values.ToArray());
        worst = 0.0;
        for (int i = 0; i < FunctionalCases; i++)
            for (int r = 0; r < Locations; r++) worst = Math.Max(worst, Math.Abs(Functional(fx0[i], fx1[i], t[r]) - functionalTruth[i * Locations + r]));
        Assert.True(worst <= 1e-12, $"the functional mean differs from S10's by {worst:G3}");
        _out.WriteLine($"functional: within {worst:G3} of S10's truth at its 60 × 20 points");

        // The borehole's monotone directions: a step of 1% of the range up, from points where it fits.
        var bore = new Xoshiro256PlusPlus(seed: 1993);
        var point = new double[8];
        (int Input, int Sign)[] directions = [(0, 1), (3, 1), (5, -1), (6, -1)];
        for (int s = 0; s < 1000; s++)
        {
            for (int v = 0; v < 8; v++) point[v] = BoreholeInputs[v].Low + 0.99 * (BoreholeInputs[v].High - BoreholeInputs[v].Low) * bore.NextDouble();
            double f = Borehole(point);
            foreach ((int input, int sign) in directions)
            {
                double[] moved = [.. point];
                moved[input] += 0.01 * (BoreholeInputs[input].High - BoreholeInputs[input].Low);
                Assert.True(sign * (Borehole(moved) - f) > 0, $"borehole: {BoreholeInputs[input].Name} moves f the wrong way at sample {s}");
            }
        }
        _out.WriteLine("borehole: increasing in r_w and H_u, decreasing in H_l and L, at 1,000 random points");

        // The designs and the holdouts.
        foreach (Benchmark b in NoiseFree)
            for (int k = 0; k < Designs; k++)
            {
                Problem design = Design(b, k), again = Design(b, k);
                Assert.Equal(Bass.Fingerprint(new BassConfig(), design.Data), Bass.Fingerprint(new BassConfig(), again.Data));
                Assert.Equal(design.Truth, again.Truth);
                if (b.Functional)
                {
                    Assert.Equal(FunctionalCases, design.Data.Response.Cases);
                    continue;
                }
                (string Name, double Low, double High)[] inputs = Inputs(b.Surface);
                for (int v = 0; v < inputs.Length; v++)
                {
                    ReadOnlySpan<double> column = ((ContinuousPredictor)design.Data.Predictors[v]).Values.Span;
                    var cells = new int[DesignSize];
                    foreach (double value in column)
                    {
                        Assert.InRange(value, inputs[v].Low, inputs[v].High);
                        cells[Math.Min(DesignSize - 1, (int)((value - inputs[v].Low) / (inputs[v].High - inputs[v].Low) * DesignSize))]++;
                    }
                    Assert.All(cells, c => Assert.Equal(1, c));
                }
                Assert.Equal(design.Truth, ((ScalarResponse)design.Data.Response).Y.ToArray());   // noise-free: y is the truth
            }
        foreach (Benchmark b in NoisyControls)
        {
            Problem noisy = Design(b, 0), clean = Design(NoiseFree.Single(n => n.Surface == b.Surface), 0);
            Assert.Equal(clean.Truth, noisy.Truth);
            double[] residual = [.. ((ScalarResponse)noisy.Data.Response).Y.ToArray().Zip(noisy.Truth, (y, f) => y - f)];
            double sd = Math.Sqrt(residual.Average(e => e * e));
            Assert.InRange(sd / b.NoiseSd, 0.85, 1.15);
            _out.WriteLine($"{b.Name}: design 0's noise sd {sd:F3} against {b.NoiseSd}; the designs' inputs are the noise-free twin's");
        }
        for (int surface = 0; surface <= FunctionalSurface; surface++)
        {
            Holdout h = HoldoutOf(surface);
            Assert.Equal(h.Hash(), HoldoutOf(surface).Hash());
            Assert.Equal(surface == FunctionalSurface ? FunctionalHoldoutCases * Locations : HoldoutSize, h.Truth.Length);
            _out.WriteLine($"holdout {surface}: {h.Cases} cases × {h.Locations}, truth sd {PopulationSd(h.Truth):G4}, hash {h.Hash()}");
        }
        _out.WriteLine($"Ishigami's noise 0.37 against its sd {Math.Sqrt(expectedVariance):F3}: {0.37 / Math.Sqrt(expectedVariance):F3}");
    }

    /// <summary>
    /// The arms are what the table says: each resolves its form, τ's prior on its τ₀ (1/N under the g-form, 1 under the
    /// ridge), the held τ where it holds one and the count prior as its row, on S10's sampling and seed; its fingerprint
    /// differs from every other arm's, and with its coefficient and count priors set back to A0's it is A0's own.
    /// </summary>
    [Fact]
    public void Arms_AreWhatTheySay()
    {
        foreach (Benchmark b in NoiseFree.Concat(NoisyControls))
        {
            Problem p = Design(b, 0);
            string a0 = Bass.Fingerprint(Config(b, 0, Arms[0]), p.Data);
            var prints = new HashSet<string>();
            foreach (Arm arm in Arms)
            {
                BassConfig config = Config(b, 0, arm);
                BassSampling s = config.Sampling;
                Assert.Equal((4, MasterSeed + 100 * b.Surface, 5000, 2000, 1), (s.Chains, s.MasterSeed, s.BurnIn, s.Samples, s.Thin));
                Assert.Null(s.Tempering);
                Assert.Null(s.Smc);
                Assert.Null(s.FlatHistogram);
                Assert.Equal(new BassConfig().Sigma2, config.Sigma2);
                Assert.Equal(GPriorCentre.Mean, config.Intercept);

                BassResolved r = BassResolved.Resolve(config, p.Data);
                double n0 = arm.Ridge ? 1.0 : r.Observations;
                Assert.Equal(arm.Ridge ? typeof(RidgeForm) : typeof(GForm), config.Coefficients.Form.GetType());
                Assert.Equal(arm.Count, config.Count);
                Assert.Equal(n0, r.Scale.UnitObservations);
                Assert.Equal(1.0 / n0, r.Scale.UnitInformation);
                if (arm.Held)
                {
                    Assert.Equal((double?)(arm.Ridge ? 1.0 : 1.0 / r.Observations), r.Scale.HeldTau);
                    Assert.Null(r.Scale.Gamma);
                }
                else
                {
                    Assert.Null(r.Scale.HeldTau);
                    Assert.Equal(0.5, r.Scale.Gamma!.Value.Shape);
                    Assert.Equal(0.5 * n0, r.Scale.Gamma.Value.Rate, 1e-12 * n0);
                }
                Assert.True(prints.Add(Bass.Fingerprint(config, p.Data)), $"{b.Name}, {arm.Name}: the fingerprint repeats another arm's.");
                Assert.Equal(a0, Bass.Fingerprint(config with { Coefficients = Arms[0].Coefficients, Count = Arms[0].Count }, p.Data));
            }
            _out.WriteLine($"{b.Name}: N = {BassResolved.Resolve(Config(b, 0, Arms[0]), p.Data).Observations}, eight distinct fingerprints, " +
                           $"A0's {a0} under each arm's priors set back");
        }
    }

    /// <summary>
    /// The metrics on a posterior set by hand: four points with truth 0, 1, 2, 3 and four draws each, placed so the mean is
    /// the truth + 1, −1, +1, −1 and the 95% band (type 7: x₀ + 0.075(x₁ − x₀), x₂ + 0.925(x₃ − x₂)) covers points 0 and 3
    /// only. NRMSE is 1/√1.25, the truth's population sd being √1.25, and the coverage ½.
    /// </summary>
    [Fact]
    public void Metrics_OnAPosteriorSetByHand()
    {
        double[] truth = [0, 1, 2, 3];
        double[][] draws =
        [
            [-1, 0.5, 1.5, 3],         // mean 1, band [−0.8875, 2.8875] ∋ 0
            [-0.5, -0.25, 0.25, 0.5],  // mean 0, band [−0.48125, 0.48125] ∌ 1
            [2.5, 2.75, 3.25, 3.5],    // mean 3, band [2.51875, 3.48125] ∌ 2
            [0, 1, 2, 5],              // mean 2, band [0.075, 4.775] ∋ 3
        ];
        var latent = new double[4, 4, 1];
        for (int d = 0; d < 4; d++)
            for (int i = 0; i < 4; i++) latent[d, i, 0] = draws[i][d];
        var mean = new double[4];
        var lower = new double[4];
        var upper = new double[4];
        Band(latent, 0, mean, lower, upper);
        Assert.Equal(new[] { 1.0, 0.0, 3.0, 2.0 }, mean);
        Assert.Equal(new[] { -0.8875, -0.48125, 2.51875, 0.075 }, lower, new ToleranceComparer(1e-15));
        Assert.Equal(new[] { 2.8875, 0.48125, 3.48125, 4.775 }, upper, new ToleranceComparer(1e-15));
        (double nrmse, double coverage) = Metrics(mean, lower, upper, truth);
        Assert.Equal(1 / Math.Sqrt(1.25), nrmse, 1e-15);
        Assert.Equal(0.5, coverage);
        _out.WriteLine($"hand-set posterior: NRMSE {nrmse:R} against 1/√1.25, coverage {coverage} against ½");
    }

    /// <summary>
    /// The holdout path reads the latent surface: a short fit's chunked prediction at its own training cases gives the
    /// posterior mean and the 95% band of the draws' fitted values (S10's <see cref="BassQualification.Fitted"/>), to
    /// rounding, for a scalar and a functional benchmark under A0 and A7.
    /// </summary>
    [Fact]
    public void HoldoutPath_ReadsTheLatentSurface()
    {
        if (!BassQualification.Enabled) return;   // fits
        foreach (Benchmark b in new[] { NoiseFree[0], NoiseFree[3] })
            foreach (Arm arm in new[] { Arms[0], Arms[7] })
            {
                Problem p = Design(b, 0);
                BassConfig config = Config(b, 0, arm, scale: 50);
                BassFit fit = Bass.Fit(config, p.Data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
                Holdout training = new([.. p.Data.Predictors.Select(q => q.Name)],
                                       [.. p.Data.Predictors.Select(q => ((ContinuousPredictor)q).Values.ToArray())], p.Truth, b.Functional ? Locations : 1);
                (double[] mean, double[] lower, double[] upper) = Surface(fit, p.Data, training);
                double[][] fitted = BassQualification.Fitted(BassModel.Create(config, p.Data), fit);
                double worst = 0.0;
                for (int k = 0; k < p.Truth.Length; k++)
                {
                    (double m, _, double lo, double hi) = BassQualification.At(fitted, k);
                    double scaleK = 1.0 + Math.Abs(m);
                    worst = Math.Max(worst, Math.Max(Math.Abs(mean[k] - m), Math.Max(Math.Abs(lower[k] - lo), Math.Abs(upper[k] - hi))) / scaleK);
                }
                Assert.True(worst <= 1e-9, $"{b.Name}, {arm.Name}: the holdout path is {worst:G3} from the draws' fitted values");
                _out.WriteLine($"{b.Name}, {arm.Name}: the holdout path at the training cases within {worst:G3} (relative) of the draws' fitted values");
            }
    }

    /// <summary>
    /// The pairing holds before the runs, on every benchmark's first design at short lengths: the eight arms read the same
    /// data and master seed, the four arms of each count prior start every chain from the same structure, the two sampled
    /// arms of a count prior from the same τ/τ₀, and the held arms from their held τ. Across count priors the starts are
    /// each prior's own exact draw from the same streams (S7), so they differ.
    /// </summary>
    [Fact]
    public void Pairing_Holds()
    {
        if (!BassQualification.Enabled) return;   // fits
        foreach (Benchmark b in NoiseFree.Concat(NoisyControls))
        {
            Problem p = Design(b, 0);
            var starts = new List<(Arm Arm, BassFit Fit)>();
            foreach (Arm arm in Arms)
            {
                BassConfig config = Config(b, 0, arm);
                config = config with { Sampling = config.Sampling with { BurnIn = 5, Samples = 10 } };
                starts.Add((arm, Bass.Fit(config, p.Data, new BassRunOptions { MaxDegreeOfParallelism = 1 })));
            }
            CheckPairing(b.Name, [.. starts.Select(s => (s.Arm, Account(s.Fit, p.Data)))]);
            _out.WriteLine($"{b.Name}: one data fingerprint and master seed; starts paired within each count prior");
        }
    }

    // A fit's pairing account: the data it read (fingerprinted under the default config, so every arm reads one value), its
    // master seed, and each chain's start, its structure and τ/τ₀.
    private sealed record PairingAccount(string Data, int MasterSeed, string[] Structures, double[] TauUnits);

    private static PairingAccount Account(BassFit fit, BassData data)
    {
        double tau0 = fit.Resolved.Scale.UnitInformation;
        return new PairingAccount(
            Bass.Fingerprint(new BassConfig(), data), fit.Identity.MasterSeed,
            [.. fit.Diagnostics.Chains.Select(c => string.Join(";", c.Start.Terms.Select(term => term.ToString())))],
            [.. fit.Diagnostics.Chains.Select(c => c.Start.Tau / tau0)]);
    }

    private static void CheckPairing(string where, IReadOnlyList<(Arm Arm, PairingAccount Account)> arms)
    {
        Assert.Single(arms.Select(a => a.Account.Data).Distinct());
        Assert.Single(arms.Select(a => a.Account.MasterSeed).Distinct());
        foreach (IGrouping<BassCountPrior, (Arm Arm, PairingAccount Account)> prior in arms.GroupBy(a => a.Arm.Count))
        {
            Assert.Equal(4, prior.Count());
            Assert.Single(prior.Select(a => string.Join("|", a.Account.Structures)).Distinct());
            (Arm Arm, PairingAccount Account)[] sampled = [.. prior.Where(a => !a.Arm.Held)];
            Assert.Equal(2, sampled.Length);
            for (int c = 0; c < sampled[0].Account.TauUnits.Length; c++)
            {
                double a = sampled[0].Account.TauUnits[c], z = sampled[1].Account.TauUnits[c];
                Assert.True(Math.Abs(a - z) <= 1e-12 * Math.Abs(a), $"{where}, chain {c}: τ/τ₀ starts {a:R} and {z:R} under {prior.Key}");
            }
            foreach ((Arm arm, PairingAccount account) in prior.Where(a => a.Arm.Held))
                Assert.All(account.TauUnits, u => Assert.Equal(1.0, u, 1e-12));
        }
    }

    // ---------------------------------------------------------------------------------------------------------------
    // The comparison.

    /// <summary>The comparison: every arm on every benchmark and design, the noisy controls, the gate, the rule and the
    /// report.</summary>
    [Fact]
    public void Comparison()
    {
        if (!BassQualification.Enabled) return;   // heavy
        int scale = Setting("BASS_FR1_SCALE", 1), designs = Math.Min(Designs, Setting("BASS_FR1_DESIGNS", Designs));
        Benchmark[] benchmarks = [.. NoiseFree, .. NoisyControls];
        ArtifactRun run = HarnessArtifacts.Create("test-runs", "BassDeterministicProfile", "Maths.Regression.Tests", new Dictionary<string, object?>
        {
            ["Scale"] = scale, ["Designs"] = designs, ["OrderSeed"] = OrderSeed, ["Configuration"] = Configuration,
            ["ServerGC"] = GCSettings.IsServerGC, ["Processors"] = Environment.ProcessorCount,
        });
        string csvPath = Path.Combine(run.RunDirectory, "deterministic-profile-results.csv");
        using var csv = new StreamWriter(csvPath) { NewLine = "\n" };
        csv.WriteLine(Reading.Header);
        var clock = Stopwatch.StartNew();

        var holdouts = new Dictionary<int, Holdout>();
        for (int surface = 0; surface <= FunctionalSurface; surface++) holdouts[surface] = HoldoutOf(surface);

        // The warm-up: every arm on every noise-free benchmark's first design, short and discarded.
        foreach (Benchmark b in NoiseFree)
        {
            Problem p = Design(b, 0);
            foreach (Arm arm in Arms) Bass.Fit(Config(b, 0, arm, Math.Max(50, scale)), p.Data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        }
        _out.WriteLine($"warm-up {clock.Elapsed.TotalSeconds:F0} s ({Configuration})");

        var readings = new List<Reading>();
        int[] orders = SeedTree.Derive(OrderSeed, Designs * benchmarks.Length);
        for (int k = 0; k < designs; k++)
            for (int bi = 0; bi < benchmarks.Length; bi++)
            {
                Benchmark b = benchmarks[bi];
                int block = k * benchmarks.Length + bi;
                Problem p = Design(b, k);
                Arm[] order = Shuffled(Arms, orders[block]);
                var held = new List<(Arm Arm, Reading Reading, PairingAccount Account)>();
                for (int position = 0; position < order.Length; position++)
                {
                    (Reading r, PairingAccount account) = Read(b, k, order[position], scale, block, position, p, holdouts[b.Surface]);
                    csv.WriteLine(r.Row());
                    csv.Flush();
                    held.Add((order[position], r, account));
                }
                CheckPairing($"{b.Name}, k = {k}", [.. held.Select(h => (h.Arm, h.Account))]);
                readings.AddRange(held.Select(h => h.Reading));
                _out.WriteLine($"{b.Name}, k = {k}: {string.Join(", ", held.Select(h => $"{h.Arm.Name} {h.Reading.Seconds:F1}+{h.Reading.ReadSeconds:F1} s"))}");
            }

        Verdict verdict = Rule(readings, [.. NoiseFree.Select(b => b.Name)], designs);
        Verdict controls = Rule(readings, [.. NoisyControls.Select(b => b.Name)], designs);
        string report = Report(verdict, controls, readings, scale, designs, clock.Elapsed);
        _out.WriteLine(report);
        _out.WriteLine($"CSV {csvPath}");
        _out.WriteLine($"report {run.WriteText("deterministic-profile-report.md", report)}");

        // Checked after the report is written, so a failure keeps it: the holdout is one per surface, whatever the arm,
        // design or noise, and every reading is produced.
        foreach (IGrouping<int, Reading> surface in readings.GroupBy(r => benchmarks.Single(b => b.Name == r.Benchmark).Surface))
            Assert.Single(surface.Select(r => r.HoldoutHash).Distinct());
        Assert.All(readings, r =>
        {
            Assert.True(double.IsFinite(r.Nrmse) && double.IsFinite(r.Coverage) && double.IsFinite(r.MedianM) && double.IsFinite(r.SigmaMedian),
                        $"{r.Benchmark}, k = {r.K}, {r.Arm}: a reading is missing");
            Assert.True(!r.Functional || double.IsFinite(r.ResidualMean), $"{r.Benchmark}, k = {r.K}, {r.Arm}: no residual check");
        });
    }

    private static int Setting(string name, int fallback)
        => int.TryParse(Environment.GetEnvironmentVariable(name), NumberStyles.Integer, CultureInfo.InvariantCulture, out int v) && v > 0 ? v : fallback;

    // The arms in an order drawn from `seed` (Fisher–Yates), as F-M2's pass shuffles them.
    private static Arm[] Shuffled(Arm[] arms, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        Arm[] order = [.. arms];
        for (int i = order.Length - 1; i > 0; i--)
        {
            int j = (int)(rng.NextDouble() * (i + 1));
            (order[i], order[j]) = (order[j], order[i]);
        }
        return order;
    }

    /// <summary>One fit's readings, a row of the CSV.</summary>
    internal sealed record Reading(
        string Benchmark, bool Functional, int K, int MasterSeed, string Arm, int Block, int Position,
        double Seconds, double CpuSeconds, double ReadSeconds,
        double Nrmse, double Coverage,
        double PredictiveRHat, string PredictiveBy, double PredictiveEss, string Verdict,
        double MedianM, double LowerM, double UpperM, int MaxM, int MaxBasis, double MassAtCap,
        double SigmaMedian, double TauUnitsMedian,
        double ResidualMean, double ResidualExpected, double ResidualMedian,
        string HoldoutHash, string DrawsHash, string Fingerprint)
    {
        /// <summary>Whether any retained draw sits at the basis cap.</summary>
        public bool AtCap => MaxM >= MaxBasis;

        public static string Header => string.Join(",",
            "benchmark", "k", "master_seed", "arm", "block", "position", "seconds", "cpu_seconds", "read_seconds", "nrmse", "coverage",
            "predictive_rhat", "predictive_rhat_by", "predictive_ess", "verdict", "median_M", "lower_M", "upper_M", "max_M", "max_basis",
            "mass_at_cap", "sigma_median", "tau_over_tau0_median", "residual_r1_mean", "residual_r1_expected", "residual_r1_median",
            "holdout_hash", "draws_hash", "fingerprint");

        public string Row() => string.Join(",",
            Benchmark, I(K), I(MasterSeed), Arm, I(Block), I(Position), F(Seconds), F(CpuSeconds), F(ReadSeconds), F(Nrmse), F(Coverage),
            F(PredictiveRHat), PredictiveBy, F(PredictiveEss), Verdict, F(MedianM), F(LowerM), F(UpperM), I(MaxM), I(MaxBasis),
            F(MassAtCap), F(SigmaMedian), F(TauUnitsMedian), F(ResidualMean), F(ResidualExpected), F(ResidualMedian),
            HoldoutHash, DrawsHash, Fingerprint);
    }

    private static string F(double value) => double.IsNaN(value) ? "" : value.ToString("R", CultureInfo.InvariantCulture);
    private static string I(int value) => value.ToString(CultureInfo.InvariantCulture);

    // One fit timed on one thread, then read: everything after Bass.Fit is untimed.
    private static (Reading Reading, PairingAccount Account) Read(Benchmark b, int k, Arm arm, int scale, int block, int position, Problem p, Holdout holdout)
    {
        BassConfig config = Config(b, k, arm, scale);
        GC.Collect();
        GC.WaitForPendingFinalizers();
        GC.Collect();
        TimeSpan cpu = Process.GetCurrentProcess().TotalProcessorTime;
        var clock = Stopwatch.StartNew();
        BassFit fit = Bass.Fit(config, p.Data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        double seconds = clock.Elapsed.TotalSeconds;
        double cpuSeconds = (Process.GetCurrentProcess().TotalProcessorTime - cpu).TotalSeconds;

        var reading = Stopwatch.StartNew();
        Assert.Null(fit.DrawWeights);
        Assert.Equal(arm.Count, fit.Config.Count);
        Assert.Equal(arm.Coefficients, fit.Config.Coefficients);
        (double[] mean, double[] lower, double[] upper) = Surface(fit, p.Data, holdout);
        (double nrmse, double coverage) = Metrics(mean, lower, upper, holdout.Truth);

        BassDiagnostics d = fit.Diagnostics;
        (double rHat, string by, double ess) = Predictive(d);
        string verdict = $"{Short(d.Verdict.Outcome)} (P {Short(d.Verdict.Group(BassDiagnostics.PredictiveGroup).Outcome)}, " +
                         $"L {Short(d.Verdict.Group(BassDiagnostics.LevelGroup).Outcome)}, S {Short(d.Verdict.Group(BassDiagnostics.StructuralGroup).Outcome)})";
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        double[] m = [.. draws.Select(draw => (double)draw.Structure.TermIds.Length).Order()];
        double[] sigma = [.. draws.Select(draw => Math.Sqrt(draw.Sigma2)).Order()];
        double[] tau = [.. draws.Select(draw => draw.Tau / fit.Resolved.Scale.UnitInformation).Order()];
        ResidualDependence? check = fit.ResidualDependence;
        if (b.Functional) Assert.Equal(ResidualStatistic.LagOne, check!.Statistic);
        double residualMedian = check is null ? double.NaN : check.Quantiles[Array.IndexOf(check.Probabilities, 0.5)];
        double readSeconds = reading.Elapsed.TotalSeconds;

        return (new Reading(b.Name, b.Functional, k, config.Sampling.MasterSeed, arm.Name, block, position, seconds, cpuSeconds, readSeconds,
                            nrmse, coverage, rHat, by, ess, verdict,
                            PosteriorInterval.Quantile(m, 0.5), PosteriorInterval.Quantile(m, 0.025), PosteriorInterval.Quantile(m, 0.975),
                            (int)m[^1], fit.Resolved.MaxBasis, d.MassAtCap,
                            PosteriorInterval.Quantile(sigma, 0.5), PosteriorInterval.Quantile(tau, 0.5),
                            check?.Mean ?? double.NaN, check?.Expected ?? double.NaN, residualMedian,
                            holdout.Hash(), BassBitIdentityTests.Hash(fit), fit.Identity.Fingerprint),
                Account(fit, p.Data));
    }

    private static string Short(ConvergenceOutcome outcome) => outcome switch
    {
        ConvergenceOutcome.Converged => "C",
        ConvergenceOutcome.NotConverged => "N",
        _ => "U",
    };

    /// <summary>
    /// The gate's reading: the predictive group's largest rank-normalized R̂ over σ² and the monitored fitted values, the
    /// quantity that set it, and the group's least ESS (the less favourable of bulk and tail). A constant quantity is left
    /// out; a frozen or unavailable one, or a group with no estimate, reads +∞.
    /// </summary>
    internal static (double RHat, string By, double Ess) Predictive(BassDiagnostics d)
    {
        (string Name, ConvergenceDiagnostic Diagnostic)[] group = [("sigma2", d.Sigma2), .. d.Fitted.Select((v, i) => ($"f{i + 1}", v.Diagnostic))];
        double worst = double.NegativeInfinity, least = double.PositiveInfinity;
        string by = "";
        foreach ((string name, ConvergenceDiagnostic q) in group)
        {
            if (q.RHat.State == DiagnosticState.Constant) continue;
            double value = q.RHat.State == DiagnosticState.Estimated ? q.RHat.Value : double.PositiveInfinity;
            if (value > worst) (worst, by) = (value, name);
            least = Math.Min(least, q.Ess.State == DiagnosticState.Estimated ? q.Ess.Value : 0.0);
        }
        return double.IsNegativeInfinity(worst) ? (double.PositiveInfinity, "none estimated", 0.0) : (worst, by, least);
    }

    /// <summary>The posterior mean and the 95% latent band at every holdout point (case-major, locations inner), the latent
    /// draws predicted a chunk of cases at a time.</summary>
    internal static (double[] Mean, double[] Lower, double[] Upper) Surface(BassFit fit, BassData data, Holdout holdout)
    {
        int points = holdout.Truth.Length, perChunk = Math.Max(1, ChunkPoints / holdout.Locations);
        var mean = new double[points];
        var lower = new double[points];
        var upper = new double[points];
        for (int start = 0; start < holdout.Cases; start += perChunk)
        {
            int count = Math.Min(perChunk, holdout.Cases - start);
            BassPrediction prediction = fit.Predict(data, holdout.Slice(start, count), new BassPredictOptions { Predictive = false });
            Assert.Equal(holdout.Locations, prediction.Locations);
            Band(prediction.Latent, start * holdout.Locations, mean, lower, upper);
        }
        return (mean, lower, upper);
    }

    /// <summary>Each point's mean over the draws of <paramref name="latent"/> (draw × case × location), and its 2.5% and
    /// 97.5% type-7 quantiles, written from <paramref name="offset"/>.</summary>
    internal static void Band(double[,,] latent, int offset, double[] mean, double[] lower, double[] upper)
    {
        int s = latent.GetLength(0), n = latent.GetLength(1), m = latent.GetLength(2);
        var values = new double[s];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++)
            {
                double sum = 0.0;
                for (int d = 0; d < s; d++) sum += values[d] = latent[d, i, l];
                Array.Sort(values);
                int k = offset + i * m + l;
                mean[k] = sum / s;
                lower[k] = PosteriorInterval.Quantile(values, 0.025);
                upper[k] = PosteriorInterval.Quantile(values, 0.975);
            }
    }

    /// <summary>NRMSE, the RMSE of <paramref name="mean"/> against <paramref name="truth"/> over the truth's population sd,
    /// and the share of points whose truth lies in [lower, upper].</summary>
    internal static (double Nrmse, double Coverage) Metrics(double[] mean, double[] lower, double[] upper, double[] truth)
    {
        double squares = 0.0;
        int covered = 0;
        for (int k = 0; k < truth.Length; k++)
        {
            squares += (mean[k] - truth[k]) * (mean[k] - truth[k]);
            if (lower[k] <= truth[k] && truth[k] <= upper[k]) covered++;
        }
        return (Math.Sqrt(squares / truth.Length) / PopulationSd(truth), covered / (double)truth.Length);
    }

    private static double PopulationSd(double[] values)
    {
        double mean = values.Average();
        return Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / values.Length);
    }

    // ---------------------------------------------------------------------------------------------------------------
    // The rule.

    /// <summary>A fit left out by the gate: its benchmark, design and arm, and its predictive R̂ and the quantity that set it.</summary>
    internal sealed record Exclusion(string Benchmark, int K, string Arm, double RHat, string By);

    /// <summary>
    /// One arm's standing: per scored benchmark ln r_f, its design-averaged coverage and A0's, the per-design log ratios
    /// and their sd; the score, whether it clears, and the score with each scored benchmark left out.
    /// </summary>
    internal sealed record Standing(string Arm, double[] LogRatio, double[] Coverage, double[] A0Coverage, double[][] DesignLogRatios,
                                    double[] Spread, double Score, bool Clears, double[] LeaveOneOut);

    /// <summary>The gate's exclusions, the scored benchmarks with their kept designs, the unscored ones, the standings and
    /// the winner; null when none clears.</summary>
    internal sealed record Verdict(Exclusion[] Exclusions, string[] Scored, int[][] Kept, string[] Unscored, Standing[] Standings, string? Winner);

    /// <summary>
    /// The gate, then the rule, as decided. The gate leaves out a (benchmark, design) for every arm when any arm's predictive
    /// R̂ exceeds 1.05 (or has no estimate); a benchmark with two or more designs left out, or none kept, is unscored. Per
    /// scored benchmark f, r_f = exp(mean over kept designs of ln NRMSE(arm) − the same for A0); the score is the geometric
    /// mean of r_f; an arm clears when its score is at most 0.9, no r_f exceeds 1.05, and its coverage averaged over kept
    /// designs is no more than 0.05 below A0's on any scored benchmark. The lowest clearing score wins.
    /// </summary>
    internal static Verdict Rule(IReadOnlyList<Reading> readings, string[] benchmarks, int designs)
    {
        Reading At(string arm, string benchmark, int k) => readings.Single(r => r.Arm == arm && r.Benchmark == benchmark && r.K == k);
        Exclusion[] exclusions =
        [
            .. readings.Where(r => benchmarks.Contains(r.Benchmark) && r.K < designs && !(r.PredictiveRHat <= Gate))
                       .OrderBy(r => Array.IndexOf(benchmarks, r.Benchmark)).ThenBy(r => r.K).ThenBy(r => r.Arm)
                       .Select(r => new Exclusion(r.Benchmark, r.K, r.Arm, r.PredictiveRHat, r.PredictiveBy)),
        ];
        var scored = new List<string>();
        var kept = new List<int[]>();
        var unscored = new List<string>();
        foreach (string f in benchmarks)
        {
            int[] keep = [.. Enumerable.Range(0, designs).Where(k => !exclusions.Any(e => e.Benchmark == f && e.K == k))];
            if (designs - keep.Length >= 2 || keep.Length == 0) unscored.Add(f);
            else
            {
                scored.Add(f);
                kept.Add(keep);
            }
        }

        string a0 = Arms[0].Name;
        var standings = new List<Standing>();
        foreach (Arm arm in Arms.Skip(1))
        {
            int n = scored.Count;
            double[] logRatio = new double[n], coverage = new double[n], a0Coverage = new double[n], spread = new double[n];
            var designRatios = new double[n][];
            for (int i = 0; i < n; i++)
            {
                string f = scored[i];
                int[] keep = kept[i];
                logRatio[i] = keep.Average(k => Math.Log(At(arm.Name, f, k).Nrmse)) - keep.Average(k => Math.Log(At(a0, f, k).Nrmse));
                coverage[i] = keep.Average(k => At(arm.Name, f, k).Coverage);
                a0Coverage[i] = keep.Average(k => At(a0, f, k).Coverage);
                designRatios[i] = [.. keep.Select(k => Math.Log(At(arm.Name, f, k).Nrmse / At(a0, f, k).Nrmse))];
                spread[i] = StandardDeviation(designRatios[i]);
            }
            double score = n == 0 ? double.NaN : Math.Exp(logRatio.Average());
            bool clears = n > 0 && score <= Clears + Slack && logRatio.All(l => Math.Exp(l) <= Ceiling + Slack)
                          && Enumerable.Range(0, n).All(i => a0Coverage[i] - coverage[i] <= CoverageMargin + Slack);
            double[] leaveOneOut = [.. Enumerable.Range(0, n).Select(i => n < 2 ? double.NaN : Math.Exp(logRatio.Where((_, j) => j != i).Average()))];
            standings.Add(new Standing(arm.Name, logRatio, coverage, a0Coverage, designRatios, spread, score, clears, leaveOneOut));
        }
        return new Verdict(exclusions, [.. scored], [.. kept], [.. unscored], [.. standings], standings.Where(s => s.Clears).MinBy(s => s.Score)?.Arm);
    }

    private static double StandardDeviation(double[] values)
    {
        if (values.Length < 2) return double.NaN;
        double mean = values.Average();
        return Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (values.Length - 1));
    }

    /// <summary>
    /// The rule's arithmetic against a hand computation, on synthetic readings of two benchmarks and three designs:
    /// <list type="bullet">
    /// <item>benchmark "a": A0's NRMSE 0.1, 0.2, 0.4; A1's half of it; A2's 0.8 of it at design 0 and twice at design 1;
    /// A3's half with coverage 0.10 below A0's. Benchmark "b": A0's 1, A1's 0.81, A2's 0.9, A3's 0.5.</item>
    /// <item>The gate: A4's predictive R̂ 1.06 at ("b", 2) leaves that design out for every arm, so "b" is scored on two;
    /// A5's frozen R̂ (+∞) at ("c", 0) and ("c", 1) leaves "c" unscored.</item>
    /// <item>A1: r_a = ½, r_b = 0.81, score √(½·0.81) = 0.6364, clears. A2: r_a = exp((ln 0.8 + ln 2 + 0)/3) = 1.1696
    /// above 1.05, fails. A3: r_a = ½, r_b = ½, score ½, fails on coverage. A1 wins.</item>
    /// </list>
    /// </summary>
    [Fact]
    public void Rule_MatchesAHandComputation()
    {
        var readings = new List<Reading>();
        void Add(string f, int k, string arm, double nrmse, double coverage = 0.95, double rHat = 1.001)
            => readings.Add(new Reading(f, false, k, 0, arm, 0, 0, 1, 1, 0, nrmse, coverage, rHat, "f1", 400, "C", 5, 3, 8, 9, 499, 0,
                                        0.01, 1, double.NaN, double.NaN, double.NaN, "h", "d", "p"));
        double[] a0 = [0.1, 0.2, 0.4];
        for (int k = 0; k < 3; k++)
        {
            Add("a", k, "A0", a0[k]);
            Add("a", k, "A1", a0[k] / 2);
            Add("a", k, "A2", a0[k] * (k == 0 ? 0.8 : k == 1 ? 2.0 : 1.0));
            Add("a", k, "A3", a0[k] / 2, coverage: 0.85);
            Add("b", k, "A0", 1.0);
            Add("b", k, "A1", 0.81);
            Add("b", k, "A2", 0.9);
            Add("b", k, "A3", 0.5);
            Add("c", k, "A0", 1.0);
            Add("c", k, "A1", 0.5);
            Add("c", k, "A2", 0.5);
            Add("c", k, "A3", 0.5);
            foreach (string f in new[] { "a", "b", "c" })
                foreach (string arm in new[] { "A4", "A5", "A6", "A7" })
                    Add(f, k, arm, 1.0, rHat: (f, k, arm) switch
                    {
                        ("b", 2, "A4") => 1.06,
                        ("c", 0 or 1, "A5") => double.PositiveInfinity,
                        _ => 1.001,
                    });
        }

        Verdict v = Rule(readings, ["a", "b", "c"], 3);
        Assert.Equal(new[] { "a", "b" }, v.Scored);
        Assert.Equal(new[] { "c" }, v.Unscored);
        Assert.Equal(new[] { 0, 1, 2 }, v.Kept[0]);
        Assert.Equal(new[] { 0, 1 }, v.Kept[1]);
        Assert.Equal(3, v.Exclusions.Length);
        Standing a1 = v.Standings.Single(s => s.Arm == "A1"), a2 = v.Standings.Single(s => s.Arm == "A2"), a3 = v.Standings.Single(s => s.Arm == "A3");
        Assert.Equal(Math.Log(0.5), a1.LogRatio[0], 1e-14);
        Assert.Equal(Math.Log(0.81), a1.LogRatio[1], 1e-14);
        Assert.Equal(Math.Sqrt(0.5 * 0.81), a1.Score, 1e-14);
        Assert.True(a1.Clears);
        Assert.Equal(Math.Exp((Math.Log(0.8) + Math.Log(2.0)) / 3), Math.Exp(a2.LogRatio[0]), 1e-14);
        Assert.False(a2.Clears);
        Assert.Equal(0.5, a3.Score, 1e-14);
        Assert.Equal(0.85, a3.Coverage[0], 1e-14);
        Assert.False(a3.Clears);
        Assert.Equal("A1", v.Winner);
        Assert.Equal(new[] { 0.81, 0.5 }, a1.LeaveOneOut, new ToleranceComparer(1e-14));
        _out.WriteLine($"hand computation: A1 {a1.Score:F4} clears, A2 r_a {Math.Exp(a2.LogRatio[0]):F4} fails, A3 {a3.Score:F4} fails on coverage; winner {v.Winner}");
    }

    private sealed class ToleranceComparer(double tolerance) : IEqualityComparer<double>
    {
        public bool Equals(double a, double b) => Math.Abs(a - b) <= tolerance;
        public int GetHashCode(double value) => 0;
    }

    // ---------------------------------------------------------------------------------------------------------------
    // The report.

    // The verdict first, then the gate, the context, the noisy controls, the per-fit table and the times, as markdown.
    private static string Report(Verdict verdict, Verdict controls, List<Reading> readings, int scale, int designs, TimeSpan elapsed)
    {
        var md = new StringBuilder();
        string Label(string arm) => $"{arm}, {Arms.Single(a => a.Name == arm).Label}";
        string BenchLabel(string name) => NoiseFree.Concat(NoisyControls).Single(b => b.Name == name).Label;
        void Row(params object[] cells) => md.Append("| ").AppendJoin(" | ", cells).Append(" |\n");
        void Head(params string[] cells)
        {
            Row(cells);
            Row([.. cells.Select(_ => "---")]);
        }

        md.Append($"Lengths 1/{scale} of S10's, {designs} design(s) per benchmark, {Configuration}, {Environment.ProcessorCount} logical processors, " +
                  $"server GC {GCSettings.IsServerGC}; {elapsed.TotalMinutes:F1} minutes in all.\n\n");

        void Standings(Verdict v, bool ruled)
        {
            Head(["arm", "score", .. v.Scored.Select(f => $"r_f {BenchLabel(f)}"), .. v.Scored.Select(f => $"coverage {BenchLabel(f)} (A0's)"), ruled ? "clears" : "loses to A0 (score > 1)"]);
            foreach (Standing s in v.Standings)
                Row([Label(s.Arm), $"{s.Score:F3}", .. s.LogRatio.Select(l => $"{Math.Exp(l):F3}"),
                     .. s.Coverage.Select((c, i) => $"{c:F3} ({s.A0Coverage[i]:F3})"), ruled ? (s.Clears ? "yes" : "no") : (s.Score > 1 ? "yes" : "no")]);
        }

        md.Append("### The verdict\n\n");
        Standings(verdict, ruled: true);
        md.Append(verdict.Winner is { } w
            ? $"\nThe winner is {Label(w)}, score {verdict.Standings.Single(s => s.Arm == w).Score:F3}.\n\n"
            : "\nNo arm clears: no profile is built, every option stays an option, and F-R1 closes as read.\n\n");

        void Gate(Verdict v, string what)
        {
            md.Append($"### The gate ({what})\n\n");
            md.Append(v.Exclusions.Length == 0 ? "No fit's predictive R̂ exceeds 1.05; every design is kept.\n" : "Fits whose predictive R̂ exceeds 1.05, each leaving its (benchmark, design) out for every arm:\n\n");
            foreach (Exclusion e in v.Exclusions) md.Append($"- {BenchLabel(e.Benchmark)}, k = {e.K}, {e.Arm}: R̂ {e.RHat:F3} ({e.By})\n");
            md.Append($"\nScored: {string.Join(", ", v.Scored.Select((f, i) => $"{BenchLabel(f)} (designs {string.Join(", ", v.Kept[i])})"))}" +
                      $"{(v.Scored.Length == 0 ? "none" : "")}. Unscored: {(v.Unscored.Length == 0 ? "none" : string.Join(", ", v.Unscored.Select(BenchLabel)))}.\n\n");
        }
        Gate(verdict, "noise-free");

        md.Append("### Context for the verdict\n\nEach arm's per-design ln(NRMSE/NRMSE_A0), and their sd:\n\n");
        Head(["arm", .. verdict.Scored.Select(f => BenchLabel(f))]);
        foreach (Standing s in verdict.Standings)
            Row([Label(s.Arm), .. s.DesignLogRatios.Select((r, i) => $"{string.Join(", ", r.Select(x => x.ToString("F3")))} (sd {s.Spread[i]:F3})")]);
        md.Append("\nEach arm's score with one benchmark left out:\n\n");
        Head(["arm", .. verdict.Scored.Select(f => $"without {BenchLabel(f)}")]);
        foreach (Standing s in verdict.Standings) Row([Label(s.Arm), .. s.LeaveOneOut.Select(x => $"{x:F3}")]);

        md.Append("\n### The noisy controls (read, not ruled)\n\n");
        Standings(controls, ruled: false);
        md.Append(verdict.Winner is { } winner && controls.Standings.SingleOrDefault(s => s.Arm == winner) is { } ws
            ? $"\nThe winner {winner} on the noisy controls: score {ws.Score:F3}, r_f {string.Join(", ", ws.LogRatio.Select((l, i) => $"{BenchLabel(controls.Scored[i])} {Math.Exp(l):F3}"))}.\n\n"
            : "\n");
        Gate(controls, "noisy controls");

        md.Append("### Per fit\n\nNRMSE and coverage on the holdout; the predictive R̂ (the quantity setting it) and least ESS; D10's verdict (run; predictive, level, structural); " +
                  "M's median [95%] and largest, against the cap; σ's median; τ/τ₀'s median; the fit's seconds; the functional residual r₁'s mean against its expectation:\n\n");
        Head("benchmark", "k", "arm", "NRMSE", "coverage", "pred. R̂", "pred. ESS", "D10", "M", "max M / cap", "σ med.", "τ/τ₀ med.", "s", "r₁");
        foreach (Reading r in readings.OrderBy(r => r.Block).ThenBy(r => r.Arm))
            Row(BenchLabel(r.Benchmark), r.K, r.Arm, $"{r.Nrmse:G4}", $"{r.Coverage:F3}", $"{r.PredictiveRHat:F3} ({r.PredictiveBy})", $"{r.PredictiveEss:F0}", r.Verdict,
                $"{r.MedianM:G4} [{r.LowerM:G4}, {r.UpperM:G4}]", $"{r.MaxM} / {r.MaxBasis}", $"{r.SigmaMedian:G3}", $"{r.TauUnitsMedian:G3}", $"{r.Seconds:F1}",
                r.Functional ? $"{r.ResidualMean:F3} ({r.ResidualExpected:F3})" : "");

        md.Append("\n### Times\n\nThe whole fit's seconds, the mean over designs, and the reading's (prediction on the holdout) beside it:\n\n");
        Head(["arm", .. NoiseFree.Concat(NoisyControls).Select(b => b.Label)]);
        foreach (Arm arm in Arms)
            Row([Label(arm.Name), .. NoiseFree.Concat(NoisyControls).Select(b =>
            {
                Reading[] c = [.. readings.Where(r => r.Arm == arm.Name && r.Benchmark == b.Name)];
                return c.Length == 0 ? "" : $"{c.Average(r => r.Seconds):F1} ({c.Average(r => r.ReadSeconds):F1})";
            })]);
        double fits = readings.Sum(r => r.Seconds), reads = readings.Sum(r => r.ReadSeconds);
        md.Append($"\nThe fits: {fits / 60:F1} minutes timed; the readings {reads / 60:F1} minutes; CPU over wall {readings.Sum(r => r.CpuSeconds) / fits:F3}.\n");
        if (scale > 1 || designs < Designs)
            md.Append($"\nProjected at full lengths and {Designs} designs: about {(fits + reads) * scale * Designs / designs / 3600:F2} hours of fits and readings " +
                      "(lengths × the scale, designs × the design count), the warm-up beside.\n");

        Reading[] capped = [.. readings.Where(r => r.AtCap)];
        md.Append($"\nFits with a draw at the basis cap: {(capped.Length == 0 ? "none" : string.Join("; ", capped.Select(r => $"{r.Benchmark} k = {r.K} {r.Arm} ({r.MassAtCap:P1} of draws)")))}.\n");
        return md.ToString();
    }
}
