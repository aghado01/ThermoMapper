using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E36's slice 2 recovery gate (the BPPR brief): functions from outside the ridge model, whose directions are known in closed form
/// (<see cref="ProjectionPursuitReferences"/>), fitted with a directional basis and read through E8 on its sampled route. The principal
/// angles between E[C]'s leading k eigenvectors and the true span are the gate, under the plain chain, under tempering and under SMC. Every
/// input declares its reference box as its domain, so E8's box metric is the inputs' own units, the directions' frame.
/// </summary>
public sealed class DirectionalRecoveryTests(ITestOutputHelper output)
{
    private static BassData Data(ProjectionPursuitReference reference, int n, double noise, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        int p = reference.Function.Dimension;
        ReferenceBox box = reference.Law;
        var x = new double[p][];
        for (int j = 0; j < p; j++) x[j] = new double[n];
        var y = new double[n];
        var point = new double[p];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < p; j++) point[j] = x[j][i] = box.Lower[j] + (box.Upper[j] - box.Lower[j]) * rng.NextDouble();
            y[i] = reference.Function.Value(point) + noise * Normal.SampleStandard(rng);
        }
        return new BassData([.. Enumerable.Range(0, p).Select(j => new ContinuousPredictor($"x{j}", x[j], new InputDomain(box.Lower[j], box.Upper[j])))],
                            new ScalarResponse(y));
    }

    private static BassConfig Config(string sampler, int seed)
    {
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 3, MaxBasis = 25, Directional = new DirectionalBasis() },
            // 2,000 burn-in sweeps leave the two ridges unsettled (a second angle of 0.30 rad, M 7.1 against 4.0 after 10,000).
            Sampling = new BassSampling { Chains = 2, MasterSeed = seed, BurnIn = 10_000, Samples = 1000, Thin = 5 },
        };
        return sampler switch
        {
            "tempering" => config with { Sampling = config.Sampling with { Tempering = new TemperingSpec { InverseTemperatures = [1.0, 0.5, 0.25, 0.1] } } },
            // Four mutation steps a stage leave the weaker ridge unsettled (0.21 rad on the two ridges); twelve give 0.018, thirty 0.006.
            "smc" => config with { Sampling = config.Sampling with { Smc = new SmcSpec { Particles = 200, Replicates = 2, MutationSteps = 12, TargetConditionalEss = 0.9 } } },
            _ => config,
        };
    }

    // The true directions orthonormalized by Gram–Schmidt, as columns.
    private static double[,] Truth(ProjectionPursuitReference reference)
    {
        int p = reference.Function.Dimension, k = reference.Directions.Length;
        var basis = new double[p, k];
        for (int c = 0; c < k; c++)
        {
            double[] v = (double[])reference.Directions[c].Clone();
            for (int b = 0; b < c; b++)
            {
                double dot = 0.0;
                for (int j = 0; j < p; j++) dot += v[j] * basis[j, b];
                for (int j = 0; j < p; j++) v[j] -= dot * basis[j, b];
            }
            double norm = Math.Sqrt(v.Sum(e => e * e));
            for (int j = 0; j < p; j++) basis[j, c] = v[j] / norm;
        }
        return basis;
    }

    [Theory]
    [InlineData("single index", "chain", 3621)]
    [InlineData("two ridges", "chain", 3622)]
    [InlineData("single index", "tempering", 3623)]
    [InlineData("two ridges", "tempering", 3625)]
    [InlineData("single index", "smc", 3624)]
    [InlineData("two ridges", "smc", 3626)]
    public void TheReferenceDirections_AreRecoveredThroughE8(string name, string sampler, int seed)
    {
        ProjectionPursuitReference reference = name == "two ridges" ? ProjectionPursuitReferences.TwoRidges() : ProjectionPursuitReferences.SingleIndex();
        BassData data = Data(reference, 300, 0.1, seed);
        BassFit fit = Bass.Fit(Config(sampler, seed), data);
        int k = reference.Directions.Length, p = reference.Function.Dimension;
        ActiveSubspaceAnalysis analysis = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Dimension = k, Points = 1024 });
        Assert.Null(analysis.Poincare);   // the sampled route, which a directional fit takes
        double[,] vectors = analysis.Posterior.MeanCVectors;
        var leading = new double[p, k];
        for (int j = 0; j < p; j++)
            for (int c = 0; c < k; c++) leading[j, c] = vectors[j, c];
        double[] angles = SubspacePosterior.PrincipalAngles(leading, Truth(reference));
        double[] eigenvalues = analysis.Posterior.MeanCEigenvalues;
        int[] counts = [.. fit.Draws.SelectMany(c => c).Select(d => d.Structure.TermIds.Length)];
        output.WriteLine($"{name}, {sampler}: principal angles to the truth {string.Join(", ", angles.Select(a => $"{a:F4}"))} rad; "
                         + $"E[C]'s eigenvalues {string.Join(", ", eigenvalues.Select(e => $"{e:G3}"))}; mean M {counts.Average():F1}");
        Assert.True(angles.Max() <= 0.1, $"the largest principal angle is {angles.Max():F4} rad");
    }

    // Several outputs (E34) share one structure, whose columns a directional basis makes as it makes a scalar fit's: nothing on that path
    // reads a term's parts. Each output's posterior mean surface against its function, which neither shares a profile with the other nor
    // comes from the model: exp(0.6x₀ + 0.8x₁), and sin(2(0.6x₀ + 0.8x₁)) + x₂.
    [Fact]
    public void SeveralOutputs_ShareTheDirectionalStructure_AndFollowEachFunction()
    {
        const int n = 200;
        var rng = new Xoshiro256PlusPlus(3627);
        var x = new double[4][];
        for (int j = 0; j < 4; j++) x[j] = new double[n];
        var y = new double[n, 2];
        var truth = new double[n, 2];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < 4; j++) x[j][i] = rng.NextDouble();
            double index = 0.6 * x[0][i] + 0.8 * x[1][i];
            truth[i, 0] = Math.Exp(index);
            truth[i, 1] = Math.Sin(2.0 * index) + x[2][i];
            for (int k = 0; k < 2; k++) y[i, k] = truth[i, k] + 0.05 * Normal.SampleStandard(rng);
        }
        var data = new BassData([.. Enumerable.Range(0, 4).Select(j => new ContinuousPredictor($"x{j}", x[j], new InputDomain(0.0, 1.0)))],
                                new MultiOutputResponse(y, ["a", "b"]));
        var config = new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 3, MaxBasis = 25, Directional = new DirectionalBasis() },
            Sampling = new BassSampling { Chains = 2, MasterSeed = 3627, BurnIn = 5000, Samples = 500, Thin = 2 },
        };
        BassFit fit = Bass.Fit(config, data);
        double[,] mean = BassMultiOutput.Predict(fit, data, PredictionCases.Of(data), new BassPredictOptions { Predictive = false, MeanSurface = true }).MeanSurface!;
        for (int k = 0; k < 2; k++)
        {
            double average = Enumerable.Range(0, n).Average(i => truth[i, k]);
            double spread = Math.Sqrt(Enumerable.Range(0, n).Sum(i => Math.Pow(truth[i, k] - average, 2)) / n);
            double error = Math.Sqrt(Enumerable.Range(0, n).Sum(i => Math.Pow(mean[i, k] - truth[i, k], 2)) / n);
            output.WriteLine($"output {k}: root mean square error against the function {error:F4}, the function's spread {spread:F3}");
            Assert.True(error <= 0.1 * spread, $"output {k}: {error:F4} against a spread of {spread:F3}");
        }
    }
}
