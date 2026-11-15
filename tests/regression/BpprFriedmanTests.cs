using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E36's slice 3 exploration on Friedman's function as the paper's §3.1 runs it (exploration, not validation; the BPPR brief): n = 300 cases
/// uniform on [0, 1]⁶, x₆ inert, noise sd 1; 20,000 sweeps, the last 2,000 kept; 2,000 points out of sample. The paper reports an
/// out-of-sample RMSE of 1.05 for the noisy response and 0.38 for f, and 95% intervals covering 94.9% of the responses and 91.1% of f, from
/// one run of <c>BayesPPR</c> at its defaults. The engine runs those defaults in its own terms: K = 4, A = 3, ω₀ = υ₀ = 1, p₀ = 2/3,
/// n_active = 20, κ = 1000, Zellner–Siow (τ = 1/g with g ~ IG(½, n/2), so τ ~ Gamma(½, rate n/2)), Jeffreys σ², the intercept inside
/// the g-prior, M capped at min(150, ⌊n/K⌋ − 2) = 73, and the count law at h₁ = 10⁷, h₂ = 10⁶, Poisson(10) to 10⁻⁶ in its ratios. One
/// run each, three seeds, so the readings are reported beside the paper's, with bounds that catch only a gross failure.
/// </summary>
public sealed class BpprFriedmanTests(ITestOutputHelper output)
{
    [Theory]
    [InlineData(3671)]
    [InlineData(3681)]
    [InlineData(3691)]
    public void Friedman_ReadBesideThePapers(int seed)
    {
        ProjectionPursuitReference reference = ProjectionPursuitReferences.Friedman();
        const int n = 300, held = 2000, p = 6;
        var rng = new Xoshiro256PlusPlus(seed);
        double[][] Draw(int count) => [.. Enumerable.Range(0, p).Select(_ => Enumerable.Range(0, count).Select(_ => rng.NextDouble()).ToArray())];
        double F(double[][] x, int i) => reference.Function.Value([.. Enumerable.Range(0, p).Select(j => x[j][i])]);
        double[][] x = Draw(n);
        double[] y = [.. Enumerable.Range(0, n).Select(i => F(x, i) + Normal.SampleStandard(rng))];
        double[][] xt = Draw(held);
        double[] ft = [.. Enumerable.Range(0, held).Select(i => F(xt, i))];
        double[] yt = [.. ft.Select(f => f + Normal.SampleStandard(rng))];
        var data = new BassData([.. Enumerable.Range(0, p).Select(j => new ContinuousPredictor($"x{j + 1}", x[j]))], new ScalarResponse(y));
        var config = new BassConfig
        {
            Structure = new BassStructureSpec
            {
                MaxInteraction = 3, MaxBasis = 73, OrderWeight = 1.0, VariableWeight = 1.0, Directional = new DirectionalBasis { Columns = 4 },
            },
            Intercept = GPriorCentre.Zero,
            Sigma2 = Sigma2Prior.Jeffreys,
            Coefficients = BassCoefficientPrior.G(TauPrior.Gamma(0.5, 0.5 * n)),
            Count = new BassCountPrior(1e7, 1e6),
            Sampling = new BassSampling { Chains = 1, MasterSeed = seed, BurnIn = 18_000, Samples = 2_000, Thin = 1 },
        };
        var clock = System.Diagnostics.Stopwatch.StartNew();
        BassFit fit = Bass.Fit(config, data);
        double seconds = clock.Elapsed.TotalSeconds;
        PredictionCases cases = PredictionCases.Of([.. Enumerable.Range(0, p).Select(j => (BassPredictor)new ContinuousPredictor($"x{j + 1}", xt[j]))]);
        BassPrediction prediction = fit.Predict(data, cases, new BassPredictOptions { Seed = seed + 1 });
        double[,,] latent = prediction.Latent, predictive = prediction.Predictive!;
        int draws = latent.GetLength(0);
        double squaredY = 0.0, squaredF = 0.0;
        int coveredY = 0, coveredF = 0;
        var column = new double[draws];
        (double Low, double High) Interval(double[,,] source, int i)
        {
            for (int d = 0; d < draws; d++) column[d] = source[d, i, 0];
            Array.Sort(column);
            return (Quantile(column, 0.025), Quantile(column, 0.975));
        }
        for (int i = 0; i < held; i++)
        {
            double mean = 0.0;
            for (int d = 0; d < draws; d++) mean += latent[d, i, 0];
            mean /= draws;
            squaredY += (mean - yt[i]) * (mean - yt[i]);
            squaredF += (mean - ft[i]) * (mean - ft[i]);
            (double lowY, double highY) = Interval(predictive, i);
            (double lowF, double highF) = Interval(latent, i);
            if (yt[i] >= lowY && yt[i] <= highY) coveredY++;
            if (ft[i] >= lowF && ft[i] <= highF) coveredF++;
        }
        double rmseY = Math.Sqrt(squaredY / held), rmseF = Math.Sqrt(squaredF / held);
        double coverY = (double)coveredY / held, coverF = (double)coveredF / held;
        BassDraw[] kept = [.. fit.Draws.SelectMany(c => c)];
        double sigma = kept.Average(d => Math.Sqrt(d.Sigma2));
        double meanM = kept.Average(d => d.Structure.TermCount(4));
        output.WriteLine($"seed {seed}: Friedman, n = {n}, 20,000 sweeps in {seconds:F1} s; mean M {meanM:F1} ridges; posterior mean σ {sigma:F3} (paper's 1.05)");
        output.WriteLine($"out of sample, {held} points: RMSE for y {rmseY:F3} (paper's 1.05), for f {rmseF:F3} (paper's 0.38)");
        output.WriteLine($"95% intervals: prediction intervals cover {coverY:P1} of y (paper's 94.9%), credible intervals {coverF:P1} of f (paper's 91.1%)");
        output.WriteLine($"moves: {string.Join(", ", fit.Diagnostics.MoveKinds.Entries.Where(e => e.Attempts > 0).Select(e => $"{e.Kind} {e.Accepted}/{e.Attempts} (refused {e.Refused})"))}; "
                         + $"refusals {string.Join(", ", fit.Diagnostics.Refusals.Entries.Select(e => $"{e.Category} {e.Count}"))}; M's R̂-free range over the kept draws "
                         + $"{kept.Min(d => d.Structure.TermCount(4))} to {kept.Max(d => d.Structure.TermCount(4))}");
        Assert.True(rmseF < 1.0, $"RMSE for f {rmseF:F3}");
        Assert.True(rmseY < 1.5, $"RMSE for y {rmseY:F3}");
        Assert.InRange(coverY, 0.85, 1.0);
        Assert.InRange(coverF, 0.6, 1.0);
    }

    // R's type-7 quantile of sorted values.
    private static double Quantile(double[] sorted, double probability)
    {
        double index = (sorted.Length - 1) * probability;
        int lo = (int)Math.Floor(index);
        double h = index - lo;
        return lo + 1 < sorted.Length ? sorted[lo] + h * (sorted[lo + 1] - sorted[lo]) : sorted[lo];
    }
}
