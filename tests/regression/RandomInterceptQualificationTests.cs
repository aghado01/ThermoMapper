using System;
using System.Diagnostics;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Random-effect intercepts (E21) recover group effects from a generator outside the model, behind <see cref="BassQualification"/>; without
/// it each returns at once. Each dataset puts n = 246 cases on twelve unbalanced levels (n_g = 4 + 3g), x uniform on [0, 1]², the surface
/// sin(2πx₀) + x₁ and Gaussian noise of σ 0.3. The levels' effects are twelve draws of 0.8·(E − 1), E ~ Exp(1), sorted so the larger
/// levels take the larger effects: skewed, tied to the levels' sizes, and not centred on the observations, where the block's prior is a
/// centred Normal independent of them. The block declares s = 2, near the generator's spread of 0.8/0.3 but not it.
/// <list type="bullet">
/// <item><b>Recovery:</b> each level's posterior mean against its true effect under the convention, u*_g − Σ_h n_h u*_h / n, as standardized
/// errors pooled over four datasets' 48 levels, and the 90% intervals' coverage of them.</item>
/// <item><b>The convention under level centring:</b> the intercept carries the observation-weighted mean of the true effects, so the
/// fit's level over the training cases sits at the surface's mean there plus Σ n_g u*_g / n; the levels' unweighted mean, which the
/// effects' tie to the sizes sets apart, is the alternative it must beat.</item>
/// <item><b>What the block buys:</b> σ against the generator's; the predictive mean at held-out cases of the observed levels against their
/// true mean, sin(2πx₀) + x₁ + u*_g, beside the same fit without the block, whose noise absorbs the effects; and the population surface,
/// which gains only in variance there, not made worse.</item>
/// </list>
/// A hinge fit (BASS's defaults: degree 1, the g-form with τ drawn) and a Legendre fit (KHAOS's: degree up to 8, order up to 2, the
/// modified g-prior, the zero centre, Jeffreys σ², inputs declared on [0, 1]), each on the same four datasets.
/// </summary>
public sealed class RandomInterceptQualificationTests(ITestOutputHelper output)
{
    private const int Levels = 12, Datasets = 4, Held = 400;
    private const double Sigma = 0.3, Scale = 2.0;

    private sealed record Generated(BassData Data, double[] Effects, double[] Surface, PredictionCases Population, PredictionCases Grouped,
                                    double[] TestSurface, double[] TestMean);

    private static double Truth(double x0, double x1) => Math.Sin(2 * Math.PI * x0) + x1;

    private static Generated Generate(int seed, bool declared)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        int[] sizes = [.. Enumerable.Range(0, Levels).Select(g => 4 + 3 * g)];
        int n = sizes.Sum();
        int[] codes = [.. sizes.SelectMany((size, g) => Enumerable.Repeat(g, size)).OrderBy(_ => rng.NextDouble())];
        double[] effects = [.. Enumerable.Range(0, Levels).Select(_ => 0.8 * (-Math.Log(1.0 - rng.NextDouble()) - 1.0)).Order()];
        var x0 = new double[n];
        var x1 = new double[n];
        var surface = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            surface[i] = Truth(x0[i], x1[i]);
            y[i] = surface[i] + effects[codes[i]] + Sigma * Maths.Distributions.Normal.SampleStandard(rng);
        }
        string[] names = [.. Enumerable.Range(0, Levels).Select(g => $"G{g}")];
        InputDomain? domain = declared ? new InputDomain(0.0, 1.0) : null;
        var data = new BassData([new ContinuousPredictor("x0", x0, domain), new ContinuousPredictor("x1", x1, domain)], new ScalarResponse(y),
                                new CategoricalPredictor("level", codes, names));
        // Held-out cases, each level equally likely, so the small levels count as much as the large.
        double[] t0 = [.. Enumerable.Range(0, Held).Select(_ => rng.NextDouble())], t1 = [.. Enumerable.Range(0, Held).Select(_ => rng.NextDouble())];
        int[] t = [.. Enumerable.Range(0, Held).Select(_ => (int)(rng.NextDouble() * Levels))];
        var population = new PredictionCases([new ContinuousColumn("x0", t0), new ContinuousColumn("x1", t1)]);
        var grouped = new PredictionCases([.. population.Columns, new CategoricalColumn("level", [.. t.Select(g => names[g])])]);
        double[] testSurface = [.. t0.Zip(t1, Truth)];
        return new Generated(data, effects, surface, population, grouped, testSurface, [.. testSurface.Select((f, i) => f + effects[t[i]])]);
    }

    private static BassConfig Config(string basis, int seed, int cases) => basis == "hinge"
        ? new BassConfig { Sampling = new BassSampling { Chains = 2, MasterSeed = seed, BurnIn = 3_000, Samples = 1_000 } }
        : KhaosReachQualificationTests.Config(null, cases, seed);

    [Theory]
    [InlineData("hinge")]
    [InlineData("legendre")]
    public void GroupEffects_AreRecovered_UnderTheConvention(string basis)
    {
        if (!BassQualification.Enabled) return;
        int covered = 0, intervals = 0;
        double z2 = 0.0, levelGap = 0.0, alternativeGap = 0.0;
        TimeSpan withTime = TimeSpan.Zero, withoutTime = TimeSpan.Zero;   // the cost row's wall clock, indicative
        for (int k = 0; k < Datasets; k++)
        {
            Generated g = Generate(9100 + k, declared: basis == "legendre");
            BassData data = g.Data;
            int n = data.Response.Cases;
            BassConfig plain = Config(basis, 300 + k, n), block = plain with { RandomIntercept = new BassRandomIntercept(Scale) };
            var watch = Stopwatch.StartNew();
            BassFit fit = Bass.Fit(block, data, new BassRunOptions { MaxDegreeOfParallelism = 2 });
            withTime += watch.Elapsed;
            BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];

            // The convention's target: the true effects less their observation-weighted mean, which the level takes.
            int[] codes = data.Grouping!.Codes.ToArray();
            double[] counts = [.. Enumerable.Range(0, Levels).Select(l => (double)codes.Count(c => c == l))];
            double weighted = Enumerable.Range(0, Levels).Sum(l => counts[l] * g.Effects[l]) / n, unweighted = g.Effects.Average();
            double[] target = [.. g.Effects.Select(u => u - weighted)];
            RandomInterceptBlock fitted = BassModel.Create(block, data).RandomIntercept!;
            var mean = new double[Levels];
            double worstCentre = 0.0, squares = 0.0;
            for (int l = 0; l < Levels; l++)
            {
                int observed = fitted.ObservedOf(l);
                double[] u = [.. draws.Select(d => d.GroupEffects![observed]).Order()];
                mean[l] = u.Average();
                double sd = Math.Sqrt(u.Sum(v => (v - mean[l]) * (v - mean[l])) / (u.Length - 1));
                z2 += Math.Pow((mean[l] - target[l]) / sd, 2);
                squares += Math.Pow(mean[l] - target[l], 2);
                double lo = u[(int)Math.Floor(0.05 * (u.Length - 1))], hi = u[(int)Math.Ceiling(0.95 * (u.Length - 1))];
                if (target[l] >= lo && target[l] <= hi) covered++;
                intervals++;
            }
            foreach (BassDraw d in draws)
                worstCentre = Math.Max(worstCentre, Math.Abs(Enumerable.Range(0, Levels).Sum(l => counts[l] * d.GroupEffects![fitted.ObservedOf(l)])));
            double rmse = Math.Sqrt(squares / Levels), correlation = Correlation(mean, target);

            // The level over the training cases: the population fit's mean there, draw by draw.
            BassPrediction population = fit.Predict(data, PredictionCases.Of(data.Predictors), new BassPredictOptions { Predictive = false });
            double[] levels = [.. Enumerable.Range(0, population.DrawCount).Select(d => Enumerable.Range(0, n).Average(i => population.Latent[d, i, 0]))];
            double level = levels.Average(), levelSd = Math.Sqrt(levels.Sum(v => (v - level) * (v - level)) / (levels.Length - 1));
            double levelTarget = g.Surface.Average() + weighted, levelAlternative = g.Surface.Average() + unweighted;
            levelGap += Math.Pow(level - levelTarget, 2);
            alternativeGap += Math.Pow(level - levelAlternative, 2);

            // What the block buys, beside the fit without it.
            BassData bare = new(data.Predictors, data.Response);
            watch.Restart();
            BassFit without = Bass.Fit(plain, bare, new BassRunOptions { MaxDegreeOfParallelism = 2 });
            withoutTime += watch.Elapsed;
            double sigmaWith = Math.Sqrt(draws.Average(d => d.Sigma2)), sigmaWithout = Math.Sqrt(without.Draws.SelectMany(c => c).Average(d => d.Sigma2));
            double groupedWith = Error(fit, data, g.Grouped, g.TestMean, 0.0), groupedWithout = Error(without, bare, g.Population, g.TestMean, 0.0);
            double surfaceWith = Error(fit, data, g.Population, g.TestSurface, weighted), surfaceWithout = Error(without, bare, g.Population, g.TestSurface, weighted);

            output.WriteLine($"{basis}, dataset {k}: the effects' RMSE {rmse:F3}, their spread {Math.Sqrt(target.Sum(t => t * t) / Levels):F3}, correlation {correlation:F4}, " +
                             $"Σ n_g u_g at most {worstCentre:E1}; level {level:F4} ± {levelSd:F4} against the convention's {levelTarget:F4} " +
                             $"(the unweighted mean's {levelAlternative:F4}); σ {sigmaWith:F3} against {Sigma} (without the block {sigmaWithout:F3}); " +
                             $"held-out cases of observed levels, RMSE {groupedWith:F4} (without the block {groupedWithout:F4}); " +
                             $"the population surface {surfaceWith:F4} (without {surfaceWithout:F4})");
            Assert.True(correlation > 0.9, $"dataset {k}: correlation {correlation:F4}");
            Assert.True(worstCentre < 1e-9, $"dataset {k}: Σ n_g u_g = {worstCentre:E2}");
            Assert.True(Math.Abs(level - levelTarget) < 4 * levelSd, $"dataset {k}: level {level:F4} ± {levelSd:F4} against {levelTarget:F4}");
            Assert.True(Math.Abs(sigmaWith / Sigma - 1.0) < 0.15, $"dataset {k}: σ {sigmaWith:F3}");
            Assert.True(groupedWith < 0.5 * groupedWithout, $"dataset {k}: held-out RMSE {groupedWith:F4} with the block, {groupedWithout:F4} without");
            Assert.True(surfaceWith < 1.1 * surfaceWithout, $"dataset {k}: the surface's RMSE {surfaceWith:F4} with the block, {surfaceWithout:F4} without");
        }
        double coverage = (double)covered / intervals;
        output.WriteLine($"{basis}: the 90% intervals cover {covered} of {intervals} levels' effects ({coverage:F3}); the standardized errors' mean square " +
                         $"{z2 / intervals:F3}; the level's squared gaps summed, {levelGap:F5} from the convention's target and {alternativeGap:F5} from the unweighted mean's; " +
                         $"the fits took {withTime.TotalSeconds:F1} s with the block and {withoutTime.TotalSeconds:F1} s without");
        Assert.True(coverage >= 0.8, $"coverage {coverage:F3}");
        Assert.True(z2 / intervals < 2.0, $"the standardized errors' mean square {z2 / intervals:F3}");
        Assert.True(levelGap < alternativeGap, $"{levelGap:F5} against {alternativeGap:F5}");
    }

    // The predictive mean's RMSE at held-out cases against their true mean plus the level the convention adds to it.
    private static double Error(BassFit fit, BassData data, PredictionCases held, double[] truth, double level)
    {
        double[,] mean = fit.Predict(data, held, new BassPredictOptions { Predictive = false }).Mean();
        return Math.Sqrt(Enumerable.Range(0, truth.Length).Sum(i => Math.Pow(mean[i, 0] - truth[i] - level, 2)) / truth.Length);
    }

    private static double Correlation(double[] a, double[] b)
    {
        double ma = a.Average(), mb = b.Average(), sab = 0.0, saa = 0.0, sbb = 0.0;
        for (int i = 0; i < a.Length; i++)
        {
            sab += (a[i] - ma) * (b[i] - mb);
            saa += (a[i] - ma) * (a[i] - ma);
            sbb += (b[i] - mb) * (b[i] - mb);
        }
        return sab / Math.Sqrt(saa * sbb);
    }
}
