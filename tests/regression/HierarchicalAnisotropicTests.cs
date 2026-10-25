using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// <see cref="HierarchicalPSpline"/> with anisotropic smoothing: the population surface's x- and y-penalties each get
/// their own smoothing variance. On replicate surfaces gentle in x and fast in y, every scheme puts τ_0x² below τ_0y²
/// and the schemes agree on both; each draw's τ_0² is the components' geometric mean; a penalty without components is
/// refused; and a resumed run continues exactly.
/// </summary>
public sealed class HierarchicalAnisotropicTests
{
    private readonly ITestOutputHelper _out;
    public HierarchicalAnisotropicTests(ITestOutputHelper output) => _out = output;

    private static double Gaussian(Xoshiro256PlusPlus rng)
    {
        double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    private static (TensorDesign Td, double[][] Replicates) Surfaces(int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        const int g = 16;
        var xs = new double[g * g];
        var ys = new double[g * g];
        for (int a = 0; a < g; a++)
            for (int b = 0; b < g; b++) { xs[a * g + b] = (a + 0.5) / g; ys[a * g + b] = (b + 0.5) / g; }
        var reps = new double[3][];
        for (int r = 0; r < 3; r++)
            reps[r] = xs.Select((x, i) => Math.Sin(1.5 * Math.PI * x) * Math.Sin(5.0 * Math.PI * ys[i]) + 0.05 * r + 0.1 * Gaussian(rng))
                        .ToArray();
        var knots = new KnotConfig(Enumerable.Range(1, 5).Select(k => k / 6.0).ToArray());
        return (TensorProductDesign.Build(new SplineBasis(3), knots, xs, new SplineBasis(3), knots, ys), reps);
    }

    [Fact]
    public void EveryScheme_SmoothsXMoreThanY_AndTheSchemesAgree()
    {
        (TensorDesign td, double[][] reps) = Surfaces(seed: 41);
        var penalty = new TensorPenalty(td.NuX, td.NuY);
        var medians = new List<(HierarchicalScheme Scheme, double X, double Y)>();
        foreach (HierarchicalScheme scheme in Enum.GetValues<HierarchicalScheme>())
        {
            var model = scheme == HierarchicalScheme.ParameterExpanded
                ? new HierarchicalPSpline(td.Design, reps, penalty, deviationPrior: new HalfTDeviationPrior(), anisotropicSmoothing: true)
                : new HierarchicalPSpline(td.Design, reps, penalty, anisotropicSmoothing: true);
            HierarchicalResult h = model.Run(chains: 2, masterSeed: 6, burn: 600, samples: 900, scheme: scheme);
            PosteriorInterval tx = h.ComponentSmoothingVarianceInterval(0), ty = h.ComponentSmoothingVarianceInterval(1);
            _out.WriteLine($"[hps-aniso] {scheme}: τ0x² {tx.Median:E2} [{tx.Lo:E2}, {tx.Hi:E2}] τ0y² {ty.Median:E2} [{ty.Lo:E2}, {ty.Hi:E2}] " +
                           $"R̂ {string.Join(", ", h.ComponentSmoothingVarianceConvergence!.Select(c => c.RHat.Value.ToString("F3")))}");
            Assert.True(tx.Hi < ty.Lo, $"{scheme}: the x-variance's interval lies below the y-variance's");
            Assert.Equal(2, h.ComponentSmoothingVarianceConvergence!.Count);
            foreach (HierarchicalDraw d in h.PosteriorDraws[0])
                Assert.Equal(Math.Sqrt(d.ComponentSmoothingVariance![0] * d.ComponentSmoothingVariance[1]), d.PopulationSmoothingVariance, 10);
            medians.Add((scheme, tx.Median, ty.Median));
        }
        double x0 = medians[0].X, y0 = medians[0].Y;
        foreach ((HierarchicalScheme scheme, double x, double y) in medians)
        {
            Assert.True(Math.Abs(Math.Log(x / x0)) < 0.4, $"{scheme}: τ0x² median {x:E2} against the centred {x0:E2}");
            Assert.True(Math.Abs(Math.Log(y / y0)) < 0.4, $"{scheme}: τ0y² median {y:E2} against the centred {y0:E2}");
        }
    }

    [Fact]
    public void AnisotropicSmoothing_RefusesAPenaltyWithoutComponents_AndResumesExactly()
    {
        double[] x = Enumerable.Range(0, 40).Select(i => (i + 0.5) / 40).ToArray();
        double[][] curves = [x.Select(v => Math.Sin(6 * v)).ToArray(), x.Select(v => Math.Sin(6 * v) + 0.1).ToArray()];
        double[,] z = new SplineBasis(3).Design(new KnotConfig([0.25, 0.5, 0.75]), x);
        Assert.Throws<ArgumentException>(() => new HierarchicalPSpline(z, curves, new DifferencePenalty(2), anisotropicSmoothing: true));
        Assert.Throws<InvalidOperationException>(() => new HierarchicalPSpline(z, curves, new DifferencePenalty(2))
            .Run(chains: 1, burn: 10, samples: 10).ComponentSmoothingVarianceInterval(0));

        (TensorDesign td, double[][] reps) = Surfaces(seed: 42);
        var model = new HierarchicalPSpline(td.Design, reps, new TensorPenalty(td.NuX, td.NuY), anisotropicSmoothing: true);
        HierarchicalResult straight = model.Run(chains: 2, masterSeed: 5, burn: 100, samples: 150);
        HierarchicalResult first = model.Run(chains: 2, masterSeed: 5, burn: 100, samples: 50, checkpoint: true);
        HierarchicalResult resumed = model.Run(chains: 2, masterSeed: 5, burn: 100, samples: 100, resumeFrom: first.Checkpoint);
        Assert.Equal(straight.PopulationCoefficients, resumed.PopulationCoefficients);
        for (int c = 0; c < 2; c++)
            for (int i = 0; i < straight.PosteriorDraws[c].Count; i++)
                Assert.Equal(straight.PosteriorDraws[c][i].ComponentSmoothingVariance, resumed.PosteriorDraws[c][i].ComponentSmoothingVariance);
    }
}
