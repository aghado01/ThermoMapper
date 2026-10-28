using System;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Evidence;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BARS's evidence (<see cref="BarsResult.Evidence"/>): stepping-stone over each chain's traced ladder, on by default for a
/// run that can give one. The truth comes from outside the ladder. With at most one knot the knot prior is a count and
/// one uniform position, restricted to the support, so ln p(y) = ln Z₁ − ln Z₀ is a quadrature over the knot's position;
/// with no knots it is the empty design's marginal.
/// </summary>
public sealed class BarsEvidenceTests
{
    private readonly ITestOutputHelper _out;
    public BarsEvidenceTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void NoKnots_TheEvidenceIsTheEmptyDesignsMarginal()
    {
        (double[] x, double[] y) = Data();
        var model = new WeightedNormalModel();
        BarsResult r = Tempered(model, new UniformPrior(0)).Run(x, y, grid: x, chains: 2, masterSeed: 5, burn: 200, samples: 300);
        double exact = model.LogMarginalLikelihood(new SplineBasis(3).Design(new KnotConfig([]), x), y, null);
        Assert.NotNull(r.Evidence);
        Assert.Equal(exact, r.Evidence!.LogEvidence, 1e-9);
        Assert.Equal(0.0, r.Evidence.MonteCarloError);
    }

    [Fact]
    public void OneKnot_TheEvidenceMatchesQuadratureOverTheKnot()
    {
        (double[] x, double[] y) = Data();
        var model = new WeightedNormalModel();
        var prior = new UniformPrior(1);
        double truth = QuadratureEvidence(model, prior, x, y, points: 20_000);
        BarsResult r = Tempered(model, prior).Run(x, y, grid: x, chains: 2, masterSeed: 11, burn: 1000, samples: 4000);
        BarsEvidence e = r.Evidence!;
        _out.WriteLine($"quadrature {truth:F5}; BARS {e.LogEvidence:F5} ± {e.MonteCarloError:F5}; chains " +
                       string.Join(", ", Array.ConvertAll([.. e.Chains], c => $"{c.LogRatio:F5} ± {c.MonteCarloError:F5}")));
        Assert.Equal(2, e.Chains.Count);
        Assert.True(e.MonteCarloError > 0.0);
        Assert.True(Math.Abs(e.LogEvidence - truth) < 4.0 * e.MonteCarloError, $"BARS {e.LogEvidence:F5} against {truth:F5}");
        foreach (LadderEstimate c in e.Chains)
            Assert.True(Math.Abs(c.LogRatio - truth) < 4.0 * c.MonteCarloError, $"chain {c.LogRatio:F5} against {truth:F5}");
    }

    [Fact]
    public void Evidence_IsReportedWhereTheRunCanGiveIt_AndTheKnobTurnsItOff()
    {
        var model = new WeightedNormalModel();
        var prior = new PoissonPrior(4.0);
        var kernel = new LocalBetaKernel(50.0);
        var basis = new SplineBasis(3);
        Assert.True(Tempered(model, prior).ReportsEvidence);
        Assert.False(Tempered(model, prior, evidence: false).ReportsEvidence);
        Assert.False(new BarsEnsemble(basis, model, prior, kernel).ReportsEvidence);                        // untempered
        Assert.False(new BarsEnsemble(basis, model, prior, kernel, temperLevels: 4).ReportsEvidence);       // no b = 0 rung
        Assert.False(new BarsEnsemble(basis, model, prior, kernel, temperLevels: 4, temperPath: TemperingPath.Powered,
                                      temperReferenceRung: true).ReportsEvidence);                          // not the prior at b = 0
        var jeffreys = new WeightedNormalModel(form: TemperedForm.Linear, sigma2Prior: Sigma2Prior.Jeffreys);
        Assert.False(Tempered(jeffreys, prior).ReportsEvidence);                                             // a constant omitted

        (double[] x, double[] y) = Data();
        Assert.Null(new BarsEnsemble(basis, model, prior, kernel).Run(x, y, grid: x, chains: 1, masterSeed: 2, burn: 50, samples: 50).Evidence);
        Assert.Null(Tempered(model, new UniformPrior(0), evidence: false).Run(x, y, grid: x, chains: 1, masterSeed: 2, burn: 50, samples: 50).Evidence);
    }

    private static BarsEnsemble Tempered(IObservationModel model, IComplexityPrior prior, bool evidence = true)
        => new(new SplineBasis(3), model, prior, new LocalBetaKernel(50.0), temperLevels: 5, temperReferenceRung: true,
               evidence: evidence);

    // ln Z₁ − ln Z₀ for a count prior on {0, 1}: Z_b = p(0)·e^{ℓ_b(∅)} + p(1)·∫ 1[t on the support] e^{ℓ_b(t)} dt, by the
    // midpoint rule in t. The target's own evaluation says which positions the support gate and admission keep.
    private static double QuadratureEvidence(IObservationModel model, IComplexityPrior prior, double[] x, double[] y, int points)
    {
        var target = new SplineTarget(new SplineBasis(3), model, prior, x, y);
        ITargetEvaluation empty = target.Evaluate(new KnotConfig([]));
        double h = 1.0 / points;
        var terms = new double[points + 1];
        terms[0] = empty.LogPrior + empty.LogLikelihood(1.0);
        double kept = 0.0;
        for (int i = 0; i < points; i++)
        {
            ITargetEvaluation e = target.Evaluate(new KnotConfig([(i + 0.5) * h]));
            double l = double.IsNegativeInfinity(e.LogPrior) ? double.NegativeInfinity : e.LogLikelihood(1.0);
            terms[i + 1] = double.IsNegativeInfinity(l) ? double.NegativeInfinity : e.LogPrior + l + Math.Log(h);
            if (!double.IsNegativeInfinity(l)) kept += h;
        }
        double max = double.NegativeInfinity;
        foreach (double t in terms) max = Math.Max(max, t);
        double sum = 0.0;
        foreach (double t in terms) sum += Math.Exp(t - max);
        double logZ0 = Math.Log(Math.Exp(prior.LogPrior(0)) + Math.Exp(prior.LogPrior(1)) * kept);
        return max + Math.Log(sum) - logZ0;
    }

    private static (double[] X, double[] Y) Data()
    {
        var rng = new Xoshiro256PlusPlus(seed: 21);
        const int n = 60;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.2 * Normal.SampleStandard(rng);
        }
        return (x, y);
    }
}
