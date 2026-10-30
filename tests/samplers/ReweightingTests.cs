using System;
using System.Linq;
using Maths.Samplers.FlatHistogram;
using Xunit;
using Xunit.Abstractions;
using static Maths.Samplers.Tests.FlatHistogramToys;

namespace Maths.Samplers.Tests;

/// <summary>
/// The reweighting and the ln Z(b) curve (N17 Phase 1, step 3), on the double well. One production run on the frozen
/// weights gives the tempered laws at every b, by the self-normalized weights e^{b·ℓ₁}/w(ℓ₁), and ln Z(b) by the ratio
/// estimator, which needs the prior normalized; both match the prior's quadrature within the production chains' errors.
/// At b = 1 the reweighted posterior holds both wells equally, which a chain at b = 1 cannot cross to.
/// </summary>
public sealed class ReweightingTests
{
    private readonly ITestOutputHelper _out;
    public ReweightingTests(ITestOutputHelper output) => _out = output;

    private const double Depth = 10.0;

    // Eight chains, so the spread-based errors have seven degrees of freedom; a 4.5-error bound fails about once in a
    // few hundred readings under Student's t.
    private const double Bound = 4.5;

    private static MulticanonicalResult<double> WellProduction(FlatHistogramWindow? window = null, int nodes = 31, int seed = 31)
    {
        FlatHistogramSpec spec = SamcLearningTests.WellSpec(64_000, 40_000, chains: 8, window: window) with { Nodes = nodes };
        return MulticanonicalRun.Run(spec, DoubleWell(Depth), new WellModel(Depth, 0.25), seed);
    }

    [Fact]
    public void Weigh_IsTheSelfNormalizedImportanceWeight_AndKishsEss()
    {
        var weights = new MulticanonicalWeights(-4.0, 0.0, 0.5, [0.0, 1.0, -1.0, 2.0, 0.5]);
        double[] levels0 = [-3.5, -0.25, -6.0], levels1 = [-1.0, 1.0, -2.0];
        var draws = new MulticanonicalDraws(weights, [levels0, levels1], [new double[3], new double[3]], 1);
        WeightedDraws w = Reweighting.Weigh(draws, 0.7);
        double[] all = [.. levels0, .. levels1];
        double[] raw = Array.ConvertAll(all, l => Math.Exp(0.7 * l - weights.LogWeight(l)));
        double total = raw.Sum();
        for (int s = 0; s < all.Length; s++) Assert.Equal(raw[s] / total, w.Weights[s], 14);
        Assert.Equal(1.0 / w.Weights.Sum(v => v * v), w.EffectiveSampleSize, 12);
        Assert.Equal(1.0, w.Weights.Sum(), 14);
        Assert.Equal(2, w.ChainEffectiveSampleSizes.Length);

        Assert.Equal(double.NegativeInfinity, Reweighting.LogImportance(double.NegativeInfinity, 0.0, 0.3));
        Assert.Equal(-2.0, Reweighting.LogImportance(double.NegativeInfinity, 2.0, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Reweighting.Weigh(draws, 1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => Reweighting.Expectation(draws, 0.5, 1));
    }

    [Fact]
    public void DoubleWell_ReweightedLaws_MatchQuadrature_AtThreeTemperatures()
    {
        MulticanonicalResult<double> run = WellProduction();
        MulticanonicalDraws draws = run.Draws!;
        foreach (double b in new[] { 0.1, 0.5, 1.0 })
        {
            (string Name, ReweightedEstimate Estimate, double Exact)[] checks =
            [
                ("E[ℓ₁]", Reweighting.Expectation(draws, b, 0), WellExpectation(Depth, b, x => WellLevel(x, Depth))),
                ("E[x²]", Reweighting.Expectation(draws, b, (_, f) => f[1] * f[1]), WellExpectation(Depth, b, x => x * x)),
                ("P(x > 0)", Reweighting.Expectation(draws, b, 2), WellExpectation(Depth, b, x => x > 0.0 ? 1.0 : 0.0)),
            ];
            WeightedDraws weighted = Reweighting.Weigh(draws, b);
            foreach ((string name, ReweightedEstimate e, double exact) in checks)
            {
                double z = (e.Mean - exact) / e.MonteCarloError;
                _out.WriteLine($"[reweight] b = {b}: {name} {e.Mean:F4} ± {e.MonteCarloError:F4} against {exact:F4} (z = {z:F2}); Kish ESS {weighted.EffectiveSampleSize:F0} of {weighted.Weights.Length}");
                Assert.True(Math.Abs(z) < Bound, $"b = {b}: {name} is {e.Mean:F4} ± {e.MonteCarloError:F4}, the quadrature's {exact:F4}");
            }
        }
        // Both wells, equally, at b = 1.
        Assert.Equal(0.5, WellExpectation(Depth, 1.0, x => x > 0.0 ? 1.0 : 0.0), 9);
    }

    [Fact]
    public void DoubleWell_LogEvidence_MatchesQuadrature_OnAGrid()
    {
        // The ratio estimator needs the prior normalized: the quadrature's masses sum to one.
        (_, _, double[] masses) = WellQuadrature(Depth, 400_000);
        Assert.Equal(0.0, LogSumExp(masses), 9);

        MulticanonicalResult<double> run = WellProduction(seed: 37);
        double[] grid = [0.0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0];
        LogEvidenceCurve curve = Reweighting.LogEvidence(run.Draws!, grid);
        Assert.Equal(0.0, curve.LogEvidence[0], 12);
        for (int g = 0; g < grid.Length; g++)
        {
            double exact = WellLogEvidence(Depth, grid[g]);
            double z = g == 0 ? 0.0 : (curve.LogEvidence[g] - exact) / curve.MonteCarloError[g];
            _out.WriteLine($"[evidence] b = {grid[g]}: ln Z {curve.LogEvidence[g]:F4} ± {curve.MonteCarloError[g]:F4} against {exact:F4} (z = {z:F2}); " +
                           $"ESS {curve.EffectiveSampleSize[g]:F0}");
            Assert.True(Math.Abs(z) < Bound, $"b = {grid[g]}: ln Z is {curve.LogEvidence[g]:F4} ± {curve.MonteCarloError[g]:F4}, the quadrature's {exact:F4}");
        }
        (double ratio, double error) = curve.LogRatio(2, grid.Length - 1);
        double exactRatio = WellLogEvidence(Depth, 1.0) - WellLogEvidence(Depth, 0.05);
        Assert.True(Math.Abs(ratio - exactRatio) < Bound * error, $"ln Z(1)/Z(0.05) is {ratio:F4} ± {error:F4}, the quadrature's {exactRatio:F4}");
        _out.WriteLine($"[evidence] denominator's Kish ESS {curve.NormalizerEffectiveSampleSize:F0} of {run.Draws!.Chains * run.Draws.Draws}");
    }

    [Fact]
    public void ANarrowWindow_StillReweighsExactly_ThroughItsExtensions()
    {
        // A window over ℓ₁ in [−12, −2] only: the extensions carry the draws below and above, and the laws stay exact.
        MulticanonicalResult<double> run = WellProduction(new FlatHistogramWindow(-12.0, -2.0), nodes: 11, seed: 41);
        foreach (double b in new[] { 0.5, 1.0 })
        {
            ReweightedEstimate e = Reweighting.Expectation(run.Draws!, b, 0);
            double exact = WellExpectation(Depth, b, x => WellLevel(x, Depth));
            _out.WriteLine($"[reweight] narrow window, b = {b}: E[ℓ₁] {e.Mean:F4} ± {e.MonteCarloError:F4} against {exact:F4}");
            Assert.True(Math.Abs(e.Mean - exact) < Bound * e.MonteCarloError);
        }
    }

    [Fact]
    public void OneChain_HasNoSpreadError()
    {
        FlatHistogramSpec spec = SamcLearningTests.WellSpec(2_000, 500, chains: 1);
        MulticanonicalResult<double> run = MulticanonicalRun.Run(spec, DoubleWell(Depth), new WellModel(Depth, 0.25), 3);
        Assert.True(double.IsNaN(Reweighting.Expectation(run.Draws!, 1.0, 0).MonteCarloError));
        LogEvidenceCurve curve = Reweighting.LogEvidence(run.Draws!, [0.5]);
        Assert.True(double.IsNaN(curve.MonteCarloError[0]));
        Assert.True(double.IsFinite(curve.LogEvidence[0]));
    }
}
