using System;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's active subspace and co-active pairs on a fit of several named outputs (E34): C = E[JᵀW_yJ] with H_ka = B_{a+1,k}, W_y K × K in
/// <see cref="MultiOutputResponse.Names"/>'s order, on the exact route and on the sampled route through E34's gradient.
/// <list type="bullet">
/// <item>Under W_y = aaᵀ for a non-axis-aligned a, C is the co-active C_fg of (a, a) on the same draws, which checks H's orientation,
/// W_y's off-diagonal and the outputs' order at once: the two exact calculations to rounding, the sampled against the exact within its
/// reported Monte Carlo error.</item>
/// <item>Two combinations of one fit through <see cref="BassFit.CoActiveSubspace(BassData, CoActiveSpec)"/> are
/// <see cref="BassMultiOutput.CoActive"/>'s, bit for bit, and each draw carries its weight once. Two fits pair by position with product
/// weights, a named combination on either side. The fixture's draw weights are unequal, so a squared same-fit weight shows.</item>
/// <item>At K = 1 under the identity, the readout is the scalar fit's on the same retained coefficients, bit for bit.</item>
/// <item>Declared locations, output patterns and a metric that is not K × K are refused, each with its reason.</item>
/// </list>
/// </summary>
public sealed class MultiOutputActiveSubspaceTests(ITestOutputHelper output)
{
    private static readonly double[] A = [0.6, -1.3];
    private static readonly double[] B = [0.3, 1.0];

    private static double[,] OuterOf(double[] a)
    {
        var w = new double[a.Length, a.Length];
        for (int i = 0; i < a.Length; i++) for (int j = 0; j < a.Length; j++) w[i, j] = a[i] * a[j];
        return w;
    }

    // A scalar fit holding the several-output fit's retained surfaces read through `functional`: its draws' coefficients are the
    // projections, its terms and weights the fit's, nothing refitted.
    private static (BassFit Fit, BassData Data) Projected(BassFit fit, BassData data, double[] functional)
    {
        var response = (MultiOutputResponse)data.Response;
        double[] y = [.. Enumerable.Range(0, response.Cases).Select(i => Enumerable.Range(0, response.Outputs).Sum(l => functional[l] * response[i, l]))];
        var scalarData = new BassData(data.Predictors, new ScalarResponse(y));
        var scalar = Bass.Fit(fit.Config with { OutputCovariance = null, Sigma2 = Sigma2Prior.InverseGamma(3, 1) }, scalarData, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        var store = new DrawStore<BassDraw>();
        foreach (var draw in fit.Draws[0])
            store.Add(draw with { Coefficients = BassMultiOutput.Project(draw, functional), OutputCoefficients = null, Sigma = null, Sigma2 = 1, OutputRates = null });
        return (scalar with { Terms = fit.Terms, Draws = [store], DrawWeights = fit.DrawWeights }, scalarData);
    }

    private static double Gap(double[,] x, double[,] y)
    {
        double scale = 0.0, gap = 0.0;
        for (int i = 0; i < x.GetLength(0); i++)
            for (int j = 0; j < x.GetLength(1); j++)
            {
                scale = Math.Max(scale, Math.Abs(y[i, j]));
                gap = Math.Max(gap, Math.Abs(x[i, j] - y[i, j]));
            }
        return gap / Math.Max(scale, double.Epsilon);
    }

    private static void SameBits(double[,] x, double[,] y) =>
        Assert.Equal(x.Cast<double>().Select(BitConverter.DoubleToInt64Bits), y.Cast<double>().Select(BitConverter.DoubleToInt64Bits));

    [Fact]
    public void UnderAaT_TheActiveMatrixIsTheCoActivePairOfAWithItself()
    {
        var (fit, data) = MultiOutputReadoutTests.Fixed();
        var metric = new OutputMetric(OuterOf(A));
        var exact = fit.ActiveSubspace(data, new ActiveSubspaceSpec { OutputMetric = metric, Route = ActiveSubspaceRoute.Exact });
        var pairs = fit.CoActiveSubspace(data, new CoActiveSpec { Output = A, OtherOutput = A, Route = ActiveSubspaceRoute.Exact });
        Assert.Equal(pairs.Posterior.Weights, exact.Posterior.Weights);
        double worst = 0.0;
        for (int d = 0; d < fit.Draws[0].Count; d++) worst = Math.Max(worst, Gap(exact.Posterior.Draws[d].C, pairs.Posterior.Pairs[d].Cfg));
        Assert.True(worst < 1e-12, $"exact: worst relative gap {worst:E2}");

        // The sampled route on shared points: against E34's sampled pair on the same points to rounding, and against the exact within
        // five of its reported standard errors.
        int[] draws = [0, 7, 39];
        var spec = new ActiveSubspaceSpec { OutputMetric = metric, Route = ActiveSubspaceRoute.Sampled, Points = 40_000, Seed = 20261219, Draws = draws };
        var sampled = fit.ActiveSubspace(data, spec);
        var sampledPairs = fit.CoActiveSubspace(data, new CoActiveSpec { Output = A, OtherOutput = A, Route = ActiveSubspaceRoute.Sampled, Points = 40_000, Seed = 20261219, Draws = draws });
        double sampledWorst = 0.0, z = 0.0;
        for (int s = 0; s < draws.Length; s++)
        {
            sampledWorst = Math.Max(sampledWorst, Gap(sampled.Posterior.Draws[s].C, sampledPairs.Posterior.Pairs[s].Cfg));
            double[,] c = sampled.Posterior.Draws[s].C, reference = pairs.Posterior.Pairs[draws[s]].Cfg, error = sampled.MonteCarloError![s];
            for (int i = 0; i < c.GetLength(0); i++)
                for (int j = 0; j < c.GetLength(1); j++)
                {
                    Assert.InRange(Math.Abs(c[i, j] - reference[i, j]), 0, 5 * error[i, j] + 1e-9);
                    if (error[i, j] > 0) z = Math.Max(z, Math.Abs(c[i, j] - reference[i, j]) / error[i, j]);
                }
        }
        Assert.True(sampledWorst < 1e-10, $"sampled: worst relative gap {sampledWorst:E2}");
        output.WriteLine($"C under aaᵀ against C_fg of (a, a): exact {worst:E2}, sampled {sampledWorst:E2}; sampled against exact, largest |z| {z:F2}");

        // Every draw at 40,000 points spans two of the gradient's draw batches (26 draws each at K = 2, two inputs); a draw reads the same
        // bits in either batch, in the active readout and in E34's pair alike.
        var whole = fit.ActiveSubspace(data, spec with { Draws = null });
        var wholePairs = fit.CoActiveSubspace(data, new CoActiveSpec { Output = A, OtherOutput = A, Route = ActiveSubspaceRoute.Sampled, Points = 40_000, Seed = 20261219 });
        for (int s = 0; s < draws.Length; s++)
        {
            SameBits(sampled.Posterior.Draws[s].C, whole.Posterior.Draws[draws[s]].C);
            SameBits(sampledPairs.Posterior.Pairs[s].Cfg, wholePairs.Posterior.Pairs[draws[s]].Cfg);
        }
    }

    [Fact]
    public void TwoCombinationsOfOneFit_AreE34sPair_AndEachDrawCarriesItsWeightOnce()
    {
        var (fit, data) = MultiOutputReadoutTests.Fixed();
        foreach (var route in new[] { ActiveSubspaceRoute.Exact, ActiveSubspaceRoute.Sampled })
        {
            var spec = new CoActiveSpec { Output = A, OtherOutput = B, Route = route, Points = 2_000, Seed = 20261114 };
            var e8 = fit.CoActiveSubspace(data, spec);
            var e34 = BassMultiOutput.CoActive(fit, data, spec);
            Assert.Equal(e34.Posterior.Weights, e8.Posterior.Weights);
            for (int d = 0; d < e34.Posterior.Pairs.Count; d++)
            {
                SameBits(e34.Posterior.Pairs[d].Cf, e8.Posterior.Pairs[d].Cf);
                SameBits(e34.Posterior.Pairs[d].Cg, e8.Posterior.Pairs[d].Cg);
                SameBits(e34.Posterior.Pairs[d].Cfg, e8.Posterior.Pairs[d].Cfg);
                if (route == ActiveSubspaceRoute.Sampled) SameBits(e34.MonteCarloError![d].Cfg, e8.MonteCarloError![d].Cfg);
            }
            SameBits(e34.Posterior.MeanCfg, e8.Posterior.MeanCfg);
            for (int d = 0; d < fit.Draws[0].Count; d++) Assert.Equal(fit.DrawWeights![d], e8.Posterior.Weights[d], 14);
        }
    }

    [Fact]
    public void TwoFits_PairByPosition_WithProductWeights_ANamedCombinationOnEitherSide()
    {
        var (fit, data) = MultiOutputReadoutTests.Fixed();
        var (scalar, scalarData) = Projected(fit, data, B);   // y_B = BᵀY's retained surfaces, draw for draw
        var oneFit = fit.CoActiveSubspace(data, new CoActiveSpec { Output = A, OtherOutput = B, Route = ActiveSubspaceRoute.Exact });
        var namedFirst = fit.CoActiveSubspace(data, scalar, scalarData, new CoActiveSpec { Output = A, Route = ActiveSubspaceRoute.Exact });
        var oneFitSwapped = fit.CoActiveSubspace(data, new CoActiveSpec { Output = B, OtherOutput = A, Route = ActiveSubspaceRoute.Exact });
        var namedSecond = scalar.CoActiveSubspace(scalarData, fit, data, new CoActiveSpec { OtherOutput = A, Route = ActiveSubspaceRoute.Exact });

        // The pairs are draw d of each fit, the same surfaces as the one fit's pair of A with B; only the weights differ: the product of
        // the two fits' weights, w_d², normalized, where one fit's draw carries w_d once.
        double[] weights = [.. fit.DrawWeights!];
        double squares = weights.Sum(w => w * w);
        for (int d = 0; d < fit.Draws[0].Count; d++)
        {
            SameBits(oneFit.Posterior.Pairs[d].Cfg, namedFirst.Posterior.Pairs[d].Cfg);
            SameBits(oneFit.Posterior.Pairs[d].Cf, namedFirst.Posterior.Pairs[d].Cf);
            SameBits(oneFitSwapped.Posterior.Pairs[d].Cfg, namedSecond.Posterior.Pairs[d].Cfg);
            SameBits(oneFitSwapped.Posterior.Pairs[d].Cg, namedSecond.Posterior.Pairs[d].Cg);
            double w = weights[d];
            Assert.Equal(w * w / squares, namedFirst.Posterior.Weights[d], 14);
            Assert.Equal(w * w / squares, namedSecond.Posterior.Weights[d], 14);
            Assert.Equal(w, oneFit.Posterior.Weights[d], 14);
        }
        Assert.NotEqual(oneFit.Posterior.Weights[0], namedFirst.Posterior.Weights[0], 6);

        // The sampled route reads the named side through E34's gradient and the scalar side through its own: the same surfaces, summed in
        // another order, so to rounding.
        int[] draws = [0, 7, 39];
        var sampledOne = fit.CoActiveSubspace(data, new CoActiveSpec { Output = A, OtherOutput = B, Route = ActiveSubspaceRoute.Sampled, Points = 2_000, Seed = 20261114, Draws = draws });
        var sampledTwo = fit.CoActiveSubspace(data, scalar, scalarData,
                                              new CoActiveSpec { Output = A, Route = ActiveSubspaceRoute.Sampled, Points = 2_000, Seed = 20261114, Draws = draws, OtherDraws = draws });
        double worst = 0.0;
        for (int s = 0; s < draws.Length; s++) worst = Math.Max(worst, Gap(sampledTwo.Posterior.Pairs[s].Cfg, sampledOne.Posterior.Pairs[s].Cfg));
        Assert.True(worst < 1e-10, $"sampled: worst relative gap {worst:E2}");
        output.WriteLine($"two fits against one, sampled: worst relative gap {worst:E2}");
    }

    [Fact]
    public void AtOneOutput_UnderTheIdentity_TheReadoutIsTheScalarFits_BitForBit()
    {
        var (fit, data) = MultiOutputReadoutTests.Fixed("ridge", outputs: 1);
        var (scalar, scalarData) = Projected(fit, data, [1.0]);   // the same retained coefficients: B·1 + 0
        var one = scalar.ActiveSubspace(scalarData, new ActiveSubspaceSpec { Route = ActiveSubspaceRoute.Exact, MatrixPoincare = true });
        foreach (OutputMetric? metric in new[] { null, new OutputMetric(new double[,] { { 1.0 } }) })
        {
            var several = fit.ActiveSubspace(data, new ActiveSubspaceSpec { OutputMetric = metric, Route = ActiveSubspaceRoute.Exact, MatrixPoincare = true });
            for (int d = 0; d < fit.Draws[0].Count; d++) SameBits(one.Posterior.Draws[d].C, several.Posterior.Draws[d].C);
            SameBits(one.Posterior.MeanC, several.Posterior.MeanC);
            SameBits(one.Poincare!.Total, several.Poincare!.Total);
            SameBits(one.Poincare.Bound, several.Poincare.Bound);
            Assert.Null(one.Poincare.JacobiRatio);   // hinge terms: neither side reads the Jacobi measure
            Assert.Equal(one.MatrixPoincare!.Ratios, several.MatrixPoincare!.Ratios);
        }
    }

    [Fact]
    public void NamedOutputs_RefuseLocationsPatternsAndAMetricNotKByK()
    {
        var (fit, data) = MultiOutputReadoutTests.Fixed(count: 2);
        var identity = new double[,] { { 1, 0 }, { 0, 1 } };
        var located = new OutputMetric(identity) { Locations = [new ContinuousColumn("t", [0.25, 0.75])] };
        var e = Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, new ActiveSubspaceSpec { OutputMetric = located }));
        Assert.Equal("metric", e.ParamName);
        Assert.Contains("named outputs, not output locations", e.Message);

        e = Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, new ActiveSubspaceSpec { OutputPatterns = 1 }));
        Assert.Equal("patterns", e.ParamName);
        Assert.Contains("E34", e.Message);

        e = Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, new ActiveSubspaceSpec { OutputMetric = new OutputMetric(new double[3, 3]) }));
        Assert.Contains("2 × 2, one row per output, in the response's Names order", e.Message);

        var (scalar, scalarData) = Projected(fit, data, B);
        e = Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, scalar, scalarData, new CoActiveSpec()));
        Assert.Equal("Output", e.ParamName);
        e = Assert.Throws<ArgumentException>(() => scalar.CoActiveSubspace(scalarData, fit, data, new CoActiveSpec { OtherOutput = [1.0] }));
        Assert.Equal("OtherOutput", e.ParamName);
    }
}
