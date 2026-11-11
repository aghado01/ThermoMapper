using System;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 6, item 1: a basis fit paired with another fit, in the terms-and-H form, on the adequacy fixture's cases (a direct-mode fit
/// and PCA basis fits of rank 1 and 2 on the same two inputs):
/// <list type="bullet">
/// <item>a direct fit with a basis fit, and two basis fits: C_f, C_g and C_fg per pair against tensor quadrature of the two gradients,
/// the direct fit's from its readout, the basis fit's reconstructed from its components';</item>
/// <item>the two routes within the sampled route's Monte Carlo error;</item>
/// <item>a basis fit paired with itself is its two outputs read together, pair for pair bit for bit;</item>
/// <item>the refusals.</item>
/// </list>
/// </summary>
public sealed class CoActiveBasisPairTests(ITestOutputHelper output)
{
    private static (BassData Data, BassFit Direct, BassBasisFit Rank1, BassBasisFit Rank2) Fixture()
    {
        (BassData data, BassFit direct, BassPredictor[] predictors, double[,] y, _, _) = RankAndOutputSideTests.AdequacyFixture.Value;
        return (data, direct, Bases.Value.Rank1, Bases.Value.Rank2);
    }

    private static readonly Lazy<(BassBasisFit Rank1, BassBasisFit Rank2)> Bases = new(() =>
    {
        (_, _, BassPredictor[] predictors, double[,] y, _, _) = RankAndOutputSideTests.AdequacyFixture.Value;
        BassBasisFit Of(int rank) => BassComponents.Fit(BassDecomposition.FromPca(y, new BassPcaSpec { Rank = rank }), predictors, BassOutputComponentsTests.Config());
        return (Of(1), Of(2));
    });

    private static double[] OutputA(int m) => [.. Enumerable.Range(0, m).Select(l => 0.2 + 0.3 * Math.Sin(1.3 * l))];

    private static double[] OutputB(int m) => [.. Enumerable.Range(0, m).Select(l => l % 3 == 0 ? 1.0 : -0.4 + 0.05 * l)];

    // Each point's gradient along a, as a p-vector per point: a direct fit's draw from its readout, or a basis fit's joint draw from its
    // reconstructed Jacobian.
    private static double[][] DirectAlong(BassFit fit, BassData data, PredictionCases cases, int draw, double[] a)
    {
        BassGradient g = fit.Gradient(data, cases, [draw]);
        double[][,] j = FunctionalActiveSubspaceTests.DirectJacobians(g, 0);
        return [.. j.Select(jr => Along(jr, a))];
    }

    private static double[][] BasisAlong(BassBasisFit fit, PredictionCases cases, int joint, double[] a)
        => [.. FunctionalActiveSubspaceTests.BasisJacobians(fit, cases, joint).Select(jr => Along(jr, a))];

    private static double[] Along(double[,] j, double[] a)
    {
        int m = j.GetLength(0), p = j.GetLength(1);
        var v = new double[p];
        for (int i = 0; i < p; i++)
            for (int l = 0; l < m; l++) v[i] += a[l] * j[l, i];
        return v;
    }

    // Σ_r w_r u_r v_rᵀ.
    private static double[,] Cross(double[] w, double[][] u, double[][] v)
    {
        int p = u[0].Length;
        var c = new double[p, p];
        for (int r = 0; r < w.Length; r++)
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++) c[i, j] += w[r] * u[r][i] * v[r][j];
        return c;
    }

    /// <summary>
    /// A direct fit with the rank-2 basis fit, and the rank-1 basis fit with the rank-2 one: each pair's C_f, C_g and C_fg against tensor
    /// quadrature of the two sides' gradients along their functionals, split at the knots, within 10⁻¹⁰ of the largest entry.
    /// </summary>
    [Theory]
    [InlineData("direct and basis")]
    [InlineData("two bases")]
    public void Pairs_AreTheQuadratureOfTheirGradients(string kind)
    {
        (BassData data, BassFit direct, BassBasisFit rank1, BassBasisFit rank2) = Fixture();
        int m = data.Response.Locations;
        double[] a = OutputA(m), b = OutputB(m);
        bool basisFirst = kind == "two bases";
        int countF = basisFirst ? rank1.JointDraws.Count : direct.Draws.Sum(c => c.Count);
        int[] drawsF = ActiveSubspaceFollowUpTests.Spread(countF, 5), drawsG = ActiveSubspaceFollowUpTests.Spread(rank2.JointDraws.Count, 5);
        var spec = new CoActiveSpec { Output = a, OtherOutput = b, Draws = drawsF, OtherDraws = drawsG };
        CoActiveAnalysis pair = basisFirst ? rank1.CoActiveSubspace(rank2, spec) : direct.CoActiveSubspace(data, rank2, spec);
        BassInputs inputs = BassInputs.Prepare(data, direct.Transform);
        (PredictionCases cases, double[] w) = FunctionalActiveSubspaceTests.Tensor(inputs, pair.Measure!);
        double worst = 0.0;
        for (int s = 0; s < drawsF.Length; s++)
        {
            double[][] f = basisFirst ? BasisAlong(rank1, cases, drawsF[s], a) : DirectAlong(direct, data, cases, drawsF[s], a);
            double[][] g = BasisAlong(rank2, cases, drawsG[s], b);
            CoActivePair got = pair.Posterior.Pairs[s];
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(got.Cf, Cross(w, f, f)));
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(got.Cg, Cross(w, g, g)));
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(got.Cfg, Cross(w, f, g)));
        }
        output.WriteLine($"{kind}: {drawsF.Length} pairs; worst gap to quadrature {worst:E2}; concordance of the means {pair.Posterior.Concordance:F4}");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
        Assert.Equal(ActiveSubspaceRoute.Exact, pair.Route);
        Assert.Equal(drawsF, pair.DrawIndex);
        Assert.Equal(drawsG, pair.OtherDrawIndex);
    }

    /// <summary>The direct fit with the rank-2 basis fit: the exact and sampled routes, each entry of each pair's three matrices within
    /// 4.5 of the sampled route's standard errors.</summary>
    [Fact]
    public void DirectAndBasis_ExactAndSampledRoutesAgree()
    {
        (BassData data, BassFit direct, _, BassBasisFit rank2) = Fixture();
        int m = data.Response.Locations;
        var spec = new CoActiveSpec
        {
            Output = OutputA(m), OtherOutput = OutputB(m), Points = 8192, Seed = 3,
            Draws = ActiveSubspaceFollowUpTests.Spread(direct.Draws.Sum(c => c.Count), 8),
            OtherDraws = ActiveSubspaceFollowUpTests.Spread(rank2.JointDraws.Count, 8),
        };
        CoActiveAnalysis exact = direct.CoActiveSubspace(data, rank2, spec with { Route = ActiveSubspaceRoute.Exact });
        CoActiveAnalysis sampled = direct.CoActiveSubspace(data, rank2, spec with { Route = ActiveSubspaceRoute.Sampled });
        double worst = 0.0;
        for (int s = 0; s < exact.Posterior.Pairs.Count; s++)
        {
            CoActivePair e = exact.Posterior.Pairs[s], q = sampled.Posterior.Pairs[s];
            CoActiveError se = sampled.MonteCarloError![s];
            foreach ((double[,] x, double[,] y, double[,] z) in new[] { (e.Cf, q.Cf, se.Cf), (e.Cg, q.Cg, se.Cg), (e.Cfg, q.Cfg, se.Cfg) })
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++)
                    {
                        double gap = Math.Abs(x[i, j] - y[i, j]);
                        if (z[i, j] > 0.0) worst = Math.Max(worst, gap / z[i, j]);
                        else Assert.True(gap <= 1e-10 * Math.Max(1.0, Math.Abs(x[i, j])), $"pair {s}: a zero error with gap {gap:E2}");
                    }
        }
        output.WriteLine($"worst |exact − sampled| / se {worst:F2} over {exact.Posterior.Pairs.Count} pairs");
        Assert.True(worst <= 4.5, $"{worst:F2}");
    }

    /// <summary>A basis fit paired with itself, draw for draw, is its two outputs read together: each pair's matrices the same bits,
    /// the posterior's means within rounding of the weights' normalization.</summary>
    [Fact]
    public void BasisWithItself_IsItsTwoOutputs()
    {
        (BassData data, _, _, BassBasisFit rank2) = Fixture();
        int m = data.Response.Locations;
        int[] draws = ActiveSubspaceFollowUpTests.Spread(rank2.JointDraws.Count, 7);
        var spec = new CoActiveSpec { Output = OutputA(m), OtherOutput = OutputB(m), Draws = draws };
        CoActiveAnalysis pair = rank2.CoActiveSubspace(rank2, spec with { OtherDraws = draws });
        CoActiveAnalysis two = rank2.CoActiveSubspace(spec);
        for (int s = 0; s < draws.Length; s++)
            foreach ((double[,] x, double[,] y) in new[] { (pair.Posterior.Pairs[s].Cf, two.Posterior.Pairs[s].Cf), (pair.Posterior.Pairs[s].Cg, two.Posterior.Pairs[s].Cg),
                                                          (pair.Posterior.Pairs[s].Cfg, two.Posterior.Pairs[s].Cfg) })
                for (int i = 0; i < 2; i++)
                    for (int j = 0; j < 2; j++) Assert.Equal(BitConverter.DoubleToInt64Bits(y[i, j]), BitConverter.DoubleToInt64Bits(x[i, j]));
        Assert.Equal(two.Posterior.Concordance, pair.Posterior.Concordance, 14);
    }

    /// <summary>Refused: a basis side without its output functional, fits on other inputs, and selections of different lengths.</summary>
    [Fact]
    public void Pairs_RefuseWhatTheyDoNotRead()
    {
        (BassData data, BassFit direct, BassBasisFit rank1, BassBasisFit rank2) = Fixture();
        int m = data.Response.Locations;
        Assert.Throws<ArgumentException>(() => direct.CoActiveSubspace(data, rank2, new CoActiveSpec { Output = OutputA(m) }));
        Assert.Throws<ArgumentException>(() => rank1.CoActiveSubspace(rank2, new CoActiveSpec { OtherOutput = OutputB(m) }));
        BassBasisFit elsewhere = BassOutputComponentsTests.FitOf();
        int me = elsewhere.Decomposition.Locations;
        Assert.ThrowsAny<ArgumentException>(() => direct.CoActiveSubspace(data, elsewhere, new CoActiveSpec { Output = OutputA(m), OtherOutput = OutputB(me) }));
        Assert.Throws<ArgumentException>(() => rank1.CoActiveSubspace(rank2, new CoActiveSpec
        {
            Output = OutputA(m), OtherOutput = OutputB(m), Draws = [0, 1, 2], OtherDraws = [0, 1],
        }));
    }
}
