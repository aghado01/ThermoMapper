using System;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 8a, item A: co-active output functionals at declared locations (<see cref="CoActiveSpec.OutputLocations"/>,
/// <see cref="CoActiveSpec.OtherOutputLocations"/>) on a direct functional fit. The training locations declared give the pairs bit
/// for bit, on both routes, for two outputs of one fit and for two fits; at new locations, each side's C and the pair's C_fg meet tensor
/// quadrature of the gradient readout read there, the two outputs at different locations; the routes agree; at f = g the pair reduces
/// to C, the readout's under W = aaᵀ at the same locations; and what does not read locations refuses them.
/// </summary>
public sealed class CoActiveOutputLocationsTests(ITestOutputHelper output)
{
    private static readonly double[] T1 = [0.04, 0.3, 0.55, 0.81, 0.97];
    private static readonly double[] T2 = [0.12, 0.2, 0.47, 0.6, 0.73, 0.88, 0.99];

    private static ContinuousColumn[] At(double[] t) => [new ContinuousColumn("t", t)];

    private static double[] OutputA(int m) => [.. Enumerable.Range(0, m).Select(l => 0.2 + 0.3 * Math.Sin(1.3 * l))];

    private static double[] OutputB(int m) => [.. Enumerable.Range(0, m).Select(l => l % 3 == 0 ? 1.0 : -0.4 + 0.05 * l)];

    private static ContinuousColumn[] TrainingLocations(BassData data)
        => [.. ((FunctionalResponse)data.Response).OutputVariables.Select(v => new ContinuousColumn(v.Name, v.Values.ToArray()))];

    /// <summary>
    /// The training locations declared reproduce the pairs bit for bit, C_f, C_g and C_fg, on both routes: for two outputs of one fit with
    /// both sides declared and with one, and for the fit paired with itself as two fits.
    /// </summary>
    [Fact]
    public void TrainingLocationsDeclared_SameBits()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int m = data.Response.Locations;
        ContinuousColumn[] training = TrainingLocations(data);
        int[] sel = Spread(fit.Draws.Sum(c => c.Count), 6);
        int compared = 0;
        foreach (ActiveSubspaceRoute route in new[] { ActiveSubspaceRoute.Exact, ActiveSubspaceRoute.Sampled })
        {
            var spec = new CoActiveSpec { Draws = sel, Output = OutputA(m), OtherOutput = OutputB(m), Route = route, Points = 512 };
            CoActiveAnalysis bare = fit.CoActiveSubspace(data, spec);
            CoActiveAnalysis both = fit.CoActiveSubspace(data, spec with { OutputLocations = training, OtherOutputLocations = training });
            CoActiveAnalysis one = fit.CoActiveSubspace(data, spec with { OtherOutputLocations = training });
            CoActiveAnalysis twoBare = fit.CoActiveSubspace(data, fit, data, spec with { OtherDraws = sel });
            CoActiveAnalysis twoBoth = fit.CoActiveSubspace(data, fit, data, spec with { OtherDraws = sel, OutputLocations = training, OtherOutputLocations = training });
            compared += SameBits(bare, both) + SameBits(bare, one) + SameBits(twoBare, twoBoth);
        }
        output.WriteLine($"{compared} pairs the same bits with the training locations declared and not, on both routes");
    }

    /// <summary>
    /// Two outputs of one fit at new locations, a over five and b over seven, none of them a training location: C_f, C_g and C_fg
    /// against tensor quadrature of aᵀJ at the first five and bᵀJ at the other seven, J from the gradient readout read there; the fit
    /// paired with itself as two fits gives the same pairs.
    /// </summary>
    [Fact]
    public void NewLocations_MeetTheGradientReadoutThere()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        double[] a = OutputA(T1.Length), b = OutputB(T2.Length);
        int[] sel = Spread(fit.Draws.Sum(c => c.Count), 6);
        var spec = new CoActiveSpec { Draws = sel, Output = a, OtherOutput = b, OutputLocations = At(T1), OtherOutputLocations = At(T2) };
        CoActiveAnalysis one = fit.CoActiveSubspace(data, spec);
        CoActiveAnalysis two = fit.CoActiveSubspace(data, fit, data, spec with { OtherDraws = sel });
        (PredictionCases cases, double[] w) = FunctionalActiveSubspaceTests.Tensor(BassInputs.Prepare(data, fit.Transform), one.Measure!);
        BassGradient gf = BassGradients.Of(fit.ModelFor(data), fit.Terms, fit.Draws, fit.DrawWeights, cases, sel, At(T1));
        BassGradient gg = BassGradients.Of(fit.ModelFor(data), fit.Terms, fit.Draws, fit.DrawWeights, cases, sel, At(T2));
        double worst = 0.0;
        for (int s = 0; s < sel.Length; s++)
        {
            int ss = s;
            Func<int, int, double> fa = (r, i) => Along(a, gf, ss, r, i), fb = (r, i) => Along(b, gg, ss, r, i);
            CoActivePair pair = one.Posterior.Pairs[s];
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(pair.Cfg, Cross(w, fa, fb)));
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(pair.Cf, Cross(w, fa, fa)));
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(pair.Cg, Cross(w, fb, fb)));
            Assert.Equal(pair.Cfg, two.Posterior.Pairs[s].Cfg);
        }
        output.WriteLine($"{sel.Length} pairs at {T1.Length} and {T2.Length} new locations: worst gap to quadrature {worst:E2}; concordance of the means {one.Posterior.Concordance:F4}");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
    }

    /// <summary>At the new locations the exact and sampled routes agree within 4.5 of the sampled route's standard errors.</summary>
    [Fact]
    public void NewLocations_ExactAndSampledRoutesAgree()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        var spec = new CoActiveSpec
        {
            Draws = Spread(fit.Draws.Sum(c => c.Count), 10), Points = 8192, Seed = 11,
            Output = OutputA(T1.Length), OtherOutput = OutputB(T2.Length), OutputLocations = At(T1), OtherOutputLocations = At(T2),
        };
        double worst = WorstInErrors(fit.CoActiveSubspace(data, spec with { Route = ActiveSubspaceRoute.Exact }),
                                     fit.CoActiveSubspace(data, spec with { Route = ActiveSubspaceRoute.Sampled }));
        output.WriteLine($"worst |exact − sampled| / se {worst:F2}");
        Assert.True(worst <= 4.5, $"{worst:F2}");
    }

    /// <summary>
    /// At f = g, one output functional at the same new locations on both sides: C_fg = C_f = C_g within rounding, concordance 1, and C_f
    /// the readout's C under W = aaᵀ with the output metric's locations the same.
    /// </summary>
    [Fact]
    public void AtFEqualsG_ReducesToC()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        double[] a = OutputA(T1.Length);
        int[] sel = Spread(fit.Draws.Sum(c => c.Count), 6);
        CoActiveAnalysis aa = fit.CoActiveSubspace(data, new CoActiveSpec { Draws = sel, Output = a, OtherOutput = a, OutputLocations = At(T1), OtherOutputLocations = At(T1) });
        ActiveSubspaceAnalysis waa = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Draws = sel, OutputMetric = new OutputMetric(Outer(a)) { Locations = At(T1) } });
        double worst = 0.0;
        for (int s = 0; s < sel.Length; s++)
        {
            CoActivePair pair = aa.Posterior.Pairs[s];
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(pair.Cfg, pair.Cf));
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(pair.Cg, pair.Cf));
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(pair.Cf, waa.Posterior.Draws[s].C));
            Assert.Equal(1.0, pair.Concordance, 12);
        }
        output.WriteLine($"a = b at {T1.Length} new locations: worst gap among C_fg, C_f, C_g and the readout's C under aaᵀ {worst:E2}");
        Assert.True(worst <= 1e-12, $"{worst:E2}");
    }

    /// <summary>
    /// Refused: locations outside t's box, a column naming no output variable, a functional of the wrong length for its locations,
    /// locations on a scalar fit's side, on a basis fit (paired with a direct fit, with another basis fit, or two of its outputs), and on a
    /// fit of several named outputs, whose outputs have no locations.
    /// </summary>
    [Fact]
    public void Refusals()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        var spec = new CoActiveSpec { Draws = [0], Output = OutputA(2), OtherOutput = OutputB(2), OutputLocations = At([0.2, 0.5]), OtherOutputLocations = At([0.3, 0.6]) };
        fit.CoActiveSubspace(data, spec);
        Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, spec with { OutputLocations = At([0.2, 1.2]) }));
        Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, spec with { OtherOutputLocations = [new ContinuousColumn("u", [0.3, 0.6])] }));
        Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, spec with { Output = OutputA(3) }));
        Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, spec with { OtherOutput = OutputB(data.Response.Locations) }));
        Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, fit, data, spec with { OtherDraws = [0], OtherOutputLocations = At([0.3, -0.1]) }));

        (BassData scalarData, BassFit scalar) = CoActiveFixture.OnFunctionalInputs();
        fit.CoActiveSubspace(data, scalar, scalarData, new CoActiveSpec { Draws = [0], OtherDraws = [0], Output = OutputA(2), OutputLocations = At([0.2, 0.5]) });
        Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, scalar, scalarData,
                                                                    new CoActiveSpec { Draws = [0], OtherDraws = [0], Output = OutputA(2), OutputLocations = At([0.2, 0.5]), OtherOutputLocations = At([0.3]) }));

        BassBasisFit basis = BassOutputComponentsTests.FitOf();
        int mb = basis.Decomposition.Locations;
        Assert.Throws<ArgumentException>(() => fit.CoActiveSubspace(data, basis, new CoActiveSpec
        {
            Draws = [0], OtherDraws = [0], Output = OutputA(data.Response.Locations), OtherOutput = OutputB(mb), OtherOutputLocations = At([0.2, 0.5]),
        }));
        Assert.Throws<ArgumentException>(() => basis.CoActiveSubspace(basis, new CoActiveSpec
        {
            Draws = [0], OtherDraws = [0], Output = OutputA(mb), OtherOutput = OutputB(mb), OutputLocations = At([0.2, 0.5]),
        }));
        Assert.Throws<ArgumentException>(() => basis.CoActiveSubspace(new CoActiveSpec { Draws = [0], Output = OutputA(mb), OtherOutput = OutputB(mb), OtherOutputLocations = At([0.2]) }));

        (BassFit named, BassData namedData) = MultiOutputReadoutTests.Fixed();
        int k = namedData.Response.Outputs;
        var pair = new CoActiveSpec { Draws = [0], Output = OutputA(k), OtherOutput = OutputB(k) };
        named.CoActiveSubspace(namedData, pair);
        Assert.Throws<ArgumentException>(() => named.CoActiveSubspace(namedData, pair with { OutputLocations = At([0.2]) }));
        Assert.Throws<ArgumentException>(() => BassMultiOutput.CoActive(named, namedData, pair with { OtherOutputLocations = At([0.2]) }));
        Assert.Throws<ArgumentException>(() => named.CoActiveSubspace(namedData, named, namedData, pair with { OtherDraws = [0], OutputLocations = At([0.2]) }));
    }

    // ── Helpers ──────────────────────────────────────────────────────────────────────────────────────────────────────

    private static int[] Spread(int count, int k) => [.. Enumerable.Range(0, k).Select(i => (int)((long)i * (count - 1) / (k - 1)))];

    // Each pair's three matrices the same bits in two readouts; the pairs compared.
    private static int SameBits(CoActiveAnalysis x, CoActiveAnalysis y)
    {
        Assert.Equal(x.Posterior.Pairs.Count, y.Posterior.Pairs.Count);
        for (int s = 0; s < x.Posterior.Pairs.Count; s++)
        {
            CoActivePair p = x.Posterior.Pairs[s], q = y.Posterior.Pairs[s];
            foreach ((double[,] u, double[,] v) in new[] { (p.Cf, q.Cf), (p.Cg, q.Cg), (p.Cfg, q.Cfg) })
                for (int i = 0; i < u.GetLength(0); i++)
                    for (int j = 0; j < u.GetLength(1); j++)
                        Assert.Equal(BitConverter.DoubleToInt64Bits(u[i, j]), BitConverter.DoubleToInt64Bits(v[i, j]));
        }
        return x.Posterior.Pairs.Count;
    }

    // Σ_r w_r ∂_i f(x_r) ∂_j g(x_r) over the rule's points, two inputs.
    private static double[,] Cross(double[] w, Func<int, int, double> gradF, Func<int, int, double> gradG)
    {
        var c = new double[2, 2];
        for (int r = 0; r < w.Length; r++)
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++) c[i, j] += w[r] * gradF(r, i) * gradG(r, j);
        return c;
    }

    // aᵀ∂y/∂x_i at point r of draw s.
    private static double Along(double[] a, BassGradient g, int s, int r, int i)
    {
        double v = 0.0;
        for (int l = 0; l < a.Length; l++) v += a[l] * g.Partials[i][s, r, l];
        return v;
    }

    private static double[,] Outer(double[] a)
    {
        var o = new double[a.Length, a.Length];
        for (int i = 0; i < a.Length; i++)
            for (int j = 0; j < a.Length; j++) o[i, j] = a[i] * a[j];
        return o;
    }

    // The worst |exact − sampled|/se over each pair's three matrices, where the error is positive.
    private static double WorstInErrors(CoActiveAnalysis exact, CoActiveAnalysis sampled)
    {
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
        return worst;
    }
}
