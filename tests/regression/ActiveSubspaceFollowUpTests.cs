using System;
using System.Collections.Generic;
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
/// E8's slice 6 gates (the active-subspace brief's follow-ups, items 2 to 6; the pairs of fits and the Legendre moments have their own
/// classes):
/// <list type="bullet">
/// <item>the co-active projection mean on pairs whose co-active subspaces and signs are known, on an opposed pair, and under weights;</item>
/// <item>output locations beyond the training ones: the training ones declared reproduce the readout bit for bit, new ones meet tensor
/// quadrature of the gradient readout there, the routes agree, and what is refused;</item>
/// <item>the matrix check's Gram cache is pure: every draw's ratio the same bits read alone or among the others;</item>
/// <item>Poincaré constants for truncated normals: the smaller of σ² and (b − a)²/π² taken, the Gaussian equality on a linear factor, a
/// seeded sweep of the two cited bounds, no violation on fitted draws, and the laws with no constant;</item>
/// <item>activity scores on known spectra and through the readout.</item>
/// </list>
/// </summary>
public sealed class ActiveSubspaceFollowUpTests(ITestOutputHelper output)
{
    // ── The co-active projection mean ────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Twenty pairs on p = 5 sharing an eigenbasis Q, each scaled by its own size: 70% have V's top two by |λ| on q₁ (+3) and q₂ (−2.5),
    /// 30% on q₁ (+3) and q₃ (+2.2), with C_f = C_g = Q|Λ|Qᵀ. Then P̄ = q₁q₁ᵀ + 0.7q₂q₂ᵀ + 0.3q₃q₃ᵀ, the mean subspace is span(q₁, q₂),
    /// the first group sits on it and the second at 90° in one angle, the selected eigenvalues keep their signs, and the directional
    /// concordance is 1 along q₁ and Σw(−2.5s)/Σw(2.5s) mixed with the second group's +0.3s along q₂.
    /// </summary>
    [Fact]
    public void CoActiveProjection_KnownSubspacesAndSigns()
    {
        const int p = 5, n = 20;
        double[,] q = RandomOrthogonal(p, 20261110);
        var rng = new Xoshiro256PlusPlus(5);
        var pairs = new CoActivePair[n];
        var scales = new double[n];
        for (int s = 0; s < n; s++)
        {
            bool first = s % 10 < 7;
            double[] lambda = first ? [3.0, -2.5, 0.4, -0.2, 0.1] : [3.0, 0.3, 2.2, -0.2, 0.1];
            scales[s] = 0.5 + rng.NextDouble();
            double[] scaled = [.. lambda.Select(l => l * scales[s])];
            pairs[s] = new CoActivePair(Spectral(q, [.. scaled.Select(Math.Abs)]), Spectral(q, [.. scaled.Select(Math.Abs)]), Spectral(q, scaled));
        }
        CoActiveProjection mean = CoActivePosterior.Of(pairs, [.. Enumerable.Repeat(1.0, n)]).ProjectionMean(2);
        double[] expected = [1.0, 0.7, 0.3, 0.0, 0.0];
        double eig = expected.Select((e, j) => Math.Abs(mean.Eigenvalues[j] - e)).Max();
        double span = SubspacePosterior.PrincipalAngles(Columns(mean.Vectors, 2), Columns(q, 2)).Max();
        double angles = 0.0, selected = 0.0;
        for (int s = 0; s < n; s++)
        {
            bool first = s % 10 < 7;
            double[] want = first ? [0.0, 0.0] : [0.0, Math.PI / 2];
            angles = Math.Max(angles, want.Select((a, j) => Math.Abs(mean.AnglesToMean[s][j] - a)).Max());
            double[] signed = first ? [3.0, -2.5] : [3.0, 2.2];
            for (int j = 0; j < 2; j++) selected = Math.Max(selected, Math.Abs(mean.Selected[s, j] - signed[j] * scales[s]) / scales[s]);
        }
        double num = 0.0, den = 0.0;
        for (int s = 0; s < n; s++)
        {
            double along = s % 10 < 7 ? -2.5 : 0.3;
            num += along * scales[s];
            den += Math.Abs(along) * scales[s];
        }
        output.WriteLine($"P̄'s eigenvalues {string.Join(", ", mean.Eigenvalues.Select(v => v.ToString("F15")))}; mean subspace to span(q₁, q₂) {span:E2} rad; " +
                         $"angles to the known ones {angles:E2}; selected λ {selected:E2}; concordances {mean.Concordances[0]:F15}, {mean.Concordances[1]:F15} " +
                         $"(want 1 and {num / den:F15})");
        Assert.True(eig <= 1e-14, $"{eig:E2}");
        Assert.True(span <= 1e-12 && angles <= 1e-12 && selected <= 1e-14, $"{span:E2}, {angles:E2}, {selected:E2}");
        Assert.True(Math.Abs(mean.Concordances[0] - 1.0) <= 1e-13, $"{mean.Concordances[0]:F15}");
        Assert.True(Math.Abs(mean.Concordances[1] - num / den) <= 1e-13, $"{mean.Concordances[1]:F15}");
    }

    /// <summary>
    /// An opposed pair, g = −f, so C_fg = −C_f: the |λ| order is C_f's own, so P̄ is the projection mean of the C_f draws
    /// (<see cref="SubspacePosterior"/>), every selected eigenvalue is negative, and every directional concordance is −1. Weighted pairs
    /// read as replicated ones.
    /// </summary>
    [Fact]
    public void CoActiveProjection_OpposedPairs_AndWeights()
    {
        const int p = 4, n = 9, k = 2;
        var rng = new Xoshiro256PlusPlus(20261111);
        var cf = new double[n][,];
        var pairs = new CoActivePair[n];
        for (int s = 0; s < n; s++)
        {
            double[] lambda = [4.0 + rng.NextDouble(), 2.0 + rng.NextDouble(), 0.5 * rng.NextDouble(), 0.1 * rng.NextDouble()];
            cf[s] = Spectral(RandomOrthogonal(p, 100 + s / 3), lambda);
            pairs[s] = new CoActivePair(cf[s], cf[s], Negated(cf[s]));
        }
        double[] weights = [.. Enumerable.Range(0, n).Select(s => 1.0 + s % 3)];
        CoActiveProjection mean = CoActivePosterior.Of(pairs, weights).ProjectionMean(k);
        double[,] reference = SubspacePosterior.Of(cf, weights, k).Projection;
        double gap = MaxGap(mean.Projection, reference);
        Assert.True(gap <= 1e-14, $"{gap:E2}");
        for (int s = 0; s < n; s++)
            for (int j = 0; j < k; j++) Assert.True(mean.Selected[s, j] < 0.0);
        Assert.All(mean.Concordances, c => Assert.True(Math.Abs(c + 1.0) <= 1e-13, $"{c:F15}"));

        // Weights as replicates: pair s taken 1 + s % 3 times.
        var replicated = new List<CoActivePair>();
        for (int s = 0; s < n; s++)
            for (int r = 0; r < 1 + s % 3; r++) replicated.Add(pairs[s]);
        double[,] again = CoActivePosterior.Of(replicated, [.. Enumerable.Repeat(1.0, replicated.Count)]).ProjectionMean(k).Projection;
        double weighted = MaxGap(mean.Projection, again);
        output.WriteLine($"opposed: P̄ against the C_f draws' projection mean {gap:E2}; weighted against replicated {weighted:E2}");
        Assert.True(weighted <= 1e-14, $"{weighted:E2}");
        Assert.Throws<ArgumentOutOfRangeException>(() => CoActivePosterior.Of(pairs, weights).ProjectionMean(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => CoActivePosterior.Of(pairs, weights).ProjectionMean(p + 1));
    }

    // ── Output locations ─────────────────────────────────────────────────────────────────────────────────────────────

    private static ContinuousColumn[] TrainingLocations(BassData data)
        => [.. ((FunctionalResponse)data.Response).OutputVariables.Select(v => new ContinuousColumn(v.Name, v.Values.ToArray()))];

    /// <summary>The training locations declared as an output metric's locations reproduce the default readout bit for bit, C and the
    /// output patterns, on both routes.</summary>
    [Fact]
    public void OutputLocations_TheTrainingOnes_ReproduceTheReadout()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int m = data.Response.Locations;
        double[,] w = FunctionalActiveSubspaceTests.RandomMetric(m, 61);
        int[] draws = Spread(fit.Draws.Sum(c => c.Count), 8);
        foreach (ActiveSubspaceRoute route in new[] { ActiveSubspaceRoute.Exact, ActiveSubspaceRoute.Sampled })
        {
            var spec = new ActiveSubspaceSpec { Draws = draws, Route = route, Points = 512, OutputMetric = new OutputMetric(w), OutputPatterns = 2 };
            ActiveSubspaceAnalysis a = fit.ActiveSubspace(data, spec);
            ActiveSubspaceAnalysis b = fit.ActiveSubspace(data, spec with { OutputMetric = new OutputMetric(w) { Locations = TrainingLocations(data) } });
            for (int d = 0; d < draws.Length; d++) AssertSame(a.Posterior.Draws[d].C, b.Posterior.Draws[d].C);
            Assert.Equal(a.Outputs!.MeanEigenvalues, b.Outputs!.MeanEigenvalues);
        }
    }

    /// <summary>
    /// At nine new locations in t's box, one of them a training location (1/6): C against tensor quadrature of the gradient readout's
    /// Jacobian read there, under the identity and a full W; the readout at the training location is its column at the training
    /// locations, bit for bit; and at W = I, tr E[G] = tr E[C] (slice 5's identity) over the new locations.
    /// </summary>
    [Theory]
    [InlineData("identity")]
    [InlineData("full")]
    public void OutputLocations_New_AreTheQuadratureOfTheJacobianThere(string metric)
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        double[] t = [0.03, 1.0 / 6.0, 0.27, 0.42, 0.58, 0.66, 0.79, 0.93, 0.995];
        ContinuousColumn[] at = [new ContinuousColumn("t", t)];
        int m = t.Length;
        double[,] w = metric == "identity" ? FunctionalActiveSubspaceTests.Identity(m) : FunctionalActiveSubspaceTests.RandomMetric(m, 67);
        int[] draws = Spread(fit.Draws.Sum(c => c.Count), 6);
        var spec = new ActiveSubspaceSpec { Draws = draws, Dimension = 2, OutputMetric = new OutputMetric(w) { Locations = at }, OutputPatterns = 2 };
        ActiveSubspaceAnalysis a = fit.ActiveSubspace(data, spec);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        (PredictionCases cases, double[] weights) = FunctionalActiveSubspaceTests.Tensor(inputs, a.Measure!);
        BassGradient g = BassGradients.Of(fit.ModelFor(data), fit.Terms, fit.Draws, fit.DrawWeights, cases, draws, at);
        double worst = 0.0;
        for (int d = 0; d < draws.Length; d++)
            worst = Math.Max(worst, FunctionalActiveSubspaceTests.RelativeGap(a.Posterior.Draws[d].C,
                                                                              FunctionalActiveSubspaceTests.Quadrature(FunctionalActiveSubspaceTests.DirectJacobians(g, d), weights, w)));
        BassGradient training = fit.Gradient(data, cases, draws);
        for (int i = 0; i < g.Partials.Count; i++)
            for (int d = 0; d < draws.Length; d++)
                for (int r = 0; r < cases.Cases; r++) Assert.Equal(training.Partials[i][d, r, 1], g.Partials[i][d, r, 1]);
        double traces = 0.0;
        if (metric == "identity")
        {
            double c = 0.0;
            for (int i = 0; i < 2; i++) c += a.Posterior.MeanC[i, i];
            traces = Math.Abs(a.Outputs!.MeanTrace - c) / c;
        }
        output.WriteLine($"{metric}: {m} new locations, {draws.Length} draws; worst gap to quadrature {worst:E2}; tr E[G] − tr E[C] {traces:E2}");
        Assert.True(worst <= 1e-10, $"{worst:E2}");
        Assert.True(traces <= 1e-12, $"{traces:E2}");
        Assert.Equal(0, a.Poincare!.Violations);
    }

    /// <summary>At the new locations, the exact and sampled routes agree within 4.5 of the sampled route's standard errors.</summary>
    [Fact]
    public void OutputLocations_New_ExactAndSampledRoutesAgree()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        double[] t = [0.05, 0.31, 0.5, 0.74, 0.88];
        var spec = new ActiveSubspaceSpec
        {
            Draws = Spread(fit.Draws.Sum(c => c.Count), 10), Points = 8192, Seed = 7,
            OutputMetric = new OutputMetric(FunctionalActiveSubspaceTests.RandomMetric(t.Length, 71)) { Locations = [new ContinuousColumn("t", t)] },
        };
        ActiveSubspaceAnalysis exact = fit.ActiveSubspace(data, spec with { Route = ActiveSubspaceRoute.Exact });
        ActiveSubspaceAnalysis sampled = fit.ActiveSubspace(data, spec with { Route = ActiveSubspaceRoute.Sampled });
        double worst = 0.0;
        for (int d = 0; d < exact.Posterior.Draws.Count; d++)
            for (int i = 0; i < 2; i++)
                for (int j = i; j < 2; j++)
                {
                    double se = sampled.MonteCarloError![d][i, j], gap = Math.Abs(exact.Posterior.Draws[d].C[i, j] - sampled.Posterior.Draws[d].C[i, j]);
                    if (se > 0.0) worst = Math.Max(worst, gap / se);
                    else Assert.True(gap <= 1e-10 * Math.Max(1.0, Math.Abs(exact.Posterior.Draws[d].C[i, j])), $"draw {d}: a zero error with gap {gap:E2}");
                }
        output.WriteLine($"worst |exact − sampled| / se {worst:F2} over {exact.Posterior.Draws.Count} draws at {t.Length} new locations");
        Assert.True(worst <= 4.5, $"{worst:F2}");
    }

    /// <summary>Refused: locations outside t's box, a metric of another size, a column that names no output variable, and locations
    /// asked of a basis fit, whose basis exists only at its training locations, or of basis adequacy.</summary>
    [Fact]
    public void OutputLocations_RefuseWhatTheyDoNotRead()
    {
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        OutputMetric At(double[] t, string name = "t") => new(FunctionalActiveSubspaceTests.Identity(t.Length)) { Locations = [new ContinuousColumn(name, t)] };
        Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, new ActiveSubspaceSpec { OutputMetric = At([0.2, 1.2]) }));
        Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, new ActiveSubspaceSpec { OutputMetric = At([0.2, 0.5], "u") }));
        Assert.Throws<ArgumentException>(() => fit.ActiveSubspace(data, new ActiveSubspaceSpec
        {
            OutputMetric = new OutputMetric(FunctionalActiveSubspaceTests.Identity(3)) { Locations = [new ContinuousColumn("t", [0.2, 0.5])] },
        }));
        BassBasisFit basis = BassOutputComponentsTests.FitOf();
        int mb = basis.Decomposition.Locations;
        var training = new OutputMetric(FunctionalActiveSubspaceTests.Identity(mb)) { Locations = [new ContinuousColumn("t", [.. Enumerable.Range(0, mb).Select(l => l / (mb - 1.0))])] };
        Assert.Throws<ArgumentException>(() => basis.ActiveSubspace(new ActiveSubspaceSpec { OutputMetric = training }));
        (BassData ad, BassFit direct, _, _, _, _) = RankAndOutputSideTests.AdequacyFixture.Value;
        Assert.Throws<ArgumentException>(() => basis.Adequacy(direct, ad, new ActiveSubspaceSpec { OutputMetric = training }));
    }

    // ── The matrix check's Gram cache ────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The cache is pure: every draw's ratio is the same bits whether the draw is read alone, where nothing is cached, or among every
    /// draw of the fit, where its structure's Grams may come from an earlier draw: a scalar fit and a direct functional fit under a full W.
    /// </summary>
    [Fact]
    public void MatrixCheck_GramCache_IsPure()
    {
        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        BassData data = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        var fullW = new OutputMetric(FunctionalActiveSubspaceTests.RandomMetric(data.Response.Locations, 73));
        foreach ((string name, BassFit model, BassData on, OutputMetric? w) in new[] { ("scalar", f, df, (OutputMetric?)null), ("functional", fit, data, fullW) })
        {
            BassDraw[] pooled = [.. model.Draws.SelectMany(c => c)];
            int[] draws = [.. Enumerable.Range(0, pooled.Length)];
            IReadOnlyList<double> all = model.ActiveSubspace(on, new ActiveSubspaceSpec { Draws = draws, MatrixPoincare = true, OutputMetric = w }).MatrixPoincare!.Ratios;
            for (int d = 0; d < draws.Length; d++)
            {
                double alone = model.ActiveSubspace(on, new ActiveSubspaceSpec { Draws = [d], MatrixPoincare = true, OutputMetric = w }).MatrixPoincare!.Ratios[0];
                Assert.Equal(BitConverter.DoubleToInt64Bits(alone), BitConverter.DoubleToInt64Bits(all[d]));
            }
            int structures = pooled.Select(p => string.Join(",", p.Structure.TermIds)).Distinct().Count();
            output.WriteLine($"{name}: {draws.Length} draws over {structures} structures, every ratio the same bits alone and among them");
            Assert.True(structures < draws.Length, $"{structures} structures in {draws.Length} draws: the cache would not be exercised");
        }
    }

    // ── Poincaré constants for truncated normals ─────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The per-input bound reads the constant chosen: for x0 a normal truncated to [0.1, 0.9] of the box, min(σ², 0.8²/π²) in box
    /// units, at σ = 0.08 (σ² binds) and at σ = 2 (the interval binds), Bound/C_vv equal to it; x1 uniform keeps 1/π². Neither check
    /// finds a violation, and the matrix check is computed.
    /// </summary>
    [Theory]
    [InlineData(0.08)]
    [InlineData(2.0)]
    public void TruncatedNormal_TakesTheSmallerConstant(double sd)
    {
        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        var scale = (RangeScale)f.Transform.Predictors[0];
        double width = scale.Max - scale.Min, lower = scale.Min + 0.1 * width, upper = scale.Min + 0.9 * width;
        var spec = new ActiveSubspaceSpec
        {
            MatrixPoincare = true,
            Measures = new Dictionary<string, InputMeasure>
            {
                [f.Transform.Predictors[0].Name] = new MixtureMeasure([new NormalComponent(1.0, scale.Min + 0.45 * width, sd * width)]) { Lower = lower, Upper = upper },
            },
        };
        ActiveSubspaceAnalysis a = f.ActiveSubspace(df, spec);
        double expected = Math.Min(sd * sd, 0.64 / (Math.PI * Math.PI)), uniform = 1.0 / (Math.PI * Math.PI);
        double worst = 0.0;
        for (int d = 0; d < a.Posterior.Draws.Count; d++)
        {
            double[,] c = a.Posterior.Draws[d].C;
            if (c[0, 0] > 0.0) worst = Math.Max(worst, Math.Abs(a.Poincare!.Bound[d, 0] / c[0, 0] - expected) / expected);
            if (c[1, 1] > 0.0) worst = Math.Max(worst, Math.Abs(a.Poincare!.Bound[d, 1] / c[1, 1] - uniform) / uniform);
        }
        output.WriteLine($"σ = {sd} of the box: constant {expected:E6} ({(sd * sd < 0.64 / (Math.PI * Math.PI) ? "σ²" : "(b − a)²/π²")}), " +
                         $"worst relative gap {worst:E2}; per-input violations {a.Poincare!.Violations}; matrix ratios up to {a.MatrixPoincare!.Ratios.Max():F4}");
        Assert.True(worst <= 1e-15, $"{worst:E2}");
        Assert.Equal(0, a.Poincare!.Violations);
        Assert.Equal(MatrixPoincareStatus.Computed, a.MatrixPoincareStatus);
        Assert.Equal(0, a.MatrixPoincare!.Violations);
    }

    /// <summary>
    /// The Gaussian equality through the moment service: a hinge (x − 0.05)₊, linear on [0.1, 0.9], under N(0.5, 0.05²) truncated
    /// there, eight sds each side, has Var/(σ²·E[φ′²]) = Var/σ² within 10⁻¹² of one and at most one; a constant 5% too small fails.
    /// </summary>
    [Fact]
    public void TruncatedNormal_LinearFactor_MeetsTheBound()
    {
        double[] x0 = [.. Enumerable.Range(0, 21).Select(i => i / 20.0)];
        double[] x1 = [.. Enumerable.Range(0, 21).Select(i => (i * 7 % 21) / 20.0)];
        var data = new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse([.. x0.Select((v, i) => v + x1[i])]));
        var transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        const double sd = 0.05;
        var measure = ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure>
        {
            ["x0"] = new MixtureMeasure([new NormalComponent(1.0, 0.5, sd)]) { Lower = 0.1, Upper = 0.9 },
        });
        ReadOnlySpan<double> knots = inputs.ContinuousKnots(0);
        int k = -1;
        for (int i = 0; i < knots.Length; i++)
            if (Math.Abs(knots[i] - 0.05) < 1e-12) k = i;
        Assert.True(k >= 0, "x0's knots hold 0.05");
        var moments = BassMoments.On(inputs, 1, measure);
        BassTerm[] term = [new BassTerm([new HingeFactor(0, 1, k)])];
        double variance = moments.CovarianceGram(term)[0, 0];
        double Ratio(double constant) => variance / moments.GradientGram(term, new double[,] { { constant, 0.0 }, { 0.0, 0.0 } })[0, 0];
        double exact = Ratio(sd * sd), short5 = Ratio(0.95 * sd * sd);
        output.WriteLine($"Var/(σ²E[φ′²]) = 1 − {1.0 - exact:E2}; with 0.95σ²: {short5:F6}");
        Assert.InRange(exact, 1.0 - 1e-12, 1.0 + 1e-15);
        Assert.True(short5 > 1.0, $"{short5:F6}");
    }

    /// <summary>
    /// The two cited bounds, checked numerically: on 300 seeded truncated normals (sds from 0.02 to 3, means within three sds of
    /// intervals inside [−0.5, 1.5]) and random smooth functions, Var(f) ≤ min(σ², (b − a)²/π²)·E[f′²] by 400-node Gauss–Legendre quadrature; and the
    /// bound is sharp at both ends: a linear f under a wide truncation, and cos(π(x − a)/(b − a)) under a nearly flat normal, each within
    /// 10⁻³ of equality.
    /// </summary>
    [Fact]
    public void TruncatedNormal_CitedBounds_SeededSweep()
    {
        var rng = new Xoshiro256PlusPlus(20261112);
        double worst = 0.0;
        for (int i = 0; i < 300; i++)
        {
            double a = -0.5 + rng.NextDouble(), b = Math.Min(1.5, a + 0.05 + 1.5 * rng.NextDouble());
            double sigma = 0.02 * Math.Pow(150.0, rng.NextDouble()), mu = a - 3.0 * sigma + (b - a + 6.0 * sigma) * rng.NextDouble();
            double[] amp = [.. Enumerable.Range(0, 3).Select(_ => Normal.SampleStandard(rng))];
            double[] freq = [.. Enumerable.Range(0, 3).Select(_ => 20.0 * rng.NextDouble())];
            double[] phase = [.. Enumerable.Range(0, 3).Select(_ => 6.3 * rng.NextDouble())];
            double cubic = Normal.SampleStandard(rng);
            double F(double x) => amp.Select((c, j) => c * Math.Sin(freq[j] * x + phase[j])).Sum() + cubic * x * x * x;
            double Slope(double x) => amp.Select((c, j) => c * freq[j] * Math.Cos(freq[j] * x + phase[j])).Sum() + 3.0 * cubic * x * x;
            double ratio = Rayleigh(mu, sigma, a, b, F, Slope);
            if (double.IsNaN(ratio)) continue;
            worst = Math.Max(worst, ratio);
        }
        double linear = Rayleigh(0.5, 0.05, 0.0, 1.0, x => x, _ => 1.0);
        double flat = Rayleigh(0.5, 200.0, 0.2, 0.7, x => Math.Cos(Math.PI * (x - 0.2) / 0.5), x => -Math.PI / 0.5 * Math.Sin(Math.PI * (x - 0.2) / 0.5));
        output.WriteLine($"worst Var/(min(σ², L²/π²)·E[f′²]) over the sweep {worst:F6}; at the ends: linear {linear:F8}, flat {flat:F8}");
        Assert.True(worst <= 1.0 + 1e-10, $"{worst:F12}");
        Assert.InRange(linear, 1.0 - 1e-3, 1.0 + 1e-10);
        Assert.InRange(flat, 1.0 - 1e-3, 1.0 + 1e-10);
    }

    // Var(f)/(min(σ², L²/π²)·E[f′²]) under N(μ, σ²) truncated to [a, b], by 400-node Gauss–Legendre on 8 panels.
    private static double Rayleigh(double mu, double sigma, double a, double b, Func<double, double> f, Func<double, double> slope)
    {
        var law = Truncated.Of(Normal.MeanSd(mu, sigma), a, b);
        GaussLegendre rule = GaussLegendre.Of(50);
        double mass = 0.0, m1 = 0.0, m2 = 0.0, energy = 0.0;
        const int panels = 8;
        double step = (b - a) / panels;
        for (int q = 0; q < panels; q++)
        {
            double mid = a + (q + 0.5) * step;
            for (int r = 0; r < rule.Count; r++)
            {
                double x = mid + 0.5 * step * rule.Nodes[r], w = 0.5 * step * rule.Weights[r] * Math.Exp(law.LogPdf(x));
                double v = f(x), s = slope(x);
                mass += w;
                m1 += w * v;
                m2 += w * v * v;
                energy += w * s * s;
            }
        }
        double mean = m1 / mass, variance = m2 / mass - mean * mean;
        double constant = Math.Min(sigma * sigma, (b - a) * (b - a) / (Math.PI * Math.PI));
        return energy > 0.0 && variance > 1e-12 * m2 / mass ? variance / (constant * energy / mass) : double.NaN;
    }

    /// <summary>Laws with no constant here leave both checks unread for that input: a mixture of two normals, a Student t component, the
    /// empirical marginal. The matrix check, asked, is not applicable, and the per-input ratio is NaN.</summary>
    [Fact]
    public void LawsWithoutAConstant_AreNotApplicable()
    {
        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        string x0 = f.Transform.Predictors[0].Name;
        InputMeasure[] laws =
        [
            new MixtureMeasure([new NormalComponent(0.5, 0.3, 0.1), new NormalComponent(0.5, 0.7, 0.1)]),
            new MixtureMeasure([new StudentComponent(1.0, 4.0, 0.5, 0.2)]),
            new EmpiricalMeasure(),
        ];
        foreach (InputMeasure law in laws)
        {
            ActiveSubspaceAnalysis a = f.ActiveSubspace(df, new ActiveSubspaceSpec { MatrixPoincare = true, Measures = new Dictionary<string, InputMeasure> { [x0] = law } });
            Assert.Equal(MatrixPoincareStatus.NotApplicable, a.MatrixPoincareStatus);
            Assert.Null(a.MatrixPoincare);
            for (int d = 0; d < a.Posterior.Draws.Count; d++) Assert.True(double.IsNaN(a.Poincare!.Ratio[d, 0]));
        }
    }

    // ── Activity scores ──────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Six draws with C_s = Q diag(λ_s) Qᵀ on a shared Q: the diagonal is C's, summing to tr C; the truncation at k is
    /// Σ_{j≤k} λ_j Q_ij², the diagonal at k = p; E[C]'s own truncation reads the weighted spectrum on Q; and under a trace-share rule each
    /// draw's truncation sits at its own k_s.
    /// </summary>
    [Fact]
    public void ActivityScores_OnKnownSpectra()
    {
        const int p = 4, n = 6;
        double[,] q = RandomOrthogonal(p, 20261113);
        var rng = new Xoshiro256PlusPlus(9);
        var lambdas = new double[n][];
        var cs = new double[n][,];
        for (int s = 0; s < n; s++)
        {
            double top = 2.0 + rng.NextDouble();
            lambdas[s] = s % 2 == 0 ? [top, 0.05, 0.02, 0.01] : [top, 1.5, 0.3, 0.01];
            cs[s] = Spectral(q, lambdas[s]);
        }
        double[] weights = [.. Enumerable.Range(0, n).Select(s => 1.0 + s)];
        for (int k = 1; k <= p; k++)
        {
            SubspacePosterior posterior = SubspacePosterior.Of(cs, weights, k);
            var rule = RankRule.TraceShare(0.95);
            ActivityScores a = ActivityScores.Of(posterior, RankPosterior.Of(posterior.Draws, weights, rule));
            double total = weights.Sum(), worst = 0.0;
            var meanSpectrum = new double[p];
            for (int s = 0; s < n; s++)
                for (int j = 0; j < p; j++) meanSpectrum[j] += weights[s] / total * lambdas[s][j];
            for (int s = 0; s < n; s++)
            {
                int ks = rule.Apply(lambdas[s]);
                double trace = 0.0;
                for (int i = 0; i < p; i++)
                {
                    Assert.Equal(cs[s][i, i], a.Draws[s, i]);
                    trace += a.Draws[s, i];
                    worst = Math.Max(worst, Math.Abs(a.Truncated[s, i] - Partial(lambdas[s], q, i, k)));
                    worst = Math.Max(worst, Math.Abs(a.RankTruncated![s, i] - Partial(lambdas[s], q, i, ks)));
                    if (k == p) worst = Math.Max(worst, Math.Abs(a.Truncated[s, i] - a.Draws[s, i]));
                }
                worst = Math.Max(worst, Math.Abs(trace - lambdas[s].Sum()));
            }
            for (int i = 0; i < p; i++)
            {
                Assert.Equal(posterior.MeanC[i, i], a.Mean[i]);
                worst = Math.Max(worst, Math.Abs(a.MeanCTruncated[i] - Partial(meanSpectrum, q, i, k)));
            }
            output.WriteLine($"k = {k}: worst gap {worst:E2}");
            Assert.True(worst <= 1e-14, $"k = {k}: {worst:E2}");
        }
    }

    /// <summary>Through the readout: the activity scores' mean is E[C]'s diagonal bit for bit, on both routes; without a rank rule
    /// there is no truncation at k_s, and with one there is.</summary>
    [Fact]
    public void ActivityScores_ThroughTheReadout()
    {
        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        foreach (ActiveSubspaceRoute route in new[] { ActiveSubspaceRoute.Exact, ActiveSubspaceRoute.Sampled })
        {
            ActiveSubspaceAnalysis a = f.ActiveSubspace(df, new ActiveSubspaceSpec { Route = route, Points = 512 });
            for (int i = 0; i < a.Inputs.Count; i++) Assert.Equal(a.Posterior.MeanC[i, i], a.Activity.Mean[i]);
            Assert.Null(a.Activity.RankTruncated);
            ActiveSubspaceAnalysis ranked = f.ActiveSubspace(df, new ActiveSubspaceSpec { Route = route, Points = 512, Rank = RankRule.TraceShare(0.9) });
            Assert.NotNull(ranked.Activity.RankTruncated);
            Assert.Equal(a.Inputs.Count, a.Activity.Mean.Length);
        }
    }

    // ── Helpers ──────────────────────────────────────────────────────────────────────────────────────────────────────

    private static double Partial(double[] lambda, double[,] q, int i, int k)
    {
        double s = 0.0;
        for (int j = 0; j < k; j++) s += lambda[j] * q[i, j] * q[i, j];
        return s;
    }

    // An orthogonal p × p matrix from a seeded normal matrix by modified Gram–Schmidt.
    internal static double[,] RandomOrthogonal(int p, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var q = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) q[i, j] = Normal.SampleStandard(rng);
        for (int j = 0; j < p; j++)
        {
            for (int k = 0; k < j; k++)
            {
                double dot = 0.0;
                for (int i = 0; i < p; i++) dot += q[i, j] * q[i, k];
                for (int i = 0; i < p; i++) q[i, j] -= dot * q[i, k];
            }
            double norm = 0.0;
            for (int i = 0; i < p; i++) norm += q[i, j] * q[i, j];
            norm = Math.Sqrt(norm);
            for (int i = 0; i < p; i++) q[i, j] /= norm;
        }
        return q;
    }

    // Q diag(λ) Qᵀ.
    internal static double[,] Spectral(double[,] q, double[] lambda)
    {
        int p = q.GetLength(0);
        var m = new double[p, p];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
            {
                double s = 0.0;
                for (int k = 0; k < lambda.Length; k++) s += q[i, k] * lambda[k] * q[j, k];
                m[i, j] = s;
            }
        return m;
    }

    private static double[,] Negated(double[,] m)
    {
        var n = (double[,])m.Clone();
        for (int i = 0; i < n.GetLength(0); i++)
            for (int j = 0; j < n.GetLength(1); j++) n[i, j] = -m[i, j];
        return n;
    }

    private static double[,] Columns(double[,] m, int k)
    {
        int p = m.GetLength(0);
        var c = new double[p, k];
        for (int i = 0; i < p; i++)
            for (int j = 0; j < k; j++) c[i, j] = m[i, j];
        return c;
    }

    private static double MaxGap(double[,] a, double[,] b)
    {
        double worst = 0.0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) worst = Math.Max(worst, Math.Abs(a[i, j] - b[i, j]));
        return worst;
    }

    private static void AssertSame(double[,] a, double[,] b)
    {
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) Assert.Equal(BitConverter.DoubleToInt64Bits(a[i, j]), BitConverter.DoubleToInt64Bits(b[i, j]));
    }

    internal static int[] Spread(int count, int k) => [.. Enumerable.Range(0, k).Select(i => (int)((long)i * (count - 1) / (k - 1)))];
}
