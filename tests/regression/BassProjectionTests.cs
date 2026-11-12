using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 7a gates (rotate-and-refit's frame and projected data, <see cref="BassProjection"/>):
/// <list type="bullet">
/// <item>the bounding box is the zonotope's: every vertex of the input box projects inside, and the extreme vertices reach each end,
/// for frames read in the box metric, a diagonal metric and a full symmetric one, where A = L⁻ᵀW₁ is checked;</item>
/// <item>the input-scale contract: a fit whose boxes differ from the analysis's, one over other levels, other inputs, and an analysis
/// without its boxes are refused;</item>
/// <item>the pushforward: under a product law its mean and covariance against Aᵀm and AᵀΣA within Monte Carlo error, the categorical
/// levels at their probabilities; under a joint sample, exactly;</item>
/// <item>projected data and cases agree bit for bit, and what is refused (d, the names, the pushforward's points);</item>
/// <item>decision 4's gate: in high p a dense direction's bounding box is several times the data's spread, the default g-form and count
/// support read that scaling not at all (draw for draw the same fit), and the column-mass support and the raw ridge read it;</item>
/// <item>slice 7d's gate: with hinges normalized over the training range, the column-mass support and the raw ridge read it not at
/// all either.</item>
/// </list>
/// </summary>
public sealed class BassProjectionTests(ITestOutputHelper output)
{
    // ── Fixtures ─────────────────────────────────────────────────────────────────────────────────────────────────────

    // p continuous inputs on [0, 1], optionally a categorical z. The box is the training range, a row of zeros and a row of ones
    // making it [0, 1]; or, declared, [0, 1] as every input's domain with the cases uniform inside, so no case sits on a corner.
    private static BassData Synthetic(int p, int n, int seed, string[]? levels = null, InputDomain? firstDomain = null, bool declared = false)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var predictors = new List<BassPredictor>();
        var y = new double[n];
        for (int i = 0; i < p; i++)
        {
            double[] x = [.. Enumerable.Range(0, n).Select(r => !declared && r == 0 ? 0.0 : !declared && r == 1 ? 1.0 : rng.NextDouble())];
            InputDomain? domain = declared ? new InputDomain(0.0, 1.0) : i == 0 ? firstDomain : null;
            predictors.Add(new ContinuousPredictor($"x{i}", x, domain));
            for (int r = 0; r < n; r++) y[r] += x[r];
        }
        if (levels is not null) predictors.Add(new CategoricalPredictor("z", [.. Enumerable.Range(0, n).Select(r => r % levels.Length)], levels));
        return new BassData(predictors, new ScalarResponse(y));
    }

    // An analysis of one draw whose C is given, read in the box metric unless a map is, under the uniform law on the boxes.
    private static ActiveSubspaceAnalysis Analysis(BassData data, double[,] c, int k, double[,]? map = null, InputMetric metric = InputMetric.Box)
    {
        var transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        var measure = ProductMeasure.Resolve(transform, inputs);
        RangeScale[] boxes = [.. transform.Predictors.OfType<RangeScale>()];
        int p = boxes.Length;
        map ??= Identity(p);
        return new ActiveSubspaceAnalysis(
            [.. boxes.Select(b => b.Name)], [.. transform.Predictors.OfType<TrainingLevels>().Select(t => t.Name)], ActiveSubspaceRoute.Exact, metric, map,
            measure, 0, SubspacePosterior.Of([c], [1.0], k), null, null, [0], Boxes: [.. boxes.Select(b => new InputDomain(b.Min, b.Max))]);
    }

    private static double[,] Identity(int p)
    {
        var m = new double[p, p];
        for (int i = 0; i < p; i++) m[i, i] = 1.0;
        return m;
    }

    // ── The bounding box ─────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// For p = 2 to 8 and d to 3, on random spectra, in the box metric, a diagonal metric and a full symmetric positive definite one:
    /// LᵀA = W₁ to 10⁻¹³; every one of the 2^p vertices projects inside each domain; the extreme vertices reach each unpadded end to
    /// 10⁻¹⁵ of Σ|a|.
    /// </summary>
    [Fact]
    public void BoundingBox_IsTheZonotopes()
    {
        var rng = new Xoshiro256PlusPlus(20261111);
        double worstMap = 0.0, worstEnd = 0.0;
        int vertices = 0;
        for (int p = 2; p <= 8; p++)
        {
            BassData data = Synthetic(p, 20, 100 + p);
            var transform = BassInputTransform.Fit(data);
            double[] lambda = [.. Enumerable.Range(0, p).Select(i => Math.Pow(0.5, i) * (1.0 + rng.NextDouble()))];
            double[,] c = ActiveSubspaceFollowUpTests.Spectral(ActiveSubspaceFollowUpTests.RandomOrthogonal(p, 200 + p), lambda);
            var diagonal = new double[p, p];
            for (int i = 0; i < p; i++) diagonal[i, i] = 0.2 + rng.NextDouble();
            double[,] q = ActiveSubspaceFollowUpTests.RandomOrthogonal(p, 300 + p);
            double[,] full = ActiveSubspaceFollowUpTests.Spectral(q, [.. Enumerable.Range(0, p).Select(_ => 0.3 + rng.NextDouble())]);
            foreach (double[,]? map in new[] { null, diagonal, full })
                for (int d = 1; d <= Math.Min(3, p); d++)
                {
                    BassProjection projection = BassProjection.Of(transform, Analysis(data, c, d, map), d);
                    double[,] a = projection.Map, w = projection.Frame, l = map ?? Identity(p);
                    for (int i = 0; i < p; i++)
                        for (int j = 0; j < d; j++)
                        {
                            double s = 0.0;
                            for (int k = 0; k < p; k++) s += l[k, i] * a[k, j];
                            worstMap = Math.Max(worstMap, Math.Abs(s - w[i, j]));
                        }
                    for (int j = 0; j < d; j++)
                    {
                        double size = Enumerable.Range(0, p).Sum(i => Math.Abs(a[i, j]));
                        double low = double.PositiveInfinity, high = double.NegativeInfinity;
                        for (int v = 0; v < 1 << p; v++, vertices++)
                        {
                            double u = 0.0;
                            for (int i = 0; i < p; i++) u += a[i, j] * ((v >> i) & 1);
                            Assert.True(projection.Domains[j].Contains(u), $"p = {p}, d = {d}: a vertex projects outside the domain");
                            low = Math.Min(low, u);
                            high = Math.Max(high, u);
                        }
                        double pad = 1e-14 * size;
                        worstEnd = Math.Max(worstEnd, Math.Max(Math.Abs(low - (projection.Domains[j].Lower + pad)), Math.Abs(high - (projection.Domains[j].Upper - pad))) / size);
                    }
                }
        }
        output.WriteLine($"{vertices} vertex projections inside their domains; LᵀA − W₁ {worstMap:E2}; the extreme vertices to each end {worstEnd:E2} of Σ|a|");
        Assert.True(worstMap <= 1e-13, $"{worstMap:E2}");
        Assert.True(worstEnd <= 1e-15, $"{worstEnd:E2}");
    }

    // ── The input-scale contract ─────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The analysis's frame holds only in its boxes. Refused: the review's case, a fit that keeps the names but declares x0's domain
    /// [−2, 2] where the analysis read [0, 1]; an analysis without its boxes; other inputs or another order; a categorical over other
    /// levels. A fit's own analysis is taken.
    /// </summary>
    [Fact]
    public void InputScaleContract_IsEnforced()
    {
        BassData data = Synthetic(3, 30, 7, ["a", "b", "c"]);
        var transform = BassInputTransform.Fit(data);
        ActiveSubspaceAnalysis analysis = Analysis(data, ActiveSubspaceFollowUpTests.Spectral(ActiveSubspaceFollowUpTests.RandomOrthogonal(3, 9), [3.0, 1.0, 0.2]), 1);
        Assert.Equal(1, BassProjection.Of(transform, analysis, 1).Dimension);
        var wider = BassInputTransform.Fit(Synthetic(3, 30, 7, ["a", "b", "c"], new InputDomain(-2.0, 2.0)));
        ArgumentException e = Assert.Throws<ArgumentException>(() => BassProjection.Of(wider, analysis, 1));
        output.WriteLine(e.Message);
        Assert.Throws<ArgumentException>(() => BassProjection.Of(transform, analysis with { Boxes = null }, 1));
        Assert.Throws<ArgumentException>(() => BassProjection.Of(BassInputTransform.Fit(Synthetic(2, 30, 7, ["a", "b", "c"])), analysis, 1));
        Assert.Throws<ArgumentException>(() => BassProjection.Of(BassInputTransform.Fit(Synthetic(3, 30, 7, ["a", "b"])), analysis, 1));

        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        BassProjection own = BassProjection.Of(f.Transform, f.ActiveSubspace(df), 1);
        Assert.Equal(f.Transform.Predictors.OfType<RangeScale>().Count(), own.Inputs.Count);
    }

    // ── The pushforward ──────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// A scalar fit with a categorical z under x0's truncated mixture and x1's uniform law, at d = 2: 200,000 pushforward points, their
    /// mean and covariance against Aᵀm and AᵀΣA from the laws' interval moments, each within 4.5 standard errors, and z's levels at their
    /// probabilities; a weighted joint sample's pushforward carries its weights, its weighted mean and covariance Aᵀm̂ and AᵀΣ̂A to
    /// 10⁻¹³.
    /// </summary>
    [Fact]
    public void Pushforward_Moments()
    {
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = SensitivityFixture.FitScalar(1);
        var measures = new Dictionary<string, InputMeasure>
        {
            ["x0"] = new MixtureMeasure([new NormalComponent(0.4, 0.3, 0.15), new StudentComponent(0.6, 3.0, 0.7, 0.2)]) { Lower = 0.05, Upper = 0.95 },
            ["z"] = new CategoricalMeasure(new Dictionary<string, double> { ["a"] = 0.2, ["b"] = 0.5, ["c"] = 0.3 }),
        };
        ActiveSubspaceAnalysis analysis = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Measures = measures, Dimension = 2 });
        BassProjection projection = BassProjection.Of(fit.Transform, analysis, 2);
        double[,] a = projection.Map;
        int p = projection.Inputs.Count;
        var m = new double[p];
        var v = new double[p];
        Span<double> moments = stackalloc double[3];
        for (int i = 0; i < p; i++)
        {
            ContinuousMeasure law = analysis.Measure!.Continuous(analysis.Measure.OfContinuous(i))!;
            double centre = 0.5 * (law.Lower + law.Upper);
            law.Measure.IntervalMoments(Math.BitDecrement(law.Lower), law.Upper, centre, 0.0, moments);
            double mean = moments[1] / moments[0];
            m[i] = centre + mean;
            v[i] = moments[2] / moments[0] - mean * mean;
        }
        const int n = 200_000;
        JointSample sample = projection.Pushforward(n, 17);
        double[][] u = [.. projection.Names.Select(name => ((ContinuousColumn)sample.Points.Columns.First(col => col.Name == name)).Values.ToArray())];
        double worst = 0.0;
        for (int j = 0; j < 2; j++)
        {
            double expected = Enumerable.Range(0, p).Sum(i => a[i, j] * m[i]), mean = u[j].Average();
            double sd = Math.Sqrt(u[j].Sum(x => (x - mean) * (x - mean)) / (n - 1));
            worst = Math.Max(worst, Math.Abs(mean - expected) / (sd / Math.Sqrt(n)));
            for (int k = j; k < 2; k++)
            {
                double expectedCov = Enumerable.Range(0, p).Sum(i => a[i, j] * v[i] * a[i, k]);
                double meanK = u[k].Average();
                double[] products = [.. Enumerable.Range(0, n).Select(r => (u[j][r] - mean) * (u[k][r] - meanK))];
                double cov = products.Sum() / (n - 1), pm = products.Average();
                double se = Math.Sqrt(products.Sum(x => (x - pm) * (x - pm)) / (n - 1) / n);
                worst = Math.Max(worst, Math.Abs(cov - expectedCov) / se);
            }
        }
        string[] levels = ((CategoricalColumn)sample.Points.Columns.First(col => col.Name == "z")).Levels.ToArray();
        foreach ((string level, double probability) in new[] { ("a", 0.2), ("b", 0.5), ("c", 0.3) })
        {
            double share = levels.Count(l => l == level) / (double)n;
            worst = Math.Max(worst, Math.Abs(share - probability) / Math.Sqrt(probability * (1 - probability) / n));
        }
        output.WriteLine($"product law, {n} points: worst |moment − expected| over its standard error {worst:F2}");
        Assert.True(worst <= 4.5, $"{worst:F2}");

        // A weighted joint sample: exact algebra.
        var rng = new Xoshiro256PlusPlus(23);
        RangeScale[] boxes = [.. fit.Transform.Predictors.OfType<RangeScale>()];
        const int k2 = 500;
        double[][] x = [.. boxes.Select(b => Enumerable.Range(0, k2).Select(_ => b.Min + rng.NextDouble() * (b.Max - b.Min)).ToArray())];
        double[] weights = [.. Enumerable.Range(0, k2).Select(_ => 0.1 + rng.NextDouble())];
        var points = new PredictionCases([new ContinuousColumn("x0", x[0]), new ContinuousColumn("x1", x[1]),
                                          new CategoricalColumn("z", [.. Enumerable.Range(0, k2).Select(r => r % 2 == 0 ? "a" : "c")])]);
        JointSample pushed = projection.Pushforward(new JointSample(points, weights));
        Assert.Equal(weights, pushed.Weights);
        double total = weights.Sum();
        double[][] box = [.. Enumerable.Range(0, p).Select(i => x[i].Select(t => boxes[i].Scale(t)).ToArray())];
        double[] mx = [.. Enumerable.Range(0, p).Select(i => Enumerable.Range(0, k2).Sum(r => weights[r] * box[i][r]) / total)];
        double[][] pu = [.. projection.Names.Select(name => ((ContinuousColumn)pushed.Points.Columns.First(col => col.Name == name)).Values.ToArray())];
        double exact = 0.0;
        for (int j = 0; j < 2; j++)
        {
            double mu = Enumerable.Range(0, k2).Sum(r => weights[r] * pu[j][r]) / total;
            exact = Math.Max(exact, Math.Abs(mu - Enumerable.Range(0, p).Sum(i => a[i, j] * mx[i])));
            for (int k = 0; k < 2; k++)
            {
                double muK = Enumerable.Range(0, k2).Sum(r => weights[r] * pu[k][r]) / total;
                double cov = Enumerable.Range(0, k2).Sum(r => weights[r] * (pu[j][r] - mu) * (pu[k][r] - muK)) / total;
                double expected = 0.0;
                for (int i = 0; i < p; i++)
                    for (int l = 0; l < p; l++)
                        expected += a[i, j] * a[l, k] * Enumerable.Range(0, k2).Sum(r => weights[r] * (box[i][r] - mx[i]) * (box[l][r] - mx[l])) / total;
                exact = Math.Max(exact, Math.Abs(cov - expected));
            }
        }
        output.WriteLine($"joint sample of {k2} weighted points: the pushforward's moments against Aᵀm̂ and AᵀΣ̂A {exact:E2}");
        Assert.True(exact <= 1e-13, $"{exact:E2}");
    }

    // ── Projected data and refusals ──────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The training data projected and its cases projected agree bit for bit; the active variables carry the bounding boxes as their
    /// domains, and z, the response and the case weights ride along. Refused: d outside 1 to p; names that are not d distinct new ones;
    /// a pushforward of fewer than two points, or of a joint-sample analysis's law; cases without a projected input, or with another
    /// continuous column.
    /// </summary>
    [Fact]
    public void ProjectedData_AndRefusals()
    {
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = SensitivityFixture.FitScalar(1);
        ActiveSubspaceAnalysis analysis = fit.ActiveSubspace(data);
        BassProjection projection = BassProjection.Of(fit.Transform, analysis, 2);
        BassData projected = projection.Project(data);
        var cases = new PredictionCases([.. data.Predictors.Select<BassPredictor, PredictionColumn>(pr => pr switch
        {
            ContinuousPredictor c => new ContinuousColumn(c.Name, c.Values.ToArray()),
            CategoricalPredictor z => new CategoricalColumn(z.Name, [.. z.Codes.ToArray().Select(code => z.Levels[code])]),
            _ => throw new InvalidOperationException(),
        })]);
        double[,] u = projection.Project(cases);
        for (int j = 0; j < 2; j++)
        {
            var column = (ContinuousPredictor)projected.Predictors[j];
            Assert.Equal(projection.Names[j], column.Name);
            Assert.Equal(projection.Domains[j], column.Domain);
            double[] values = column.Values.ToArray();
            for (int r = 0; r < values.Length; r++) Assert.Equal(BitConverter.DoubleToInt64Bits(u[r, j]), BitConverter.DoubleToInt64Bits(values[r]));
        }
        Assert.Same(data.Response, projected.Response);
        Assert.Contains(projected.Predictors, pr => pr is CategoricalPredictor { Name: "z" });

        int p = projection.Inputs.Count;
        Assert.Throws<ArgumentOutOfRangeException>(() => BassProjection.Of(fit.Transform, analysis, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => BassProjection.Of(fit.Transform, analysis, p + 1));
        Assert.Throws<ArgumentException>(() => BassProjection.Of(fit.Transform, analysis, 2, ["u"]));
        Assert.Throws<ArgumentException>(() => BassProjection.Of(fit.Transform, analysis, 2, ["u", "u"]));
        Assert.Throws<ArgumentException>(() => BassProjection.Of(fit.Transform, analysis, 2, ["u", "z"]));
        Assert.Throws<ArgumentOutOfRangeException>(() => projection.Pushforward(1, 3));
        Assert.Throws<InvalidOperationException>(() => BassProjection.Of(fit.Transform, analysis with { Measure = null }, 1).Pushforward(100, 3));
        Assert.Throws<ArgumentException>(() => projection.Project(new PredictionCases([new ContinuousColumn("x0", [0.5])])));
        Assert.Throws<ArgumentException>(() => projection.Project(new PredictionCases([new ContinuousColumn("x0", [0.5]), new ContinuousColumn("x1", [0.5]),
                                                                                     new ContinuousColumn("w", [0.5])])));
    }

    // ── Decision 4's gate ────────────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// A dense direction w = 1/√p, n = 300 uniform cases: the bounding box's width over the projected data's range grows with p (about
    /// 1, 2 and 4 at p = 2, 10 and 40). At p = 40 the projected data, y = sin(3u) plus noise of sd 0.05, are fitted twice, with the
    /// bounding box as u's domain and with its training range: under the default g-form and count support the two fits are draw for
    /// draw the same structures and their posterior means agree to 10⁻⁹, so nothing in them reads the box's scale. The column-mass
    /// support and the raw ridge read it: their two fits part.
    /// </summary>
    [Fact]
    public void HighDimensionalDomain_HowTheFitReadsItsScale()
    {
        var spreads = new List<(int P, double Ratio)>();
        BassData? data40 = null;
        BassProjection? projection40 = null;
        foreach (int p in new[] { 2, 10, 40 })
        {
            (BassData data, BassProjection projection) = DenseDirection(p);
            spreads.Add((p, projection.Spreads(data)[0]));
            if (p == 40) (data40, projection40) = (data, projection);
        }
        output.WriteLine("bounding box width over the projected range: " + string.Join(", ", spreads.Select(s => $"p = {s.P}: {s.Ratio:F2}")));
        Assert.True(spreads[0].Ratio < spreads[1].Ratio && spreads[1].Ratio < spreads[2].Ratio && spreads[2].Ratio > 3.0, "the ratio grows with p");

        // y = sin(3u) + noise on the projected inputs, with the bounding box and with the training range.
        BassData box = projection40!.Project(data40!);
        double[] u = ((ContinuousPredictor)box.Predictors[0]).Values.ToArray();
        var rng = new Xoshiro256PlusPlus(41);
        double[] y = [.. u.Select(t => Math.Sin(3.0 * t) + 0.05 * Normal.SampleStandard(rng))];
        var withBox = new BassData([new ContinuousPredictor("u1", u, projection40.Domains[0])], new ScalarResponse(y));
        var withRange = new BassData([new ContinuousPredictor("u1", u)], new ScalarResponse(y));
        var cases = new PredictionCases([new ContinuousColumn("u1", u)]);
        BassConfig defaults = SensitivityFixture.Config(maxBasis: 8);
        foreach ((string name, BassConfig config, bool reads) in new[]
        {
            ("g-form, count support", defaults, false),
            ("column-mass support", defaults with { Structure = defaults.Structure with { Support = new SupportRule(SupportKind.ColumnMass) } }, true),
            ("raw ridge", defaults with { Coefficients = BassCoefficientPrior.RawRidge(TauPrior.UnitInformationGamma()) }, true),
        })
        {
            BassFit a = Bass.Fit(config, withBox), b = Bass.Fit(config, withRange);
            BassDraw[] da = [.. a.Draws.SelectMany(ch => ch)], db = [.. b.Draws.SelectMany(ch => ch)];
            int differ = 0;
            for (int s = 0; s < da.Length; s++)
            {
                BassTerm[] ta = [.. da[s].Structure.TermIds.Select(id => a.Terms[id])], tb = [.. db[s].Structure.TermIds.Select(id => b.Terms[id])];
                if (!ta.SequenceEqual(tb)) differ++;
            }
            double[,] ma = a.Predict(withBox, cases).Mean(), mb = b.Predict(withRange, cases).Mean();
            double scale = 0.0, gap = 0.0;
            for (int r = 0; r < u.Length; r++)
            {
                scale = Math.Max(scale, Math.Abs(ma[r, 0]));
                gap = Math.Max(gap, Math.Abs(ma[r, 0] - mb[r, 0]));
            }
            output.WriteLine($"{name}: {differ} of {da.Length} draws with other structures; posterior means apart by {gap / scale:E2} of their size");
            if (reads) Assert.True(differ > 0 || gap / scale > 1e-6, $"{name} reads the box's scale, yet the fits agree");
            else
            {
                Assert.Equal(0, differ);
                Assert.True(gap / scale <= 1e-9, $"{name}: {gap / scale:E2}");
            }
        }
    }

    // Decision 4's fixture: p inputs declared on [0, 1] with 300 uniform cases, the frame a dense direction w = 1/√p, d = 1.
    private static (BassData Data, BassProjection Projection) DenseDirection(int p)
    {
        BassData data = Synthetic(p, 300, 1000 + p, declared: true);
        var transform = BassInputTransform.Fit(data);
        var w = new double[p, 1];
        for (int i = 0; i < p; i++) w[i, 0] = 1.0 / Math.Sqrt(p);
        double[,] c = ActiveSubspaceFollowUpTests.Spectral(Complete(w), [.. Enumerable.Range(0, p).Select(i => i == 0 ? 1.0 : 1e-3 / (i + 1))]);
        return (data, BassProjection.Of(transform, Analysis(data, c, 1), 1));
    }

    /// <summary>
    /// Slice 7d's gate on decision 4's p = 40 data. With hinges normalized over the training range, the fit with the bounding box as
    /// u's domain is the fit with the training range as its domain, under the column-mass support and the raw ridge as under the
    /// g-form with count support: draw for draw the same structures, posterior means to 10⁻⁹ of their size. On the training-range
    /// data the training range is the box, so the two rules there are the same bits.
    /// </summary>
    [Fact]
    public void HighDimensionalDomain_TrainingRangeHinges_ReadTheDataScale()
    {
        (BassData data, BassProjection projection) = DenseDirection(40);
        BassData box = projection.Project(data);
        double[] u = ((ContinuousPredictor)box.Predictors[0]).Values.ToArray();
        var rng = new Xoshiro256PlusPlus(41);
        double[] y = [.. u.Select(t => Math.Sin(3.0 * t) + 0.05 * Normal.SampleStandard(rng))];
        var withBox = new BassData([new ContinuousPredictor("u1", u, projection.Domains[0])], new ScalarResponse(y));
        var withRange = new BassData([new ContinuousPredictor("u1", u)], new ScalarResponse(y));
        var cases = new PredictionCases([new ContinuousColumn("u1", u)]);
        BassConfig defaults = SensitivityFixture.Config(maxBasis: 8);
        static BassConfig Trained(BassConfig c) => c with { Structure = c.Structure with { HingeNormalization = HingeNormalization.TrainingRange } };
        foreach ((string name, BassConfig config) in new[]
        {
            ("g-form, count support", defaults),
            ("column-mass support", defaults with { Structure = defaults.Structure with { Support = new SupportRule(SupportKind.ColumnMass) } }),
            ("raw ridge", defaults with { Coefficients = BassCoefficientPrior.RawRidge(TauPrior.UnitInformationGamma()) }),
        })
        {
            BassFit wide = Bass.Fit(Trained(config), withBox), reference = Bass.Fit(config, withRange), trained = Bass.Fit(Trained(config), withRange);
            BassDraw[] dw = [.. wide.Draws.SelectMany(ch => ch)], dr = [.. reference.Draws.SelectMany(ch => ch)], dt = [.. trained.Draws.SelectMany(ch => ch)];
            int differ = 0;
            for (int s = 0; s < dw.Length; s++)
            {
                BassTerm[] tw = [.. dw[s].Structure.TermIds.Select(id => wide.Terms[id])], tr = [.. dr[s].Structure.TermIds.Select(id => reference.Terms[id])];
                if (!tw.SequenceEqual(tr)) differ++;
            }
            double[,] mw = wide.Predict(withBox, cases).Mean(), mr = reference.Predict(withRange, cases).Mean(), mt = trained.Predict(withRange, cases).Mean();
            double scale = 0.0, gap = 0.0;
            for (int r = 0; r < u.Length; r++)
            {
                scale = Math.Max(scale, Math.Abs(mr[r, 0]));
                gap = Math.Max(gap, Math.Abs(mw[r, 0] - mr[r, 0]));
            }
            output.WriteLine($"{name}, training-range hinges: {differ} of {dw.Length} draws with other structures than the training-range domain's; " +
                             $"posterior means apart by {gap / scale:E2} of their size");
            Assert.Equal(0, differ);
            Assert.True(gap / scale <= 1e-9, $"{name}: {gap / scale:E2}");
            Assert.Equal(dr.Length, dt.Length);
            for (int s = 0; s < dr.Length; s++) Assert.Equal(dr[s].Structure.TermIds.Select(id => reference.Terms[id]), dt[s].Structure.TermIds.Select(id => trained.Terms[id]));
            Assert.Equal(mr, mt);
        }
    }

    // A p × p orthonormal basis whose first column is w's direction, by Gram–Schmidt on w then the axes.
    private static double[,] Complete(double[,] w)
    {
        int p = w.GetLength(0);
        var q = new double[p, p];
        var columns = new List<double[]>();
        IEnumerable<double[]> candidates = new[] { Enumerable.Range(0, p).Select(i => w[i, 0]).ToArray() }
            .Concat(Enumerable.Range(0, p).Select(e => Enumerable.Range(0, p).Select(i => i == e ? 1.0 : 0.0).ToArray()));
        foreach (double[] v in candidates)
        {
            if (columns.Count == p) break;
            double[] r = (double[])v.Clone();
            foreach (double[] c in columns)
            {
                double dot = r.Zip(c, (x, y) => x * y).Sum();
                for (int i = 0; i < p; i++) r[i] -= dot * c[i];
            }
            double norm = Math.Sqrt(r.Sum(x => x * x));
            if (norm < 1e-8) continue;
            columns.Add([.. r.Select(x => x / norm)]);
        }
        for (int j = 0; j < p; j++)
            for (int i = 0; i < p; i++) q[i, j] = columns[j][i];
        return q;
    }
}
