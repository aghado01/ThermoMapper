using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 7d: hinges normalized over each variable's training range (<see cref="HingeNormalization.TrainingRange"/>), beside the
/// box, R's rule and the default.
/// <list type="bullet">
/// <item>Without a declared domain the training range is the box, and the two rules are the same bits: columns, derivatives, κ.</item>
/// <item>Under declared domains wider than the data, the training-range columns are the same function of the raw inputs as the box
/// rule's columns on the same data undeclared, at the training cases and beyond them to the domains' ends; their derivatives in raw
/// units too. κ itself moves with the box, and the box rule's columns read it.</item>
/// <item>The moments take the rule: means, cross moments and the gradient Gram are the box rule's times κ's ratios, and a mean
/// against quadrature of the column.</item>
/// <item>The prefix tables take it: a border from tables against the dots, under declared domains.</item>
/// <item>A hinge facing out at its training end: zero on the training cases, c = 0 left at 1, refused by the support.</item>
/// <item>The setting: a Legendre basis refuses it, null stays out of the serialized form, a stated one round-trips, two fits paired by
/// a co-active readout must share it, and a projected fit resolves null to the training range and keeps a stated one.</item>
/// </list>
/// </summary>
public sealed class HingeNormalizationTests(ITestOutputHelper output)
{
    private static readonly InputDomain[] Domains = [new(0.0, 10.0), new(-3.0, 6.0), new(1.5, 9.0)];

    // Three inputs uniform on [2, 5]: once undeclared, so each box is its training range, and once in Domains, wider on both sides.
    private static (BassData Range, BassData Wide) Data(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[][] x = [.. Enumerable.Range(0, 3).Select(_ => Enumerable.Range(0, n).Select(_ => 2.0 + 3.0 * rng.NextDouble()).ToArray())];
        double[] y = [.. Enumerable.Range(0, n).Select(r => Math.Sin(2.0 * x[0][r]) + x[1][r] * x[2][r] / 10.0 + 0.05 * Normal.SampleStandard(rng))];
        var range = new BassData([.. Enumerable.Range(0, 3).Select(i => new ContinuousPredictor($"x{i}", x[i]))], new ScalarResponse(y));
        var wide = new BassData([.. Enumerable.Range(0, 3).Select(i => new ContinuousPredictor($"x{i}", x[i], Domains[i]))], new ScalarResponse(y));
        return (range, wide);
    }

    private static (BassBasis Basis, BassInputTransform Transform) Basis(BassData data, int degree, HingeNormalization normalization)
    {
        BassInputTransform transform = BassInputTransform.Fit(data);
        return (new BassBasis(BassInputs.Prepare(data, transform), degree, normalization), transform);
    }

    // Every one-way hinge, each knot and sign, then 300 random two- and three-way products.
    private static List<BassTerm> Terms(BassInputs inputs, int seed)
    {
        var terms = new List<BassTerm>();
        for (int v = 0; v < 3; v++)
            for (int k = 0; k < inputs.ContinuousKnots(v).Length; k++)
                foreach (int sign in new[] { 1, -1 }) terms.Add(new BassTerm([new HingeFactor(v, sign, k)]));
        var rng = new Xoshiro256PlusPlus(seed);
        for (int t = 0; t < 300; t++)
        {
            int order = 2 + (int)(rng.NextUInt64() % 2);
            int[] variables = [.. Enumerable.Range(0, 3).OrderBy(_ => rng.NextDouble()).Take(order).Order()];
            terms.Add(new BassTerm([.. variables.Select(v => new HingeFactor(v, rng.NextDouble() < 0.5 ? 1 : -1, (int)(rng.NextUInt64() % (ulong)inputs.ContinuousKnots(v).Length)))]));
        }
        return terms;
    }

    // Raw cases uniform over each declared domain, so most lie beyond the training range [2, 5], in each transform's box units.
    private static (CaseInputs Range, CaseInputs Wide, int Count) Fresh(BassInputTransform range, BassInputTransform wide, int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[][] raw = [.. Domains.Select(d => Enumerable.Range(0, n).Select(_ => d.Lower + (d.Upper - d.Lower) * rng.NextDouble()).ToArray())];
        CaseInputs At(BassInputTransform t) => new([.. raw.Select((column, v) => column.Select(((RangeScale)t.Predictors[v]).Scale).ToArray())], []);
        return (At(range), At(wide), n);
    }

    private static double[] Column(BassBasis basis, BassTerm term, CaseInputs at, int count, int variable = -1)
    {
        var column = new double[count];
        if (variable < 0) basis.CaseColumn(term, at, column);
        else basis.CaseDerivative(term, variable, at, column);
        return column;
    }

    private static double Width(BassInputTransform t, int v) => ((RangeScale)t.Predictors[v]).Max - ((RangeScale)t.Predictors[v]).Min;

    /// <summary>On data without a declared domain the training range is the box: every term's column, derivatives and κ, the same bits.</summary>
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void Undeclared_TrainingRangeIsTheBox_SameBits(int degree)
    {
        (BassData range, _) = Data(150, 3);
        (BassBasis box, BassInputTransform t) = Basis(range, degree, HingeNormalization.Box);
        (BassBasis trained, _) = Basis(range, degree, HingeNormalization.TrainingRange);
        (CaseInputs fresh, _, int m) = Fresh(t, t, 200, 4);
        const int n = 150;
        foreach (BassTerm term in Terms(box.Inputs, 5))
        {
            Assert.Equal(box.Kernels.Case[0].Scale(term), trained.Kernels.Case[0].Scale(term));
            Assert.Equal(Column(box, term, box.Inputs.Cases, n), Column(trained, term, trained.Inputs.Cases, n));
            Assert.Equal(Column(box, term, fresh, m), Column(trained, term, fresh, m));
            foreach (HingeFactor f in term.Continuous)
                Assert.Equal(Column(box, term, fresh, m, f.Variable), Column(trained, term, fresh, m, f.Variable));
        }
    }

    /// <summary>
    /// Under declared domains wider than the data, at degrees 1 and 2: the training-range columns against the box rule's on the same
    /// data undeclared, at the 200 training cases and at 300 cases over the domains, most beyond the training range, to 10⁻¹² of the
    /// column's size; their derivatives in raw units to 10⁻¹¹. κ moves with the box: (W/w)^q apart for a one-way hinge, W the domain's
    /// width and w the training range's. The box rule's columns on the declared data are the undeclared ones over that ratio. Terms with
    /// a factor facing out at its training end are counted apart: peak 0 under both rules, zero on the training cases.
    /// </summary>
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void Declared_TrainingRangeColumns_AreTheDataScales(int degree)
    {
        (BassData range, BassData wide) = Data(200, 11);
        (BassBasis reference, BassInputTransform tRange) = Basis(range, degree, HingeNormalization.Box);
        (BassBasis trained, BassInputTransform tWide) = Basis(wide, degree, HingeNormalization.TrainingRange);
        (BassBasis boxed, _) = Basis(wide, degree, HingeNormalization.Box);
        (CaseInputs freshRange, CaseInputs freshWide, int m) = Fresh(tRange, tWide, 300, 12);
        double worstColumn = 0.0, worstDerivative = 0.0, largest = 0.0;
        int terms = 0, guarded = 0;
        foreach (BassTerm term in Terms(reference.Inputs, 13))
        {
            // A factor facing out at its training end has peak 0: c is left at 1, the column is zero on the data, and beyond them it is
            // in each box's own units (OutwardHingeAtTheTrainingEnd_ZeroSupport_Refused). The support never admits it.
            if (reference.Kernels.ContinuousPeak(term) == 0.0)
            {
                Assert.Equal(0.0, trained.Kernels.ContinuousPeak(term));
                Assert.All(Column(trained, term, trained.Inputs.Cases, 200), v => Assert.Equal(0.0, v));
                guarded++;
                continue;
            }
            foreach ((CaseInputs a, CaseInputs b, int count) in new[] { (trained.Inputs.Cases, reference.Inputs.Cases, 200), (freshWide, freshRange, m) })
            {
                double[] x = Column(trained, term, a, count), y = Column(reference, term, b, count);
                double size = Math.Max(1.0, y.Max(Math.Abs));
                largest = Math.Max(largest, y.Max(Math.Abs));
                worstColumn = Math.Max(worstColumn, x.Zip(y, (p, q) => Math.Abs(p - q)).Max() / size);
                foreach (HingeFactor f in term.Continuous)
                {
                    // ∂/∂x_raw = ∂/∂x_box / width, each in its own box.
                    double[] dx = Column(trained, term, a, count, f.Variable), dy = Column(reference, term, b, count, f.Variable);
                    double wa = Width(tWide, f.Variable), wb = Width(tRange, f.Variable);
                    double dsize = Math.Max(1.0, dy.Max(v => Math.Abs(v / wb)));
                    worstDerivative = Math.Max(worstDerivative, dx.Zip(dy, (p, q) => Math.Abs(p / wa - q / wb)).Max() / dsize);
                }
            }
            terms++;
        }
        output.WriteLine($"degree {degree}: {terms} terms, {guarded} more facing out at a training end; columns at the training cases and over the domains (largest {largest:F1}) apart by at most " +
                         $"{worstColumn:E2} of their size, raw derivatives by {worstDerivative:E2}");
        Assert.True(worstColumn <= 1e-12, $"{worstColumn:E2}");
        Assert.True(worstDerivative <= 1e-11, $"{worstDerivative:E2}");
        Assert.True(largest > 1.5, "the fresh cases reach beyond the training range, where a training-range column exceeds one");

        // κ moves with the box, and the box rule's columns read it: one-way hinges at an interior knot.
        for (int v = 0; v < 3; v++)
        {
            int k = reference.Inputs.ContinuousKnots(v).Length / 2;
            var term = new BassTerm([new HingeFactor(v, 1, k)]);
            double ratio = Math.Pow(Width(tWide, v) / Width(tRange, v), degree);
            double kt = trained.Kernels.Case[0].Scale(term), kr = reference.Kernels.Case[0].Scale(term);
            Assert.True(Math.Abs(kt / kr - ratio) <= 1e-12 * ratio, $"κ's ratio {kt / kr} against (W/w)^q = {ratio}");
            double[] xb = Column(boxed, term, boxed.Inputs.Cases, 200), xr = Column(reference, term, reference.Inputs.Cases, 200);
            double shrink = xr.Max() / xb.Max();
            Assert.True(shrink > 1.2, $"the box rule's column on the declared domain is not smaller: {shrink}");
        }
    }

    /// <summary>
    /// The moments under the training range, with declared domains, the measure uniform on them: each term's mean, each pair's cross
    /// moment and the gradient Gram are the box rule's times κ's ratios, to 10⁻¹²; and a one-way hinge's mean against the midpoint rule
    /// on its column, 2·10⁵ points, to 10⁻⁸.
    /// </summary>
    [Fact]
    public void Moments_TakeTheRule()
    {
        (_, BassData wide) = Data(200, 21);
        BassInputTransform transform = BassInputTransform.Fit(wide);
        BassInputs inputs = BassInputs.Prepare(wide, transform);
        ProductMeasure measure = ProductMeasure.Resolve(transform, inputs);
        BassMoments box = BassMoments.On(inputs, 1, measure), trained = BassMoments.On(inputs, 1, measure, HingeNormalization.TrainingRange);
        BassTerm[] terms = [.. Terms(inputs, 22).Where((_, i) => i % 37 == 0)];
        double[] r = [.. terms.Select(t => trained.Scale(t) / box.Scale(t))];
        double worst = 0.0;
        for (int a = 0; a < terms.Length; a++)
        {
            worst = Math.Max(worst, Rel(trained.Mean(terms[a]), box.Mean(terms[a]) * r[a]));
            for (int b = a; b < terms.Length; b += 5) worst = Math.Max(worst, Rel(trained.Cross(terms[a], terms[b]), box.Cross(terms[a], terms[b]) * r[a] * r[b]));
        }
        double[,] gt = trained.GradientGram(terms, Identity(3)), gb = box.GradientGram(terms, Identity(3));
        for (int a = 0; a < terms.Length; a++)
            for (int b = 0; b < terms.Length; b++) worst = Math.Max(worst, Rel(gt[a, b], gb[a, b] * r[a] * r[b]));
        Assert.True(worst <= 1e-12, $"{worst:E2}");

        // A one-way hinge's mean against quadrature of its training-range column over the domain.
        var basis = new BassBasis(inputs, 1, HingeNormalization.TrainingRange);
        const int points = 200_000;
        var at = new CaseInputs([[.. Enumerable.Range(0, points).Select(i => (i + 0.5) / points)], new double[points], new double[points]], []);
        foreach (BassTerm term in new[] { new BassTerm([new HingeFactor(0, 1, 60)]), new BassTerm([new HingeFactor(0, -1, 140)]) })
        {
            double quadrature = Column(basis, term, at, points).Average();
            output.WriteLine($"{term}: mean {trained.Mean(term):F10} against quadrature {quadrature:F10}");
            Assert.True(Math.Abs(trained.Mean(term) - quadrature) <= 1e-8 * Math.Max(1.0, Math.Abs(quadrature)));
        }
    }

    private static double Rel(double a, double b) => a.Equals(b) ? 0.0 : Math.Abs(a - b) / Math.Max(1e-300, Math.Max(Math.Abs(a), Math.Abs(b)));

    private static double[,] Identity(int p)
    {
        var m = new double[p, p];
        for (int i = 0; i < p; i++) m[i, i] = 1.0;
        return m;
    }

    // The rich surface's 1000 cases on eight inputs uniform on [0, 1], declared on [−1, 2.5], the training range the hinges' interval.
    private static BassModel WideModel(HingeNormalization normalization)
    {
        BassCost.Point point = BassCost.Base with { N = 1000, P = 8, Cap = 20 };
        BassData data = BassCost.Data(point);
        var wide = new BassData([.. data.Predictors.Select(p => (BassPredictor)new ContinuousPredictor(p.Name, ((ContinuousPredictor)p).Values.ToArray(), new InputDomain(-1.0, 2.5)))],
                                data.Response);
        BassConfig config = BassCost.Config(point);
        return BassModel.Create(config with { Structure = config.Structure with { HingeNormalization = normalization } }, wide);
    }

    /// <summary>
    /// The prefix tables divide by the kernel's c: under declared domains and the training range, a border from tables against the dots
    /// over 400 random one-way hinges, within lever 1's bound (BassPrefixTableTests), where the box's c would put them (W/w)^q apart.
    /// </summary>
    [Fact]
    public void PrefixTables_TakeTheRule()
    {
        BassModel model = WideModel(HingeNormalization.TrainingRange);
        BassDesign design = model.NewDesign();
        design.TablesOn = true;
        var rng = new Xoshiro256PlusPlus(seed: 5301);
        for (int step = 0; step < 600 && design.Count < 20; step++)
        {
            var counts = new BassTermCounts(model.Prior);
            foreach (BassTerm t in design.Terms) counts.Add(t);
            if (model.Proposal.DrawBirth(rng, counts) is { } born && design.ProposeBirth(born).Admitted) design.Accept();
        }
        Assert.True(design.Count >= 10);
        double worst = 0.0;
        int compared = 0;
        for (int trial = 0; trial < 400; trial++)
        {
            var one = new BassTerm([new HingeFactor((int)(rng.NextUInt64() % 8), rng.NextDouble() < 0.5 ? 1 : -1, (int)(rng.NextUInt64() % 1000))]);
            if (design.BorderComparison(one) is not var (tabled, dotted)) continue;
            for (int j = 0; j < tabled.Length; j++)
            {
                worst = Math.Max(worst, Math.Abs(tabled[j] - dotted[j]) / Math.Max(Math.Sqrt(((BorderedCholesky)design.Factor).Entry(j, j)), 1e-300));
                compared++;
            }
        }
        Assert.True(compared > 0);
        output.WriteLine($"{compared} border entries under the training range; the worst difference is {worst:E2} of √G_jj");
        Assert.True(worst <= 4000 * 1.1102230246251565e-16 * Math.Sqrt(1000));
    }

    /// <summary>
    /// A hinge facing out at its variable's training end, under declared domains: its peak over the training range is 0, so c is left
    /// at 1 (κ = 1); its column is zero at every training case, finite and positive beyond the end inside the domain; the support
    /// refuses its birth. On the box rule its c is not 0, and it is refused the same way.
    /// </summary>
    [Fact]
    public void OutwardHingeAtTheTrainingEnd_ZeroSupport_Refused()
    {
        foreach (HingeNormalization normalization in new[] { HingeNormalization.TrainingRange, HingeNormalization.Box })
        {
            BassModel model = WideModel(normalization);
            BassBasis basis = model.Basis;
            int last = basis.Inputs.ContinuousKnots(2).Length - 1;
            foreach (BassTerm term in new[] { new BassTerm([new HingeFactor(2, 1, last)]), new BassTerm([new HingeFactor(2, -1, 0)]) })
            {
                double kappa = basis.Kernels.Case[0].Scale(term);
                if (normalization == HingeNormalization.TrainingRange) Assert.Equal(1.0, kappa);
                else Assert.NotEqual(1.0, kappa);
                double[] column = Column(basis, term, basis.Inputs.Cases, 1000);
                Assert.All(column, v => Assert.Equal(0.0, v));
                double beyond = term.Continuous[0].Sign > 0 ? 0.95 : 0.05;   // inside the declared [−1, 2.5], past the data's [0, 1]
                var at = new CaseInputs([.. Enumerable.Range(0, 8).Select(v => new[] { v == 2 ? beyond : 0.5 })], []);
                double value = Column(basis, term, at, 1)[0];
                Assert.True(value > 0.0 && double.IsFinite(value), $"{term} beyond the training end: {value}");
                BassDesign design = model.NewDesign();
                Assert.False(design.ProposeBirth(term).Admitted);
            }
        }
    }

    /// <summary>
    /// The setting: a Legendre basis refuses it; an unknown value is refused; null stays out of the serialized form; a stated one
    /// round-trips. Two fits a co-active readout pairs must share it.
    /// </summary>
    [Fact]
    public void Setting_RefusedWhereItReadsNothing_AndSerialized()
    {
        BassConfig legendre = KhaosReadoutTests.Config(2);
        ArgumentException refused = Assert.Throws<ArgumentException>(() =>
            (legendre with { Structure = legendre.Structure with { HingeNormalization = HingeNormalization.Box } }).Validate());
        Assert.Contains("Structure.HingeNormalization", refused.Message);
        BassConfig hinge = SensitivityFixture.Config();
        Assert.Throws<ArgumentException>(() => (hinge with { Structure = hinge.Structure with { HingeNormalization = (HingeNormalization)7 } }).Validate());

        Assert.DoesNotContain("HingeNormalization", hinge.ToJson());
        BassConfig trained = hinge with { Structure = hinge.Structure with { HingeNormalization = HingeNormalization.TrainingRange } };
        string json = trained.ToJson();
        Assert.Contains("\"HingeNormalization\":\"TrainingRange\"", json);
        Assert.Equal(HingeNormalization.TrainingRange, BassConfig.FromJson(json).Structure.HingeNormalization);
        Assert.Equal(json, BassConfig.FromJson(json).ToJson());
        Assert.Equal(HingeNormalization.Box, hinge.Structure.ResolvedHingeNormalization);

        (_, BassData wide) = Data(120, 31);
        BassFit a = Bass.Fit(hinge, wide), b = Bass.Fit(trained, wide);
        ArgumentException paired = Assert.Throws<ArgumentException>(() => a.CoActiveSubspace(wide, b, wide));
        Assert.Contains("normalized over", paired.Message);
    }

    /// <summary>
    /// A projected fit resolves an unset normalization to the training range and records it in the refit's config; a stated box is kept.
    /// </summary>
    [Fact]
    public void ProjectedFit_DefaultsToTheTrainingRange()
    {
        BassProjectedFitTests.RidgeFits f = BassProjectedFitTests.Fixture.Value;
        Assert.Null(f.Native.Config.Structure.HingeNormalization);
        Assert.Equal(HingeNormalization.TrainingRange, f.Refit.Fit.Config.Structure.HingeNormalization);
        BassConfig config = f.Refit.Fit.Config with { Structure = f.Refit.Fit.Config.Structure with { HingeNormalization = HingeNormalization.Box } };
        BassProjectedFit stated = BassProjectedFit.Of(f.Projection, f.Data, config);
        Assert.Equal(HingeNormalization.Box, stated.Fit.Config.Structure.HingeNormalization);
    }
}
