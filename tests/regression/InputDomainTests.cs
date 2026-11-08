using System;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// A continuous input's declared domain (<see cref="InputDomain"/>, E19's sixth decision): what it refuses; the box it scales the
/// input by, with the training range kept for the extrapolation flags; the fingerprint, which takes it only when declared; the
/// sensitivity measure's default, uniform over it; the resolution's log; and a hinge fit's law, which it leaves alone, since the
/// hinges' knots are the data's and the g-form reads each column up to its scale.
/// </summary>
public sealed class InputDomainTests(ITestOutputHelper output)
{
    // BassSweepTests' small grid, each predictor declaring `domain` (null: none).
    private static BassData Declared(BassData data, InputDomain? domain)
        => new([.. data.Predictors.Cast<ContinuousPredictor>().Select(c => new ContinuousPredictor(c.Name, c.Values.ToArray(), domain))], data.Response);

    [Fact]
    public void Domain_RefusesBadEnds_AndAValueOutsideIt()
    {
        Assert.Throws<ArgumentException>(() => new InputDomain(1.0, 1.0));
        Assert.Throws<ArgumentException>(() => new InputDomain(2.0, 1.0));
        Assert.Throws<ArgumentException>(() => new InputDomain(double.NaN, 1.0));
        Assert.Throws<ArgumentException>(() => new InputDomain(0.0, double.PositiveInfinity));
        var domain = new InputDomain(-1.0, 2.0);
        Assert.True(domain.Contains(-1.0) && domain.Contains(2.0) && !domain.Contains(2.0000001));
        _ = new ContinuousPredictor("x", [-1.0, 0.5, 2.0], domain);   // the ends belong to it
        Assert.Throws<ArgumentException>(() => new ContinuousPredictor("x", [-1.0, 0.5, 2.5], domain));
        Assert.Throws<ArgumentException>(() => new RangeScale("x", -2.0, 1.0, domain));
    }

    [Fact]
    public void Scale_ReadsTheDomain_AndFlagsReadTheTrainingRange()
    {
        BassData data = Declared(BassSweepTests.Data(), new InputDomain(-1.0, 3.0));   // the inputs on {0, ½, 1}
        BassInputTransform transform = BassInputTransform.Fit(data);
        var x0 = (RangeScale)transform.Predictors[0];
        Assert.Equal((-1.0, 3.0, 0.0, 1.0), (x0.Min, x0.Max, x0.TrainingMin, x0.TrainingMax));
        Assert.Equal((0.0, 1.0, 0.25), (x0.Scale(-1.0), x0.Scale(3.0), x0.Scale(0.0)));
        Assert.True(x0.Contains(0.0) && x0.Contains(1.0) && !x0.Contains(2.0) && !x0.Contains(-0.5));

        BassInputs inputs = BassInputs.Prepare(data, transform);
        Assert.Equal(new[] { 0.25, 0.375, 0.5 }, inputs.Cases.Continuous(0).ToArray().Distinct().Order());

        // A case inside the domain but past the data is flagged: the hinges extrapolate past the last knot either way.
        var cases = new PredictionCases([new ContinuousColumn("x0", [0.5, 2.0]), new ContinuousColumn("x1", [0.5, 0.5])]);
        (CaseInputs scaled, bool[] outside) = cases.Prepare(transform, inputs, default);
        Assert.Equal(new[] { false, true }, outside);
        Assert.Equal(0.75, scaled.Continuous(0)[1]);

        // Undeclared, the box is the training range, as before.
        var plain = (RangeScale)BassInputTransform.Fit(BassSweepTests.Data()).Predictors[0];
        Assert.Equal((0.0, 1.0), (plain.Min, plain.Max));
        Assert.Null(plain.Domain);
    }

    [Fact]
    public void Fingerprint_TakesADomain_OnlyWhenDeclared()
    {
        BassData data = BassSweepTests.Data();
        string none = Bass.Fingerprint(new BassConfig(), data);
        Assert.Equal(none, Bass.Fingerprint(new BassConfig(), Declared(data, null)));
        string unit = Bass.Fingerprint(new BassConfig(), Declared(data, new InputDomain(0.0, 1.0)));
        Assert.NotEqual(none, unit);   // the same box as the training range, declared, is still a declaration
        Assert.NotEqual(unit, Bass.Fingerprint(new BassConfig(), Declared(data, new InputDomain(0.0, 2.0))));
        Assert.Equal(unit, Bass.Fingerprint(new BassConfig(), Declared(data, new InputDomain(0.0, 1.0))));
    }

    [Fact]
    public void Measure_IsUniformOverTheDomain_AndSaysItExtrapolates()
    {
        BassData data = Declared(BassSweepTests.Data(), new InputDomain(-1.0, 3.0));
        BassInputTransform transform = BassInputTransform.Fit(data);
        ProductMeasure measure = ProductMeasure.Resolve(transform, BassInputs.Prepare(data, transform));
        ContinuousMeasure x0 = measure.Continuous(0)!;
        Assert.Equal((0.0, 1.0, 0.25, 0.5), (x0.Lower, x0.Upper, x0.TrainingLower, x0.TrainingUpper));
        Assert.True(x0.Extrapolates);

        // A measure on the training range stays inside it.
        ProductMeasure onData = ProductMeasure.Resolve(transform, BassInputs.Prepare(data, transform),
                                                       new System.Collections.Generic.Dictionary<string, InputMeasure> { ["x0"] = new UniformMeasure(0.0, 1.0) });
        Assert.False(onData.Continuous(0)!.Extrapolates);

        BassData plain = BassSweepTests.Data();
        BassInputTransform unscaled = BassInputTransform.Fit(plain);
        ContinuousMeasure x = ProductMeasure.Resolve(unscaled, BassInputs.Prepare(plain, unscaled)).Continuous(0)!;
        Assert.Equal((0.0, 1.0, 0.0, 1.0, false), (x.Lower, x.Upper, x.TrainingLower, x.TrainingUpper, x.Extrapolates));
    }

    [Fact]
    public void Resolution_LogsEachDomain()
    {
        BassModel model = BassModel.Create(new BassConfig(), Declared(BassSweepTests.Data(), new InputDomain(-1.0, 3.0)));
        BassAutoPick[] domains = [.. model.Resolved.AutoPicks.Where(p => p.Setting.StartsWith("Domain[", StringComparison.Ordinal))];
        Assert.Equal(new[] { "Domain[x0]", "Domain[x1]" }, domains.Select(p => p.Setting));
        Assert.Equal("[-1, 3]", domains[0].Value);
        Assert.DoesNotContain(BassModel.Create(new BassConfig(), BassSweepTests.Data()).Resolved.AutoPicks, p => p.Setting.StartsWith("Domain[", StringComparison.Ordinal));
    }

    [Fact]
    public void HingeLaw_IsTheTrainingRangesUnderADomain()
    {
        // The enumerated law at a held τ, with and without a domain wider than the data: the knots are the data's, so the same
        // structures are admitted, and the g-form reads each column up to its scale, which the domain alone changes.
        BassConfig config = BassSweepTests.Config() with { Coefficients = BassCoefficientPrior.G(TauPrior.Held()) };
        BassKernelEnumeration.Space plain = BassKernelEnumeration.Enumerate(BassModel.Create(config, BassSweepTests.Data()));
        BassKernelEnumeration.Space wide = BassKernelEnumeration.Enumerate(BassModel.Create(config, Declared(BassSweepTests.Data(), new InputDomain(-0.5, 2.0))));
        Assert.Equal(plain.States.Count, wide.States.Count);
        double[] p = plain.Posterior(), q = wide.Posterior();
        double worst = 0.0;
        for (int k = 0; k < p.Length; k++)
        {
            int at = wide.Find(plain.States[k].Terms);
            Assert.True(at >= 0, "the same structures are admitted");
            worst = Math.Max(worst, Math.Abs(p[k] - q[at]));
        }
        output.WriteLine($"{p.Length} structures; worst |Δπ| under the domain [−0.5, 2] {worst:E2}");
        Assert.True(worst < 1e-10, $"the hinge law moves with the domain by {worst:E2}");
    }
}
