using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using Maths.Regression.Sensitivity;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class SampledVarianceFormTests(ITestOutputHelper output)
{
    public static TheoryData<string> Readouts() => ["scalar", "directional", "named", "integrated", "conditioned", "partly conditioned", "basis"];

    private static SobolDecomposition[] Read(string kind, SensitivitySpec spec)
    {
        if (kind == "basis")
        {
            BasisSobolAnalysis a = BassOutputComponentsTests.FitOf().Sensitivity(spec);
            return [a.AtLocation(0), a.AtLocation(7), a.Integrated];
        }
        if (kind == "named")
        {
            var (fit, data) = MultiOutputReadoutTests.Fixed();
            return [.. BassMultiOutput.Sensitivity(fit, data, new double[] { 1.0, -0.4 }, spec).Decompositions];
        }
        if (kind == "directional")
        {
            BassData data = SampledSensitivityTests.DirectionalData();
            return [.. SampledSensitivityTests.DirectionalFit(data, 1).Sensitivity(data, spec).Decompositions];
        }
        var (f, d) = kind switch
        {
            "scalar" => (SensitivityFixture.FitScalar(), SensitivityFixture.Scalar()),
            "partly conditioned" => (SensitivityFixture.FitTwoOutputs(), SensitivityFixture.TwoOutputs()),
            _ => (SensitivityFixture.FitFunctional(), SensitivityFixture.Functional()),
        };
        if (kind == "conditioned") spec = spec with { Functional = FunctionalForm.Condition };
        if (kind == "partly conditioned") spec = spec with { Functional = FunctionalForm.Condition, ConditionOn = ["t0"] };
        return [.. f.Sensitivity(d, spec).Decompositions];
    }

    // All reported numeric tables, indices and status bits, in a stable binary order.
    private static string Fingerprint(IEnumerable<SobolDecomposition> decompositions)
    {
        using var stream = new MemoryStream();
        using var writer = new BinaryWriter(stream);
        foreach (SobolDecomposition d in decompositions)
        {
            writer.Write(d.Count);
            foreach (string input in d.Inputs) writer.Write(input);
            foreach (int[] subset in d.Subsets) { writer.Write(subset.Length); foreach (int v in subset) writer.Write(v); }
            foreach (Array? table in new Array?[] { d.Variance, d.MainVariance, d.TotalVariance, d.OtherVariance,
                d.VarianceError, d.MainVarianceError, d.TotalVarianceError, d.MainIndexError, d.TotalIndexError })
            {
                writer.Write(table is not null);
                if (table is not null) foreach (double value in table) writer.Write(value);
            }
            for (int j = 0; j < d.Count; j++)
            {
                writer.Write((int)d.Status[j]);
                for (int k = 0; k < d.Subsets.Count; k++) writer.Write(d.MainIndex(j, k));
                for (int k = 0; k < d.Inputs.Count; k++) writer.Write(d.TotalIndex(j, k));
            }
        }
        return Convert.ToHexString(SHA256.HashData(stream.ToArray()));
    }

    [Theory]
    [MemberData(nameof(Readouts))]
    [Trait("Tier", "Pinned")]
    public void Saltelli_ReproducesThePriorPublicTables(string kind)
    {
        // Captured from the public readouts at a83db0c, before SensitivitySpec exposed the form.
        string expected = kind switch
        {
            "scalar" => "A7654FB6D07BA3D6BC1E3EF97A80E066D0DB8FADAA46B8464910975E27F2AE04",
            "directional" => "99E4E0F6AE8A40C5D0CA56482C0194DFD58EDCF3CA736A57E81A60D04C8D8145",
            "named" => "E174D499BEA71779907D35EEA6BD3A201F12B1CA8480BF161DF4B8976AD12BEB",
            "integrated" => "5E7D328785E783D930ED1EC984D568E878595BFA6684A859D002F92C6F42B6C0",
            "conditioned" => "14B911C8A79DD5571E4DC393652E354A21FE183765E9D1B916190C3FBEFCD4FE",
            "partly conditioned" => "793CFD7E908D1920320AC7CA879601D74B5654A4894B1B636738CEB5C02A4CDA",
            "basis" => "4600AB9DE4DB99335BFE8492BCA5DD8F50A2428B11ACA270A39CE1A8AE2EE128",
            _ => throw new ArgumentException(kind),
        };
        var spec = new SensitivitySpec { Route = SensitivityRoute.Sampled, Points = 256, Seed = 739, Draws = [0, 1], MainVarianceForm = MainVarianceForm.Saltelli };
        string actual = Fingerprint(Read(kind, spec));
        output.WriteLine($"{kind}: {actual}");
        Assert.Equal(expected, actual);
    }

    [Theory]
    [MemberData(nameof(Readouts))]
    public void OmittedForm_IsTwoDifference_WithMatchedPoints(string kind)
    {
        var spec = new SensitivitySpec { Route = SensitivityRoute.Sampled, Points = 256, Seed = 739, Draws = [0, 1] };
        Assert.Equal(MainVarianceForm.TwoDifference, spec.MainVarianceForm);
        Assert.Equal(Fingerprint(Read(kind, spec with { MainVarianceForm = MainVarianceForm.TwoDifference })), Fingerprint(Read(kind, spec)));
        // Auto reaches the same sampled design for terms without factors.
        if (kind == "directional")
            foreach (MainVarianceForm form in new[] { MainVarianceForm.Saltelli, MainVarianceForm.TwoDifference })
                Assert.Equal(Fingerprint(Read(kind, spec with { MainVarianceForm = form })),
                             Fingerprint(Read(kind, spec with { Route = SensitivityRoute.Auto, MainVarianceForm = form })));
    }

    [Theory]
    [InlineData("scalar")]
    [InlineData("named")]
    [InlineData("integrated")]
    [InlineData("conditioned")]
    [InlineData("partly conditioned")]
    [InlineData("basis")]
    public void AnalyticTables_IgnoreTheSampledForm(string kind)
    {
        var spec = new SensitivitySpec { Route = SensitivityRoute.Exact, Draws = [0, 1] };
        string exact = Fingerprint(Read(kind, spec));
        Assert.Equal(exact, Fingerprint(Read(kind, spec with { MainVarianceForm = MainVarianceForm.Saltelli })));
        Assert.Equal(exact, Fingerprint(Read(kind, spec with { Route = SensitivityRoute.Auto })));
        Assert.Equal(exact, Fingerprint(Read(kind, spec with { Route = SensitivityRoute.Auto, MainVarianceForm = MainVarianceForm.Saltelli })));
    }

    [Theory]
    [InlineData(MainVarianceForm.Saltelli)]
    [InlineData(MainVarianceForm.TwoDifference)]
    public void SampledForms_RespectTheSubsetBudget(MainVarianceForm form)
    {
        var spec = new SensitivitySpec { Route = SensitivityRoute.Sampled, MainVarianceForm = form, Points = 16, SubsetBudget = 1 };
        Assert.Throws<ArgumentException>(() => Read("directional", spec));
    }

    [Fact]
    public void UnknownSampledForm_IsRefused()
    {
        var spec = new SensitivitySpec { Route = SensitivityRoute.Sampled, MainVarianceForm = (MainVarianceForm)99, Points = 16, Draws = [0] };
        foreach (string kind in new[] { "scalar", "named", "integrated", "basis" })
            Assert.Throws<ArgumentOutOfRangeException>(() => Read(kind, spec));
    }
}
