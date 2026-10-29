using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Presentation;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Sensitivity summaries (stage 9, O9; P11; R's <c>plot.bassSob</c>): each index's mean and quantiles over the defined draws
/// against <see cref="DrawSample"/>, the zero-variance draws left out, under equal and unequal weights, the means the
/// decomposition's own bit for bit; the variances over every draw; a direct functional fit's conditioned form per
/// coordinate; and a basis fit's per location and integrated.
/// </summary>
public sealed class BassPresentationSensitivityTests
{
    private readonly ITestOutputHelper _out;
    public BassPresentationSensitivityTests(ITestOutputHelper output) => _out = output;

    // A summary of `decomposition` against DrawSample over the draws `include` keeps, under the decomposition's weights.
    private static void AssertAgainstDraws(SobolSummary s, SobolDecomposition decomposition)
    {
        int n = decomposition.Count;
        double Weight(int d) => decomposition.Weights?[d] ?? 1.0;
        PosteriorSummary Reference(Func<int, double> value, bool definedOnly)
        {
            var sample = new DrawSample();
            double sum = 0.0, mass = 0.0;
            for (int d = 0; d < n; d++)
            {
                if (definedOnly && !decomposition.IsDefined(d)) continue;
                sample.Add(value(d), Weight(d));
                sum += Weight(d) * value(d);
                mass += Weight(d);
            }
            return new PosteriorSummary(sum / mass, double.NaN, double.NaN, double.NaN, s.Variance.Probabilities,
                                        [.. s.Variance.Probabilities.Select(sample.Quantile)], sample.Count);
        }
        void Same(PosteriorSummary expected, PosteriorSummary actual)
        {
            Assert.Equal(expected.Mean, actual.Mean, 14);
            Assert.Equal(expected.Quantiles, actual.Quantiles);
        }

        double[] meanMain = decomposition.MeanMainIndices(), meanTotal = decomposition.MeanTotalIndices();
        for (int u = 0; u < decomposition.Subsets.Count; u++)
        {
            Same(Reference(d => decomposition.MainIndex(d, u), true), s.MainIndex[u]);
            Assert.Equal(meanMain[u], s.MainIndex[u].Mean);   // the view's own mean, bit for bit
            Same(Reference(d => decomposition.MainVariance[d, u], false), s.MainVariance[u]);
        }
        for (int v = 0; v < decomposition.Inputs.Count; v++)
        {
            Same(Reference(d => decomposition.TotalIndex(d, v), true), s.TotalIndex[v]);
            Assert.Equal(meanTotal[v], s.TotalIndex[v].Mean);
            Same(Reference(d => decomposition.TotalVariance[d, v], false), s.TotalVariance[v]);
        }
        Same(Reference(d => decomposition.Variance[d], false), s.Variance);
        Assert.Equal(n, s.Draws);
        Assert.Equal(Enumerable.Range(0, n).Count(decomposition.IsDefined), s.Defined);
        Assert.Equal(decomposition.Subsets, s.Subsets);
        Assert.Equal(decomposition.Inputs, s.Inputs);
    }

    [Fact]
    public void DirectSummary_AgreesWithDrawSample_OverTheDefinedDraws_UnderEqualAndUnequalWeights()
    {
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        foreach (double[]? weights in new[] { null, BassPresentationTests.Weights(fit.Draws.Sum(c => c.Count)) })
        {
            BassFit weighted = fit with { DrawWeights = weights };
            SobolAnalysis analysis = weighted.Sensitivity(data);
            SobolAnalysisSummary summary = weighted.SensitivitySummary(data);
            SobolDecomposition decomposition = Assert.Single(analysis.Decompositions);
            AssertAgainstDraws(Assert.Single(summary.Decompositions), decomposition);
            Assert.Null(summary.Form);
            Assert.Null(summary.Coordinates);
            Assert.Empty(summary.OutputVariables);
            Assert.Equal(analysis.MaxOrder, summary.MaxOrder);
        }
        SobolSummary s = fit.SensitivitySummary(data).Decompositions[0];
        _out.WriteLine(string.Join("; ", s.Subsets.Select((u, k) =>
            $"{{{string.Join(",", u.Select(v => s.Inputs[v]))}}} {s.MainIndex[k].Mean:F3} [{s.MainIndex[k].Quantiles[0]:F3}, {s.MainIndex[k].Quantiles[^1]:F3}]")));
    }

    [Fact]
    public void DirectSummary_LeavesTheZeroVarianceDrawsOut_OfTheIndicesAlone()
    {
        // Two intercept-only draws beside the fit's: their surfaces are constant, so their indices are undefined, while
        // their variances, all zero, are values.
        BassFit fit = SensitivityFixture.FitScalar();
        BassData data = SensitivityFixture.Scalar();
        var store = fit.Draws[0].Copy();
        store.Add(new BassDraw(new BassStructure([]), 1.0, 0.01, [1.0], 1.0));
        store.Add(new BassDraw(new BassStructure([]), 1.0, 0.01, [2.0], 1.0));
        BassFit mixed = fit with { Draws = [store, fit.Draws[1]] };
        SobolDecomposition decomposition = mixed.Sensitivity(data).Decompositions[0];
        SobolSummary s = mixed.SensitivitySummary(data, probabilities: [0.0, 0.5, 1.0]).Decompositions[0];
        Assert.Equal(decomposition.Count - 2, s.Defined);
        Assert.Equal(s.Defined, s.MainIndex[0].Draws);
        Assert.Equal(s.Draws, s.Variance.Draws);
        Assert.Equal(0.0, s.Variance.Min);
        AssertAgainstDraws(s, decomposition);
    }

    [Fact]
    public void ConditionedSummary_IsOnePerCoordinate()
    {
        BassFit fit = SensitivityFixture.FitFunctional();
        BassData data = SensitivityFixture.Functional();
        var spec = new SensitivitySpec { Functional = FunctionalForm.Condition, Coordinates = [[0.0], [0.5], [1.2]] };
        SobolAnalysis analysis = fit.Sensitivity(data, spec);
        SobolAnalysisSummary summary = fit.SensitivitySummary(data, spec);
        Assert.Equal(FunctionalForm.Condition, summary.Form);
        Assert.Equal(new[] { "t" }, summary.OutputVariables);
        Assert.Equal(3, summary.Decompositions.Count);
        Assert.Equal(new[] { false, false, true }, summary.CoordinatesOutside);
        for (int c = 0; c < 3; c++)
        {
            Assert.Equal(analysis.Coordinates![c], summary.Coordinates![c]);
            AssertAgainstDraws(summary.Decompositions[c], analysis.Decompositions[c]);
        }
        SobolSummary integrated = Assert.Single(fit.SensitivitySummary(data).Decompositions);
        Assert.Equal(new[] { "x0", "x1", "t" }, integrated.Inputs);
    }

    [Fact]
    public void BasisSummary_IsPerLocationAndIntegrated()
    {
        BassBasisFit fit = BassOutputComponentsTests.FitOf();
        BasisSobolAnalysis analysis = fit.Sensitivity();
        BasisSobolSummary summary = fit.SensitivitySummary();
        Assert.Equal(analysis.Locations, summary.PerLocation.Count);
        foreach (int l in new[] { 0, 7, analysis.Locations - 1 }) AssertAgainstDraws(summary.PerLocation[l], analysis.AtLocation(l));
        AssertAgainstDraws(summary.Integrated, analysis.Integrated);
        Assert.Equal(BasisSobolAnalysis.LocationInput, summary.Integrated.Inputs[^1]);
        Assert.Equal((analysis.MaxOrder, analysis.IntegratedMaxOrder), (summary.MaxOrder, summary.IntegratedMaxOrder));
        IReadOnlyList<double> curve = [.. summary.PerLocation.Select(s => s.MainIndex[0].Mean)];
        _out.WriteLine($"S(z0) across the locations: {string.Join(", ", curve.Select(v => v.ToString("F2")))}");
    }
}
