using System;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class MultiOutputTemperingReferenceTests(ITestOutputHelper output)
{
    [Theory]
    [InlineData(GPriorCentre.Zero)]
    [InlineData(GPriorCentre.Mean)]
    [InlineData(GPriorCentre.Flat)]
    public void NormalizedRungLawAndEvidenceAgreeWithDenseQuadrature(GPriorCentre centre)
    {
        var model = BassModel.Create(MultiOutputTemperingReference.Config(centre), MultiOutputTemperingReference.Data());
        double[] rungs = centre == GPriorCentre.Flat ? [1, 0.4, 0.05] : [1, 0.4, 0.05, 0];
        var dense = MultiOutputTemperingReference.Exact(model, 6, rungs);
        var refined = MultiOutputTemperingReference.Exact(model, 6, rungs, 0.25);
        var engine = BassTemperingTests.Exact(model, rungs);
        for (int r = 0; r < rungs.Length; r++)
        {
            output.WriteLine($"{centre}, b={rungs[r]}: ln Z={dense[r].LogZ:G12}, E[ln τ]={dense[r].LogTauMean:G9}");
            Assert.InRange(Math.Abs(dense[r].LogZ - refined[r].LogZ), 0, 1e-9);
            Assert.InRange(Math.Abs(dense[r].LogZ - engine[r].LogZ), 0, 2e-8);
            Assert.InRange(Math.Abs(dense[r].LogTauMean - engine[r].LogTauMean), 0, 2e-8);
            Assert.InRange(Math.Abs(dense[r].TauMean - engine[r].TauMean), 0, 2e-8);
            for (int m = 0; m < dense[r].Count.Length; m++)
                Assert.InRange(Math.Abs(dense[r].Count[m] - engine[r].Count[m]), 0, 1e-9);
        }
        if (centre != GPriorCentre.Flat) Assert.InRange(Math.Abs(dense[^1].LogZ), 0, 1e-12);
        Assert.Equal(centre != GPriorCentre.Flat, model.Target.HasCompleteValues);
    }

    [Theory]
    [InlineData(GPriorCentre.Zero, 0, true)]
    [InlineData(GPriorCentre.Mean, 0, true)]
    [InlineData(GPriorCentre.Flat, 0, false)]
    [InlineData(GPriorCentre.Flat, 0.03125, false)]
    [InlineData(GPriorCentre.Flat, 0.0625, false)]
    [InlineData(GPriorCentre.Flat, 0.125, true)]
    public void PropernessUsesNuPlusNbMinusTheFlatInterceptCount(GPriorCentre centre, double b, bool proper)
    {
        var model = BassModel.Create(MultiOutputTemperingReference.Config(centre, 1.5), MultiOutputTemperingReference.Data());
        Assert.Equal(proper, model.Target.Defines(b));
        var design = model.NewDesign();
        var fit = BassStructureFit.Of(design, 0);
        if (proper) Assert.True(double.IsFinite(model.Target.Evaluate(fit, -3).LogLikelihood(b)));
        else Assert.Throws<ArgumentOutOfRangeException>(() => model.Target.Evaluate(fit, -3).LogLikelihood(b));
    }
}
