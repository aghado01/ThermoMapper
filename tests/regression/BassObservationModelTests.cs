using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The observation-model seam (F-X1): the Gaussian core fits a working response y − offset with weights W, and Gaussian
/// BASS hands it the response itself and unit weights.
/// <list type="bullet">
/// <item>Gaussian hands the response through as it is: no copy, no offset, no weights, and a model built on it holds that.</item>
/// <item>An offset is the same as shifting the data, bit for bit: a design on y − offset has the statistics of a design on
/// shifted data, the resolved record reads the working response's level and spread, and a chain of a model with the offset
/// walks the shifted data's, draw for draw. Scalar and direct functional.</item>
/// <item>The contract checks its inputs, and the core refuses weights it has no Gram for.</item>
/// </list>
/// </summary>
public sealed class BassObservationModelTests
{
    private sealed class Shifted(double[] offset) : IBassObservationModel
    {
        public WorkingResponse Working(BassResponse response) => new(response, offset);

        public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot) => [];
    }

    private sealed class Weighted(double[] weights) : IBassObservationModel
    {
        public WorkingResponse Working(BassResponse response) => new(response, weights: weights);

        public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot) => [];
    }

    private static double[] Offset(int n) => [.. Enumerable.Range(0, n).Select(i => 2.5 + 0.3 * Math.Sin(0.7 * i) + (i % 5 == 0 ? 1.0 : 0.0))];

    // The data with its response less the offset: what a shifted response is, by construction.
    private static BassData Moved(BassData data, double[] offset)
    {
        ReadOnlySpan<double> y = data.Response.Values;
        BassResponse response = data.Response switch
        {
            ScalarResponse => new ScalarResponse([.. y.ToArray().Select((v, i) => v - offset[i])]),
            FunctionalResponse f => Functional(f, y, offset),
            _ => throw new InvalidOperationException(),
        };
        return new BassData(data.Predictors, response);
    }

    private static FunctionalResponse Functional(FunctionalResponse f, ReadOnlySpan<double> y, double[] offset)
    {
        var z = new double[f.Cases, f.Locations];
        for (int i = 0; i < f.Cases; i++)
            for (int r = 0; r < f.Locations; r++) z[i, r] = y[i * f.Locations + r] - offset[i * f.Locations + r];
        return new FunctionalResponse(z, f.OutputVariables);
    }

    [Fact]
    public void Gaussian_HandsTheResponseThrough_AsItIs()
    {
        BassData data = SensitivityFixture.Scalar();
        WorkingResponse working = GaussianObservation.Instance.Working(data.Response);
        Assert.Same(data.Response, working.Response);
        Assert.True(working.IsUnweighted);
        Assert.True(working.Offset.IsEmpty && working.Weights.IsEmpty);
        Assert.True(working.Values.SequenceEqual(data.Response.Values));

        BassModel model = BassModel.Create(SensitivityFixture.Config(), data);
        Assert.Same(GaussianObservation.Instance, model.Observation);
        Assert.Same(data.Response, model.Working.Response);
        Assert.True(model.Working.IsUnweighted);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void AnOffset_IsTheShiftedData_InTheDesignBitForBit(bool functional)
    {
        BassData data = functional ? SensitivityFixture.Functional() : SensitivityFixture.Scalar();
        double[] offset = Offset(data.Response.Observations);
        BassData moved = Moved(data, offset);
        BassModel gaussian = BassModel.Create(SensitivityFixture.Config(), moved);
        BassModel shifted = BassModel.Create(SensitivityFixture.Config(), data, new Shifted(offset));

        // The resolved record reads the working response: its level, its spread, and the σ² prior scaled to it.
        Assert.Equal(gaussian.Resolved.ResponseMean, shifted.Resolved.ResponseMean);
        Assert.Equal(gaussian.Resolved.ResponseMeanSquare, shifted.Resolved.ResponseMeanSquare);
        Assert.Equal(gaussian.Resolved.Sigma2?.Shape, shifted.Resolved.Sigma2?.Shape);
        Assert.Equal(gaussian.Resolved.Sigma2?.Scale, shifted.Resolved.Sigma2?.Scale);

        BassTerm[] terms = functional
            ?
            [
                new([new HingeFactor(0, 1, 2)]),
                new(functional: [new HingeFactor(0, -1, 3)]),
                new([new HingeFactor(1, -1, 4)], functional: [new HingeFactor(0, 1, 1)]),
            ]
            :
            [
                new([new HingeFactor(0, 1, 3)]),
                new([new HingeFactor(1, -1, 5)], [new CategoricalFactor(0, [1])]),
                new([new HingeFactor(0, -1, 6), new HingeFactor(1, 1, 2)]),
            ];
        BassDesign a = gaussian.NewDesign(), b = shifted.NewDesign();
        Assert.Equal(a.ResponseMean, b.ResponseMean);
        Assert.Equal(a.ResponseSquares, b.ResponseSquares);
        foreach (BassTerm term in terms)
        {
            DesignProposal pa = a.ProposeBirth(term), pb = b.ProposeBirth(term);
            Assert.Equal(pa, pb);
            Assert.True(pa.Admitted);
            a.Accept();
            b.Accept();
        }
        Assert.Equal(a.ResidualSquares, b.ResidualSquares);
        Assert.Equal(a.FittedSquares, b.FittedSquares);
        Assert.True(a.Coefficients.SequenceEqual(b.Coefficients));
        Assert.True(a.Residual.SequenceEqual(b.Residual));
    }

    [Fact]
    public void AModelWithAnOffset_WalksTheShiftedDatasChain_DrawForDraw()
    {
        BassData data = SensitivityFixture.Scalar();
        double[] offset = Offset(data.Response.Observations);
        BassModel gaussian = BassModel.Create(SensitivityFixture.Config(), Moved(data, offset));
        BassModel shifted = BassModel.Create(SensitivityFixture.Config(), data, new Shifted(offset));
        var a = new BassChain(gaussian, [], Math.Log(0.05), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        var b = new BassChain(shifted, [], Math.Log(0.05), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        for (int step = 0; step < 150; step++)
        {
            BassState sa = a.Step(), sb = b.Step();
            Assert.Equal(sa.Count, sb.Count);
            Assert.Equal(sa.LogTau, sb.LogTau);
            Assert.Equal(sa.Fit.ResidualSquares, sb.Fit.ResidualSquares);
            Assert.Equal(sa.Fit.FittedSquares, sb.Fit.FittedSquares);
        }
        Assert.True(a.Current.Count > 0, "the chain should have taken births");
    }

    [Fact]
    public void TheContract_ChecksItsInputs_AndTheCoreRefusesWeights()
    {
        BassData data = SensitivityFixture.Scalar();
        int n = data.Response.Observations;
        Assert.Throws<ArgumentException>(() => new WorkingResponse(data.Response, new double[n - 1]));
        Assert.Throws<ArgumentException>(() => new WorkingResponse(data.Response, [.. Enumerable.Repeat(double.NaN, n)]));
        Assert.Throws<ArgumentException>(() => new WorkingResponse(data.Response, weights: new double[n - 1]));
        Assert.Throws<ArgumentException>(() => new WorkingResponse(data.Response, weights: new double[n]));   // zero weights
        Assert.Throws<ArgumentException>(() => new WorkingResponse(data.Response, weights: [.. Enumerable.Repeat(-1.0, n)]));

        var weighted = new WorkingResponse(data.Response, weights: [.. Enumerable.Repeat(2.0, n)]);
        Assert.False(weighted.IsUnweighted);
        Assert.Equal(n, weighted.Weights.Length);
        BassModel gaussian = BassModel.Create(SensitivityFixture.Config(), data);
        Assert.Throws<NotSupportedException>(() => new BassDesign(gaussian.Basis, data, gaussian.Admission, working: weighted));
        Assert.Throws<NotSupportedException>(() => BassModel.Create(SensitivityFixture.Config(), data, new Weighted([.. Enumerable.Repeat(2.0, n)])).NewDesign());

        // A working response is of its own response.
        BassData other = SensitivityFixture.Scalar(seed: 8);
        Assert.Throws<ArgumentException>(() => new BassDesign(gaussian.Basis, other, gaussian.Admission, working: new WorkingResponse(data.Response)));
    }
}
