using System;
using System.Collections.Generic;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// The retained BAPS draws (<see cref="BapsResult.PosteriorDraws"/>): every draw of every chain is kept with its
/// coefficients, λ, σ², τ² and (adaptive mode) η, and the reported summaries are exactly what those draws give.
/// </summary>
public sealed class BapsDrawTests
{
    private static (double[,] Z, double[] Y) Data()
    {
        var rng = new Xoshiro256PlusPlus(seed: 31);
        int n = 60;
        var x = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = (i + 0.5) / n;
            double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
            y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.1 * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
        var knots = new double[12];
        for (int k = 0; k < knots.Length; k++) knots[k] = (k + 1.0) / (knots.Length + 1);
        return (new SplineBasis(3).Design(new KnotConfig(knots), x), y);
    }

    public static TheoryData<BapsLambdaUpdate> Modes => new()
    {
        BapsLambdaUpdate.Gibbs, BapsLambdaUpdate.MarginalEvidence, BapsLambdaUpdate.AdaptiveGibbs,
    };

    [Theory]
    [MemberData(nameof(Modes))]
    public void PosteriorDraws_KeepEveryDraw_AndGiveTheSummaries(BapsLambdaUpdate mode)
    {
        (double[,] z, double[] y) = Data();
        int nu = z.GetLength(1), nLoc = nu - 2;
        BapsResult r = new BapsSampler(z, y, new DifferencePenalty(2)).Run(mode, chains: 3, masterSeed: 7, burn: 200, samples: 400);

        Assert.Equal(3, r.PosteriorDraws.Count);
        var coef = new double[nu];
        var local = new double[nLoc];
        var lambdas = new List<double>();
        double sigma = 0.0;
        foreach (IReadOnlyList<BapsDraw> chain in r.PosteriorDraws)
        {
            Assert.Equal(400, chain.Count);
            foreach (BapsDraw d in chain)
            {
                Assert.Equal(nu, d.Coefficients.Length);
                for (int p = 0; p < nu; p++) coef[p] += d.Coefficients[p];
                lambdas.Add(d.Lambda);
                sigma += Math.Sqrt(d.Sigma2);
                if (mode == BapsLambdaUpdate.MarginalEvidence)
                {
                    Assert.True(double.IsNaN(d.Tau2));   // λ is sampled alone
                }
                else
                {
                    Assert.Equal(d.Sigma2 / d.Tau2, d.Lambda);
                }
                if (mode == BapsLambdaUpdate.AdaptiveGibbs)
                {
                    Assert.Equal(nLoc, d.Eta!.Length);
                    for (int k = 0; k < nLoc; k++) local[k] += d.Lambda * d.Eta[k];
                }
                else
                {
                    Assert.Null(d.Eta);
                }
            }
        }
        // Distinct draws are distinct arrays, not one reused buffer.
        Assert.NotSame(r.PosteriorDraws[0][0].Coefficients, r.PosteriorDraws[0][1].Coefficients);

        for (int p = 0; p < nu; p++) Assert.Equal(coef[p] / r.Draws, r.Coefficients[p], 12);
        Assert.Equal(sigma / r.Draws, r.SigmaMean, 12);
        lambdas.Sort();
        Assert.Equal(lambdas[(lambdas.Count - 1) / 2] + 0.5 * (lambdas[lambdas.Count / 2] - lambdas[(lambdas.Count - 1) / 2]),
                     r.LambdaMedian);
        if (mode == BapsLambdaUpdate.AdaptiveGibbs)
            for (int k = 0; k < nLoc; k++) Assert.Equal(local[k] / r.Draws, r.LocalSmoothing![k], 12);
        else
            Assert.Null(r.LocalSmoothing);
    }
}
