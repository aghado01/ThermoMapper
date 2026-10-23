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
                double[] mean = d.ConditionalMean ?? d.Coefficients;   // the coefficient mean is Rao–Blackwellized where it can be
                for (int p = 0; p < nu; p++) coef[p] += mean[p];
                lambdas.Add(d.Lambda);
                sigma += Math.Sqrt(d.Sigma2);
                if (mode == BapsLambdaUpdate.MarginalEvidence)
                {
                    // λ is the chain's state; σ² and β are conditional draws given it, and τ² = σ²/λ.
                    Assert.Equal(nu, d.ConditionalMean!.Length);
                    Assert.NotEqual(d.ConditionalMean[0], d.Coefficients[0]);
                    Assert.Equal(d.Lambda, d.Sigma2 / d.Tau2, 1e-15 * d.Lambda);
                }
                else
                {
                    Assert.Null(d.ConditionalMean);
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

    [Theory]
    [MemberData(nameof(Modes))]
    public void IntervalsAndFitBand_AreReadFromTheDraws(BapsLambdaUpdate mode)
    {
        (double[,] z, double[] y) = Data();
        int n = z.GetLength(0), nu = z.GetLength(1);
        BapsResult r = new BapsSampler(z, y, new DifferencePenalty(2)).Run(mode, chains: 3, masterSeed: 7, burn: 200, samples: 400);

        var sigma2 = new List<double>();
        var tau2 = new List<double>();
        var fits = new List<double>[n];
        var conditionalFits = new List<double>[n];
        for (int i = 0; i < n; i++) { fits[i] = new List<double>(); conditionalFits[i] = new List<double>(); }
        foreach (IReadOnlyList<BapsDraw> chain in r.PosteriorDraws)
            foreach (BapsDraw d in chain)
            {
                sigma2.Add(d.Sigma2);
                tau2.Add(d.Tau2);
                for (int i = 0; i < n; i++)
                {
                    double f = 0.0, c = 0.0;
                    for (int j = 0; j < nu; j++)
                    {
                        f += z[i, j] * d.Coefficients[j];
                        c += z[i, j] * (d.ConditionalMean ?? d.Coefficients)[j];
                    }
                    fits[i].Add(f);
                    conditionalFits[i].Add(c);
                }
            }
        sigma2.Sort();
        tau2.Sort();
        Assert.Equal(new PosteriorInterval(Type7(sigma2, 0.5), Type7(sigma2, 0.025), Type7(sigma2, 0.975)), r.Sigma2Interval);
        Assert.Equal(new PosteriorInterval(Type7(tau2, 0.5), Type7(tau2, 0.025), Type7(tau2, 0.975)), r.Tau2Interval);

        PointwiseBand band = r.FitBand(z);
        double drawSpread = 0.0, conditionalSpread = 0.0;
        for (int i = 0; i < n; i++)
        {
            double mean = 0.0;
            for (int j = 0; j < nu; j++) mean += z[i, j] * r.Coefficients[j];
            Assert.Equal(mean, band.Mean[i], 12);
            fits[i].Sort();
            Assert.Equal(Type7(fits[i], 0.025), band.Lo[i]);
            Assert.Equal(Type7(fits[i], 0.975), band.Hi[i]);
            Assert.True(band.Lo[i] < band.Mean[i] && band.Mean[i] < band.Hi[i], $"row {i}: the band must hold the mean");
            drawSpread += band.Sd[i];
            conditionalSpread += Sd(conditionalFits[i]);
        }
        // Draws carry the coefficients' own uncertainty, which conditional means leave out: in the marginal-evidence
        // mode the band is wider than the spread of β̂(λ) alone.
        if (mode == BapsLambdaUpdate.MarginalEvidence)
            Assert.True(drawSpread > 2.0 * conditionalSpread, $"draw sd {drawSpread / n:E3} against conditional-mean sd {conditionalSpread / n:E3}");
        Assert.Throws<ArgumentException>(() => r.FitBand(new double[3, nu + 1]));
    }

    private static double Type7(List<double> sorted, double p)
    {
        double pos = p * (sorted.Count - 1);
        int lo = (int)Math.Floor(pos);
        int hi = Math.Min(lo + 1, sorted.Count - 1);
        return sorted[lo] + (pos - lo) * (sorted[hi] - sorted[lo]);
    }

    private static double Sd(List<double> values)
    {
        double mean = 0.0;
        foreach (double v in values) mean += v;
        mean /= values.Count;
        double ss = 0.0;
        foreach (double v in values) ss += (v - mean) * (v - mean);
        return Math.Sqrt(ss / values.Count);
    }
}
