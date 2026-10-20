using System;
using System.Threading.Tasks;
using Maths.Regression.Bgp;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// <see cref="GpRegression.Fit"/> must stay reentrant: an ensemble drives many bandwidth chains concurrently over
/// one shared instance. Concurrent fits and predictions over immutable inputs and a stateless kernel must equal the
/// serial ones bit for bit. This protects the local-scratch repair of the historical shared-buffer defect (R06).
/// </summary>
public sealed class BgpConcurrencyTests
{
    [Fact]
    public void ConcurrentFits_MatchSerialFits_BitForBit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 606);
        const int n = 160, dim = 2, bandwidths = 24, repeats = 8;
        var x = new double[n, dim];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = rng.NextDouble();
            x[i, 1] = rng.NextDouble();
            y[i] = Math.Sin(6.0 * x[i, 0]) * Math.Cos(4.0 * x[i, 1]) + 0.05 * (rng.NextDouble() - 0.5);
        }
        var xCopy = (double[,])x.Clone();
        var yCopy = (double[])y.Clone();
        var xTest = new double[,] { { 0.1, 0.9 }, { 0.5, 0.5 }, { 0.75, 0.2 } };

        var gp = new GpRegression(x, y, sigma2: 0.01, new SquaredExponentialKernel());
        var ts = new double[bandwidths];
        for (int b = 0; b < bandwidths; b++) ts[b] = 0.005 * Math.Pow(100.0, b / (bandwidths - 1.0));

        var serial = new GpFit[bandwidths];
        var serialPred = new double[bandwidths][];
        for (int b = 0; b < bandwidths; b++)
        {
            serial[b] = gp.Fit(ts[b]);
            serialPred[b] = gp.PredictMean(serial[b], xTest);
        }

        var concurrent = new GpFit[repeats * bandwidths];
        var concurrentPred = new double[repeats * bandwidths][];
        Parallel.For(0, repeats * bandwidths, new ParallelOptions { MaxDegreeOfParallelism = Math.Max(2, Environment.ProcessorCount) },
            k =>
            {
                GpFit fit = gp.Fit(ts[k % bandwidths]);
                concurrent[k] = fit;
                concurrentPred[k] = gp.PredictMean(fit, xTest);
            });

        for (int k = 0; k < repeats * bandwidths; k++)
        {
            int b = k % bandwidths;
            Assert.Equal(serial[b].Bandwidth, concurrent[k].Bandwidth);
            Assert.Equal(serial[b].LogMarginal, concurrent[k].LogMarginal);
            Assert.Equal(serial[b].Alpha, concurrent[k].Alpha);
            Assert.Equal(serialPred[b], concurrentPred[k]);
        }

        // The inputs are read, never written.
        Assert.Equal(xCopy, x);
        Assert.Equal(yCopy, y);
    }
}
