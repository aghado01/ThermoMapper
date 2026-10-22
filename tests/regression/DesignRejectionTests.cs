using System;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Bars;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Every refused knot configuration is counted by why it was refused (step 3.4), so the truncation that the support
/// rules and admission impose stays visible: off the simplex, a zero count prior, too few observations, an
/// inadmissible factor, an IRLS fit that did not settle, a zero response.
/// </summary>
public sealed class DesignRejectionTests
{
    private static double[] Even(int n)
    {
        var x = new double[n];
        for (int i = 0; i < n; i++) x[i] = (i + 0.5) / n;
        return x;
    }

    private static DesignRejection OutcomeOf(IObservationModel model, IComplexityPrior prior, double[] x, double[] y, double[] knots)
    {
        var tally = new Tally<DesignRejection>();
        ITargetEvaluation e = new SplineTarget(new SplineBasis(3), model, prior, x, y, outcomes: tally).Evaluate(new KnotConfig(knots));
        Assert.Equal(1, tally.Total);
        foreach ((DesignRejection reason, long count) in tally.Entries)
            if (count == 1)
            {
                Assert.Equal(reason == DesignRejection.None, double.IsFinite(Temperature.LogPosterior(e)));
                return reason;
            }
        throw new InvalidOperationException("No outcome was counted.");
    }

    [Fact]
    public void EachReason_IsCounted()
    {
        double[] x = Even(40);
        var y = new double[x.Length];
        for (int i = 0; i < y.Length; i++) y[i] = Math.Sin(2.0 * Math.PI * x[i]);
        var normal = new WeightedNormalModel();
        var poisson = new PoissonPrior(4.0);

        Assert.Equal(DesignRejection.None, OutcomeOf(normal, poisson, x, y, [0.3, 0.6]));
        Assert.Equal(DesignRejection.OutsideSimplex, OutcomeOf(normal, poisson, x, y, [0.6, 0.3]));
        Assert.Equal(DesignRejection.ZeroCountPrior, OutcomeOf(normal, new UniformPrior(1), x, y, [0.3, 0.6]));
        Assert.Equal(DesignRejection.TooFewObservations, OutcomeOf(normal, poisson, Even(6), new double[6], [0.3, 0.6]));
        Assert.Equal(DesignRejection.ZeroResponse, OutcomeOf(normal, poisson, x, new double[x.Length], [0.3, 0.6]));

        // Knots packed into a gap in the data leave cubic B-splines with nothing under them.
        var gapped = new double[40];
        for (int i = 0; i < 20; i++) { gapped[i] = 0.4 * (i + 0.5) / 20; gapped[20 + i] = 0.6 + 0.4 * (i + 0.5) / 20; }
        var gy = new double[40];
        for (int i = 0; i < 40; i++) gy[i] = Math.Sin(2.0 * Math.PI * gapped[i]);
        Assert.Equal(DesignRejection.Inadmissible,
                     OutcomeOf(normal, poisson, gapped, gy, [0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58]));

        var counts = new double[x.Length];
        for (int i = 0; i < counts.Length; i++) counts[i] = Math.Round(6.0 + 4.0 * Math.Sin(2.0 * Math.PI * x[i]));
        Assert.Equal(DesignRejection.NotConverged, OutcomeOf(new PoissonModel(maxIterations: 1), poisson, x, counts, [0.3, 0.6]));
        Assert.Equal(DesignRejection.None, OutcomeOf(new PoissonModel(), poisson, x, counts, [0.3, 0.6]));
    }

    [Fact]
    public void BarsResult_ReportsEveryEvaluationsOutcome()
    {
        double[] x = Even(80);
        var y = new double[x.Length];
        for (int i = 0; i < y.Length; i++) y[i] = Math.Sin(2.0 * Math.PI * x[i]) + 0.1 * Math.Cos(17.0 * i);

        // A cap of two knots: births beyond it are proposed and refused by the count prior.
        BarsResult r = Bars.Run(new BarsConfig { Prior = new UniformPrior(2), Chains = 2, MasterSeed = 3, BurnIn = 200, MaxSamples = 300 },
                                x, y, x);
        Tally<DesignRejection> e = r.Evaluations;
        Assert.True(e.Count(DesignRejection.None) > 0);
        Assert.True(e.Count(DesignRejection.ZeroCountPrior) > 0);
        Assert.InRange(e.Total, 2 * (200 + 300) / 2, 2 * (200 + 300) + 2);   // at most one per step, plus each start
    }
}
