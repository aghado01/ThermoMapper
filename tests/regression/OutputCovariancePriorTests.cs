using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// Σ's prior for several outputs (E34's decision 2 and its options, <see cref="OutputCovariancePrior"/>): its JSON form, which a config
/// without it leaves out, so every existing fingerprint holds; its checks; its resolution per kind; and the target under a declared prior
/// against the matrix-t reference.
/// </summary>
public sealed class OutputCovariancePriorTests
{
    private static readonly double[,] Declared = { { 2.0, 0.4, -0.1 }, { 0.4, 1.5, 0.2 }, { -0.1, 0.2, 0.8 } };

    [Fact]
    public void Config_RoundTripsEachKind_AndLeavesItOutWhenNull()
    {
        Assert.DoesNotContain("OutputCovariance", new BassConfig().ToJson());
        foreach (OutputCovariancePrior prior in new[]
                 {
                     OutputCovariancePrior.ScaledToResponse(), OutputCovariancePrior.ScaledToResponse(0.5), OutputCovariancePrior.Standardized(0.01),
                     OutputCovariancePrior.InverseWishart(5.0, Declared),
                 })
        {
            var config = new BassConfig { OutputCovariance = prior };
            BassConfig back = BassConfig.FromJson(config.ToJson());
            Assert.Equal(prior, back.OutputCovariance);
            Assert.Equal(config.ToJson(), back.ToJson());
        }
        Assert.Equal(OutputCovariancePrior.InverseWishart(5.0, Declared), OutputCovariancePrior.InverseWishart(5.0, (double[,])Declared.Clone()));
        Assert.NotEqual(OutputCovariancePrior.InverseWishart(5.0, Declared), OutputCovariancePrior.InverseWishart(5.5, Declared));
    }

    [Fact]
    public void Prior_RefusesParametersItsKindDoesNotTake()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OutputCovariancePrior.ScaledToResponse(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => OutputCovariancePrior.Standardized(double.PositiveInfinity));
        Assert.Throws<ArgumentException>(() => new OutputCovariancePrior(OutputCovariancePriorKind.ScaledToResponse, 1e-3, degreesOfFreedom: 4.0));
        Assert.Throws<ArgumentException>(() => new OutputCovariancePrior(OutputCovariancePriorKind.InverseWishart, shape: 1e-3));
        Assert.Throws<ArgumentOutOfRangeException>(() => OutputCovariancePrior.InverseWishart(2.0, Declared));   // ν₀ ≤ K − 1
        Assert.Throws<ArgumentException>(() => OutputCovariancePrior.InverseWishart(5.0, new[,] { { 1.0, 0.3 }, { 0.2, 1.0 } }));
        Assert.Throws<ArgumentException>(() => OutputCovariancePrior.InverseWishart(5.0, new[,] { { 1.0, 2.0 }, { 2.0, 1.0 } }));
        Assert.Throws<ArgumentException>(() => new OutputCovariancePrior(OutputCovariancePriorKind.InverseWishart, degreesOfFreedom: 5.0, scale: [[1.0, 0.0], [0.0]]));
    }

    [Fact]
    public void Resolution_ReadsEachKind_AndRefusesWhatDoesNotFit()
    {
        const int n = 50, k = 3;
        (List<BassPredictor> predictors, double[,] y) = MultiOutputCoreTests.Data(n, k, 20261113);
        var data = new BassData(predictors, new MultiOutputResponse(y, ["a", "b", "c"]));
        double[,] squares = BassModel.Create(new BassConfig(), data).Resolved.Outputs!.MeanSquares;

        InverseWishart standardized = BassModel.Create(new BassConfig { OutputCovariance = OutputCovariancePrior.Standardized(0.25) }, data).Resolved.Outputs!.Prior;
        Assert.Equal(k - 1 + 0.5, standardized.DegreesOfFreedom);
        double[,] psi = standardized.Scale;
        for (int a = 0; a < k; a++)
            for (int d = 0; d < k; d++) Assert.Equal(a == d ? 0.5 * squares[a, a] : 0.0, psi[a, d]);

        BassResolved declared = BassModel.Create(new BassConfig { OutputCovariance = OutputCovariancePrior.InverseWishart(5.0, Declared), Sigma2 = Sigma2Prior.Jeffreys }, data).Resolved;
        Assert.Equal(5.0, declared.Outputs!.Prior.DegreesOfFreedom);
        Assert.Equal(Declared, declared.Outputs.Prior.Scale);
        Assert.Contains(declared.AutoPicks, p => p.Setting == "Sigma2" && p.Value == "not read");

        // Another K, OutputCovariance on one response, and collinear outputs: refused under decision 2's form alone.
        Assert.Throws<ArgumentException>(() => BassModel.Create(new BassConfig { OutputCovariance = OutputCovariancePrior.InverseWishart(5.0, new[,] { { 1.0, 0.0 }, { 0.0, 1.0 } }) }, data));
        Assert.Throws<ArgumentException>(() => BassModel.Create(new BassConfig { OutputCovariance = OutputCovariancePrior.Standardized() },
                                                                new BassData(predictors, new ScalarResponse([.. Enumerable.Range(0, n).Select(i => y[i, 0])]))));
        var collinear = (double[,])y.Clone();
        for (int i = 0; i < n; i++) collinear[i, 2] = 2.0 * y[i, 0] - 0.5 * y[i, 1] + 7.0;
        var collinearData = new BassData(predictors, new MultiOutputResponse(collinear, ["a", "b", "c"]));
        Assert.Throws<ArgumentException>(() => BassModel.Create(new BassConfig(), collinearData));
        Assert.NotNull(BassModel.Create(new BassConfig { OutputCovariance = OutputCovariancePrior.Standardized() }, collinearData).Resolved.Outputs);
        Assert.NotNull(BassModel.Create(new BassConfig { OutputCovariance = OutputCovariancePrior.InverseWishart(5.0, Declared) }, collinearData).Resolved.Outputs);
    }

    [Fact]
    public void Target_UnderADeclaredPrior_IsTheMatrixTReference()
    {
        const int n = 50, k = 3;
        (List<BassPredictor> predictors, double[,] y) = MultiOutputCoreTests.Data(n, k, 20261114);
        BassModel model = BassModel.Create(new BassConfig { OutputCovariance = OutputCovariancePrior.InverseWishart(5.0, Declared) },
                                           new BassData(predictors, new MultiOutputResponse(y, ["a", "b", "c"])));
        BassDesign design = model.NewDesign();
        var rng = new Xoshiro256PlusPlus(20261114);
        double worst = 0.0;
        int compared = 0;
        for (int index = 0; index < 40; index++)
        {
            var terms = new List<BassTerm>();
            while (terms.Count < 1 + rng.NextInt(4))
            {
                int v = rng.NextInt(2);
                var term = new BassTerm([new HingeFactor(v, rng.NextDouble() < 0.5 ? -1 : 1, rng.NextInt(model.Inputs.ContinuousKnots(v).Length))], [], []);
                if (!terms.Contains(term)) terms.Add(term);
            }
            if (design.Load(terms) is not null) continue;
            double logTau = Math.Log(0.02 + rng.NextDouble());
            double engine = model.Target.Evaluate(BassStructureFit.Of(design, 0.0), logTau).LogMarginal;
            var xc = new double[n, design.Count];
            var column = new double[n];
            for (int j = 0; j < design.Count; j++)
            {
                model.Basis.CaseColumn(design.Terms[j], model.Inputs.Cases, column);
                double mean = column.Average();
                for (int i = 0; i < n; i++) xc[i, j] = column[i] - mean;
            }
            double reference = MatrixTReference.LogMarginal(MatrixTReference.Collapse(MultiOutputCase.GForm(y, xc, GPriorCentre.Mean, Math.Exp(logTau), 1.0)),
                                                            new ReferenceSigmaPrior(5.0, Declared));
            worst = Math.Max(worst, Math.Abs(engine - reference) / Math.Max(1.0, Math.Abs(reference)));
            compared++;
        }
        Assert.True(compared > 20, $"only {compared} structures loaded");
        Assert.True(worst < 1e-12, $"worst relative gap {worst:E2}");
    }
}
