using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class MultiOutputCostTests(ITestOutputHelper output)
{
    [Fact]
    public void SeparatePencilSetupProjectionsAndRepeatedMarginals()
    {
        if (!BassQualification.Enabled) return;
        const int n = 256;
        var rng = new Xoshiro256PlusPlus(20261210);
        double[][] x = [.. Enumerable.Range(0, 4).Select(_ => Enumerable.Range(0, n).Select(_ => rng.NextDouble()).ToArray())];
        foreach (var (m, k) in new[] { (8, 2), (24, 5) })
        {
            var y = new double[n, k];
            var psi = new double[k, k];
            for (int l = 0; l < k; l++)
            {
                psi[l, l] = 1;
                for (int i = 0; i < n; i++) y[i, l] = Math.Sin((l + 2) * x[l % 4][i]) + x[(l + 1) % 4][i] + 0.15 * Normal.SampleStandard(rng);
            }
            var data = new BassData([.. x.Select((v, j) => (BassPredictor)new ContinuousPredictor($"x{j}", v))],
                                    new MultiOutputResponse(y, [.. Enumerable.Range(0, k).Select(l => $"y{l}")]));
            BassTerm[] terms = [.. Enumerable.Range(0, m).Select(j => new BassTerm([new HingeFactor(j % 4, 1, 24 + 27 * (j / 4))]))];
            var replacement = new BassTerm([new HingeFactor(0, 1, 25)]);
            foreach (bool ridge in new[] { false, true })
            {
                var config = new BassConfig { Intercept = GPriorCentre.Zero,
                    Coefficients = ridge ? BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()) : BassCoefficientPrior.G(),
                    Structure = new BassStructureSpec { MaxBasis = m + 1, MaxInteraction = 1 },
                    OutputCovariance = OutputCovariancePrior.InverseWishart(k + 4, psi) };
                var model = BassModel.Create(config, data);
                var design = model.NewDesign();
                Assert.Null(design.Load(terms));
                var fit = BassStructureFit.Of(design, terms.Sum(model.Prior.LogPrior));
                const double tau = 0.07, b = 0.4;
                double sink = 0;
                object? retained = null;
                var timers = new List<(string Name, Func<int, double> Batch)>();
                void Add(string name, Action action) => timers.Add((name, BassCost.Batch(BassCost.Meter.Time, action, batch: 200)));
                Add("load and setup", () => {
                    var fresh = model.NewDesign();
                    if (fresh.Load(terms) is not null) throw new InvalidOperationException("The benchmark structure was refused.");
                    retained = fresh.OutputForms;
                });
                Add("change proposal", () => {
                    var proposal = design.ProposeChange(0, replacement);
                    if (!proposal.Admitted) throw new InvalidOperationException("The benchmark proposal was refused.");
                    retained = proposal.Outputs;
                });
                Add("cached shrinkage", () => retained = model.Target.OutputsMarginal(fit, tau, b));
                Add("full marginal", () => sink = model.Target.Evaluate(fit, Math.Log(tau)).LogLikelihood(b));
                Add("full slope", () => sink = model.Target.Evaluate(fit, Math.Log(tau)).LogLikelihoodSlope(b));
                if (ridge)
                {
                    var c = MultiOutputCoefficientTests.DenseCase(model, design, 1, 1);
                    double[,] gram = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(c.X), c.X);
                    var zero = new double[m + 1];
                    var basePencil = design.Pencil!;
                    double[,] xc = MultiOutputCoreTests.Columns(model, design);
                    var crosses = new double[k][];
                    for (int l = 0; l < k; l++)
                    {
                        crosses[l] = new double[m];
                        for (int j = 0; j < m; j++)
                            for (int i = 0; i < n; i++) crosses[l][j] += xc[i, j] * y[i, l];
                    }
                    double[] means = model.Resolved.Outputs!.Means;
                    Add("pencil decomposition", () => retained = PriorPencil.Of(gram, c.PriorPrecision, zero, m, true));
                    Add("K projections", () => retained = OutputPencil.Of(basePencil, crosses, means, GPriorCentre.Zero, n));
                    var accepted = model.NewDesign();
                    Assert.Null(accepted.Load(terms));
                    var candidate = accepted.ProposeChange(0, replacement);
                    var scored = candidate.Outputs!.Prior;
                    accepted.Accept();
                    Assert.Same(scored, accepted.OutputForms!.Prior);
                    Add("accepted form read", () => retained = accepted.OutputForms);
                }
                var readings = BassCost.Interleave(timers, repeats: 5);
                foreach (var timer in timers)
                    output.WriteLine($"M={m}, K={k}, {(ridge ? "normalized ridge" : "g-form")}, {timer.Name}: " +
                                     $"{BassCost.Median(readings[timer.Name]):F3} μs, IQR {BassCost.Iqr(readings[timer.Name]):F3}");
                output.WriteLine($"clock={(BassCost.Cycles ? "thread cycles" : "wall")}, indicative; other test hosts={string.Join(", ", BassCost.OtherTestHosts())}");
                GC.KeepAlive(sink);
                GC.KeepAlive(retained);
            }
        }
    }
}
