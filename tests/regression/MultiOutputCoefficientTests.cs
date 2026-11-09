using System;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class MultiOutputCoefficientTests(ITestOutputHelper output)
{
    internal static BassConfig Config(string form, GPriorCentre centre = GPriorCentre.Zero, int outputs = 2)
    {
        CoefficientForm covariance = form switch { "ridge" => new RidgeForm(), "raw" => new RawRidgeForm(), _ => new ModifiedGForm(0.8) };
        var psi = new double[outputs, outputs];
        for (int a = 0; a < outputs; a++)
            for (int d = 0; d < outputs; d++) psi[a, d] = a == d ? 1.0 + a : 0.2;
        return new BassConfig {
            Structure = new BassStructureSpec { MaxBasis = 4, MaxInteraction = 2, Legendre = form == "modified" ? new LegendreBasis { MaxDegree = 3 } : null },
            Coefficients = new BassCoefficientPrior(covariance, TauPrior.UnitInformationGamma()),
            Intercept = centre, OutputCovariance = OutputCovariancePrior.InverseWishart(outputs + 5, psi)
        };
    }

    internal static BassData Data(int outputs)
    {
        const int n = 40;
        var rng = new Xoshiro256PlusPlus(20261128);
        double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        double[] z = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        var y = new double[n, outputs];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < outputs; l++) y[i, l] = (l + 1) * (0.5 + Math.Sin(2 * x[i]) + z[i] * x[i]) + 0.15 * Normal.SampleStandard(rng);
        return new BassData([new ContinuousPredictor("x", x), new ContinuousPredictor("z", z)],
                            new MultiOutputResponse(y, [.. Enumerable.Range(0, outputs).Select(l => $"y{l}")]));
    }

    internal static BassTerm[] Terms(string form) => form == "modified"
        ? [new(legendre: [new LegendreFactor(0, 1)]), new(legendre: [new LegendreFactor(1, 2)]),
           new(legendre: [new LegendreFactor(0, 1), new LegendreFactor(1, 1)])]
        : [new([new HingeFactor(0, 1, 5)]), new([new HingeFactor(1, -1, 12)]),
           new([new HingeFactor(0, 1, 8), new HingeFactor(1, 1, 6)])];

    internal static MultiOutputCase DenseCase(BassModel model, BassDesign design, double tau, double b)
    {
        var response = (MultiOutputResponse)model.Data.Response;
        int n = response.Cases, k = response.Outputs, m = design.Count;
        bool flat = model.Target.Intercept == GPriorCentre.Flat;
        int offset = flat ? 0 : 1, p = m + offset;
        double[,] xc = MultiOutputCoreTests.Columns(model, design);
        var means = new double[m];
        var raw = new double[n];
        for (int j = 0; j < m; j++)
        {
            model.Basis.CaseColumn(design.Terms[j], model.Inputs.Cases, raw);
            means[j] = raw.Average();
        }
        double[,] gram = ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(xc), xc);
        var x = new double[n, p];
        var prior = new double[p, p];
        var centre = new double[p, k];
        var shifts = new double[m];
        for (int j = 0; j < m; j++)
        {
            double scale = model.Config.Coefficients.Form is ModifiedGForm modified ? 1.0 / modified.Factor(design.Terms[j]) : 1.0;
            shifts[j] = model.Config.Coefficients.Form is RawRidgeForm ? -means[j] : (scale - 1) * means[j];
            for (int q = 0; q < m; q++)
            {
                double other = model.Config.Coefficients.Form is ModifiedGForm g ? 1.0 / g.Factor(design.Terms[q]) : 1.0;
                double d = model.Config.Coefficients.Form switch {
                    RidgeForm => j == q ? gram[j, j] / n : 0,
                    RawRidgeForm => j == q ? 1 : 0,
                    _ => scale * other * gram[j, q]
                };
                prior[j + offset, q + offset] = d;
            }
        }
        if (!flat)
        {
            double mass = model.Config.Coefficients.Form is ModifiedGForm ? n : 1;
            prior[0, 0] = mass;
            for (int j = 0; j < m; j++)
            {
                prior[0, j + 1] = prior[j + 1, 0] = mass * shifts[j];
                for (int q = 0; q < m; q++) prior[j + 1, q + 1] += mass * shifts[j] * shifts[q];
            }
            for (int l = 0; l < k; l++) centre[0, l] = model.Target.Intercept == GPriorCentre.Mean ? model.Resolved.Outputs!.Means[l] : 0;
        }
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++) prior[i, j] *= tau;
        var y = new double[n, k];
        for (int i = 0; i < n; i++)
        {
            if (!flat) x[i, 0] = 1;
            for (int j = 0; j < m; j++) x[i, j + offset] = xc[i, j];
            for (int l = 0; l < k; l++) y[i, l] = response.Values[i * k + l];
        }
        return new MultiOutputCase(y, x, prior, centre, b, flat);
    }

    [Theory]
    [InlineData("ridge")]
    [InlineData("raw")]
    [InlineData("modified")]
    public void EveryCoefficientFormMatchesTheDenseMatrixTAndScalarReduction(string form)
    {
        double worst = 0;
        int compared = 0;
        foreach (var centre in new[] { GPriorCentre.Zero, GPriorCentre.Mean, GPriorCentre.Flat })
            foreach (int k in new[] { 1, 2, 3 })
            {
                var data = Data(k);
                var model = BassModel.Create(Config(form, centre, k), data);
                var design = model.NewDesign();
                Assert.Null(design.Load(Terms(form)));
                var fit = BassStructureFit.Of(design, design.Terms.Sum(model.Prior.LogPrior));
                var prior = model.Resolved.Outputs!.Prior;
                var reference = new ReferenceSigmaPrior(prior.DegreesOfFreedom, prior.Scale);
                var evaluation = model.Target.Evaluate(fit, Math.Log(0.03));
                const double h = 1e-5, temperature = 0.4;
                double derivative = (evaluation.LogLikelihood(temperature + h) - evaluation.LogLikelihood(temperature - h)) / (2 * h);
                Assert.InRange(Math.Abs(derivative - evaluation.LogLikelihoodSlope(temperature)), 0, 2e-7);
                foreach (double tau in new[] { 1e-8, 1e-2, 1.0, 1e4 })
                    foreach (double b in new[] { 1e-6, 0.2, 1.0 })
                    {
                        double expected = MatrixTReference.LogMarginal(MatrixTReference.Collapse(DenseCase(model, design, tau, b)), reference);
                        double actual = model.Target.Evaluate(fit, Math.Log(tau)).LogLikelihood(b);
                        worst = Math.Max(worst, Math.Abs(actual - expected) / Math.Max(1, Math.Abs(expected)));
                        compared++;
                    }
                if (k == 1)
                {
                    double[] y = [.. data.Response.Values];
                    var scalar = BassModel.Create(model.Config with { OutputCovariance = null,
                        Sigma2 = Sigma2Prior.InverseGamma(prior.DegreesOfFreedom / 2, prior.Scale[0, 0] / 2) },
                        new BassData(data.Predictors, new ScalarResponse(y)));
                    var one = scalar.NewDesign();
                    Assert.Null(one.Load(Terms(form)));
                    var singleFit = BassStructureFit.Of(one, one.Terms.Sum(scalar.Prior.LogPrior));
                    foreach (double tau in new[] { 1e-5, 0.03, 5.0 })
                    {
                        double a = model.Target.Evaluate(fit, Math.Log(tau)).LogMarginal;
                        double d = scalar.Target.Evaluate(singleFit, Math.Log(tau)).LogMarginal;
                        Assert.InRange(Math.Abs(a - d) / Math.Max(1, Math.Abs(d)), 0, 2e-10);
                    }
                }
            }
        output.WriteLine($"{form}: {compared} dense cases, maximum relative marginal error {worst:G4}");
        Assert.True(worst < 2e-9);
    }

    [Theory]
    [InlineData("ridge")]
    [InlineData("raw")]
    [InlineData("modified")]
    public void ProposedCoefficientBlocksMatchFreshDenseFits(string form)
    {
        var model = BassModel.Create(Config(form), Data(2));
        var rng = new Xoshiro256PlusPlus(20261208);
        var prior = new BassPriorSampler(model);
        var iw = model.Resolved.Outputs!.Prior;
        var reference = new ReferenceSigmaPrior(iw.DegreesOfFreedom, iw.Scale);
        var fresh = model.NewDesign();
        int compared = 0, reused = 0;
        double worst = 0;
        for (int i = 0; i < 200; i++)
        {
            BassTerm[] terms = i % 2 == 0 ? Terms(form) : [.. prior.DrawState(rng).State.Terms];
            var moves = new BassStructureMoves(model, terms, Math.Log(0.03));
            foreach (var move in moves.Moves)
            {
                if (move.Weight(moves.Start) == 0 || move.Propose(moves.Start, rng) is not { } proposal) continue;
                var state = proposal.Candidate;
                var refusal = fresh.Load(state.Terms);
                Assert.Equal(refusal is not null, model.Target.Evaluate(state).OffSupport);
                if (refusal is not null) continue;
                double expected = MatrixTReference.LogMarginal(MatrixTReference.Collapse(DenseCase(model, fresh, 0.03, 0.4)), reference);
                double actual = model.Target.Evaluate(state).LogLikelihood(0.4);
                worst = Math.Max(worst, Math.Abs(expected - actual) / Math.Max(1, Math.Abs(expected)));
                bool pending = moves.Holds(state) && !state.SharesStructure(moves.Start);
                var accepted = moves.DesignAt(state);
                if (pending)
                {
                    Assert.Same(state.Fit.Outputs!.Prior, accepted.OutputForms!.Prior);
                    reused++;
                }
                compared++;
            }
        }
        output.WriteLine($"{form}: {compared} proposal blocks, {reused} accepted pencils reused, maximum relative matrix-t error {worst:G4}");
        Assert.True(compared > 100 && reused > 100 && worst < 2e-9);
    }

    [Theory]
    [InlineData("ridge", GPriorCentre.Zero)]
    [InlineData("raw", GPriorCentre.Mean)]
    [InlineData("modified", GPriorCentre.Flat)]
    public void CoefficientDrawsHaveTheDenseConditionalMatrixNormalLaw(string form, GPriorCentre centre)
    {
        var model = BassModel.Create(Config(form, centre), Data(2));
        var design = model.NewDesign();
        Assert.Null(design.Load(Terms(form)));
        var fit = BassStructureFit.Of(design, design.Terms.Sum(model.Prior.LogPrior));
        const double tau = 0.3;
        var c = DenseCase(model, design, tau, 1);
        var precision = ReferenceAlgebra.Sum(c.PriorPrecision,
            ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(c.X), c.X));
        var pRoot = ReferenceAlgebra.Cholesky(precision)!;
        var rhs = ReferenceAlgebra.Sum(ReferenceAlgebra.Multiply(ReferenceAlgebra.Transpose(c.X), c.Y),
                                      ReferenceAlgebra.Multiply(c.PriorPrecision, c.PriorMean!));
        var mean = ReferenceAlgebra.Solve(pRoot, rhs);
        int p = c.Columns, offset = centre == GPriorCentre.Flat ? 1 : 0, dimension = 2 * p;
        var first = new RunningMoments[dimension];
        var second = new RunningMoments[dimension];
        var cross = new RunningMoments();
        var rng = new Xoshiro256PlusPlus(20261129);
        for (int sample = 0; sample < 20_000; sample++)
        {
            var (sigma, coefficients) = BassOutputDraws.Sample(model.Target, design, fit, tau, rng);
            double l11 = Math.Sqrt(sigma[0, 0]), l21 = sigma[1, 0] / l11, l22 = Math.Sqrt(sigma[1, 1] - l21 * l21);
            var whitened = new double[p, 2];
            for (int j = 0; j < p; j++)
            {
                whitened[j, 0] = (coefficients[j + offset, 0] - mean[j, 0]) / l11;
                whitened[j, 1] = (coefficients[j + offset, 1] - mean[j, 1] - l21 * whitened[j, 0]) / l22;
            }
            var z = new double[dimension];
            for (int j = 0; j < p; j++)
                for (int l = 0; l < 2; l++)
                {
                    double value = 0;
                    for (int r = j; r < p; r++) value += pRoot[r, j] * whitened[r, l];
                    z[2 * j + l] = value;
                    first[2 * j + l].Add(value);
                    second[2 * j + l].Add(value * value);
                }
            cross.Add(z[0] * z[^1]);
        }
        double worst = 0;
        foreach (var moments in first) worst = Math.Max(worst, Math.Abs(moments.Mean) / Math.Sqrt(moments.Variance / moments.Count));
        foreach (var moments in second) worst = Math.Max(worst, Math.Abs(moments.Mean - 1) / Math.Sqrt(moments.Variance / moments.Count));
        worst = Math.Max(worst, Math.Abs(cross.Mean) / Math.Sqrt(cross.Variance / cross.Count));
        output.WriteLine($"{form}, {centre}: {dimension} standardized coefficient coordinates, worst moment |z|={worst:F3}");
        Assert.True(worst < 4.5);
    }
}
