using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// Priors that stop reading a moving response (the options note's steps 1 and 2): the frozen reference, the frozen Gram, the
/// coefficient prior's constant and σ²'s hierarchical prior. A design's working response cannot move yet, so a moving response is
/// a model per snapshot here, each chain segment starting where the last ended.
/// <list type="bullet">
/// <item><b>Nothing changes unasked.</b> A reference that restates today's read is today's target, chain for chain.</item>
/// <item><b>Dense algebra.</b> S, the determinant term and the conditional mean against matrices built here from the columns, with
/// the level's energy, for each form and each Gram; and the reductions among the constant, the frozen Gram and τ.</item>
/// <item><b>The target stays one target.</b> Across snapshots a frozen read keeps its priors where a re-read moves them, and the
/// frozen and invariant routes are equivariant under y → ky + m, deterministically and through chains, where a declared prior is
/// not.</item>
/// <item><b>The hierarchical family.</b> Against the half-Cauchy read independently, against both its limits, against enumeration
/// through a chain's visits, with σ² integrated and with σ² drawn apart, and against an evidence at the run's level.</item>
/// </list>
/// </summary>
public sealed class BassPriorReferenceTests
{
    private readonly ITestOutputHelper _out;
    public BassPriorReferenceTests(ITestOutputHelper output) => _out = output;

    private sealed class Snapshot(Func<BassResponse, WorkingResponse> working) : IBassObservationModel
    {
        public WorkingResponse Working(BassResponse response) => working(response);

        public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot) => [];
    }

    private static void Close(double expected, double actual, double relative, string what) => BassWeightedCoreTests.Close(expected, actual, relative, what);

    private static double[] Unit(int n) => [.. Enumerable.Repeat(1.0, n)];

    // A working response's offset at one snapshot: what y − β√w·v is less y, with a mean of its own.
    private static double[] Offset(int n, int snapshot)
    {
        var rng = new Xoshiro256PlusPlus(100 + snapshot);
        double shift = 0.3 * (snapshot + 1);
        return [.. Enumerable.Range(0, n).Select(_ => shift + 0.4 * (rng.NextDouble() - 0.5))];
    }

    private static BassData Rescaled(BassData data, double k, double m)
        => new(data.Predictors, new ScalarResponse([.. data.Response.Values.ToArray().Select(v => k * v + m)]));

    // ---- Nothing changes unasked ----

    [Theory]
    [InlineData("g", GPriorCentre.Mean)]
    [InlineData("g", GPriorCentre.Flat)]
    [InlineData("g", GPriorCentre.Zero)]
    [InlineData("ridge", GPriorCentre.Mean)]
    [InlineData("ridge", GPriorCentre.Zero)]
    public void AReferenceThatRestatesTheRead_IsTodaysTarget_ChainForChain(string form, GPriorCentre intercept)
    {
        BassData data = BassWeightedCoreTests.Scattered(60, 11);
        BassConfig config = BassSweepTests.Config(maxBasis: 6) with { Intercept = intercept, Sigma2 = Sigma2Prior.ScaledToResponse(0.5) };
        if (form == "ridge") config = config with { Coefficients = BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()) };
        BassModel plain = BassModel.Create(config, data);
        BassModel restated = BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, reference:
            intercept == GPriorCentre.Zero ? BassPriorReference.About(r.Values, 0.0) : BassPriorReference.Moments(r.Values))));
        Assert.Equal(plain.Resolved.CentreLevel, restated.Resolved.CentreLevel);
        Assert.Equal(plain.Resolved.ResponseMeanSquare, restated.Resolved.ResponseMeanSquare);
        Assert.Equal(plain.Resolved.Sigma2!.Value.Scale, restated.Resolved.Sigma2!.Value.Scale);
        var a = new BassChain(plain, [], Math.Log(0.05), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        var b = new BassChain(restated, [], Math.Log(0.05), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        for (int s = 0; s < 200; s++)
        {
            BassState x = a.Step(), y = b.Step();
            Assert.Equal(x.Terms, y.Terms);
            Assert.Equal(x.LogTau, y.LogTau);
        }
        Assert.True(a.Current.Count > 0, "the chain should have taken births");
    }

    // ---- Dense algebra ----

    [Theory]
    [InlineData("g", "unit")]
    [InlineData("ridge", "unit")]
    [InlineData("g", "apart")]
    [InlineData("ridge", "apart")]
    [InlineData("g", "weighted")]
    [InlineData("ridge", "weighted")]
    [InlineData("g", "frozen")]
    [InlineData("ridge", "frozen")]
    [InlineData("g", "frozen, unit likelihood")]
    [InlineData("ridge", "frozen, unit likelihood")]
    public void AFrozenReference_MatchesDenseAlgebra(string form, string weighting)
    {
        bool ridge = form == "ridge";
        BassData data = BassSweepTests.Data(large: true);
        double[] observed = data.Response.Values.ToArray();
        int n = observed.Length;
        double[] offset = Offset(n, 2);
        double[] z = [.. observed.Zip(offset, (y, o) => y - o)];
        double[]? w = weighting is "unit" or "frozen, unit likelihood" ? null : BassWeightedCoreTests.Weights(n);
        double[] w0 = BassWeightedCoreTests.Weights(n, seed: 9);
        bool frozen = weighting.StartsWith("frozen");
        BassPriorReference reference = frozen ? BassPriorReference.FrozenGram(observed, w0) : BassPriorReference.Moments(observed);
        WeightedPriorGram gram = frozen ? WeightedPriorGram.Frozen : weighting == "weighted" ? WeightedPriorGram.Weighted : WeightedPriorGram.Unweighted;
        BassConfig config = BassWeightedCoreTests.Config(ridge ? new RidgeForm() : new GForm(), GPriorCentre.Mean, Sigma2Prior.InverseGamma(2.0, 0.3), 0.2);
        BassModel model = BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, offset, w, priorGram: gram, reference: reference)));
        Assert.Equal(reference.Level, model.Resolved.CentreLevel);
        Assert.Equal(reference.Level, model.Resolved.PriorLevel);
        double[] likelihood = w ?? Unit(n);
        double[] priorWeight = gram switch { WeightedPriorGram.Weighted => likelihood, WeightedPriorGram.Frozen => w0, _ => Unit(n) };
        double sumLog = likelihood.Sum(Math.Log);
        double mean = z.Zip(likelihood, (v, k) => v * k).Sum() / likelihood.Sum();
        Assert.True(Math.Abs(mean - reference.Level!.Value) > 0.1, $"the fixture's ȳ {mean} should stand away from the level {reference.Level}");

        foreach (int count in new[] { 0, 1, 3, 5 })
        {
            BassTerm[] terms = BassWeightedCoreTests.Terms[..count];
            BassState state = BassMixing.Loaded(model, terms, Math.Log(0.2), out BassRefusal? refusal) ?? throw new InvalidOperationException($"{refusal}");
            double[][] columns = BassWeightedCoreTests.Columns(model, state.Terms);
            foreach ((double tau, double b) in new[] { (0.2, 1.0), (0.2, 0.3), (0.013, 1.0), (3.0, 0.7) })
            {
                BassWeightedCoreTests.Dense dense = BassWeightedCoreTests.Reference(z, likelihood, columns, priorWeight, reference.Level, GPriorCentre.Mean, ridge, tau, b);
                GPriorMarginal marginal = model.Target.Marginal(state.Fit, tau, b);
                string at = $"M = {count}, τ = {tau}, b = {b}";
                Close(dense.S, marginal.S, 1e-9, $"S, {at}");
                Close(dense.LogDeterminant + 0.5 * b * sumLog, marginal.LogDeterminantTerm, 1e-9, $"the determinant term, {at}");

                // The slope in b carries the level's energy too.
                const double h = 1e-5;
                InverseGamma prior = model.Target.Sigma2!.Value;
                double numeric = (model.Target.Marginal(state.Fit, tau, b + (b == 1.0 ? 0.0 : h)).LogMarginal(prior) - model.Target.Marginal(state.Fit, tau, b - h).LogMarginal(prior)) / (b == 1.0 ? h : 2.0 * h);
                Close(numeric, marginal.LogMarginalSlope(prior), 2e-4, $"the slope, {at}");
            }

            BassDesign design = model.NewDesign();
            Assert.Null(design.Load(terms));
            GPriorConditional conditional = model.Target.CoefficientBlock(design)(0.2, 0.8);
            BassWeightedCoreTests.Dense reading = BassWeightedCoreTests.Reference(z, likelihood, columns, priorWeight, reference.Level, GPriorCentre.Mean, ridge, 0.2, 0.8);
            for (int i = 0; i < n; i += 7)
            {
                double f = conditional.Mean[0];
                for (int j = 0; j < count; j++) f += conditional.Mean[j + 1] * (columns[j][i] - columns[j].Zip(likelihood, (v, k) => v * k).Sum() / likelihood.Sum());
                Close(reading.FittedMean[i], f, 1e-9, $"the mean's fit at row {i}, M = {count}");
            }
        }
    }

    // Centred on ȳ the level's energy is exactly none; at a frozen level it is the zero centre's form about that level.
    [Fact]
    public void TheLevelsEnergy_IsNoneAtTheMean_AndTheZeroCentresFormElsewhere()
    {
        const double rss = 4.2, q = 7.9, mean = 1.3, tau = 0.5, b = 0.6;
        GPriorMarginal plain = CollapsedGPrior.CentredMarginalAt(3, rss, q, 40, mean, GPriorCentre.Mean, tau, b);
        GPriorMarginal atMean = CollapsedGPrior.CentredMarginalAt(3, rss, q, 40, mean, GPriorCentre.Mean, tau, b, level: mean);
        Assert.Equal(plain.S, atMean.S);
        Assert.Equal(plain.LogDeterminantTerm, atMean.LogDeterminantTerm);
        GPriorMarginal zero = CollapsedGPrior.CentredMarginalAt(3, rss, q, 40, mean, GPriorCentre.Zero, tau, b);
        GPriorMarginal atZero = CollapsedGPrior.CentredMarginalAt(3, rss, q, 40, mean, GPriorCentre.Mean, tau, b, level: 0.0);
        Assert.Equal(zero.S, atZero.S);
        Assert.Equal(zero.LogDeterminantTerm, atZero.LogDeterminantTerm);
        Assert.Throws<ArgumentException>(() => CollapsedGPrior.CentredMarginalAt(3, rss, q, 40, mean, GPriorCentre.Zero, tau, b, level: 0.0));
        Assert.Throws<ArgumentException>(() => CollapsedGPrior.CentredMarginalAt(3, rss, q, 40, mean, GPriorCentre.Flat, tau, b, level: 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollapsedGPrior.CentredMarginalAt(3, rss, q, 40, mean, GPriorCentre.Mean, tau, b, level: double.NaN));
    }

    // The collapsed predictive under a frozen read, against the joint normal it comes from: z | σ² ~ N(level·1, σ²(W⁻¹ + X(τD)⁻¹Xᵀ)), X
    // and D built here from the columns alone, each row conditioned on the rest by dense algebra.
    [Theory]
    [InlineData("g", "unit")]
    [InlineData("ridge", "unit")]
    [InlineData("g", "apart")]
    [InlineData("ridge", "apart")]
    [InlineData("g", "frozen")]
    [InlineData("ridge", "frozen")]
    public void TheCollapsedPredictive_MatchesTheJointNormal(string form, string weighting)
    {
        const double sigma2 = 0.4, tau = 0.2;
        bool ridge = form == "ridge", frozen = weighting == "frozen";
        BassData data = BassSweepTests.Data(large: true);
        double[] observed = data.Response.Values.ToArray();
        int n = observed.Length;
        double[] offset = Offset(n, 2);
        double[] z = [.. observed.Zip(offset, (y, o) => y - o)];
        double[]? w = weighting == "unit" ? null : BassWeightedCoreTests.Weights(n);
        double[] w0 = BassWeightedCoreTests.Weights(n, seed: 9);
        BassPriorReference reference = frozen ? BassPriorReference.FrozenGram(observed, w0) : BassPriorReference.Moments(observed);
        BassConfig config = BassWeightedCoreTests.Config(ridge ? new RidgeForm() : new GForm(), GPriorCentre.Mean, Sigma2Prior.Held(sigma2), tau);
        BassModel model = BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, offset, w, priorGram: frozen ? WeightedPriorGram.Frozen : WeightedPriorGram.Unweighted, reference: reference)));
        var table = new BassTermTable();
        BassTerm[] terms = BassWeightedCoreTests.Terms[..3];
        BassDraw draw = new BassChain(model, terms, Math.Log(tau), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2)).Compose(new Xoshiro256PlusPlus(3), table);
        PointwiseLogDensity reading = BassPointwise.LogDensity(model, table, [(IReadOnlyList<BassDraw>)[draw]], UnitPartition.EachObservation(n), PointwiseContent.CollapsedPredictive);

        double[] likelihood = w ?? Unit(n), priorWeight = frozen ? w0 : Unit(n);
        double priorMass = priorWeight.Sum();
        double[][] columns = BassWeightedCoreTests.Columns(model, terms);
        int p = terms.Length + 1;
        var x = new double[n, p];
        for (int i = 0; i < n; i++) x[i, 0] = 1.0;
        for (int j = 1; j < p; j++)
        {
            double centre = columns[j - 1].Zip(priorWeight, (v, k) => v * k).Sum() / priorMass;
            for (int i = 0; i < n; i++) x[i, j] = columns[j - 1][i] - centre;
        }
        var d = new double[p, p];
        d[0, 0] = ridge ? 1.0 : priorMass;
        for (int j = 1; j < p; j++)
            for (int k = 1; k < p; k++)
            {
                double own = 0.0;
                for (int i = 0; i < n; i++) own += priorWeight[i] * x[i, j] * x[i, k];
                d[j, k] = !ridge ? own : j == k ? own / priorMass : 0.0;
            }
        var factorD = new DenseCholesky(p);
        factorD.Decompose(d);
        var covariance = new double[n, n];
        var column = new double[p];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < p; j++) column[j] = x[i, j];
            double[] solved = factorD.Solve(column);
            for (int u = 0; u <= i; u++)
            {
                double q = 0.0;
                for (int j = 0; j < p; j++) q += x[u, j] * solved[j];
                covariance[i, u] = covariance[u, i] = q / tau;
            }
            covariance[i, i] += 1.0 / likelihood[i];
        }
        var factor = new DenseCholesky(n);
        factor.Decompose(covariance);
        Assert.True(factor.Status.PositiveDefinite);
        double[] e = [.. z.Select(v => v - reference.Level!.Value)];
        double[] solvedE = factor.Solve(e);
        var unit = new double[n];
        for (int i = 0; i < n; i++)
        {
            Array.Clear(unit);
            unit[i] = 1.0;
            double precision = factor.Solve(unit)[i];
            double expected = 0.5 * Math.Log(precision / (2.0 * Math.PI * sigma2)) - solvedE[i] * solvedE[i] / (2.0 * sigma2 * precision);
            Close(expected, reading.Values[0, i], 1e-9, $"row {i}");
        }
    }

    // ---- The constant and the frozen Gram ----

    // κ reads τ as τ/κ whatever the form, for the marginal and the coefficients' block alike. τ's own prior does not move.
    [Theory]
    [InlineData("g")]
    [InlineData("ridge")]
    [InlineData("apart")]
    public void TheConstant_ReadsTauAsTauOverC(string form)
    {
        const double c = 4.0;
        BassData data = BassSweepTests.Data(large: true);
        int n = data.Response.Observations;
        double[]? w = form == "apart" ? BassWeightedCoreTests.Weights(n) : null;
        BassConfig config = BassWeightedCoreTests.Config(form == "ridge" ? new RidgeForm() : new GForm(), GPriorCentre.Flat, Sigma2Prior.Held(1.0), 0.2);
        BassModel one = BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, weights: w)));
        BassModel constant = BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, weights: w, reference: new BassPriorReference(coefficientConstant: c))));
        Assert.Equal(c, constant.Target.CoefficientConstant);
        Assert.Equal(1.0, one.Target.CoefficientConstant);
        Assert.Equal(one.Target.Scale, constant.Target.Scale);
        BassTerm[] terms = BassWeightedCoreTests.Terms[..3];
        BassState x = BassMixing.Loaded(one, terms, 0.0, out _)!, y = BassMixing.Loaded(constant, terms, 0.0, out _)!;
        GPriorMarginal scaled = constant.Target.Marginal(y.Fit, 0.2, 0.7), direct = one.Target.Marginal(x.Fit, 0.2 / c, 0.7);
        Assert.Equal(direct.S, scaled.S);
        Assert.Equal(direct.LogDeterminantTerm, scaled.LogDeterminantTerm);
        BassDesign d0 = one.NewDesign(), d1 = constant.NewDesign();
        Assert.Null(d0.Load(terms));
        Assert.Null(d1.Load(terms));
        GPriorConditional a = constant.Target.CoefficientBlock(d1)(0.2, 0.7), b = one.Target.CoefficientBlock(d0)(0.2 / c, 0.7);
        Assert.True(a.Mean.SequenceEqual(b.Mean));
        Assert.Equal(a.CovarianceScale, b.CovarianceScale);

        // A composed draw reads it too: at a vanishing σ² the coefficients are the block's mean.
        BassConfig quiet = config with { Sigma2 = Sigma2Prior.Held(1e-30) };
        BassModel held = BassModel.Create(quiet, data, new Snapshot(r => new WorkingResponse(r, weights: w, reference: new BassPriorReference(coefficientConstant: c))));
        BassDraw drawn = new BassChain(held, terms, Math.Log(0.2), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2)).Compose(new Xoshiro256PlusPlus(3), new BassTermTable());
        BassDraw atOver = new BassChain(BassModel.Create(BassWeightedCoreTests.Config(form == "ridge" ? new RidgeForm() : new GForm(), GPriorCentre.Flat, Sigma2Prior.Held(1e-30), 0.2 / c), data, new Snapshot(r => new WorkingResponse(r, weights: w))),
                                        terms, Math.Log(0.2 / c), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2)).Compose(new Xoshiro256PlusPlus(3), new BassTermTable());
        for (int j = 0; j < drawn.Coefficients.Length; j++) Close(atOver.Coefficients[j], drawn.Coefficients[j], 1e-9, $"a_{j}");
    }

    // Sabanés Bové and Held's constant from a link's variance and slope at 0, their Table 1's values.
    [Fact]
    public void TheLinkConstants_AreTheTables()
    {
        Assert.Equal(4.0, BassPriorReference.Link(0.25, 0.25).CoefficientConstant, 12);
        Assert.Equal(Math.PI / 2.0, BassPriorReference.Link(0.25, 1.0 / Math.Sqrt(2.0 * Math.PI)).CoefficientConstant, 12);
        Assert.Equal(Math.PI * Math.PI / 4.0, BassPriorReference.Link(0.25, 1.0 / Math.PI).CoefficientConstant, 12);
        Assert.Equal(1.0, BassPriorReference.Link(1.0, 1.0).CoefficientConstant, 12);
        Assert.Null(BassPriorReference.Link(0.25, 0.25).Level);
        Assert.Null(BassPriorReference.Link(0.25, 0.25).MeanSquare);
    }

    // The ladder's automatic hot end reads the constant: the b where τ₀/κ halves a unit-information column's shrinkage.
    [Fact]
    public void TheAutoHotEnd_ReadsTheConstant()
    {
        var spec = new TemperingSpec();
        Assert.Equal(1.0 / 100.0, BassResolved.AutoHotEnd(spec, GPriorCentre.Mean, 100.0, null));
        Assert.Equal(1.0 / 400.0, BassResolved.AutoHotEnd(spec, GPriorCentre.Mean, 100.0, null, 4.0));
        Assert.Equal(2.0 / 100.0, BassResolved.AutoHotEnd(spec, GPriorCentre.Flat, 100.0, null, 4.0));   // the pole stays
        ResolvedTauPrior held = TauPrior.Held(0.002).For(100.0);
        Assert.Equal(0.002 / 4.0, BassResolved.HalvingTau(new GForm(), held, 100.0, 4.0));
        Assert.Equal(0.002 / 400.0, BassResolved.HalvingTau(new RidgeForm(), held, 100.0, 4.0));
        BassData data = BassSweepTests.Data(large: true);
        BassConfig config = BassSweepTests.Config(maxBasis: 3) with
        {
            Sigma2 = Sigma2Prior.Held(1.0),
            Sampling = new BassSampling { Moves = new EqualOverAvailable(), Tempering = new TemperingSpec { Levels = 4 } },
        };
        BassModel model = BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, reference: new BassPriorReference(level: 0.0, coefficientConstant: 4.0))));
        Assert.Equal(1.0 / (4.0 * data.Response.Observations), model.Resolved.HotEnd!.Value, 15);
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        var ladder = new BassLadder(model, config.Sampling.Tempering!, start.Terms, start.LogTau, 2);
        Assert.Equal(model.Resolved.HotEnd!.Value, ladder.Tempering.InverseTemperatures.Min(), 12);
    }

    // W₀ reads in information per observation: under the g-form W₀ = (1/κ)·1 is the constant κ on the unweighted columns, Li and Clyde's
    // information at β = 0 for a link. A unit W₀ is the unweighted prior under both forms, and W₀ = W the weighted prior at a held τ.
    [Theory]
    [InlineData("g")]
    [InlineData("ridge")]
    public void TheFrozenGram_ReducesToTheOtherPriors(string form)
    {
        bool ridge = form == "ridge";
        BassData data = BassSweepTests.Data(large: true);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        double[] w = BassWeightedCoreTests.Weights(n);
        BassConfig config = BassWeightedCoreTests.Config(ridge ? new RidgeForm() : new GForm(), GPriorCentre.Mean, Sigma2Prior.InverseGamma(2.0, 0.3), 0.2);
        BassTerm[] terms = BassWeightedCoreTests.Terms[..4];
        BassModel Built(WeightedPriorGram gram, BassPriorReference? reference) => BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, weights: w, priorGram: gram, reference: reference)));
        void Same(BassModel a, BassModel b, string what)
        {
            BassState x = BassMixing.Loaded(a, terms, 0.0, out _)!, u = BassMixing.Loaded(b, terms, 0.0, out _)!;
            foreach (double t in new[] { 1.0, 0.4 })
            {
                GPriorMarginal p = a.Target.Marginal(x.Fit, 0.2, t), q = b.Target.Marginal(u.Fit, 0.2, t);
                Close(p.S, q.S, 1e-9, $"{what}: S at b = {t}");
                Close(p.LogDeterminantTerm, q.LogDeterminantTerm, 1e-9, $"{what}: the determinant term at b = {t}");
            }
        }
        BassPriorReference moments = BassPriorReference.Moments(y);
        BassModel unweighted = Built(WeightedPriorGram.Unweighted, moments);
        Same(Built(WeightedPriorGram.Frozen, BassPriorReference.FrozenGram(y, Unit(n))), unweighted, "a unit W₀ against the unweighted prior");
        Same(Built(WeightedPriorGram.Frozen, BassPriorReference.FrozenGram(y, w)), Built(WeightedPriorGram.Weighted, new BassPriorReference(y.Zip(w, (v, k) => v * k).Sum() / w.Sum(), 1.0)), "W₀ = W against the weighted prior");
        if (ridge) return;
        const double c = 4.0;
        double[] oneOverC = [.. Enumerable.Repeat(1.0 / c, n)];
        BassModel constant = Built(WeightedPriorGram.Unweighted, BassPriorReference.Moments(y, c));
        BassModel information = Built(WeightedPriorGram.Frozen, BassPriorReference.FrozenGram(y, oneOverC));
        Same(information, constant, "W₀ = (1/κ)·1 against the constant κ");

        // Unit-information τ on the two sides: both stated on N, so the same law for the same covariance.
        Assert.Equal(unweighted.Resolved.Scale.UnitObservations, information.Resolved.Scale.UnitObservations);
    }

    // A design holding the frozen Gram carries it through proposals as through a load, and hands it across a ladder's replicas.
    [Theory]
    [InlineData("g")]
    [InlineData("ridge")]
    public void TheFrozenGram_CarriesThroughAProposalAsThroughALoad(string form)
    {
        BassData data = BassSweepTests.Data(large: true);
        double[] y = data.Response.Values.ToArray();
        double[] w = BassWeightedCoreTests.Weights(y.Length), w0 = BassWeightedCoreTests.Weights(y.Length, seed: 9);
        BassConfig config = BassWeightedCoreTests.Config(form == "ridge" ? new RidgeForm() : new GForm(), GPriorCentre.Mean, Sigma2Prior.InverseGamma(2.0, 0.3), 0.2);
        BassModel model = BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, weights: w, priorGram: WeightedPriorGram.Frozen, reference: BassPriorReference.FrozenGram(r.Values, w0))));
        BassDesign built = model.NewDesign();
        foreach (BassTerm term in BassWeightedCoreTests.Terms[..4])
        {
            Assert.True(built.ProposeBirth(term).Admitted);
            built.Accept();
        }
        BassDesign loaded = model.NewDesign();
        Assert.Null(loaded.Load(built.Terms));
        GPriorMarginal a = model.Target.Marginal(BassStructureFit.Of(built, 0.0), 0.2, 0.6), b = model.Target.Marginal(BassStructureFit.Of(loaded, 0.0), 0.2, 0.6);
        Close(a.S, b.S, 1e-10, "S");
        Close(a.LogDeterminantTerm, b.LogDeterminantTerm, 1e-10, "the determinant term");

        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        var ladder = new BassLadder(model, new TemperingSpec { Levels = 4, HotEnd = 0.1 }, start.Terms, start.LogTau, 2);
        for (int s = 0; s < 30; s++) ladder.Step();
        BassState cold = ladder.Current;
        BassState fresh = BassMixing.Loaded(model, cold.Terms, cold.LogTau, out BassRefusal? refusal)!;
        Assert.Null(refusal);
        Close(fresh.Fit.ResidualSquares, cold.Fit.ResidualSquares, 1e-9, "the ladder's residual");
    }

    // ---- The rules ----

    [Fact]
    public void TheRules_ReadWhatTheyName_AndAreEquivariant()
    {
        var rng = new Xoshiro256PlusPlus(5);
        double[] y = [.. Enumerable.Range(0, 201).Select(_ => 2.0 + BassQualification.Gaussian(rng))];
        const double k = 37.0, m = -11.0;
        double[] moved = [.. y.Select(v => k * v + m)];
        double[] w0 = BassWeightedCoreTests.Weights(y.Length, seed: 4);
        (string Name, Func<double[], BassPriorReference> Rule)[] rules =
        [
            ("moments", v => BassPriorReference.Moments(v)),
            ("robust", v => BassPriorReference.Robust(v)),
            ("quantile", v => BassPriorReference.Quantile(v, 0.3)),
            ("frozen", v => BassPriorReference.FrozenGram(v, w0)),
        ];
        foreach ((string name, Func<double[], BassPriorReference> rule) in rules)
        {
            BassPriorReference a = rule(y), b = rule(moved);
            Close(k * a.Level!.Value + m, b.Level!.Value, 1e-12, $"{name}: the level");
            Close(k * k * a.MeanSquare!.Value, b.MeanSquare!.Value, 1e-12, $"{name}: the mean square");
        }

        // Moments, against the sums.
        double mean = y.Average();
        Close(mean, BassPriorReference.Moments(y).Level!.Value, 1e-13, "the mean");
        Close(y.Average(v => (v - mean) * (v - mean)), BassPriorReference.Moments(y).MeanSquare!.Value, 1e-13, "the mean square");
        Close(y.Average(v => v * v), BassPriorReference.About(y, 0.0).MeanSquare!.Value, 1e-13, "the mean square about 0");

        // The median and MAD are the outlier's no more than one rank's: the mean moves with it.
        double[] outlier = [.. y];
        outlier[0] = 1e6;
        BassPriorReference robust = BassPriorReference.Robust(y), hit = BassPriorReference.Robust(outlier);
        Assert.True(Math.Abs(robust.Level!.Value - hit.Level!.Value) < 0.1, "the median stays near");
        Assert.True(Math.Abs(robust.MeanSquare!.Value - hit.MeanSquare!.Value) < 0.5, "the MAD stays near");
        Assert.True(BassPriorReference.Moments(outlier).MeanSquare > 1e6, "the mean square does not");
        Close(1.0, BassPriorReference.Robust(y).MeanSquare!.Value, 0.4, "the normal-consistent scale of unit noise");

        // When over half the values tie the MAD is 0, and the mean square about the median stands in.
        double[] tied = [1, 1, 1, 1, 1, 2, 3];
        Assert.Equal(1.0, BassPriorReference.Robust(tied).Level);
        Close(5.0 / 7.0, BassPriorReference.Robust(tied).MeanSquare!.Value, 1e-13, "the tied scale");

        // The sample quantile minimizes the check loss over the sample, and the scale is the asymmetric Laplace's maximum likelihood.
        foreach (double q in new[] { 0.1, 0.3, 0.5, 0.9 })
        {
            double Loss(double level) => y.Average(v => v < level ? (q - 1.0) * (v - level) : q * (v - level));
            BassPriorReference quantile = BassPriorReference.Quantile(y, q);
            double level = quantile.Level!.Value, best = y.Min(Loss);
            Assert.Contains(level, y);
            Close(best, Loss(level), 1e-13, $"the check loss at q = {q}");
            double scale = Math.Sqrt(quantile.MeanSquare!.Value);
            double LogLikelihood(double s) => y.Sum(v => Math.Log(q * (1.0 - q) / s) - (v < level ? (q - 1.0) * (v - level) : q * (v - level)) / s);
            Assert.True(LogLikelihood(scale) >= LogLikelihood(scale * 1.01) && LogLikelihood(scale) >= LogLikelihood(scale * 0.99), $"the scale at q = {q}");
        }

        Assert.Throws<ArgumentOutOfRangeException>(() => BassPriorReference.Quantile(y, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => BassPriorReference.Quantile(y, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new BassPriorReference(double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => new BassPriorReference(0.0, -1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new BassPriorReference(coefficientConstant: 0.0));
        Assert.Throws<ArgumentException>(() => new BassPriorReference(gramWeights: [1.0, 0.0]));
        Assert.Throws<ArgumentException>(() => BassPriorReference.FrozenGram(y, [1.0]));
    }

    // ---- The contract's refusals ----

    [Fact]
    public void TheRefusals_ReadTheReference()
    {
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        double[] toConstant = [.. y.Select(v => v - 3.0)];   // y − offset = 3 everywhere
        BassConfig scaled = BassSweepTests.Config() with { Sigma2 = Sigma2Prior.ScaledToResponse(0.5) };
        Snapshot With(BassPriorReference? reference, WeightedPriorGram gram = WeightedPriorGram.Unweighted) => new(r => new WorkingResponse(r, toConstant, priorGram: gram, reference: reference));

        // A working response that happens to be constant is refused only when the priors read it.
        Assert.Throws<ArgumentException>(() => BassModel.Create(scaled, data, With(null)));
        BassModel admitted = BassModel.Create(scaled, data, With(BassPriorReference.Moments(y)));
        Assert.Equal(y.Average(), admitted.Resolved.CentreLevel, 12);
        Assert.Null(admitted.NewDesign().Load([]));

        // Improper σ² still has no posterior on it, reference or not.
        Assert.Throws<ArgumentException>(() => BassModel.Create(scaled with { Sigma2 = Sigma2Prior.Jeffreys }, data, With(BassPriorReference.Moments(y))));

        // The reference's own spread is the gate; a held σ² lifts it.
        BassPriorReference spreadless = new(1.0, 0.0);
        Assert.Throws<ArgumentException>(() => BassModel.Create(scaled, data, With(spreadless)));
        Assert.Throws<ArgumentException>(() => BassModel.Create(scaled with { Intercept = GPriorCentre.Zero }, data, new Snapshot(r => new WorkingResponse(r, reference: new BassPriorReference(0.0, 0.0)))));
        Assert.NotNull(BassModel.Create(scaled with { Sigma2 = Sigma2Prior.Held(1.0) }, data, With(spreadless)));

        // What the config reads the reference must state; what it does not read, it need not.
        Assert.Throws<ArgumentException>(() => BassModel.Create(scaled, data, With(new BassPriorReference(meanSquare: 1.0))));
        Assert.Throws<ArgumentException>(() => BassModel.Create(scaled, data, new Snapshot(r => new WorkingResponse(r, reference: new BassPriorReference(level: 1.0)))));
        Assert.Throws<ArgumentException>(() => BassModel.Create(scaled with { Sigma2 = Sigma2Prior.HalfCauchy() }, data, new Snapshot(r => new WorkingResponse(r, reference: new BassPriorReference(level: 1.0)))));
        Assert.NotNull(BassModel.Create(scaled with { Intercept = GPriorCentre.Flat, Sigma2 = Sigma2Prior.Held(1.0) }, data, new Snapshot(r => new WorkingResponse(r, reference: BassPriorReference.Link(0.25, 0.25)))));

        // The zero centre's level is 0, and a reference that says otherwise contradicts it.
        Assert.Throws<ArgumentException>(() => BassModel.Create(scaled with { Intercept = GPriorCentre.Zero }, data, new Snapshot(r => new WorkingResponse(r, reference: new BassPriorReference(1.0, 1.0)))));
        Assert.NotNull(BassModel.Create(scaled with { Intercept = GPriorCentre.Zero }, data, new Snapshot(r => new WorkingResponse(r, reference: BassPriorReference.About(r.Values, 0.0)))));

        // The frozen Gram and its weights come together, at the response's length.
        Assert.Throws<ArgumentException>(() => new WorkingResponse(data.Response, priorGram: WeightedPriorGram.Frozen));
        Assert.Throws<ArgumentException>(() => new WorkingResponse(data.Response, priorGram: WeightedPriorGram.Frozen, reference: BassPriorReference.Moments(y)));
        Assert.Throws<ArgumentException>(() => new WorkingResponse(data.Response, reference: BassPriorReference.FrozenGram(y, Unit(n))));
        Assert.Throws<ArgumentException>(() => new WorkingResponse(data.Response, priorGram: WeightedPriorGram.Frozen, reference: new BassPriorReference(gramWeights: [1.0, 1.0])));
        var unitLikelihood = new WorkingResponse(data.Response, priorGram: WeightedPriorGram.Frozen, reference: BassPriorReference.FrozenGram(y, Unit(n)));
        Assert.True(unitLikelihood.PriorStandsApart);
        Assert.Equal(new GPriorWeights(n, n, 0.0), unitLikelihood.Weighting);
    }

    // The exact-fit floor reads the reference's spread, so F does not move with a working response.
    [Fact]
    public void TheExactFitFloor_ReadsTheReference()
    {
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        const double s2 = 9.0;
        BassConfig config = BassSweepTests.Config() with { Intercept = GPriorCentre.Flat, Sigma2 = Sigma2Prior.Jeffreys };
        BassModel Built(double[] offset) => BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, offset, reference: new BassPriorReference(meanSquare: s2))));
        double floor = 1e3 * 2.220446049250313e-16;
        double exact = floor * floor * y.Length * s2;
        foreach (double[] offset in new[] { Offset(y.Length, 0), Offset(y.Length, 3) })
        {
            BassTarget target = Built(offset).Target;
            Assert.Equal(BassRefusal.ExactFit, target.Refusal(new BassStructureFit(0, 0.0, 0.5 * exact, 0.0)));
            Assert.Null(target.Refusal(new BassStructureFit(0, 0.0, 2.0 * exact, 0.0)));
        }
    }

    // ---- The target stays one target ----

    // Across snapshots a frozen read's priors do not move, and a re-read's do.
    [Fact]
    public void AcrossSnapshots_AFrozenReadHoldsThePriors_WhereARereadMovesThem()
    {
        BassData data = BassSweepTests.Data(large: true);
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        BassConfig config = BassSweepTests.Config() with { Sigma2 = Sigma2Prior.ScaledToResponse(0.5) };
        BassPriorReference reference = BassPriorReference.Moments(y);
        var frozen = new List<BassResolved>();
        var reread = new List<BassResolved>();
        for (int t = 0; t < 3; t++)
        {
            double[] offset = Offset(n, t);
            frozen.Add(BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, offset, reference: reference))).Resolved);
            reread.Add(BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, offset))).Resolved);
        }
        foreach (BassResolved r in frozen)
        {
            Assert.Equal(reference.Level!.Value, r.CentreLevel);
            Assert.Equal(reference.MeanSquare!.Value, r.ResponseMeanSquare);
            Assert.Equal(0.5 * reference.MeanSquare.Value, r.Sigma2!.Value.Scale);
        }
        Assert.NotEqual(reread[0].CentreLevel, reread[1].CentreLevel);
        Assert.NotEqual(reread[0].ResponseMeanSquare, reread[1].ResponseMeanSquare);
        Assert.NotEqual(frozen[0].ResponseMean, frozen[1].ResponseMean);   // ȳ is a statistic, and moves
    }

    // The routes, stated on one response: the config and the rule.
    private static (BassConfig Config, Func<double[], BassPriorReference?> Rule, WeightedPriorGram Gram, double[]? W0) Route(string route, bool ridge, int n)
    {
        BassConfig Of(GPriorCentre intercept, Sigma2Prior sigma2) => BassWeightedCoreTests.Config(ridge ? new RidgeForm() : new GForm(), intercept, sigma2, 0.2, 6);
        double[] w0 = BassWeightedCoreTests.Weights(n, seed: 9);
        return route switch
        {
            "moments" => (Of(GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5)), v => BassPriorReference.Moments(v), WeightedPriorGram.Unweighted, null),
            "robust" => (Of(GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5)), v => BassPriorReference.Robust(v), WeightedPriorGram.Unweighted, null),
            "quantile" => (Of(GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5)), v => BassPriorReference.Quantile(v, 0.3), WeightedPriorGram.Unweighted, null),
            "frozen" => (Of(GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5)), v => BassPriorReference.FrozenGram(v, w0), WeightedPriorGram.Frozen, w0),
            "invariant" => (Of(GPriorCentre.Flat, Sigma2Prior.Jeffreys), _ => null, WeightedPriorGram.Unweighted, null),
            "hierarchical" => (Of(GPriorCentre.Mean, Sigma2Prior.HalfCauchy()), v => BassPriorReference.Moments(v), WeightedPriorGram.Unweighted, null),
            "declared" => (Of(GPriorCentre.Mean, Sigma2Prior.InverseGamma(2.0, 0.3)), v => BassPriorReference.Moments(v), WeightedPriorGram.Unweighted, null),
            _ => throw new ArgumentOutOfRangeException(nameof(route)),
        };
    }

    // y → ky + m, offsets k·offset: the frozen, invariant and hierarchical routes' log posterior differences do not move, whatever
    // the snapshot, on each form and through weights; a declared prior's do. The trilemma, on the targets themselves.
    [Theory]
    [InlineData("moments", false, true)]
    [InlineData("moments", true, true)]
    [InlineData("robust", false, true)]
    [InlineData("quantile", false, true)]
    [InlineData("frozen", false, true)]
    [InlineData("frozen", true, true)]
    [InlineData("invariant", false, true)]
    [InlineData("invariant", true, true)]
    [InlineData("hierarchical", false, true)]
    [InlineData("hierarchical", true, true)]
    [InlineData("declared", false, false)]
    public void TheRoutes_AreEquivariant_AsTheTrilemmaAllows(string route, bool ridge, bool equivariant)
    {
        const double k = 250.0, m = -40.0;
        BassData data = BassSweepTests.Data(large: true), moved = Rescaled(data, k, m);
        int n = data.Response.Observations;
        (BassConfig config, Func<double[], BassPriorReference?> rule, WeightedPriorGram gram, _) = Route(route, ridge, n);
        double[] w = BassWeightedCoreTests.Weights(n);
        double largest = 0.0;
        for (int t = 0; t < 3; t++)
        {
            double[] offset = Offset(n, t), scaled = [.. offset.Select(o => k * o)];
            BassModel a = BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, offset, w, priorGram: gram, reference: rule(r.Values.ToArray()))));
            BassModel b = BassModel.Create(config, moved, new Snapshot(r => new WorkingResponse(r, scaled, w, priorGram: gram, reference: rule(r.Values.ToArray()))));
            foreach (double v in new[] { -0.5, 0.0, 1.5 })
            {
                double Gap(BassModel model)
                {
                    double Value(int count)
                    {
                        BassState state = BassMixing.Loaded(model, BassWeightedCoreTests.Terms[..count], 0.0, out BassRefusal? refusal) ?? throw new InvalidOperationException($"{refusal}");
                        if (model.Target.Sigma2Rate is { } rate) state = state.WithLogRate(rate.CentreLogRate + v);
                        return model.Target.Evaluate(state).LogMarginal;
                    }
                    return Value(1) - Value(4);
                }
                double x = Gap(a), y = Gap(b);
                if (equivariant) Close(x, y, 1e-8, $"{route}, snapshot {t}");
                largest = Math.Max(largest, Math.Abs(x - y));
            }
        }
        if (!equivariant) Assert.True(largest > 1e-2, $"A declared prior should not be equivariant: the gap moved {largest}.");
    }

    // A chain whose working response moves: one segment per snapshot, the streams and the state carried across. On y and on ky + m
    // the chains visit the same structures, step by step.
    [Theory]
    [InlineData("moments", false)]
    [InlineData("moments", true)]
    [InlineData("robust", false)]
    [InlineData("frozen", false)]
    [InlineData("invariant", false)]
    [InlineData("hierarchical", false)]
    [InlineData("hierarchical", true)]
    public void AMovingResponse_WalksTheSameStructures_OnKYPlusM(string route, bool ridge)
    {
        const double k = 250.0, m = -40.0;
        BassData data = BassWeightedCoreTests.Scattered(60, 11), moved = Rescaled(data, k, m);
        int n = data.Response.Observations;
        (BassConfig config, Func<double[], BassPriorReference?> rule, WeightedPriorGram gram, _) = Route(route, ridge, n);
        BassCoefficientPrior sampled = new(config.Coefficients.Form, TauPrior.UnitInformationGamma());
        config = config with { Coefficients = sampled };
        var rngA = (new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        var rngB = (new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        IReadOnlyList<BassTerm> termsA = [], termsB = [];
        double tauA = Math.Log(0.05), tauB = tauA, rateA = double.NaN, rateB = double.NaN;
        int changes = 0;
        BassTerm[] last = [];
        for (int t = 0; t < 6; t++)
        {
            double[] offset = Offset(n, t), scaled = [.. offset.Select(o => k * o)];
            double[] w = BassWeightedCoreTests.Weights(n, seed: 30 + t);
            BassModel a = BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, offset, w, priorGram: gram, reference: rule(r.Values.ToArray()))));
            BassModel b = BassModel.Create(config, moved, new Snapshot(r => new WorkingResponse(r, scaled, w, priorGram: gram, reference: rule(r.Values.ToArray()))));
            var chainA = new BassChain(a, termsA, tauA, rngA.Item1, rngA.Item2, logRate: rateA);
            var chainB = new BassChain(b, termsB, tauB, rngB.Item1, rngB.Item2, logRate: rateB);
            for (int s = 0; s < 40; s++)
            {
                BassState x = chainA.Step(), y = chainB.Step();
                Assert.Equal(x.Terms, y.Terms);
                Close(x.LogTau, y.LogTau, 1e-6, $"snapshot {t}, sweep {s}: ln τ");
                if (a.Target.Sigma2Rate is not null) Close(x.LogRate + 2.0 * Math.Log(k), y.LogRate, 1e-6, $"snapshot {t}, sweep {s}: ln r");
                if (!x.Terms.SequenceEqual(last)) changes++;
                last = [.. x.Terms];
            }
            termsA = chainA.Current.Terms;
            termsB = chainB.Current.Terms;
            tauA = chainA.Current.LogTau;
            tauB = chainB.Current.LogTau;
            rateA = chainA.Current.LogRate;
            rateB = chainB.Current.LogRate;
        }
        Assert.True(changes > 15, $"The chain changed structure {changes} times.");
    }

    // ---- Enumeration ----

    // Detailed balance and stationarity on every structure of F under a frozen level and frozen Gram, with the posterior against
    // dense algebra and a chain's visits.
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Enumeration_AFrozenReference(bool frozenGram)
    {
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        double[] offset = Offset(n, 1);
        double[] z = [.. y.Zip(offset, (a, b) => a - b)];
        double[] w = BassWeightedCoreTests.Weights(n, seed: 21), w0 = BassWeightedCoreTests.Weights(n, seed: 22);
        BassPriorReference reference = frozenGram ? BassPriorReference.FrozenGram(y, w0) : BassPriorReference.Moments(y);
        BassConfig config = BassWeightedCoreTests.EnumerationConfig(Sigma2Prior.ScaledToResponse(0.5));
        BassModel model = BassModel.Create(config, data, new Snapshot(r => new WorkingResponse(r, offset, w, priorGram: frozenGram ? WeightedPriorGram.Frozen : WeightedPriorGram.Unweighted, reference: reference)));
        string label = frozenGram ? "the frozen Gram" : "a frozen level";
        BassWeightedCoreTests.CheckKernel(model, label, _out);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        BassWeightedCoreTests.DensePosterior(space, z, w, frozenGram ? w0 : Unit(n), reference.Level, 0.1, null, model.Target.Sigma2, w.Sum(Math.Log));
        BassWeightedCoreTests.CheckVisits(space, space.Posterior(), c => new BassChain(model, [], Math.Log(0.1), new Xoshiro256PlusPlus(400 + c), new Xoshiro256PlusPlus(500 + c)), _out, label);
    }

    // ---- The hierarchical family ----

    // ln ∫ Gamma(r) exp(logValue(IG(a, r))) dr, by the trapezoid on v = ln r about r's mean.
    private static double OverRate(ResolvedSigma2Rate rate, Func<InverseGamma, double> logValue)
    {
        const int points = 24_001;
        const double low = -80.0, high = 40.0;
        double step = (high - low) / (points - 1), top = double.NegativeInfinity;
        var values = new double[points];
        for (int i = 0; i < points; i++)
        {
            double v = rate.CentreLogRate + low + i * step, density = rate.LogDensity(v);
            values[i] = double.IsNegativeInfinity(density) ? double.NegativeInfinity : density + logValue(rate.PriorAt(v));
            top = Math.Max(top, values[i]);
        }
        double sum = 0.0;
        foreach (double value in values) sum += Math.Exp(value - top);
        return top + Math.Log(sum * step);
    }

    [Fact]
    public void TheRate_IsItsGamma()
    {
        ResolvedSigma2Rate rate = Sigma2Prior.Hierarchical(0.5, 2.0).RateFor(9.0)!;
        Assert.Equal(4.5, rate.Centre);
        Assert.Equal(2.0, rate.Hyper.Shape);
        Assert.Equal(2.0 / 4.5, rate.Hyper.Rate, 14);
        Assert.Equal(Math.Log(4.5), rate.CentreLogRate, 14);
        const int points = 20_001;
        double total = 0.0, step = 80.0 / (points - 1);
        for (int i = 0; i < points; i++) total += Math.Exp(rate.LogDensity(-40.0 + i * step)) * step;
        Close(1.0, total, 1e-10, "the density's integral on v");
        Close(rate.Hyper.LogPdf(3.0) + Math.Log(3.0), rate.LogDensity(Math.Log(3.0)), 1e-12, "the Gamma with its Jacobian");
        Assert.Equal(double.NegativeInfinity, rate.LogDensity(-800.0));
        Assert.Equal(double.NegativeInfinity, rate.LogDensity(800.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => rate.LogDensity(double.NaN));
        Close(InverseGamma.ShapeScale(0.5, 3.0).LogPdf(0.7), rate.LogConditional(0.7, Math.Log(3.0)), 1e-12, "σ²'s conditional prior");
        InverseGamma at = rate.PriorAt(Math.Log(3.0));
        Assert.Equal(0.5, at.Shape);
        Close(3.0, at.Scale, 1e-14, "the inverse gamma's rate");

        var rng = new Xoshiro256PlusPlus(3);
        double[] drawn = [.. Enumerable.Range(0, 200_000).Select(_ => Math.Exp(rate.SampleLogRate(rng)))];
        Close(4.5, drawn.Average(), 0.01, "r's mean");

        // A tiny shape's draws underflow, and are drawn again.
        ResolvedSigma2Rate tiny = Sigma2Prior.Hierarchical(0.5, 1e-3).RateFor(1.0)!;
        for (int i = 0; i < 2_000; i++) Assert.False(double.IsNegativeInfinity(tiny.LogDensity(tiny.SampleLogRate(rng))));
    }

    // At a = h = ½ the marginal over r is the marginal given σ² integrated against σ ~ half-Cauchy(s): the half-Cauchy read from its own
    // distribution, on the structure's own S, tempered too. Nothing of the mixing is read on that side.
    [Theory]
    [InlineData(1.0)]
    [InlineData(0.3)]
    public void AtOneHalf_TheFamilyIsTheHalfCauchy(double b)
    {
        BassData data = BassSweepTests.Data(large: true);
        BassModel model = BassModel.Create(BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.HalfCauchy(), 0.2), data);
        ResolvedSigma2Rate rate = model.Target.Sigma2Rate!;
        double s = Math.Sqrt(model.Resolved.ResponseMeanSquare);
        var halfCauchy = HalfCauchy.WithScale(s);
        foreach (int count in new[] { 0, 2, 5 })
        {
            BassState state = BassMixing.Loaded(model, BassWeightedCoreTests.Terms[..count], 0.0, out BassRefusal? refusal) ?? throw new InvalidOperationException($"{refusal}");
            GPriorMarginal marginal = model.Target.Marginal(state.Fit, 0.2, b);
            double mixed = OverRate(rate, prior => marginal.LogMarginal(prior));

            // On u = ln σ: exp(ln p(y | σ²)) · half-Cauchy(σ) · σ du.
            const int points = 40_001;
            const double low = -30.0, high = 30.0;
            double step = (high - low) / (points - 1), top = double.NegativeInfinity;
            var values = new double[points];
            for (int i = 0; i < points; i++)
            {
                double u = low + i * step, sigma = Math.Exp(u);
                values[i] = marginal.LogMarginalGiven(sigma * sigma) + halfCauchy.LogPdf(sigma) + u;
                top = Math.Max(top, values[i]);
            }
            double sum = 0.0;
            foreach (double value in values) sum += Math.Exp(value - top);
            Close(top + Math.Log(sum * step), mixed, 1e-8, $"M = {count}");
        }
    }

    // h → ∞ is the scaled prior, and h → 0 is the Jeffreys prior for any a, across structures: the family's two ends.
    [Theory]
    [InlineData(0.5)]
    [InlineData(2.0)]
    public void TheFamily_RunsFromThePlugInToTheInvariant(double a)
    {
        BassData data = BassSweepTests.Data(large: true);
        double Gap(BassConfig config, double h)
        {
            BassModel model = BassModel.Create(config with { Sigma2 = Sigma2Prior.Hierarchical(a, h) }, data);
            double Value(int count)
            {
                BassState state = BassMixing.Loaded(model, BassWeightedCoreTests.Terms[..count], 0.0, out BassRefusal? refusal) ?? throw new InvalidOperationException($"{refusal}");
                GPriorMarginal marginal = model.Target.Marginal(state.Fit, 0.2);
                return OverRate(model.Target.Sigma2Rate!, prior => marginal.LogMarginal(prior));
            }
            return Value(1) - Value(5);
        }
        double Plain(BassConfig config, Sigma2Prior sigma2)
        {
            BassModel model = BassModel.Create(config with { Sigma2 = sigma2 }, data);
            double Value(int count)
            {
                BassState state = BassMixing.Loaded(model, BassWeightedCoreTests.Terms[..count], 0.0, out BassRefusal? refusal) ?? throw new InvalidOperationException($"{refusal}");
                return model.Target.Marginal(state.Fit, 0.2).LogMarginal(model.Target.Sigma2);
            }
            return Value(1) - Value(5);
        }
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(a), 0.2);
        double scaled = Plain(config, Sigma2Prior.ScaledToResponse(a)), jeffreys = Plain(config, Sigma2Prior.Jeffreys);
        _out.WriteLine($"a = {a}: plug-in {scaled:F5}, h = 10⁴ {Gap(config, 1e4):F5}; Jeffreys {jeffreys:F5}, h = 10⁻⁵ {Gap(config, 1e-5):F5}");
        Close(scaled, Gap(config, 1e4), 1e-3, "h → ∞");
        Close(jeffreys, Gap(config, 1e-5), 1e-3, "h → 0");
        Assert.True(Math.Abs(scaled - jeffreys) > 1e-3, "The fixture should tell the two ends apart.");
    }

    // The posterior over the structures with r integrated, by quadrature on each structure's own marginal.
    private static double[] OverRatePosterior(BassKernelEnumeration.Space space, BassModel hierarchical, double tau)
    {
        var logs = new double[space.States.Count];
        for (int i = 0; i < logs.Length; i++)
        {
            BassKernelEnumeration.State state = space.States[i];
            GPriorMarginal marginal = hierarchical.Target.Marginal(state.Loaded.Fit, tau);
            logs[i] = hierarchical.Prior.LogCollectionWeight(state.Terms.Length) + state.Terms.Sum(hierarchical.Prior.LogPrior)
                      + OverRate(hierarchical.Target.Sigma2Rate!, prior => marginal.LogMarginal(prior));
        }
        double top = logs.Max();
        double[] p = [.. logs.Select(v => Math.Exp(v - top))];
        double total = p.Sum();
        return [.. p.Select(v => v / total)];
    }

    // σ² integrated: the chain's slice step on v, and the moves given it, visit the posterior with r integrated out.
    [Fact]
    public void Enumeration_TheRatesBlock_WalksThePosterior()
    {
        BassData data = BassSweepTests.Data();
        BassModel plain = BassModel.Create(BassWeightedCoreTests.EnumerationConfig(Sigma2Prior.ScaledToResponse(0.5)), data);
        BassModel model = BassModel.Create(BassWeightedCoreTests.EnumerationConfig(Sigma2Prior.HalfCauchy()), data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(plain);
        double[] pi = OverRatePosterior(space, model, 0.1);
        double apart = pi.Zip(space.Posterior(), (x, y) => Math.Abs(x - y)).Max();
        _out.WriteLine($"the largest gap from the plug-in's posterior {apart:E2}");
        BassWeightedCoreTests.CheckVisits(space, pi, c => new BassChain(model, [], Math.Log(0.1), new Xoshiro256PlusPlus(610 + c), new Xoshiro256PlusPlus(620 + c)), _out, "the hierarchical prior");
    }

    private sealed class Sigma2Block(BassTarget target) : IBlock<ReversibleJumpChain<BassState>>
    {
        public string Name => "sigma2";

        public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
        {
            double tau = target.HeldTau!.Value;
            double sigma2 = target.Marginal(state.Current.Fit, tau).Sigma2Posterior(target.Sigma2At(state.Current)).Sample(rng);
            state.ReplaceTarget(target.Given(sigma2));
            return true;
        }
    }

    // σ² drawn apart, the target given it, its rate reading it: the form E1's w reads. The same posterior over structures.
    [Fact]
    public void Enumeration_TheRatesBlock_ReadsASigma2DrawnApart()
    {
        BassData data = BassSweepTests.Data();
        BassModel plain = BassModel.Create(BassWeightedCoreTests.EnumerationConfig(Sigma2Prior.ScaledToResponse(0.5)), data);
        BassModel model = BassModel.Create(BassWeightedCoreTests.EnumerationConfig(Sigma2Prior.HalfCauchy()), data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(plain);
        double[] pi = OverRatePosterior(space, model, 0.1);
        const int chains = 4, burn = 2_000, sweeps = 30_000;
        var visits = new double[pi.Length];
        for (int c = 0; c < chains; c++)
        {
            var moves = new BassStructureMoves(model, [], Math.Log(0.1));
            var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target.Given(0.1), moves.Start, new Xoshiro256PlusPlus(800 + c));
            var scan = new BlockScan<ReversibleJumpChain<BassState>>([new MoveStep<BassState>(), new Sigma2Block(model.Target), new BassSigma2RateBlock()], chain, new Xoshiro256PlusPlus(900 + c), x => x);
            for (int s = 0; s < burn; s++) scan.Step();
            for (int s = 0; s < sweeps; s++)
            {
                scan.Step();
                visits[space.Find(chain.Current.Terms)] += 1.0;
            }
        }
        double total = chains * sweeps;
        foreach (int i in Enumerable.Range(0, pi.Length).OrderByDescending(i => pi[i]).Take(10))
            Assert.True(Math.Abs(pi[i] - visits[i] / total) <= 0.02, $"state {i}: π {pi[i]:F4}, visits {visits[i] / total:F4}");
        _out.WriteLine($"largest π {pi.Max():F4}, its visits {visits[Array.IndexOf(pi, pi.Max())] / total:F4}");
    }

    // At a fixed structure, the block's draws of v follow its conditional with σ² integrated, by quadrature, whether σ² is integrated
    // or drawn apart beside it: a rate that read no σ² when given it would follow its prior instead.
    [Theory]
    [InlineData(0.5, 0.5, false)]
    [InlineData(2.0, 5.0, false)]
    [InlineData(0.5, 0.5, true)]
    [InlineData(2.0, 5.0, true)]
    public void TheRatesBlock_FollowsItsConditional(double a, double h, bool drawnApart)
    {
        BassData data = BassSweepTests.Data(large: true);
        BassModel model = BassModel.Create(BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.Hierarchical(a, h), 0.2), data);
        ResolvedSigma2Rate rate = model.Target.Sigma2Rate!;
        var moves = new BassStructureMoves(model, BassWeightedCoreTests.Terms[..2], Math.Log(0.2));
        var chain = new ReversibleJumpChain<BassState>([new BassCoefficientPriorTests.Stay()], drawnApart ? model.Target.Given(1.0) : model.Target, moves.Start, new Xoshiro256PlusPlus(7));
        IBlock<ReversibleJumpChain<BassState>>[] blocks = drawnApart ? [new Sigma2Block(model.Target), new BassSigma2RateBlock()] : [new BassSigma2RateBlock()];
        var scan = new BlockScan<ReversibleJumpChain<BassState>>(blocks, chain, new Xoshiro256PlusPlus(8), c => c);
        const int sweeps = 60_000;
        for (int i = 0; i < 1_000; i++) scan.Step();
        double[] trace = [.. Enumerable.Range(0, sweeps).Select(_ => scan.Step().Current.LogRate)];
        GPriorMarginal marginal = model.Target.Marginal(moves.Start.Fit, 0.2);
        const int points = 24_001;
        double step = 120.0 / (points - 1), top = double.NegativeInfinity;
        var logs = new double[points];
        for (int i = 0; i < points; i++)
        {
            double v = rate.CentreLogRate - 80.0 + i * step, density = rate.LogDensity(v);
            logs[i] = double.IsNegativeInfinity(density) ? density : density + marginal.LogMarginal(rate.PriorAt(v));
            top = Math.Max(top, logs[i]);
        }
        double mass = 0.0, first = 0.0;
        for (int i = 0; i < points; i++)
        {
            double weight = Math.Exp(logs[i] - top);
            mass += weight;
            first += weight * (rate.CentreLogRate - 80.0 + i * step);
        }
        double expected = first / mass, error = BassCoefficientPriorTests.BatchError(trace);
        double priorMean = SpecialFunctions.Digamma(h) - Math.Log(rate.Hyper.Rate);
        _out.WriteLine($"E[v] {expected:F4}, the chain's {trace.Average():F4} ± {error:F4}; the prior's {priorMean:F4}");
        Assert.True(Math.Abs(expected - priorMean) > 0.1, "The fixture should tell the conditional from the prior.");
        Assert.True(Math.Abs(expected - trace.Average()) <= 4.0 * error + 1e-3, $"E[v] {expected}, the chain's {trace.Average()} ± {error}");
    }

    // τ's block and a chain's composition read the target the chain holds now, not the one it was built with.
    [Fact]
    public void TheTauBlockAndTheComposition_ReadTheChainsTarget()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2) with { Coefficients = BassCoefficientPrior.G(TauPrior.UnitInformationGamma()) };
        BassModel model = BassModel.Create(config, data);
        Assert.True(model.Target.SamplesTau);
        var chain = new BassChain(model, BassWeightedCoreTests.Terms[..2], Math.Log(0.02), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2));
        BassTarget given = model.Target.Given(0.37);
        chain.Moves.ReplaceTarget(given);
        Assert.Same(given, chain.Moves.Target);
        new BassTauBlock().Update(chain.Moves, new Xoshiro256PlusPlus(3));
        Close(given.Evaluate(chain.Current).LogMarginal, ((BassEvaluation)chain.Moves.CurrentEvaluation).LogMarginal, 1e-13, "τ's block's evaluation");
        Assert.Equal(0.37, chain.Compose(new Xoshiro256PlusPlus(4), new BassTermTable()).Sigma2);
    }

    // The state's rate reaches the target's values, the composition and a retained draw; the evaluation's own prior says so.
    [Fact]
    public void TheStatesRate_ReachesWhatReadsIt()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.Hierarchical(0.5, 2.0), 0.2);
        BassModel model = BassModel.Create(config, data);
        ResolvedSigma2Rate rate = model.Target.Sigma2Rate!;
        Assert.Equal(model.Resolved.Sigma2Rate, rate);
        Assert.Equal(rate.Centre, model.Target.Sigma2!.Value.Scale, 12);
        var moves = new BassStructureMoves(model, BassWeightedCoreTests.Terms[..2], Math.Log(0.2));
        Assert.Equal(rate.CentreLogRate, moves.Start.LogRate);
        BassState state = moves.Start.WithLogRate(rate.CentreLogRate + 1.3);
        Assert.Equal(moves.Start.LogTau, state.LogTau);
        Assert.Equal(state.LogRate, state.WithLogTau(0.1).LogRate);
        BassEvaluation evaluation = model.Target.Evaluate(state);
        InverseGamma atRate = rate.PriorAt(state.LogRate);
        Assert.Equal(model.Target.Sigma2At(state)!.Value.Scale, atRate.Scale);
        Assert.Equal(atRate.Scale, evaluation.Sigma2!.Value.Scale);
        Close(model.Target.Marginal(state.Fit, 0.2).LogMarginal(atRate), evaluation.LogMarginal, 1e-13, "the marginal at the state's rate");
        Close(model.Target.LogStructurePrior(state.Fit) + rate.LogDensity(state.LogRate), evaluation.LogPrior, 1e-13, "the prior with the rate's density");
        Assert.NotEqual(evaluation.LogMarginal, model.Target.Evaluate(moves.Start).LogMarginal);

        // Given σ², the prior carries σ²'s conditional, and the marginal does not read the rate.
        BassTarget given = model.Target.Given(0.4);
        BassEvaluation conditioned = given.Evaluate(state);
        Close(given.LogStructurePrior(state.Fit) + rate.LogDensity(state.LogRate) + rate.LogConditional(0.4, state.LogRate), conditioned.LogPrior, 1e-13, "the prior given σ²");
        Assert.Equal(given.Evaluate(moves.Start).LogMarginal, conditioned.LogMarginal);

        // A rate off the doubles lies outside the support, with no value to compute.
        Assert.True(model.Target.Evaluate(state.WithLogRate(800.0)).OffSupport);
        Assert.True(model.Target.Evaluate(state.WithLogRate(-800.0)).OffSupport);
        Assert.Equal(double.NegativeInfinity, model.Target.Evaluate(state.WithLogRate(800.0)).LogMarginal);

        // The fit's overloads have no rate to read.
        Assert.Throws<InvalidOperationException>(() => model.Target.Evaluate(state.Fit));

        // A retained draw records its rate, and a draw's σ² reads it: at a huge rate σ² is huge.
        BassChain chain = new(model, BassWeightedCoreTests.Terms[..2], Math.Log(0.2), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2), logRate: rate.CentreLogRate + 200.0);
        BassDraw draw = chain.Compose(new Xoshiro256PlusPlus(3), new BassTermTable());
        Assert.Equal(Math.Exp(rate.CentreLogRate + 200.0), draw.Sigma2Rate!.Value);
        Assert.True(draw.Sigma2 > 1e70, $"σ² = {draw.Sigma2}");

        // Under every other prior the state has no rate, and a draw records none.
        BassModel scaled = BassModel.Create(config with { Sigma2 = Sigma2Prior.ScaledToResponse(0.5) }, data);
        var other = new BassChain(scaled, BassWeightedCoreTests.Terms[..2], Math.Log(0.2), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2), logRate: 4.0);
        Assert.True(double.IsNaN(other.Current.LogRate));
        Assert.Null(other.Compose(new Xoshiro256PlusPlus(3), new BassTermTable()).Sigma2Rate);
    }

    // The prior sampler draws r from its prior, on the stream after τ's; a prior without a rate spends nothing on it.
    [Fact]
    public void ThePriorSampler_DrawsTheRate()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassModel model = BassModel.Create(BassSweepTests.Config(maxBasis: 3) with { Sigma2 = Sigma2Prior.Hierarchical(0.5, 2.0) }, data);
        var sampler = new BassPriorSampler(model);
        var rng = new Xoshiro256PlusPlus(2);
        double[] rates = [.. Enumerable.Range(0, 20_000).Select(_ => Math.Exp(sampler.DrawState(rng).State.LogRate))];
        Close(model.Target.Sigma2Rate!.Centre, rates.Average(), 0.03, "r's prior mean");
        BassModel scaled = BassModel.Create(BassSweepTests.Config(maxBasis: 3) with { Sigma2 = Sigma2Prior.ScaledToResponse(0.5) }, data);
        Assert.True(double.IsNaN(new BassPriorSampler(scaled).DrawState(new Xoshiro256PlusPlus(2)).State.LogRate));
    }

    // The block sits last in the scale slot, on every sampler, behind what the observation model draws.
    [Fact]
    public void TheRatesBlock_IsLastInTheScaleSlot()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassModel model = BassModel.Create(BassSweepTests.Config(maxBasis: 3) with { Sigma2 = Sigma2Prior.HalfCauchy() }, data);
        Assert.Equal(["move-step", "tau", "sigma2 rate"], BassSweepSlots.Blocks(model).Select(b => b is MoveStep<BassState> ? "move-step" : b.Name));
        Assert.Equal(["tau", "sigma2 rate"], BassSweepSlots.Blocks(model, BassSweepSlot.Tau).Select(b => b.Name));
        BassModel gaussian = BassModel.Create(BassSweepTests.Config(maxBasis: 3), data);
        Assert.DoesNotContain(BassSweepSlots.Blocks(gaussian), b => b is BassSigma2RateBlock);
    }

    // A run's evidence against the enumeration's, r integrated by quadrature: the rate through a ladder down to its exact reference
    // rung, and through SMC from its prior draws.
    [Theory]
    [InlineData("ladder")]
    [InlineData("linear ladder")]
    [InlineData("smc")]
    public void AHierarchicalPrior_ReportsTheEvidence_ThatMatchesTheEnumeration(string sampler)
    {
        BassData data = BassSweepTests.Data();
        BassSampling sampling = sampler.EndsWith("ladder")
            ? new BassSampling
            {
                Chains = 2, BurnIn = 1_500, Samples = 3_000, Moves = new EqualOverAvailable(), ReportEvidence = true,
                Tempering = new TemperingSpec { Levels = sampler == "ladder" ? 8 : 16, HotEnd = 0.05, ReferenceRung = true, Form = sampler == "ladder" ? TemperedForm.Integrated : TemperedForm.Linear },
            }
            : new BassSampling { Moves = new EqualOverAvailable(), ReportEvidence = true, Smc = new SmcSpec { Particles = 4_000, Replicates = 4 } };
        BassConfig config = BassWeightedCoreTests.EnumerationConfig(Sigma2Prior.HalfCauchy()) with { Sampling = sampling };
        BassFit fit = Bass.Fit(config, data);
        BassEvidence evidence = fit.Evidence ?? throw new InvalidOperationException("The run reports no evidence.");
        Assert.Equal(Sigma2Prior.HalfCauchy(), evidence.Priors.Sigma2Prior);
        Assert.NotNull(evidence.Priors.Sigma2Rate);
        BassModel model = BassModel.Create(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(BassModel.Create(BassWeightedCoreTests.EnumerationConfig(Sigma2Prior.ScaledToResponse(0.5)), data));
        var logPosterior = new double[space.States.Count];
        var logPrior = new double[space.States.Count];
        for (int i = 0; i < logPrior.Length; i++)
        {
            BassKernelEnumeration.State state = space.States[i];
            logPrior[i] = model.Prior.LogCollectionWeight(state.Terms.Length) + state.Terms.Sum(model.Prior.LogPrior);
            GPriorMarginal marginal = model.Target.Marginal(state.Loaded.Fit, 0.1);
            logPosterior[i] = logPrior[i] + OverRate(model.Target.Sigma2Rate!, prior => marginal.LogMarginal(prior));
        }
        double top = logPosterior.Max(), topPrior = logPrior.Max();
        double exact = top + Math.Log(logPosterior.Sum(v => Math.Exp(v - top))) - topPrior - Math.Log(logPrior.Sum(v => Math.Exp(v - topPrior)));
        _out.WriteLine($"{sampler}: ln Z {evidence.LogEvidence:F3}, the enumeration's {exact:F3}");
        Assert.True(Math.Abs(evidence.LogEvidence - exact) <= 0.15, $"ln Z {evidence.LogEvidence:F3}, the enumeration's {exact:F3}");
        Assert.All(fit.Draws.SelectMany(c => c), draw => Assert.NotNull(draw.Sigma2Rate));
    }

    // A flat-histogram run carries the rate through its pilot's linear ladder and its multicanonical chains.
    [Fact]
    public void AFlatHistogramRun_CarriesTheRate()
    {
        BassConfig config = BassFlatHistogramTests.WithFlatHistogram(BassWeightedCoreTests.EnumerationConfig(Sigma2Prior.HalfCauchy()), BassFlatHistogramTests.Spec());
        BassFit fit = Bass.Fit(config, BassSweepTests.Data());
        Assert.NotEmpty(fit.Draws.SelectMany(c => c));
        Assert.All(fit.Draws.SelectMany(c => c), draw => Assert.True(draw.Sigma2Rate > 0.0, $"r = {draw.Sigma2Rate}"));
        Assert.True(fit.Draws.SelectMany(c => c).Select(d => d.Sigma2Rate).Distinct().Count() > 10, "The rate should move.");
    }

    // A draw's collapsed predictive reads its own rate: the same as a declared IG(a, r) at that r. Each r is computed once per τ.
    [Fact]
    public void ThePointwiseDensities_ReadTheDrawsRate()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig config = BassWeightedCoreTests.Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.Hierarchical(0.5, 2.0), 0.2);
        BassModel model = BassModel.Create(config, data);
        double[] rates = [0.4, 1.7];
        var table = new BassTermTable();
        BassDraw composed = new BassChain(model, BassWeightedCoreTests.Terms[..3], Math.Log(0.2), new Xoshiro256PlusPlus(1), new Xoshiro256PlusPlus(2)).Compose(new Xoshiro256PlusPlus(3), table);
        IReadOnlyList<BassDraw>[] draws = [[composed with { Sigma2Rate = rates[0] }, composed with { Sigma2Rate = rates[1] }, composed with { Sigma2Rate = rates[0] }]];
        PointwiseLogDensity reading = BassPointwise.LogDensity(model, table, draws, content: PointwiseContent.CollapsedPredictive);
        for (int r = 0; r < rates.Length; r++)
        {
            BassConfig declared = config with { Sigma2 = Sigma2Prior.InverseGamma(0.5, rates[r]) };
            BassModel plain = BassModel.Create(declared, data);
            IReadOnlyList<BassDraw>[] one = [[composed with { Sigma2Rate = null }]];
            PointwiseLogDensity expected = BassPointwise.LogDensity(plain, table, one, content: PointwiseContent.CollapsedPredictive);
            for (int u = 0; u < expected.Values.GetLength(1); u++) Close(expected.Values[0, u], reading.Values[r, u], 1e-12, $"rate {rates[r]}, unit {u}");
        }
        for (int u = 0; u < reading.Values.GetLength(1); u++) Assert.Equal(reading.Values[0, u], reading.Values[2, u]);
    }

    [Fact]
    public void AHierarchicalPrior_IsAPriorKind_ThatChangesNoOtherKindsSerializedForm()
    {
        BassConfig config = new();
        Assert.DoesNotContain("\"HyperShape\"", config.ToJson());
        Assert.Contains("\"HyperShape\":0.5", (config with { Sigma2 = Sigma2Prior.HalfCauchy() }).ToJson());
        BassConfig back = BassConfig.FromJson((config with { Sigma2 = Sigma2Prior.Hierarchical(0.5, 2.0) }).ToJson());
        Assert.Equal(Sigma2Prior.Hierarchical(0.5, 2.0), back.Sigma2);
        Assert.Equal("Hierarchical(0.5, 2)", back.Sigma2.ToString());
        Assert.True(back.Sigma2.IsProper);
        Assert.True(back.Sigma2.ReadsResponse);
        Assert.True(back.Sigma2.IsHierarchical);
        Assert.False(back.Sigma2.IsHeld);
        Assert.Equal(0.5 * 9.0, back.Sigma2.For(9.0)!.Value.Scale);
        Assert.Null(back.Sigma2.For(0.0));
        Assert.Null(back.Sigma2.RateFor(0.0));
        Assert.Null(Sigma2Prior.ScaledToResponse().RateFor(1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sigma2Prior.Hierarchical(0.0, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sigma2Prior.Hierarchical(1.0, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sigma2Prior.Hierarchical(1.0, double.PositiveInfinity));
        Assert.Throws<ArgumentException>(() => new Sigma2Prior(Sigma2PriorKind.ScaledToResponse, 1.0, hyperShape: 1.0));
        Assert.Throws<ArgumentException>(() => new Sigma2Prior(Sigma2PriorKind.Hierarchical, 1.0, 1.0, hyperShape: 1.0));
        Assert.Throws<ArgumentException>(() => new WeightedNormalModel(sigma2Prior: Sigma2Prior.HalfCauchy()));
    }
}
