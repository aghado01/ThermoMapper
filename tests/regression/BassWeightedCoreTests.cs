using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline;
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
/// The weighted core and the marginal given σ²: the design under <see cref="WorkingResponse.Weights"/>, its two readings of a
/// weight and its two coefficient priors, and the structure step with σ² held.
/// <list type="bullet">
/// <item><b>Integer weights are duplicated rows</b> under <see cref="WeightMeaning.Replication"/> with the prior on the weighted
/// columns: through the design, the marginal and a chain's states and draws. Under <see cref="WeightMeaning.Precision"/> the
/// statistics and the marginal given σ² agree up to a constant free of the structure, and σ²'s integral does not.</item>
/// <item><b>Dense algebra.</b> S, the determinant term, the conditional mean and variances, and the leave-unit-out predictive,
/// against matrices built here from the columns, for each form, convention and prior.</item>
/// <item><b>σ²'s integral.</b> The marginal given σ² integrated over σ²'s prior by quadrature is the collapsed marginal.</item>
/// <item><b>Enumeration.</b> Detailed balance and stationarity on every structure of F, for a weighted target and for one given
/// σ², with the posterior checked against dense algebra and against a chain's visits.</item>
/// </list>
/// </summary>
public sealed class BassWeightedCoreTests
{
    private readonly ITestOutputHelper _out;
    public BassWeightedCoreTests(ITestOutputHelper output) => _out = output;

    private sealed class Observation(Func<BassResponse, WorkingResponse> working) : IBassObservationModel
    {
        public WorkingResponse Working(BassResponse response) => working(response);

        public IReadOnlyList<IBlock<ReversibleJumpChain<BassState>>> Blocks(BassSweepSlot slot) => [];
    }

    private static Observation Weighted(double[] weights, WeightMeaning meaning, WeightedPriorGram gram)
        => new(r => new WorkingResponse(r, weights: weights, meaning: meaning, priorGram: gram));

    private static double[] Weights(int n, int seed = 3)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        return [.. Enumerable.Range(0, n).Select(_ => Math.Exp((2.0 * rng.NextDouble() - 1.0) * Math.Log(4.0)))];
    }

    private static int[] Counts(int n) => [.. Enumerable.Range(0, n).Select(i => 1 + (i * 7 + i / 3) % 3)];

    private static double[] AsDoubles(int[] counts) => [.. counts.Select(c => (double)c)];

    // The rows repeated, each its count's times, in order.
    private static BassData Duplicated(BassData data, int[] counts)
    {
        IEnumerable<T> Repeat<T>(IEnumerable<T> values) => values.SelectMany((v, i) => Enumerable.Repeat(v, counts[i]));
        BassPredictor[] predictors = [.. data.Predictors.Select(p => (BassPredictor)new ContinuousPredictor(p.Name, [.. Repeat(((ContinuousPredictor)p).Values.ToArray())]))];
        return new BassData(predictors, new ScalarResponse([.. Repeat(data.Response.Values.ToArray())]));
    }

    private static BassTerm Hinge(int variable, int sign, int knot) => new(continuous: [new HingeFactor(variable, sign, knot)]);

    private static BassTerm Product(int knot0, int knot1) => new(continuous: [new HingeFactor(0, 1, knot0), new HingeFactor(1, -1, knot1)]);

    private static readonly BassTerm[] Terms = [Hinge(0, 1, 0), Hinge(1, -1, 2), Hinge(0, 1, 1), Product(0, 2), Hinge(1, 1, 1)];

    private static BassConfig Config(CoefficientForm form, GPriorCentre intercept, Sigma2Prior sigma2, double tau, int maxBasis = 5)
        => BassSweepTests.Config(maxBasis: maxBasis) with
        {
            Intercept = intercept,
            Sigma2 = sigma2,
            Coefficients = new BassCoefficientPrior(form, TauPrior.Held(tau)),
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = maxBasis, Support = new SupportRule(SupportKind.Count, 1.0) },
        };

    private static void Close(double expected, double actual, double relative, string what)
    {
        double tolerance = relative * Math.Max(1.0, Math.Abs(expected));
        Assert.True(Math.Abs(expected - actual) <= tolerance, $"{what}: {expected:R} against {actual:R}, tolerance {tolerance:E2}");
    }

    // ---- Unit weights ----

    [Theory]
    [InlineData(WeightMeaning.Precision, WeightedPriorGram.Weighted)]
    [InlineData(WeightMeaning.Precision, WeightedPriorGram.Unweighted)]
    [InlineData(WeightMeaning.Replication, WeightedPriorGram.Weighted)]
    [InlineData(WeightMeaning.Replication, WeightedPriorGram.Unweighted)]
    public void DeclaredUnitWeights_ReduceToGaussian_UnderEveryReading(WeightMeaning meaning, WeightedPriorGram gram)
    {
        BassData data = BassSweepTests.Data();
        double[] unit = [.. Enumerable.Repeat(1.0, data.Response.Observations)];
        foreach (CoefficientForm form in new CoefficientForm[] { new GForm(), new RidgeForm() })
            foreach (GPriorCentre intercept in new[] { GPriorCentre.Mean, GPriorCentre.Zero, GPriorCentre.Flat })
            {
                BassConfig config = Config(form, intercept, Sigma2Prior.ScaledToResponse(0.5), 0.2);
                BassModel gaussian = BassModel.Create(config, data), weighted = BassModel.Create(config, data, Weighted(unit, meaning, gram));
                Close(gaussian.Resolved.CentreLevel, weighted.Resolved.CentreLevel, 1e-13, "the level");
                Close(gaussian.Resolved.ResponseMeanSquare, weighted.Resolved.ResponseMeanSquare, 1e-13, "s²");
                BassState a = BassMixing.Loaded(gaussian, Terms[..3], Math.Log(0.2), out _)!, b = BassMixing.Loaded(weighted, Terms[..3], Math.Log(0.2), out _)!;
                foreach (double t in new[] { 1.0, 0.6 })
                {
                    GPriorMarginal m = gaussian.Target.Marginal(a.Fit, 0.2, t), w = weighted.Target.Marginal(b.Fit, 0.2, t);
                    Close(m.S, w.S, 1e-10, $"S, {form}, {intercept}");
                    Close(m.LogDeterminantTerm, w.LogDeterminantTerm, 1e-10, $"the determinant term, {form}, {intercept}");
                    Close(m.LogMarginal(gaussian.Target.Sigma2), w.LogMarginal(weighted.Target.Sigma2), 1e-10, $"the marginal, {form}, {intercept}");
                    Close(m.LogMarginalSlope(gaussian.Target.Sigma2), w.LogMarginalSlope(weighted.Target.Sigma2), 1e-9, $"its slope, {form}, {intercept}");
                }
            }
    }

    // ---- Integer weights and duplicated rows ----

    // Two designs walked through births, a refused repeat, a change, a rebuild and a death, compared at every step.
    [Theory]
    [InlineData(WeightMeaning.Replication, WeightedPriorGram.Weighted)]
    [InlineData(WeightMeaning.Replication, WeightedPriorGram.Unweighted)]
    [InlineData(WeightMeaning.Precision, WeightedPriorGram.Weighted)]
    [InlineData(WeightMeaning.Precision, WeightedPriorGram.Unweighted)]
    public void IntegerWeights_EqualDuplicatedRows_ThroughTheDesign(WeightMeaning meaning, WeightedPriorGram gram)
    {
        BassData data = BassSweepTests.Data(large: true);
        int[] counts = Counts(data.Response.Observations);
        BassData repeated = Duplicated(data, counts);
        BassConfig config = Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.ScaledToResponse(0.5), 0.2, maxBasis: 6);
        BassDesign a = BassModel.Create(config, data, Weighted(AsDoubles(counts), meaning, gram)).NewDesign();
        BassDesign b = BassModel.Create(config, repeated).NewDesign();
        void Same(string step)
        {
            Assert.Equal(b.Count, a.Count);
            Close(b.ResponseMean, a.ResponseMean, 1e-13, $"{step}: ȳ");
            Close(b.ResponseSquares, a.ResponseSquares, 1e-12, $"{step}: ‖ỹ‖²");
            Close(b.ResidualSquares, a.ResidualSquares, 1e-11, $"{step}: ‖r‖²");
            Close(b.FittedSquares, a.FittedSquares, 1e-11, $"{step}: Q");
            if (a.Count > 0) Close(b.MinAngle, a.MinAngle, 1e-10, $"{step}: the smallest angle");
            for (int j = 0; j < a.Count; j++) Close(b.Coefficients[j], a.Coefficients[j], 1e-10, $"{step}: β̂_{j}");
        }
        void Proposed(DesignProposal p, DesignProposal q, string step)
        {
            Assert.Equal(q.Refusal, p.Refusal);
            if (q.Refusal is not null) return;
            Assert.Equal(q.Columns, p.Columns);
            Close(q.ResidualSquares, p.ResidualSquares, 1e-11, $"{step}: the proposal's ‖r‖²");
            Close(q.FittedSquares, p.FittedSquares, 1e-11, $"{step}: the proposal's Q");
            Close(q.MinAngle, p.MinAngle, 1e-10, $"{step}: the proposal's angle");
        }
        Same("intercept");
        foreach (BassTerm term in Terms)
        {
            Proposed(a.ProposeBirth(term), b.ProposeBirth(term), $"birth of {term}");
            a.Accept();
            b.Accept();
            Same($"birth of {term}");
        }
        Assert.Equal(BassRefusal.Numerical, a.ProposeBirth(Terms[0]).Refusal);
        Proposed(a.ProposeChange(1, Hinge(0, -1, 3)), b.ProposeChange(1, Hinge(0, -1, 3)), "change");
        a.Accept();
        b.Accept();
        Same("change");
        a.Rebuild();
        b.Rebuild();
        Same("rebuild");
        Proposed(a.ProposeDeath(0), b.ProposeDeath(0), "death");
        a.Accept();
        b.Accept();
        Same("death");
        BassDesign c = BassModel.Create(config, data, Weighted(AsDoubles(counts), meaning, gram)).NewDesign();
        Assert.Null(c.Load(a.Terms));
        Close(a.ResidualSquares, c.ResidualSquares, 1e-11, "a load's ‖r‖²");
        Close(a.FittedSquares, c.FittedSquares, 1e-11, "a load's Q");
    }

    [Theory]
    [InlineData("g", GPriorCentre.Mean)]
    [InlineData("g", GPriorCentre.Zero)]
    [InlineData("g", GPriorCentre.Flat)]
    [InlineData("ridge", GPriorCentre.Mean)]
    [InlineData("ridge", GPriorCentre.Zero)]
    [InlineData("ridge", GPriorCentre.Flat)]
    public void IntegerWeights_EqualDuplicatedRows_ThroughTheMarginal(string form, GPriorCentre intercept)
    {
        BassData data = BassSweepTests.Data(large: true);
        int[] counts = Counts(data.Response.Observations);
        BassData repeated = Duplicated(data, counts);
        CoefficientForm coefficient = form == "g" ? new GForm() : new RidgeForm();
        foreach (Sigma2Prior sigma2 in new[] { Sigma2Prior.ScaledToResponse(0.5), Sigma2Prior.Jeffreys, Sigma2Prior.InverseGamma(2.0, 0.3), Sigma2Prior.Held(0.4) })
        {
            BassConfig config = Config(coefficient, intercept, sigma2, 0.2);
            BassModel a = BassModel.Create(config, data, Weighted(AsDoubles(counts), WeightMeaning.Replication, WeightedPriorGram.Weighted));
            BassModel b = BassModel.Create(config, repeated);
            Assert.Equal(b.Resolved.CentreLevel, a.Resolved.CentreLevel, 1e-12);
            Assert.Equal(b.Resolved.ResponseMeanSquare, a.Resolved.ResponseMeanSquare, 1e-12);
            foreach (int terms in new[] { 0, 2, 5 })
            {
                BassState x = BassMixing.Loaded(a, Terms[..terms], Math.Log(0.2), out _)!, y = BassMixing.Loaded(b, Terms[..terms], Math.Log(0.2), out _)!;
                foreach (double t in new[] { 1.0, 0.5, 0.05 })
                {
                    if (!a.Target.Defines(t)) continue;
                    ITargetEvaluation p = a.Target.Evaluate(x), q = b.Target.Evaluate(y);
                    Close(q.LogLikelihood(t), p.LogLikelihood(t), 1e-10, $"ℓ_{t}, {sigma2}, M = {terms}");
                    Close(((BassEvaluation)q).LogLikelihoodSlope(t), ((BassEvaluation)p).LogLikelihoodSlope(t), 1e-9, $"the slope at {t}, {sigma2}, M = {terms}");
                    Assert.Equal(b.Target.Marginal(y.Fit, 0.2, t).DegreesOfFreedom, a.Target.Marginal(x.Fit, 0.2, t).DegreesOfFreedom, 1e-10);
                }
            }
        }
    }

    // Inputs scattered, so no two terms share a column's span by the grid: a tied move's log ratio is rounding, on whose sign a
    // chain draws its uniform or not, and its stream then leaves the other's.
    private static BassData Scattered(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = rng.NextDouble();
            x1[i] = rng.NextDouble();
            y[i] = 1 + 2 * x0[i] + 3 * Math.Max(x0[i] - 0.5, 0) * Math.Max(0.5 - x1[i], 0) + 0.3 * BassQualification.Gaussian(rng);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
    }

    // A chain with these weights walks the duplicated rows' chain, state by state, on one seed: τ sampled on the unit information
    // that Σw counts, and each retained draw's σ² and coefficients.
    [Theory]
    [InlineData("g")]
    [InlineData("ridge")]
    public void IntegerWeights_EqualDuplicatedRows_ThroughAChain(string form)
    {
        BassData data = Scattered(60, 11);
        int[] counts = Counts(data.Response.Observations);
        BassData repeated = Duplicated(data, counts);
        BassConfig config = BassSweepTests.Config(maxBasis: 6);
        if (form == "ridge") config = config with { Coefficients = BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()) };
        BassModel weighted = BassModel.Create(config, data, Weighted(AsDoubles(counts), WeightMeaning.Replication, WeightedPriorGram.Weighted));
        BassModel plain = BassModel.Create(config, repeated);
        Assert.Equal(plain.Resolved.SupportMinimum, weighted.Resolved.SupportMinimum, 1e-12);
        Assert.Equal(plain.Resolved.Scale.UnitObservations, weighted.Resolved.Scale.UnitObservations, 1e-12);
        var a = new BassChain(weighted, [], Math.Log(0.02), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        var b = new BassChain(plain, [], Math.Log(0.02), new Xoshiro256PlusPlus(5), new Xoshiro256PlusPlus(6));
        var composeA = new Xoshiro256PlusPlus(7);
        var composeB = new Xoshiro256PlusPlus(7);
        var tableA = new BassTermTable();
        var tableB = new BassTermTable();
        int moved = 0;
        BassTerm[] last = [];
        for (int s = 0; s < 400; s++)
        {
            BassState x = a.Step(), y = b.Step();
            Assert.Equal(y.Terms, x.Terms);
            Close(y.LogTau, x.LogTau, 1e-8, $"sweep {s}: ln τ");
            if (!x.Terms.SequenceEqual(last)) moved++;
            last = [.. x.Terms];
            if (s % 40 != 39) continue;
            BassDraw p = a.Compose(composeA, tableA), q = b.Compose(composeB, tableB);
            Close(q.Sigma2, p.Sigma2, 1e-8, $"sweep {s}: σ²");
            Assert.Equal(q.Coefficients.Length, p.Coefficients.Length);
            for (int j = 0; j < p.Coefficients.Length; j++) Close(q.Coefficients[j], p.Coefficients[j], 1e-7, $"sweep {s}: a_{j}");
        }
        Assert.True(moved > 20, $"The chain moved {moved} times.");
    }

    // Under precision weights the likelihood is one observation at each row, so given σ² the marginal differs from the duplicates'
    // by a constant free of the structure and τ; integrated over σ² it does not, since ν counts N.
    [Fact]
    public void PrecisionWeights_DifferFromDuplicates_ByAConstantGivenSigma2()
    {
        BassData data = BassSweepTests.Data(large: true);
        int[] counts = Counts(data.Response.Observations);
        BassData repeated = Duplicated(data, counts);
        BassConfig config = Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.InverseGamma(2.0, 0.3), 0.2);
        BassModel a = BassModel.Create(config, data, Weighted(AsDoubles(counts), WeightMeaning.Precision, WeightedPriorGram.Weighted));
        BassModel b = BassModel.Create(config, repeated);
        const double sigma2 = 0.4;
        int n = counts.Length, total = counts.Sum();
        double sumLog = counts.Sum(c => Math.Log(c));
        var integrated = new List<double>();
        foreach (int terms in new[] { 0, 2, 5 })
        {
            BassState x = BassMixing.Loaded(a, Terms[..terms], Math.Log(0.2), out _)!, y = BassMixing.Loaded(b, Terms[..terms], Math.Log(0.2), out _)!;
            foreach (double t in new[] { 1.0, 0.5 })
            {
                double gap = b.Target.Marginal(y.Fit, 0.2, t).LogMarginalGiven(sigma2) - a.Target.Marginal(x.Fit, 0.2, t).LogMarginalGiven(sigma2);
                Close(t * (-0.5 * (total - n) * Math.Log(2.0 * Math.PI * sigma2) - 0.5 * sumLog), gap, 1e-10, $"the gap at M = {terms}, b = {t}");
            }
            integrated.Add(b.Target.Marginal(y.Fit, 0.2).LogMarginal(b.Target.Sigma2) - a.Target.Marginal(x.Fit, 0.2).LogMarginal(a.Target.Sigma2));
        }
        Assert.True(Math.Abs(integrated[2] - integrated[0]) > 1.0, $"σ²'s integral separates them: {string.Join(", ", integrated)}");
    }

    // ---- Dense algebra ----

    private sealed record Dense(double S, double LogDeterminant, double[] FittedMean, double[,] Covariance, double[,] X);

    private static double[] Solve(double[,] a, double[] b)
    {
        var factor = new DenseCholesky(a.GetLength(0));
        factor.Decompose(a);
        Assert.True(factor.Status.PositiveDefinite);
        return factor.Solve(b);
    }

    private static double LogDet(double[,] a)
    {
        var factor = new DenseCholesky(a.GetLength(0));
        factor.Decompose(a);
        Assert.True(factor.Status.PositiveDefinite);
        return factor.LogDet;
    }

    // The block from the columns alone: the prior on (α, β) in the prior's own centring, its likelihood at weights w and power b.
    // α flat drops its prior row and spends the observation in ν, which the production form states apart.
    private static Dense Reference(double[] y, double[] w, double[][] z, bool priorOnWeighted, GPriorCentre intercept, bool ridge, double tau, double b)
    {
        int n = y.Length, m = z.Length, p = m + 1;
        double mass = w.Sum();
        double[] priorWeight = priorOnWeighted ? w : [.. Enumerable.Repeat(1.0, n)];
        double priorMass = priorWeight.Sum();
        double level = intercept == GPriorCentre.Zero ? 0.0 : y.Zip(priorWeight, (v, k) => v * k).Sum() / priorMass;
        var x = new double[n, p];
        for (int i = 0; i < n; i++) x[i, 0] = 1.0;
        for (int j = 0; j < m; j++)
        {
            double centre = z[j].Zip(priorWeight, (v, k) => v * k).Sum() / priorMass;
            for (int i = 0; i < n; i++) x[i, j + 1] = z[j][i] - centre;
        }
        var d = new double[p, p];
        var a = new double[p, p];
        for (int j = 0; j < p; j++)
            for (int k = 0; k < p; k++)
            {
                double like = 0.0, own = 0.0;
                for (int i = 0; i < n; i++)
                {
                    like += w[i] * x[i, j] * x[i, k];
                    own += priorWeight[i] * x[i, j] * x[i, k];
                }
                a[j, k] = like;
                d[j, k] = own;
            }
        if (ridge)
        {
            for (int j = 1; j < p; j++)
                for (int k = 1; k < p; k++) d[j, k] = j == k ? d[j, j] / priorMass : 0.0;
            d[0, 0] = 1.0;
        }
        for (int k = 1; k < p; k++) d[0, k] = d[k, 0] = 0.0;
        bool flat = intercept == GPriorCentre.Flat;
        var precision = new double[p, p];
        for (int j = 0; j < p; j++)
            for (int k = 0; k < p; k++) precision[j, k] = b * a[j, k] + (flat && j == 0 ? 0.0 : tau * d[j, k]);
        var e = new double[n];
        for (int i = 0; i < n; i++) e[i] = y[i] - (flat ? 0.0 : level);
        var cross = new double[p];
        for (int j = 0; j < p; j++)
            for (int i = 0; i < n; i++) cross[j] += b * w[i] * x[i, j] * e[i];
        double[] theta = Solve(precision, cross);
        double s = 0.0;
        var fitted = new double[n];
        for (int i = 0; i < n; i++)
        {
            double f = 0.0;
            for (int j = 0; j < p; j++) f += x[i, j] * theta[j];
            s += b * w[i] * (e[i] - f) * (e[i] - f);
            fitted[i] = f + (flat ? 0.0 : level);
        }
        for (int j = flat ? 1 : 0; j < p; j++)
            for (int k = flat ? 1 : 0; k < p; k++) s += tau * d[j, k] * theta[j] * theta[k];
        int start = flat ? 1 : 0, q = p - start;
        var prior = new double[q, q];
        for (int j = 0; j < q; j++)
            for (int k = 0; k < q; k++) prior[j, k] = tau * d[j + start, k + start];
        double logDet = 0.5 * (q > 0 ? LogDet(prior) : 0.0) - 0.5 * LogDet(precision);
        var covariance = new double[p, p];
        for (int k = 0; k < p; k++)
        {
            var unit = new double[p];
            unit[k] = 1.0;
            double[] c = Solve(precision, unit);
            for (int j = 0; j < p; j++) covariance[j, k] = c[j];
        }
        return new Dense(s, logDet, fitted, covariance, x);
    }

    private static double[][] Columns(BassModel model, IEnumerable<BassTerm> terms) => [.. terms.Select(t => BassPointwiseTests.Column(model, t))];

    [Theory]
    [InlineData("g", GPriorCentre.Mean, WeightedPriorGram.Unweighted)]
    [InlineData("g", GPriorCentre.Zero, WeightedPriorGram.Unweighted)]
    [InlineData("g", GPriorCentre.Flat, WeightedPriorGram.Unweighted)]
    [InlineData("ridge", GPriorCentre.Mean, WeightedPriorGram.Unweighted)]
    [InlineData("ridge", GPriorCentre.Zero, WeightedPriorGram.Unweighted)]
    [InlineData("ridge", GPriorCentre.Flat, WeightedPriorGram.Unweighted)]
    [InlineData("g", GPriorCentre.Mean, WeightedPriorGram.Weighted)]
    [InlineData("g", GPriorCentre.Zero, WeightedPriorGram.Weighted)]
    [InlineData("g", GPriorCentre.Flat, WeightedPriorGram.Weighted)]
    [InlineData("ridge", GPriorCentre.Mean, WeightedPriorGram.Weighted)]
    [InlineData("ridge", GPriorCentre.Zero, WeightedPriorGram.Weighted)]
    [InlineData("ridge", GPriorCentre.Flat, WeightedPriorGram.Weighted)]
    public void TheBlock_MatchesDenseAlgebra(string form, GPriorCentre intercept, WeightedPriorGram gram)
    {
        bool ridge = form == "ridge";
        BassData data = BassSweepTests.Data(large: true);
        ReadOnlySpan<double> yValues = data.Response.Values;
        double[] y = yValues.ToArray();
        double[] w = Weights(y.Length);
        BassConfig config = Config(ridge ? new RidgeForm() : new GForm(), intercept, Sigma2Prior.InverseGamma(2.0, 0.3), 0.2);
        BassModel model = BassModel.Create(config, data, Weighted(w, WeightMeaning.Precision, gram));
        double sumLog = w.Sum(Math.Log);
        foreach (int count in new[] { 0, 1, 3, 5 })
        {
            BassTerm[] terms = Terms[..count];
            BassState state = BassMixing.Loaded(model, terms, Math.Log(0.2), out BassRefusal? refusal) ?? throw new InvalidOperationException($"{refusal}");
            double[][] z = Columns(model, state.Terms);
            foreach ((double tau, double b) in new[] { (0.2, 1.0), (0.2, 0.3), (0.013, 1.0), (3.0, 0.7) })
            {
                Dense dense = Reference(y, w, z, gram == WeightedPriorGram.Weighted, intercept, ridge, tau, b);
                GPriorMarginal marginal = model.Target.Marginal(state.Fit, tau, b);
                string at = $"M = {count}, τ = {tau}, b = {b}";
                Close(dense.S, marginal.S, 1e-9, $"S, {at}");
                Close(dense.LogDeterminant + 0.5 * b * sumLog, marginal.LogDeterminantTerm, 1e-9, $"the determinant term, {at}");
            }

            BassDesign design = model.NewDesign();
            Assert.Null(design.Load(terms));
            GPriorConditional conditional = model.Target.CoefficientBlock(design)(0.2, 0.8);
            Dense reference = Reference(y, w, z, gram == WeightedPriorGram.Weighted, intercept, ridge, 0.2, 0.8);
            // The fit from the conditional mean, in the design's coordinates, at each row: the columns centred under W.
            for (int i = 0; i < y.Length; i += 7)
            {
                double f = conditional.Mean[0];
                for (int j = 0; j < count; j++) f += conditional.Mean[j + 1] * (z[j][i] - z[j].Zip(w, (v, k) => v * k).Sum() / w.Sum());
                Close(reference.FittedMean[i], f, 1e-9, $"the mean's fit at row {i}, M = {count}");
            }
            // The variances of the fit at those rows, from the covariance's forms.
            var rows = new double[3, count + 1];
            for (int t = 0; t < 3; t++)
            {
                rows[t, 0] = 1.0;
                for (int j = 0; j < count; j++) rows[t, j + 1] = z[j][t * 11] - z[j].Zip(w, (v, k) => v * k).Sum() / w.Sum();
            }
            double[] variances = conditional.ConditionalVariance(model.Target.Sigma2, rows);
            double expected = conditional.Sigma2Posterior(model.Target.Sigma2).Mean;
            for (int t = 0; t < 3; t++)
            {
                double v = 0.0;
                for (int j = 0; j <= count; j++)
                    for (int k = 0; k <= count; k++) v += reference.X[t * 11, j] * reference.Covariance[j, k] * reference.X[t * 11, k];
                Close(expected * v, variances[t], 1e-8, $"the fit's variance at row {t * 11}, M = {count}");
            }
        }
    }

    [Theory]
    [InlineData("g", WeightedPriorGram.Unweighted)]
    [InlineData("ridge", WeightedPriorGram.Unweighted)]
    [InlineData("g", WeightedPriorGram.Weighted)]
    public void ColumnsCarryWeights_ThroughAProposalAsThroughALoad(string form, WeightedPriorGram gram)
    {
        BassData data = BassSweepTests.Data(large: true);
        double[] w = Weights(data.Response.Observations);
        BassConfig config = Config(form == "g" ? new GForm() : new RidgeForm(), GPriorCentre.Mean, Sigma2Prior.InverseGamma(2.0, 0.3), 0.2);
        BassModel model = BassModel.Create(config, data, Weighted(w, WeightMeaning.Precision, gram));
        BassDesign walked = model.NewDesign();
        foreach (BassTerm term in Terms[..4])
        {
            walked.ProposeBirth(term);
            walked.Accept();
        }
        void Same(DesignProposal p, IReadOnlyList<BassTerm> terms)
        {
            BassState loaded = BassMixing.Loaded(model, terms, Math.Log(0.2), out _)!;
            BassState proposed = BassMixing.Proposed(model, [.. terms], p, Math.Log(0.2), out _)!;
            foreach (double b in new[] { 1.0, 0.4 })
            {
                GPriorMarginal x = model.Target.Marginal(loaded.Fit, 0.2, b), y = model.Target.Marginal(proposed.Fit, 0.2, b);
                Close(x.S, y.S, 1e-10, "S");
                Close(x.LogDeterminantTerm, y.LogDeterminantTerm, 1e-10, "the determinant term");
            }
        }
        Same(walked.ProposeBirth(Terms[4]), Terms);
        Same(walked.ProposeDeath(1), [Terms[0], Terms[2], Terms[3]]);
        Same(walked.ProposeChange(2, Hinge(0, -1, 3)), [Terms[0], Terms[1], Terms[3], Hinge(0, -1, 3)]);
    }

    // A direct response's per-row weights, against explicit weighted least squares over the stacked rows, with a snapshot's
    // restored design continuing bit for bit.
    [Theory]
    [InlineData(WeightedPriorGram.Weighted)]
    [InlineData(WeightedPriorGram.Unweighted)]
    public void DirectFunctional_HoldsPerRowWeights(WeightedPriorGram gram)
    {
        BassData data = BassKernelEnumeration.FunctionalData();
        double[] y = data.Response.Values.ToArray();
        double[] w = Weights(y.Length, seed: 9);
        BassConfig config = Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.InverseGamma(2.0, 0.3), 0.2) with
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 5, Support = new SupportRule(SupportKind.Count, 1.0, 1.0) },
        };
        BassModel model = BassModel.Create(config, data, Weighted(w, WeightMeaning.Precision, gram));
        BassTerm[] terms =
        [
            new(continuous: [new HingeFactor(0, 1, 0)]),
            new(continuous: [new HingeFactor(0, 1, 1)], functional: [new HingeFactor(0, 1, 1)]),
            new(functional: [new HingeFactor(1, 1, 0)]),
        ];
        BassDesign design = model.NewDesign();
        foreach (BassTerm term in terms)
        {
            Assert.True(design.ProposeBirth(term).Admitted, $"{term}");
            design.Accept();
        }
        double[][] z = Columns(model, design.Terms);
        double mass = w.Sum(), mean = y.Zip(w, (v, k) => v * k).Sum() / mass;
        Dense dense = Reference(y, w, z, gram == WeightedPriorGram.Weighted, GPriorCentre.Mean, false, 0.2, 1.0);
        GPriorMarginal marginal = model.Target.Marginal(BassStructureFit.Of(design, 0.0), 0.2);
        Close(dense.S, marginal.S, 1e-9, "S");
        Close(dense.LogDeterminant + 0.5 * w.Sum(Math.Log), marginal.LogDeterminantTerm, 1e-9, "the determinant term");
        double squares = y.Select((v, i) => w[i] * (v - mean) * (v - mean)).Sum();
        Close(squares, design.ResponseSquares, 1e-11, "‖ỹ‖²_W");
        Close(squares, design.ResidualSquares + design.FittedSquares, 1e-10, "‖r‖²_W + Q");

        var restored = model.NewDesign();
        restored.Restore(design.Snapshot());
        var handed = model.NewDesign();
        handed.Restore(design.Snapshot(shareColumns: true));
        BassTerm next = new(continuous: [new HingeFactor(0, -1, 2)], functional: [new HingeFactor(1, -1, 1)]);
        DesignProposal a = design.ProposeBirth(next), b = restored.ProposeBirth(next), c = handed.ProposeBirth(next);
        Assert.Equal(a.ResidualSquares, b.ResidualSquares);
        Assert.Equal(a.FittedSquares, b.FittedSquares);
        Assert.Equal(a.ResidualSquares, c.ResidualSquares);
        design.Accept();
        restored.Accept();
        Assert.True(design.Residual.SequenceEqual(restored.Residual));
        Assert.Equal(design.PriorTerms is null, restored.PriorTerms is null);
    }

    // The gate reads the weighted columns: two columns that differ only where the weights nearly vanish are one column to the fit.
    [Fact]
    public void TheGate_ReadsTheWeights()
    {
        BassData data = BassSweepTests.Data(large: true);
        ReadOnlySpan<double> x0 = ((ContinuousPredictor)data.Predictors[0]).Values.Span;
        var w = new double[x0.Length];
        for (int i = 0; i < w.Length; i++) w[i] = x0[i] < 0.5 ? 1e-22 : 1.0;   // the two hinges below differ by a constant from ½ up
        BassTerm first = Hinge(0, 1, 1), second = Hinge(0, 1, 2);
        BassConfig config = Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.InverseGamma(2.0, 0.3), 0.2);
        BassDesign plain = BassModel.Create(config, data).NewDesign();
        Assert.True(plain.ProposeBirth(first).Admitted);
        plain.Accept();
        Assert.True(plain.ProposeBirth(second).Admitted);
        BassDesign weighted = BassModel.Create(config, data, Weighted(w, WeightMeaning.Precision, WeightedPriorGram.Weighted)).NewDesign();
        Assert.True(weighted.ProposeBirth(first).Admitted);
        weighted.Accept();
        Assert.Equal(BassRefusal.Numerical, weighted.ProposeBirth(second).Refusal);
    }

    // Near the gate a strict factorization decides, and an admitted candidate is accepted from it: its account and the state after
    // are a fresh load's, under each prior; just past the gate it is refused.
    [Theory]
    [InlineData("g", WeightedPriorGram.Weighted)]
    [InlineData("ridge", WeightedPriorGram.Weighted)]
    [InlineData("g", WeightedPriorGram.Unweighted)]
    [InlineData("ridge", WeightedPriorGram.Unweighted)]
    public void NearTheGate_TheStrictFactorizationDecides(string form, WeightedPriorGram gram)
    {
        BassData data = BassSweepTests.Data(large: true);
        double[] w = Weights(data.Response.Observations);
        Observation reading = Weighted(w, WeightMeaning.Precision, gram);
        BassConfig loose = Config(form == "g" ? new GForm() : new RidgeForm(), GPriorCentre.Mean, Sigma2Prior.InverseGamma(2.0, 0.3), 0.2);
        BassDesign Walked(BassConfig config)
        {
            BassDesign d = BassModel.Create(config, data, reading).NewDesign();
            foreach (BassTerm term in Terms[..2])
            {
                Assert.True(d.ProposeBirth(term).Admitted);
                d.Accept();
            }
            return d;
        }
        double angle = Walked(loose).ProposeBirth(Terms[2]).MinAngle;
        BassConfig Floor(double floor) => loose with { Structure = loose.Structure with { AngleFloor = floor } };

        BassConfig near = Floor(angle / 1.5);
        BassModel model = BassModel.Create(near, data, reading);
        BassDesign design = Walked(near);
        DesignProposal p = design.ProposeBirth(Terms[2]);
        Assert.True(p.Admitted);
        Close(angle, p.MinAngle, 1e-10, "the angle");
        BassState proposed = BassMixing.Proposed(model, Terms[..3], p, Math.Log(0.2), out _)!;
        BassState loaded = BassMixing.Loaded(model, Terms[..3], Math.Log(0.2), out _)!;
        Close(model.Target.Marginal(loaded.Fit, 0.2).S, model.Target.Marginal(proposed.Fit, 0.2).S, 1e-10, "S");
        design.Accept();
        BassDesign fresh = model.NewDesign();
        Assert.Null(fresh.Load(Terms[..3]));
        Close(fresh.ResidualSquares, design.ResidualSquares, 1e-10, "‖r‖² after the acceptance");
        Close(model.Target.Marginal(BassStructureFit.Of(fresh, 0.0), 0.2).LogDeterminantTerm, model.Target.Marginal(BassStructureFit.Of(design, 0.0), 0.2).LogDeterminantTerm, 1e-10, "the determinant term after it");

        Assert.Equal(BassRefusal.Numerical, Walked(Floor(angle * 1.5)).ProposeBirth(Terms[2]).Refusal);
    }

    // An exact fit under W: the move that explains it forms its residual explicitly, and the Jeffreys floor reads ‖·‖²_W.
    [Theory]
    [InlineData(WeightedPriorGram.Weighted)]
    [InlineData(WeightedPriorGram.Unweighted)]
    public void AnExactFit_ReadsTheWeightedNorm(WeightedPriorGram gram)
    {
        BassData data = BassSweepTests.Data(exact: true);
        double[] w = Weights(data.Response.Observations);
        BassTerm term = Hinge(0, 1, 1);
        foreach ((Sigma2Prior sigma2, BassRefusal? expected) in new (Sigma2Prior, BassRefusal?)[] { (Sigma2Prior.Jeffreys, BassRefusal.ExactFit), (Sigma2Prior.InverseGamma(2.0, 0.3), null), (Sigma2Prior.Held(0.5), null) })
        {
            BassModel model = BassModel.Create(Config(new GForm(), GPriorCentre.Mean, sigma2, 0.2), data, Weighted(w, WeightMeaning.Precision, gram));
            BassDesign design = model.NewDesign();
            DesignProposal p = design.ProposeBirth(term);
            Assert.True(p.ResidualSquares <= 1e-24 * design.ResponseSquares, $"‖r′‖² = {p.ResidualSquares:R} of ‖ỹ‖² = {design.ResponseSquares:R}");
            Assert.Equal(expected, model.Target.Refusal(BassStructureFit.Of(p, 0.0)));
        }
    }

    // The support threshold counts rows, or the rows weights stand for, by what the weight states.
    [Fact]
    public void TheSupportThreshold_CountsRowsAsTheWeightsState()
    {
        BassData data = BassSweepTests.Data(large: true);
        int n = data.Response.Observations;
        ReadOnlySpan<double> x0 = ((ContinuousPredictor)data.Predictors[0]).Values.Span;
        double[] w = [.. Enumerable.Repeat(3.0, n)];
        BassTerm term = Hinge(0, 1, 3);   // nonzero on the top level of x₀, a fifth of the rows
        int nonzero = 0;
        foreach (double v in x0)
            if (v > 0.75) nonzero++;
        BassConfig config = Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.InverseGamma(2.0, 0.3), 0.2) with
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 5, Support = new SupportRule(SupportKind.Count, nonzero + 1.0) },
        };
        Assert.Equal(BassRefusal.Support, BassModel.Create(config, data, Weighted(w, WeightMeaning.Precision, WeightedPriorGram.Weighted)).NewDesign().ProposeBirth(term).Refusal);
        Assert.True(BassModel.Create(config, data, Weighted(w, WeightMeaning.Replication, WeightedPriorGram.Weighted)).NewDesign().ProposeBirth(term).Admitted);
    }

    // ---- σ²'s integral ----

    // ln ∫ exp(ln p(y | σ²)) π(σ²) dσ², by the trapezoid on u = ln σ², which converges fast on a smooth integrand that vanishes at both ends.
    private static double Integrate(Func<double, double> logGiven, InverseGamma? prior)
    {
        const int points = 40_001;
        const double low = -40.0, high = 40.0;
        double step = (high - low) / (points - 1);
        var values = new double[points];
        double top = double.NegativeInfinity;
        for (int k = 0; k < points; k++)
        {
            double u = low + k * step, sigma2 = Math.Exp(u);
            values[k] = logGiven(sigma2) + (prior is { } p ? p.LogPdf(sigma2) + u : 0.0);   // Jeffreys: dσ²/σ² is du
            top = Math.Max(top, values[k]);
        }
        double sum = 0.0;
        foreach (double v in values) sum += Math.Exp(v - top);
        return top + Math.Log(sum * step);
    }

    [Theory]
    [InlineData("g", GPriorCentre.Mean, WeightedPriorGram.Weighted)]
    [InlineData("g", GPriorCentre.Zero, WeightedPriorGram.Weighted)]
    [InlineData("g", GPriorCentre.Flat, WeightedPriorGram.Weighted)]
    [InlineData("ridge", GPriorCentre.Mean, WeightedPriorGram.Weighted)]
    [InlineData("ridge", GPriorCentre.Flat, WeightedPriorGram.Weighted)]
    [InlineData("g", GPriorCentre.Mean, WeightedPriorGram.Unweighted)]
    [InlineData("g", GPriorCentre.Zero, WeightedPriorGram.Unweighted)]
    [InlineData("g", GPriorCentre.Flat, WeightedPriorGram.Unweighted)]
    [InlineData("ridge", GPriorCentre.Zero, WeightedPriorGram.Unweighted)]
    public void TheMarginalGivenSigma2_IntegratesToTheCollapsedMarginal(string form, GPriorCentre intercept, WeightedPriorGram gram)
    {
        BassData data = BassSweepTests.Data(large: true);
        double[] w = Weights(data.Response.Observations);
        foreach (Sigma2Prior sigma2 in new[] { Sigma2Prior.ScaledToResponse(0.5), Sigma2Prior.InverseGamma(2.0, 0.3), Sigma2Prior.Jeffreys })
        {
            BassConfig config = Config(form == "g" ? new GForm() : new RidgeForm(), intercept, sigma2, 0.2);
            BassModel model = BassModel.Create(config, data, Weighted(w, WeightMeaning.Precision, gram));
            foreach (int count in new[] { 0, 2, 5 })
            {
                BassState state = BassMixing.Loaded(model, Terms[..count], Math.Log(0.2), out _)!;
                foreach ((double tau, double b) in new[] { (0.2, 1.0), (0.2, 0.4), (0.03, 1.0) })
                {
                    GPriorMarginal m = model.Target.Marginal(state.Fit, tau, b);
                    double numeric = Integrate(m.LogMarginalGiven, model.Target.Sigma2);
                    Close(m.LogMarginal(model.Target.Sigma2), numeric, 1e-9, $"{sigma2}, M = {count}, τ = {tau}, b = {b}");

                    // Its slope in b, by a central difference of the form itself.
                    const double h = 1e-5, s2 = 0.37;
                    double up = model.Target.Marginal(state.Fit, tau, b - h).LogMarginalGiven(s2), down = model.Target.Marginal(state.Fit, tau, b - 2 * h).LogMarginalGiven(s2);
                    double top = model.Target.Marginal(state.Fit, tau, b).LogMarginalGiven(s2);
                    double slope = (3.0 * top - 4.0 * up + down) / (2.0 * h);   // one-sided: b may sit at 1
                    Close(slope, m.LogMarginalGivenSlope(s2), 1e-5, $"the slope, {sigma2}, M = {count}, τ = {tau}, b = {b}");
                }
            }
        }
    }

    [Fact]
    public void TheMarginalGivenSigma2_IsFiniteOnAnExactFit_AndZeroAtTheReferenceRung()
    {
        BassData data = BassSweepTests.Data(exact: true);
        BassConfig config = Config(new GForm(), GPriorCentre.Mean, Sigma2Prior.Held(0.5), 0.2);
        BassModel model = BassModel.Create(config, data);
        BassState state = BassMixing.Loaded(model, [Hinge(0, 1, 1)], Math.Log(0.2), out BassRefusal? refusal)!;
        Assert.Null(refusal);
        Assert.Null(model.Target.Refusal(state.Fit));
        BassEvaluation evaluation = model.Target.Evaluate(state);
        Assert.True(double.IsFinite(evaluation.LogLikelihood(1.0)));
        Assert.Equal(0.0, evaluation.LogLikelihood(0.0));
        Assert.True(model.Target.Defines(0.0));
        Assert.True(model.Target.HasCompleteValues);
    }

    // ---- The form's selection ----

    [Fact]
    public void AHeldSigma2_IsAPriorKind_ThatChangesNoOtherKindsSerializedForm()
    {
        BassConfig config = new();
        Assert.DoesNotContain("\"Value\"", config.ToJson());
        Assert.Contains("\"Value\":0.25", (config with { Sigma2 = Sigma2Prior.Held(0.25) }).ToJson());
        BassConfig back = BassConfig.FromJson((config with { Sigma2 = Sigma2Prior.Held(0.25) }).ToJson());
        Assert.Equal(Sigma2Prior.Held(0.25), back.Sigma2);
        Assert.Equal("Held(0.25)", back.Sigma2.ToString());
        Assert.True(back.Sigma2.IsProper);
        Assert.False(back.Sigma2.ReadsResponse);
        Assert.Null(back.Sigma2.For(1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sigma2Prior.Held(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sigma2Prior.Held(double.PositiveInfinity));
        Assert.Throws<ArgumentException>(() => new Sigma2Prior(Sigma2PriorKind.Jeffreys, value: 1.0));
        Assert.Throws<ArgumentException>(() => new Sigma2Prior(Sigma2PriorKind.Held, 1.0, value: 1.0));
        Assert.Throws<ArgumentException>(() => new WeightedNormalModel(sigma2Prior: Sigma2Prior.Held(1.0)));

        // Given σ², every b is defined but the flat intercept's pole.
        Assert.True(BassTarget.Defines(TemperedForm.Integrated, GPriorCentre.Mean, null, 10, 0.0, givenSigma2: true));
        Assert.False(BassTarget.Defines(TemperedForm.Integrated, GPriorCentre.Mean, null, 10, 0.0));
        Assert.False(BassTarget.Defines(TemperedForm.Integrated, GPriorCentre.Flat, null, 10, 0.0, givenSigma2: true));
        Assert.True(BassTarget.Defines(TemperedForm.Integrated, GPriorCentre.Flat, null, 10, 0.01, givenSigma2: true));
    }

    [Fact]
    public void AHeldSigma2_LiftsTheResponsesRefusals_AndARunRetainsIt()
    {
        BassData flat = new([new ContinuousPredictor("x", [.. Enumerable.Range(0, 20).Select(i => i / 19.0)])], new ScalarResponse([.. Enumerable.Repeat(2.0, 20)]));
        Assert.Throws<ArgumentException>(() => BassResolved.Resolve(BassSweepTests.Config(), flat));
        BassResolved resolved = BassResolved.Resolve(BassSweepTests.Config() with { Sigma2 = Sigma2Prior.Held(1.0) }, flat);
        Assert.Equal(1.0, resolved.GivenSigma2);
        Assert.Null(resolved.Sigma2);

        BassData data = BassSweepTests.Data();
        BassConfig config = BassSweepTests.Config() with
        {
            Sigma2 = Sigma2Prior.Held(0.09),
            Sampling = new BassSampling { Chains = 2, BurnIn = 50, Samples = 20, Moves = new EqualOverAvailable() },
        };
        BassFit fit = Bass.Fit(config, data);
        foreach (DrawStore<BassDraw> chain in fit.Draws)
            foreach (BassDraw draw in chain) Assert.Equal(0.09, draw.Sigma2);
    }

    // σ² from its collapsed conditional given the structure, and the chain handed its target given that σ².
    private sealed class Sigma2Block(BassTarget target) : IBlock<ReversibleJumpChain<BassState>>
    {
        public string Name => "sigma2";

        public bool Update(ReversibleJumpChain<BassState> state, Xoshiro256PlusPlus rng)
        {
            double tau = target.HeldTau!.Value;
            double sigma2 = target.Marginal(state.Current.Fit, tau).Sigma2Posterior(target.Sigma2).Sample(rng);
            state.ReplaceTarget(target.Given(sigma2));
            return true;
        }
    }

    // ---- Enumeration ----

    private static void CheckKernel(BassModel model, string label, ITestOutputHelper log)
    {
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        var flows = BassKernelEnumeration.Kernel(space, (key, s) => BassKernelEnumeration.D01(space, key, s, null));
        (Dictionary<string, int> transitions, double worst) = BassKernelEnumeration.Check(space, flows);
        log.WriteLine($"{label}: {space.States.Count} states, {transitions.Values.Sum()} transitions, worst {worst:E2}");
        Assert.True(worst <= BassKernelEnumeration.Tolerance, $"{label}: worst {worst:E2}");
    }

    // π(S) by dense algebra at the held τ: the term measure, plus the marginal from the columns alone.
    private static double[] DensePosterior(BassKernelEnumeration.Space space, double[] y, double[] w, bool priorOnWeighted, double tau, double? given, InverseGamma? prior, double sumLog)
    {
        BassModel model = space.Model;
        var logs = new double[space.States.Count];
        for (int s = 0; s < logs.Length; s++)
        {
            BassKernelEnumeration.State state = space.States[s];
            Dense dense = Reference(y, w, Columns(model, state.Terms), priorOnWeighted, GPriorCentre.Mean, false, tau, 1.0);
            double logPrior = model.Prior.LogCollectionWeight(state.Terms.Length) + state.Terms.Sum(model.Prior.LogPrior);
            double common = dense.LogDeterminant + 0.5 * sumLog;
            int n = y.Length;
            double likelihood = given is double sigma2
                ? -0.5 * n * Math.Log(2.0 * Math.PI * sigma2) + common - 0.5 * dense.S / sigma2
                : -0.5 * n * Math.Log(2.0 * Math.PI) + common + prior!.Value.Shape * Math.Log(prior.Value.Scale)
                  + SpecialFunctions.LogGamma(prior.Value.Shape + 0.5 * n) - SpecialFunctions.LogGamma(prior.Value.Shape)
                  - (prior.Value.Shape + 0.5 * n) * Math.Log(prior.Value.Scale + 0.5 * dense.S);
            logs[s] = logPrior + likelihood;
            Close(logs[s], state.LogPi, 1e-9, $"state {s}: {string.Join(", ", state.Terms.Select(t => t.ToString()))}");
        }
        double top = logs.Max();
        double[] p = [.. logs.Select(v => Math.Exp(v - top))];
        double total = p.Sum();
        return [.. p.Select(v => v / total)];
    }

    // A chain's visits against π: by M and on the most probable structures, within the tolerance its length earns.
    private static void CheckVisits(BassKernelEnumeration.Space space, Func<int, BassChain> start, ITestOutputHelper log, string label, Action<BassChain>? each = null)
    {
        double[] pi = space.Posterior();
        const int chains = 4, burn = 2_000, sweeps = 30_000;
        var visits = new double[pi.Length];
        for (int c = 0; c < chains; c++)
        {
            BassChain chain = start(c);
            for (int s = 0; s < burn; s++) chain.Step();
            for (int s = 0; s < sweeps; s++)
            {
                chain.Step();
                int found = space.Find(chain.Current.Terms);
                Assert.True(found >= 0);
                visits[found] += 1.0;
            }
        }
        double total = chains * sweeps;
        int cap = space.States.Max(s => s.Terms.Length);
        var byCountPi = new double[cap + 1];
        var byCountVisits = new double[cap + 1];
        for (int i = 0; i < pi.Length; i++)
        {
            byCountPi[space.States[i].Terms.Length] += pi[i];
            byCountVisits[space.States[i].Terms.Length] += visits[i] / total;
        }
        for (int m = 0; m <= cap; m++)
            Assert.True(Math.Abs(byCountPi[m] - byCountVisits[m]) <= 0.03, $"{label}: M = {m}: π {byCountPi[m]:F4}, visits {byCountVisits[m]:F4}");
        foreach (int i in Enumerable.Range(0, pi.Length).OrderByDescending(i => pi[i]).Take(10))
            Assert.True(Math.Abs(pi[i] - visits[i] / total) <= 0.02, $"{label}: state {i}: π {pi[i]:F4}, visits {visits[i] / total:F4}");
        log.WriteLine($"{label}: M's law {string.Join(", ", byCountPi.Select(v => v.ToString("F3")))} against {string.Join(", ", byCountVisits.Select(v => v.ToString("F3")))}");
    }

    private static BassConfig EnumerationConfig(Sigma2Prior sigma2, double tau = 0.1)
        => BassSweepTests.Config(maxBasis: 3) with { Sigma2 = sigma2, Coefficients = BassCoefficientPrior.G(TauPrior.Held(tau)) };

    [Theory]
    [InlineData(WeightedPriorGram.Weighted)]
    [InlineData(WeightedPriorGram.Unweighted)]
    public void Enumeration_AWeightedTarget(WeightedPriorGram gram)
    {
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        double[] w = Weights(y.Length, seed: 21);
        BassConfig config = EnumerationConfig(Sigma2Prior.InverseGamma(2.0, 0.3));
        BassModel model = BassModel.Create(config, data, Weighted(w, WeightMeaning.Precision, gram));
        CheckKernel(model, $"weighted, prior on the {gram} columns", _out);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        DensePosterior(space, y, w, gram == WeightedPriorGram.Weighted, 0.1, null, model.Target.Sigma2, w.Sum(Math.Log));
        CheckVisits(space, c => new BassChain(model, [], Math.Log(0.1), new Xoshiro256PlusPlus(400 + c), new Xoshiro256PlusPlus(500 + c)), _out, $"weighted, {gram}");
    }

    [Theory]
    [InlineData(false, WeightedPriorGram.Weighted)]
    [InlineData(true, WeightedPriorGram.Weighted)]
    [InlineData(true, WeightedPriorGram.Unweighted)]
    public void Enumeration_ATargetGivenSigma2(bool weighted, WeightedPriorGram gram)
    {
        const double sigma2 = 0.25;
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        double[] w = weighted ? Weights(y.Length, seed: 21) : [.. Enumerable.Repeat(1.0, y.Length)];
        BassConfig config = EnumerationConfig(Sigma2Prior.Held(sigma2));
        BassModel model = weighted ? BassModel.Create(config, data, Weighted(w, WeightMeaning.Precision, gram)) : BassModel.Create(config, data);
        CheckKernel(model, $"given σ², {(weighted ? "weighted, " + gram : "unit weights")}", _out);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        DensePosterior(space, y, w, gram == WeightedPriorGram.Weighted, 0.1, sigma2, null, w.Sum(Math.Log));
        CheckVisits(space, c => new BassChain(model, [], Math.Log(0.1), new Xoshiro256PlusPlus(600 + c), new Xoshiro256PlusPlus(700 + c)), _out, "given σ²");
    }

    // σ² drawn apart, the target given it, at the sweep's Scale slot: over the structures the walk is the collapsed posterior's,
    // σ² integrated by its own draws. The form E1's σ² block and E32's outside σ² read.
    [Fact]
    public void Sigma2DrawnApart_WalksTheCollapsedPosterior()
    {
        BassData data = BassSweepTests.Data();
        BassConfig config = EnumerationConfig(Sigma2Prior.InverseGamma(2.0, 0.3));
        BassModel model = BassModel.Create(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        double[] pi = space.Posterior();
        const int chains = 4, burn = 2_000, sweeps = 30_000;
        var visits = new double[pi.Length];
        for (int c = 0; c < chains; c++)
        {
            var moves = new BassStructureMoves(model, [], Math.Log(0.1));
            var chain = new ReversibleJumpChain<BassState>(moves.Moves, model.Target.Given(0.1), moves.Start, new Xoshiro256PlusPlus(800 + c));
            var scan = new BlockScan<ReversibleJumpChain<BassState>>([new MoveStep<BassState>(), new Sigma2Block(model.Target)], chain, new Xoshiro256PlusPlus(900 + c), x => x);
            for (int s = 0; s < burn; s++) scan.Step();
            for (int s = 0; s < sweeps; s++)
            {
                scan.Step();
                visits[space.Find(chain.Current.Terms)] += 1.0;
            }
        }
        double total = chains * sweeps;
        int cap = space.States.Max(s => s.Terms.Length);
        var byCountPi = new double[cap + 1];
        var byCountVisits = new double[cap + 1];
        for (int i = 0; i < pi.Length; i++)
        {
            byCountPi[space.States[i].Terms.Length] += pi[i];
            byCountVisits[space.States[i].Terms.Length] += visits[i] / total;
        }
        _out.WriteLine($"M's law {string.Join(", ", byCountPi.Select(v => v.ToString("F3")))} against {string.Join(", ", byCountVisits.Select(v => v.ToString("F3")))}");
        for (int m = 0; m <= cap; m++)
            Assert.True(Math.Abs(byCountPi[m] - byCountVisits[m]) <= 0.03, $"M = {m}: π {byCountPi[m]:F4}, visits {byCountVisits[m]:F4}");
        foreach (int i in Enumerable.Range(0, pi.Length).OrderByDescending(i => pi[i]).Take(10))
            Assert.True(Math.Abs(pi[i] - visits[i] / total) <= 0.02, $"state {i}: π {pi[i]:F4}, visits {visits[i] / total:F4}");
    }

    // ---- What reads the design ----

    // A retained density under precision weights is N(y | f, σ²/w); under weights that count rows, w copies of the row's. The
    // collapsed predictive of one row against the joint Gaussian, given σ².
    [Theory]
    [InlineData(WeightedPriorGram.Weighted)]
    [InlineData(WeightedPriorGram.Unweighted)]
    public void ThePointwiseDensities_ReadTheWeights(WeightedPriorGram gram)
    {
        BassData data = BassSweepTests.Data();
        double[] y = data.Response.Values.ToArray();
        int n = y.Length;
        double[] w = Weights(n, seed: 21);
        const double sigma2 = 0.25, tau = 0.1;
        BassConfig config = EnumerationConfig(Sigma2Prior.Held(sigma2), tau);
        BassModel model = BassModel.Create(config, data, Weighted(w, WeightMeaning.Precision, gram));
        var chain = new BassChain(model, [], Math.Log(tau), new Xoshiro256PlusPlus(31), new Xoshiro256PlusPlus(32));
        var table = new BassTermTable();
        var composer = new Xoshiro256PlusPlus(33);
        for (int s = 0; s < 200; s++) chain.Step();
        BassDraw draw = chain.Compose(composer, table);
        Assert.Equal(sigma2, draw.Sigma2);

        PointwiseLogDensity retained = BassPointwise.LogDensity(model, table, [[draw]], null, PointwiseContent.RetainedParameters);
        BassTerm[] terms = [.. draw.Structure.TermIds.Select(id => table[id])];
        double[][] z = Columns(model, terms);
        for (int i = 0; i < n; i++)
        {
            double f = draw.Coefficients[0];
            for (int j = 0; j < z.Length; j++) f += draw.Coefficients[j + 1] * z[j][i];
            double expected = -0.5 * Math.Log(2.0 * Math.PI * sigma2 / w[i]) - 0.5 * w[i] * (y[i] - f) * (y[i] - f) / sigma2;
            Close(expected, retained.Values[0, i], 1e-11, $"row {i}");
        }

        PointwiseLogDensity collapsed = BassPointwise.LogDensity(model, table, [[draw]], null, PointwiseContent.CollapsedPredictive);
        // y | σ² ~ N(Xμ₀, σ²(W⁻¹ + XD⁻¹Xᵀ/τ)), with D the prior's in its own centring; a row's predictive by conditioning.
        bool onWeighted = gram == WeightedPriorGram.Weighted;
        double[] pw = onWeighted ? w : [.. Enumerable.Repeat(1.0, n)];
        Dense dense = Reference(y, pw, z, onWeighted, GPriorCentre.Mean, false, tau, 1.0);
        double level = y.Zip(pw, (v, k) => v * k).Sum() / pw.Sum();
        int p = z.Length + 1;
        var d = new double[p, p];
        for (int j = 0; j < p; j++)
            for (int k = 0; k < p; k++)
            {
                double own = 0.0;
                for (int i = 0; i < n; i++) own += pw[i] * dense.X[i, j] * dense.X[i, k];
                d[j, k] = own;
            }
        for (int k = 1; k < p; k++) d[0, k] = d[k, 0] = 0.0;
        var covariance = new double[n, n];
        for (int k = 0; k < p; k++)
        {
            var unit = new double[p];
            unit[k] = 1.0;
            double[] column = Solve(d, unit);
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++) covariance[a, b] += dense.X[a, k] * (Enumerable.Range(0, p).Sum(j => column[j] * dense.X[b, j])) / tau;
        }
        for (int a = 0; a < n; a++) covariance[a, a] += 1.0 / w[a];
        double[] e = [.. y.Select(v => v - level)];
        foreach (int i in new[] { 0, 5, 11, 17 })
        {
            int[] rest = [.. Enumerable.Range(0, n).Where(r => r != i)];
            var block = new double[n - 1, n - 1];
            var cross = new double[n - 1];
            var eRest = new double[n - 1];
            for (int a = 0; a < n - 1; a++)
            {
                cross[a] = covariance[i, rest[a]];
                eRest[a] = e[rest[a]];
                for (int b = 0; b < n - 1; b++) block[a, b] = covariance[rest[a], rest[b]];
            }
            double[] solved = Solve(block, cross);
            double mean = solved.Zip(eRest, (s, v) => s * v).Sum(), variance = covariance[i, i] - solved.Zip(cross, (s, c) => s * c).Sum();
            double expected = -0.5 * Math.Log(2.0 * Math.PI * sigma2 * variance) - 0.5 * (e[i] - mean) * (e[i] - mean) / (sigma2 * variance);
            Close(expected, collapsed.Values[0, i], 1e-9, $"the predictive of row {i}");
        }
    }

    [Fact]
    public void ThePointwiseDensities_UnderWeightsThatCountRows_AreTheirCopies()
    {
        BassData data = BassSweepTests.Data();
        int[] counts = Counts(data.Response.Observations);
        BassData repeated = Duplicated(data, counts);
        BassConfig config = EnumerationConfig(Sigma2Prior.InverseGamma(2.0, 0.3));
        BassModel a = BassModel.Create(config, data, Weighted(AsDoubles(counts), WeightMeaning.Replication, WeightedPriorGram.Weighted));
        BassModel b = BassModel.Create(config, repeated);
        var chain = new BassChain(a, [], Math.Log(0.1), new Xoshiro256PlusPlus(41), new Xoshiro256PlusPlus(42));
        var table = new BassTermTable();
        for (int s = 0; s < 200; s++) chain.Step();
        BassDraw draw = chain.Compose(new Xoshiro256PlusPlus(43), table);
        PointwiseLogDensity weighted = BassPointwise.LogDensity(a, table, [[draw]], null, PointwiseContent.RetainedParameters);
        PointwiseLogDensity plain = BassPointwise.LogDensity(b, table, [[draw]], null, PointwiseContent.RetainedParameters);
        int row = 0;
        for (int i = 0; i < counts.Length; i++)
        {
            double sum = 0.0;
            for (int c = 0; c < counts[i]; c++) sum += plain.Values[0, row++];
            Close(sum, weighted.Values[0, i], 1e-11, $"row {i}");
        }
        Assert.Throws<NotSupportedException>(() => BassPointwise.LogDensity(a, table, [[draw]], null, PointwiseContent.CollapsedPredictive));
    }

    // A ladder's swaps hand a replica's design over, its pencil and prior side included: each replica's state then reads as a
    // fresh load of its structure does.
    [Theory]
    [InlineData(WeightedPriorGram.Unweighted)]
    [InlineData(WeightedPriorGram.Weighted)]
    public void ALadder_HandsItsWeightedDesignsOver(WeightedPriorGram gram)
    {
        BassData data = BassSweepTests.Data(large: true);
        double[] w = Weights(data.Response.Observations);
        BassConfig config = BassSweepTests.Config(maxBasis: 5);
        BassModel model = BassModel.Create(config, data, Weighted(w, WeightMeaning.Precision, gram));
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(1)).State;
        var ladder = new BassLadder(model, new TemperingSpec { Levels = 4, HotEnd = 0.1 }, start.Terms, start.LogTau, 2);
        for (int step = 0; step < 150; step++)
        {
            ladder.Step();
            BassState cold = ladder.Current;
            BassState loaded = BassMixing.Loaded(model, cold.Terms, cold.LogTau, out BassRefusal? refusal)!;
            Assert.Null(refusal);
            GPriorMarginal a = model.Target.Marginal(cold.Fit, cold.Tau), b = model.Target.Marginal(loaded.Fit, cold.Tau);
            Close(b.S, a.S, 1e-8, $"step {step}: S");
            Close(b.LogDeterminantTerm, a.LogDeterminantTerm, 1e-8, $"step {step}: the determinant term");
        }
    }

    // A held σ² keeps complete values, so a tempered run down to the reference rung reports an evidence, against the enumeration's.
    [Fact]
    public void AHeldSigma2_ReportsAnEvidence_ThatMatchesTheEnumeration()
    {
        BassData data = BassSweepTests.Data();
        BassConfig config = EnumerationConfig(Sigma2Prior.Held(0.25)) with
        {
            Sampling = new BassSampling
            {
                Chains = 2, BurnIn = 1_500, Samples = 3_000, Moves = new EqualOverAvailable(), ReportEvidence = true,
                Tempering = new TemperingSpec { Levels = 8, HotEnd = 0.05, ReferenceRung = true },
            },
        };
        BassFit fit = Bass.Fit(config, data);
        BassEvidence evidence = fit.Evidence ?? throw new InvalidOperationException("The run reports no evidence.");
        Assert.Equal(Sigma2Prior.Held(0.25), evidence.Priors.Sigma2Prior);
        Assert.Null(evidence.Priors.Sigma2);
        BassModel model = BassModel.Create(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        BassKernelEnumeration.Space prior = BassKernelEnumeration.Enumerate(model, 0.0);
        double top = space.States.Max(s => s.LogPi), topPrior = prior.States.Max(s => s.LogPi);
        double exact = top + Math.Log(space.States.Sum(s => Math.Exp(s.LogPi - top))) - topPrior - Math.Log(prior.States.Sum(s => Math.Exp(s.LogPi - topPrior)));
        _out.WriteLine($"ln Z {evidence.LogEvidence:F3}, the enumeration's {exact:F3}");
        Assert.True(Math.Abs(evidence.LogEvidence - exact) <= 0.5, $"ln Z {evidence.LogEvidence:F3}, the enumeration's {exact:F3}");
    }

    // Under weights the conditional mean's surface reads the centring under W, so a draw's uncentred coefficients are its own.
    [Fact]
    public void AComposedDraw_UncentresUnderW()
    {
        BassData data = BassSweepTests.Data();
        double[] w = Weights(data.Response.Observations, seed: 21);
        BassConfig config = EnumerationConfig(Sigma2Prior.Held(1e-30));   // σ² ≈ 0: each draw is its conditional mean
        foreach (WeightedPriorGram gram in new[] { WeightedPriorGram.Weighted, WeightedPriorGram.Unweighted })
        {
            BassModel model = BassModel.Create(config, data, Weighted(w, WeightMeaning.Precision, gram));
            var chain = new BassChain(model, [Terms[0], Terms[1]], Math.Log(0.1), new Xoshiro256PlusPlus(51), new Xoshiro256PlusPlus(52));
            var table = new BassTermTable();
            BassDraw draw = chain.Compose(new Xoshiro256PlusPlus(53), table);
            BassTerm[] terms = [.. draw.Structure.TermIds.Select(id => table[id])];
            double[][] z = Columns(model, terms);
            double[] y = data.Response.Values.ToArray();
            Dense dense = Reference(y, w, z, gram == WeightedPriorGram.Weighted, GPriorCentre.Mean, false, 0.1, 1.0);
            for (int i = 0; i < y.Length; i += 3)
            {
                double f = draw.Coefficients[0];
                for (int j = 0; j < z.Length; j++) f += draw.Coefficients[j + 1] * z[j][i];
                Close(dense.FittedMean[i], f, 1e-9, $"{gram}: the fit at row {i}");
            }
        }
    }
}
