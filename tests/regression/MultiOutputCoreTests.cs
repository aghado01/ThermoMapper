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

/// <summary>
/// E34's core (the multi-output brief, slice 1): K outputs on one basis, Σ integrated.
/// <list type="bullet">
/// <item>The design's K-output side against the scalar engine: a random move sequence mirrored on one scalar design per output and one
/// per pair of outputs, whose sums of squares give the diagonal of the K × K forms and, by polarization, the rest; every proposal and
/// every accepted state.</item>
/// <item>The target against the matrix-t reference at random structures, τ and b, under the three intercept conventions.</item>
/// <item>K = 1 against scalar BASS: the resolution's s² and prior, the design's statistics bit for bit, the target to rounding, and the
/// draw of Σ and the coefficients to rounding on the same stream.</item>
/// <item>Equivariance under decision 2: Y → YL shifts every structure's marginal by −N ln|det L| at every τ.</item>
/// <item>The refusals of what later slices bring.</item>
/// </list>
/// </summary>
public sealed class MultiOutputCoreTests(ITestOutputHelper output)
{
    private const int Seed = 20261110;

    // ── Data ───────────────────────────────────────────────────────────────────────────────────────────────────────────────

    // n cases of two continuous inputs and a categorical one; K outputs from surfaces and noise outside the model, the noise correlated.
    internal static (List<BassPredictor> Predictors, double[,] Y) Data(int n, int k, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x1 = [.. Enumerable.Range(0, n).Select(_ => Math.Round(rng.NextDouble() * 30) / 30)];
        double[] x2 = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        int[] z = [.. Enumerable.Range(0, n).Select(i => i % 3)];
        var y = new double[n, k];
        for (int i = 0; i < n; i++)
        {
            double shared = 0.2 * Normal.SampleStandard(rng);
            for (int l = 0; l < k; l++)
            {
                double surface = l switch
                {
                    0 => Math.Sin(3.0 * x1[i]) + x1[i] * x2[i],
                    1 => 2.0 * Math.Exp(-4.0 * (x1[i] - 0.5) * (x1[i] - 0.5)) - x2[i],
                    2 => Math.Abs(x2[i] - 0.3) + (z[i] == 1 ? 0.5 : 0.0),
                    _ => Math.Cos(2.0 * l * x2[i]) * x1[i],
                };
                y[i, l] = (l + 1) * (surface + shared * (l % 2 == 0 ? 1.0 : -0.6) + 0.1 * Normal.SampleStandard(rng)) + 3.0 * l;
            }
        }
        return ([new ContinuousPredictor("x1", x1), new ContinuousPredictor("x2", x2), new CategoricalPredictor("z", z, ["a", "b", "c"])], y);
    }

    private static string[] Names(int k) => [.. Enumerable.Range(0, k).Select(l => $"y{l}")];

    private static BassModel Model(List<BassPredictor> predictors, BassResponse response, GPriorCentre intercept = GPriorCentre.Mean, int maxBasis = 12)
        => BassModel.Create(new BassConfig { Intercept = intercept, Structure = new BassStructureSpec { MaxBasis = maxBasis } }, new BassData(predictors, response));

    private static double[] Column(double[,] y, int l) => [.. Enumerable.Range(0, y.GetLength(0)).Select(i => y[i, l])];

    private static BassTerm RandomTerm(BassBasis basis, Xoshiro256PlusPlus rng)
    {
        BassInputs inputs = basis.Inputs;
        while (true)
        {
            var hinges = new List<HingeFactor>();
            for (int v = 0; v < inputs.Cases.ContinuousCount; v++)
                if (rng.NextDouble() < 0.45)
                    hinges.Add(new HingeFactor(v, rng.NextDouble() < 0.5 ? -1 : 1, rng.NextInt(inputs.ContinuousKnots(v).Length)));
            var subsets = new List<CategoricalFactor>();
            if (rng.NextDouble() < 0.3) subsets.Add(new CategoricalFactor(0, [rng.NextInt(3)]));
            if (hinges.Count + subsets.Count > 0) return new BassTerm(hinges, subsets, []);
        }
    }

    // ── The design against the scalar engine ───────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Design_FollowsTheScalarDesigns_OutputByOutput_AndInPairs()
    {
        const int n = 80, k = 3;
        (List<BassPredictor> predictors, double[,] y) = Data(n, k, Seed);
        BassModel multi = Model(predictors, new MultiOutputResponse(y, Names(k)));
        BassDesign design = multi.NewDesign();
        var singles = new BassDesign[k];
        for (int l = 0; l < k; l++) singles[l] = Model(predictors, new ScalarResponse(Column(y, l))).NewDesign();
        var pairs = new Dictionary<(int, int), BassDesign>();
        for (int a = 0; a < k; a++)
            for (int d = 0; d < a; d++)
            {
                double[] sum = [.. Enumerable.Range(0, n).Select(i => y[i, a] + y[i, d])];
                pairs[(a, d)] = Model(predictors, new ScalarResponse(sum)).NewDesign();
            }
        IEnumerable<BassDesign> scalars = singles.Concat(pairs.Values);

        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 1)[0]);
        double worst = 0.0;
        int proposals = 0, accepted = 0;
        for (int step = 0; step < 300; step++)
        {
            int m = design.Count;
            double u = rng.NextDouble();
            BassTerm? term = null;
            int index = -1;
            if (m == 0 || u < 0.45) term = RandomTerm(multi.Basis, rng);
            else if (u < 0.7) index = rng.NextInt(m);
            else { index = rng.NextInt(m); term = RandomTerm(multi.Basis, rng); }

            DesignProposal Propose(BassDesign d) => term is null ? d.ProposeDeath(index) : index < 0 ? d.ProposeBirth(term) : d.ProposeChange(index, term);
            DesignProposal several = Propose(design);
            DesignProposal[] one = [.. singles.Select(Propose)];
            var two = pairs.ToDictionary(p => p.Key, p => Propose(p.Value));
            Assert.All(one.Concat(two.Values), p => Assert.Equal(several.Refusal, p.Refusal));
            if (!several.Admitted) continue;
            proposals++;
            worst = Math.Max(worst, Gap(several.Outputs!, one, two, k));
            if (rng.NextDouble() < 0.6)
            {
                design.Accept();
                foreach (BassDesign d in scalars) d.Accept();
                accepted++;
                OutputForms current = design.OutputForms!;
                DesignProposal[] heldOne = [.. singles.Select(d => new DesignProposal(null, d.Count, d.ResidualSquares, d.FittedSquares, 0, 0))];
                var heldTwo = pairs.ToDictionary(p => p.Key, p => new DesignProposal(null, p.Value.Count, p.Value.ResidualSquares, p.Value.FittedSquares, 0, 0));
                worst = Math.Max(worst, Gap(current, heldOne, heldTwo, k));
            }
        }
        output.WriteLine($"{proposals} admitted proposals, {accepted} accepted; worst relative gap {worst:E2}");
        Assert.True(accepted > 30, $"only {accepted} moves accepted");
        Assert.True(worst < 1e-9, $"worst relative gap {worst:E2}");
    }

    // The forms' gap from the scalar designs': R_kk and F_kk against output k's, R_kl and F_kl against (pair − R_kk − R_ll)/2, each relative to
    // the pair's own scale.
    private static double Gap(OutputForms forms, DesignProposal[] one, Dictionary<(int, int), DesignProposal> two, int k)
    {
        double worst = 0.0;
        for (int a = 0; a < k; a++)
        {
            worst = Math.Max(worst, Relative(one[a].ResidualSquares, forms.Residual[a, a], one[a].ResidualSquares + one[a].FittedSquares));
            worst = Math.Max(worst, Relative(one[a].FittedSquares, forms.Fitted[a, a], one[a].ResidualSquares + one[a].FittedSquares));
            for (int d = 0; d < a; d++)
            {
                DesignProposal p = two[(a, d)];
                double scale = one[a].ResidualSquares + one[a].FittedSquares + one[d].ResidualSquares + one[d].FittedSquares;
                double r = 0.5 * (p.ResidualSquares - one[a].ResidualSquares - one[d].ResidualSquares);
                double f = 0.5 * (p.FittedSquares - one[a].FittedSquares - one[d].FittedSquares);
                worst = Math.Max(worst, Relative(r, forms.Residual[a, d], scale));
                worst = Math.Max(worst, Relative(f, forms.Fitted[a, d], scale));
                Assert.Equal(forms.Residual[a, d], forms.Residual[d, a]);
                Assert.Equal(forms.Fitted[a, d], forms.Fitted[d, a]);
            }
        }
        return worst;
    }

    private static double Relative(double expected, double actual, double scale) => Math.Abs(actual - expected) / Math.Max(scale, 1e-300);

    // ── The target against the matrix-t reference ─────────────────────────────────────────────────────────────────────────

    public static TheoryData<GPriorCentre> Centres => new() { GPriorCentre.Mean, GPriorCentre.Zero, GPriorCentre.Flat };

    [Theory]
    [MemberData(nameof(Centres))]
    public void Target_IsTheMatrixTReference_AtRandomStructures(GPriorCentre centre)
    {
        const int n = 60, k = 3;
        (List<BassPredictor> predictors, double[,] y) = Data(n, k, Seed + 1);
        BassModel model = Model(predictors, new MultiOutputResponse(y, Names(k)), centre);
        BassDesign design = model.NewDesign();
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 2)[1]);
        double[] levels = centre == GPriorCentre.Zero ? new double[k] : [.. Enumerable.Range(0, k).Select(l => Column(y, l).Average())];
        double worst = 0.0;
        int compared = 0;
        for (int index = 0; index < 150; index++)
        {
            int m = rng.NextInt(6);
            var terms = new List<BassTerm>();
            while (terms.Count < m)
            {
                BassTerm t = RandomTerm(model.Basis, rng);
                if (!terms.Contains(t)) terms.Add(t);
            }
            if (design.Load(terms) is not null) continue;
            BassStructureFit fit = BassStructureFit.Of(design, 0.0);
            double logTau = Math.Log(Math.Pow(10.0, -3.0 + 5.0 * rng.NextDouble())), b = rng.NextDouble() < 0.5 ? 1.0 : 0.05 + 0.95 * rng.NextDouble();
            double engine = model.Target.Evaluate(fit, logTau).LogLikelihood(b);

            double[,] xc = Columns(model, design);
            ReferenceSigmaPrior prior = ReferenceSigmaPrior.ScaledToResponse(y, levels);
            double reference = MatrixTReference.LogMarginal(MatrixTReference.Collapse(MultiOutputCase.GForm(y, xc, centre, Math.Exp(logTau), b)), prior);
            worst = Math.Max(worst, Math.Abs(engine - reference) / Math.Max(1.0, Math.Abs(reference)));
            compared++;
        }
        output.WriteLine($"{centre}: {compared} structures, worst relative gap {worst:E2}");
        Assert.True(compared > 40, $"only {compared} structures loaded");
        Assert.True(worst < 1e-10, $"worst relative gap {worst:E2}");
    }

    // The design's terms' columns as evaluated, centred over the cases: X_c, which the reference reads.
    internal static double[,] Columns(BassModel model, BassDesign design)
    {
        int n = model.Data.Response.Cases, m = design.Count;
        var xc = new double[n, m];
        var column = new double[n];
        for (int j = 0; j < m; j++)
        {
            model.Basis.CaseColumn(design.Terms[j], model.Inputs.Cases, column);
            double mean = column.Average();
            for (int i = 0; i < n; i++) xc[i, j] = column[i] - mean;
        }
        return xc;
    }

    // ── K = 1 against scalar BASS ─────────────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [MemberData(nameof(Centres))]
    public void AtOneOutput_TheCoreIsScalarBass(GPriorCentre centre)
    {
        const int n = 70;
        (List<BassPredictor> predictors, double[,] y) = Data(n, 1, Seed + 2);
        BassModel multi = Model(predictors, new MultiOutputResponse(y, ["y"]), centre);
        BassModel scalar = Model(predictors, new ScalarResponse(Column(y, 0)), centre);

        // The resolution: S_Y is s², bit for bit, and IW(2c, 2c·s²) is IG(c, c·s²).
        BassResolvedOutputs outputs = multi.Resolved.Outputs!;
        Assert.Equal(scalar.Resolved.ResponseMeanSquare, outputs.MeanSquares[0, 0]);
        Assert.Equal(2e-3, outputs.Prior.DegreesOfFreedom, 15);
        Close(InverseGamma.ShapeScale(1e-3, 1e-3 * scalar.Resolved.ResponseMeanSquare).LogPdf(0.37), outputs.Prior.LogPdf(new[,] { { 0.37 } }), 1e-13);

        BassDesign several = multi.NewDesign(), one = scalar.NewDesign();
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 3)[2]);
        double worstTarget = 0.0, worstDraw = 0.0;
        int compared = 0;
        for (int index = 0; index < 120; index++)
        {
            var terms = new List<BassTerm>();
            int m = rng.NextInt(6);
            while (terms.Count < m)
            {
                BassTerm t = RandomTerm(multi.Basis, rng);
                if (!terms.Contains(t)) terms.Add(t);
            }
            BassRefusal? refused = several.Load(terms);
            Assert.Equal(refused, one.Load(terms));
            if (refused is not null) continue;
            compared++;
            OutputForms forms = several.OutputForms!;
            Assert.Equal(one.ResidualSquares, forms.Residual[0, 0]);   // the same arithmetic on a load
            Assert.Equal(one.FittedSquares, forms.Fitted[0, 0]);

            double logTau = Math.Log(Math.Pow(10.0, -3.0 + 5.0 * rng.NextDouble())), b = rng.NextDouble() < 0.5 ? 1.0 : 0.05 + 0.95 * rng.NextDouble();
            double expected = scalar.Target.Evaluate(BassStructureFit.Of(one, 0.0), logTau).LogLikelihood(b);
            double actual = multi.Target.Evaluate(BassStructureFit.Of(several, 0.0), logTau).LogLikelihood(b);
            worstTarget = Math.Max(worstTarget, Math.Abs(actual - expected) / Math.Max(1.0, Math.Abs(expected)));

            // The draw: σ² and the centred coefficients (α, β_c) on the same stream.
            int seed = SeedTree.Derive(Seed + index, 1)[0];
            var centred = new double[one.Count + 1];
            double sigma2 = scalar.Target.CoefficientBlock(one)(Math.Exp(logTau), 1.0).Sample(scalar.Target.Sigma2, new Xoshiro256PlusPlus(seed), centred);
            (double[,] sigma, double[,] coefficients) = BassOutputDraws.Sample(multi.Target, several, BassStructureFit.Of(several, 0.0), Math.Exp(logTau), new Xoshiro256PlusPlus(seed));
            worstDraw = Math.Max(worstDraw, Math.Abs(sigma[0, 0] - sigma2) / sigma2);
            for (int j = 0; j <= one.Count; j++)
                worstDraw = Math.Max(worstDraw, Math.Abs(coefficients[j, 0] - centred[j]) / Math.Max(1.0, Math.Abs(centred[j])));
        }
        output.WriteLine($"{centre}: {compared} structures; the target's worst relative gap {worstTarget:E2}, the draw's {worstDraw:E2}");
        Assert.True(compared > 30, $"only {compared} structures loaded");
        Assert.True(worstTarget < 1e-12, $"the target's worst relative gap {worstTarget:E2}");
        Assert.True(worstDraw < 1e-12, $"the draw's worst relative gap {worstDraw:E2}");
    }

    // ── Equivariance ──────────────────────────────────────────────────────────────────────────────────────────────────────

    [Theory]
    [MemberData(nameof(Centres))]
    public void Marginal_UnderDecisionTwo_ShiftsEveryStructureAlike_WhenTheOutputsAreMixed(GPriorCentre centre)
    {
        // Y → YL takes S_Y to LᵀS_YL and the levels with it, so ln p(YL | s, τ) = ln p(Y | s, τ) − (N − f) ln|det L| at every structure and τ:
        // the structure posterior and τ's are unchanged. L rescales by 10⁻³ to 10³, mixes, permutes and flips a sign.
        const int n = 60, k = 3;
        (List<BassPredictor> predictors, double[,] y) = Data(n, k, Seed + 3);
        double[,] l = { { 0.0, 1e3, 0.2 }, { -0.7, 0.4, 1.1 }, { 1e-3, 0.0, 0.9 } };
        var yl = new double[n, k];
        for (int i = 0; i < n; i++)
            for (int d = 0; d < k; d++)
                for (int a = 0; a < k; a++) yl[i, d] += y[i, a] * l[a, d];
        double logDet = Math.Log(Math.Abs(l[0, 0] * (l[1, 1] * l[2, 2] - l[1, 2] * l[2, 1]) - l[0, 1] * (l[1, 0] * l[2, 2] - l[1, 2] * l[2, 0])
                                          + l[0, 2] * (l[1, 0] * l[2, 1] - l[1, 1] * l[2, 0])));
        BassModel plain = Model(predictors, new MultiOutputResponse(y, Names(k)), centre);
        BassModel mixed = Model(predictors, new MultiOutputResponse(yl, Names(k)), centre);
        BassDesign a0 = plain.NewDesign(), a1 = mixed.NewDesign();
        double expected = -(n - (centre == GPriorCentre.Flat ? 1.0 : 0.0)) * logDet, worst = 0.0;
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(Seed, 4)[3]);
        for (int index = 0; index < 40; index++)
        {
            var terms = new List<BassTerm>();
            int m = rng.NextInt(6);
            while (terms.Count < m)
            {
                BassTerm t = RandomTerm(plain.Basis, rng);
                if (!terms.Contains(t)) terms.Add(t);
            }
            if (a0.Load(terms) is not null) continue;
            Assert.Null(a1.Load(terms));
            double logTau = Math.Log(Math.Pow(10.0, -3.0 + 5.0 * rng.NextDouble()));
            double shift = mixed.Target.Evaluate(BassStructureFit.Of(a1, 0.0), logTau).LogMarginal - plain.Target.Evaluate(BassStructureFit.Of(a0, 0.0), logTau).LogMarginal;
            worst = Math.Max(worst, Math.Abs(shift - expected) / Math.Abs(expected));
        }
        output.WriteLine($"{centre}: shift −(N − f) ln|det L| = {expected:R}, worst relative gap {worst:E2}");
        Assert.True(worst < 1e-10, $"worst relative gap {worst:E2}");
    }

    // ── Refusals ──────────────────────────────────────────────────────────────────────────────────────────────────────────

    [Fact]
    public void Refuses_WhatLaterSlicesBring_AndOutputsThatCannotBeFitted()
    {
        (List<BassPredictor> predictors, double[,] y) = Data(40, 3, Seed + 4);
        var data = new BassData(predictors, new MultiOutputResponse(y, Names(3)));
        void Refused(BassConfig config) => Assert.Throws<ArgumentException>(() => BassModel.Create(config, data));
        Refused(new BassConfig { Noise = BassNoise.StudentT(4.0) });
        Refused(new BassConfig { Sampling = new BassSampling { Tempering = new TemperingSpec() } });
        Refused(new BassConfig { Sampling = new BassSampling { Smc = new SmcSpec() } });
        Refused(new BassConfig { Coefficients = new BassCoefficientPrior(new RidgeForm(), TauPrior.UnitInformationGamma()) });
        Refused(new BassConfig { Sigma2 = Sigma2Prior.Jeffreys });
        Refused(new BassConfig { Sigma2 = Sigma2Prior.InverseGamma(1.0, 1.0) });
        Refused(new BassConfig { Sigma2 = Sigma2Prior.HalfCauchy() });

        // A constant output, and collinear outputs: S_Y singular.
        var constant = (double[,])y.Clone();
        for (int i = 0; i < 40; i++) constant[i, 1] = 2.5;
        Assert.Throws<ArgumentException>(() => BassModel.Create(new BassConfig(), new BassData(predictors, new MultiOutputResponse(constant, Names(3)))));
        var collinear = (double[,])y.Clone();
        for (int i = 0; i < 40; i++) collinear[i, 2] = 2.0 * y[i, 0] - 0.5 * y[i, 1] + 7.0;
        ArgumentException e = Assert.Throws<ArgumentException>(() => BassModel.Create(new BassConfig(), new BassData(predictors, new MultiOutputResponse(collinear, Names(3)))));
        Assert.Contains("collinear", e.Message);

        // Missing outputs, and the names.
        var missing = (double[,])y.Clone();
        missing[3, 1] = double.NaN;
        Assert.Throws<ArgumentException>(() => new MultiOutputResponse(missing, Names(3)));
        Assert.Throws<ArgumentException>(() => new MultiOutputResponse(y, ["a", "b"]));
        Assert.Throws<ArgumentException>(() => new MultiOutputResponse(y, ["a", "b", "a"]));
        Assert.Throws<ArgumentException>(() => new MultiOutputResponse(y, ["a", "", "c"]));

        // A working response that moves, and the scalar marginal.
        BassModel model = BassModel.Create(new BassConfig(), data);
        Assert.Throws<NotSupportedException>(() => new WorkingResponse(data.Response, weights: Enumerable.Repeat(1.5, 40).ToArray()));
        Assert.Throws<InvalidOperationException>(() => model.Working.ResponseMean);
        Assert.Throws<InvalidOperationException>(() => model.Target.Marginal(BassStructureFit.Of(model.NewDesign(), 0.0), 1.0));
    }

    private static void Close(double expected, double actual, double rel)
        => Assert.True(Math.Abs(actual - expected) <= rel * Math.Max(1.0, Math.Abs(expected)), $"expected {expected:R}, got {actual:R}");
}
