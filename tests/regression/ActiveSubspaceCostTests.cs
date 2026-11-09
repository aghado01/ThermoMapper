using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's cost row (the active-subspace brief, slice 1), run with the qualification switch (<c>BASS_QUALIFICATION</c> = 1) on the clock
/// <c>BASS_COST_CLOCK</c> names, on the sensitivity fixture (two continuous inputs and a categorical one):
/// <list type="bullet">
/// <item>per hinge factor pair on x0, its derivative moments (<see cref="DesignFactor.SlopePair"/>) beside its value moments
/// (<see cref="DesignFactor.Pair"/>), under the uniform law, a normal–Student mixture and the empirical marginal of 60 atoms;</item>
/// <item>per term pair, <see cref="BassMoments.GradientPair"/> from a cold service;</item>
/// <item>per draw, C on the exact route (<see cref="BassMoments.GradientOuterProduct"/>) cold and with the caches warm, on the degree-one fit's
/// 60 draws;</item>
/// <item>the readout (<see cref="BassFit.ActiveSubspace"/>, slice 2) per draw: the exact route at 3,000 draws (the fit's 60, cycled),
/// with its consistency checks and without them as the summary alone; the sampled route at 256 draws and 4,096 points.</item>
/// <item>slice 3, per draw or pair: a direct functional fit's C over its 7 locations under the identity and a full W, a basis fit's over
/// its joint draws, two scalar fits' co-active pairs on both routes, and two outputs of the functional fit.</item>
/// <item>slice 4: a scalar fit's readout under a joint sample of 4,096 and 16,384 points, box and whitened; the basis fit's sampled
/// route per joint draw; and whitening itself, <see cref="Maths.Regression.Subspaces.Whitening"/>'s moments and root, at 16,384 points
/// of 2 and 10 inputs.</item>
/// <item>slice 5, per draw at 3,000 draws: the exact readout with and without the matrix check on a scalar, a functional and a basis
/// fit; with per-draw ranks; a functional fit's output side on both routes; and basis adequacy per draw of the direct fit.</item>
/// </list>
/// Indicative: the best of five on a shared machine.
/// </summary>
public sealed class ActiveSubspaceCostTests(ITestOutputHelper output)
{
    private const int Repeats = 5;

    [Fact]
    public void ReadoutCostRow()
    {
        if (!BassQualification.Enabled) return;
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = SensitivityFixture.FitScalar(1);
        int pooled = fit.Draws.Sum(c => c.Count);
        int[] many = [.. Enumerable.Range(0, 3000).Select(d => d % pooled)];
        var spec = new ActiveSubspaceSpec { Dimension = 2, Draws = many };
        double readout = Best(() => fit.ActiveSubspace(data, spec), many.Length);
        ActiveSubspaceAnalysis a = fit.ActiveSubspace(data, spec);
        double[][,] c = [.. a.Posterior.Draws.Select(d => d.C)];
        double[] w = [.. a.Posterior.Weights];
        double summary = Best(() => Maths.Regression.Subspaces.SubspacePosterior.Of(c, w, 2), many.Length);
        int[] some = [.. Enumerable.Range(0, 256).Select(d => d % pooled)];
        double sampled = Best(() => fit.ActiveSubspace(data, spec with { Draws = some, Route = ActiveSubspaceRoute.Sampled }), some.Length);
        output.WriteLine($"clock: {(BassCost.Cycles ? "thread cycles" : "wall")}; per draw: exact readout at {many.Length} draws {readout:F1} µs " +
                         $"(the summary alone {summary:F1} µs), sampled readout at {some.Length} draws and 4,096 points {sampled:F0} µs");
    }

    [Fact]
    public void FunctionalAndCoActiveCostRow()
    {
        if (!BassQualification.Enabled) return;
        BassData functional = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int m = functional.Response.Locations, pooled = fit.Draws.Sum(c => c.Count);
        int[] many = [.. Enumerable.Range(0, 3000).Select(d => d % pooled)];
        var spec = new ActiveSubspaceSpec { Draws = many };
        double identity = Best(() => fit.ActiveSubspace(functional, spec), many.Length);
        double full = Best(() => fit.ActiveSubspace(functional, spec with { OutputMetric = new OutputMetric(FunctionalActiveSubspaceTests.RandomMetric(m, 3)) }), many.Length);
        double[] a = [.. Enumerable.Range(0, m).Select(l => 1.0 + 0.1 * l)], b = [.. Enumerable.Range(0, m).Select(l => l % 2 == 0 ? 1.0 : -0.5)];
        double outputs = Best(() => fit.CoActiveSubspace(functional, new CoActiveSpec { Draws = many, Output = a, OtherOutput = b }), many.Length);

        BassBasisFit basis = BassOutputComponentsTests.FitOf();
        int[] joint = [.. Enumerable.Range(0, 3000).Select(d => d % basis.JointDraws.Count)];
        double basisC = Best(() => basis.ActiveSubspace(new ActiveSubspaceSpec { Draws = joint }), joint.Length);

        (BassData df, BassData dg) = CoActiveFixture.Data();
        (BassFit f, BassFit g) = CoActiveFixture.Fits();
        int[] pairs = [.. Enumerable.Range(0, 3000).Select(d => d % f.Draws.Sum(c => c.Count))];
        double exact = Best(() => f.CoActiveSubspace(df, g, dg, new CoActiveSpec { Draws = pairs, OtherDraws = pairs }), pairs.Length);
        int[] some = pairs[..256];
        double sampled = Best(() => f.CoActiveSubspace(df, g, dg, new CoActiveSpec { Draws = some, OtherDraws = some, Route = ActiveSubspaceRoute.Sampled }), some.Length);
        output.WriteLine($"clock: {(BassCost.Cycles ? "thread cycles" : "wall")}; per draw at 3,000 draws: functional C over {m} locations {identity:F1} µs " +
                         $"(identity W), {full:F1} µs (full W); two outputs {outputs:F1} µs; basis C over {basis.Decomposition.Locations} locations and " +
                         $"{basis.Components.Count} components {basisC:F1} µs; per pair: two scalar fits exact {exact:F1} µs, sampled at 4,096 points {sampled:F0} µs");
    }

    [Fact]
    public void JointSampleCostRow()
    {
        if (!BassQualification.Enabled) return;
        (BassData data, _) = CoActiveFixture.Data();
        (BassFit fit, _) = CoActiveFixture.Fits();
        int pooled = fit.Draws.Sum(c => c.Count);
        int[] some = [.. Enumerable.Range(0, 256).Select(d => d % pooled)];
        var rng = new Maths.Rng.Xoshiro256PlusPlus(5);
        var rows = new List<string>();
        foreach (int n in new[] { 4096, 16384 })
        {
            double[][] x = JointSampleActiveSubspaceTests.Triangle(rng, n);
            var sample = new JointSample(new PredictionCases([new ContinuousColumn("x0", x[0]), new ContinuousColumn("x1", x[1])]));
            var spec = new ActiveSubspaceSpec { Draws = some, Sample = sample };
            double box = Best(() => fit.ActiveSubspace(data, spec), some.Length);
            double white = Best(() => fit.ActiveSubspace(data, spec with { Metric = InputMetric.Whitened }), some.Length);
            rows.Add($"{n} points: {box:F0} µs box, {white:F0} µs whitened");
        }
        BassBasisFit basis = BassOutputComponentsTests.FitOf();
        int[] joint = [.. Enumerable.Range(0, 128).Select(d => d % basis.JointDraws.Count)];
        double basisSampled = Best(() => basis.ActiveSubspace(new ActiveSubspaceSpec { Draws = joint, Route = ActiveSubspaceRoute.Sampled }), joint.Length);
        var whitening = new List<string>();
        foreach (int p in new[] { 2, 10 })
        {
            double[][] cols = [.. Enumerable.Range(0, p).Select(_ => Enumerable.Range(0, 16384).Select(_ => rng.NextDouble()).ToArray())];
            double cost = Best(() => Maths.Regression.Subspaces.Whitening.SymmetricRoot(Maths.Regression.Subspaces.Whitening.Moments(cols).Covariance), 1);
            whitening.Add($"p {p}: {cost:F0} µs");
        }
        output.WriteLine($"clock: {(BassCost.Cycles ? "thread cycles" : "wall")}; per draw under a joint sample, {string.Join("; ", rows)}; " +
                         $"the basis fit's sampled route per joint draw at 4,096 points {basisSampled:F0} µs; whitening at 16,384 points, {string.Join(", ", whitening)}");
    }

    [Fact]
    public void SettledItemsCostRow()
    {
        if (!BassQualification.Enabled) return;
        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        int[] many = [.. Enumerable.Range(0, 3000).Select(d => d % f.Draws.Sum(c => c.Count))];
        var spec = new ActiveSubspaceSpec { Draws = many };
        double scalarOff = Best(() => f.ActiveSubspace(df, spec), many.Length);
        double plain = Best(() => f.ActiveSubspace(df, spec with { MatrixPoincare = true }), many.Length);
        double ranked = Best(() => f.ActiveSubspace(df, spec with { Rank = Maths.Regression.Subspaces.RankRule.TraceShare(0.95) }), many.Length) - scalarOff;
        BassData functional = SensitivityFixture.Functional();
        BassFit fit = SensitivityFixture.FitFunctional();
        int[] fdraws = [.. Enumerable.Range(0, 3000).Select(d => d % fit.Draws.Sum(c => c.Count))];
        var fspec = new ActiveSubspaceSpec { Draws = fdraws };
        double without = Best(() => fit.ActiveSubspace(functional, fspec), fdraws.Length);
        double functionalOn = Best(() => fit.ActiveSubspace(functional, fspec with { MatrixPoincare = true }), fdraws.Length);
        BassBasisFit basisFit = BassOutputComponentsTests.FitOf();
        int[] joint = [.. Enumerable.Range(0, 3000).Select(d => d % basisFit.JointDraws.Count)];
        double basisOff = Best(() => basisFit.ActiveSubspace(new ActiveSubspaceSpec { Draws = joint }), joint.Length);
        double basisOn = Best(() => basisFit.ActiveSubspace(new ActiveSubspaceSpec { Draws = joint, MatrixPoincare = true }), joint.Length);
        double exactG = Best(() => fit.ActiveSubspace(functional, fspec with { OutputPatterns = 2 }), fdraws.Length);
        int[] some = fdraws[..256];
        double sampledG = Best(() => fit.ActiveSubspace(functional, fspec with { Draws = some, OutputPatterns = 2, Route = ActiveSubspaceRoute.Sampled }), some.Length)
                          - Best(() => fit.ActiveSubspace(functional, fspec with { Draws = some, Route = ActiveSubspaceRoute.Sampled }), some.Length);
        (BassData ad, BassFit direct, BassPredictor[] predictors, double[,] y, _, _) = RankAndOutputSideTests.AdequacyFixture.Value;
        BassBasisFit basis = BassComponents.Fit(BassDecomposition.FromPca(y, new BassPcaSpec { Rank = 1 }), predictors, BassOutputComponentsTests.Config());
        int[] adraws = [.. Enumerable.Range(0, 600).Select(d => d % direct.Draws.Sum(c => c.Count))];
        double adequacy = Best(() => basis.Adequacy(direct, ad, new ActiveSubspaceSpec { Draws = adraws }), adraws.Length);
        output.WriteLine($"clock: {(BassCost.Cycles ? "thread cycles" : "wall")}; per draw at 3,000 draws, the exact readout with the matrix check and " +
                         $"without: scalar {plain:F1} and {scalarOff:F1} µs, functional {functionalOn:F1} and {without:F1} µs, basis {basisOn:F1} and {basisOff:F1} µs; " +
                         $"per-draw ranks add {ranked:F1} µs to the scalar; a functional fit's exact readout {without:F1} µs, with two output patterns " +
                         $"{exactG:F1} µs; the sampled output side's own cost at 4,096 points {sampledG:F0} µs; basis adequacy per direct draw (12 locations) {adequacy:F1} µs");
    }

    private static double Best(Action call, int calls)
    {
        call();
        double best = double.PositiveInfinity;
        for (int r = 0; r < Repeats; r++)
        {
            long t0 = BassCost.Tick();
            call();
            best = Math.Min(best, (BassCost.Tick() - t0) / BassCost.TicksPerSecond * 1e6 / calls);
        }
        return best;   // µs per call
    }

    [Fact]
    public void CostRow()
    {
        if (!BassQualification.Enabled) return;
        BassData data = SensitivityFixture.Scalar();
        BassFit fit = SensitivityFixture.FitScalar(1);
        BassInputs inputs = BassInputs.Prepare(data, fit.Transform);
        var laws = new (string Name, Dictionary<string, InputMeasure> Measures)[]
        {
            ("uniform", []),
            ("mixture", new() { ["x0"] = new MixtureMeasure([new NormalComponent(0.4, 0.3, 0.15), new StudentComponent(0.6, 2.0, 0.7, 0.2)]) }),
            ("empirical", new() { ["x0"] = new EmpiricalMeasure() }),
        };
        ReadOnlySpan<double> knots = inputs.ContinuousKnots(0);
        var factors = new List<HingeDesignFactor>();
        for (int k = 0; k < knots.Length; k++)
            foreach (int s in new[] { 1, -1 }) factors.Add(new HingeDesignFactor(0, s, k, knots[k], 1));
        int pairs = factors.Count * factors.Count;
        BassDraw[] draws = [.. fit.Draws.SelectMany(c => c)];
        BassTerm[][] drawTerms = [.. draws.Select(d => d.Structure.TermIds.Select(id => fit.Terms[id]).ToArray())];
        BassTerm[] termSet = [.. draws.SelectMany(d => d.Structure.TermIds).Distinct().Select(id => fit.Terms[id])];
        output.WriteLine($"clock: {(BassCost.Cycles ? "thread cycles" : "wall")}; {factors.Count} hinge factors on x0, {termSet.Length} distinct terms, {draws.Length} draws " +
                         $"of up to {drawTerms.Max(t => t.Length)} terms");
        foreach ((string name, Dictionary<string, InputMeasure> measures) in laws)
        {
            var measure = ProductMeasure.Resolve(fit.Transform, inputs, measures);
            ContinuousMeasure law = measure.Continuous(0)!;
            double value = Best(() =>
            {
                foreach (HingeDesignFactor a in factors)
                    foreach (HingeDesignFactor b in factors) a.Pair(b, 0.5, 0.5, measure, out _);
            }, pairs);
            double slope = Best(() =>
            {
                foreach (HingeDesignFactor a in factors)
                    foreach (HingeDesignFactor b in factors) a.SlopePair(b, law, out _);
            }, pairs);
            int termPairs = termSet.Length * termSet.Length;
            double gradientPair = Best(() =>
            {
                BassMoments cold = BassMoments.On(inputs, 1, measure);
                foreach (BassTerm a in termSet)
                    foreach (BassTerm b in termSet) cold.GradientPair(a, b);
            }, termPairs);
            double coldDraw = Best(() =>
            {
                BassMoments cold = BassMoments.On(inputs, 1, measure);
                for (int d = 0; d < draws.Length; d++) cold.GradientOuterProduct(drawTerms[d], draws[d].Coefficients.AsSpan(1));
            }, draws.Length);
            BassMoments warm = BassMoments.On(inputs, 1, measure);
            double warmDraw = Best(() =>
            {
                for (int d = 0; d < draws.Length; d++) warm.GradientOuterProduct(drawTerms[d], draws[d].Coefficients.AsSpan(1));
            }, draws.Length);
            output.WriteLine($"{name}: factor pair, value moments {value:F2} µs, derivative moments {slope:F2} µs; term pair cold {gradientPair:F2} µs; " +
                             $"C per draw cold {coldDraw:F1} µs, warm {warmDraw:F2} µs");
        }
    }
}
