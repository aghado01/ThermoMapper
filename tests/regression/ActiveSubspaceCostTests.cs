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
/// 60 draws.</item>
/// </list>
/// Indicative: the best of five on a shared machine.
/// </summary>
public sealed class ActiveSubspaceCostTests(ITestOutputHelper output)
{
    private const int Repeats = 5;

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
