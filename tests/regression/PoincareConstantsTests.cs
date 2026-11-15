using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Regression.Subspaces;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 8a, item B: Poincaré constants for the remaining laws. Mixtures of several components and Student t components take the
/// finite-element constant (<see cref="IntervalPoincare"/>; its own gates are <c>IntervalPoincareTests</c>), which both checks read:
/// the per-input bound is that constant times C_vv, neither check finds a violation on fitted draws under such laws on either route,
/// and the matrix check is computed. Through the moment service, independent of the elements' quadrature, the matrix check's pencil over
/// a dense span of hinges reaches the constant from below and fires when it is scaled down. The uniform law keeps its exact constant,
/// the single truncated normal its bound (slice 6), and the empirical marginal has none, which the readout states.
/// </summary>
public sealed class PoincareConstantsTests(ITestOutputHelper output)
{
    // Laws in x's own units over its box [min, min + width]: a bimodal mixture, a Student t with three degrees of freedom, and a normal
    // mixed with a t.
    private static IEnumerable<(string Name, InputMeasure Law)> Laws(double min, double width) =>
    [
        ("bimodal", new MixtureMeasure([new NormalComponent(0.5, min + 0.25 * width, 0.06 * width), new NormalComponent(0.5, min + 0.75 * width, 0.06 * width)])),
        ("t3", new MixtureMeasure([new StudentComponent(1.0, 3.0, min + 0.4 * width, 0.15 * width)])),
        ("normal and t", new MixtureMeasure([new NormalComponent(0.7, min + 0.3 * width, 0.1 * width), new StudentComponent(0.3, 2.0, min + 0.8 * width, 0.05 * width)])),
    ];

    /// <summary>
    /// A scalar fit with x0 under each law and x1 under the bimodal one: both constants numerical, the per-input bound each constant times
    /// C_vv, each constant the finite-element one of the law in box units; no violation of either check on the exact or sampled route, and
    /// the matrix check computed.
    /// </summary>
    [Fact]
    public void NumericalConstants_AreWhatTheChecksRead()
    {
        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        var s0 = (RangeScale)f.Transform.Predictors[0];
        var s1 = (RangeScale)f.Transform.Predictors[1];
        (string _, InputMeasure bimodal) = Laws(s1.Min, s1.Max - s1.Min).First();
        double worstRead = 0.0, worstRatio = 0.0, worstMatrix = 0.0;
        foreach ((string name, InputMeasure law) in Laws(s0.Min, s0.Max - s0.Min))
        {
            var measures = new Dictionary<string, InputMeasure> { [s0.Name] = law, [s1.Name] = bimodal };
            foreach (ActiveSubspaceRoute route in new[] { ActiveSubspaceRoute.Exact, ActiveSubspaceRoute.Sampled })
            {
                ActiveSubspaceAnalysis a = f.ActiveSubspace(df, new ActiveSubspaceSpec { Measures = measures, MatrixPoincare = true, Route = route, Points = 4096 });
                PoincareCheck check = a.Poincare!;
                for (int i = 0; i < 2; i++)
                {
                    InputPoincareConstant c = check.Constants![i];
                    Assert.Equal(PoincareConstantKind.Numerical, c.Kind);
                    Assert.True(c.Elements >= IntervalPoincare.FirstGrid * 2 && c.RelativeError <= 1e-7, $"{name}, input {i}: {c.Elements} elements, {c.RelativeError:E2}");
                    ContinuousMeasure resolved = a.Measure!.Continuous(a.Measure.OfContinuous(i))!;
                    IntervalPoincare direct = IntervalPoincare.Of(((IContinuousDistribution)resolved.Measure).LogPdf, resolved.Lower, resolved.Upper);
                    Assert.Equal(direct.Constant, c.Value);
                    for (int d = 0; d < a.Posterior.Draws.Count; d++)
                    {
                        double cvv = a.Posterior.Draws[d].C[i, i];
                        if (cvv > 0.0) worstRead = Math.Max(worstRead, Math.Abs(check.Bound[d, i] / cvv - c.Value) / c.Value);
                        worstRatio = Math.Max(worstRatio, check.Ratio[d, i]);
                    }
                }
                Assert.Equal(0, check.Violations);
                Assert.Equal(MatrixPoincareStatus.Computed, a.MatrixPoincareStatus);
                Assert.Equal(0, a.MatrixPoincare!.Violations);
                worstMatrix = Math.Max(worstMatrix, a.MatrixPoincare.Ratios.Max());
                if (route == ActiveSubspaceRoute.Exact)
                    output.WriteLine($"{name}: constants {check.Constants![0].Value:E4} and {check.Constants[1].Value:E4} in box units, " +
                                     $"{check.Constants[0].Elements} and {check.Constants[1].Elements} elements");
            }
        }
        output.WriteLine($"worst |Bound/C_vv − constant|/constant {worstRead:E2}; largest per-input ratio {worstRatio:F4}; largest matrix ratio {worstMatrix:F4}");
        // Indicative cost: the exact readout with both inputs under the bimodal law, against both uniform.
        var both = new Dictionary<string, InputMeasure> { [s0.Name] = Laws(s0.Min, s0.Max - s0.Min).First().Law, [s1.Name] = bimodal };
        double Seconds(ActiveSubspaceSpec spec) { long t = Stopwatch.GetTimestamp(); f.ActiveSubspace(df, spec); return Stopwatch.GetElapsedTime(t).TotalSeconds; }
        Seconds(new ActiveSubspaceSpec());
        double uniformTime = Seconds(new ActiveSubspaceSpec { MatrixPoincare = true }), mixtureTime = Seconds(new ActiveSubspaceSpec { Measures = both, MatrixPoincare = true });
        output.WriteLine($"indicative: the exact readout with the matrix check, {uniformTime:F3} s under uniform laws, {mixtureTime:F3} s under two bimodal ones");
        Assert.True(worstRead <= 1e-15, $"{worstRead:E2}");
    }

    /// <summary>
    /// Each kind as the readout states it: a uniform law's exact (b − a)²/π², a single truncated normal's bound min(σ², (b − a)²/π²)
    /// (slice 6's readings unchanged), and the empirical marginal's none, with its ratios NaN and the matrix check not applicable.
    /// </summary>
    [Fact]
    public void EachKind_AsStated()
    {
        (BassData df, _) = CoActiveFixture.Data();
        (BassFit f, _) = CoActiveFixture.Fits();
        var s0 = (RangeScale)f.Transform.Predictors[0];
        double width = s0.Max - s0.Min;
        ActiveSubspaceAnalysis closed = f.ActiveSubspace(df, new ActiveSubspaceSpec
        {
            MatrixPoincare = true,
            Measures = new Dictionary<string, InputMeasure> { [s0.Name] = new MixtureMeasure([new NormalComponent(1.0, s0.Min + 0.45 * width, 0.08 * width)]) },
        });
        InputPoincareConstant normal = closed.Poincare!.Constants![0], uniform = closed.Poincare.Constants[1];
        Assert.Equal(PoincareConstantKind.Bound, normal.Kind);
        Assert.Equal(Math.Min(0.08 * 0.08, 1.0 / (Math.PI * Math.PI)), normal.Value, 15);
        Assert.Equal(PoincareConstantKind.Exact, uniform.Kind);
        Assert.Equal(1.0 / (Math.PI * Math.PI), uniform.Value, 15);
        Assert.Equal((0, 0.0), (normal.Elements, normal.RelativeError));
        Assert.Equal(MatrixPoincareStatus.Computed, closed.MatrixPoincareStatus);

        ActiveSubspaceAnalysis empirical = f.ActiveSubspace(df, new ActiveSubspaceSpec
        {
            MatrixPoincare = true, Measures = new Dictionary<string, InputMeasure> { [s0.Name] = new EmpiricalMeasure() },
        });
        InputPoincareConstant none = empirical.Poincare!.Constants![0];
        Assert.Equal(PoincareConstantKind.None, none.Kind);
        Assert.True(double.IsNaN(none.Value));
        for (int d = 0; d < empirical.Posterior.Draws.Count; d++) Assert.True(double.IsNaN(empirical.Poincare.Ratio[d, 0]));
        Assert.Equal(MatrixPoincareStatus.NotApplicable, empirical.MatrixPoincareStatus);
        Assert.Null(empirical.MatrixPoincare);
    }

    /// <summary>
    /// Through the moment service: on training values of x0 in [0, 1], every 0.04 and every 0.003 across [0.44, 0.56], the hinges
    /// (x0 − k)₊ at every knot below 1 span the piecewise linear functions on that grid. The fine run is where the bimodal law's first
    /// eigenfunction turns: its derivative goes as 1/ρ, which across the gap is about sech(69(x − ½)), so a uniform grid at 0.02 reaches
    /// only 0.935 of the constant. Under each law on x0, the largest Var(f)/(C·E[f′²]) over the span, the matrix check's pencil with
    /// D = diag(C, 0), is at most 1 within C's stated error and above 0.97: the constant is reached from below, not loose. With C scaled
    /// by 0.95 the ratio exceeds 1, so the check fires.
    /// </summary>
    [Fact]
    public void HingeSpan_ReachesTheConstant_AndAScaledOneFires()
    {
        double[] x0 = [.. Enumerable.Range(0, 26).Select(i => i / 25.0).Concat(Enumerable.Range(0, 41).Select(j => 0.44 + 0.003 * j))
                                    .Select(x => Math.Round(x, 9)).Distinct().Order()];
        double[] x1 = [.. Enumerable.Range(0, x0.Length).Select(i => (i * 7 % x0.Length) / (x0.Length - 1.0))];
        var data = new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse([.. x0.Select((v, i) => v + x1[i])]));
        var transform = BassInputTransform.Fit(data);
        BassInputs inputs = BassInputs.Prepare(data, transform);
        double[] knots = inputs.ContinuousKnots(0).ToArray();
        BassTerm[] terms = [.. Enumerable.Range(0, knots.Length).Where(k => knots[k] < 1.0).Select(k => new BassTerm([new HingeFactor(0, 1, k)]))];
        foreach ((string name, InputMeasure law) in Laws(0.0, 1.0))
        {
            var measure = ProductMeasure.Resolve(transform, inputs, new Dictionary<string, InputMeasure> { ["x0"] = law });
            ContinuousMeasure resolved = measure.Continuous(0)!;
            IntervalPoincare constant = IntervalPoincare.Of(((IContinuousDistribution)resolved.Measure).LogPdf, resolved.Lower, resolved.Upper);
            var moments = BassMoments.On(inputs, 1, measure);
            double[,] covariance = moments.CovarianceGram(terms);
            double Ratio(double c) => MatrixPoincare.PencilRatio(covariance, moments.GradientGram(terms, new double[,] { { c, 0.0 }, { 0.0, 0.0 } }));
            double reached = Ratio(constant.Constant), scaled = Ratio(0.95 * constant.Constant);
            output.WriteLine($"{name}: constant {constant.Constant:E4}; the span's largest ratio 1 − {1.0 - reached:E2}; with 0.95C {scaled:F4}");
            Assert.InRange(reached, 0.97, 1.0 + constant.RelativeError + 1e-9);
            Assert.True(scaled > 1.0, $"{name}: {scaled:F6}");
        }
    }
}
