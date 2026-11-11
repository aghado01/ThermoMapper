using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 6, item 7: a Legendre factor's moments under a uniform law without the monomial expansion, so high degrees keep full
/// precision (<see cref="LegendreDesignFactor"/>):
/// <list type="bullet">
/// <item>on [0, 1] to degree 40, the closed forms: E[ψ_a] = δ_a0, E[ψ_aψ_b] = δ_ab, the centred pair, and the slope moments against
/// the test's own closed forms (<see cref="LegendreExpansion"/>) and against Gauss–Legendre quadrature of the recurrence;</item>
/// <item>on other intervals, the exact rule against a rule of twice the nodes, and at low degree against the expansion, which a
/// truncated wider uniform still takes;</item>
/// <item>C of hand-built draws of degree up to 25 against tensor quadrature of the polynomials' own derivatives, within 10⁻¹⁴.</item>
/// </list>
/// </summary>
public sealed class LegendreUniformMomentTests(ITestOutputHelper output)
{
    private const int Top = 40;

    // One continuous predictor x on [0, 1] (its box), and its law resolved from the declared measure.
    private static (BassInputs Inputs, ProductMeasure Measure) OneInput(InputMeasure? law = null, int inputs = 1)
    {
        double[] x = [.. Enumerable.Range(0, 11).Select(i => i / 10.0)];
        var predictors = new List<BassPredictor> { new ContinuousPredictor("x0", x) };
        if (inputs == 2) predictors.Add(new ContinuousPredictor("x1", [.. x.Select(v => 1.0 - v)]));
        var data = new BassData(predictors, new ScalarResponse(x));
        var transform = BassInputTransform.Fit(data);
        BassInputs prepared = BassInputs.Prepare(data, transform);
        var measures = law is null ? null : new Dictionary<string, InputMeasure> { ["x0"] = law };
        return (prepared, ProductMeasure.Resolve(transform, prepared, measures));
    }

    // E[g] under the uniform law on [lo, hi] by an n-point Gauss–Legendre rule.
    private static double Quadrature(double lo, double hi, int n, Func<double, double> g)
    {
        GaussLegendre rule = GaussLegendre.Of(n);
        double half = 0.5 * (hi - lo), centre = 0.5 * (hi + lo), s = 0.0;
        for (int k = 0; k < rule.Count; k++) s += rule.Weights[k] * g(centre + half * rule.Nodes[k]);
        return 0.5 * s;
    }

    /// <summary>
    /// On [0, 1] to degree 40: the means and pairs exactly as orthonormality gives them, the centred pair with the means given, and the
    /// slope moments within 2·10⁻¹⁶ of the test's closed forms and 10⁻¹³ of 64-node quadrature of the recurrence, both relative to the
    /// slope Gram's scale.
    /// </summary>
    [Fact]
    public void OnTheBox_ClosedForms_ToDegreeForty()
    {
        (_, ProductMeasure measure) = OneInput();
        ContinuousMeasure law = measure.Continuous(0)!;
        Assert.Equal(0.0, law.Lower);
        Assert.Equal(1.0, law.Upper);
        double closed = 0.0, quadrature = 0.0;
        for (int a = 0; a <= Top; a++)
        {
            var fa = new LegendreDesignFactor(0, a);
            Assert.Equal(a == 0 ? 1.0 : 0.0, fa.Mean(measure, out double error));
            Assert.Equal(0.0, error);
            double scaleA = Math.Max(1.0, LegendreExpansion.SlopeGram(a, a));
            double slopeMean = fa.SlopeMean(law, out _);
            closed = Math.Max(closed, Math.Abs(slopeMean - (a % 2 == 1 ? 2.0 * Math.Sqrt(2.0 * a + 1.0) : 0.0)) / scaleA);
            quadrature = Math.Max(quadrature, Math.Abs(slopeMean - Quadrature(0.0, 1.0, 64, x => ShiftedLegendre.Derivative(a, x))) / scaleA);
            for (int b = 0; b <= Top; b++)
            {
                var fb = new LegendreDesignFactor(0, b);
                (double cross, double covariance) = fa.Pair(fb, 0.3, -0.2, measure, out _);
                Assert.Equal(a == b ? 1.0 : 0.0, cross);
                double centred = (a == b ? 1.0 : 0.0) - 0.3 * (b == 0 ? 1.0 : 0.0) + 0.2 * (a == 0 ? 1.0 : 0.0) - 0.06;
                Assert.Equal(centred, covariance, 15);
                (double sc, double ss) = fa.SlopePair(fb, law, out _);
                double scale = Math.Max(1.0, LegendreExpansion.SlopeGram(Math.Max(a, b), Math.Max(a, b)));
                closed = Math.Max(closed, Math.Max(Math.Abs(sc - LegendreExpansion.SlopeProduct(a, b)), Math.Abs(ss - LegendreExpansion.SlopeGram(a, b))) / scale);
                int aa = a, bb = b;
                double qc = Quadrature(0.0, 1.0, 64, x => ShiftedLegendre.Derivative(aa, x) * ShiftedLegendre.Value(bb, x));
                double qs = Quadrature(0.0, 1.0, 64, x => ShiftedLegendre.Derivative(aa, x) * ShiftedLegendre.Derivative(bb, x));
                quadrature = Math.Max(quadrature, Math.Max(Math.Abs(sc - qc), Math.Abs(ss - qs)) / scale);
            }
        }
        output.WriteLine($"to degree {Top} on [0, 1]: slope moments against the closed forms {closed:E2}, against quadrature of the recurrence {quadrature:E2}");
        Assert.True(closed <= 2e-16, $"{closed:E2}");
        Assert.True(quadrature <= 1e-13, $"{quadrature:E2}");
    }

    /// <summary>
    /// On [0.1, 0.85], [−0.2, 1.3] and [0.4, 0.45], to degree 30: every moment against a rule of twice the nodes plus eight, within
    /// 10⁻¹² of the integrand's scale; and to degree six against the monomial expansion, which a uniform truncated from a wider one
    /// still takes, within 10⁻¹² of that scale, the gap at degree 30 printed beside it.
    /// </summary>
    [Fact]
    public void OnOtherIntervals_TheExactRule()
    {
        (double Lo, double Hi)[] intervals = [(0.1, 0.85), (-0.2, 1.3), (0.4, 0.45)];
        foreach ((double lo, double hi) in intervals)
        {
            (_, ProductMeasure measure) = OneInput(new UniformMeasure(lo, hi));
            ContinuousMeasure law = measure.Continuous(0)!;
            var expansion = new ContinuousMeasure(Truncated.Of(Uniform.On(lo - 1.0, hi + 1.0), lo, hi), lo, hi);
            double rule = 0.0, low = 0.0, high = 0.0;
            for (int a = 0; a <= 30; a++)
                for (int b = 0; b <= 30; b++)
                {
                    var fa = new LegendreDesignFactor(0, a);
                    var fb = new LegendreDesignFactor(0, b);
                    int aa = a, bb = b, n = (a + b + 2) / 2 * 2 + 8;
                    double size = Quadrature(lo, hi, n, x => Math.Abs(ShiftedLegendre.Derivative(aa, x) * ShiftedLegendre.Derivative(bb, x))
                                                            + Math.Abs(ShiftedLegendre.Value(aa, x) * ShiftedLegendre.Value(bb, x))) + 1.0;
                    (double cross, _) = fa.Pair(fb, 0.0, 0.0, measure, out _);
                    (double sc, double ss) = fa.SlopePair(fb, law, out _);
                    double qv = Quadrature(lo, hi, n, x => ShiftedLegendre.Value(aa, x) * ShiftedLegendre.Value(bb, x));
                    double qc = Quadrature(lo, hi, n, x => ShiftedLegendre.Derivative(aa, x) * ShiftedLegendre.Value(bb, x));
                    double qs = Quadrature(lo, hi, n, x => ShiftedLegendre.Derivative(aa, x) * ShiftedLegendre.Derivative(bb, x));
                    rule = Math.Max(rule, Math.Max(Math.Abs(cross - qv), Math.Max(Math.Abs(sc - qc), Math.Abs(ss - qs))) / size);
                    if (Math.Max(a, b) > 6 && Math.Max(a, b) < 30) continue;
                    (double ec, double es) = fa.SlopePair(fb, expansion, out _);
                    double gap = Math.Max(Math.Abs(sc - ec), Math.Abs(ss - es)) / size;
                    if (Math.Max(a, b) <= 6) low = Math.Max(low, gap);
                    else high = Math.Max(high, gap);
                }
            output.WriteLine($"[{lo}, {hi}] to degree 30: against the wider rule {rule:E2}; the expansion to degree six {low:E2}, at degree 30 {high:E2}");
            Assert.True(rule <= 1e-12, $"[{lo}, {hi}]: {rule:E2}");
            Assert.True(low <= 1e-12, $"[{lo}, {hi}]: {low:E2}");
        }
    }

    /// <summary>
    /// C = Σ c_mc_n E[∇B_m∇B_nᵀ] for three Legendre terms of degree up to 25 on two inputs, under the uniform law on the box, against a
    /// 30 × 30 tensor Gauss–Legendre rule of the polynomials' own derivatives by their recurrence (exact for the integrand's degree):
    /// within 10⁻¹⁴ of C's largest entry.
    /// </summary>
    [Fact]
    public void HighDegreeDraws_C_IsTheTensorQuadrature()
    {
        (BassInputs inputs, ProductMeasure measure) = OneInput(inputs: 2);
        var moments = BassMoments.On(inputs, 1, measure);
        (int A, int B)[] degrees = [(25, 3), (18, 0), (7, 22)];
        double[] c = [0.7, -1.3, 0.45];
        BassTerm[] terms = [.. degrees.Select(d => new BassTerm(legendre: d.B == 0
            ? [new LegendreFactor(0, d.A)]
            : [new LegendreFactor(0, d.A), new LegendreFactor(1, d.B)]))];
        double[,] exact = moments.GradientOuterProduct(terms, c);
        GaussLegendre rule = GaussLegendre.Of(30);
        var reference = new double[2, 2];
        for (int i = 0; i < rule.Count; i++)
            for (int j = 0; j < rule.Count; j++)
            {
                double x = 0.5 + 0.5 * rule.Nodes[i], y = 0.5 + 0.5 * rule.Nodes[j], w = 0.25 * rule.Weights[i] * rule.Weights[j];
                double gx = 0.0, gy = 0.0;
                for (int m = 0; m < terms.Length; m++)
                {
                    (int a, int b) = degrees[m];
                    gx += c[m] * ShiftedLegendre.Derivative(a, x) * ShiftedLegendre.Value(b, y);
                    gy += c[m] * ShiftedLegendre.Value(a, x) * ShiftedLegendre.Derivative(b, y);
                }
                reference[0, 0] += w * gx * gx;
                reference[0, 1] += w * gx * gy;
                reference[1, 0] += w * gy * gx;
                reference[1, 1] += w * gy * gy;
            }
        double scale = 0.0, worst = 0.0;
        foreach (double v in reference) scale = Math.Max(scale, Math.Abs(v));
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++) worst = Math.Max(worst, Math.Abs(exact[i, j] - reference[i, j]) / scale);
        output.WriteLine($"degrees to 25: C against tensor quadrature {worst:E2} of its largest entry {scale:E3}; the service's error estimate {moments.LargestError}");
        Assert.True(worst <= 1e-14, $"{worst:E2}");
        Assert.Equal(0.0, moments.LargestError);
    }
}
