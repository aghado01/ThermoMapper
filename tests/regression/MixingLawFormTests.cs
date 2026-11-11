using System;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The mixing laws' step for several outputs sharing one v (E34's decision 6; <see cref="IMixingLaw.StepConditionalForm"/>):
/// ∝ p(v)·v^(−bK/2)·exp(−b·E_r/(2v) − b·E_γ·v/2), E_r = rᵀΣ⁻¹r/c and E_γ = γᵀΣ⁻¹γ/c.
/// <list type="bullet">
/// <item><b>K = 1 is the scalar step</b> on the same stream, at E_r = r²/s² and E_γ = m²/s², for t, the Normal–Wald's GIG, a declared GIG
/// and the horseshoe's GBP.</item>
/// <item><b>K = 3 against quadrature:</b> E[ln v] and P(v &lt; 1) of the stated conditional on u = ln v, against exact draws (a chain for the
/// GBP, whose step is a slice step), with and without a drift, at b = 1 and 0.4.</item>
/// </list>
/// </summary>
public sealed class MixingLawFormTests(ITestOutputHelper output)
{
    private static IMixingLaw[] Laws() =>
    [
        new StudentTMixing(4.0),
        new GigMixing(new GigParameters(-0.5, 0.8 * 0.8, 1.0)),
        new GigMixing(new GigParameters(1.3, 0.7, 0.4)),
        new GbpMixing(new GbpParameters(1.0, 0.5, 0.5)),
    ];

    [Fact]
    public void KOne_IsTheScalarStep()
    {
        double worst = 0.0;
        foreach (IMixingLaw law in Laws())
            foreach ((double r, double s2, double m, double b) in new[] { (0.7, 0.5, 0.0, 1.0), (-1.3, 2.0, 0.4, 0.6), (0.05, 0.3, -0.8, 1.0) })
            {
                if (law is StudentTMixing && m != 0.0) continue;   // the skewed t is the GIG's, checked there
                var a = new Xoshiro256PlusPlus(7);
                var c = new Xoshiro256PlusPlus(7);
                double va = 1.3, vc = 1.3;
                for (int t = 0; t < 200; t++)
                {
                    va = law.StepConditional(va, r, s2, m, b, a);
                    vc = law.StepConditionalForm(vc, r * r / s2, m * m / s2, 1, b, c);
                    worst = Math.Max(worst, Math.Abs(va - vc) / va);
                }
            }
        output.WriteLine($"K = 1 against the scalar step: {worst:G3}");
        Assert.True(worst < 1e-10, $"{worst:G3}");
    }

    [Fact]
    public void KThree_MeetsItsConditional_ByQuadrature()
    {
        double worst = 0.0;
        foreach (IMixingLaw law in Laws())
            foreach ((double er, double eg, double b) in new[] { (2.4, 0.0, 1.0), (0.9, 0.6, 1.0), (4.0, 1.5, 0.4) })
            {
                if (law is StudentTMixing && eg != 0.0) continue;
                const int k = 3, draws = 60_000, batches = 40;
                (double meanLog, double below) = Quadrature(law, er, eg, k, b);
                var rng = new Xoshiro256PlusPlus(11);
                var logs = new double[draws];
                double v = 1.0;
                int count = 0;
                for (int t = 0; t < draws; t++)
                {
                    v = law.StepConditionalForm(v, er, eg, k, b, rng);
                    logs[t] = Math.Log(v);
                    if (v < 1.0) count++;
                }
                double mean = logs.Average(), se = BatchError(logs, batches);
                double z = (mean - meanLog) / se;
                double p = (double)count / draws, zp = (p - below) / Math.Sqrt(below * (1 - below) / draws * (law is GbpMixing ? 10.0 : 1.0));
                output.WriteLine($"{law.GetType().Name,-15} E_r {er} E_γ {eg} b {b}: E[ln v] {mean:F4} against {meanLog:F4} (z {z:F2}); P(v<1) {p:F4} against {below:F4} (z {zp:F2})");
                worst = Math.Max(worst, Math.Max(Math.Abs(z), Math.Abs(zp)));
            }
        Assert.True(worst < 4.0, $"{worst:F2}");
    }

    // E[ln v] and P(v < 1) of the conditional on u = ln v, ∝ p(eᵘ)·eᵘ·e^(−bKu/2)·exp(−b·E_r·e⁻ᵘ/2 − b·E_γ·eᵘ/2), by the trapezoid on a fine grid.
    private static (double MeanLog, double Below) Quadrature(IMixingLaw law, double er, double eg, int k, double b)
    {
        const double lo = -40.0, hi = 40.0;
        const int steps = 400_000;
        double h = (hi - lo) / steps, peak = double.NegativeInfinity;
        var log = new double[steps + 1];
        for (int i = 0; i <= steps; i++)
        {
            double u = lo + i * h, v = Math.Exp(u);
            double value = law.LogPrior(v) + u - 0.5 * b * k * u - 0.5 * b * er / v - 0.5 * b * eg * v;
            log[i] = double.IsFinite(value) ? value : double.NegativeInfinity;
            peak = Math.Max(peak, log[i]);
        }
        double mass = 0.0, first = 0.0, below = 0.0;
        for (int i = 0; i <= steps; i++)
        {
            double u = lo + i * h, w = Math.Exp(log[i] - peak) * (i == 0 || i == steps ? 0.5 : 1.0);
            mass += w;
            first += w * u;
            if (u < 0.0) below += w;
            else if (u == 0.0) below += 0.5 * w;
        }
        return (first / mass, below / mass);
    }

    private static double BatchError(double[] x, int batches)
    {
        int size = x.Length / batches;
        var means = new double[batches];
        for (int b = 0; b < batches; b++) means[b] = x.Skip(b * size).Take(size).Average();
        double m = means.Average(), v = means.Sum(y => (y - m) * (y - m)) / (batches - 1);
        return Math.Sqrt(v / batches);
    }
}
