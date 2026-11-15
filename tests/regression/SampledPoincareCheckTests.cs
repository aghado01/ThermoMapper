using System;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E8's slice 8a′, the sampled-sensitivity brief's decision 7: a scalar fit holding a term without one-variable factors reads
/// decision 6's per-input check from Sobol's sampled first-order totals. The totals are Sobol's own at the readout's N and a seed derived
/// apart from C's, bit for bit, with their errors; the bound is the input's constant times C_vv; the hand-built directional fit meets it
/// on both degrees; draws read alone or repeated keep their rows; and on a fit near the uniform law's extremal cos(πu) the ratio is
/// within a few percent of one, where a constant 5% short would be counted.
/// </summary>
public sealed class SampledPoincareCheckTests(ITestOutputHelper output)
{
    /// <summary>
    /// The sampled-sensitivity chip's hand-built directional fit, degrees 1 and 2: the check's totals and their errors are Sobol's sampled
    /// route's at first order, N = 8192 and the seed SeedTree.Derive(seed, 1)[0], bit for bit; each input's constant is the uniform law's,
    /// exact, and Bound/C_vv equals it; no violation.
    /// </summary>
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    public void DirectionalFit_TotalsAreSobolsSampled_AndMeetTheBound(int degree)
    {
        BassData data = SampledSensitivityTests.DirectionalData();
        BassFit fit = SampledSensitivityTests.DirectionalFit(data, degree);
        ActiveSubspaceAnalysis a = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Points = 8192, Seed = 6 });
        Assert.Equal(ActiveSubspaceRoute.Sampled, a.Route);
        PoincareCheck check = a.Poincare!;
        SobolDecomposition s = fit.Sensitivity(data, new SensitivitySpec
        {
            Route = SensitivityRoute.Sampled, MaxOrder = 1, Points = 8192, Seed = SeedTree.Derive(6, 1)[0],
        }).Decompositions[0];
        double worstRead = 0.0, worstRatio = 0.0;
        for (int i = 0; i < a.Inputs.Count; i++)
        {
            Assert.Equal(PoincareConstantKind.Exact, check.Constants![i].Kind);
            for (int d = 0; d < a.Posterior.Draws.Count; d++)
            {
                Assert.Equal(BitConverter.DoubleToInt64Bits(s.TotalVariance[d, i]), BitConverter.DoubleToInt64Bits(check.Total[d, i]));
                Assert.Equal(BitConverter.DoubleToInt64Bits(s.TotalVarianceError![d, i]), BitConverter.DoubleToInt64Bits(check.TotalError![d, i]));
                double cvv = a.Posterior.Draws[d].C[i, i];
                if (cvv > 0.0) worstRead = Math.Max(worstRead, Math.Abs(check.Bound[d, i] / cvv - check.Constants[i].Value) / check.Constants[i].Value);
                worstRatio = Math.Max(worstRatio, check.Ratio[d, i]);
            }
        }
        output.WriteLine($"degree {degree}: totals Sobol's bit for bit; worst |Bound/C_vv − constant|/constant {worstRead:E2}; largest ratio {worstRatio:F3}; violations {check.Violations}");
        Assert.True(worstRead <= 1e-15, $"{worstRead:E2}");
        Assert.Null(check.JacobiRatio);
        Assert.Equal(0, check.Violations);
    }

    /// <summary>Draws asked as [1, 0, 1]: the repeated draw's rows are the same bits, and a draw's row is what the whole readout gives
    /// it, its totals and its errors.</summary>
    [Fact]
    public void Draws_ReadAloneOrRepeated_KeepTheirRows()
    {
        BassData data = SampledSensitivityTests.DirectionalData();
        BassFit fit = SampledSensitivityTests.DirectionalFit(data, 1);
        var spec = new ActiveSubspaceSpec { Points = 2048, Seed = 9 };
        PoincareCheck whole = fit.ActiveSubspace(data, spec).Poincare!;
        PoincareCheck asked = fit.ActiveSubspace(data, spec with { Draws = [1, 0, 1] }).Poincare!;
        for (int i = 0; i < whole.Total.GetLength(1); i++)
        {
            Assert.Equal(asked.Total[0, i], asked.Total[2, i]);
            Assert.Equal(whole.Total[1, i], asked.Total[0, i]);
            Assert.Equal(whole.Total[0, i], asked.Total[1, i]);
            Assert.Equal(whole.TotalError![1, i], asked.TotalError![0, i]);
            Assert.Equal(whole.Bound[0, i], asked.Bound[1, i]);
        }
    }

    /// <summary>
    /// The check's power: one draw of 16 single-variable directional hinges on x1, the piecewise linear interpolant of cos(πu) on x1's
    /// box, u its box coordinate, the uniform law's extremal, under the default uniform law. At N = 65,536 its x1 ratio is within 3% of
    /// one and no violation is counted; with the constant 5% short the same total, bound and errors would be counted.
    /// </summary>
    [Fact]
    public void NearTheExtremal_TheRatioIsNearOne_AndAShortConstantWouldFire()
    {
        BassData data = SampledSensitivityTests.DirectionalData();
        BassFit fit = CosineFit(data, 16);
        ActiveSubspaceAnalysis a = fit.ActiveSubspace(data, new ActiveSubspaceSpec { Points = 65536, Seed = 3 });
        PoincareCheck check = a.Poincare!;
        const int x1 = 1;
        double t = check.Total[0, x1], bound = check.Bound[0, x1], constant = check.Constants![x1].Value;
        double toBound = bound / a.Posterior.Draws[0].C[x1, x1];
        double seBound = toBound * a.MonteCarloError![0][x1, x1], seTotal = check.TotalError![0, x1];
        double ratio = t / bound, slack = 3.0 * Math.Sqrt(seTotal * seTotal + seBound * seBound) / bound;
        double shortBound = 0.95 * bound, shortSlack = 3.0 * Math.Sqrt(seTotal * seTotal + 0.95 * 0.95 * seBound * seBound);
        output.WriteLine($"x1: total/bound {ratio:F4}, three errors {slack:P2} of the bound; with the constant 5% short, total/bound {t / shortBound:F4} " +
                         $"against 1 + {shortSlack / shortBound:P2}; the other inputs' totals {check.Total[0, 0]:E2}, {check.Total[0, 2]:E2}; constant {constant:F6}");
        Assert.InRange(ratio, 0.97, 1.0 + slack);
        Assert.Equal(0, check.Violations);
        Assert.True(t > shortBound + shortSlack, "a constant 5% short is not counted");
    }

    // One draw: Σ_j Δs_j (u − u_j)₊ for the interpolant of cos(πu) at u_j = j/J, each (u − u_j)₊ a single-variable directional hinge on x1's
    // standardized projection z = (x1 − mean)/sd, (u − u_j)₊ = (z − z_j)₊·sd/W, W x1's box width, the hinge over its normalizer c_j.
    private static BassFit CosineFit(BassData data, int pieces)
    {
        BassFit run = SampledSensitivityTests.DirectionalFit(data, 1);
        var scale = (RangeScale)run.Transform.Predictors[1];
        double width = scale.Max - scale.Min;
        (double mean, double deviation) = DirectionalHingeTests.Moments(data)[1];
        var direction = new ProjectionDirection([1], [1.0]);
        double[][] probe = [[0.5], [0.0], [11.0]];
        var table = new BassTermTable();
        var ids = new int[pieces];
        var beta = new double[pieces + 1];
        beta[0] = 1.0;
        double previous = 0.0;
        for (int j = 0; j < pieces; j++)
        {
            double u = (double)j / pieces, slope = (Math.Cos(Math.PI * (j + 1) / pieces) - Math.Cos(Math.PI * u)) * pieces;
            double knot = (scale.Min + u * width - mean) / deviation;
            ids[j] = table.Intern(new BassTerm(directional: new DirectionalHinge(direction, knot)));
            double c = DirectionalHingeTests.Hand(data, probe, direction, knot, 1, HingeNormalization.Box).C;
            beta[j + 1] = (slope - previous) * deviation / width * c;
            previous = slope;
        }
        var store = new DrawStore<BassDraw>();
        store.Add(new BassDraw(new BassStructure(ids), 1.0, 1.0, beta, 1.0));
        return run with { Terms = table, Draws = [store], DrawWeights = null };
    }
}
