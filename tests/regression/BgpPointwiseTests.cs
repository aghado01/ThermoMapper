using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.Bgp;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BGP's pointwise density (OQ17, the collapsed content): each term is the leave-unit-out predictive with f, the constant
/// mean and s² integrated out, so it equals the difference of two closed-form evidences, the full data's and the data's
/// without the unit. The sampler scores every draw by its (t, η), PSIS-LOO reads the result and WAIC refuses it.
/// </summary>
public sealed class BgpPointwiseTests
{
    private readonly ITestOutputHelper _out;
    public BgpPointwiseTests(ITestOutputHelper output) => _out = output;

    private static (double[,] X, double[] Y) Data(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[n, 1];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i, 0] = (i + 0.5) / n;
            double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
            y[i] = 3.0 + Math.Sin(2.0 * Math.PI * x[i, 0]) + 0.2 * Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }
        return (x, y);
    }

    private static GpRegression Gp(double[,] x, double[] y, bool unknownScale) => unknownScale
        ? GpRegression.WithUnknownScale(x, y, new SquaredExponentialKernel(), GpInputScaling.None)
        : new GpRegression(x, y, 0.04, new SquaredExponentialKernel(), GpInputScaling.None);

    private static GpFit FitAt(GpRegression gp, double t, double eta) =>
        gp.NoiseModel == GpNoiseModel.Known ? gp.Fit(t) : gp.Fit(t, eta);

    public static IEnumerable<object[]> Partitions()
    {
        foreach (bool unknownScale in new[] { false, true })
        {
            yield return new object[] { unknownScale, "each" };
            yield return new object[] { unknownScale, "blocks" };
            yield return new object[] { unknownScale, "interleaved" };
        }
    }

    private static UnitPartition Partition(string kind, int n) => kind switch
    {
        "each" => UnitPartition.EachObservation(n),
        "blocks" => UnitPartition.ConsecutiveBlocks(n / 4, 4),
        _ => new UnitPartition(Enumerable.Range(0, n).Select(i => i % 6).ToArray()),
    };

    [Theory]
    [MemberData(nameof(Partitions))]
    public void Terms_AreEvidenceDifferences(bool unknownScale, string kind)
    {
        (double[,] x, double[] y) = Data(36, seed: 7);
        GpRegression gp = Gp(x, y, unknownScale);
        const double t = 0.02, eta = 0.08;
        GpFit fit = FitAt(gp, t, eta);
        UnitPartition units = Partition(kind, 36);
        double[] terms = gp.LeaveOutLogDensities(fit, units);
        Assert.Equal(units.Units, terms.Length);

        for (int u = 0; u < units.Units; u++)
        {
            int[] rest = Enumerable.Range(0, 36).Except(units.Members[u]).ToArray();
            var xRest = new double[rest.Length, 1];
            for (int k = 0; k < rest.Length; k++) xRest[k, 0] = x[rest[k], 0];
            GpRegression gpRest = Gp(xRest, rest.Select(i => y[i]).ToArray(), unknownScale);
            double expected = fit.LogMarginal - FitAt(gpRest, t, eta).LogMarginal;
            Assert.True(Math.Abs(terms[u] - expected) < 1e-8 * Math.Max(1.0, Math.Abs(expected)),
                        $"unit {u}: {terms[u]:R} against the evidence difference {expected:R}");
        }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Sampler_ScoresEachDrawByItsFit_ForPsisLooOnly(bool unknownScale)
    {
        (double[,] x, double[] y) = Data(40, seed: 8);
        GpRegression gp = Gp(x, y, unknownScale);
        var sampler = new BgpSampler(gp, new EmpiricalBayesBandwidthPrior(gp));
        BgpResult r = sampler.Run(x, draws: 400, burn: 150, seed: 4, chains: 2);
        PointwiseLogDensity density = sampler.PointwiseLogDensity(r);
        Assert.Equal(r.Draws, density.Draws);
        Assert.Equal(40, density.Units.Units);
        Assert.True(density.LeavesUnitsOut);
        Assert.Equal(unknownScale ? 5 : 3, density.Latents.Count);

        int row = 0;
        foreach (IReadOnlyList<BgpDraw> chain in r.BandwidthDraws)
            foreach (BgpDraw d in chain)
            {
                if (row % 53 == 0)
                {
                    double[] terms = gp.LeaveOutLogDensities(d.Fit, density.Units);
                    for (int i = 0; i < 40; i++) Assert.Equal(terms[i], density.Values[row, i]);
                }
                row++;
            }

        Assert.Throws<ArgumentException>(() => PredictiveChecks.Waic(density));
        PredictiveCheck loo = PredictiveChecks.PsisLoo(density);
        Assert.True(double.IsNaN(loo.EffectiveParameters));
        _out.WriteLine($"elpd {loo.Elpd:F3} (se {loo.StandardError:F3}); max k̂ {loo.ParetoK!.Max():F3}");
        // A noise sd of 0.2 puts a good fit's per-unit predictive density near −½log(2π·0.04) − ½ ≈ 0.19.
        Assert.InRange(loo.Elpd / 40, -0.3, 0.6);
    }

    /// <summary>
    /// With a known σ², t is the only sampled parameter, so exact LOO is a one-dimensional quadrature:
    /// log p(yᵢ | y₋ᵢ) = log ∫p(Y | t)p(t)dt − log ∫p(Y₋ᵢ | t)p(t)dt, independent of the chain. PSIS-LOO on the collapsed
    /// terms matches it per unit, and that includes the units its k̂ flags. Those flags come from clumping: the terms vary
    /// only with t, the chain's slow part, so the largest ratios come from a few excursions and the rejections' repeats.
    /// The Monte Carlo errors hold too. With the chains' relative efficiency the errors against exact LOO are about one
    /// reported error each (rms 0.5–1.3 over five runs); read as independent draws, the reported errors are 2.3–3 times
    /// smaller.
    /// </summary>
    [Fact]
    public void PsisLoo_MatchesExactLoo_ByQuadratureOverTheBandwidth()
    {
        const int n = 40;
        (double[,] x, double[] y) = Data(n, seed: 8);
        GpRegression gp = Gp(x, y, unknownScale: false);
        var prior = new EmpiricalBayesBandwidthPrior(gp);
        var sampler = new BgpSampler(gp, prior);
        PointwiseLogDensity density = sampler.PointwiseLogDensity(sampler.Run(x, draws: 1000, burn: 200, seed: 4, chains: 2));
        PredictiveCheck loo = PredictiveChecks.PsisLoo(density);
        PredictiveCheck asIndependent = PredictiveChecks.PsisLoo(density with { ChainLengths = null });

        // The midpoint rule in log t over the prior's support, with the Jacobian t; the step cancels in each difference.
        const int grid = 800;
        double lo = Math.Log(prior.LowerBound);
        double[] logT = Enumerable.Range(0, grid).Select(g => lo * (1.0 - (g + 0.5) / grid)).ToArray();
        double[] logPrior = logT.Select(u => prior.LogDensity(Math.Exp(u)) + u).ToArray();
        double LogEvidence(GpRegression g) =>
            PredictiveChecks.LogSumExp(logT.Select((u, k) => g.Fit(Math.Exp(u)).LogMarginal + logPrior[k]).ToArray());
        double full = LogEvidence(gp);
        var delta = new double[n];
        for (int i = 0; i < n; i++)
        {
            int[] rest = Enumerable.Range(0, n).Where(j => j != i).ToArray();
            var xRest = new double[n - 1, 1];
            for (int k = 0; k < n - 1; k++) xRest[k, 0] = x[rest[k], 0];
            double exact = full - LogEvidence(Gp(xRest, rest.Select(j => y[j]).ToArray(), unknownScale: false));
            delta[i] = loo.PerUnit[i] - exact;
        }
        double[] flagged = Enumerable.Range(0, n).Where(i => loo.ParetoK![i] > 0.7).Select(i => delta[i]).ToArray();
        _out.WriteLine($"Σ Δ {delta.Sum():F4}; max |Δᵢ| {delta.Max(Math.Abs):F5}; {flagged.Length} units with k̂ > 0.7, |Δ| {string.Join(", ", flagged.Select(d => Math.Abs(d).ToString("F5")))}");
        Assert.All(delta, d => Assert.True(Math.Abs(d) < 0.03, $"Δᵢ {d:F5}"));
        Assert.True(Math.Abs(delta.Sum()) < 0.15, $"Σ Δ {delta.Sum():F4}");
        Assert.All(flagged, d => Assert.True(Math.Abs(d) < 0.01, $"a flagged unit's Δ {d:F5}"));

        double Z(PredictiveCheck c) => Math.Sqrt(delta.Zip(c.MonteCarloErrors!, (d, e) => d / e * (d / e)).Average());
        _out.WriteLine($"rms error over the reported error: {Z(loo):F2} with r_eff (mean {loo.RelativeEfficiency!.Average():F2}), {Z(asIndependent):F2} as independent");
        Assert.InRange(Z(loo), 0.3, 1.5);
        Assert.True(Z(asIndependent) > 2.0 * Z(loo), $"{Z(asIndependent):F2} against {Z(loo):F2}");
    }

    [Fact]
    public void Units_MustLeaveEnoughOutside()
    {
        (double[,] x, double[] y) = Data(12, seed: 9);
        GpRegression known = Gp(x, y, unknownScale: false), unknown = Gp(x, y, unknownScale: true);
        GpFit knownFit = known.Fit(0.05), unknownFit = unknown.Fit(0.05, 0.1);
        Assert.Throws<ArgumentException>(() => known.LeaveOutLogDensities(knownFit, UnitPartition.EachObservation(11)));
        Assert.Throws<ArgumentException>(() => known.LeaveOutLogDensities(knownFit, UnitPartition.ConsecutiveBlocks(1, 12)));
        // A known σ² needs one observation outside the unit; an unknown scale needs two.
        var allButOne = new UnitPartition(Enumerable.Range(0, 12).Select(i => i < 11 ? 0 : 1).ToArray());
        Assert.All(known.LeaveOutLogDensities(knownFit, allButOne), v => Assert.True(double.IsFinite(v)));
        Assert.Throws<ArgumentException>(() => unknown.LeaveOutLogDensities(unknownFit, allButOne));
    }
}
