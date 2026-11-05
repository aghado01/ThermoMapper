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
/// The efficiency pass's lever 5, the marginal's constants (<c>efficiency-pass-brief.md</c>): <see cref="GPriorMarginal.LogMarginal(InverseGamma?, Sigma2PriorTerms)"/>
/// reads <c>g₁ln g₂</c>, <c>lnΓ(g₁ + ν/2)</c> and <c>lnΓ(g₁)</c> from a memo instead of computing them at every evaluation, in the same
/// expression and order as <see cref="GPriorMarginal.LogMarginal(InverseGamma?)"/>, so the same bits. Checked against it over chosen cases and a
/// seeded random parameter sweep, far ranges included, with one memo shared across all of them so that its slots are overwritten, its
/// entries reused across a moving rate, and its misses and hits both read.
/// </summary>
public sealed class BassMarginalTermsTests
{
    private readonly ITestOutputHelper _out;
    public BassMarginalTermsTests(ITestOutputHelper output) => _out = output;

    private static double LogUniform(Xoshiro256PlusPlus rng, double low, double high) => Math.Exp(Math.Log(low) + (Math.Log(high) - Math.Log(low)) * rng.NextDouble());

    private static bool SameBits(double a, double b) => (double.IsNaN(a) && double.IsNaN(b)) || BitConverter.DoubleToInt64Bits(a) == BitConverter.DoubleToInt64Bits(b);

    [Fact]
    public void LogMarginal_WithTheMemo_IsTheDirectFormBitForBit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 9200);
        var terms = new Sigma2PriorTerms();
        var shapes = new List<double> { 0.5, 1.0, 2.0, 3.0, 1e-3, 1e3, 1e8 };
        var scales = new List<double> { 0.3, 1.0, 1e-6, 1e6 };
        int proper = 0, improper = 0, rate = 0;
        for (int trial = 0; trial < 60000; trial++)
        {
            int rows = rng.NextDouble() < 0.5 ? 1 + (int)(rng.NextUInt64() % 2000) : 1 + (int)(rng.NextUInt64() % 5_000_000);
            double b = rng.NextDouble() < 0.5 ? 1.0 : (rng.NextDouble() < 0.1 ? 1e-12 : 1e-3 + 0.999 * rng.NextDouble());
            int spent = rng.NextDouble() < 0.2 ? 1 : 0;
            double s = LogUniform(rng, 1e-9, 1e9), logDet = 100.0 * (rng.NextDouble() - 0.5);
            double tau = LogUniform(rng, 1e-6, 1e6);
            var marginal = new GPriorMarginal(s, logDet, rows, tau, b, spent);
            InverseGamma? prior = null;
            if (rng.NextDouble() < 0.85)
            {
                double shape = rng.NextDouble() < 0.5 ? shapes[(int)(rng.NextUInt64() % (ulong)shapes.Count)] : LogUniform(rng, 1e-3, 1e4);
                double scale = rng.NextDouble() < 0.5 ? scales[(int)(rng.NextUInt64() % (ulong)scales.Count)] : LogUniform(rng, 1e-8, 1e8);
                prior = InverseGamma.ShapeScale(shape, scale);
                if (trial % 3 == 0) rate++;
            }
            double want, got;
            try { want = marginal.LogMarginal(prior); }
            catch (InvalidOperationException)
            {
                Assert.Throws<InvalidOperationException>(() => marginal.LogMarginal(prior, terms));
                continue;
            }
            got = marginal.LogMarginal(prior, terms);
            Assert.True(SameBits(want, got), $"trial {trial}: {want:R} against {got:R} (rows {rows}, b {b:R}, shape {prior?.Shape:R}, scale {prior?.Scale:R})");
            // Asked again, from the memo's hit.
            Assert.True(SameBits(want, marginal.LogMarginal(prior, terms)));
            if (prior is null) improper++;
            else proper++;
        }
        _out.WriteLine($"{proper} proper and {improper} improper evaluations bit-identical, {rate} of them under a moving scale");
        Assert.True(proper > 40000 && improper > 3000);
    }

    [Fact]
    public void CollectionWeight_WithTheMemo_IsTheDirectFormBitForBit()
    {
        var rng = new Xoshiro256PlusPlus(seed: 9300);
        var data = BassCost.Data(BassCost.Base with { N = 200, P = 4 });
        foreach ((double shape, double rate) in new[] { (10.0, 10.0), (1.0, 1.0), (0.5, 3.0), (3.0, 0.01), (1e3, 50.0), (2.5, 1e6), (LogUniform(rng, 0.1, 100), LogUniform(rng, 0.1, 100)) })
        {
            BassConfig config = BassCost.Config(BassCost.Base with { N = 200, P = 4 }) with { Count = new BassCountPrior(shape, rate) };
            BassStructurePrior prior = BassModel.Create(config, data).Prior;
            for (int pass = 0; pass < 2; pass++)   // the second pass reads what the first kept
                foreach (int terms in Enumerable.Range(0, 1100).Concat([2000, 100000, 1023, 1024, 1025]).OrderBy(_ => rng.NextUInt64()))
                {
                    double want = SpecialFunctions.LogGamma(shape + terms) - SpecialFunctions.LogGamma(shape) - terms * Math.Log(rate + 1.0);
                    Assert.True(SameBits(want, prior.LogCollectionWeight(terms)), $"shape {shape}, rate {rate}, M {terms}");
                }
            Assert.Throws<ArgumentOutOfRangeException>(() => prior.LogCollectionWeight(-1));
        }
    }

    [Fact]
    public void LogMarginal_WithTheMemo_RefusesWhereTheDirectFormDoes()
    {
        var terms = new Sigma2PriorTerms();
        var exact = new GPriorMarginal(0.0, 0.0, 10, 1.0, 1.0);
        Assert.Throws<InvalidOperationException>(() => exact.LogMarginal(null));
        Assert.Throws<InvalidOperationException>(() => exact.LogMarginal(null, terms));
        var empty = new GPriorMarginal(1.0, 0.0, 10, 1.0, 0.0);
        Assert.Equal(0.0, empty.LogMarginal(InverseGamma.ShapeScale(2.0, 0.5)));
        Assert.Equal(0.0, empty.LogMarginal(InverseGamma.ShapeScale(2.0, 0.5), terms));
    }
}
