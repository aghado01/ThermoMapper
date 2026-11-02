using System;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;

namespace Maths.Distributions.Tests;

/// <summary>
/// Truncated draws (<see cref="IContinuousDistribution.SampleOn"/>) against the truncated law's own density and interval
/// moments, which the draws do not read: a chi-square over bins whose masses are a composite Gauss–Legendre rule on the log
/// density, relative to the anchor, so a narrow or far-tail bin keeps its digits. Each case lies in the far tail, on a narrow
/// interval or far from the mode. A generator stepped alongside counts the uniforms, so the acceptance rate is measured
/// without a hook in the sampler. R's exact truncated moments are in the oracle project.
/// </summary>
public sealed class TruncatedDrawTests
{
    private const int Draws = 200_000;
    private const int Bins = 20;
    private const double Inf = double.PositiveInfinity;

    // qchisq(1 − 1e-6, 19). Seventy-odd cases at this level leave a false alarm near 1 in 14,000.
    private const double ChiSquareCut = 63.68;

    // The hat's floor as measured over a grid of each family, 86.8% at Gamma(½, 1) on [½, ∞); it proves more than half.
    private const double HatFloor = 0.85;

    public static TheoryData<double, double, double, double> NormalCases => new()
    {
        { 0.0, 1.0, 10.0, Inf },                // the brief's case: 1 − Φ(10) is below ε
        { 0.0, 1.0, 40.0, Inf },
        { 0.0, 1.0, 1000.0, Inf },
        { 0.0, 1.0, double.NegativeInfinity, -40.0 },
        { 0.0, 1.0, 10.0, 10.1 },               // a tenth of a scale width in the tail
        { 0.0, 1.0, 10.0, 10.0001 },            // narrow: flat in a law whose scale is 0.1
        { 0.0, 1.0, 40.0, 40.001 },
        { 0.0, 1.0, -1e-6, 1e-6 },              // narrow around the mean
        { 0.0, 1.0, 0.0, 1e-3 },
        { 0.0, 1.0, -0.5, 0.5 },
        { 0.0, 1.0, -3.0, 4.0 },
        { 0.0, 1.0, 0.0, Inf },
        { 0.0, 1.0, double.NegativeInfinity, Inf },
        { 5.0, 0.1, 0.0, 1.0 },                 // 40–50 sds below the mean
        { -3.0, 0.01, 0.0, Inf },               // a probit latent with the wrong sign, 300 sds out
        { 1e6, 1.0, 1e6 + 1.0, 1e6 + 1.0 + 1e-6 },
    };

    [Theory]
    [MemberData(nameof(NormalCases))]
    public void Normal_MatchesItsDensityAndIntervalMoments(double mean_, double sd_, double lower, double upper)
    {
        Normal parent = Normal.MeanSd(mean_, sd_);
        double[] draws = Check(parent, lower, upper, HatFloor);

        // The interval moments: a series in the density, sharing nothing with the hat.
        var truncated = Truncated.Of(parent, lower, upper);
        double centre = double.IsFinite(lower) ? lower : double.IsFinite(upper) ? upper : mean_;
        var moments = new double[3];
        truncated.IntervalMoments(double.NegativeInfinity, double.PositiveInfinity, centre, 0.0, moments);
        double mean = draws.Sum(x => x - centre) / draws.Length;
        double second = draws.Sum(x => (x - centre) * (x - centre)) / draws.Length;
        double sd = Math.Sqrt(Math.Max(moments[2] - moments[1] * moments[1], 0.0));
        // The first-moment bound holds on the scale of the law's spread; the second on its own.
        Assert.True(Math.Abs(mean - moments[1]) <= 5.0 * sd / Math.Sqrt(draws.Length) + 4e-16 * Math.Abs(centre),
                    $"first moment {mean:R} against {moments[1]:R}, sd {sd:R}");
        double secondSd = Math.Sqrt(draws.Sum(x => Math.Pow((x - centre) * (x - centre) - second, 2.0)) / draws.Length);
        Assert.True(Math.Abs(second - moments[2]) <= 6.0 * secondSd / Math.Sqrt(draws.Length) + 4e-15 * centre * centre,
                    $"second moment {second:R} against {moments[2]:R}");
    }

    [Theory]
    [InlineData(2.0, 1.0, 1000.0, Inf)]         // far above the mode: S(1000) underflows
    [InlineData(2.0, 1.0, 40.0, 41.0)]
    [InlineData(1000.0, 1.0, 0.0, 1.0)]         // far below it: xᵅ⁻¹ alone
    [InlineData(1000.0, 1.0, 800.0, 801.0)]
    [InlineData(1000.0, 1.0, 1000.0, 1000.001)] // narrow at the mode
    [InlineData(1e6, 1.0, 1e6, 1e6 + 1.0)]
    [InlineData(1e6, 1.0, 1.1e6, Inf)]
    [InlineData(0.5, 1.0, 0.0, 1e-8)]           // below 1 the density diverges at 0
    [InlineData(0.01, 1.0, 5.0, 6.0)]
    [InlineData(0.01, 1.0, 0.0, 1e-50)]
    [InlineData(1e-3, 1.0, 1e-200, 1e-100)]
    [InlineData(1e-3, 1.0, 1.0, Inf)]
    [InlineData(1.0, 1.0, 50.0, Inf)]           // the exponential: memoryless
    [InlineData(3.5, 0.01, 2000.0, 2100.0)]     // a small rate: scale 100
    [InlineData(3.5, 100.0, 0.5, Inf)]
    public void Gamma_MatchesItsDensity(double shape, double rate, double lower, double upper)
        => Check(Gamma.ShapeRate(shape, rate), lower, upper, HatFloor);

    [Theory]
    [InlineData(50.0, 1e-3, 1.0, Inf)]          // a variance floor far above the mode: F-R1's s2.lower trap
    [InlineData(2.5, 1e-4, 10.0, Inf)]
    [InlineData(500.0, 1.0, 1.0, Inf)]
    [InlineData(3.0, 2.0, 0.0, 0.01)]           // far below the mode: e^(−β/x) falls fast
    [InlineData(3.0, 2.0, 0.01, 0.0101)]
    [InlineData(0.5, 50.0, 1e-6, 1e-5)]         // Zellner–Siow's g
    [InlineData(0.5, 50.0, 1e6, Inf)]
    [InlineData(1000.0, 1000.0, 1.0, 1.0 + 1e-3)]
    [InlineData(1e-3, 1.0, 1e50, 1e100)]
    [InlineData(1e-3, 1.0, 0.0, 1.0)]
    public void InverseGamma_MatchesItsDensity(double shape, double scale, double lower, double upper)
        => Check(InverseGamma.ShapeScale(shape, scale), lower, upper, HatFloor);

    [Theory]
    [InlineData(2.0, 5.0, 0.999, 1.0)]
    [InlineData(2.0, 5.0, 0.0, 1e-6)]
    [InlineData(2.0, 5.0, 0.2, 0.2001)]
    [InlineData(0.01, 0.01, 0.4, 0.6)]
    [InlineData(0.5, 0.5, 0.0, 1.0)]
    [InlineData(1000.0, 1000.0, 0.9, 1.0)]
    [InlineData(5.0, 0.1, 0.0, 0.5)]
    [InlineData(1e-3, 3.0, 0.5, 1.0)]
    public void Beta_MatchesItsDensity(double a, double b, double lower, double upper)
        => Check(Beta.Shapes(a, b), lower, upper, HatFloor);

    [Theory]
    [InlineData(3.0, 1e6, Inf)]                 // a power law far out
    [InlineData(3.0, 1e6, 1.0001e6)]
    [InlineData(1.0, 100.0, Inf)]               // the Cauchy has no mean
    [InlineData(0.1, 0.0, Inf)]
    [InlineData(0.1, 1e30, 1e40)]
    [InlineData(1.0, -1e-3, 1e-3)]
    [InlineData(30.0, -20.0, -10.0)]
    [InlineData(1e4, 5.0, 6.0)]
    [InlineData(1e8, 8.0, Inf)]                 // nearly normal, where the tails lose digits and the draws read none
    public void StudentT_MatchesItsDensity(double nu, double lower, double upper)
        => Check(StudentT.Standard(nu), lower, upper, HatFloor);

    [Theory]
    [InlineData(5.0, 1000.0)]                   // far above the location, in units of its scale
    [InlineData(1.9, 2.1)]                      // around it
    [InlineData(-100.0, 1.0)]                   // far below, on the negative side of the location
    public void StudentT_MovesWithItsLocationAndScale(double lower, double upper)
        => Check(StudentT.LocationScale(3.0, 2.0, 0.5), lower, upper, HatFloor);

    [Theory]
    [InlineData(0.5, 2.0, 100.0, Inf)]
    [InlineData(0.5, 2.0, 0.0, 1e-6)]
    [InlineData(3.0, 1.0, 3.0, 3.001)]
    [InlineData(1.0, 1.0, 50.0, 60.0)]
    public void Weibull_MatchesItsDensity_InClosedForm(double shape, double scale, double lower, double upper)
        => Check(Weibull.ShapeScale(shape, scale), lower, upper, 1.0);

    [Theory]
    [InlineData(1.0, 0.0, Inf)]
    [InlineData(1.0, 1e6, Inf)]
    [InlineData(1.0, 1e6, 1.001e6)]
    [InlineData(2.0, 0.5, 3.0)]
    [InlineData(2.0, 0.0, 1e-3)]
    public void HalfCauchy_MatchesItsDensity_InClosedForm(double scale, double lower, double upper)
        => Check(HalfCauchy.WithScale(scale), lower, upper, 1.0);

    [Fact]
    public void Uniform_DrawsItsClippedInterval()
    {
        Uniform u = Uniform.On(1.0, 3.0);
        Check(u, 2.0, Inf, 1.0);
        var rng = new Xoshiro256PlusPlus(3);
        double[] all = [.. Enumerable.Range(0, 20000).Select(_ => u.Sample(rng))];
        Assert.InRange(all.Min(), 1.0, 1.01);
        Assert.InRange(all.Max(), 2.99, 3.0);
        Assert.Equal(2.0, all.Average(), 0.03);
        Assert.Throws<ArgumentException>(() => u.SampleOn(rng, 4.0, 5.0));
    }

    [Fact]
    public void Truncated_DrawsFromItsParent_AndIntersectsNestedBounds()
    {
        var truncated = Truncated.Of(Normal.Standard, 10.0, Inf);
        Xoshiro256PlusPlus a = new(5), b = new(5);
        for (int i = 0; i < 100; i++) Assert.Equal(Normal.Standard.SampleOn(b, 10.0, Inf), truncated.Sample(a));
        for (int i = 0; i < 100; i++) Assert.Equal(Normal.Standard.SampleOn(b, 11.0, 12.0), truncated.SampleOn(a, 11.0, 12.0));
        for (int i = 0; i < 100; i++) Assert.InRange(truncated.SampleOn(a, 5.0, 10.5), 10.0, 10.5);
        Assert.Throws<ArgumentException>(() => truncated.SampleOn(a, 0.0, 9.0));
        // Twice truncated, through the interface: the law is the intersection's.
        Check(Truncated.Of(Gamma.ShapeRate(2.0, 1.0), 1.0, 10.0), 5.0, Inf, HatFloor);
    }

    [Fact]
    public void Mixture_PicksByWeightAndMass()
    {
        // On [2, 3]: N(0, 1) has mass 0.0214 at weight 0.9, N(3, 1) has 0.3413 at weight 0.1, so the second holds 0.639.
        Normal a = Normal.Standard, b = Normal.MeanSd(3.0, 1.0);
        Mixture mix = Mixture.Of([a, b], [0.9, 0.1]);
        var rng = new Xoshiro256PlusPlus(11);
        double[] on = [.. Enumerable.Range(0, Draws).Select(_ => mix.SampleOn(rng, 2.0, 3.0))];
        double m0 = a.IntervalProbability(2.0, 3.0), m1 = b.IntervalProbability(2.0, 3.0);
        var ma = new double[2];
        var mb = new double[2];
        a.IntervalMoments(2.0, 3.0, 0.0, 0.0, ma);
        b.IntervalMoments(2.0, 3.0, 0.0, 0.0, mb);
        double expected = (0.9 * ma[1] + 0.1 * mb[1]) / (0.9 * m0 + 0.1 * m1);
        Assert.Equal(expected, on.Average(), 5.0 * 0.29 / Math.Sqrt(Draws));
        double[] all = [.. Enumerable.Range(0, Draws).Select(_ => mix.Sample(rng))];
        Assert.Equal(0.3, all.Average(), 5.0 * 1.7 / Math.Sqrt(Draws));
        Check(Truncated.Of(mix, 2.0, Inf), 2.0, 3.0, 1.0);   // a pick ahead of each draw, so no count of proposals
        Assert.Throws<ArgumentException>(() => Mixture.Of([Uniform.On(0.0, 1.0)], [1.0]).SampleOn(rng, 2.0, 3.0));
    }

    [Fact]
    public void Draws_FollowTheOneStream_Deterministically()
    {
        // Pinned: regenerating these is a deliberate act.
        var rng = new Xoshiro256PlusPlus(20261102);
        double[] pinned = [10.006999867083126, 1004.9956235364975, 1.0525064754544284];
        double[] drawn =
        [
            Normal.Standard.SampleOn(rng, 10.0, Inf),
            Gamma.ShapeRate(2.0, 1.0).SampleOn(rng, 1000.0, Inf),
            InverseGamma.ShapeScale(50.0, 1e-3).SampleOn(rng, 1.0, Inf),
        ];
        for (int i = 0; i < pinned.Length; i++)
            Assert.True(Math.Abs(pinned[i] - drawn[i]) <= 1e-12 * Math.Abs(pinned[i]), $"drew {string.Join(", ", drawn.Select(v => v.ToString("R")))}");

        Xoshiro256PlusPlus a = new(7), b = new(7);
        var cases = new (IContinuousDistribution D, double L, double U)[]
        {
            (Normal.Standard, 3.0, Inf), (Gamma.ShapeRate(0.5, 1.0), 2.0, 3.0), (InverseGamma.ShapeScale(2.0, 1.0), 1.0, 2.0),
            (Beta.Shapes(2.0, 3.0), 0.5, 0.6), (StudentT.Standard(3.0), 2.0, Inf), (Weibull.ShapeScale(2.0, 1.0), 1.0, 2.0), (HalfCauchy.WithScale(1.0), 1.0, 2.0),
        };
        foreach ((IContinuousDistribution d, double l, double u) in cases)
            for (int i = 0; i < 200; i++)
                Assert.Equal(d.SampleOn(a, l, u), d.SampleOn(b, l, u));
        Assert.Equal(a.SaveState(), b.SaveState());
    }

    [Fact]
    public void Draws_ReadNoTail_SoTheyHoldWhereTheTailsCannot()
    {
        // At these sizes the incomplete gamma's series gives up and ln B cancels, so Truncated.Of cannot be built; the draw
        // needs neither. Gamma(10¹², 1) above its mean is a half-normal of sd 10⁶, and likewise the t.
        var rng = new Xoshiro256PlusPlus(20261102);
        Assert.Throws<ArithmeticException>(() => Gamma.ShapeRate(1e12, 1.0).Quantile(0.5));
        double[] gamma = [.. Enumerable.Range(0, Draws).Select(_ => Gamma.ShapeRate(1e12, 1.0).SampleOn(rng, 1e12, Inf) - 1e12)];
        Assert.All(gamma, x => Assert.True(x >= 0.0));
        double halfNormal = 1e6 * Math.Sqrt(2.0 / Math.PI);
        Assert.Equal(halfNormal, gamma.Average(), 5.0 * 1e6 * Math.Sqrt(1.0 - 2.0 / Math.PI) / Math.Sqrt(Draws));

        double[] t = [.. Enumerable.Range(0, Draws).Select(_ => StudentT.Standard(1e14).SampleOn(rng, 0.0, Inf))];
        Assert.Equal(Math.Sqrt(2.0 / Math.PI), t.Average(), 5.0 * Math.Sqrt(1.0 - 2.0 / Math.PI) / Math.Sqrt(Draws));
    }

    [Fact]
    public void Refusals()
    {
        var rng = new Xoshiro256PlusPlus(1);
        foreach (IContinuousDistribution d in new IContinuousDistribution[]
                 {
                     Normal.Standard, Gamma.ShapeRate(2.0, 1.0), InverseGamma.ShapeScale(2.0, 1.0), Beta.Shapes(2.0, 3.0), StudentT.Standard(3.0),
                     Weibull.ShapeScale(2.0, 1.0), HalfCauchy.WithScale(1.0), Uniform.On(0.0, 1.0),
                 })
        {
            Assert.Throws<ArgumentException>(() => d.SampleOn(rng, 1.0, 1.0));
            Assert.Throws<ArgumentException>(() => d.SampleOn(rng, 2.0, 1.0));
            Assert.Throws<ArgumentException>(() => d.SampleOn(rng, double.NaN, 1.0));
            Assert.Throws<ArgumentNullException>(() => d.SampleOn(null!, 0.25, 0.75));
        }
        Assert.Throws<ArgumentException>(() => Gamma.ShapeRate(2.0, 1.0).SampleOn(rng, -2.0, -1.0));
        Assert.Throws<ArgumentException>(() => InverseGamma.ShapeScale(2.0, 1.0).SampleOn(rng, -2.0, 0.0));
        Assert.Throws<ArgumentException>(() => Beta.Shapes(2.0, 3.0).SampleOn(rng, 1.0, 2.0));
    }

    // Draws, then the bounds, the chi-square over bins from the density and the acceptance rate. Returns the draws.
    private static double[] Check<T>(T parent, double lower, double upper, double acceptance) where T : struct, IContinuousDistribution
    {
        var rng = new Xoshiro256PlusPlus(20261102);
        var counter = Xoshiro256PlusPlus.Restore(rng.SaveState().s0, rng.SaveState().s1, rng.SaveState().s2, rng.SaveState().s3);
        var draws = new double[Draws];
        for (int i = 0; i < Draws; i++) draws[i] = parent.SampleOn(rng, lower, upper);

        (double lo, double hi) = Support(parent, lower, upper);
        Assert.All(draws, x => Assert.InRange(x, lo, hi));

        var truncated = Truncated.Of(parent, lower, upper);
        double[] edges = Edges(truncated, lo, hi);
        double[] expected = BinMasses(truncated, edges, parent.Quantile(0.0));
        var counts = new int[Bins];
        foreach (double x in draws)
        {
            int bin = Array.BinarySearch(edges, x);
            if (bin < 0) bin = ~bin - 1;
            counts[Math.Clamp(bin, 0, Bins - 1)]++;
        }
        double chi = 0.0;
        for (int k = 0; k < Bins; k++)
        {
            double e = expected[k] * Draws;
            if (e > 0.0) chi += (counts[k] - e) * (counts[k] - e) / e;
        }
        Assert.True(chi <= ChiSquareCut, $"chi-square {chi:F1} over {Bins} bins on [{lower}, {upper}]; counts {string.Join(",", counts)}");

        // A generator stepped until it meets the sampler's counts its uniforms: three per proposal under the hat.
        if (acceptance < 1.0)
        {
            long uniforms = 0;
            var state = rng.SaveState();
            while (counter.SaveState() != state)
            {
                counter.NextUInt64();
                uniforms++;
            }
            Assert.Equal(0L, uniforms % 3);
            double rate = 3.0 * Draws / uniforms;
            Assert.True(rate >= acceptance - 0.01, $"accepted {rate:P1} on [{lower}, {upper}]");
        }
        return draws;
    }

    private static (double Lo, double Hi) Support<T>(T parent, double lower, double upper) where T : struct, IContinuousDistribution
    {
        double lo = Math.Max(lower, parent.Quantile(0.0)), hi = Math.Min(upper, parent.Quantile(1.0));
        return (lo, hi);
    }

    // Bin edges at the truncated law's quantiles; the last bin may run to an infinite end. Edges need not be exact, since each
    // bin's mass is computed from the density, so a far-tail or narrow interval that the quantile cannot resolve falls back
    // on even spacing.
    private static double[] Edges<T>(Truncated<T> truncated, double lo, double hi) where T : struct, IContinuousDistribution
    {
        var edges = new double[Bins + 1];
        edges[0] = lo;
        edges[Bins] = hi;
        for (int k = 1; k < Bins; k++) edges[k] = truncated.Quantile((double)k / Bins);
        bool resolved = true;
        for (int k = 1; k <= Bins; k++) resolved &= edges[k] > edges[k - 1];
        if (!resolved && double.IsFinite(hi))
            for (int k = 1; k < Bins; k++) edges[k] = lo + (hi - lo) * k / Bins;
        return edges;
    }

    // Each bin's mass: an adaptive rule on the density, relative to its largest value on the interval, in ln x over a bin of
    // decades. A bin that runs to an infinite end, or to the 0 of a support, takes its mass from the log tails, which is the
    // density's integral that a finite rule cannot reach; the two must agree on the normalizer, which is checked.
    private static double[] BinMasses<T>(Truncated<T> truncated, double[] edges, double supportLower) where T : struct, IContinuousDistribution
    {
        var masses = new double[Bins];
        var fromTails = new bool[Bins];
        double reference = double.NegativeInfinity;
        for (int k = 0; k <= Bins; k++)
            if (double.IsFinite(edges[k])) reference = Math.Max(reference, truncated.LogPdf(edges[k]));
        double tails = 0.0, density = 0.0;
        for (int k = 0; k < Bins; k++)
        {
            double a = edges[k], b = edges[k + 1];
            fromTails[k] = !double.IsFinite(a) || !double.IsFinite(b) || (a == 0.0 && supportLower == 0.0);
            if (fromTails[k]) tails += masses[k] = truncated.IntervalProbability(a, b);
            else density += masses[k] = Adaptive(truncated, a, b, reference);
        }
        double scale = (1.0 - tails) / density;
        Assert.True(Math.Abs(scale * Math.Exp(-reference) - 1.0) <= 1e-6, $"the density integrates to {1.0 / (scale * Math.Exp(-reference)):R} of the normalizer");
        for (int k = 0; k < Bins; k++)
            if (!fromTails[k]) masses[k] *= scale;
        return masses;
    }

    private static double Adaptive<T>(Truncated<T> truncated, double a, double b, double reference) where T : struct, IContinuousDistribution
    {
        bool logarithmic = a > 0.0 && b / a > 2.0;
        double f(double u)
        {
            double x = logarithmic ? Math.Exp(u) : u;
            return Math.Exp(truncated.LogPdf(x) - reference) * (logarithmic ? x : 1.0);
        }
        return Panel(f, logarithmic ? Math.Log(a) : a, logarithmic ? Math.Log(b) : b, 80);
    }

    // Splits while the 10- and 20-node rules differ, or while a panel's end stands above what the rules saw, where a spike
    // at the end falls between their nodes.
    private static double Panel(Func<double, double> f, double a, double b, int depth)
    {
        double coarse = Rule(f, a, b, 10), fine = Rule(f, a, b, 20);
        double end = Math.Max(f(a), f(b)) * (b - a);
        if (depth == 0 || (Math.Abs(fine - coarse) <= 1e-10 * fine && end <= 4.0 * fine) || (fine == 0.0 && end == 0.0) || end < 1e-20) return fine;
        double mid = 0.5 * (a + b);
        return Panel(f, a, mid, depth - 1) + Panel(f, mid, b, depth - 1);
    }

    private static double Rule(Func<double, double> f, double a, double b, int nodes)
    {
        GaussLegendre rule = GaussLegendre.Of(nodes);
        double mid = 0.5 * (a + b), half = 0.5 * (b - a), sum = 0.0;
        for (int i = 0; i < rule.Count; i++) sum += rule.Weights[i] * f(mid + half * rule.Nodes[i]);
        return half * sum;
    }
}
