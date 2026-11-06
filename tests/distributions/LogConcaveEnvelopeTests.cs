using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Distributions.Tests;

/// <summary>
/// The tangent hat behind every truncated draw, at its design points. Each family's own log ratio, at parameters over many orders of
/// magnitude, is drawn with reaches from the root of h = −1 out to far past it and to +∞, against its law by quadrature, and the hat's
/// acceptance read at each; and the families' draws are pinned, so a change to the envelope that moves any draw shows. The seeds are
/// fixed, so a failure repeats.
/// </summary>
public sealed class LogConcaveEnvelopeTests
{
    private const int CasesPerFamily = 30;
    private const int Draws = 2_000;
    private const double KolmogorovCut = 2.6;   // P(√n·D > 2.6) ≈ 3 in a million per case

    // The reach handed to the envelope, as a multiple of the root of h = −1 on each side.
    private static readonly double[] Multiples = [1.0, 10.0, 1e3, 1e6, 1e12, 1e100, double.PositiveInfinity];

    private static readonly string[] Families =
        ["Normal", "StudentT", "Beta", "Gamma", "GIG", "MHN linear", "MHN edge", "MHN log", "GIG off the peak"];

    private readonly ITestOutputHelper _out;
    public LogConcaveEnvelopeTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void FarReaches_DrawEachFamilysLaw()
    {
        var rng = new Xoshiro256PlusPlus(20261106);
        var tally = new Dictionary<(string, int), Tally>();
        var failures = new List<string>();
        for (int c = 0; c < CasesPerFamily; c++)
            foreach (string family in Families) Case(rng, family, tally, failures);

        _out.WriteLine("family, reach multiple: cases, throws, failures, acceptance min / median");
        foreach (string family in Families)
            for (int m = 0; m < Multiples.Length; m++)
            {
                Tally t = tally[(family, m)];
                double[] a = [.. t.Acceptance.Order()];
                string read = a.Length == 0 ? "-" : $"{a[0]:F3} / {a[a.Length / 2]:F3}";
                _out.WriteLine($"{family}, ×{Multiples[m]:G3}: {t.Cases}, {t.Throws}, {t.Failures}, {read}");
            }
        Assert.True(failures.Count == 0, string.Join("\n", failures.Take(20)));
    }

    // One case of a family's log ratio on spans that keep its peak at δ = 0: a side the ratio rises toward is cut there.
    private static void Case(Xoshiro256PlusPlus rng, string family, Dictionary<(string, int), Tally> tally, List<string> failures)
    {
        double Long() => rng.NextDouble() < 0.6 ? double.PositiveInfinity : LogUniform(rng, 1e-4, 1e4);
        (double, double) Spans(double slope) => slope < 0.0 ? (0.0, Long()) : slope > 0.0 ? (Long(), 0.0) : (Long(), Long());
        switch (family)
        {
            case "Normal":
            {
                double z = rng.NextDouble() < 0.3 ? 0.0 : Signed(rng) * LogUniform(rng, 1e-3, 1e8);
                var h = new Normal.LogRatio(z);
                (double left, double right) = Spans(h.Slope(0.0));
                Run($"{family} z {z:R}", family, h, left, right, rng, tally, failures);
                break;
            }
            case "StudentT":
            {
                double nu = LogUniform(rng, 0.1, 1e6), small = rng.NextDouble() < 0.3 ? 0.5 : LogUniform(rng, 1e-15, 0.5);
                bool above = rng.NextDouble() < 0.5;
                double up = above ? small : 1.0 - small, down = above ? 1.0 - small : small;
                var h = new StudentT.AsinhRatio(nu, down - up, up, down);
                (double left, double right) = Spans(h.Slope(0.0));
                Run($"{family} ν {nu:R}, (1 − |τ|)/2 {small:R}", family, h, left, right, rng, tally, failures);
                break;
            }
            case "Beta":
            {
                double a = LogUniform(rng, 0.02, 1e5), b = LogUniform(rng, 0.02, 1e5), x0 = a / (a + b), c0 = b / (a + b);
                if (rng.NextDouble() < 0.7)
                {
                    double y = Signed(rng) * LogUniform(rng, 1e-3, 30.0);
                    x0 = 1.0 / (1.0 + Math.Exp(-y));
                    c0 = 1.0 / (1.0 + Math.Exp(y));
                }
                var h = new Beta.LogitRatio(a, b, x0, c0);
                (double left, double right) = Spans(h.Slope(0.0));
                Run($"{family}({a:R}, {b:R}) at {x0:R}", family, h, left, right, rng, tally, failures);
                break;
            }
            case "Gamma":
            {
                double shape = LogUniform(rng, 1e-6, 1e6), g0 = rng.NextDouble() < 0.3 ? shape : LogUniform(rng, 1e-6, 1e6);
                var h = new Gamma.LogGammaRatio(shape, g0);
                (double left, double right) = Spans(h.Slope(0.0));
                Run($"{family} shape {shape:R}, g₀ {g0:R}", family, h, left, right, rng, tally, failures);
                break;
            }
            case "GIG":
            case "GIG off the peak":
            {
                double up = LogUniform(rng, 1e-12, 1e6), down = LogUniform(rng, 1e-12, 1e6);
                double d = family == "GIG" ? 0.0 : Signed(rng) * LogUniform(rng, 1e-3, 30.0);
                var h = new GeneralizedInverseGaussian.LogRatio(up, down, d);
                (double left, double right) = Spans(h.Slope(0.0));
                Run($"{family} A {up:R}, B {down:R}, d {d:R}", family, h, left, right, rng, tally, failures);
                break;
            }
            case "MHN linear":
            {
                double power = LogUniform(rng, 1e-6, 1e4), beta = LogUniform(rng, 1e-6, 1e6), mode = LogUniform(rng, 1e-6, 1e6);
                Run($"{family} α − 1 {power:R}, β {beta:R}, x* {mode:R}", family, new ModifiedHalfNormal.LinearRatio(power, beta, mode),
                    mode, double.PositiveInfinity, rng, tally, failures);
                break;
            }
            case "MHN edge":
            {
                double beta = LogUniform(rng, 1e-6, 1e6), gamma = rng.NextDouble() < 0.2 ? 0.0 : -LogUniform(rng, 1e-6, 1e6);
                Run($"{family} β {beta:R}, γ {gamma:R}", family, new ModifiedHalfNormal.EdgeRatio(beta, gamma),
                    0.0, double.PositiveInfinity, rng, tally, failures);
                break;
            }
            default:
            {
                double q = LogUniform(rng, 1e-8, 1e6), l = rng.NextDouble() < 0.2 ? 0.0 : -LogUniform(rng, 1e-8, 1e6);
                Run($"{family} q {q:R}, l {l:R}", family, new ModifiedHalfNormal.LogRatio(q, l),
                    double.PositiveInfinity, double.PositiveInfinity, rng, tally, failures);
                break;
            }
        }
    }

    private static void Run<T>(string name, string family, T h, double left, double right, Xoshiro256PlusPlus rng,
                               Dictionary<(string, int), Tally> tally, List<string> failures) where T : struct, IConcaveLogRatio
    {
        var law = new Law<T>(h, left, right);
        var deltas = new double[Draws];
        for (int m = 0; m < Multiples.Length; m++)
        {
            if (!tally.TryGetValue((family, m), out Tally? t)) tally[(family, m)] = t = new Tally();
            t.Cases++;
            double leftReach = left > 0.0 ? law.LeftRoot * Multiples[m] : 0.0, rightReach = right > 0.0 ? law.RightRoot * Multiples[m] : 0.0;
            string at = $"{name} on [-{left:R}, {right:R}], reach ×{Multiples[m]:G3}";
            // Over half of the proposals accepted is the envelope's own guarantee, which a design point found holds at any reach.
            try
            {
                double acceptance = law.Mass / LogConcaveEnvelope.HatArea(h, left, right, leftReach, rightReach);
                t.Acceptance.Add(acceptance);
                if (!(acceptance > 0.5)) failures.Add($"acceptance {acceptance:G3}: {at}");
                for (int i = 0; i < Draws; i++)
                {
                    deltas[i] = LogConcaveEnvelope.Draw(rng, h, left, right, leftReach, rightReach);
                    if (!(deltas[i] >= -left && deltas[i] <= right)) throw new ArithmeticException($"δ = {deltas[i]:R} off the interval");
                }
            }
            catch (ArithmeticException e)
            {
                t.Throws++;
                failures.Add($"{e.Message}: {at}");
                continue;
            }
            Array.Sort(deltas);
            double gap = 0.0;
            for (int i = 0; i < Draws; i++)
            {
                double f = law.Cdf(deltas[i]);
                gap = Math.Max(gap, Math.Max(Math.Abs(f - (double)i / Draws), Math.Abs(f - (i + 1.0) / Draws)));
            }
            double statistic = Math.Sqrt(Draws) * gap;
            if (statistic > KolmogorovCut)
            {
                t.Failures++;
                failures.Add($"√n·D = {statistic:F2}: {at}");
            }
        }
    }

    private sealed class Tally
    {
        public int Cases, Throws, Failures;
        public readonly List<double> Acceptance = [];
    }

    // e^h on [−left, right], each side by five-point Gauss–Legendre on even panels out to where concavity puts h below −60: past the
    // root s of h = −1, h(t) ≤ −t/s.
    private sealed class Law<T> where T : struct, IConcaveLogRatio
    {
        private readonly SideLaw<T> _below, _above;

        public Law(T h, double left, double right)
        {
            _below = new SideLaw<T>(h, -1.0, left);
            _above = new SideLaw<T>(h, 1.0, right);
        }

        public double LeftRoot => _below.Root;
        public double RightRoot => _above.Root;
        public double Mass => _below.Mass + _above.Mass;

        public double Cdf(double delta)
            => (delta < 0.0 ? _below.Mass - _below.Within(-delta) : _below.Mass + _above.Within(delta)) / Mass;
    }

    private sealed class SideLaw<T> where T : struct, IConcaveLogRatio
    {
        private const int Panels = 800;
        private static readonly double[] Nodes = [0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640];
        private static readonly double[] Weights = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];
        private readonly T _h;
        private readonly double _sign, _width;
        private readonly double[] _cumulative = new double[Panels + 1];

        public SideLaw(T h, double sign, double end)
        {
            _h = h;
            _sign = sign;
            Root = RootOf(h, sign, end);
            double reach = Root < end ? Math.Min(end, 60.0 * Root) : end;
            _width = reach / Panels;
            for (int k = 0; k < Panels; k++) _cumulative[k + 1] = _cumulative[k] + Integral(k * _width, (k + 1) * _width);
        }

        public double Root { get; }
        public double Mass => _cumulative[Panels];

        public double Within(double t)
        {
            if (!(_width > 0.0)) return 0.0;
            int k = Math.Min((int)(t / _width), Panels - 1);
            return t >= Panels * _width ? Mass : _cumulative[k] + Integral(k * _width, t);
        }

        private double Integral(double a, double b)
        {
            double mid = 0.5 * (a + b), half = 0.5 * (b - a), sum = 0.0;
            for (int i = 0; i < Nodes.Length; i++) sum += Weights[i] * Math.Exp(_h.Value(_sign * (mid + half * Nodes[i])));
            return sum * half;
        }

        // Bisection for h = −1, or the end where h never falls that far; an h that overflows (NaN) counts as past it.
        private static double RootOf(T h, double sign, double end)
        {
            if (!(end > 0.0)) return 0.0;
            bool Past(double t) => !(h.Value(sign * t) > -1.0);
            if (double.IsFinite(end) && !Past(end)) return end;
            double lo = 0.0, hi = Math.Min(1.0, end);
            while (!Past(hi))
            {
                lo = hi;
                hi = Math.Min(2.0 * hi, end);
            }
            for (int i = 0; i < 200; i++)
            {
                double mid = 0.5 * (lo + hi);
                if (Past(mid)) hi = mid;
                else lo = mid;
            }
            return hi;
        }
    }

    // The families' draws hashed: their truncated draws on the sweep's intervals, the GIG's on and off its tails, the modified
    // half-normal's on its branches. Each case draws from its own stream, so a change that moves one case's draws moves that case's
    // hash alone, and the readout names it. A hash moves only when a change moves a draw; regenerating it is a deliberate act.
    // Baselined with the design point's guard (2026-11-05): before it the hash was cc7b23455512c2f9, and the guard moved 31 of the
    // 1,311 cases, all the modified half-normal's, each one whose own reach left Newton's search out of steps short of h = −1: 25 in
    // its branch in x at α from just above 1 to 2.2, where the search met x = 0, and 6 in its ln x branch at α near 10⁻³, from a
    // reach on the exponential side. No other family's draw moved.
    [Fact]
    public void Draws_ArePinned()
    {
        var rng = new Xoshiro256PlusPlus(20261107);
        var cases = new List<Func<Xoshiro256PlusPlus, double>>();
        for (int c = 0; c < 400; c++)
        {
            (_, IContinuousDistribution d) = DistributionSweepTests.Family(rng);
            if (DistributionSweepTests.Interval(rng, d) is { } box && double.IsFinite(d.LogIntervalProbability(box.Lower, box.Upper)))
                cases.Add(r => d.SampleOn(r, box.Lower, box.Upper));
        }
        for (int c = 0; c < 200; c++)
        {
            var gig = GeneralizedInverseGaussian.Of(Signed(rng) * LogUniform(rng, 1e-3, 100.0), LogUniform(rng, 1e-8, 1e4), LogUniform(rng, 1e-8, 1e4));
            double x = gig.Mode, above = x * Math.Exp(3.0 * rng.NextDouble()), below = x * Math.Exp(-3.0 * rng.NextDouble());
            cases.Add(r => gig.Sample(r));
            cases.Add(r => gig.SampleOn(r, above, double.PositiveInfinity));
            cases.Add(r => gig.SampleOn(r, 0.0, below));
            var mhn = ModifiedHalfNormal.Of(LogUniform(rng, 1e-3, 1e3), LogUniform(rng, 1e-3, 1e3), Signed(rng) * LogUniform(rng, 1e-3, 100.0));
            cases.Add(r => mhn.Sample(r));
            // α just above 1 with γ < 0: a mode near 0, where h stays above −1 until x = 0 itself.
            var near = ModifiedHalfNormal.Of(1.0 + LogUniform(rng, 1e-8, 1e-2), LogUniform(rng, 1e-3, 1e3), -LogUniform(rng, 1e-3, 100.0));
            cases.Add(r => near.Sample(r));
        }
        ulong total = Offset;
        for (int c = 0; c < cases.Count; c++)
        {
            var stream = new Xoshiro256PlusPlus(20261108 + c);
            ulong hash = Offset;
            try
            {
                for (int i = 0; i < 25; i++) hash = Add(hash, cases[c](stream));
            }
            catch (ArithmeticException)
            {
                hash = 0;
            }
            _out.WriteLine($"case {c}: {hash:x16}");
            total = Add(total, BitConverter.Int64BitsToDouble((long)hash));
        }
        _out.WriteLine($"draws' hash {total:x16}");
        Assert.Equal("ce13801abaff15d3", total.ToString("x16"));
    }

    // Found by a sweep of the GIG's own truncated draws, 4 of 4,000 cases at small |p|, a and b on far upper tails: the GIG's reach there
    // lands where h is astronomical, and before the guard the hat had no finite area or accepted nothing.
    [Theory]
    [InlineData(-0.003589629176923744, 8.802241158702192E-07, 9.183297216618438E-06, 257.3726870546604)]
    [InlineData(-0.011834994023155954, 8.157473382205912E-07, 4.928133071695612E-07, 2.6209812051973156)]
    [InlineData(-0.0056879131793294005, 3.203977344011842E-08, 0.059015684952721616, 305.6792403399156)]
    public void Gig_OnAFarUpperTail_DrawsItsLaw(double p, double a, double b, double lower)
    {
        var gig = GeneralizedInverseGaussian.Of(p, a, b);
        var law = new Law<GigFromEnd>(new GigFromEnd(p, a * lower, b / lower), 0.0, double.PositiveInfinity);
        var rng = new Xoshiro256PlusPlus(20261109);
        var deltas = new double[Draws];
        for (int i = 0; i < Draws; i++)
        {
            double x = gig.SampleOn(rng, lower, double.PositiveInfinity);
            Assert.InRange(x, lower, double.PositiveInfinity);
            deltas[i] = Math.Log(x / lower);
        }
        Array.Sort(deltas);
        double gap = 0.0;
        for (int i = 0; i < Draws; i++)
        {
            double f = law.Cdf(deltas[i]);
            gap = Math.Max(gap, Math.Max(Math.Abs(f - (double)i / Draws), Math.Abs(f - (i + 1.0) / Draws)));
        }
        Assert.True(Math.Sqrt(Draws) * gap <= KolmogorovCut, $"√n·D = {Math.Sqrt(Draws) * gap:F2}");
    }

    // The GIG's log density in δ = ln(x/x₀) above x₀, from its value there: p·δ − (a·x₀(e^δ − 1) + (b/x₀)(e^(−δ) − 1))/2.
    private readonly struct GigFromEnd(double p, double ax, double bx) : IConcaveLogRatio
    {
        public double Value(double delta) => p * delta - 0.5 * (ax * SpecialFunctions.Expm1(delta) + bx * SpecialFunctions.Expm1(-delta));
        public double Slope(double delta) => p - 0.5 * (ax * Math.Exp(delta) - bx * Math.Exp(-delta));
    }

    private const ulong Offset = 14695981039346656037UL;

    // FNV-1a over the draw's bits.
    private static ulong Add(ulong hash, double x)
    {
        ulong bits = (ulong)BitConverter.DoubleToInt64Bits(x);
        for (int i = 0; i < 8; i++, bits >>= 8) hash = (hash ^ (bits & 0xff)) * 1099511628211UL;
        return hash;
    }

    private static double LogUniform(Xoshiro256PlusPlus rng, double lo, double hi) => DistributionSweepTests.LogUniform(rng, lo, hi);

    private static double Signed(Xoshiro256PlusPlus rng) => rng.NextDouble() < 0.5 ? -1.0 : 1.0;
}
