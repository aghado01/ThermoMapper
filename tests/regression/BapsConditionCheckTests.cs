using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.RegularExpressions;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The samplers' condition check as a setting of three values (<see cref="ConditionCheck"/>; shared-code defects brief,
/// slice 6).
/// <list type="bullet">
/// <item>The three modes draw the same chain bit for bit over the cost measurement's 19 cases.</item>
/// <item>The screen admits only what the estimate admits: on every sampler system, across designs from well identified to
/// nearly unidentified and coefficient ratios up to the refusal and past it, chosen and seeded at random. The slack the
/// sweep sees against the bound is written to the output.</item>
/// <item>Every mode refuses an unidentified design before the first draw. Every factorization and screened refuse a
/// coefficient drawn far out at the same step, which reference only admits.</item>
/// <item>A resumed run matches the uninterrupted one in each mode, its records included.</item>
/// </list>
/// </summary>
public sealed class BapsConditionCheckTests(ITestOutputHelper output)
{
    private static readonly ConditionCheck[] Modes = [ConditionCheck.EveryFactorization, ConditionCheck.Screened, ConditionCheck.ReferenceOnly];

    // ── Fixtures, as the cost measurement's probe builds them ──────────────────────────────────────────────────────

    private static double[,] Design1D(int nu, double[] x) =>
        new SplineBasis(3).Design(new KnotConfig([.. Enumerable.Range(1, nu - 4).Select(k => k / (nu - 3.0))]), x);

    private static double[] Grid(int n) => [.. Enumerable.Range(0, n).Select(i => (i + 0.5) / n)];

    // x confined to a width δ around 0.5: as δ → 0 the design stops identifying the linear part of the penalty's null space.
    private static double[] Confined(int n, double width) => [.. Grid(n).Select(u => 0.5 + width * (u - 0.5))];

    private static double[] Response(double[] x, int seed, double shift = 0.0)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        return [.. x.Select(v => Math.Sin(2 * Math.PI * v) + 0.3 * Math.Cos(7 * v) + shift + 0.2 * Normal.SampleStandard(rng))];
    }

    // The spectral-tensor report's surface; y's spread around 0.5 scaled by `width`, 0 leaving y's linear part unidentified.
    private static (TensorDesign Td, double[] Y) Tensor(int nx, int ny, double width = 1.0)
    {
        int gx = Math.Max(20, 2 * nx), gy = Math.Max(20, 2 * ny);
        var xs = new List<double>();
        var vs = new List<double>();
        var y = new List<double>();
        for (int a = 0; a < gx; a++)
            for (int b = 0; b < gy; b++)
            {
                double u = (a + 0.5) / gx, v = 0.5 + width * ((b + 0.5) / gy - 0.5);
                xs.Add(u);
                vs.Add(v);
                y.Add(Math.Sin(1.5 * Math.PI * u) * Math.Sin(5 * Math.PI * v) + 0.1 * Math.Cos(37 * u + 53 * v) + 0.05 * Math.Sin(89 * u - 71 * v));
            }
        KnotConfig Knots(int k) => new([.. Enumerable.Range(1, k - 4).Select(j => j / (k - 3.0))]);
        return (TensorProductDesign.Build(new SplineBasis(3), Knots(nx), [.. xs], new SplineBasis(3), Knots(ny), [.. vs]), [.. y]);
    }

    private static double[][] Replicates(double[] x, int reps, int seed) =>
        [.. Enumerable.Range(0, reps).Select(r => Response(x, seed + 31 * r, 0.1 * (r - reps / 2)))];

    private static IEnumerable<double> Values(BapsResult r)
    {
        foreach (IReadOnlyList<BapsDraw> chain in r.PosteriorDraws)
            foreach (BapsDraw d in chain)
            {
                foreach (double v in d.Coefficients) yield return v;
                if (d.ConditionalMean is { } mean) foreach (double v in mean) yield return v;
                yield return d.Lambda;
                yield return d.Sigma2;
                yield return d.Tau2;
                if (d.Eta is { } eta) foreach (double v in eta) yield return v;
                if (d.ComponentTau2 is { } c) foreach (double v in c) yield return v;
            }
    }

    private static IEnumerable<double> Values(HierarchicalResult r)
    {
        foreach (IReadOnlyList<HierarchicalDraw> chain in r.PosteriorDraws)
            foreach (HierarchicalDraw d in chain)
            {
                foreach (double v in d.Population) yield return v;
                foreach (double v in d.PopulationConditionalMean) yield return v;
                if (d.Replicates is { } reps) foreach (double[] rep in reps) foreach (double v in rep) yield return v;
                yield return d.NoiseVariance;
                yield return d.DeviationVariance;
                yield return d.PopulationSmoothingVariance;
            }
    }

    private static long[] Bits(IEnumerable<double> values) => [.. values.Select(BitConverter.DoubleToInt64Bits)];

    // The cost measurement's 19 cases, each a run under a mode: 20 burn-in steps and 30 draws on one chain.
    public static TheoryData<int> Cases
    {
        get
        {
            var cases = new TheoryData<int>();
            for (int i = 0; i < 19; i++) cases.Add(i);
            return cases;
        }
    }

    private static (string Name, Func<ConditionCheck, (long[] Bits, IReadOnlyList<ConditionRecord> Records)> Run) Case(int index)
    {
        int i = 0;
        foreach (int nu in new[] { 12, 48, 192 })
            foreach (int ratio in new[] { 2, 10, 50 })
                if (i++ == index)
                {
                    double[] x = Grid(nu * ratio);
                    double[,] z = Design1D(nu, x);
                    double[] y = Response(x, 20261114 + 100 * nu + ratio);
                    return ($"Gibbs nu={nu} n/nu={ratio}", m => RunBaps(new Baps(z, y, new DifferencePenalty(2), conditionCheck: m), BapsScheme.Gibbs));
                }
        {
            double[] x = Grid(480);
            double[,] z = Design1D(48, x);
            double[] y = Response(x, 20261115);
            if (i++ == index) return ("LocallyAdaptiveGibbs nu=48", m => RunBaps(new Baps(z, y, new DifferencePenalty(2), conditionCheck: m), BapsScheme.LocallyAdaptiveGibbs));
            if (i++ == index)
                return ("Collapsed nu=48", m => RunBaps(new Baps(z, y, new DifferencePenalty(2), noisePrior: NoiseVariancePrior.Jeffreys, conditionCheck: m), BapsScheme.Collapsed));
        }
        foreach (int nu in new[] { 24, 96 })
        {
            double[] x = Grid(5 * nu);
            double[,] z = Design1D(nu, x);
            double[][] ys = Replicates(x, 8, 20261116 + nu);
            foreach (HierarchicalScheme scheme in new[] { HierarchicalScheme.Centred, HierarchicalScheme.CollapsedPopulation, HierarchicalScheme.ParameterExpanded })
                if (i++ == index)
                    return ($"Hierarchical {scheme} nu={nu}", m =>
                    {
                        var model = new HierarchicalPSpline(z, ys, new DifferencePenalty(2),
                            deviationPrior: scheme == HierarchicalScheme.ParameterExpanded ? new HalfTDeviationPrior() : null, conditionCheck: m);
                        HierarchicalResult r = model.Run(chains: 1, masterSeed: 7, burn: 20, samples: 30, scheme: scheme);
                        return (Bits(Values(r)), r.ConditionRecords!);
                    });
        }
        foreach ((int nx, int ny) in new[] { (8, 8), (16, 16) })
            if (i++ == index)
            {
                (TensorDesign td, double[] y) = Tensor(nx, ny);
                return ($"AnisotropicGibbs {nx}x{ny}", m => RunBaps(new Baps(td.Design, y, new TensorPenalty(nx, ny, 2, 2), conditionCheck: m), BapsScheme.AnisotropicGibbs));
            }
        throw new ArgumentOutOfRangeException(nameof(index));

        static (long[], IReadOnlyList<ConditionRecord>) RunBaps(Baps model, BapsScheme scheme)
        {
            BapsResult r = model.Run(scheme, chains: 1, masterSeed: 7, burn: 20, samples: 30);
            return (Bits(Values(r)), r.ConditionRecords!);
        }
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void TheThreeModes_DrawTheSameChain_AndRecordWhatTheirChecksRan(int index)
    {
        (string name, var run) = Case(index);
        var (everyBits, every) = run(ConditionCheck.EveryFactorization);
        var (screenedBits, screened) = run(ConditionCheck.Screened);
        var (referenceBits, reference) = run(ConditionCheck.ReferenceOnly);
        Assert.Equal(everyBits, screenedBits);
        Assert.Equal(everyBits, referenceBits);
        Assert.Equal(every.Count, screened.Count);
        Assert.Equal(every.Count, reference.Count);
        for (int s = 0; s < every.Count; s++)
        {
            ConditionRecord e = every[s], c = screened[s], r = reference[s];
            Assert.Equal(e.System, c.System);
            Assert.Equal(e.Factorizations, c.Factorizations);
            Assert.Equal(e.Factorizations, r.Factorizations);
            Assert.True(e.Factorizations > 0);
            Assert.Equal(e.WidestRatio, c.WidestRatio);   // the ratio is read in every mode, on the same coefficients
            Assert.Equal(e.WidestRatio, r.WidestRatio);

            Assert.True(double.IsNaN(e.ReferenceCondition));
            Assert.Equal(e.Factorizations, e.Estimates);
            Assert.True(e.LargestCondition is >= 1.0 and <= BapsAdmission.MaxCondition);

            Assert.True(c.ReferenceCondition is >= 1.0 and <= BapsAdmission.MaxCondition);
            Assert.Equal(c.ReferenceCondition, r.ReferenceCondition);
            Assert.InRange(c.Estimates, 0, c.Factorizations);
            if (c.Estimates > 0) Assert.True(c.LargestCondition <= e.LargestCondition);

            Assert.Equal(0, r.Estimates);
            Assert.True(double.IsNaN(r.LargestCondition));
            output.WriteLine(string.Create(CultureInfo.InvariantCulture,
                $"{name} | {e.System}: κ_ref {c.ReferenceCondition:G3}, widest ratio {e.WidestRatio:G3}, largest κ̂ {e.LargestCondition:G3}, screened estimates {c.Estimates}/{c.Factorizations}"));
        }
    }

    // ── The invariant: the screen admits only what the estimate admits ─────────────────────────────────────────────

    // What one system's sweep saw.
    private sealed class Sweep(string form)
    {
        public long Steps, Cleared, FellThrough, Refused, NotPositiveDefinite;
        public double Slack, AdmittedSlack, LeastS = double.PositiveInfinity, MostS;

        public void Visit(SamplerAdmission admission, ISpdFactor factor, double spread)
        {
            Steps++;
            bool cleared = !admission.Estimates(spread);
            if (!factor.Status.PositiveDefinite)
            {
                Assert.False(cleared, $"{form}: the screen cleared a step at ratio {spread:G3} whose factor is not positive definite.");
                NotPositiveDefinite++;
                return;
            }
            double condition = factor.EquilibratedCondition();
            bool read = double.IsFinite(admission.ReferenceCondition);   // a reference that is not positive definite clears nothing
            double slack = read ? condition / (admission.ReferenceCondition * spread) : double.NaN;
            if (double.IsFinite(slack)) Slack = Math.Max(Slack, slack);
            if (cleared)
            {
                Cleared++;
                Assert.True(condition <= BapsAdmission.MaxCondition,
                    $"{form}: the screen cleared ratio {spread:G3} at κ_ref {admission.ReferenceCondition:G3}, and the estimate reads {condition:G3}.");
            }
            else if (condition <= BapsAdmission.MaxCondition) FellThrough++;
            else Refused++;
            if (condition <= BapsAdmission.MaxCondition && double.IsFinite(slack))
            {
                AdmittedSlack = Math.Max(AdmittedSlack, slack);
                double s = admission.ScreenScale / admission.ReferenceCondition;
                LeastS = Math.Min(LeastS, s);
                MostS = Math.Max(MostS, s);
                Assert.True(slack <= s, $"{form}: slack {slack:G3} past S = {s:G3}.");
            }
        }

        public override string ToString() => string.Create(CultureInfo.InvariantCulture,
            $"{form}: {Steps} systems, cleared {Cleared}, fell through and admitted {FellThrough}, refused {Refused}, not PD {NotPositiveDefinite}; " +
            $"largest κ̂/(κ_ref·ratio) {Slack:G3}, on admitted {AdmittedSlack:G3}; S from {LeastS:G3} to {MostS:G3}");
    }

    private static double[,] GramBand(double[,] z, int bandwidth)
    {
        int n = z.GetLength(0), nu = z.GetLength(1);
        var band = new double[bandwidth + 1, nu];
        for (int j = 0; j < nu; j++)
            for (int d = 0; d <= bandwidth && j + d < nu; d++)
            {
                double s = 0.0;
                for (int i = 0; i < n; i++) s += z[i, j + d] * z[i, j];
                band[d, j] = s;
            }
        return band;
    }

    private static int GramBandwidth(double[,] z)
    {
        int n = z.GetLength(0), nu = z.GetLength(1), most = 0;
        for (int i = 0; i < n; i++)
        {
            int first = -1, last = -1;
            for (int j = 0; j < nu; j++)
                if (z[i, j] != 0.0) { if (first < 0) first = j; last = j; }
            if (first >= 0) most = Math.Max(most, last - first);
        }
        return most;
    }

    private static BandCholesky Factor(double[,] band, int nu, int bandwidth)
    {
        var factor = new BandCholesky(nu, bandwidth, BandFactorization.Ldlt);
        factor.DecomposeBanded(band);
        return factor;
    }

    private static readonly double[] Widths = [1.0, 0.3, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6];

    // A seeded draw of log-uniform values: the random half of each sweep.
    private static double LogUniform(Xoshiro256PlusPlus rng, double lowExponent, double highExponent) =>
        Math.Pow(10.0, lowExponent + (highExponent - lowExponent) * rng.NextDouble());

    [Fact]
    public void ZtZPlusLambdaP_TheScreenAdmitsOnlyWhatTheEstimateAdmits()
    {
        var sweep = new Sweep("ZᵀZ + λP");
        var rng = new Xoshiro256PlusPlus(20261114);
        IEnumerable<(int Nu, double Width, double[] Ratios)> Systems()
        {
            foreach (int nu in new[] { 12, 48 })
                foreach (double width in Widths)
                    yield return (nu, width, [.. Enumerable.Range(-12, 33).Select(k => Math.Pow(10.0, k))]);
            for (int t = 0; t < 40; t++)
                yield return (8 + (int)(rng.NextDouble() * 60), LogUniform(rng, -6, 0), [.. Enumerable.Range(0, 8).Select(_ => LogUniform(rng, -12, 20))]);
        }
        foreach ((int nu, double width, double[] ratios) in Systems())
        {
            double[] x = Confined(10 * nu, width);
            double[,] z = Design1D(nu, x);
            var penalty = new DifferencePenalty(2);
            SamplerAdmission admission = new Baps(z, Response(x, nu), penalty, conditionCheck: ConditionCheck.Screened).Admission(BapsScheme.Gibbs);
            int bw = Math.Max(GramBandwidth(z), penalty.Bandwidth);
            double[,] gram = GramBand(z, bw);
            double lambda0 = admission.ReferenceWeights[1];
            foreach (double ratio in ratios)
            {
                var band = (double[,])gram.Clone();
                penalty.AccumulateInto(band, nu, ratio * lambda0);
                sweep.Visit(admission, Factor(band, nu, bw), admission.Spread([1.0, ratio * lambda0]));
            }
        }
        output.WriteLine(sweep.ToString());
    }

    [Fact]
    public void LocallyAdaptive_TheScreenAdmitsOnlyWhatTheEstimateAdmits()
    {
        var sweep = new Sweep("ZᵀZ + Dᵀdiag(λᵢ)D");
        var rng = new Xoshiro256PlusPlus(20261115);
        const int nu = 24;
        var penalty = new DifferencePenalty(2);
        int terms = penalty.DifferenceCount(nu);
        foreach (double width in new[] { 1.0, 1e-2, 1e-4, 1e-6 })
        {
            double[] x = Confined(10 * nu, width);
            double[,] z = Design1D(nu, x);
            SamplerAdmission admission = new Baps(z, Response(x, 3), penalty, conditionCheck: ConditionCheck.Screened)
                .Admission(BapsScheme.LocallyAdaptiveGibbs);
            int bw = Math.Max(GramBandwidth(z), penalty.Bandwidth);
            double[,] gram = GramBand(z, bw);
            double lambda0 = admission.ReferenceWeights[1];
            foreach (double spreadOfEta in new[] { 0.0, 1.0, 4.0, 12.0 })
                for (int k = -12; k <= 20; k += 2)
                    for (int draw = 0; draw < 3; draw++)
                    {
                        var weights = new double[terms];
                        double scale = Math.Pow(10.0, k) * lambda0;
                        for (int i = 0; i < terms; i++) weights[i] = scale * Math.Exp(spreadOfEta * Normal.SampleStandard(rng));
                        var band = (double[,])gram.Clone();
                        penalty.AccumulateInto(band, nu, weights);
                        sweep.Visit(admission, Factor(band, nu, bw), admission.Spread([1.0, .. weights]));
                    }
        }
        output.WriteLine(sweep.ToString());
    }

    [Fact]
    public void Anisotropic_TheScreenAdmitsOnlyWhatTheEstimateAdmits()
    {
        var sweep = new Sweep("ZᵀZ + Σⱼ λⱼPⱼ");
        var rng = new Xoshiro256PlusPlus(20261116);
        int[] exponents = [-8, -4, 0, 4, 8, 12, 16];
        foreach (double width in new[] { 1.0, 1e-1, 1e-2, 1e-3, 1e-4, 0.0 })
        {
            (TensorDesign td, double[] y) = Tensor(8, 8, width);
            var penalty = new TensorPenalty(8, 8, 2, 2);
            const int nu = 64;
            SamplerAdmission admission = new Baps(td.Design, y, penalty, conditionCheck: ConditionCheck.Screened).Admission(BapsScheme.AnisotropicGibbs);
            int bw = Math.Max(GramBandwidth(td.Design), penalty.Bandwidth);
            double[,] gram = GramBand(td.Design, bw);
            double rx = admission.ReferenceWeights[1], ry = admission.ReferenceWeights[2];
            var pairs = new List<(double, double)>();
            foreach (int a in exponents) foreach (int b in exponents) pairs.Add((Math.Pow(10.0, a), Math.Pow(10.0, b)));
            for (int t = 0; t < 30; t++) pairs.Add((LogUniform(rng, -12, 20), LogUniform(rng, -12, 20)));
            foreach ((double ax, double ay) in pairs)
            {
                var band = (double[,])gram.Clone();
                penalty.AccumulateComponents(band, nu, [ax * rx, ay * ry]);
                sweep.Visit(admission, Factor(band, nu, bw), admission.Spread([1.0, ax * rx, ay * ry]));
            }
        }
        output.WriteLine(sweep.ToString());
    }

    [Fact]
    public void HierarchicalSystems_TheScreenAdmitsOnlyWhatTheEstimateAdmits()
    {
        const int nu = 24, reps = 4;
        var penalty = new DifferencePenalty(2);
        var rng = new Xoshiro256PlusPlus(20261117);
        var replicate = new Sweep("ZᵀZ/σ² + I/τ_u²");
        var population = new Sweep("(R/τ_u²)I + P/τ_0²");
        var collapsed = new Sweep("β_0's collapsed precision");
        var deviation = new Sweep("ξ²ZᵀZ/σ² + I/σ_η²");
        var expanded = new Sweep("R·ZᵀZ/σ² + P/τ_0²");
        double[] exponents = [.. Enumerable.Range(-10, 11).Select(k => 2.0 * k)];
        foreach (double width in new[] { 1.0, 1e-2, 1e-4, 1e-6 })
        {
            double[] x = Confined(5 * nu, width);
            double[,] z = Design1D(nu, x);
            double[][] ys = Replicates(x, reps, 5);
            var centred = new HierarchicalPSpline(z, ys, penalty, conditionCheck: ConditionCheck.Screened);
            var halfT = new HierarchicalPSpline(z, ys, penalty, deviationPrior: new HalfTDeviationPrior(), conditionCheck: ConditionCheck.Screened);
            SamplerAdmission?[] c = centred.Admissions(HierarchicalScheme.Centred);
            SamplerAdmission?[] k = centred.Admissions(HierarchicalScheme.CollapsedPopulation);
            SamplerAdmission?[] e = halfT.Admissions(HierarchicalScheme.ParameterExpanded);
            int gbw = GramBandwidth(z), bw = Math.Max(gbw, penalty.Bandwidth);
            double[,] gram = GramBand(z, bw);
            var p = new double[bw + 1, nu];
            penalty.AccumulateInto(p, nu, 1.0);

            // Two coefficients each: the ratio's whole range, chosen, and seeded draws.
            IEnumerable<(double A, double B)> Pairs()
            {
                foreach (double a in exponents) yield return (1.0, Math.Pow(10.0, a));
                for (int t = 0; t < 20; t++) yield return (LogUniform(rng, -10, 10), LogUniform(rng, -20, 20));
            }
            double[,] Sum(double ca, double[,] a, double cb, double[,] b)
            {
                var band = new double[bw + 1, nu];
                for (int d = 0; d <= bw; d++)
                    for (int j = 0; j < nu; j++) band[d, j] = ca * a[d, j] + cb * b[d, j];
                return band;
            }
            var identity = new double[bw + 1, nu];
            for (int j = 0; j < nu; j++) identity[0, j] = 1.0;
            foreach ((double a, double b) in Pairs())
            {
                // ZᵀZ/σ² + I/τ_u² at σ² = 1/a, τ_u² = 1/b; ξ²ZᵀZ/σ² + I/σ_η² the same form.
                replicate.Visit(c[HierarchicalPSpline.ReplicateSite]!, Factor(Sum(a, gram, b, identity), nu, bw), c[HierarchicalPSpline.ReplicateSite]!.Spread([a, b]));
                deviation.Visit(e[HierarchicalPSpline.ExpandedDeviationSite]!, Factor(Sum(a, gram, b, identity), nu, bw), e[HierarchicalPSpline.ExpandedDeviationSite]!.Spread([a, b]));
                // (R/τ_u²)I + P/τ_0² and R·ZᵀZ/σ² + P/τ_0².
                population.Visit(c[HierarchicalPSpline.PopulationSite]!, Factor(Sum(a, identity, b, p), nu, bw), c[HierarchicalPSpline.PopulationSite]!.Spread([a, b]));
                expanded.Visit(e[HierarchicalPSpline.ExpandedPopulationSite]!, Factor(Sum(a, gram, b, p), nu, bw), e[HierarchicalPSpline.ExpandedPopulationSite]!.Spread([a, b]));
            }
            // ξ = 0: ZᵀZ's coefficient vanishes and the screen hands the step to the estimate.
            deviation.Visit(e[HierarchicalPSpline.ExpandedDeviationSite]!, Factor(Sum(0.0, gram, 1.0, identity), nu, bw),
                            e[HierarchicalPSpline.ExpandedDeviationSite]!.Spread([0.0, 1.0]));
            Assert.True(e[HierarchicalPSpline.ExpandedDeviationSite]!.Estimates(e[HierarchicalPSpline.ExpandedDeviationSite]!.Spread([0.0, 1.0])));

            // The collapsed precision P/τ_0² + Σₖ R·λₖ/(σ² + τ_u²λₖ)·vₖvₖᵀ, dense, over (σ², τ_u², τ_0²).
            SamplerAdmission dense = k[HierarchicalPSpline.CollapsedSite]!;
            HierarchicalPSpline.CollapsedReference reference = centred.Collapsed();
            var gramDense = new double[nu, nu];
            for (int d = 0; d <= gbw; d++)
                for (int j = 0; j + d < nu; j++) gramDense[j + d, j] = gramDense[j, j + d] = gram[d, j];
            EigenResult eigen = DenseEigen.DecomposeSymmetric(gramDense);
            var pDense = new double[nu, nu];
            for (int d = 0; d <= penalty.Bandwidth; d++)
                for (int j = 0; j + d < nu; j++) pDense[j + d, j] = pDense[j, j + d] = p[d, j];
            var triples = new List<(double, double, double)>();
            foreach (double s2 in new[] { 1e-6, 1.0, 1e6 })
                foreach (double tu in new[] { 1e-6, 1.0, 1e6 })
                    for (int t0 = -12; t0 <= 12; t0 += 4) triples.Add((s2, tu, Math.Pow(10.0, t0)));
            for (int t = 0; t < 20; t++) triples.Add((LogUniform(rng, -8, 8), LogUniform(rng, -8, 8), LogUniform(rng, -16, 16)));
            foreach ((double sigma2, double tauU2, double tau02) in triples)
            {
                var a = new double[nu, nu];
                for (int i = 0; i < nu; i++)
                    for (int j = 0; j < nu; j++) a[i, j] = pDense[i, j] / tau02;
                for (int m = 0; m < nu; m++)
                {
                    double lambda = Math.Max(0.0, eigen.Eigenvalues[m]), w = reps * lambda / (sigma2 + tauU2 * lambda);
                    double[] v = eigen.Eigenvectors[m];
                    for (int i = 0; i < nu; i++)
                        for (int j = 0; j < nu; j++) a[i, j] += w * v[i] * v[j];
                }
                var factor = new DenseCholesky(nu);
                factor.Decompose(a);
                collapsed.Visit(dense, factor, HierarchicalPSpline.CollapsedSpread(reference, sigma2, tauU2, [1.0 / tau02]));
            }
        }
        foreach (Sweep s in new[] { replicate, population, collapsed, deviation, expanded }) output.WriteLine(s.ToString());
    }

    // ── Refusals ──────────────────────────────────────────────────────────────────────────────────────────────────

    // A run's refusal: a chain's arrives through the ensemble's parallel loop, the reference point's directly.
    private static InvalidOperationException Refusal(Action run)
    {
        Exception e = Assert.ThrowsAny<Exception>(run);
        if (e is AggregateException { InnerExceptions: [InvalidOperationException inner] }) return inner;
        return Assert.IsType<InvalidOperationException>(e);
    }

    private static InvalidOperationException? RefusalOrNone(Action run)
    {
        try { run(); return null; }
        catch (AggregateException e) when (e.InnerExceptions is [InvalidOperationException inner]) { return inner; }
        catch (InvalidOperationException e) { return e; }
    }

    [Fact]
    public void EveryMode_RefusesAnUnidentifiedDesign_BeforeTheFirstDraw()
    {
        double[] x = [.. Enumerable.Repeat(0.5, 60)];   // one distinct x: the linear part of the null space is unidentified
        double[,] z = Design1D(12, x);
        double[] y = Response(x, 11);
        (TensorDesign td, double[] ty) = Tensor(8, 8, width: 0.0);
        double[][] ys = Replicates(x, 4, 13);
        foreach (ConditionCheck mode in Modes)
        {
            var refusals = new List<InvalidOperationException>
            {
                Refusal(() =>
                    new Baps(z, y, new DifferencePenalty(2), conditionCheck: mode).Run(BapsScheme.Gibbs, chains: 1, burn: 0, samples: 2)),
                Refusal(() =>
                    new Baps(td.Design, ty, new TensorPenalty(8, 8, 2, 2), conditionCheck: mode).Run(BapsScheme.AnisotropicGibbs, chains: 1, burn: 0, samples: 2)),
                Refusal(() =>
                    new HierarchicalPSpline(z, ys, new DifferencePenalty(2), conditionCheck: mode)
                        .Run(chains: 1, burn: 0, samples: 2, scheme: HierarchicalScheme.CollapsedPopulation)),
                Refusal(() =>
                    new HierarchicalPSpline(z, ys, new DifferencePenalty(2), deviationPrior: new HalfTDeviationPrior(), conditionCheck: mode)
                        .Run(chains: 1, burn: 0, samples: 2, scheme: HierarchicalScheme.ParameterExpanded)),
            };
            if (mode == ConditionCheck.ReferenceOnly)
                foreach (InvalidOperationException refusal in refusals) Assert.Contains("trace-balance reference point", refusal.Message);
        }
    }

    // λ far out: positive definite, with κ(E) past 10¹⁴. Every factorization and screened refuse it, reference only admits it.
    [Fact]
    public void ACoefficientFarOut_IsRefusedByEveryFactorizationAndScreened_AndAdmittedByReferenceOnly()
    {
        double[] x = Grid(120);
        double[,] z = Design1D(12, x);
        double[] y = Response(x, 17);
        var penalty = new DifferencePenalty(2);
        SamplerAdmission[] admissions = [.. Modes.Select(m => new Baps(z, y, penalty, conditionCheck: m).Admission(BapsScheme.Gibbs))];
        double lambda0 = admissions[1].ReferenceWeights[1];
        double far = double.NaN;
        for (int k = 8; k <= 22 && double.IsNaN(far); k++)
        {
            var probe = GramBand(z, 3);
            penalty.AccumulateInto(probe, 12, Math.Pow(10.0, k) * lambda0);
            BandCholesky factor = Factor(probe, 12, 3);
            if (factor.Status.PositiveDefinite && factor.EquilibratedCondition() > BapsAdmission.MaxCondition) far = Math.Pow(10.0, k) * lambda0;
        }
        Assert.False(double.IsNaN(far), "no λ found where the factor is positive definite and κ(E) passes 10¹⁴");
        var tallies = new ConditionTally[3];
        string?[] reasons = new string?[3];
        for (int m = 0; m < 3; m++)
        {
            tallies[m] = new ConditionTally();
            var band = GramBand(z, 3);
            penalty.AccumulateInto(band, 12, far);
            reasons[m] = admissions[m].Admit(Factor(band, 12, 3), admissions[m].Spread([1.0, far]), tallies[m]);
        }
        Assert.Equal(SamplerAdmission.Refused, reasons[0]);
        Assert.NotNull(reasons[1]);
        Assert.Contains("the reference point's is", reasons[1]);
        Assert.Contains("the coefficients' ratio to it", reasons[1]);
        Assert.Null(reasons[2]);
        output.WriteLine($"λ/λ₀ = {far / lambda0:G3}: screened says '{reasons[1]}'");

        // In a chain: the collapsed scheme's wide proposals reach λ far out, and the first two modes refuse the same step.
        var lambdaInMessage = new Regex(@"at λ = (\S+) \(");
        int refusedSeeds = 0;
        for (int seed = 1; seed <= 12; seed++)
        {
            InvalidOperationException? Run(ConditionCheck mode) => RefusalOrNone(() =>
                new Baps(z, y, penalty, mhLogStep: 12.0, noisePrior: NoiseVariancePrior.Jeffreys, conditionCheck: mode)
                    .Run(BapsScheme.Collapsed, chains: 1, masterSeed: seed, burn: 200, samples: 2));
            InvalidOperationException? every = Run(ConditionCheck.EveryFactorization), screened = Run(ConditionCheck.Screened);
            Assert.Equal(every is null, screened is null);
            if (every is null) continue;
            refusedSeeds++;
            string atEvery = lambdaInMessage.Match(every.Message).Groups[1].Value;
            Assert.NotEmpty(atEvery);
            Assert.Equal(atEvery, lambdaInMessage.Match(screened!.Message).Groups[1].Value);
            InvalidOperationException? reference = Run(ConditionCheck.ReferenceOnly);
            if (reference is not null) Assert.Contains("not positive definite at λ", reference.Message);
        }
        Assert.True(refusedSeeds > 0, "no seed's chain reached a λ the gate refuses");
    }

    // ── Resume, and the setting's own checks ──────────────────────────────────────────────────────────────────────

    [Fact]
    public void AResumedRun_MatchesTheUninterruptedOne_InEachMode_RecordsIncluded()
    {
        double[] x = Grid(240);
        double[,] z = Design1D(24, x);
        double[] y = Response(x, 19);
        double[][] ys = Replicates(Grid(120), 4, 23);
        double[,] hz = Design1D(24, Grid(120));
        foreach (ConditionCheck mode in Modes)
        {
            var baps = new Baps(z, y, new DifferencePenalty(2), conditionCheck: mode);
            BapsResult whole = baps.Run(BapsScheme.Gibbs, chains: 2, masterSeed: 5, burn: 15, samples: 40);
            BapsResult first = baps.Run(BapsScheme.Gibbs, chains: 2, masterSeed: 5, burn: 15, samples: 20, checkpoint: true);
            BapsResult resumed = baps.Run(BapsScheme.Gibbs, chains: 2, masterSeed: 5, burn: 15, samples: 20, resumeFrom: first.Checkpoint);
            Assert.Equal(Bits(Values(whole)), Bits(Values(resumed)));
            Assert.Equal(whole.ConditionRecords, resumed.ConditionRecords);
            Assert.NotEqual(first.ConditionRecords![0].Factorizations, whole.ConditionRecords![0].Factorizations);

            var hierarchical = new HierarchicalPSpline(hz, ys, new DifferencePenalty(2), conditionCheck: mode);
            HierarchicalResult h = hierarchical.Run(chains: 2, masterSeed: 5, burn: 15, samples: 40);
            HierarchicalResult hFirst = hierarchical.Run(chains: 2, masterSeed: 5, burn: 15, samples: 20, checkpoint: true);
            HierarchicalResult hResumed = hierarchical.Run(chains: 2, masterSeed: 5, burn: 15, samples: 20, resumeFrom: hFirst.Checkpoint);
            Assert.Equal(Bits(Values(h)), Bits(Values(hResumed)));
            Assert.Equal(h.ConditionRecords, hResumed.ConditionRecords);
        }

        // A checkpoint resumes only under the mode that started it (slice 7): one taken under every factorization, the default
        // before slice 7, is refused by the screened default and resumes when that mode is named; a screened one resumes under
        // the default.
        var pinned = new Baps(z, y, new DifferencePenalty(2), conditionCheck: ConditionCheck.EveryFactorization);
        BapsResult before = pinned.Run(BapsScheme.Gibbs, chains: 2, masterSeed: 5, burn: 15, samples: 20, checkpoint: true);
        var byDefault = new Baps(z, y, new DifferencePenalty(2));
        Assert.Throws<ArgumentException>(() => byDefault.Run(BapsScheme.Gibbs, chains: 2, masterSeed: 5, burn: 15, samples: 20, resumeFrom: before.Checkpoint));
        BapsResult continued = pinned.Run(BapsScheme.Gibbs, chains: 2, masterSeed: 5, burn: 15, samples: 20, resumeFrom: before.Checkpoint);
        BapsResult screened = byDefault.Run(BapsScheme.Gibbs, chains: 2, masterSeed: 5, burn: 15, samples: 20, checkpoint: true);
        BapsResult screenedOn = new Baps(z, y, new DifferencePenalty(2), conditionCheck: ConditionCheck.Screened)
            .Run(BapsScheme.Gibbs, chains: 2, masterSeed: 5, burn: 15, samples: 20, resumeFrom: screened.Checkpoint);
        Assert.Equal(Bits(Values(continued)), Bits(Values(screenedOn)));
    }

    // Slice 7: screened is the default on both samplers, and a default run draws what every factorization draws.
    [Fact]
    public void TheDefault_IsScreened_AndDrawsTheEveryFactorizationChain()
    {
        double[] x = Grid(240);
        double[,] z = Design1D(24, x);
        double[] y = Response(x, 29);
        BapsResult byDefault = new Baps(z, y, new DifferencePenalty(2)).Run(BapsScheme.Gibbs, chains: 2, masterSeed: 3, burn: 20, samples: 30);
        BapsResult every = new Baps(z, y, new DifferencePenalty(2), conditionCheck: ConditionCheck.EveryFactorization)
            .Run(BapsScheme.Gibbs, chains: 2, masterSeed: 3, burn: 20, samples: 30);
        Assert.Equal(Bits(Values(every)), Bits(Values(byDefault)));
        ConditionRecord record = Assert.Single(byDefault.ConditionRecords!);
        Assert.Equal(ConditionCheck.Screened, record.Check);
        Assert.True(double.IsFinite(record.ReferenceCondition));

        double[][] ys = Replicates(Grid(120), 4, 31);
        double[,] hz = Design1D(24, Grid(120));
        HierarchicalResult h = new HierarchicalPSpline(hz, ys, new DifferencePenalty(2)).Run(chains: 2, masterSeed: 3, burn: 20, samples: 30);
        HierarchicalResult hEvery = new HierarchicalPSpline(hz, ys, new DifferencePenalty(2), conditionCheck: ConditionCheck.EveryFactorization)
            .Run(chains: 2, masterSeed: 3, burn: 20, samples: 30);
        Assert.Equal(Bits(Values(hEvery)), Bits(Values(h)));
        Assert.All(h.ConditionRecords!, r => Assert.Equal(ConditionCheck.Screened, r.Check));
    }

    [Fact]
    public void AnUndefinedMode_IsRefused()
    {
        double[] x = Grid(60);
        double[,] z = Design1D(12, x);
        Assert.Throws<ArgumentOutOfRangeException>(() => new Baps(z, Response(x, 1), new DifferencePenalty(2), conditionCheck: (ConditionCheck)3));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new HierarchicalPSpline(z, Replicates(x, 2, 1), new DifferencePenalty(2), conditionCheck: (ConditionCheck)3));
    }
}
