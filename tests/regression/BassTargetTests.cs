using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Xunit;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's collapsed target (stage 4, S5; D03, D05). The shared block's statistics-only marginal and its slope in b;
/// τ's unit-information prior on u = ln τ (F-O8's identities); the prior's pieces and the two modes; the integrated and
/// linear forms and their endpoints; the exact-fit exclusion under Jeffreys; shift invariance under level centring and
/// Liang, and its failure under the zero centre (F-O9); and scale invariance under the scaled σ² prior. The marginals
/// themselves are checked against the stage-4 oracle's quadrature in <c>BassTargetParityTests</c>.
/// </summary>
public sealed class BassTargetTests
{
    // x₀ and x₁ on a five-value grid, three times (N = 75); y = 1 + 2x₀ + 3(x₀ − ½)₊(½ − x₁)₊ + noise, or exactly
    // 1 + 2(x₀ − ½)₊, which one hinge and the intercept fit.
    private static BassData Data(double offset = 0.0, double scale = 1.0, bool exact = false, bool constant = false)
    {
        var rng = new Xoshiro256PlusPlus(seed: 17);
        var x0 = new List<double>();
        var x1 = new List<double>();
        var y = new List<double>();
        for (int rep = 0; rep < 3; rep++)
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++)
                {
                    double a = i * 0.25, b = j * 0.25;
                    x0.Add(a);
                    x1.Add(b);
                    double noise = 0.3 * Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());
                    double value = constant ? 4.0
                        : exact ? 1 + 2 * Math.Max(a - 0.5, 0)
                        : 1 + 2 * a + 3 * Math.Max(a - 0.5, 0) * Math.Max(0.5 - b, 0) + noise;
                    y.Add(offset + scale * value);
                }
        return new BassData([new ContinuousPredictor("x0", [.. x0]), new ContinuousPredictor("x1", [.. x1])], new ScalarResponse([.. y]));
    }

    private sealed record Built(BassTarget Target, BassStructurePrior Prior, Func<BassTerm[], BassStructureFit> Fit, BassResolved Resolved);

    private static Built Build(BassData data, BassConfig config, TemperedForm form = TemperedForm.Integrated)
    {
        BassResolved resolved = BassResolved.Resolve(config, data);
        BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
        var prior = BassStructurePrior.From(config, resolved, inputs);
        var basis = new BassBasis(inputs, config.Structure.Degree);
        BassStructureFit Fit(BassTerm[] terms)
        {
            var design = new BassDesign(basis, data, BassAdmission.From(config, resolved));
            return design.Load(terms) is { } refusal
                ? BassStructureFit.Refused(refusal)
                : BassStructureFit.Of(design, terms.Sum(prior.LogPrior));
        }
        return new(new BassTarget(config, resolved, prior, form), prior, Fit, resolved);
    }

    // Knots index the grid {0, ¼, ½, ¾, 1}.
    private static readonly BassTerm[][] Structures =
    [
        [],
        [new BassTerm([new HingeFactor(0, 1, 2)])],
        [new BassTerm([new HingeFactor(0, 1, 2)]), new BassTerm([new HingeFactor(0, 1, 2), new HingeFactor(1, -1, 2)])],
        [new BassTerm([new HingeFactor(0, -1, 1)]), new BassTerm([new HingeFactor(1, 1, 3)]), new BassTerm([new HingeFactor(0, 1, 3)])],
    ];

    private static BassConfig Config(GPriorCentre centre = GPriorCentre.Mean, Sigma2Prior? sigma2 = null, BassCoefficientPrior? coefficients = null)
        => new() { Intercept = centre, Sigma2 = sigma2 ?? Sigma2Prior.ScaledToResponse(1e-3), Coefficients = coefficients ?? new SampledG() };

    public static TheoryData<GPriorCentre> Centres() => new() { GPriorCentre.Mean, GPriorCentre.Flat, GPriorCentre.Zero };

    [Theory]
    [MemberData(nameof(Centres))]
    public void Slope_IsTheDerivativeOfTheMarginal(GPriorCentre centre)
    {
        InverseGamma?[] priors = [InverseGamma.ShapeScale(1e-3, 1e-3 * 2.4), null, InverseGamma.ShapeScale(2.0, 3.0)];
        foreach (InverseGamma? sigma in priors)
            foreach (double tau in new[] { 0.01, 1.0 })
                foreach (double b in new[] { 0.05, 0.4, 1.0 })
                {
                    double h = 1e-6 * b;
                    double L(double at) => CollapsedGPrior.CentredMarginalAt(3, 4.2, 7.9, 40, 1.3, centre, tau, at).LogMarginal(sigma);
                    double difference = (L(b + h) - L(b - h)) / (2 * h);
                    double slope = CollapsedGPrior.CentredMarginalAt(3, 4.2, 7.9, 40, 1.3, centre, tau, b).LogMarginalSlope(sigma);
                    Assert.True(Math.Abs(slope - difference) <= 1e-6 * Math.Max(1.0, Math.Abs(slope)),
                                $"{centre}, σ² {sigma?.ToString() ?? "Jeffreys"}, τ {tau}, b {b}: {slope:R} against {difference:R}");
                }

        // At b = 0 under a proper prior the marginal is exactly 0, and the slope is the one-sided derivative there,
        // carrying −(k/2)/τ.
        if (centre == GPriorCentre.Flat) return;
        var proper = InverseGamma.ShapeScale(2.0, 3.0);
        double L0(double at) => CollapsedGPrior.CentredMarginalAt(3, 4.2, 7.9, 40, 1.3, centre, 0.5, at).LogMarginal(proper);
        double step = 1e-5;
        double forward = (-3 * L0(0) + 4 * L0(step) - L0(2 * step)) / (2 * step);
        double at0 = CollapsedGPrior.CentredMarginalAt(3, 4.2, 7.9, 40, 1.3, centre, 0.5, 0.0).LogMarginalSlope(proper);
        Assert.Equal(0.0, L0(0));
        Assert.True(Math.Abs(at0 - forward) <= 1e-5 * Math.Abs(at0), $"{centre} at b = 0: {at0:R} against {forward:R}");
    }

    [Theory]
    [MemberData(nameof(Centres))]
    public void CentredMarginal_IsTheBlocksMarginal(GPriorCentre centre)
    {
        BassData data = Data();
        BassConfig config = Config(centre);
        BassResolved resolved = BassResolved.Resolve(config, data);
        BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
        var design = new BassDesign(new BassBasis(inputs, 1), data, BassAdmission.From(config, resolved));
        Assert.Null(design.Load(Structures[3]));
        CollapsedGPrior block = CollapsedGPrior.FromCentredLeastSquares(design.Factor, design.Coefficients.ToArray(), design.ResidualSquares,
                                                                        design.FittedSquares, design.Observations, design.ResponseMean, centre);
        foreach (double tau in new[] { 0.02, 3.0 })
            foreach (double b in new[] { 0.3, 1.0 })
            {
                GPriorMarginal expected = block.MarginalAt(tau, b);
                GPriorMarginal actual = CollapsedGPrior.CentredMarginalAt(design.Count, design.ResidualSquares, design.FittedSquares,
                                                                          design.Observations, design.ResponseMean, centre, tau, b);
                Assert.Equal(expected.S, actual.S);
                Assert.Equal(expected.LogDeterminantTerm, actual.LogDeterminantTerm);
                Assert.Equal(expected.DegreesOfFreedom, actual.DegreesOfFreedom);
                Assert.Equal(expected.LogMarginal(resolved.Sigma2), actual.LogMarginal(resolved.Sigma2));
            }
    }

    [Theory]
    [InlineData(0.5)]
    [InlineData(2.0)]
    public void TauPrior_IsTheUnitInformationGamma(double shape)
    {
        // F-O8's identities on u = ln τ: the density integrates to one, E[τ] = 1/N, and E[ln τ] = ψ(a) − ln(aN).
        Built built = Build(Data(), Config(coefficients: new SampledG(shape)));
        int n = built.Target.Observations;
        double centre = Math.Log(1.0 / n);
        double mass = 0, meanTau = 0, meanLog = 0;
        foreach ((double u, double w) in Nodes(centre - 90 / shape, centre + 8))
        {
            double p = w * Math.Exp(built.Target.LogTauDensity(u));
            mass += p;
            meanTau += p * Math.Exp(u);
            meanLog += p * u;
        }
        Assert.True(Math.Abs(mass - 1.0) <= 1e-13, $"mass {mass:R}");
        Assert.True(Math.Abs(meanTau * n - 1.0) <= 1e-12, $"E[τ]·N = {meanTau * n:R}");
        double expected = SpecialFunctions.Digamma(shape) - Math.Log(shape * n);
        Assert.True(Math.Abs(meanLog - expected) <= 1e-11 * Math.Abs(expected), $"E[ln τ] = {meanLog:R} against {expected:R}");
        Assert.Contains(built.Resolved.AutoPicks, p => p.Setting == "TauPrior");
    }

    [Fact]
    public void Prior_ReadsTheStructureAndTau()
    {
        Built g = Build(Data(), Config());
        BassStructureFit fit = g.Fit(Structures[2]);
        double structure = g.Prior.LogCollectionWeight(2) + Structures[2].Sum(g.Prior.LogPrior);
        BassEvaluation ev = g.Target.Evaluate(fit, Math.Log(0.03));
        Assert.Equal(structure + g.Target.LogTauDensity(Math.Log(0.03)), ev.LogPrior, 12);
        Assert.Equal(0.03, ev.Tau, 15);
        Assert.Throws<InvalidOperationException>(() => g.Target.Evaluate(fit));

        Built fixedG = Build(Data(), Config(coefficients: new FixedG()));
        BassEvaluation held = fixedG.Target.Evaluate(fixedG.Fit(Structures[2]));
        Assert.Equal(structure, held.LogPrior, 12);
        Assert.Equal(1.0 / 75, held.Tau);
        Assert.False(fixedG.Target.SamplesTau);
        Assert.Throws<InvalidOperationException>(() => fixedG.Target.Evaluate(fit, 0.0));
        Assert.Throws<InvalidOperationException>(() => fixedG.Target.LogTauDensity(0.0));
    }

    [Fact]
    public void Forms_AndTheirEndpoints()
    {
        BassData data = Data();
        double u = Math.Log(0.02);
        Built integrated = Build(data, Config());
        Built linear = Build(data, Config(), TemperedForm.Linear);
        BassEvaluation a = integrated.Target.Evaluate(integrated.Fit(Structures[2]), u);
        BassEvaluation l = linear.Target.Evaluate(linear.Fit(Structures[2]), u);

        // b = 1 is one target on both forms; the linear form is b·ℓ₁ with slope ℓ₁.
        Assert.Equal(a.LogLikelihood(1.0), l.LogLikelihood(1.0));
        Assert.Equal(a.LogMarginal, a.LogLikelihood(1.0));
        foreach (double b in new[] { 0.0, 0.2, 0.7 })
        {
            Assert.Equal(b * l.LogMarginal, l.LogLikelihood(b));
            Assert.Equal(l.LogMarginal, l.LogLikelihoodSlope(b));
        }
        Assert.NotEqual(l.LogLikelihood(0.2), a.LogLikelihood(0.2));

        // The proper defaults reach b = 0 on the integrated form: the prior, ℓ₀ = 0.
        Assert.Equal(0.0, a.LogLikelihood(0.0));
        Assert.True(integrated.Target.Defines(0.0));

        // Jeffreys and the flat intercept do not; the linear form does.
        Built jeffreys = Build(data, Config(sigma2: Sigma2Prior.Jeffreys));
        BassEvaluation j = jeffreys.Target.Evaluate(jeffreys.Fit(Structures[2]), u);
        Assert.Throws<ArgumentOutOfRangeException>(() => j.LogLikelihood(0.0));
        Assert.True(double.IsFinite(j.LogLikelihood(1e-3)));
        Built flat = Build(data, Config(GPriorCentre.Flat));
        BassEvaluation f = flat.Target.Evaluate(flat.Fit(Structures[2]), u);
        Assert.Throws<ArgumentOutOfRangeException>(() => f.LogLikelihood(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => f.LogLikelihood(0.5 / 75));   // g₁ + (Nb − 1)/2 < 0
        Assert.True(double.IsFinite(f.LogLikelihood(2.0 / 75)));
        Built flatLinear = Build(data, Config(GPriorCentre.Flat, Sigma2Prior.Jeffreys), TemperedForm.Linear);
        Assert.Equal(0.0, flatLinear.Target.Evaluate(flatLinear.Fit(Structures[2]), u).LogLikelihood(0.0));

        // The powered path powers prior and likelihood alike, τ in its log coordinate.
        var powered = new Temperature(0.3, TemperingPath.Powered);
        Assert.Equal(0.3 * (a.LogPrior + a.LogMarginal), powered.LogDensity(a), 12);
        Assert.Equal(a.LogPrior + a.LogMarginal, Temperature.LogPosterior(a), 12);
    }

    [Fact]
    public void ExactFits_AreRefusedUnderJeffreysOnly()
    {
        BassData exact = Data(exact: true);
        Built jeffreys = Build(exact, Config(sigma2: Sigma2Prior.Jeffreys));
        BassStructureFit fit = jeffreys.Fit(Structures[1]);
        Assert.True(fit.ResidualSquares < 1e-25, $"‖r‖² = {fit.ResidualSquares:E2}");
        BassEvaluation refused = jeffreys.Target.Evaluate(fit, Math.Log(0.02));
        Assert.Equal(BassRefusal.ExactFit, refused.Refusal);
        Assert.Equal(double.NegativeInfinity, refused.LogPrior);
        Assert.Equal(double.NegativeInfinity, refused.LogLikelihood(1.0));
        Assert.Null(jeffreys.Target.Evaluate(jeffreys.Fit(Structures[0]), 0.0).Refusal);   // the intercept alone leaves a residual

        // The proper default keeps it, with a finite marginal.
        Built scaled = Build(exact, Config());
        BassEvaluation kept = scaled.Target.Evaluate(scaled.Fit(Structures[1]), Math.Log(0.02));
        Assert.Null(kept.Refusal);
        Assert.True(double.IsFinite(kept.LogMarginal));

        // The design's refusals carry through.
        Assert.Equal(BassRefusal.Numerical, scaled.Target.Evaluate(BassStructureFit.Refused(BassRefusal.Numerical), 0.0).Refusal);

        // Under the zero centre a constant response is fitted exactly by the intercept: refused under Jeffreys at entry.
        BassData constant = Data(constant: true);
        var ex = Assert.Throws<ArgumentException>(() => BassResolved.Resolve(Config(GPriorCentre.Zero, Sigma2Prior.Jeffreys), constant));
        Assert.Contains("fitted exactly by the intercept", ex.Message);
        BassResolved.Resolve(Config(GPriorCentre.Zero), constant);
    }

    [Theory]
    [MemberData(nameof(Centres))]
    public void Shift_MovesOnlyTheZeroCentre(GPriorCentre centre)
    {
        // F-O9: level centring and Liang read y − ȳ, so a shift of y moves no marginal; R's zero centre reads y.
        double[] offsets = [0.0, 10.0, -250.0];
        var marginals = offsets.Select(o =>
        {
            Built built = Build(Data(offset: o), Config(centre));
            return Structures.Select(s => built.Target.Evaluate(built.Fit(s), Math.Log(0.05)).LogMarginal).ToArray();
        }).ToArray();
        double worst = 0, moved = 0;
        for (int o = 1; o < offsets.Length; o++)
            for (int s = 0; s < Structures.Length; s++)
            {
                worst = Math.Max(worst, Math.Abs(marginals[o][s] - marginals[0][s]));
                moved = Math.Max(moved, Math.Abs(marginals[o][s] - marginals[o][0] - (marginals[0][s] - marginals[0][0])));
            }
        if (centre == GPriorCentre.Zero) Assert.True(moved > 1.0, $"the zero centre's Bayes factors moved by only {moved:F3}");
        else Assert.True(worst <= 1e-9, $"{centre}: a shift moved a marginal by {worst:E2}");
    }

    [Fact]
    public void Scale_MovesEveryStructureAlike()
    {
        // Under the scaled σ² prior, y × c moves every marginal by −N ln c (D03; OQ19's route).
        Built unit = Build(Data(), Config());
        double[] reference = Structures.Select(s => unit.Target.Evaluate(unit.Fit(s), Math.Log(0.05)).LogMarginal).ToArray();
        foreach (double c in new[] { 1e-6, 1e3, 1e6 })
        {
            Built built = Build(Data(scale: c), Config());
            for (int s = 0; s < Structures.Length; s++)
            {
                double shift = built.Target.Evaluate(built.Fit(Structures[s]), Math.Log(0.05)).LogMarginal - reference[s];
                double expected = -75 * Math.Log(c);
                Assert.True(Math.Abs(shift - expected) <= 1e-9 * Math.Abs(expected), $"c = {c}, structure {s}: {shift:R} against {expected:R}");
            }
        }
    }

    [Fact]
    public void DeclaredInverseGamma_ReachesTheTarget()
    {
        // A declared IG(g₁, g₂) is held as given, reading nothing from y, and the marginal is the block's under it.
        BassData data = Data();
        Built built = Build(data, Config(sigma2: Sigma2Prior.InverseGamma(2.0, 3.0)));
        InverseGamma sigma = Assert.NotNull(built.Target.Sigma2);
        Assert.Equal(2.0, sigma.Shape);
        Assert.Equal(3.0, sigma.Scale);
        BassStructureFit fit = built.Fit(Structures[2]);
        double tau = 0.04;
        double expected = CollapsedGPrior.CentredMarginalAt(fit.Terms, fit.ResidualSquares, fit.FittedSquares, 75, built.Resolved.ResponseMean,
                                                            GPriorCentre.Mean, tau).LogMarginal(InverseGamma.ShapeScale(2.0, 3.0));
        Assert.Equal(expected, built.Target.Evaluate(fit, Math.Log(tau)).LogMarginal);
        Assert.True(built.Target.Defines(0.0));

        // Unlike the scaled prior it is in y's units, so rescaling y moves structures unequally.
        Built scaled = Build(Data(scale: 1e3), Config(sigma2: Sigma2Prior.InverseGamma(2.0, 3.0)));
        double Gap(Built b) => b.Target.Evaluate(b.Fit(Structures[2]), Math.Log(tau)).LogMarginal - b.Target.Evaluate(b.Fit(Structures[0]), Math.Log(tau)).LogMarginal;
        Assert.True(Math.Abs(Gap(scaled) - Gap(built)) > 1.0);
    }

    [Fact]
    public void TauOutsideTheDoubles_IsOffTheSupport()
    {
        // e^u rounds to 0 below u ≈ −745 and overflows above ≈ 709: −∞ for the prior and the likelihood alike, as
        // ITargetEvaluation requires off the support, rather than a throw from the marginal.
        Built built = Build(Data(), Config());
        BassStructureFit fit = built.Fit(Structures[2]);
        foreach (double u in new[] { -800.0, 800.0, double.NegativeInfinity })
        {
            BassEvaluation ev = built.Target.Evaluate(fit, u);
            Assert.True(ev.OffSupport, $"u = {u}");
            Assert.Equal(double.NegativeInfinity, ev.LogPrior);
            Assert.Equal(double.NegativeInfinity, ev.LogLikelihood(1.0));
            Assert.Equal(double.NegativeInfinity, ev.LogLikelihood(0.5));
            Assert.Equal(double.NegativeInfinity, Temperature.LogPosterior(ev));
            Assert.Throws<InvalidOperationException>(() => ev.LogLikelihoodSlope(0.5));
        }
        Assert.False(built.Target.Evaluate(fit, -700.0).OffSupport);
    }

    // Composite 12-point Gauss–Legendre on panels of width at most `width` (½ by default) over [lo, hi].
    internal static IEnumerable<(double U, double W)> Nodes(double lo, double hi, double width = 0.5)
    {
        (double[] x, double[] w) = GaussLegendre(12);
        int panels = (int)Math.Ceiling((hi - lo) / width);
        double step = (hi - lo) / panels;
        for (int p = 0; p < panels; p++)
            for (int i = 0; i < x.Length; i++)
                yield return (lo + step * (p + 0.5 * (x[i] + 1)), 0.5 * step * w[i]);
    }

    // Nodes and weights of the n-point rule on [−1, 1], by Newton's method on Pₙ.
    internal static (double[] X, double[] W) GaussLegendre(int n)
    {
        var x = new double[n];
        var w = new double[n];
        for (int i = 0; i < n; i++)
        {
            double z = Math.Cos(Math.PI * (i + 0.75) / (n + 0.5));
            for (int iteration = 0; iteration < 100; iteration++)
            {
                (double p, double dp) = Legendre(n, z);
                double next = z - p / dp;
                bool done = Math.Abs(next - z) <= 1e-16;
                z = next;
                if (done) break;
            }
            x[i] = z;
            double slope = Legendre(n, z).Slope;
            w[i] = 2 / ((1 - z * z) * slope * slope);
        }
        return (x, w);
    }

    private static (double Value, double Slope) Legendre(int n, double z)
    {
        double p0 = 1, p1 = z;
        for (int k = 2; k <= n; k++) (p0, p1) = (p1, ((2 * k - 1) * z * p1 - (k - 1) * p0) / k);
        return (p1, n * (z * p1 - p0) / (z * z - 1));
    }
}
