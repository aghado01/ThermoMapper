using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The normalized ridge's marginal (D03's ridge, decided 2026-10-29): β | σ², τ ~ N(0, σ²(τD)⁻¹), D = diag(G_c)/N, with
/// the intercept a unit-RMS column inside it (or flat). Checked against a dense computation on the explicit design
/// [1, Z_c] at b = 1 and two tempered b under every intercept convention; for independence of the columns' scale; on
/// orthogonal columns, where it is the g-form at τ/N; and, through the design layer, each bordered candidate against a
/// fresh decomposition of the candidate over randomized birth, death and change sequences, near-collinear columns
/// included.
/// </summary>
public sealed class BassRidgeTests
{
    private readonly ITestOutputHelper _out;
    public BassRidgeTests(ITestOutputHelper output) => _out = output;

    internal static readonly BassTerm[] ThreeTerms =
        [.. BassCoefficientPriorTests.TwoTerms, new BassTerm([new HingeFactor(1, 1, 1)])];

    internal static BassModel Model(BassData? data = null, BassCoefficientPrior? coefficients = null, GPriorCentre centre = GPriorCentre.Mean,
                                    Sigma2Prior? sigma2 = null, int maxBasis = 10, double angleFloor = 1e-10)
    {
        BassConfig config = BassSweepTests.Config(maxBasis: maxBasis);
        config = config with
        {
            Structure = config.Structure with { AngleFloor = angleFloor },
            Coefficients = coefficients ?? BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()),
            Intercept = centre,
            Sigma2 = sigma2 ?? Sigma2Prior.ScaledToResponse(1e-3),
        };
        return BassModel.Create(config, data ?? BassSweepTests.Data(large: true));
    }

    // The explicit centred columns of `terms` on a scalar response.
    internal static double[][] CentredColumns(BassModel model, IReadOnlyList<BassTerm> terms)
    {
        int n = model.Data.Response.Cases;
        return [.. terms.Select(t =>
        {
            var a = new double[n];
            model.Basis.CaseColumn(t, model.Basis.Inputs.Cases, a);
            double mean = a.Average();
            return a.Select(v => v - mean).ToArray();
        })];
    }

    // The marginal's S, determinant term and ν on the explicit design, by a dense factorization: under level centring and
    // the zero centre X = [1, Z_c] with y as given, prior precision τ·diag(1, d) about (level, 0); under the flat
    // intercept Z_c with ỹ, prior precision τ·diag(d), and its −½ln(Nb). S = min_β b‖y − Xβ‖² + (β − μ₀)ᵀP₀(β − μ₀).
    internal static (double S, double LogDet, double Nu) Dense(BassModel model, IReadOnlyList<BassTerm> terms, double tau, double b)
    {
        double[] y = model.Data.Response.Values.ToArray();
        int n = y.Length;
        double mean = y.Average();
        GPriorCentre centre = model.Config.Intercept;
        bool flat = centre == GPriorCentre.Flat;
        double[][] z = CentredColumns(model, terms);
        var columns = new List<double[]>();
        var precision = new List<double>();
        if (!flat)
        {
            columns.Add([.. Enumerable.Repeat(1.0, n)]);
            precision.Add(tau);
        }
        foreach (double[] c in z)
        {
            columns.Add(c);
            precision.Add(tau * c.Sum(v => v * v) / n);
        }
        double[] response = flat ? [.. y.Select(v => v - mean)] : y;
        double level = centre == GPriorCentre.Mean ? mean : 0.0;
        int p = columns.Count;
        var mu = new double[p];
        if (!flat) mu[0] = level;
        var a = new double[p, p];
        var rhs = new double[p];
        for (int i = 0; i < p; i++)
        {
            for (int j = 0; j < p; j++) a[i, j] = b * columns[i].Zip(columns[j], (u, v) => u * v).Sum() + (i == j ? precision[i] : 0.0);
            rhs[i] = b * columns[i].Zip(response, (u, v) => u * v).Sum() + precision[i] * mu[i];
        }
        double logDet = 0.0, s;
        foreach (double d in precision) logDet += 0.5 * Math.Log(d);
        if (p > 0)
        {
            var factor = new DenseCholesky(p);
            factor.Decompose(a);
            double[] m = factor.Solve(rhs);
            logDet -= 0.5 * factor.LogDet;
            double residual = 0.0;
            for (int r = 0; r < n; r++)
            {
                double fit = 0.0;
                for (int j = 0; j < p; j++) fit += columns[j][r] * m[j];
                residual += (response[r] - fit) * (response[r] - fit);
            }
            double penalty = 0.0;
            for (int j = 0; j < p; j++) penalty += precision[j] * (m[j] - mu[j]) * (m[j] - mu[j]);
            s = b * residual + penalty;
        }
        else s = b * response.Sum(v => v * v);
        if (flat) logDet -= 0.5 * Math.Log(n * b);
        return (s, logDet, n * b - (flat ? 1 : 0));
    }

    public static TheoryData<GPriorCentre, string> Conventions => new()
    {
        { GPriorCentre.Mean, "scaled" }, { GPriorCentre.Zero, "scaled" }, { GPriorCentre.Flat, "scaled" },
        { GPriorCentre.Mean, "jeffreys" }, { GPriorCentre.Zero, "jeffreys" }, { GPriorCentre.Flat, "jeffreys" },
    };

    [Theory]
    [MemberData(nameof(Conventions))]
    public void Marginal_MatchesTheDenseComputation(GPriorCentre centre, string sigma2)
    {
        BassModel model = Model(centre: centre, sigma2: sigma2 == "jeffreys" ? Sigma2Prior.Jeffreys : null);
        foreach (BassTerm[] terms in new[] { Array.Empty<BassTerm>(), ThreeTerms[..1], ThreeTerms })
        {
            BassDesign design = model.NewDesign();
            Assert.Null(design.Load(terms));
            BassStructureFit fit = BassStructureFit.Of(design, 0.0);
            foreach (double tau in new[] { 1e-4, 0.7, 40.0 })
                foreach (double b in new[] { 1.0, 0.3, 0.02 })
                {
                    GPriorMarginal m = model.Target.Marginal(fit, tau, b);
                    (double s, double logDet, double nu) = Dense(model, terms, tau, b);
                    Assert.True(Math.Abs(m.S - s) <= 1e-11 * s, $"M = {terms.Length}, τ {tau}, b {b}: S {m.S:R} against {s:R}");
                    Assert.True(Math.Abs(m.LogDeterminantTerm - logDet) <= 1e-11 * Math.Max(1.0, Math.Abs(logDet)),
                                $"M = {terms.Length}, τ {tau}, b {b}: ln det {m.LogDeterminantTerm:R} against {logDet:R}");
                    Assert.Equal(nu, m.DegreesOfFreedom, 12);
                }
        }
    }

    [Fact]
    public void Marginal_DoesNotDependOnTheColumnsScale()
    {
        // The prior is stated on the columns at unit root mean square, so rescaling the raw columns moves neither the
        // spectrum's terms nor the dense marginal.
        BassModel model = Model();
        double[][] z = CentredColumns(model, ThreeTerms);
        double[] y = model.Data.Response.Values.ToArray();
        double mean = y.Average();
        int n = y.Length;
        RidgeSpectrum Spectrum(double[] scales)
        {
            int m = z.Length;
            var g = new double[m, m];
            var c = new double[m];
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < m; j++) g[i, j] = scales[i] * scales[j] * z[i].Zip(z[j], (u, v) => u * v).Sum();
                c[i] = scales[i] * z[i].Zip(y, (u, v) => u * (v - mean)).Sum();
            }
            return RidgeSpectrum.Of(g, c, n);
        }
        RidgeSpectrum unit = Spectrum([1, 1, 1]), scaled = Spectrum([1e-3, 7, 1e4]);
        foreach (double tau in new[] { 1e-3, 1.0, 30.0 })
            foreach (double b in new[] { 1.0, 0.1 })
            {
                RidgeTerms u = unit.At(tau, b), s = scaled.At(tau, b);
                Assert.Equal(u.LogDeterminant, s.LogDeterminant, 12);
                Assert.True(Math.Abs(u.Shrinkage - s.Shrinkage) <= 1e-11 * u.Shrinkage, $"{u.Shrinkage:R} against {s.Shrinkage:R}");
                Assert.True(Math.Abs(u.ShrinkageSlope - s.ShrinkageSlope) <= 1e-11 * u.ShrinkageSlope);
                Assert.Equal(u.LogDeterminantSlope, s.LogDeterminantSlope, 10);
            }
        Assert.Equal(unit.Eigenvalues.ToArray().Sum(), 3.0 * n, 9);   // tr G̃ = N·M
    }

    [Fact]
    public void OnOrthogonalColumns_TheRidgeIsTheGFormAtTauOverN()
    {
        // Orthogonal centred columns make G̃ = N·I, so τD = (τ/N)G_c: the ridge at τ is the g-form at τ/N, column for
        // column and the intercept with them.
        const int n = 8;
        double[][] h = [[1, -1, 1, -1, 1, -1, 1, -1], [1, 1, -1, -1, 1, 1, -1, -1], [1, -1, -1, 1, 1, -1, -1, 1]];
        double[] scale = [0.3, 2.0, 11.0];
        double[] y = [3.1, 0.4, -1.2, 2.2, 0.9, 1.7, -0.3, 0.8];
        double mean = y.Average();
        var g = new double[3, 3];
        var c = new double[3];
        double fitted = 0.0, total = y.Sum(v => (v - mean) * (v - mean));
        for (int i = 0; i < 3; i++)
        {
            g[i, i] = scale[i] * scale[i] * n;
            c[i] = scale[i] * h[i].Zip(y, (u, v) => u * (v - mean)).Sum();
            fitted += c[i] * c[i] / g[i, i];
        }
        RidgeSpectrum spectrum = RidgeSpectrum.Of(g, c, n);
        foreach (GPriorCentre centre in new[] { GPriorCentre.Mean, GPriorCentre.Zero, GPriorCentre.Flat })
            foreach (double tau in new[] { 0.05, 1.0, 9.0 })
                foreach (double b in new[] { 1.0, 0.25 })
                {
                    GPriorMarginal ridge = CollapsedGPrior.CentredRidgeMarginalAt(spectrum, total - fitted, n, mean, centre, tau, b);
                    GPriorMarginal gForm = CollapsedGPrior.CentredMarginalAt(3, total - fitted, fitted, n, mean, centre, tau / n, b);
                    Assert.True(Math.Abs(ridge.S - gForm.S) <= 1e-12 * gForm.S, $"{centre}, τ {tau}, b {b}: S {ridge.S:R} against {gForm.S:R}");
                    Assert.Equal(gForm.LogDeterminantTerm, ridge.LogDeterminantTerm, 12);
                    Assert.Equal(gForm.DegreesOfFreedom, ridge.DegreesOfFreedom);
                    Assert.Equal(gForm.LogMarginal(null), ridge.LogMarginal(null), 11);
                    Assert.Equal(gForm.LogMarginalSlope(null), ridge.LogMarginalSlope(null), 10);
                }
    }

    // A fixture with many close knots, whose adjacent hinges are nearly collinear: x₀ on 160 distinct values, x₁ on 40.
    internal static BassData DenseKnots(int seed = 3)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        const int n = 160;
        var x0 = new double[n];
        var x1 = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x0[i] = i / (n - 1.0);
            x1[i] = (i * 37 % 40) / 39.0;
            y[i] = Math.Sin(5 * x0[i]) + x1[i] * x0[i] + 0.05 * (rng.NextDouble() - 0.5);
        }
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1)], new ScalarResponse(y));
    }

    [Theory]
    [InlineData(1e-10)]
    [InlineData(2e-3)]
    public void Admission_IsTheGForms(double floor)
    {
        // D06's gates read the Gram alone, so a ridge design admits and refuses exactly as a g-form one: the same proposals
        // give the same accounts, bit for bit, the ridge's terms beside them. F is unchanged.
        BassData data = DenseKnots();
        BassModel g = Model(data, BassCoefficientPrior.G(), maxBasis: 12, angleFloor: floor), ridge = Model(data, maxBasis: 12, angleFloor: floor);
        BassDesign a = g.NewDesign(), b = ridge.NewDesign();
        Assert.False(a.IsRidge);
        Assert.True(b.IsRidge);
        Assert.Throws<InvalidOperationException>(() => a.Spectrum);
        int[] knots = [g.Inputs.ContinuousKnots(0).Length, g.Inputs.ContinuousKnots(1).Length];
        var rng = new Xoshiro256PlusPlus(23);
        int refused = 0;
        for (int step = 0; step < 500; step++)
        {
            var term = new BassTerm([new HingeFactor(0, rng.NextDouble() < 0.5 ? 1 : -1, BassDraws.Below(rng, knots[0]))]);
            int count = a.Count;
            double u = rng.NextDouble();
            int index = count == 0 ? -1 : BassDraws.Below(rng, count);
            (DesignProposal pa, DesignProposal pb) = count == 0 || u < 0.5 ? (a.ProposeBirth(term), b.ProposeBirth(term))
                : u < 0.75 ? (a.ProposeDeath(index), b.ProposeDeath(index)) : (a.ProposeChange(index, term), b.ProposeChange(index, term));
            Assert.Equal(pa with { Ridge = null }, pb with { Ridge = null });
            Assert.Equal(pb.Admitted, pb.Ridge is not null);
            if (!pa.Admitted) { refused++; continue; }
            if (rng.NextDouble() < 0.6) { a.Accept(); b.Accept(); }
        }
        Assert.Equal(a.Terms, b.Terms);
        _out.WriteLine($"s_min {floor}: {refused} refusals, shared by both forms");
        Assert.True(refused > 0);
    }

    public static TheoryData<string, double, int> MoveCases => new()
    {
        { "grid", 1e-10, 11 }, { "close knots", 1e-10, 12 }, { "close knots, strict", 2e-3, 13 },
    };

    [Theory]
    [MemberData(nameof(MoveCases))]
    public void BorderedCandidates_MatchFreshDecompositions(string fixture, double floor, int seed)
    {
        BassData data = fixture == "grid" ? BassSweepTests.Data(large: true) : DenseKnots();
        BassModel model = Model(data, maxBasis: 12, angleFloor: floor);
        int[] knots = [model.Inputs.ContinuousKnots(0).Length, model.Inputs.ContinuousKnots(1).Length];
        var rng = new Xoshiro256PlusPlus(seed);
        BassTerm Draw()
        {
            int order = rng.NextDouble() < 0.6 ? 1 : 2;
            var factors = new List<HingeFactor>();
            foreach (int v in order == 1 ? new[] { BassDraws.Below(rng, 2) } : new[] { 0, 1 })
                factors.Add(new HingeFactor(v, rng.NextDouble() < 0.5 ? 1 : -1, BassDraws.Below(rng, knots[v])));
            return new BassTerm(factors);
        }
        BassDesign design = model.NewDesign();
        int bordered = 0, fresh = 0;
        double worst = 0.0;
        (double Tau, double B)[] points = [(1e-3, 1.0), (0.5, 1.0), (20.0, 1.0), (0.5, 0.05), (3.0, 0.0)];
        for (int step = 0; step < 600; step++)
        {
            List<BassTerm> current = [.. design.Terms];
            double u = rng.NextDouble();
            DesignProposal proposal;
            List<BassTerm> candidate;
            if (current.Count == 0 || (u < 0.45 && current.Count < 12))
            {
                BassTerm term = Draw();
                proposal = design.ProposeBirth(term);
                candidate = [.. current, term];
            }
            else if (u < 0.7)
            {
                int index = BassDraws.Below(rng, current.Count);
                proposal = design.ProposeDeath(index);
                candidate = [.. current.Where((_, i) => i != index)];
            }
            else
            {
                int index = BassDraws.Below(rng, current.Count);
                BassTerm term = Draw();
                proposal = design.ProposeChange(index, term);
                candidate = term.Equals(current[index]) ? current : [.. current.Where((_, i) => i != index), term];
            }
            if (!proposal.Admitted) continue;
            RidgeColumns view = proposal.Ridge!;
            if (view is RidgeSpectrum) fresh++;
            else bordered++;
            BassDesign reference = model.NewDesign();
            Assert.Null(reference.Load(candidate));
            RidgeSpectrum exact = reference.Spectrum;
            Assert.Equal(exact.Count, view.Count);
            foreach ((double tau, double b) in points)
            {
                RidgeTerms got = view.At(tau, b), want = exact.At(tau, b);
                double scale = Math.Max(want.Shrinkage, 1e-300);
                double error = Math.Max(Math.Max(Math.Abs(got.LogDeterminant - want.LogDeterminant) / Math.Max(1.0, Math.Abs(want.LogDeterminant)),
                                                 Math.Abs(got.Shrinkage - want.Shrinkage) / scale),
                                        Math.Max(Math.Abs(got.LogDeterminantSlope - want.LogDeterminantSlope) / Math.Max(1.0, Math.Abs(want.LogDeterminantSlope)),
                                                 Math.Abs(got.ShrinkageSlope - want.ShrinkageSlope) / Math.Max(want.ShrinkageSlope, 1e-300)));
                worst = Math.Max(worst, error);
                Assert.True(error <= 1e-8, $"step {step}, M {candidate.Count}, τ {tau}, b {b}: {got} against {want}");
            }
            if (rng.NextDouble() < 0.5) design.Accept();
        }
        _out.WriteLine($"{fixture}, s_min {floor}: {bordered} bordered and {fresh} fresh candidates, worst relative error {worst:G3}");
        Assert.True(bordered > 50);
        if (floor > 1e-6) Assert.True(fresh > 0, "the strict path was not exercised");
    }
}
