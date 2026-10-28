using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// BASS's design layer and admission (stage 4, S4). Every proposal and every accepted state is checked against a
/// fresh reference that shares nothing with the bordered factor: the explicitly centred design, its least-squares fit
/// and its leave-one-out angles s_j = 1/(G_jj (G⁻¹)_jj), all from a thin SVD. Random move sequences on scalar and
/// direct functional data, at the default angle floor and at one that binds; F-O2's order-freedom with hinges at
/// neighbouring knots; one constructed case per refusal reason; wholesale loads; the near-exact-fit path; and the
/// residual and factor after a long run.
/// </summary>
public sealed class BassDesignTests
{
    private readonly ITestOutputHelper _out;
    public BassDesignTests(ITestOutputHelper output) => _out = output;

    private sealed record Setup(BassData Data, BassBasis Basis);

    private static Setup Scalar(int n, int seed, Func<double, double, double>? surface = null)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x1 = Enumerable.Range(0, n).Select(_ => Math.Round(rng.NextDouble() * 40) / 40).ToArray();
        double[] x2 = Enumerable.Range(0, n).Select(_ => rng.NextDouble()).ToArray();
        int[] z = Enumerable.Range(0, n).Select(i => i % 4).ToArray();
        surface ??= (a, b) => Math.Sin(3 * a) + a * b;
        double[] y = Enumerable.Range(0, n).Select(i => surface(x1[i], x2[i]) + 0.1 * (rng.NextDouble() - 0.5) + (z[i] == 2 ? 0.3 : 0)).ToArray();
        var data = new BassData([new ContinuousPredictor("x1", x1), new ContinuousPredictor("x2", x2),
                                 new CategoricalPredictor("z", z, ["a", "b", "c", "d"])], new ScalarResponse(y));
        return new Setup(data, new BassBasis(BassInputs.Prepare(data, BassInputTransform.Fit(data)), 1));
    }

    private static Setup Functional(int n, int m, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        double[] x = Enumerable.Range(0, n).Select(_ => rng.NextDouble()).ToArray();
        int[] z = Enumerable.Range(0, n).Select(i => i % 3).ToArray();
        double[] t = Enumerable.Range(0, m).Select(r => r / (m - 1.0)).ToArray();
        var y = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int r = 0; r < m; r++) y[i, r] = Math.Sin(4 * x[i] * t[r]) + (z[i] == 1 ? t[r] : 0) + 0.05 * rng.NextDouble();
        var data = new BassData([new ContinuousPredictor("x", x), new CategoricalPredictor("z", z, ["a", "b", "c"])],
                                new FunctionalResponse(y, [new ContinuousPredictor("t", t)]));
        return new Setup(data, new BassBasis(BassInputs.Prepare(data, BassInputTransform.Fit(data)), 1));
    }

    private static BassTerm RandomTerm(BassBasis basis, Xoshiro256PlusPlus rng)
    {
        BassInputs inputs = basis.Inputs;
        while (true)
        {
            var cont = new List<HingeFactor>();
            for (int v = 0; v < inputs.Cases.ContinuousCount; v++)
                if (rng.NextDouble() < 0.45)
                    cont.Add(new HingeFactor(v, rng.NextDouble() < 0.5 ? -1 : 1, (int)(rng.NextDouble() * inputs.ContinuousKnots(v).Length)));
            var cat = new List<CategoricalFactor>();
            for (int v = 0; v < inputs.Cases.CategoricalCount; v++)
                if (rng.NextDouble() < 0.3)
                {
                    int levels = inputs.LevelCount(v);
                    int[] subset = [.. Enumerable.Range(0, levels).Where(_ => rng.NextDouble() < 0.5)];
                    if (subset.Length > 0 && subset.Length < levels) cat.Add(new CategoricalFactor(v, subset));
                }
            var func = new List<HingeFactor>();
            if (inputs.Outputs is { } outputs && rng.NextDouble() < 0.5)
                func.Add(new HingeFactor(0, rng.NextDouble() < 0.5 ? -1 : 1, (int)(rng.NextDouble() * inputs.FunctionalKnots(0).Length)));
            if (cont.Count + cat.Count + func.Count > 0) return new BassTerm(cont, cat, func);
        }
    }

    // The reference: the explicitly centred design, its least-squares fit and its leave-one-out angles, by thin SVD.
    private sealed record Reference(BassRefusal? Refusal, double ResidualSquares, double FittedSquares, double MinAngle, double[] Residual);

    private static Reference Fresh(Setup s, IReadOnlyList<BassTerm> terms, BassAdmission rules)
    {
        BassBasis basis = s.Basis;
        int n = s.Data.Response.Cases, m = s.Data.Response.Locations, rows = n * m;
        double[] y = s.Data.Response switch
        {
            ScalarResponse sr => sr.Y.ToArray(),
            FunctionalResponse fr => [.. Enumerable.Range(0, rows).Select(k => fr[k / m, k % m])],
            _ => throw new InvalidOperationException(),
        };
        double mean = y.Average();
        double[] yc = [.. y.Select(v => v - mean)];
        if (terms.Count > rules.MaxBasis) return new(BassRefusal.Cap, double.NaN, double.NaN, double.NaN, []);
        if (terms.Distinct().Count() < terms.Count) return new(BassRefusal.Numerical, double.NaN, double.NaN, double.NaN, []);
        var x = new double[rows, terms.Count];
        for (int j = 0; j < terms.Count; j++)
        {
            var a = new double[n];
            basis.CaseColumn(terms[j], basis.Inputs.Cases, a);
            var h = Enumerable.Repeat(1.0, m).ToArray();
            if (basis.Inputs.Outputs is { } outputs) basis.OutputColumn(terms[j], outputs, h);
            if (a.Count(v => v != 0) < rules.SupportMinimum) return new(BassRefusal.Support, double.NaN, double.NaN, double.NaN, []);
            if (rules.SupportMinimumFunctional is { } f && h.Count(v => v != 0) < f) return new(BassRefusal.Support, double.NaN, double.NaN, double.NaN, []);
            for (int k = 0; k < rows; k++) x[k, j] = a[k / m] * h[k % m];
            double cm = Enumerable.Range(0, rows).Average(k => x[k, j]);
            for (int k = 0; k < rows; k++) x[k, j] -= cm;
        }
        if (terms.Count == 0) return new(null, yc.Sum(v => v * v), 0.0, double.PositiveInfinity, yc);
        ThinSvdResult svd = ThinSvd.Decompose(x);
        double minAngle = double.PositiveInfinity;
        for (int j = 0; j < terms.Count; j++)
        {
            double gjj = Enumerable.Range(0, rows).Sum(k => x[k, j] * x[k, j]), inverse = 0.0;
            for (int c = 0; c < terms.Count; c++) inverse += svd.RightVectors[c][j] * svd.RightVectors[c][j] / (svd.SingularValues[c] * svd.SingularValues[c]);
            minAngle = Math.Min(minAngle, 1.0 / (gjj * inverse));
        }
        if (svd.NumericalRank < terms.Count || minAngle < rules.AngleFloor) return new(BassRefusal.Numerical, double.NaN, double.NaN, minAngle, []);
        if (rules.MaxVarianceInflation is { } cap && minAngle * cap < 1) return new(BassRefusal.VarianceInflation, double.NaN, double.NaN, minAngle, []);
        double fitted = 0.0;
        var residual = (double[])yc.Clone();
        for (int c = 0; c < terms.Count; c++)
        {
            double coefficient = Enumerable.Range(0, rows).Sum(k => svd.LeftVectors[c][k] * yc[k]);
            fitted += coefficient * coefficient;
            for (int k = 0; k < rows; k++) residual[k] -= coefficient * svd.LeftVectors[c][k];
        }
        return new(null, residual.Sum(v => v * v), fitted, minAngle, residual);
    }

    private static void Close(double expected, double actual, double relative, double scale, string what)
        => Assert.True(Math.Abs(expected - actual) <= relative * Math.Max(Math.Abs(expected), scale), $"{what}: {actual:R} against {expected:R}");

    [Theory]
    [InlineData(false, 1e-10)]
    [InlineData(false, 3e-3)]
    [InlineData(true, 1e-10)]
    [InlineData(true, 3e-3)]
    public void RandomMoves_MatchAFreshReference(bool functional, double floor)
    {
        Setup s = functional ? Functional(14, 9, seed: 5) : Scalar(70, seed: 3);
        var rules = new BassAdmission(8, SupportKind.Count, 3, functional ? 2 : null, floor, null);
        var design = new BassDesign(s.Basis, s.Data, rules);
        var rng = new Xoshiro256PlusPlus(seed: functional ? 101 : 202);
        double scale = design.ResponseSquares;
        int compared = 0, refused = 0, accepted = 0, nearAdmitted = 0;
        var reasons = Enum.GetValues<BassRefusal>().ToDictionary(r => r, _ => 0);
        for (int step = 0; step < 400; step++)
        {
            List<BassTerm> current = [.. design.Terms];
            double u = rng.NextDouble();
            List<BassTerm> candidate;
            DesignProposal p;
            if (current.Count == 0 || u < 0.4)
            {
                BassTerm t = RandomTerm(s.Basis, rng);
                candidate = [.. current, t];
                p = design.ProposeBirth(t);
            }
            else if (u < 0.65)
            {
                int index = (int)(rng.NextDouble() * current.Count);
                candidate = [.. current.Where((_, j) => j != index)];
                p = design.ProposeDeath(index);
            }
            else
            {
                int index = (int)(rng.NextDouble() * current.Count);
                BassTerm t = RandomTerm(s.Basis, rng);
                candidate = [.. current.Where((_, j) => j != index), t];
                if (t.Equals(current[index])) candidate = current;
                p = design.ProposeChange(index, t);
            }

            Reference r = Fresh(s, candidate, rules);
            bool nearFloor = r.MinAngle > floor / 1.05 && r.MinAngle < floor * 1.05;
            if (!nearFloor) Assert.Equal(r.Refusal, p.Refusal);
            if (p.Admitted && p.MinAngle < 2 * floor) nearAdmitted++;
            if (!p.Admitted) { refused++; reasons[p.Refusal!.Value]++; continue; }
            if (r.Refusal is null)
            {
                compared++;
                Close(r.ResidualSquares, p.ResidualSquares, 1e-8, 1e-12 * scale, $"step {step}: proposed ‖r‖²");
                Close(r.FittedSquares, p.FittedSquares, 1e-8, 1e-12 * scale, $"step {step}: proposed Q");
                if (double.IsFinite(r.MinAngle) && r.MinAngle > 1e-8) Close(r.MinAngle, p.MinAngle, 1e-6, 0, $"step {step}: proposed minimum angle");
                else if (!double.IsFinite(r.MinAngle)) Assert.Equal(r.MinAngle, p.MinAngle);   // no columns: +∞
            }
            if (rng.NextDouble() < 0.6)
            {
                design.Accept();
                accepted++;
                Assert.Equal(candidate.ToHashSet(), design.Terms.ToHashSet());
                Reference now = Fresh(s, design.Terms, rules with { AngleFloor = 0, MaxBasis = int.MaxValue });
                Close(now.ResidualSquares, design.ResidualSquares, 1e-9, 1e-12 * scale, $"step {step}: ‖r‖²");
                Close(now.FittedSquares, design.FittedSquares, 1e-9, 1e-12 * scale, $"step {step}: Q");
                for (int k = 0; k < now.Residual.Length; k++) Close(now.Residual[k], design.Residual[k], 1e-8, 1e-8 * Math.Sqrt(scale), $"step {step}: r[{k}]");
                Close(design.ResponseSquares, design.ResidualSquares + design.FittedSquares, 1e-10, 0, $"step {step}: ‖r‖² + Q = ‖ỹ‖²");
            }
        }
        _out.WriteLine($"[design] functional={functional} floor={floor:G}: {compared} compared, {accepted} accepted, {nearAdmitted} admitted within 2× of the floor, " +
                       $"refused {string.Join(", ", reasons.Select(kv => $"{kv.Key} {kv.Value}"))}; final M = {design.Count}");
        Assert.True(compared > 100 && accepted > 40 && refused > 10);
        if (floor > 1e-6) Assert.True(reasons[BassRefusal.Numerical] > 5, "a binding floor must refuse on the angle gate");
    }

    [Fact]
    public void TheGate_IsOrderFree_WhereANewestColumnRuleIsNot()
    {
        // Hinges at three neighbouring knots: the middle one is nearly the mean of the other two, so its angle to their
        // span is the smallest, and a newest-column rule would admit the set whenever the middle one is not last.
        double[] x = [.. Enumerable.Range(0, 200).Select(i => i / 199.0)];
        var data = new BassData([new ContinuousPredictor("x", x)], new ScalarResponse([.. x.Select(v => Math.Sin(5 * v))]));
        var s = new Setup(data, new BassBasis(BassInputs.Prepare(data, BassInputTransform.Fit(data)), 1));
        BassTerm h1 = new([new HingeFactor(0, 1, 100)]), h2 = new([new HingeFactor(0, 1, 104)]), h3 = new([new HingeFactor(0, 1, 108)]);

        var loose = new BassAdmission(10, SupportKind.Count, 1, null, 1e-14, null);
        double middle = Fresh(s, [h1, h2, h3], loose).MinAngle;
        var probe = new BassDesign(s.Basis, data, loose);
        probe.ProposeBirth(h1); probe.Accept();
        double pair = probe.ProposeBirth(h2).NewestAngle; probe.Accept();
        double newest = probe.ProposeBirth(h3).NewestAngle;
        double floor = Math.Sqrt(middle * newest);
        _out.WriteLine($"[order] min angle {middle:E3}, h2 after h1 {pair:E3}, h3 last {newest:E3}, floor {floor:E3}");
        Assert.True(middle < floor && floor < newest && floor < pair, "the fixture must separate the two rules");

        var rules = loose with { AngleFloor = floor };
        foreach (BassTerm[] order in new[] { new[] { h1, h2, h3 }, [h1, h3, h2], [h2, h1, h3], [h2, h3, h1], [h3, h1, h2], [h3, h2, h1] })
        {
            var design = new BassDesign(s.Basis, data, rules);
            Assert.True(design.ProposeBirth(order[0]).Admitted);
            design.Accept();
            Assert.True(design.ProposeBirth(order[1]).Admitted);
            design.Accept();
            Assert.Equal(BassRefusal.Numerical, design.ProposeBirth(order[2]).Refusal);
            Assert.Equal(BassRefusal.Numerical, design.Load([h1, h2, h3]));
        }
    }

    [Fact]
    public void EachRefusalReason_HasItsCase()
    {
        Setup s = Scalar(80, seed: 9);
        BassInputs inputs = s.Basis.Inputs;
        int top = inputs.ContinuousKnots(0).Length - 2;
        BassTerm a = new([new HingeFactor(0, 1, 5)]), b = new([new HingeFactor(1, -1, 20)]);
        var rules = new BassAdmission(2, SupportKind.Count, 10, null, 1e-10, null);

        var design = new BassDesign(s.Basis, s.Data, rules);
        Assert.True(design.ProposeBirth(a).Admitted); design.Accept();
        Assert.Equal(BassRefusal.Numerical, design.ProposeBirth(a).Refusal);                                   // a repeated term
        Assert.Equal(BassRefusal.Support, design.ProposeBirth(new BassTerm([new HingeFactor(0, 1, top)])).Refusal);
        Assert.True(design.ProposeBirth(b).Admitted); design.Accept();
        Assert.Equal(BassRefusal.Cap, design.ProposeBirth(new BassTerm([new HingeFactor(1, 1, 3)])).Refusal);
        Assert.Equal(BassRefusal.Numerical, design.ProposeChange(0, b).Refusal);                              // a change onto a term already there

        // A categorical subset and its complement: with the intercept they span one space, an exact collinearity.
        var factors = new BassDesign(s.Basis, s.Data, rules with { SupportMinimum = 1 });
        Assert.True(factors.ProposeBirth(new BassTerm(categorical: [new CategoricalFactor(0, [0, 1])])).Admitted); factors.Accept();
        Assert.Equal(BassRefusal.Numerical, factors.ProposeBirth(new BassTerm(categorical: [new CategoricalFactor(0, [2, 3])])).Refusal);

        // Two hinges at nearby knots on one variable, correlated past a variance inflation cap of 2.
        var capped = new BassDesign(s.Basis, s.Data, rules with { MaxVarianceInflation = 2 });
        Assert.True(capped.ProposeBirth(new BassTerm([new HingeFactor(1, 1, 10)])).Admitted); capped.Accept();
        Assert.Equal(BassRefusal.VarianceInflation, capped.ProposeBirth(new BassTerm([new HingeFactor(1, 1, 14)])).Refusal);

        // The functional part has its own threshold.
        Setup f = Functional(10, 12, seed: 4);
        var functional = new BassDesign(f.Basis, f.Data, new BassAdmission(5, SupportKind.Count, 1, 5, 1e-10, null));
        Assert.Equal(BassRefusal.Support, functional.ProposeBirth(new BassTerm(functional: [new HingeFactor(0, 1, 9)])).Refusal);
        Assert.True(functional.ProposeBirth(new BassTerm(functional: [new HingeFactor(0, 1, 3)])).Admitted);
    }

    [Fact]
    public void ColumnMass_IsTheOptionalSupportMeasure()
    {
        Setup s = Scalar(80, seed: 9);
        int knots = s.Basis.Inputs.ContinuousKnots(0).Length;
        BassTerm late = new([new HingeFactor(0, 1, knots - 6)]);   // nonzero values each below one, so its mass is below its count
        TermSupport support = s.Basis.Support(late);
        double threshold = support.CaseMass + 0.5 * (support.CaseNonzero - support.CaseMass);
        Assert.True(support.CaseNonzero - support.CaseMass > 1, $"fixture: {support}");
        var count = new BassDesign(s.Basis, s.Data, new BassAdmission(5, SupportKind.Count, threshold, null, 1e-10, null));
        var mass = new BassDesign(s.Basis, s.Data, new BassAdmission(5, SupportKind.ColumnMass, threshold, null, 1e-10, null));
        Assert.True(count.ProposeBirth(late).Admitted);
        Assert.Equal(BassRefusal.Support, mass.ProposeBirth(late).Refusal);
    }

    [Fact]
    public void Load_MatchesTheIncrementalPath_AndARefusalLeavesTheDesign()
    {
        Setup s = Scalar(70, seed: 3);
        var rules = new BassAdmission(10, SupportKind.Count, 3, null, 1e-10, null);
        var incremental = new BassDesign(s.Basis, s.Data, rules);
        var rng = new Xoshiro256PlusPlus(seed: 77);
        while (incremental.Count < 6)
            if (incremental.ProposeBirth(RandomTerm(s.Basis, rng)).Admitted) incremental.Accept();
        var loaded = new BassDesign(s.Basis, s.Data, rules);
        Assert.Null(loaded.Load([.. incremental.Terms.Reverse()]));
        Close(incremental.ResidualSquares, loaded.ResidualSquares, 1e-11, 0, "‖r‖²");
        Close(incremental.FittedSquares, loaded.FittedSquares, 1e-11, 0, "Q");

        double before = loaded.ResidualSquares;
        Assert.Equal(BassRefusal.Numerical, loaded.Load([incremental.Terms[0], incremental.Terms[0]]));
        Assert.Equal(6, loaded.Count);
        Assert.Equal(before, loaded.ResidualSquares);
    }

    [Fact]
    public void ANearExactFit_FormsTheCandidatesResidualExplicitly()
    {
        // y is a hinge plus a level, so the birth of that hinge leaves no residual at all.
        double[] x = [.. Enumerable.Range(0, 50).Select(i => i / 49.0)];
        var data = new BassData([new ContinuousPredictor("x", x)], new ScalarResponse([.. x.Select(v => 3 * Math.Max(v - x[20], 0) / (1 - x[20]) + 1)]));
        var s = new Setup(data, new BassBasis(BassInputs.Prepare(data, BassInputTransform.Fit(data)), 1));
        var design = new BassDesign(s.Basis, data, new BassAdmission(5, SupportKind.Count, 1, null, 1e-10, null));
        DesignProposal p = design.ProposeBirth(new BassTerm([new HingeFactor(0, 1, 20)]));
        Assert.True(p.Admitted);
        Assert.True(p.ResidualSquares <= 1e-24 * design.ResponseSquares, $"‖r′‖² = {p.ResidualSquares:R} of ‖ỹ‖² = {design.ResponseSquares:R}");
        design.Accept();
        Assert.True(design.ResidualSquares <= 1e-24 * design.ResponseSquares);
    }

    [Fact]
    public void ALongRun_KeepsTheResidualAndTheFactor()
    {
        Setup s = Functional(16, 11, seed: 8);
        var rules = new BassAdmission(12, SupportKind.Count, 2, 2, 1e-10, null);
        var design = new BassDesign(s.Basis, s.Data, rules);
        var rng = new Xoshiro256PlusPlus(seed: 303);
        int accepted = 0;
        while (accepted < 1500)
        {
            double u = rng.NextDouble();
            DesignProposal p = design.Count == 0 || u < 0.5 ? design.ProposeBirth(RandomTerm(s.Basis, rng))
                : u < 0.6 ? design.ProposeDeath((int)(rng.NextDouble() * design.Count))
                : design.ProposeChange((int)(rng.NextDouble() * design.Count), RandomTerm(s.Basis, rng));
            if (p.Admitted) { design.Accept(); accepted++; }
        }
        Assert.True(design.Count >= 8, $"the run should end near its cap of 12, at M = {design.Count}");
        Reference r = Fresh(s, design.Terms, rules);
        Close(r.ResidualSquares, design.ResidualSquares, 1e-9, 1e-12 * design.ResponseSquares, "‖r‖² after 1,500 moves");
        FactorDrift drift = design.Rebuild();
        _out.WriteLine($"[drift] M = {design.Count}: factor {drift.Factor:E2}, inverse diagonal {drift.InverseDiagonal:E2}");
        Assert.True(drift.PositiveDefinite && drift.Factor < 1e-10 && drift.InverseDiagonal < 1e-8);
    }
}
