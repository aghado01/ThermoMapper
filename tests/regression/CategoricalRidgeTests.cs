using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E36's slice 5, categorical inputs (decision 7, approved with 5a, 5b and 5c at their recommendations): BASS's subset factors gating a
/// direction or a ridge block, Π_v 1{c_v ∈ S_v}·b(θᵀz | t), the knot's range and a block's knots read over the term's cells; terms of
/// categorical factors alone; the determinant reading directions alone, its normalizer and cap at M_θ.
/// </summary>
public sealed class CategoricalRidgeTests(ITestOutputHelper output)
{
    private static readonly double[] Lattice = [0.6180339887498949, 0.4142135623730950, 0.7320508075688772];

    // Three continuous inputs on [0, 1] on golden-ratio lattices and a categorical input c of three levels, case i at level i mod 3; y follows
    // exp(0.6x₀ + 0.8x₁) where c ≠ 1 and 0.5 where c = 1.
    internal static BassData Data(int n = 160, int seed = 7)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var x = new double[3][];
        for (int j = 0; j < 3; j++) x[j] = [.. Enumerable.Range(0, n).Select(i => (0.5 + i * Lattice[j]) % 1.0)];
        int[] c = [.. Enumerable.Range(0, n).Select(i => i % 3)];
        double[] y = [.. Enumerable.Range(0, n).Select(i => (c[i] != 1 ? Math.Exp(0.6 * x[0][i] + 0.8 * x[1][i]) : 0.5) + 0.05 * Normal.SampleStandard(rng))];
        return new BassData([.. Enumerable.Range(0, 3).Select(j => (BassPredictor)new ContinuousPredictor($"x{j}", x[j], new InputDomain(0.0, 1.0))),
                             new CategoricalPredictor("c", c, ["a", "b", "c"])], new ScalarResponse(y));
    }

    internal static BassConfig Config(int columns = 1, DirectionPrior directions = DirectionPrior.Independent, int maxBasis = 8, int maxCategorical = 1)
    {
        BassConfig config = DirectionalChainTests.Config(maxBasis);
        return config with
        {
            Structure = config.Structure with
            {
                MaxInteraction = 2, MaxCategoricalInteraction = maxCategorical, Directional = new DirectionalBasis { Columns = columns, Directions = directions },
            },
        };
    }

    // The data's continuous predictors alone, which slice 1's and 2's hand columns read, index for index.
    internal static BassData Continuous(BassData data) => new([.. data.Predictors.OfType<ContinuousPredictor>()], data.Response);

    internal static int[][] Levels(BassData data) => [.. data.Predictors.OfType<CategoricalPredictor>().Select(p => p.Codes.ToArray())];

    // The cells of `cells`: 1 where every factor holds the case's level.
    internal static double[] Mask(int[][] levels, IReadOnlyList<CategoricalFactor> cells)
        => [.. Enumerable.Range(0, levels.Length == 0 ? 0 : levels[0].Length).Select(i => cells.All(c => c.Contains(levels[c.Variable][i])) ? 1.0 : 0.0)];

    // R's type-7 quantile of sorted values.
    private static double Quantile(double[] sorted, double probability)
    {
        double index = (sorted.Length - 1) * probability;
        int lo = (int)Math.Floor(index);
        double h = index - lo;
        return lo + 1 < sorted.Length ? sorted[lo] + h * (sorted[lo + 1] - sorted[lo]) : sorted[lo];
    }

    // The term's column by hand: its cells' mask times slice 1's hand hinge, or the mask alone for categorical factors alone.
    internal static double[] HandColumn(BassData data, BassTerm term, int degree = 1)
    {
        double[] mask = Mask(Levels(data), term.Categorical);
        if (term.Directional is not { } h) return mask;
        BassData continuous = Continuous(data);
        double[] hinge = DirectionalHingeTests.Hand(continuous, DirectionalHingeTests.Raw(continuous), h.Direction, h.Knot, degree, HingeNormalization.Box).Values;
        return [.. hinge.Select((v, i) => v * mask[i])];
    }

    private static double Dense(BassModel model, BassData data, IEnumerable<BassTerm> terms)
    {
        double[] y = ((ScalarResponse)data.Response).Values.ToArray();
        return DirectionalChainTests.DenseLogMarginal(y, [.. terms.Select(t => HandColumn(data, t))], model.Target.HeldTau!.Value, 3.0, 2.0);
    }

    private static BassTerm Gated(ProjectionDirection direction, double knot, params CategoricalFactor[] cells)
        => new(categorical: cells, directional: new DirectionalHinge(direction, knot));

    private static BassTerm Within(DirectionalTermSpace space, ProjectionDirection direction, double at, params CategoricalFactor[] cells)
    {
        KnotRange r = space.KnotRangeOf(direction, cells)!.Value;
        return Gated(direction, r.Lower + at * r.Width, cells);
    }

    private static readonly ProjectionDirection Pair = new([0, 1], [0.6, 0.8]);
    private static readonly ProjectionDirection Across = new([1, 2], [-0.8, 0.6]);
    private static readonly ProjectionDirection E0 = new([0], [1.0]);
    private static readonly ProjectionDirection MinusE2 = new([2], [-1.0]);

    // ---- columns, knots and the cells ------------------------------------------------------------------------------------------------

    [Fact]
    public void GatedColumnsAndKnots_AreTheMaskTimesTheHandColumns_OverTheCells()
    {
        BassData data = Data();
        BassModel model = BassModel.Create(Config(), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        BassData continuous = Continuous(data);
        int[][] levels = Levels(data);
        double active = space.ActiveCases, worstRange = 0.0, worstColumn = 0.0, worstKnot = 0.0;
        Assert.Equal(16.0, active);   // min(20, n/10) at n = 160
        foreach (ProjectionDirection d in new[] { Pair, Across, E0, MinusE2 })
            foreach (CategoricalFactor[] cells in new[] { new[] { new CategoricalFactor(0, [0, 2]) }, [new CategoricalFactor(0, [1])] })
            {
                // The range over the cells by hand: U the type-7 quantile at 1 − n_active/n_S of their projections, L = U − (U − Q₀)/p₀.
                double[] projection = DirectionalChainTests.Projections(continuous, d), mask = Mask(levels, cells);
                double[] reached = [.. projection.Where((_, i) => mask[i] != 0.0).Order()];
                double upper = Quantile(reached, 1.0 - active / reached.Length), lower = upper - (upper - reached[0]) / space.LocalProbability;
                KnotRange range = space.KnotRangeOf(d, cells)!.Value;
                worstRange = Math.Max(worstRange, Math.Max(Math.Abs(range.Upper - upper), Math.Abs(range.Lower - lower)) / Math.Max(1.0, Math.Abs(upper)));
                Assert.Equal(-Math.Log(range.Width), space.LogKnot(new DirectionalHinge(d, range.Lower + 0.3 * range.Width), cells));

                // The gated hinge's column: the mask times slice 1's hand hinge.
                BassTerm term = Within(space, d, 0.4, cells);
                var column = new double[data.Response.Cases];
                model.Basis.CaseColumn(term, model.Inputs.Cases, column);
                double[] hand = HandColumn(data, term);
                for (int i = 0; i < column.Length; i++) worstColumn = Math.Max(worstColumn, Math.Abs(column[i] - hand[i]) / Math.Max(1.0, Math.Abs(hand[i])));

                // A ridge block's knots and columns over the cells: the type-7 quantiles at ℓ/K of the cells' projections above t₀, and the
                // mask times the paper's columns (slice 3's exact ones).
                const int K = 4;
                double t0 = term.Directional!.Knot;
                double[] above = [.. projection.Where((p, i) => p > t0 && mask[i] != 0.0).Order()];
                double[] knots = space.RidgeKnots(term.Directional, cells, K)!;
                for (int l = 0; l <= K; l++) worstKnot = Math.Max(worstKnot, Math.Abs(knots[l + 1] - Quantile(above, (double)l / K)) / Math.Max(1.0, Math.Abs(knots[l + 1])));
                BassModel ridges = BassModel.Create(Config(columns: K), data);
                IReadOnlyList<BassTerm> block = ridges.Expansion!.Expand(term)!;
                for (int l = 1; l <= K; l++)
                {
                    Assert.Equal(cells, block[l - 1].Categorical);
                    ridges.Basis.CaseColumn(block[l - 1], ridges.Inputs.Cases, column);
                    for (int i = 0; i < column.Length; i++)
                    {
                        double paper = mask[i] * RidgeColumnTests.Paper(knots, l, projection[i]).Value;
                        worstColumn = Math.Max(worstColumn, Math.Abs(column[i] - paper) / Math.Max(1.0, Math.Abs(paper)));
                    }
                }
            }
        // With no cells, the range is the ungated one, bit for bit.
        Assert.Equal(space.KnotRangeOf(Pair), space.KnotRangeOf(Pair, []));
        output.WriteLine($"ranges over the cells against the hand quantile: {worstRange:E2}; columns: {worstColumn:E2}; block knots: {worstKnot:E2}");
        Assert.True(worstRange <= 1e-14, $"{worstRange:E2}");
        Assert.True(worstColumn <= 1e-12, $"{worstColumn:E2}");
        Assert.True(worstKnot <= 1e-14, $"{worstKnot:E2}");
    }

    // x₀, x₁ lattices; x₂ 0.3 on every case of c₁'s level 1; c₁ at level 2 on 10 cases only; c₂ never 1 where c₁ = 0.
    private static BassData EdgeData(int n = 160)
    {
        var rng = new Xoshiro256PlusPlus(11);
        int[] c1 = [.. Enumerable.Range(0, n).Select(i => i % 16 == 0 ? 2 : i % 2)];
        int[] c2 = [.. Enumerable.Range(0, n).Select(i => c1[i] == 0 ? 0 : (i / 2) % 2)];
        double[][] x = [[.. Enumerable.Range(0, n).Select(i => (0.5 + i * Lattice[0]) % 1.0)], [.. Enumerable.Range(0, n).Select(i => (0.5 + i * Lattice[1]) % 1.0)],
                        [.. Enumerable.Range(0, n).Select(i => c1[i] == 1 ? 0.3 : (0.5 + i * Lattice[2]) % 1.0)]];
        double[] y = [.. Enumerable.Range(0, n).Select(i => x[0][i] + 0.05 * Normal.SampleStandard(rng))];
        return new BassData([.. Enumerable.Range(0, 3).Select(j => (BassPredictor)new ContinuousPredictor($"x{j}", x[j], new InputDomain(0.0, 1.0))),
                             new CategoricalPredictor("c1", c1, ["a", "b", "c"]), new CategoricalPredictor("c2", c2, ["d", "e"])], new ScalarResponse(y));
    }

    [Fact]
    public void ACellThatIsEmpty_TooSmall_OrTied_GivesNoRange_AndItsTermLiesOutsideTheSupport()
    {
        BassData data = EdgeData();
        BassModel model = BassModel.Create(Config(maxCategorical: 2), data), ridges = BassModel.Create(Config(columns: 4, maxCategorical: 2), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        var small = new[] { new CategoricalFactor(0, [2]) };                               // 10 cases, n_active 16
        var empty = new[] { new CategoricalFactor(0, [0]), new CategoricalFactor(1, [1]) };  // no case
        var tied = new[] { new CategoricalFactor(0, [1]) };                                 // x₂ is 0.3 on every case
        var e2 = new ProjectionDirection([2], [1.0]);
        // Each gives no range for t₀, so its term has log prior −∞ under either basis: outside the support, which the term measure, the prior
        // sampler and the moves meet as a null draw, never redrawn. A ridge block's knots need only two cases above t₀: the empty and the tied
        // cells leave none or a tie at the top, so their block is undefined too and the design refuses it for support; the small cell's block
        // exists, and its term is outside the support through its prior alone.
        foreach ((ProjectionDirection d, CategoricalFactor[] cells, string what, bool block) in
                 new[] { (E0, small, "too small", true), (E0, empty, "empty", false), (e2, tied, "tied", false) })
        {
            Assert.Null(space.KnotRangeOf(d, cells));
            BassTerm term = Gated(d, 0.0, cells);
            Assert.Equal(double.NegativeInfinity, model.Prior.LogPrior(term));
            Assert.Equal(double.NegativeInfinity, ridges.Prior.LogPrior(term));
            Assert.Equal(block, ridges.Expansion!.Expand(term) is not null);
            if (!block) Assert.Equal(BassRefusal.Support, ridges.NewDesign().Load([term]));
            output.WriteLine($"{what}: no range, log prior −∞{(block ? "; its block exists" : "; no block, a support refusal")}");
        }
        Assert.NotNull(space.KnotRangeOf(e2, [new CategoricalFactor(0, [0])]));   // the same direction over other cells has one
    }

    // ---- the term measure ------------------------------------------------------------------------------------------------------------

    [Fact]
    public void TheTermMeasure_IsBothPartsByHand_AndItsDrawsFollowIt()
    {
        BassData data = Data();
        BassModel model = BassModel.Create(Config(), data);
        BassStructurePrior prior = model.Prior;
        DirectionalTermSpace space = prior.Directional!;
        // ln π(term) by hand: a uniform on {0, 1, 2}, 𝒥 uniform among C(3, a), θ_𝒥 on its subsphere, t uniform on the cells' range; k uniform
        // on {0, 1}, the one variable, its subset by R's size-uniform measure: size uniform on {1, 2}, then a subset of that size.
        double LogSubsphere(int a) => SpecialFunctions.LogGamma(0.5 * a) - Math.Log(2.0) - 0.5 * a * Math.Log(Math.PI);
        double Hand(BassTerm t)
        {
            double lp = -Math.Log(3.0) - Math.Log(2.0);
            if (t.Directional is { } h)
                lp += -Math.Log(new[] { 1, 3, 3 }[h.Direction.Count]) + LogSubsphere(h.Direction.Count) - Math.Log(space.KnotRangeOf(h.Direction, t.Categorical)!.Value.Width);
            foreach (CategoricalFactor c in t.Categorical) lp += -Math.Log(2.0) - Math.Log(new[] { 1, 3, 3 }[c.Levels.Count]);
            return lp;
        }
        double worst = 0.0;
        foreach (BassTerm t in new[] { new BassTerm(categorical: [new CategoricalFactor(0, [1])]), Within(space, Pair, 0.3, new CategoricalFactor(0, [0, 2])),
                                       Within(space, E0, 0.7), Within(space, MinusE2, 0.5, new CategoricalFactor(0, [2])) })
            worst = Math.Max(worst, Math.Abs(prior.LogPrior(t) - Hand(t)));
        output.WriteLine($"ln π(term) against the hand measure: worst |Δ| {worst:E2}");
        Assert.True(worst <= 1e-13, $"{worst:E2}");

        // 60,000 draws: (a, k) uniform over the six pairs, (0, 0) the intercept (a null draw, no other kind of null at these cells); the subset's
        // size uniform on {1, 2}.
        var rng = new Xoshiro256PlusPlus(3721);
        const int draws = 60_000;
        var seen = new int[3, 2];
        int sized = 0, single = 0;
        for (int k = 0; k < draws; k++)
        {
            if (prior.DrawTerm(rng) is not { } t)
            {
                seen[0, 0]++;
                continue;
            }
            seen[t.Directional?.Direction.Count ?? 0, t.Categorical.Count]++;
            foreach (CategoricalFactor c in t.Categorical)
            {
                sized++;
                if (c.Levels.Count == 1) single++;
            }
        }
        double chi = 0.0, expected = draws / 6.0;
        foreach (int count in seen) chi += (count - expected) * (count - expected) / expected;
        double z = (single - 0.5 * sized) / Math.Sqrt(0.25 * sized);
        output.WriteLine($"(a, k) over 60,000 draws: χ² {chi:F2} on 5; subsets of one level {single} of {sized}, z {z:F2}");
        Assert.True(chi < 20.5, $"χ² {chi:F2}");   // 5 degrees of freedom: p ≈ 10⁻³
        Assert.True(Math.Abs(z) < 4.0, $"{z:F2}");
    }

    // ---- each move's ratio -----------------------------------------------------------------------------------------------------------

    [Fact]
    public void EachMovesRatio_WithGates_MatchesADirectEvaluation()
    {
        BassData data = Data();
        BassModel model = BassModel.Create(Config(), data);
        BassStructurePrior prior = model.Prior;
        DirectionalTermSpace space = prior.Directional!;
        IBassProposalKernel kernel = model.Kernel;
        BassCountPrior count = model.Config.Count;
        double w1 = model.Config.Structure.OrderWeight, w2 = model.Config.Structure.VariableWeight;
        double Target(BassTerm[] s)
        {
            BassEvaluation e = DirectionalChainTests.Evaluate(model, s);
            return e.LogPrior + e.LogMarginal;
        }

        BassTerm[] s = [Within(space, Pair, 0.6, new CategoricalFactor(0, [0, 2])), new BassTerm(categorical: [new CategoricalFactor(0, [1])])];
        var counts = new BassTermCounts(prior);
        foreach (BassTerm t in s) counts.Add(t);
        // The birth's law by hand over both parts: q(a) = (w₁ + c_a)/(3w₁ + M), 𝒥 uniform at a = 1 and Wallenius' (by quadrature) at a = 2 with
        // z_u = w₂ + (terms using u); q(k) = (w₁ + c_k)/(2w₁ + M), the one variable surely. The prior's: 1/3, 1/C(3, a), 1/2, 1.
        double LogQ(BassTerm t)
        {
            int a = t.Directional?.Direction.Count ?? 0, k = t.Categorical.Count, m = counts.Terms;
            double lq = Math.Log(w1 + counts.Order(BassPart.Directional, a)) - Math.Log(3.0 * w1 + m)
                      + Math.Log(w1 + counts.Order(BassPart.Categorical, k)) - Math.Log(2.0 * w1 + m);
            if (a == 1) lq -= Math.Log(3.0);
            else if (a == 2)
            {
                double[] z = [.. Enumerable.Range(0, 3).Select(u => w2 + counts.Appearances(BassPart.Directional, u))];
                lq += Math.Log(DirectionalChainTests.WalleniusByQuadrature(z, [.. t.Directional!.Direction.Variables]));
            }
            return lq;
        }
        double LogPi(BassTerm t) => -Math.Log(3.0) - Math.Log(2.0) - (t.Directional is { } h ? Math.Log(new[] { 1, 3, 3 }[h.Direction.Count]) : 0.0);
        double worst = 0.0;
        foreach (BassTerm born in new[] { Within(space, E0, 0.55, new CategoricalFactor(0, [0])), Within(space, Across, 0.5),
                                          new BassTerm(categorical: [new CategoricalFactor(0, [2])]), Within(space, MinusE2, 0.45, new CategoricalFactor(0, [1, 2])) })
        {
            BassTerm[] grown = [.. s, born];
            int m = s.Length;
            double engine = Target(grown) - Target(s) + kernel.LogBirthRatio(born, counts);
            double direct = Dense(model, data, grown) - Dense(model, data, s) + Math.Log((m + count.Shape) / ((m + 1.0) * (count.Rate + 1.0)))
                          + LogPi(born) - LogQ(born);
            BassTermCounts after = counts.Clone();
            after.Add(born);
            double death = Target(s) - Target(grown) + kernel.LogDeathRatio(born, after);
            worst = Math.Max(worst, Math.Max(Math.Abs(engine - direct), Math.Abs(death + direct)));
            output.WriteLine($"birth of {born}: engine {engine:F12}, direct {direct:F12}; death {death:F12}");
        }

        // The changes: a gated term's direction or subset, each with probability ½; a categorical term's subset; each on the likelihood alone.
        var rng = new Xoshiro256PlusPlus(3722);
        int directions = 0, subsets = 0;
        foreach (int index in new[] { 0, 1 })
            for (int k = 0; k < 8; k++)
            {
                BassTermMove move = kernel.ProposeChange(rng, s[index], counts);
                if (move.Term is not { } moved) continue;
                if (move.Kind == BassMoveKind.DirectionChange) directions++;
                else subsets++;
                Assert.Equal(s[index].Directional is null, moved.Directional is null);
                BassTerm[] changed = [.. s];
                changed[index] = moved;
                double engine = Target(changed) - Target(s) + move.LogProposalRatio, direct = Dense(model, data, changed) - Dense(model, data, s);
                worst = Math.Max(worst, Math.Abs(engine - direct));
                output.WriteLine($"{move.Kind} of {s[index]} to {moved}: engine {engine:F12}, direct {direct:F12}");
            }
        output.WriteLine($"worst |Δ| over the moves {worst:E2}; direction changes {directions}, subset changes {subsets}");
        Assert.True(directions > 0 && subsets > 0);
        Assert.True(worst <= 1e-9, $"{worst:E2}");
    }

    // ---- Geweke's test with the prior alone ------------------------------------------------------------------------------------------

    private static double[] Functions(BassState state, DirectionalTermSpace space)
    {
        var f = new double[13];
        f[0] = state.Count;
        foreach (BassTerm t in state.Terms)
        {
            if (t.Directional is { } h)
            {
                f[h.Direction.Count]++;
                foreach (int v in h.Direction.Variables) f[3 + v]++;
                if (h.Direction.Weights[0] > 0.0) f[6]++;
                KnotRange r = space.KnotRangeOf(h.Direction, t.Categorical)!.Value;
                f[7] += (r.Upper - h.Knot) / r.Width;
            }
            else f[8]++;
            foreach (CategoricalFactor c in t.Categorical)
            {
                f[9]++;
                if (c.Levels.Count == 1) f[10]++;
                if (c.Contains(0)) f[11]++;
                if (c.Contains(2)) f[12]++;
            }
        }
        return f;
    }

    private static readonly string[] Statistics =
        ["M", "a = 1", "a = 2", "uses x0", "uses x1", "uses x2", "θ₁ > 0", "knot position", "categorical alone", "k = 1", "one level", "holds a", "holds c"];

    private double Geweke(BassModel model, int seed, Func<BassState, double[]> functions, string[] names)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var sampler = new BassPriorSampler(model);
        const int draws = 40_000, burn = 5_000, steps = 200_000, batches = 40;
        var direct = new double[draws][];
        for (int k = 0; k < draws; k++) direct[k] = functions(sampler.DrawState(rng).State);
        var chain = new BassChain(model, [], Math.Log(model.Target.HeldTau!.Value), new Xoshiro256PlusPlus(seed + 1), new Xoshiro256PlusPlus(seed + 2));
        chain.Moves.SetTemperature(new Temperature(0.0, TemperingPath.Likelihood));
        for (int t = 0; t < burn; t++) chain.Step();
        var successive = new double[steps][];
        for (int t = 0; t < steps; t++) successive[t] = functions(chain.Step());
        double worst = 0.0;
        for (int g = 0; g < names.Length; g++)
        {
            double[] a = [.. direct.Select(r => r[g])], b = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1);
            double meanB = b.Average(), se = Math.Sqrt(varA / a.Length + Math.Pow(DirectionalChainTests.BatchError(b, batches), 2));
            double z = se > 0.0 ? (meanA - meanB) / se : 0.0;
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{names[g],-18} prior {meanA,9:F5}  chain {meanB,9:F5}  z {z,6:F2}");
        }
        output.WriteLine($"worst |z| {worst:F2}; moves {string.Join(", ", chain.Structure.Kinds.Entries.Where(e => e.Attempts > 0).Select(e => $"{e.Kind} {e.Accepted}/{e.Attempts}"))}; "
                         + $"prior sampler acceptance {sampler.AcceptanceRate:F3}, determinant rejections {sampler.DeterminantRejections}");
        return worst;
    }

    [Fact]
    public void Geweke_WithGates_ThePriorChainLeavesThePriorInvariant()
    {
        BassModel model = BassModel.Create(Config(maxBasis: 4), Data(80));
        DirectionalTermSpace space = model.Prior.Directional!;
        double worst = Geweke(model, 3723, s => Functions(s, space), Statistics);
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ---- the determinant counts directions --------------------------------------------------------------------------------------------

    // One continuous input and the categorical one: a ridge and a categorical contrast are a valid two-term structure under the determinant.
    private static BassData OneInput(int n = 80)
    {
        var rng = new Xoshiro256PlusPlus(13);
        double[] x = [.. Enumerable.Range(0, n).Select(i => (0.5 + i * Lattice[0]) % 1.0)];
        int[] c = [.. Enumerable.Range(0, n).Select(i => i % 3)];
        double[] y = [.. Enumerable.Range(0, n).Select(i => x[i] * (c[i] == 1 ? 1.0 : 0.0) + 0.05 * Normal.SampleStandard(rng))];
        return new BassData([new ContinuousPredictor("x0", x, new InputDomain(0.0, 1.0)), new CategoricalPredictor("c", c, ["a", "b", "c"])], new ScalarResponse(y));
    }

    [Fact]
    public void TheDeterminant_CountsDirectionBearingTerms_ForItsNormalizerAndCap()
    {
        BassModel model = BassModel.Create(Config(directions: DirectionPrior.Determinant, maxBasis: 4), OneInput());
        DirectionalTermSpace space = model.Prior.Directional!;
        // M keeps its cap: terms without a direction take no part in det, so only M_θ ≤ p = 1 binds, through det.
        Assert.Equal(4, model.MaxBasis);
        Assert.Contains(model.Resolved.AutoPicks, a => a.Setting == "Directions" && a.Reason.Contains("M_θ ≤ p = 1"));
        BassTerm ridge = Within(space, E0, 0.5), contrast = new(categorical: [new CategoricalFactor(0, [1])]);
        Assert.Equal(0.0, DirectionDeterminant.LogCoupling([ridge, contrast], 1));   // ln 1 − ln c₁
        Assert.Equal(model.Prior.LogPrior(ridge) + model.Prior.LogPrior(contrast), model.Prior.LogTermsPrior([ridge, contrast]));
        Assert.Equal(double.NegativeInfinity, model.Prior.LogTermsPrior([ridge, Within(space, new ProjectionDirection([0], [-1.0]), 0.5), contrast]));

        string[] names = ["M", "M_θ", "M ≥ 2", "categorical alone", "θ₁ > 0", "knot position"];
        double[] Functions(BassState state)
        {
            var f = new double[names.Length];
            f[0] = state.Count;
            foreach (BassTerm t in state.Terms)
            {
                if (t.Directional is { } h)
                {
                    f[1]++;
                    if (h.Direction.Weights[0] > 0.0) f[4]++;
                    KnotRange r = space.KnotRangeOf(h.Direction, t.Categorical)!.Value;
                    f[5] += (r.Upper - h.Knot) / r.Width;
                }
                else f[3]++;
            }
            f[2] = state.Count >= 2 ? 1.0 : 0.0;
            Assert.True(f[1] <= 1.0);
            return f;
        }
        double worst = Geweke(model, 3724, Functions, names);
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ---- blocks of varying width ------------------------------------------------------------------------------------------------------

    [Fact]
    public void ARidgeAndACategoricalTerm_LoadInEitherOrder_ReloadFromTheirColumns_AndCountTwo()
    {
        BassData data = Data();
        BassModel model = BassModel.Create(Config(columns: 4, maxBasis: 2), data);
        DirectionalTermSpace space = model.Prior.Directional!;
        BassTerm ridge = Within(space, Pair, 0.4), contrast = new(categorical: [new CategoricalFactor(0, [1])]);
        double? residual = null;
        foreach (BassTerm[] units in new[] { new[] { ridge, contrast }, [contrast, ridge] })
        {
            BassDesign design = model.NewDesign();
            Assert.Null(design.Load(units));
            Assert.Equal(units, design.Terms);
            Assert.Equal(5, design.ColumnTerms.Count);
            Assert.Equal(2, design.Units);
            // Reloaded from its columns, as a draw's restore does: the same terms, columns and fit, bit for bit.
            BassDesign again = model.NewDesign();
            Assert.Null(again.Load([.. design.ColumnTerms]));
            Assert.Equal(units, again.Terms);
            Assert.Equal(design.ColumnTerms, again.ColumnTerms);
            Assert.Equal(BitConverter.DoubleToInt64Bits(design.ResidualSquares), BitConverter.DoubleToInt64Bits(again.ResidualSquares));
            // M through the table the ids index: 2, not 5/4.
            var table = new BassTermTable();
            int[] ids = [.. design.ColumnTerms.Select(table.Intern)];
            Assert.Equal(2, new BassStructure(ids).TermCount(table));
            // The fit's order does not move it beyond rounding.
            residual ??= design.ResidualSquares;
            Assert.True(Math.Abs(design.ResidualSquares - residual.Value) <= 1e-9 * residual.Value, $"{design.ResidualSquares} against {residual}");
            // M's cap on structure terms is 2, apart from the columns' 8: a third term is a cap refusal though its column would fit.
            Assert.Equal(BassRefusal.Cap, design.ProposeBirth(new BassTerm(categorical: [new CategoricalFactor(0, [2])])).Refusal);
        }
    }
}
