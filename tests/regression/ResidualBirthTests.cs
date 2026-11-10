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
using static Maths.Regression.Tests.BassKernelEnumeration;
using static Maths.Regression.Tests.BassTemperingTests;

namespace Maths.Regression.Tests;

/// <summary>
/// E4's residual birth (the proposal-kernels brief's slice 2): a parent uniformly over P(S), then a new factor ∝ exp(½·b·Δℓ), Δℓ its
/// column's profile-likelihood gain on the structure's least-squares residual, in the birth mixture beside the fit's own birth.
/// <list type="bullet">
/// <item><b>The statement, written out here:</b> the residual by modified Gram–Schmidt on dense columns, each candidate's weight from
/// its dense column, and q_R(t | S) summed over the parents the statement names. Production's tables, residual and probability are
/// read against it.</item>
/// <item><b>Primitives:</b> every candidate's score from the cumulative-sum pass against its dense column, on ties at knots, hinge
/// degrees 1 to 3, a far input range, columns that vanish on the cases, an exact fit, and a seeded sweep of structures.</item>
/// <item><b>The kernel's own law, detailed balance exactly, the chain's law and Geweke,</b> on hinge and Legendre fits, at b = 1 and
/// on a tempered rung.</item>
/// <item><b>The reach:</b> the ladder's rung laws (each rung's b), SMC's evidence (each stage's b), the flat histogram, resume, and
/// the refusals.</item>
/// </list>
/// </summary>
public sealed class ResidualBirthTests(ITestOutputHelper output)
{
    private static readonly ResidualBirth Spec = new(Share: 0.3);

    internal static BassConfig WithResidual(BassConfig config, ResidualBirth? spec = null, HistoryBirth? history = null)
        => config with { Sampling = config.Sampling with { Births = new BirthMixture { Residual = spec ?? Spec, History = history } } };

    private static BassModel Fixture(string name, bool births = false)
    {
        (BassConfig config, BassData data) = name switch
        {
            "legendre" => (KhaosChainTests.Config(raw: true), KhaosChainTests.Data()),
            "three inputs" or "quadratic" or "cubic" or "far" => (new BassConfig
            {
                Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 1, Degree = name == "quadratic" ? 2 : name == "cubic" ? 3 : 1 },
                Coefficients = BassCoefficientPrior.G(TauPrior.Held()),
            }, ThreeData(name == "far" ? 1e6 : 0.0)),
            _ => BassLocalStepTests.Fixture(name, "equal", s => s),
        };
        return BassModel.Create(births ? WithResidual(config) : config, data);
    }

    // Three continuous inputs on six cases, each a permutation of the grid, shifted by `offset`: y = 1 + 2x₀ + x₁x₂ + noise.
    private static BassData ThreeData(double offset)
    {
        double[] x0 = [0, 1, 2, 3, 4, 5], x1 = [3, 0, 4, 1, 5, 2], x2 = [5, 3, 1, 4, 2, 0];
        for (int i = 0; i < 6; i++) (x0[i], x1[i], x2[i]) = (x0[i] / 5, x1[i] / 5, x2[i] / 5);
        double[] y = [.. Enumerable.Range(0, 6).Select(i => 1 + 2 * x0[i] + x1[i] * x2[i] + 0.3 * Math.Sin(7 * i))];
        return new BassData([new ContinuousPredictor("x0", [.. x0.Select(v => v + offset)]), new ContinuousPredictor("x1", [.. x1.Select(v => v + offset)]),
                             new ContinuousPredictor("x2", [.. x2.Select(v => v + offset)])], new ScalarResponse(y));
    }

    private static double LogTau(BassModel model) => Math.Log(model.Target.HeldTau!.Value);

    private static BassPart Main(BassModel model) => model.Prior.Has(BassPart.Legendre) ? BassPart.Legendre : BassPart.Continuous;

    // ---- the statement, written out ------------------------------------------------------------------------------------------------

    private static double[] Column(BassModel model, BassTerm? term)
    {
        int n = model.Data.Response.Cases;
        var c = new double[n];
        if (term is null) Array.Fill(c, 1.0);
        else model.Basis.CaseColumn(term, model.Inputs.Cases, c);
        return c;
    }

    // r(S): the raw response less its projection on [1, B_S], by modified Gram–Schmidt.
    internal static double[] DenseResidual(BassModel model, IReadOnlyList<BassTerm> terms)
    {
        double[] r = model.Data.Response.Values.ToArray();
        var basis = new List<double[]>();
        foreach (double[] raw in new[] { Column(model, null) }.Concat(terms.Select(t => Column(model, t))))
        {
            double[] v = (double[])raw.Clone();
            foreach (double[] q in basis)
            {
                double d = 0.0;
                for (int i = 0; i < v.Length; i++) d += q[i] * v[i];
                for (int i = 0; i < v.Length; i++) v[i] -= d * q[i];
            }
            double norm = Math.Sqrt(v.Sum(x => x * x));
            for (int i = 0; i < v.Length; i++) v[i] /= norm;
            basis.Add(v);
        }
        foreach (double[] q in basis)
        {
            double d = 0.0;
            for (int i = 0; i < r.Length; i++) d += q[i] * r[i];
            for (int i = 0; i < r.Length; i++) r[i] -= d * q[i];
        }
        return r;
    }

    // ½·b·Δℓ for column c on r: Δℓ = −(n/2)ln(1 − (c̃ᵀr)²/(‖c̃‖²‖r‖²)), the share capped below 1, 0 for a column constant on the cases
    // and for an exact fit (‖r‖² at most 10⁻²⁴ of the response's spread).
    internal static double StatementLogWeight(BassModel model, double[] c, double[] r, double b)
    {
        int n = c.Length;
        double[] y = model.Data.Response.Values.ToArray();
        double ym = y.Average(), spread = y.Sum(v => (v - ym) * (v - ym));
        double rss = r.Sum(v => v * v);
        if (!(rss > 1e-24 * spread)) return 0.0;
        double cm = c.Average();
        double cross = 0.0, cc = 0.0, raw = 0.0;
        for (int i = 0; i < n; i++)
        {
            cross += (c[i] - cm) * r[i];
            cc += (c[i] - cm) * (c[i] - cm);
            raw += c[i] * c[i];
        }
        if (!(cc > 1e-13 * Math.Max(raw, 1e-300))) return 0.0;
        double share = Math.Min(cross * cross / cc / rss, 1.0 - 1e-12);
        return 0.5 * b * (-0.5 * n * Math.Log(1.0 - share));
    }

    private static int MainOrder(BassModel model, BassTerm t) => Main(model) == BassPart.Legendre ? t.Legendre.Count : t.Continuous.Count;

    // P(S) as the statement names it: the intercept, and the terms that can take one more main-part factor.
    internal static List<BassTerm?> StatementParents(BassModel model, IReadOnlyList<BassTerm> terms)
    {
        BassPart main = Main(model);
        int cap = model.Prior.Cap(main), vars = model.Prior.Variables(main), dmax = model.Prior.Legendre?.MaxDegree ?? int.MaxValue;
        var parents = new List<BassTerm?>();
        if (cap >= 1) parents.Add(null);
        foreach (BassTerm t in terms)
        {
            int order = MainOrder(model, t);
            if (order < cap && order < vars && (main != BassPart.Legendre || t.LegendreDegree < dmax)) parents.Add(t);
        }
        return parents;
    }

    // The statement's candidates of `parent`: every variable it leaves, with every sign and knot, or every degree within d_max − d_p.
    internal static List<BassTerm> Candidates(BassModel model, BassTerm? parent)
    {
        BassPart main = Main(model);
        var list = new List<BassTerm>();
        for (int u = 0; u < model.Prior.Variables(main); u++)
        {
            if (parent is not null && BassTermCounts.Variables(parent, main).Contains(u)) continue;
            if (main == BassPart.Legendre)
            {
                int budget = model.Prior.Legendre!.MaxDegree - (parent?.LegendreDegree ?? 0);
                for (int d = 1; d <= budget; d++)
                    list.Add(new BassTerm(legendre: parent is null ? [new LegendreFactor(u, d)] : [.. parent.Legendre, new LegendreFactor(u, d)]));
                continue;
            }
            for (int s = -1; s <= 1; s += 2)
                for (int k = 0; k < model.Inputs.ContinuousKnots(u).Length; k++)
                {
                    var h = new HingeFactor(u, s, k);
                    list.Add(parent is null ? new BassTerm(continuous: [h]) : new BassTerm([.. parent.Continuous, h], parent.Categorical, parent.Functional));
                }
        }
        return list;
    }

    // q_R(· | S) as the statement gives it, over every child of P(S): (1/|P(S)|)·Σ_p w_p(t)/W_p, the weights from dense columns on the
    // dense residual.
    internal static Dictionary<BassTerm, double> StatementLaw(BassModel model, IReadOnlyList<BassTerm> terms, double b)
    {
        var law = new Dictionary<BassTerm, double>();
        List<BassTerm?> parents = StatementParents(model, terms);
        if (parents.Count == 0) return law;
        double[] r = DenseResidual(model, terms);
        foreach (BassTerm? p in parents)
        {
            List<BassTerm> candidates = Candidates(model, p);
            double[] lw = [.. candidates.Select(c => StatementLogWeight(model, Column(model, c), r, b))];
            double top = lw.Max(), total = lw.Sum(v => Math.Exp(v - top));
            for (int k = 0; k < candidates.Count; k++)
                law[candidates[k]] = law.GetValueOrDefault(candidates[k]) + Math.Exp(lw[k] - top) / total / parents.Count;
        }
        return law;
    }

    private static double LogOf(Dictionary<BassTerm, double> law, BassTerm t) => law.TryGetValue(t, out double q) ? Math.Log(q) : double.NegativeInfinity;

    // ---- primitives ----------------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData("continuous")]
    [InlineData("three inputs")]
    [InlineData("quadratic")]
    [InlineData("cubic")]
    [InlineData("far")]
    [InlineData("legendre")]
    [InlineData("categorical, R's measure")]
    public void Scores_AreTheDenseColumns_AndTheResidualIsTheProjection(string fixture)
    {
        BassModel model = Fixture(fixture);
        var residual = new BassResidualBirth(model);
        List<BassTerm> terms = AllTerms(model);
        var rng = new Xoshiro256PlusPlus(20261201);
        var structures = new List<BassTerm[]> { Array.Empty<BassTerm>() };
        for (int k = 0; structures.Count < 12 && k < 400; k++)
        {
            int m = 1 + rng.NextInt(model.MaxBasis);
            BassTerm[] s = [.. Enumerable.Range(0, m).Select(_ => terms[rng.NextInt(terms.Count)]).Distinct()];
            if (BassMixing.Loaded(model, s, LogTau(model), out _) is not null) structures.Add(s);
        }
        double worstResidual = 0.0, worstScore = 0.0;
        int scored = 0, zero = 0;
        foreach (BassTerm[] s in structures)
        {
            double[] dense = DenseResidual(model, s), production = residual.Residual(s);
            double scale = Math.Sqrt(dense.Sum(v => v * v)) + 1e-300;
            worstResidual = Math.Max(worstResidual, dense.Zip(production, (a, c) => Math.Abs(a - c)).Max() / scale);
            // The residual is the structure's, whatever its terms' order, bit for bit; P(S) is the statement's.
            Assert.Equal(production, new BassResidualBirth(model).Residual([.. s.Reverse()]));
            Assert.Equal(StatementParents(model, s), residual.ParentsOf(s));
            foreach (BassTerm? p in residual.ParentsOf(s))
                foreach (double b in new[] { 1.0, 0.3 })
                {
                    double[] lw = residual.LogWeights(s, p, b);
                    List<BassTerm> candidates = Candidates(model, p);
                    var seen = new HashSet<BassTerm>();
                    for (int slot = 0; slot < lw.Length; slot++)
                    {
                        if (double.IsNegativeInfinity(lw[slot])) continue;
                        BassTerm child = residual.ChildOf(p, slot);
                        Assert.True(seen.Add(child), $"slot {slot} repeats {child}");
                        double expected = StatementLogWeight(model, Column(model, child), dense, b);
                        worstScore = Math.Max(worstScore, Math.Abs(lw[slot] - expected) / (1.0 + Math.Abs(expected)));
                        scored++;
                        if (expected == 0.0) zero++;
                    }
                    Assert.True(seen.SetEquals(candidates), $"the slots of {p?.ToString() ?? "the intercept"} are not the statement's candidates");
                }
        }
        output.WriteLine($"{fixture}: {structures.Count} structures, {scored} scores ({zero} zero); residual worst {worstResidual:E2} of ‖r‖; score worst {worstScore:E2} relative");
        Assert.True(worstResidual <= 1e-9, $"the residual sits {worstResidual:E2} from the projection");
        Assert.True(worstScore <= 1e-7, $"a score sits {worstScore:E2} from its dense column's");
    }

    [Fact]
    public void AnExactFit_WeighsEveryCandidateAlike()
    {
        BassModel probe = Fixture("three inputs");
        BassTerm t = AllTerms(probe).First(x => x.Continuous.Count == 1 && x.Continuous[0].Sign == 1 && x.Continuous[0].Knot == 1);
        double[] column = Column(probe, t);
        double[] y = [.. column.Select(v => 1.0 + 2.0 * v)];
        BassData data = new(probe.Data.Predictors, new ScalarResponse(y));
        BassModel model = BassModel.Create(probe.Config, data);
        var residual = new BassResidualBirth(model);
        foreach (BassTerm? p in residual.ParentsOf([t]))
            Assert.All(residual.LogWeights([t], p, 1.0).Where(double.IsFinite), v => Assert.Equal(0.0, v));
    }

    // ---- the kernel's own law ------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData("three inputs", 20261210)]
    [InlineData("continuous", 20261211)]
    [InlineData("categorical, R's measure", 20261212)]
    [InlineData("legendre", 20261213)]
    public void Law_SumsToOne_IsTheStatement_AndTheDrawsFollowIt(string fixture, int seed)
    {
        BassModel model = Fixture(fixture);
        var residual = new BassResidualBirth(model);
        var mixture = new BassBirthMixture(model.Kernel, null, null, new BassResidualBirth(model), Spec);
        List<BassTerm> terms = AllTerms(model);
        // Structures of F (the design admits them, so the statement's projection is defined): the empty one, the first and last
        // single terms that load, and a pair.
        bool InF(BassTerm[] s) => BassMixing.Loaded(model, s, LogTau(model), out _) is not null;
        var structures = new List<BassTerm[]> { Array.Empty<BassTerm>() };
        structures.Add(terms.Select(t => new[] { t }).First(InF));
        structures.Add(terms.Select(t => new[] { t }).Last(InF));
        if (model.MaxBasis >= 2) structures.Add(terms.SelectMany(a => terms.Where(c => !c.Equals(a)).Select(c => new[] { a, c })).Skip(7).First(InF));
        double worst = 0.0;
        foreach (BassTerm[] s in structures)
            foreach (double b in new[] { 1.0, 0.4 })
            {
                var at = new BassBirthContext(s, b);
                BassTermCounts counts = BassMixing.Counts(model, s);
                Dictionary<BassTerm, double> law = StatementLaw(model, s, b);
                double empty = model.Prior.Legendre is null ? Math.Exp(model.Proposal.LogEmptyBirth(counts)) : 0.0;
                double alone = 0.0, mixed = (1.0 - Spec.Share) * empty;
                foreach (BassTerm t in terms)
                {
                    double lr = residual.LogProbability(t, at);
                    alone += Math.Exp(lr);
                    mixed += Math.Exp(mixture.LogBirth(t, counts, at));
                    double statement = LogOf(law, t);
                    if (double.IsNegativeInfinity(statement)) Assert.True(double.IsNegativeInfinity(lr), $"{t} is no child of P(S), yet q_R = {Math.Exp(lr):G4}");
                    else worst = Math.Max(worst, Math.Abs(lr - statement));
                }
                output.WriteLine($"{fixture}, M = {s.Length}, b = {b}: Σq_R {alone:R}; the mixture's Σq with its empty draw {mixed:R}");
                Assert.Equal(1.0, alone, 1e-10);
                Assert.Equal(1.0, mixed, 1e-10);
            }
        output.WriteLine($"{fixture}: worst |ln q_R − the statement| {worst:E2}");
        Assert.True(worst <= 1e-8, $"ln q_R sits {worst:E2} from the statement");

        // The draws against q_R at the last structure, by χ², cells below five expected pooled.
        var last = new BassBirthContext(structures[^1], 1.0);
        const int draws = 100_000;
        var rng = new Xoshiro256PlusPlus(seed);
        var tally = new Dictionary<BassTerm, int>();
        for (int k = 0; k < draws; k++)
        {
            BassTerm t = residual.Draw(rng, last)!;
            tally[t] = tally.GetValueOrDefault(t) + 1;
        }
        List<(double Expected, int Observed)> cells = [.. terms.Select(t => (draws * Math.Exp(residual.LogProbability(t, last)), tally.GetValueOrDefault(t)))];
        Assert.Equal(draws, cells.Sum(c => c.Observed));
        Assert.Equal(0, cells.Where(c => c.Expected == 0.0).Sum(c => c.Observed));
        List<(double Expected, int Observed)> kept = [.. cells.Where(c => c.Expected >= 5)];
        double pooledExpected = cells.Where(c => c.Expected < 5).Sum(c => c.Expected);
        if (pooledExpected > 0) kept.Add((pooledExpected, cells.Where(c => c.Expected < 5).Sum(c => c.Observed)));
        double chi2 = kept.Sum(c => (c.Observed - c.Expected) * (c.Observed - c.Expected) / c.Expected);
        int df = kept.Count - 1;
        double z = (chi2 - df) / Math.Sqrt(2.0 * df);
        output.WriteLine($"{fixture}: {draws} draws over {kept.Count} cells, χ² {chi2:F1} on {df} df, z {z:F2}");
        Assert.True(z < 4.0, $"the draws sit z = {z:F2} from q_R");
    }

    // ---- detailed balance, exactly -------------------------------------------------------------------------------------------------

    // Births from the statement (the mixture written out: the base law, the history's when declared, the residual birth's from dense
    // columns) with the production mixture's ratio, each from its context; deaths uniform; changes the base kernel's.
    private static IEnumerable<Enumerated> Mixed(Space space, BassBirthMixture mixture, BassInclusionHistory laws, HistoryBirth? historySpec,
                                                 string key, State s)
    {
        BassModel model = space.Model;
        double b = space.InverseTemperature;
        int m = s.Terms.Length;
        double wh = historySpec?.Share ?? 0.0, wr = Spec.Share;
        switch (key)
        {
            case "birth":
                Dictionary<BassTerm, double> residual = StatementLaw(model, s.Terms, b);
                foreach (BassTerm t in space.Terms)
                {
                    if (s.Terms.Contains(t)) continue;
                    (double lb, double lh) = HistoryBirthTests.Laws(model, laws, t, s.Counts);
                    double q = (1.0 - wh - wr) * Math.Exp(lb) + wh * Math.Exp(lh) + wr * residual.GetValueOrDefault(t);
                    yield return new Enumerated(BirthDeath, [.. s.Terms, t], Math.Log(q), mixture.LogBirthRatio(t, s.Counts, new BassBirthContext(s.Terms, b)));
                }
                break;
            case "death":
                for (int k = 0; k < m; k++)
                {
                    BassTerm[] reduced = [.. s.Terms.Where((_, i) => i != k)];
                    yield return new Enumerated(BirthDeath, reduced, -Math.Log(m), mixture.LogDeathRatio(s.Terms[k], s.Counts, new BassBirthContext(reduced, b)));
                }
                break;
            default:
                foreach (Enumerated e in model.Prior.Legendre is null ? All(space, key, s, null) : KhaosChainTests.Khaos(space, key, s)) yield return e;
                break;
        }
    }

    [Theory]
    [InlineData("continuous", 1.0, false)]
    [InlineData("continuous", 0.4, false)]
    [InlineData("continuous", 1.0, true)]
    [InlineData("categorical, R's measure", 1.0, false)]
    [InlineData("three inputs", 1.0, false)]
    [InlineData("quadratic", 1.0, false)]
    [InlineData("legendre", 1.0, false)]
    [InlineData("legendre", 1.0, true)]
    public void DetailedBalance_Enumerated(string fixture, double inverseTemperature, bool withHistory)
    {
        BassModel model = Fixture(fixture);
        var historySpec = withHistory ? new HistoryBirth(Share: 0.2) : null;
        BassInclusionHistory laws = HistoryBirthTests.Skewed(model, historySpec ?? new HistoryBirth());
        var mixture = new BassBirthMixture(model.Kernel, historySpec, withHistory ? laws : null, new BassResidualBirth(model), Spec);
        Space space = model.Prior.Legendre is null ? Enumerate(model, inverseTemperature) : KhaosChainTests.Enumerate(model);
        var flows = Kernel(space, (key, s) => Mixed(space, mixture, laws, historySpec, key, s));
        (Dictionary<string, int> transitions, double worst) = Check(space, flows);
        output.WriteLine($"{fixture}, b = {inverseTemperature}{(withHistory ? ", with the history birth" : "")}: |F| = {space.States.Count}; " +
                         $"{string.Join(", ", transitions.Select(t => $"{t.Key} {t.Value}"))}; worst gap {worst:E2}");
        Assert.True(transitions[BirthDeath] > 0);
    }

    // The check's power: a residual birth whose death reads the full structure's residual, not the reduced one's, breaks balance.
    [Fact]
    public void DetailedBalance_Enumerated_RefusesTheFullStructuresResidualInADeath()
    {
        BassModel model = Fixture("continuous");
        var mixture = new BassBirthMixture(model.Kernel, null, null, new BassResidualBirth(model), Spec);
        Space space = Enumerate(model);
        BassInclusionHistory laws = HistoryBirthTests.Skewed(model, new HistoryBirth());
        IEnumerable<Enumerated> Wrong(string key, State s)
        {
            if (key != "death")
            {
                foreach (Enumerated e in Mixed(space, mixture, laws, null, key, s)) yield return e;
                yield break;
            }
            for (int k = 0; k < s.Terms.Length; k++)
                yield return new Enumerated(BirthDeath, [.. s.Terms.Where((_, i) => i != k)], -Math.Log(s.Terms.Length),
                                            mixture.LogDeathRatio(s.Terms[k], s.Counts, new BassBirthContext(s.Terms, 1.0)));
        }
        var flows = Kernel(space, Wrong);
        Assert.ThrowsAny<Xunit.Sdk.XunitException>(() => Check(space, flows));
    }

    // ---- the chain's law -----------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData("continuous", 20261220)]
    [InlineData("legendre", 20261221)]
    public void Chain_MeetsTheExactLaw(string fixture, int seed)
    {
        BassModel model = Fixture(fixture, births: true);
        Space space = model.Prior.Legendre is null ? Enumerate(model) : KhaosChainTests.Enumerate(model);
        double[] pi = space.Posterior();
        var chain = new BassChain(model, [], LogTau(model), new Xoshiro256PlusPlus(seed), new Xoshiro256PlusPlus(seed + 1));
        for (int t = 0; t < 5_000; t++) chain.Step();
        const int steps = 300_000, batches = 40;
        var visits = new int[steps];
        for (int t = 0; t < steps; t++)
        {
            BassState s = chain.Step();
            int i = space.Find(s.Terms);
            Assert.True(i >= 0, $"the chain holds a structure outside the enumeration: {string.Join(" + ", s.Terms)}");
            visits[t] = i;
        }
        BassPart main = Main(model);
        var readings = new List<(string Name, Func<State, bool> Holds)>();
        for (int mm = 0; mm <= model.MaxBasis; mm++) { int m = mm; readings.Add(($"P(M = {m})", st => st.Terms.Length == m)); }
        for (int v = 0; v < model.Prior.Variables(main); v++) { int j = v; readings.Add(($"P(x{j} in)", st => st.Terms.Any(t => BassTermCounts.Variables(t, main).Contains(j)))); }
        foreach (int top in Enumerable.Range(0, pi.Length).OrderByDescending(i => pi[i]).Take(6))
        {
            int k = top;
            readings.Add(($"P({(space.States[k].Terms.Length == 0 ? "∅" : string.Join(" + ", space.States[k].Terms))})", st => ReferenceEquals(st, space.States[k])));
        }
        double worst = 0.0;
        foreach ((string name, Func<State, bool> holds) in readings)
        {
            bool[] indicator = [.. space.States.Select(holds)];
            double exact = Enumerable.Range(0, pi.Length).Where(i => indicator[i]).Sum(i => pi[i]);
            double[] series = [.. visits.Select(i => indicator[i] ? 1.0 : 0.0)];
            double mean = series.Average(), se = KhaosChainTests.BatchError(series, batches);
            double z = se > 0.0 ? (mean - exact) / se : (Math.Abs(mean - exact) < 1e-3 ? 0.0 : double.PositiveInfinity);
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{name,-40} chain {mean:F4}  exact {exact:G6}  z {z,6:F2}");
        }
        BassKindCount births = chain.Structure.Kinds[BassMoveKind.Birth];
        output.WriteLine($"{fixture}: {steps} steps, worst |z| {worst:F2}; births {births.Accepted}/{births.Attempts}");
        Assert.True(worst < 4.5, $"the chain's worst |z| {worst:F2} against the exact law");
    }

    // ---- Geweke (qualification) ----------------------------------------------------------------------------------------------------

    private const int GewekeN = 12;

    private static BassConfig GewekeConfig(bool legendre) => WithResidual(legendre
        ? KhaosChainTests.Config(raw: true, maxBasis: 3, sigma2: Sigma2Prior.InverseGamma(3, 2)) with
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 3, Legendre = new LegendreBasis { MaxDegree = 3, DegreeMutation = 0.5 } },
        }
        : new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 3, Support = new SupportRule(SupportKind.Count, 2) },
            Coefficients = BassCoefficientPrior.G(TauPrior.Held()),
            Intercept = GPriorCentre.Zero,
            Sigma2 = Sigma2Prior.InverseGamma(3, 2),
            Sampling = new BassSampling { Chains = 1, MasterSeed = 1, BurnIn = 0, Samples = 1 },
        }, new ResidualBirth(Share: 0.5));

    private static BassData GewekeData(bool legendre, double[] y)
    {
        BassData d = KhaosBasisTests.UnitData(GewekeN, legendre ? 4 : 3, 77);
        return new BassData(d.Predictors, new ScalarResponse(y));
    }

    private static readonly string[] GewekeNames = ["M", "ln sigma2", "atan f0", "atan f6", "x0 in", "x2 in", "orders", "z0^2"];

    private static double[] Functions(BassModel model, BassState state, double sigma2, double[] f, double[] y)
    {
        BassPart main = Main(model);
        bool In(int v) => state.Terms.Any(t => BassTermCounts.Variables(t, main).Contains(v));
        return [state.Count, Math.Log(sigma2), Math.Atan(f[0]), Math.Atan(f[GewekeN / 2]), In(0) ? 1 : 0, In(2) ? 1 : 0,
                state.Terms.Sum(t => BassTermCounts.Variables(t, main).Count()), (y[0] - f[0]) * (y[0] - f[0]) / sigma2];
    }

    [Theory]
    [InlineData(false, 20261230)]
    [InlineData(true, 20261231)]
    public void Geweke_TheRedrawnResponseMovesTheResidual(bool legendre, int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int draws = 200_000, steps = 300_000, batches = 50;
        BassConfig config = GewekeConfig(legendre);
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(config, GewekeData(legendre, [.. Enumerable.Range(0, GewekeN).Select(i => Math.Sin(i))]));
        var sampler = new BassPriorSampler(prior);
        double logTau = LogTau(prior);
        (BassState S, double Sigma2, double[] F, double[] Y) Draw()
        {
            BassState state = sampler.DrawState(rng).State;
            var chain = new BassChain(prior, state.Terms, logTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
            BassSweepFit fit = chain.Structure.DrawSweepFit(chain.Current, rng, prior.Target, power: 0.0);
            double[] y = [.. Enumerable.Range(0, GewekeN).Select(i => fit.Fitted[i] + Math.Sqrt(fit.Sigma2) * Normal.SampleStandard(rng))];
            return (chain.Current, fit.Sigma2, fit.Fitted, y);
        }
        var marginal = new double[draws][];
        for (int k = 0; k < draws; k++)
        {
            (BassState s, double sigma2, double[] f, double[] y) = Draw();
            marginal[k] = Functions(prior, s, sigma2, f, y);
        }
        (BassState S, double Sigma2, double[] F, double[] Y) theta = Draw();
        var successive = new double[steps][];
        long births = 0, accepted = 0;
        for (int t = 0; t < steps; t++)
        {
            double[] y = [.. Enumerable.Range(0, GewekeN).Select(i => theta.F[i] + Math.Sqrt(theta.Sigma2) * Normal.SampleStandard(rng))];
            BassModel model = BassModel.Create(config, GewekeData(legendre, y));
            var chain = new BassChain(model, theta.S.Terms, logTau, new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)), new Xoshiro256PlusPlus(rng.NextInt(int.MaxValue)));
            chain.Step();
            BassSweepFit fit = chain.Structure.DrawSweepFit(chain.Current, rng, model.Target);
            successive[t] = Functions(model, chain.Current, fit.Sigma2, fit.Fitted, y);
            theta = (chain.Current, fit.Sigma2, fit.Fitted, y);
            chain.Structure.DesignAt(chain.Current);
            BassKindCount e = chain.Structure.Kinds[BassMoveKind.Birth];
            births += e.Attempts;
            accepted += e.Accepted;
        }
        double worst = 0;
        for (int g = 0; g < GewekeNames.Length; g++)
        {
            double[] a = [.. marginal.Select(r => r[g])], b = [.. successive.Select(r => r[g])];
            double meanA = a.Average(), varA = a.Sum(x => (x - meanA) * (x - meanA)) / (a.Length - 1);
            double meanB = b.Average(), se = Math.Sqrt(varA / a.Length + Math.Pow(KhaosChainTests.BatchError(b, batches), 2));
            double z = (meanA - meanB) / se;
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{GewekeNames[g],-10} marginal-conditional {meanA,10:F5}  successive {meanB,10:F5}  z {z,6:F2}");
        }
        output.WriteLine($"{(legendre ? "Legendre" : "hinge")}, seed {seed}: births accepted {accepted}/{births}; worst |z| {worst:F2}");
        Assert.True(worst < 4.0, $"Geweke's worst |z| {worst:F2}");
    }

    // ---- the reach -----------------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData("hinge", new[] { 1.0, 0.5, 0.2 })]
    [InlineData("legendre", new[] { 1.0, 0.4, 0.1, 0.0 })]
    public void Ladder_EachRungSamplesItsTemperedLaw_AtItsOwnB(string basis, double[] rungs)
    {
        bool reference = rungs[^1] == 0.0;
        var spec = new TemperingSpec { Levels = rungs.Length, HotEnd = reference ? rungs[^2] : rungs[^1], ReferenceRung = reference };
        BassConfig config = basis == "legendre" ? KhaosReachTests.Config("raw ridge") : BassTemperingTests.Config(spec);
        BassData data = basis == "legendre" ? KhaosChainTests.Data() : BassSweepTests.Data();
        BassModel model = BassModel.Create(WithResidual(config with { Sampling = config.Sampling with { Tempering = spec } }), data);
        const int ladders = 4, batches = 20;
        (double[][,] count, double[][] logTau) = RungBatches(model, spec, rungs, ladders, batches);
        RungLaw[] exact = Exact(model, rungs);
        for (int r = 0; r < rungs.Length; r++)
        {
            for (int m = 0; m <= model.MaxBasis; m++)
            {
                (double mean, double error) = Summary([.. Enumerable.Range(0, ladders * batches).Select(k => count[r][m, k])]);
                output.WriteLine($"{basis}, b = {rungs[r]:G3}, M = {m}: {mean:F4} against {exact[r].Count[m]:F4} (± {error:F4})");
                Assert.True(Math.Abs(mean - exact[r].Count[m]) <= 5 * error + 0.004, $"b = {rungs[r]}, M = {m}: {mean:F4} against {exact[r].Count[m]:F4} (± {error:F4})");
            }
            if (!model.Target.SamplesTau) continue;
            (double u, double uError) = Summary(logTau[r]);
            Assert.True(Math.Abs(u - exact[r].LogTauMean) <= 5 * uError + 0.01, $"b = {rungs[r]}, E[ln τ] {u:F4} against {exact[r].LogTauMean:F4} (± {uError:F4})");
        }
    }

    [Theory]
    [InlineData("raw ridge")]
    [InlineData("modified g")]
    public void Smc_AtEachStagesB_MeetsTheEnumeratedEvidence(string prior)
    {
        BassData data = KhaosChainTests.Data();
        BassConfig config = WithResidual(KhaosReachTests.Config(prior));
        config = config with { Sampling = config.Sampling with { MasterSeed = 3, ReportEvidence = true, Smc = new SmcSpec { Particles = 2000, Replicates = 4, MutationSteps = 4 } } };
        BassModel model = BassModel.Create(config, data);
        RungLaw exact = Exact(model, 1.0)[0];
        BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        BassEvidence evidence = fit.Evidence!;
        output.WriteLine($"{prior}, SMC with the residual birth: ln p(y) {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against {exact.LogZ:F4}");
        Assert.True(Math.Abs(evidence.LogEvidence - exact.LogZ) <= 4 * evidence.MonteCarloError + 0.01,
                    $"{prior}: {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against {exact.LogZ:F4}");
        Assert.Empty(fit.Diagnostics.Smc!.Flags);
    }

    [Fact]
    public void FlatHistogram_MatchesTheEnumeration()
        => BassFlatHistogramTests.MatchesTheEnumeration(BassModel.Create(WithResidual(KhaosReachTests.Config("raw ridge")), KhaosChainTests.Data()), 0.0,
                                                        "KHAOS with the residual birth, raw ridge", output, 6);

    [Theory]
    [InlineData("chain")]
    [InlineData("ladder")]
    [InlineData("both components")]
    public void Resume_IsBitForBit(string sampler)
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig Config(int samples)
        {
            BassConfig c = WithResidual(BassRunTests.Config(chains: 2, burn: 250, samples: samples, seed: 31),
                                        history: sampler == "both components" ? new HistoryBirth(Share: 0.2, Refresh: 50) : null);
            return c with { Sampling = c.Sampling with { Tempering = sampler == "ladder" ? new TemperingSpec { Levels = 3, ReferenceRung = true } : null } };
        }
        BassFit whole = Bass.Fit(Config(300), data);
        BassFit first = Bass.Fit(Config(150), data, new BassRunOptions { Checkpoint = true });
        BassFit resumed = Bass.Resume(first.Checkpoint!, Config(150), data);
        BassRunTests.AssertSameFit(whole, resumed);
        Assert.Equal(BassBitIdentityTests.Hash(whole), BassBitIdentityTests.Hash(resumed));
        BassKindCount births = whole.Diagnostics.MoveKinds[BassMoveKind.Birth];
        output.WriteLine($"{sampler}: births {births.Accepted}/{births.Attempts}");
    }

    [Fact]
    public void Smc_ResumesBitForBit()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig Config(int maxStages)
        {
            BassConfig c = WithResidual(BassSweepTests.Config(maxBasis: 4));
            return c with { Sampling = c.Sampling with { MasterSeed = 33, Smc = new SmcSpec { Particles = 150, Replicates = 2, MutationSteps = 2, TargetConditionalEss = 0.9, MaxStages = maxStages } } };
        }
        BassFit whole = Bass.Fit(Config(5000), data);
        var run = new BassRunOptions { Checkpoint = true };
        BassFit resumed = Bass.Fit(Config(3), data, run);
        int resumes = 0;
        while (!resumed.Diagnostics.Smc!.ReachedPosterior)
        {
            resumed = Bass.Resume(resumed.Checkpoint!, Config(3), data, run);
            Assert.True(++resumes <= 200, "the resumes never reached b = 1");
        }
        Assert.True(resumes > 0);
        Assert.Equal(BassBitIdentityTests.Hash(whole), BassBitIdentityTests.Hash(resumed));
    }

    // ---- the config and the refusals -----------------------------------------------------------------------------------------------

    [Fact]
    public void Config_AndData_RefuseWhatTheResidualBirthCannotRead()
    {
        BassConfig good = WithResidual(new BassConfig());
        good.Validate();
        BassConfig With(ResidualBirth? r, HistoryBirth? h) => good with { Sampling = good.Sampling with { Births = new BirthMixture { Residual = r, History = h } } };
        With(null, null).Validate();   // empty: the fit's own birth alone, a default turned off
        Assert.Contains("\"Births\":{}", With(null, null).ToJson());
        Assert.Throws<ArgumentException>(() => With(new ResidualBirth(0.0), null).Validate());
        Assert.Throws<ArgumentException>(() => With(new ResidualBirth(1.0), null).Validate());
        Assert.Throws<ArgumentException>(() => With(new ResidualBirth(0.6), new HistoryBirth(Share: 0.4)).Validate());   // no share left
        With(new ResidualBirth(0.5), new HistoryBirth(Share: 0.4)).Validate();

        string json = good.ToJson();
        Assert.Contains("\"Births\":{\"Residual\":{\"Share\":0.3}}", json);
        Assert.Equal(good, BassConfig.FromJson(json));
        Assert.Contains("\"Births\":{\"History\":{\"Share\":0.25,\"Floor\":0.01,\"Refresh\":100}}",
                        (new BassConfig() with { Sampling = new BassSampling { Births = new BirthMixture { History = new HistoryBirth() } } }).ToJson());

        // A direct functional response's residual is a matrix; a fit without continuous predictors has no factor to add.
        Assert.Throws<ArgumentException>(() => BassModel.Create(WithResidual(BassSweepTests.Config()), BassRunTests.FunctionalData()));
        BassData categorical = new([new CategoricalPredictor("z", [0, 1, 2, 0, 1, 2], ["a", "b", "c"])], new ScalarResponse([1, 2, 3, 1.5, 2.5, 3.5]));
        Assert.Throws<ArgumentException>(() => BassModel.Create(WithResidual(new BassConfig()), categorical));
        Assert.Throws<InvalidOperationException>(() =>
        {
            BassModel model = Fixture("continuous");
            new BassBirthMixture(model.Kernel, null, null, new BassResidualBirth(model), Spec).LogBirth(AllTerms(model)[0], BassMixing.Counts(model, Array.Empty<BassTerm>()));
        });
    }

    // ---- the default (the user's decision, 2026-11-09) ------------------------------------------------------------------------------

    // KHAOS's raw-ridge fits of a scalar response take the residual birth at ω_R = ¼ when they declare no mixture, logged as an auto-pick:
    // their draws are a declared residual birth's, bit for bit, and not the fit's own birth's. An empty mixture turns it off; a declared
    // one runs as declared; the modified g-prior, the normalized ridge and hinge fits keep the fit's own birth.
    [Fact]
    public void KhaosRawRidgeFits_TakeTheResidualBirthByDefault_AndAnEmptyMixtureTurnsItOff()
    {
        BassData data = KhaosChainTests.Data();
        BassConfig raw = KhaosChainTests.Config(raw: true) with { Sampling = new BassSampling { Chains = 2, MasterSeed = 61, BurnIn = 200, Samples = 150 } };
        BassConfig With(BirthMixture? births) => raw with { Sampling = raw.Sampling with { Births = births } };

        BassResolved resolved = BassModel.Create(raw, data).Resolved;
        Assert.Equal(new BirthMixture { Residual = BassResolved.DefaultResidualBirth }, resolved.Births);
        Assert.Equal(0.25, BassResolved.DefaultResidualBirth.Share);
        Assert.Contains(resolved.AutoPicks, p => p.Setting == "Births");
        Assert.Null(BassModel.Create(With(new BirthMixture()), data).Resolved.Births);
        Assert.DoesNotContain(BassModel.Create(With(new BirthMixture()), data).Resolved.AutoPicks, p => p.Setting == "Births");
        var historyOnly = new BirthMixture { History = new HistoryBirth() };
        Assert.Equal(historyOnly, BassModel.Create(With(historyOnly), data).Resolved.Births);
        Assert.Null(BassModel.Create(KhaosReachTests.Config("modified g"), data).Resolved.Births);
        Assert.Null(BassModel.Create(KhaosChainTests.Config(raw: false), data).Resolved.Births);
        Assert.Null(BassModel.Create(BassSweepTests.Config(), BassSweepTests.Data()).Resolved.Births);
        Assert.Null(BassModel.Create(BassSweepTests.Config() with { Coefficients = BassCoefficientPrior.RawRidge(TauPrior.Held(0.1)) }, BassSweepTests.Data()).Resolved.Births);

        string byDefault = BassBitIdentityTests.Hash(Bass.Fit(raw, data));
        string declared = BassBitIdentityTests.Hash(Bass.Fit(With(new BirthMixture { Residual = new ResidualBirth(0.25) }), data));
        string own = BassBitIdentityTests.Hash(Bass.Fit(With(new BirthMixture()), data));
        output.WriteLine($"by default {byDefault}, declared {declared}, the fit's own birth {own}");
        Assert.Equal(declared, byDefault);
        Assert.NotEqual(own, byDefault);
        Assert.DoesNotContain("Births", raw.ToJson());
    }
}
