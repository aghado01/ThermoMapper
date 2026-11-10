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
/// E4's history birth (the proposal-kernels brief's slice 1; Nott and Kohn 2005): the birth mixture
/// q = (1 − ω_H)·q_base + ω_H·q_H through F-X4, with q_H the fit's own birth under the weights max(ĥ_u, δ).
/// <list type="bullet">
/// <item><b>The history:</b> its records, refreshes, freeze, restart state and population reading, against counts kept here.</item>
/// <item><b>The kernel's own law:</b> on enumerated term spaces, Σ_t q(t | S) and the empty draw make 1, each ln q against the
/// statement written out here (the base law with each part's variable law replaced by Wallenius' over the weights, or
/// <c>make_weights</c>' passes over them), and the draws against q by χ².</item>
/// <item><b>Detailed balance, exactly:</b> the sweep's flows against the enumerated law, each birth's forward probability from the
/// statement and the production mixture's ratios, with a deliberately skewed ĥ and with one a chain learned.</item>
/// <item><b>The chain's law</b> by enumeration after a learned freeze; <b>Geweke</b> with a skewed ĥ (qualification).</item>
/// <item><b>The engine's reach:</b> the ladder's rungs share the b = 1 history and each samples its law; SMC reads each stage's
/// population and meets the enumerated evidence; the flat histogram learns with its weights; every resume is bit for bit.</item>
/// </list>
/// The fixtures: the hinge enumeration's (two continuous inputs; categorical; functional), three continuous inputs at M_max 1 for
/// Wallenius' law at K = 2, and KHAOS's exact-law fixture (p 4, degree 3, order 2).
/// </summary>
public sealed class HistoryBirthTests(ITestOutputHelper output)
{
    private static readonly HistoryBirth Spec = new(Share: 0.3, Floor: 0.01, Refresh: 7);

    internal static BassConfig WithBirths(BassConfig config, HistoryBirth? spec = null)
        => config with { Sampling = config.Sampling with { Births = new BirthMixture { History = spec ?? Spec } } };

    private static BassModel Fixture(string name, bool births = false)
    {
        (BassConfig config, BassData data) = name switch
        {
            "legendre" => (KhaosChainTests.Config(raw: true), KhaosChainTests.Data()),
            "three inputs" => (new BassConfig
            {
                Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 1 },
                Coefficients = BassCoefficientPrior.G(TauPrior.Held()),
            }, ThreeData()),
            _ => BassLocalStepTests.Fixture(name, "equal", s => s),
        };
        return BassModel.Create(births ? WithBirths(config) : config, data);
    }

    // Three continuous inputs on six cases, each a permutation of the grid: y = 1 + 2x₀ + x₁x₂ + noise.
    private static BassData ThreeData()
    {
        double[] x0 = [0, 1, 2, 3, 4, 5], x1 = [3, 0, 4, 1, 5, 2], x2 = [5, 3, 1, 4, 2, 0];
        for (int i = 0; i < 6; i++) (x0[i], x1[i], x2[i]) = (x0[i] / 5, x1[i] / 5, x2[i] / 5);
        double[] y = [.. Enumerable.Range(0, 6).Select(i => 1 + 2 * x0[i] + x1[i] * x2[i] + 0.3 * Math.Sin(7 * i))];
        return new BassData([new ContinuousPredictor("x0", x0), new ContinuousPredictor("x1", x1), new ContinuousPredictor("x2", x2)], new ScalarResponse(y));
    }

    private static double LogTau(BassModel model) => Math.Log(model.Target.HeldTau!.Value);

    private static BassPart Main(BassModel model) => model.Prior.Has(BassPart.Legendre) ? BassPart.Legendre : BassPart.Continuous;

    // A skewed ĥ that leaves the main part's first variable at zero, so its floor holds: a population of single-term structures
    // over every term that does not use it, weighted 2^(k mod 5).
    internal static BassInclusionHistory Skewed(BassModel model, HistoryBirth spec)
    {
        BassPart main = Main(model);
        List<BassTerm> terms = [.. AllTerms(model).Where(t => !BassTermCounts.Variables(t, main).Contains(0))];
        var history = new BassInclusionHistory(model.Prior, spec);
        history.SetFromPopulation([.. terms.Select(t => (IReadOnlyList<BassTerm>)[t])], [.. terms.Select((_, k) => Math.Pow(2, k % 5))]);
        return history;
    }

    // A history a run's burn-in learns: a chain of the model with the history birth, recorded each sweep, then frozen.
    private static BassInclusionHistory Learned(BassModel births, int sweeps, int seed)
    {
        var chain = new BassChain(births, [], LogTau(births), new Xoshiro256PlusPlus(seed), new Xoshiro256PlusPlus(seed + 1));
        BassInclusionHistory history = chain.Structure.History!;
        history.Learn();
        for (int t = 0; t < sweeps; t++) history.Record(chain.Step().Terms);
        history.Freeze();
        return history;
    }

    private static double[] FloorWeights(BassInclusionHistory history, BassPart part)
        => [.. history.Frequencies(part)!.Select(h => Math.Max(h, history.Spec.Floor))];

    // ---- the statement, written out ------------------------------------------------------------------------------------------------

    // ln q_base(t | S) and ln q_H(t | S). Hinge: the base law is the production proposal's (D01's, checked against Wallenius and the
    // enumeration elsewhere); q_H replaces each part's variable law, uniform at K = 1 or Wallenius' over w₂ + a, by Wallenius' over the
    // floored ĥ. Legendre: the inclusion law Σ_{q₀} J(q₀)P(χ | q₀)/P(valid | q₀), each validity enumerated, from make_weights' η over
    // (a^α + ε) or over the floored ĥ, with the degree and composition laws.
    internal static (double Base, double History) Laws(BassModel model, BassInclusionHistory history, BassTerm t, BassTermCounts counts)
    {
        if (model.Prior.Legendre is null)
        {
            BassTermProposal proposal = model.Proposal;
            double lb = proposal.LogBirth(t, counts), lh = lb;
            foreach (BassPart part in Enum.GetValues<BassPart>().Where(model.Prior.Has))
            {
                int[] set = [.. BassTermCounts.Variables(t, part)];
                if (set.Length == 0) continue;
                int p = model.Prior.Variables(part);
                int[] a = counts.AppearanceCounts(part);
                double baseLaw = set.Length == 1 ? -Math.Log(p) : Wallenius.LogProbability([.. Enumerable.Range(0, p).Select(u => proposal.VariableWeight + a[u])], set);
                lh += Wallenius.LogProbability(FloorWeights(history, part), set) - baseLaw;
            }
            return (lb, lh);
        }
        var legendre = (LegendreTermProposal)model.Kernel;
        LegendreBasis settings = legendre.Settings;
        int vars = legendre.Space.Variables, top = legendre.Space.MaxOrder;
        int[] appearances = counts.AppearanceCounts(BassPart.Legendre);
        double[] raw = FloorWeights(history, BassPart.Legendre);
        bool[] chi = new bool[vars];
        foreach (LegendreFactor f in t.Legendre) chi[f.Variable] = true;
        double Inclusion(Func<int, double[]> etaOf)
        {
            double[] logJ = [.. Enumerable.Range(1, top).Select(q0 => -settings.OrderPenalty * Math.Log(q0))];
            double h = logJ.Sum(Math.Exp), total = 0.0;
            for (int q0 = 1; q0 <= top; q0++)
            {
                double[] eta = etaOf(q0);
                double product = 1.0;
                for (int j = 0; j < vars; j++) product *= chi[j] ? eta[j] : 1.0 - eta[j];
                total += Math.Exp(logJ[q0 - 1]) / h * product / KhaosChainTests.Validity(eta, top);
            }
            return Math.Log(total);
        }
        double[] EtaFromCounts(int q0) { var e = new double[vars]; legendre.InclusionProbabilities(appearances, q0, e); return e; }
        double[] EtaFromHistory(int q0) { var e = new double[vars]; legendre.InclusionProbabilities(raw, q0, e); return e; }
        int q = t.Legendre.Count, d = t.LegendreDegree;
        double rest = KhaosChainTests.LogDegree(settings, q, d) + LegendreTermSpace.LogComposition(d, q);
        return (Inclusion(EtaFromCounts) + rest, Inclusion(EtaFromHistory) + rest);
    }

    private static double Mix(double share, (double Base, double History) laws)
        => Math.Log((1.0 - share) * Math.Exp(laws.Base) + share * Math.Exp(laws.History));

    // ---- the history ---------------------------------------------------------------------------------------------------------------

    [Fact]
    public void History_RecordsRefreshesAndFreezes_AsCountedHere_AndRestores()
    {
        BassModel model = Fixture("categorical, R's measure");
        var spec = new HistoryBirth(Share: 0.4, Floor: 0.05, Refresh: 3);
        var history = new BassInclusionHistory(model.Prior, spec);
        List<BassTerm> terms = AllTerms(model);
        var rng = new Xoshiro256PlusPlus(5);
        var structures = new List<BassTerm[]>();
        for (int k = 0; k < 11; k++) structures.Add([.. Enumerable.Range(0, rng.NextInt(3)).Select(_ => terms[rng.NextInt(terms.Count)])]);

        history.Record(structures[0]);   // not learning: nothing
        Assert.Equal(0, history.Recorded);
        history.Learn();
        for (int k = 0; k < 11; k++)
        {
            history.Record(structures[k]);
            int recorded = k + 1, refreshed = recorded / 3 * 3;
            Assert.Equal(recorded, history.Recorded);
            if (refreshed == 0) { Assert.Null(history.Weights); Assert.Null(history.Reading()); continue; }
            foreach (BassPart part in Enum.GetValues<BassPart>())
            {
                if (!model.Prior.Has(part)) { Assert.Null(history.Frequencies(part)); continue; }
                double[] expected = [.. Enumerable.Range(0, model.Prior.Variables(part)).Select(v =>
                    structures.Take(refreshed).Count(s => s.Any(t => BassTermCounts.Variables(t, part).Contains(v))) / (double)refreshed)];
                Assert.Equal(expected, history.Frequencies(part)!);
                Assert.Equal(expected.Select(h => Math.Max(h, 0.05)).ToArray(), history.Weights![part].ToArray());
            }
        }
        // A restart mid-learning continues as the original does.
        BassInclusionHistory copy = history.Copy();
        history.Freeze();
        copy.Freeze();
        Assert.True(history.Frozen && !history.Learning);
        Assert.Equal(11, history.Recorded);
        foreach (BassPart part in Enum.GetValues<BassPart>().Where(model.Prior.Has))
        {
            Assert.Equal(history.Frequencies(part), copy.Frequencies(part));
            double[] all = [.. Enumerable.Range(0, model.Prior.Variables(part)).Select(v =>
                structures.Count(s => s.Any(t => BassTermCounts.Variables(t, part).Contains(v))) / 11.0)];
            Assert.Equal(all, history.Frequencies(part)!);   // the freeze refreshes from every record
        }
        BassInclusionReading reading = history.Reading()!;
        Assert.Equal(11, reading.Recorded);
        Assert.Equal(history.Frequencies(BassPart.Categorical), reading[BassPart.Categorical]);
        history.Record(structures[0]);
        Assert.Equal(11, history.Recorded);   // frozen: nothing
        Assert.Throws<InvalidOperationException>(history.Learn);
        Assert.Throws<InvalidOperationException>(history.Freeze);

        // A freeze before any record leaves the base law; a population sets the weighted share.
        var empty = new BassInclusionHistory(model.Prior, spec);
        empty.Learn();
        empty.Freeze();
        Assert.Null(empty.Weights);
        var population = new BassInclusionHistory(model.Prior, spec);
        double[] w = [0.5, 0.0, 1.5, 2.0];
        population.SetFromPopulation([.. structures.Take(4)], w);
        foreach (BassPart part in Enum.GetValues<BassPart>().Where(model.Prior.Has))
        {
            double[] expected = [.. Enumerable.Range(0, model.Prior.Variables(part)).Select(v =>
                Enumerable.Range(0, 4).Where(i => structures[i].Any(t => BassTermCounts.Variables(t, part).Contains(v))).Sum(i => w[i]) / 4.0)];
            Assert.Equal(expected, population.Frequencies(part)!);
        }
        Assert.Throws<ArgumentException>(() => population.SetFromPopulation([.. structures.Take(2)], [0.0, 0.0]));
    }

    // ---- the kernel's own law ------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData("three inputs", 20261130)]
    [InlineData("categorical, R's measure", 20261131)]
    [InlineData("functional", 20261132)]
    [InlineData("legendre", 20261133)]
    public void Law_SumsToOne_IsTheStatement_AndTheDrawsFollowIt(string fixture, int seed)
    {
        BassModel model = Fixture(fixture);
        BassInclusionHistory history = Skewed(model, Spec);
        Assert.Equal(0.0, history.Frequencies(Main(model))![0]);
        var mixture = new BassBirthMixture(model.Kernel, Spec, history);
        var weighted = (IBassWeightedBirth)model.Kernel;
        List<BassTerm> terms = AllTerms(model);
        var structures = new List<BassTerm[]> { Array.Empty<BassTerm>(), new[] { terms[0] }, new[] { terms[^1] } };
        if (model.MaxBasis >= 2) structures.Add([terms[1], terms[^2]]);
        double worst = 0.0;
        foreach (BassTerm[] s in structures)
        {
            BassTermCounts counts = BassMixing.Counts(model, s);
            double empty = model.Prior.Legendre is null ? Math.Exp(model.Proposal.LogEmptyBirth(counts)) : 0.0;
            double mixed = empty, alone = empty;
            foreach (BassTerm t in terms)
            {
                double lq = mixture.LogBirth(t, counts);
                mixed += Math.Exp(lq);
                alone += Math.Exp(weighted.LogBirth(t, counts, 0, history.Weights!));
                worst = Math.Max(worst, Math.Abs(lq - Mix(Spec.Share, Laws(model, history, t, counts))));
            }
            output.WriteLine($"{fixture}, M = {s.Length}: Σq {mixed:R} (the history component alone {alone:R}); empty draw {empty:G4}");
            Assert.Equal(1.0, mixed, 1e-10);
            Assert.Equal(1.0, alone, 1e-10);
        }
        output.WriteLine($"{fixture}: {terms.Count} terms; worst |ln q − the statement| {worst:E2}");
        Assert.True(worst <= 1e-10, $"ln q sits {worst:E2} from the statement");

        // The draws against q at the last structure, by χ², cells below five expected pooled.
        BassTermCounts at = BassMixing.Counts(model, structures[^1]);
        const int draws = 200_000;
        var rng = new Xoshiro256PlusPlus(seed);
        var tally = new Dictionary<BassTerm, int>();
        int nulls = 0;
        for (int k = 0; k < draws; k++)
        {
            if (mixture.DrawBirth(rng, at) is { } t) tally[t] = tally.GetValueOrDefault(t) + 1;
            else nulls++;
        }
        var cells = new List<(double Expected, int Observed)> { (draws * (model.Prior.Legendre is null ? Math.Exp(model.Proposal.LogEmptyBirth(at)) : 0.0), nulls) };
        cells.AddRange(terms.Select(t => (draws * Math.Exp(mixture.LogBirth(t, at)), tally.GetValueOrDefault(t))));
        Assert.Equal(draws, cells.Sum(c => c.Observed));
        Assert.Equal(0, cells.Where(c => c.Expected == 0.0).Sum(c => c.Observed));   // nothing q excludes is drawn
        double pooledExpected = cells.Where(c => c.Expected < 5).Sum(c => c.Expected);
        int pooledObserved = cells.Where(c => c.Expected < 5).Sum(c => c.Observed);
        List<(double Expected, int Observed)> kept = [.. cells.Where(c => c.Expected >= 5)];
        if (pooledExpected > 0) kept.Add((pooledExpected, pooledObserved));
        double chi2 = kept.Sum(c => (c.Observed - c.Expected) * (c.Observed - c.Expected) / c.Expected);
        int df = kept.Count - 1;
        double z = (chi2 - df) / Math.Sqrt(2.0 * df);
        output.WriteLine($"{fixture}: {draws} draws over {kept.Count} cells, χ² {chi2:F1} on {df} df, z {z:F2}");
        Assert.True(z < 4.0, $"the draws sit z = {z:F2} from q");
    }

    // ---- detailed balance, exactly ------------------------------------------------------------------------------------------------

    // Births from the statement (the mixture written out) with the production mixture's ratio; deaths uniform with its ratio; the
    // changes the base kernel's, from their own enumerations.
    private static IEnumerable<Enumerated> Mixed(Space space, BassBirthMixture mixture, BassInclusionHistory history, string key, State s)
    {
        BassModel model = space.Model;
        int m = s.Terms.Length;
        switch (key)
        {
            case "birth":
                foreach (BassTerm t in space.Terms)
                {
                    if (s.Terms.Contains(t)) continue;
                    yield return new Enumerated(BirthDeath, [.. s.Terms, t], Mix(Spec.Share, Laws(model, history, t, s.Counts)), mixture.LogBirthRatio(t, s.Counts));
                }
                break;
            case "death":
                for (int k = 0; k < m; k++)
                    yield return new Enumerated(BirthDeath, [.. s.Terms.Where((_, i) => i != k)], -Math.Log(m), mixture.LogDeathRatio(s.Terms[k], s.Counts));
                break;
            default:
                foreach (Enumerated e in model.Prior.Legendre is null ? All(space, key, s, null) : KhaosChainTests.Khaos(space, key, s)) yield return e;
                break;
        }
    }

    [Theory]
    [InlineData("continuous", false)]
    [InlineData("continuous", true)]
    [InlineData("categorical, R's measure", false)]
    [InlineData("functional", false)]
    [InlineData("three inputs", false)]
    [InlineData("legendre", false)]
    [InlineData("legendre", true)]
    public void DetailedBalance_Enumerated(string fixture, bool learned)
    {
        BassModel model = Fixture(fixture);
        BassInclusionHistory history = learned ? Learned(Fixture(fixture, births: true), 3_000, 20261140) : Skewed(model, Spec);
        var mixture = new BassBirthMixture(model.Kernel, Spec, history);
        Space space = model.Prior.Legendre is null ? Enumerate(model) : KhaosChainTests.Enumerate(model);
        var flows = Kernel(space, (key, s) => Mixed(space, mixture, history, key, s));
        (Dictionary<string, int> transitions, double worst) = Check(space, flows);
        output.WriteLine($"{fixture}, {(learned ? $"learned over {history.Recorded} sweeps" : "skewed")} ĥ " +
                         $"[{string.Join(", ", history.Frequencies(Main(model))!.Select(h => h.ToString("F3")))}]: |F| = {space.States.Count}; " +
                         $"{string.Join(", ", transitions.Select(t => $"{t.Key} {t.Value}"))}; worst gap {worst:E2}");
        Assert.True(transitions[BirthDeath] > 0);
    }

    // The check's power: a mixture whose ratio reads another share than its draws breaks detailed balance.
    [Fact]
    public void DetailedBalance_Enumerated_RefusesAWrongShare()
    {
        BassModel model = Fixture("continuous");
        BassInclusionHistory history = Skewed(model, Spec);
        var wrong = new BassBirthMixture(model.Kernel, Spec with { Share = 0.5 }, history);
        Space space = Enumerate(model);
        var flows = Kernel(space, (key, s) => Mixed(space, wrong, history, key, s));
        Assert.ThrowsAny<Xunit.Sdk.XunitException>(() => Check(space, flows));
    }

    // ---- the chain's law after a learned freeze ------------------------------------------------------------------------------------

    [Theory]
    [InlineData("continuous", 20261150)]
    [InlineData("legendre", 20261151)]
    public void Chain_MeetsTheExactLaw_AfterALearnedFreeze(string fixture, int seed)
    {
        BassModel model = Fixture(fixture, births: true);
        Space space = model.Prior.Legendre is null ? Enumerate(model) : KhaosChainTests.Enumerate(model);
        double[] pi = space.Posterior();
        var chain = new BassChain(model, [], LogTau(model), new Xoshiro256PlusPlus(seed), new Xoshiro256PlusPlus(seed + 1));
        BassInclusionHistory history = chain.Structure.History!;
        history.Learn();
        for (int t = 0; t < 5_000; t++) history.Record(chain.Step().Terms);
        history.Freeze();
        Assert.Equal(5_000, history.Recorded);
        const int steps = 400_000, batches = 40;
        var visits = new int[steps];
        for (int t = 0; t < steps; t++)
        {
            BassState s = chain.Step();
            history.Record(s.Terms);   // frozen: nothing
            int i = space.Find(s.Terms);
            Assert.True(i >= 0, $"the chain holds a structure outside the enumeration: {string.Join(" + ", s.Terms)}");
            visits[t] = i;
        }
        Assert.Equal(5_000, history.Recorded);
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
            // A reading the chain never varied on, such as M = 0 here: it agrees when the law puts nearly all or none of its mass there.
            double z = se > 0.0 ? (mean - exact) / se : (Math.Abs(mean - exact) < 1e-3 ? 0.0 : double.PositiveInfinity);
            worst = Math.Max(worst, Math.Abs(z));
            output.WriteLine($"{name,-40} chain {mean:F4}  exact {exact:G6}  z {z,6:F2}");
        }
        output.WriteLine($"{fixture}: ĥ [{string.Join(", ", history.Frequencies(main)!.Select(h => h.ToString("F3")))}], {steps} steps: worst |z| {worst:F2}; " +
                         $"births {chain.Structure.Kinds.Entries.Where(e => e.Kind == BassMoveKind.Birth).Select(e => $"{e.Accepted}/{e.Attempts}").FirstOrDefault()}");
        Assert.True(worst < 4.5, $"the chain's worst |z| {worst:F2} against the exact law");
    }

    // ---- Geweke (qualification) ----------------------------------------------------------------------------------------------------

    private const int GewekeN = 12;

    // Proper priors that read no response: the zero centre, σ² ~ IG(3, 2), τ held. Hinge: three inputs at order 2, so the history's
    // weights enter Wallienius' law at K = 1 and 2; Legendre: four inputs, the swap on.
    private static BassConfig GewekeConfig(bool legendre) => WithBirths(legendre
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
        });

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
    [InlineData(false, 20261160)]
    [InlineData(true, 20261161)]
    public void Geweke_WithASkewedHistory(bool legendre, int seed)
    {
        if (!BassQualification.Enabled) return;   // heavy
        const int draws = 200_000, steps = 300_000, batches = 50;
        BassConfig config = GewekeConfig(legendre);
        var rng = new Xoshiro256PlusPlus(seed);
        BassModel prior = BassModel.Create(config, GewekeData(legendre, [.. Enumerable.Range(0, GewekeN).Select(i => Math.Sin(i))]));
        BassInclusionHistory skewed = Skewed(prior, Spec);
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
            chain.Structure.ShareHistory(skewed);
            chain.Step();
            BassSweepFit fit = chain.Structure.DrawSweepFit(chain.Current, rng, model.Target);
            successive[t] = Functions(model, chain.Current, fit.Sigma2, fit.Fitted, y);
            theta = (chain.Current, fit.Sigma2, fit.Fitted, y);
            chain.Structure.DesignAt(chain.Current);
            foreach (BassKindCount e in chain.Structure.Kinds.Entries.Where(e => e.Kind == BassMoveKind.Birth)) { births += e.Attempts; accepted += e.Accepted; }
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

    // ---- the engine's reach --------------------------------------------------------------------------------------------------------

    [Theory]
    [InlineData("hinge", new[] { 1.0, 0.5, 0.2 })]
    [InlineData("legendre", new[] { 1.0, 0.4, 0.1, 0.0 })]
    public void Ladder_SharesTheColdHistory_AndEachRungSamplesItsTemperedLaw(string basis, double[] rungs)
    {
        bool reference = rungs[^1] == 0.0;
        var spec = new TemperingSpec { Levels = rungs.Length, HotEnd = reference ? rungs[^2] : rungs[^1], ReferenceRung = reference };
        BassConfig config = basis == "legendre" ? KhaosReachTests.Config("raw ridge") : BassTemperingTests.Config(spec);
        BassData data = basis == "legendre" ? KhaosChainTests.Data() : BassSweepTests.Data();
        BassModel model = BassModel.Create(WithBirths(config with { Sampling = config.Sampling with { Tempering = spec } }), data);

        // Every rung reads the b = 1 rung's history, at the start and after a self-sizing ladder's resize.
        BassState start = new BassPriorSampler(model).DrawState(new Xoshiro256PlusPlus(3)).State;
        var sizing = new TemperingSpec { ReferenceRung = reference };
        var ladder = new BassLadder(BassModel.Create(model.Config with { Sampling = model.Config.Sampling with { Tempering = sizing } }, data), sizing,
                                    start.Terms, start.LogTau, seed: 4, start.StartLogRate);
        Assert.True(ladder.SizesItself);
        BassInclusionHistory shared = ladder.Structure.History!;
        Assert.All(ladder.Structures, s => Assert.Same(shared, s.History));
        int before = ladder.Structures.Count;
        for (int round = 1; round <= 64; round *= 2)
        {
            for (int t = 0; t < round; t++) ladder.Step();
            ladder.Resize();
            Assert.All(ladder.Structures, s => Assert.Same(shared, s.History));
        }
        output.WriteLine($"{basis}: the self-sizing ladder went from {before} to {ladder.Structures.Count} rungs, every one on the b = 1 history");

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
            output.WriteLine($"{basis}, b = {rungs[r]:G3}, E[ln τ]: {u:F4} against {exact[r].LogTauMean:F4} (± {uError:F4})");
            Assert.True(Math.Abs(u - exact[r].LogTauMean) <= 5 * uError + 0.01, $"b = {rungs[r]}, E[ln τ] {u:F4} against {exact[r].LogTauMean:F4} (± {uError:F4})");
        }
    }

    [Theory]
    [InlineData("raw ridge")]
    [InlineData("modified g")]
    public void Smc_ReadsEachStagesPopulation_AndMeetsTheEnumeratedEvidence(string prior)
    {
        BassData data = KhaosChainTests.Data();
        BassConfig config = WithBirths(KhaosReachTests.Config(prior));
        config = config with { Sampling = config.Sampling with { MasterSeed = 3, ReportEvidence = true, Smc = new SmcSpec { Particles = 2000, Replicates = 4, MutationSteps = 4 } } };
        BassModel model = BassModel.Create(config, data);
        RungLaw exact = Exact(model, 1.0)[0];

        // The kernel's stage weights are the population's weighted share.
        var kernel = new BassSmcKernel(model);
        var particles = new BassParticles(model);
        var rng = new Xoshiro256PlusPlus(9);
        BassParticle[] population = [.. Enumerable.Range(0, 50).Select(_ => particles.Draw(rng))];
        double[] w = [.. Enumerable.Range(0, 50).Select(i => (i % 7) / 21.0)];
        kernel.Prepare(default, population, w);
        double[] expected = [.. Enumerable.Range(0, 4).Select(v =>
            Enumerable.Range(0, 50).Where(i => population[i].State.Terms.Any(t => t.Legendre.Any(f => f.Variable == v))).Sum(i => w[i]) / w.Sum())];
        IReadOnlyList<double> read = kernel.History!.Frequencies(BassPart.Legendre)!;
        for (int v = 0; v < 4; v++) Assert.Equal(expected[v], read[v], 1e-14);

        BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 4 });
        BassEvidence evidence = fit.Evidence!;
        output.WriteLine($"{prior}, SMC with the history birth: ln p(y) {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against {exact.LogZ:F4}");
        Assert.True(Math.Abs(evidence.LogEvidence - exact.LogZ) <= 4 * evidence.MonteCarloError + 0.01,
                    $"{prior}: {evidence.LogEvidence:F4} ± {evidence.MonteCarloError:F4} against {exact.LogZ:F4}");
        Assert.Empty(fit.Diagnostics.Smc!.Flags);
    }

    [Theory]
    [InlineData("raw ridge", 6)]
    public void FlatHistogram_LearnsWithItsWeights_AndMatchesTheEnumeration(string prior, int nodes)
        => BassFlatHistogramTests.MatchesTheEnumeration(BassModel.Create(WithBirths(KhaosReachTests.Config(prior)), KhaosChainTests.Data()), 0.0,
                                                        $"KHAOS with the history birth, {prior}", output, nodes);

    // ---- a run: the readout, resume, and the draws without the mixture --------------------------------------------------------------

    private static BassConfig RunConfig(int samples, TemperingSpec? tempering = null)
    {
        BassConfig c = WithBirths(BassRunTests.Config(chains: 2, burn: 250, samples: samples, seed: 21));
        return c with { Sampling = c.Sampling with { Tempering = tempering } };
    }

    [Fact]
    public void Fit_ReportsEachChainsFrozenReading()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassFit fit = Bass.Fit(RunConfig(100), data);
        Assert.All(fit.Diagnostics.Chains, account =>
        {
            BassInclusionReading reading = account.Inclusion!;
            Assert.Equal(250, reading.Recorded);
            Assert.All(reading[BassPart.Continuous]!, h => Assert.InRange(h, 0.0, 1.0));
            Assert.Null(reading[BassPart.Categorical]);
        });
        output.WriteLine(string.Join("; ", fit.Diagnostics.Chains.Select((a, c) => $"chain {c}: ĥ [{string.Join(", ", a.Inclusion![BassPart.Continuous]!.Select(h => h.ToString("F3")))}]")));
        Assert.All(Bass.Fit(BassRunTests.Config(chains: 2, burn: 250, samples: 100, seed: 21), data).Diagnostics.Chains, a => Assert.Null(a.Inclusion));
    }

    [Theory]
    [InlineData("chain")]
    [InlineData("ladder")]
    public void Resume_IsBitForBit(string sampler)
    {
        BassData data = BassSweepTests.Data(large: true);
        TemperingSpec? tempering = sampler == "ladder" ? new TemperingSpec { Levels = 3, ReferenceRung = true } : null;
        BassFit whole = Bass.Fit(RunConfig(300, tempering), data);
        BassConfig half = RunConfig(150, tempering);
        BassFit first = Bass.Fit(half, data, new BassRunOptions { Checkpoint = true });
        BassFit resumed = Bass.Resume(first.Checkpoint!, half, data);
        BassRunTests.AssertSameFit(whole, resumed);
        Assert.Equal(BassBitIdentityTests.Hash(whole), BassBitIdentityTests.Hash(resumed));
        for (int c = 0; c < whole.Diagnostics.Chains.Count; c++)
            Assert.Equal(whole.Diagnostics.Chains[c].Inclusion!.ByPart[0], resumed.Diagnostics.Chains[c].Inclusion!.ByPart[0]);
    }

    [Fact]
    public void Smc_ResumesBitForBit()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig Config(int maxStages)
        {
            BassConfig c = WithBirths(BassSweepTests.Config(maxBasis: 4));
            return c with { Sampling = c.Sampling with { MasterSeed = 23, Smc = new SmcSpec { Particles = 150, Replicates = 2, MutationSteps = 2, TargetConditionalEss = 0.9, MaxStages = maxStages } } };
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
        Assert.Equal(whole.Evidence!.LogEvidence, resumed.Evidence!.LogEvidence);
    }

    [Fact]
    public void FlatHistogram_ResumesBitForBit_AcrossLearningAndProduction()
    {
        BassData data = BassSweepTests.Data(large: true);
        BassConfig plain = WithBirths(new BassConfig
        {
            Structure = new BassStructureSpec { MaxInteraction = 2, MaxBasis = 8 },
            Sampling = new BassSampling { MasterSeed = 25, Thin = 2 },
        });
        BassModel model = BassModel.Create(plain, data);
        BassFlatHistogramSpec spec = new() { Run = BassFlatHistogramTests.DeclaredRun(BassFlatHistogramTests.RoughShape(model, 8), 300, 200, chains: 3) };
        BassFit whole = Bass.Fit(BassFlatHistogramTests.WithFlatHistogram(plain, spec), data);
        Assert.All(whole.Diagnostics.Chains, a => Assert.Equal(600, a.Inclusion!.Recorded));   // 300 learning iterations of two sweeps
        BassConfig cut = BassFlatHistogramTests.WithFlatHistogram(plain, spec with { MaxIterations = 170 });
        var run = new BassRunOptions { Checkpoint = true };
        BassFit first = Bass.Fit(cut, data, run);
        BassFit second = Bass.Resume(first.Checkpoint!, cut, data, run);
        BassFit third = Bass.Resume(second.Checkpoint!, cut, data, run);
        Assert.True(third.Diagnostics.FlatHistogram!.Complete);
        Assert.Equal(BassBitIdentityTests.Hash(whole), BassBitIdentityTests.Hash(third));
        for (int c = 0; c < 3; c++)
            Assert.Equal(whole.Diagnostics.Chains[c].Inclusion!.ByPart[0], third.Diagnostics.Chains[c].Inclusion!.ByPart[0]);
    }

    // ---- the config ----------------------------------------------------------------------------------------------------------------

    [Fact]
    public void Config_RefusesWhatTheMixtureCannotRun_AndAConfigWithoutItStandsAsItWas()
    {
        BassConfig good = WithBirths(new BassConfig());
        good.Validate();
        BassConfig With(HistoryBirth h) => good with { Sampling = good.Sampling with { Births = new BirthMixture { History = h } } };
        (good with { Sampling = good.Sampling with { Births = new BirthMixture() } }).Validate();   // empty: the fit's own birth alone
        Assert.Throws<ArgumentException>(() => With(Spec with { Share = 0.0 }).Validate());
        Assert.Throws<ArgumentException>(() => With(Spec with { Share = 1.0 }).Validate());
        Assert.Throws<ArgumentException>(() => With(Spec with { Floor = 0.0 }).Validate());
        Assert.Throws<ArgumentException>(() => With(Spec with { Floor = 1.5 }).Validate());
        Assert.Throws<ArgumentException>(() => With(Spec with { Refresh = 0 }).Validate());
        With(Spec with { Floor = 1.0, Share = 0.999 }).Validate();

        Assert.DoesNotContain("Births", new BassConfig().ToJson());
        string json = good.ToJson();
        Assert.Contains("\"Births\":{\"History\":{\"Share\":0.3,\"Floor\":0.01,\"Refresh\":7}}", json);
        Assert.Equal(good, BassConfig.FromJson(json));
        BassData data = BassSweepTests.Data();
        Assert.NotEqual(Bass.Fingerprint(new BassConfig(), data), Bass.Fingerprint(good, data));
    }
}
