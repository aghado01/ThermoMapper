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
/// Whether the readout fixture's chains hold the posterior or a shared trap (the slice 1 check, behind
/// <see cref="BassQualification"/>). On slice 0's Ishigami readout fixture under KHAOS's raw ridge (degree 8, order 3, khaos's
/// τ² = 10⁵, h₁ = 4, h₂ = 40/n, σ² Jeffreys, the zero centre), in two boxes: the raw inputs with [0, 1] declared, where both
/// samplers read M about 10.1 and σ² about 0.23, and the inputs rescaled to [0, 1] exactly (no domain, the training range's box),
/// where both read M about 6.2 and σ² about 0.95. In each, the collapsed log posterior (the log marginal plus the structure prior,
/// λ integrated) of a greedy forward selection and of its best-improvement descent by deaths, swaps and births; in the rescaled box
/// also the raw structure's image under the affine map between the boxes, the closure of its multi-indices, which spans the raw
/// columns exactly, and that image's descent. Against them, chains of the readout's length from M = 0 and from the richest
/// state, traced in log π, M and E[σ² | structure]. The reading (2026-11-08): the box effect is the posterior's, not a trap. The
/// image fits better (E[σ² | S] 0.16) but scores 133 nats below the sparse mode, each term costing about 14 nats under the raw
/// ridge (ln|𝒜| ≈ 6 and ½ln(nτ²) ≈ 8); its descent returns to the mode, and the chains agree from either start. Two boxes of one
/// width, centred on ½ and shifted by 0.01, separate Ishigami's parity about ½ as the cause.
/// </summary>
public sealed class KhaosRichStartTests(ITestOutputHelper output)
{
    private const int Inputs = 4, Cases = 120, Burn = 10000, Kept = 10000;

    // Slice 0's readout fixture (KhaosReadoutParityTests.Data, n 120, its seed): Ishigami's function of the first three inputs on
    // [0, 1]³ mapped to [−π, π]³, a fourth inert, noise 0.5·N(0, 1).
    private static BassData Fixture(bool declared) => Fixture(declared ? new InputDomain(0.0, 1.0) : null, rescale: !declared);

    // The fixture's raw inputs in the box `domain`, or rescaled onto [0, 1] exactly.
    private static BassData Fixture(InputDomain? domain, bool rescale)
    {
        var rng = new Xoshiro256PlusPlus(20261107);
        var x = new double[Inputs][];
        for (int j = 0; j < Inputs; j++) x[j] = new double[Cases];
        var y = new double[Cases];
        for (int i = 0; i < Cases; i++)
        {
            for (int j = 0; j < Inputs; j++) x[j][i] = rng.NextDouble();
            double z1 = Math.PI * (2 * x[0][i] - 1), z2 = Math.PI * (2 * x[1][i] - 1), z3 = Math.PI * (2 * x[2][i] - 1);
            y[i] = Math.Sin(z1) + 7 * Math.Pow(Math.Sin(z2), 2) + 0.1 * Math.Pow(z3, 4) * Math.Sin(z1) + 0.5 * Normal.SampleStandard(rng);
        }
        if (rescale)
            for (int j = 0; j < Inputs; j++)
            {
                double min = x[j].Min(), max = x[j].Max();
                for (int i = 0; i < Cases; i++) x[j][i] = (x[j][i] - min) / (max - min);
            }
        return new BassData([.. Enumerable.Range(0, Inputs).Select(j => new ContinuousPredictor($"x{j}", x[j], domain))], new ScalarResponse(y));
    }

    private static BassConfig Config() => new()
    {
        Structure = new BassStructureSpec { MaxInteraction = 3, Legendre = new LegendreBasis { MaxDegree = 8 } },
        Count = new BassCountPrior(4, 40.0 / Cases),
        Coefficients = BassCoefficientPrior.RawRidge(TauPrior.Held(1e-5)),
        Intercept = GPriorCentre.Zero,
        Sigma2 = Sigma2Prior.Jeffreys,
        Sampling = new BassSampling { Chains = 1, MasterSeed = 1 },
    };

    private sealed record Trace(string Start, int Seed, double[] LogPosterior, double[] M, double[] Sigma2);

    // One model's structure scores at its held τ.
    private sealed class Scorer(BassModel model)
    {
        public readonly double Tau = model.Target.HeldTau!.Value;
        public readonly BassTerm[] Space = [.. model.Prior.Legendre!.Terms()];
        public BassModel Model => model;
        public int Calls;

        public double Score(IReadOnlyList<BassTerm> terms)
        {
            Calls++;
            return BassMixing.Loaded(model, terms, Math.Log(Tau), out _) is { } state ? BassMixing.LogPosterior(model, state) : double.NegativeInfinity;
        }

        public double Sigma2(BassState state) => model.Target.Marginal(state.Fit, Tau).Sigma2Posterior(model.Target.Sigma2At(state)).Mean;

        public double Sigma2(IReadOnlyList<BassTerm> terms)
            => BassMixing.Loaded(model, terms, Math.Log(Tau), out _) is { } state ? Sigma2(state) : double.NaN;

        // Greedy forward selection from `start`: the best single birth while one raises log π.
        public (List<BassTerm> Terms, double Score) Forward(IEnumerable<BassTerm> start)
        {
            List<BassTerm> terms = [.. start];
            double best = Score(terms);
            while (true)
            {
                (BassTerm add, double score) = Space.Where(t => !terms.Contains(t)).Select(t => (t, Score([.. terms, t]))).MaxBy(c => c.Item2);
                if (!(score > best)) return (terms, best);
                terms.Add(add);
                best = score;
            }
        }

        // Best-improvement descent from `start`: deaths alone while one raises log π (or leaves F's outside), then passes over every
        // death, swap and birth until none raises it.
        public (List<BassTerm> Terms, double Score, int Passes) Descend(IEnumerable<BassTerm> start)
        {
            List<BassTerm> terms = [.. start];
            double best = Score(terms);
            int passes = 0;
            for (bool improved = true; improved && terms.Count > 0; passes++)
            {
                improved = false;
                (int at, double score) = Enumerable.Range(0, terms.Count).Select(i => (i, Score([.. terms.Where((_, k) => k != i)]))).MaxBy(c => c.Item2);
                if (score > best || (double.IsNegativeInfinity(best) && double.IsFinite(score)) || (double.IsNegativeInfinity(best) && double.IsNegativeInfinity(score)))
                {
                    terms.RemoveAt(at);
                    best = score;
                    improved = true;
                }
            }
            for (bool improved = true; improved && passes < 400; passes++)
            {
                improved = false;
                List<BassTerm>? next = null;
                double target = best;
                for (int i = 0; i < terms.Count; i++)
                {
                    List<BassTerm> without = [.. terms.Where((_, k) => k != i)];
                    double s = Score(without);
                    if (s > target) (next, target) = (without, s);
                    foreach (BassTerm t in Space)
                    {
                        if (terms.Contains(t)) continue;
                        List<BassTerm> swapped = [.. without, t];
                        double w = Score(swapped);
                        if (w > target) (next, target) = (swapped, w);
                    }
                }
                foreach (BassTerm t in Space)
                {
                    if (terms.Contains(t)) continue;
                    List<BassTerm> grown = [.. terms, t];
                    double g = Score(grown);
                    if (g > target) (next, target) = (grown, g);
                }
                if (next is not null)
                {
                    terms = next;
                    best = target;
                    improved = true;
                }
            }
            return (terms, best, passes);
        }
    }

    // The span of `terms` under any affine map of each input: ψ_α(ax + b) is a combination of the ψ_β with β ≤ α, so a term's
    // image lies in the span of every multi-index below its own, the constant left out.
    private static List<BassTerm> Closure(IEnumerable<BassTerm> terms)
    {
        var closed = new HashSet<BassTerm>();
        foreach (BassTerm t in terms)
        {
            LegendreFactor[] f = [.. t.Legendre];
            int total = f.Aggregate(1, (acc, x) => acc * (x.Degree + 1));
            for (int code = 1; code < total; code++)
            {
                var factors = new List<LegendreFactor>();
                int c = code;
                foreach (LegendreFactor x in f)
                {
                    int b = c % (x.Degree + 1);
                    c /= x.Degree + 1;
                    if (b > 0) factors.Add(new LegendreFactor(x.Variable, b));
                }
                closed.Add(new BassTerm(legendre: factors));
            }
        }
        return [.. closed];
    }

    private static string Show(IEnumerable<BassTerm> terms)
        => string.Join(" + ", terms.Select(t => string.Concat(t.Legendre.Select(f => $"ψ{f.Degree}(x{f.Variable})"))).Order(StringComparer.Ordinal));

    // Chains of the readout's length from each start, traced in log π, M and E[σ² | S].
    private List<Trace> Chains(Scorer scorer, (string Name, IReadOnlyList<BassTerm> Terms)[] starts)
    {
        var traces = new List<Trace>();
        double logTau = Math.Log(scorer.Tau);
        foreach ((string name, IReadOnlyList<BassTerm> terms) in starts)
        {
            for (int seed = 0; seed < 4; seed++)
            {
                var chain = new BassChain(scorer.Model, terms, logTau, new Xoshiro256PlusPlus(1000 + seed), new Xoshiro256PlusPlus(2000 + seed));
                var early = new double[Burn];
                for (int s = 0; s < Burn; s++) early[s] = BassMixing.LogPosterior(scorer.Model, chain.Step());
                var lp = new double[Kept];
                var m = new double[Kept];
                var s2 = new double[Kept];
                for (int s = 0; s < Kept; s++)
                {
                    BassState st = chain.Step();
                    lp[s] = BassMixing.LogPosterior(scorer.Model, st);
                    m[s] = st.Terms.Count;
                    s2[s] = scorer.Sigma2(st);
                }
                traces.Add(new Trace(name, seed, lp, m, s2));
                output.WriteLine($"  from {name}, seed {seed}: log π during burn-in at steps 1, 100, 1,000, 10,000: {early[0]:F1}, {early[99]:F1}, {early[999]:F1}, {early[^1]:F1}; " +
                                 $"kept: median {BassMixing.Quantile(lp, 0.5):F2}, max {lp.Max():F2}; M {m.Average():F2} ({m.Min()}–{m.Max()}); E[σ²] {s2.Average():F4}");
            }
            Trace[] group = [.. traces.Where(t => t.Start == name)];
            double[] all = [.. group.SelectMany(t => t.LogPosterior)];
            output.WriteLine($"from {name}, 4 chains: log π median {BassMixing.Quantile(all, 0.5):F2}, 99% {BassMixing.Quantile(all, 0.99):F2}; " +
                             $"M {group.Average(t => t.M.Average()):F2}; E[σ²] {group.Average(t => t.Sigma2.Average()):F4}");
        }
        return traces;
    }

    [Fact]
    public void RichStart_AgainstTheChains()
    {
        if (!BassQualification.Enabled) return;   // heavy
        var clock = System.Diagnostics.Stopwatch.StartNew();

        // The raw inputs with [0, 1] declared: greedy, then its descent.
        var raw = new Scorer(BassModel.Create(Config(), Fixture(declared: true)));
        (List<BassTerm> rawGreedy, double rawForward) = raw.Forward([]);
        (List<BassTerm> rawLocal, double rawScore, int rawPasses) = raw.Descend(rawGreedy);
        output.WriteLine($"raw inputs, [0, 1] declared; |𝒜| {raw.Space.Length}: greedy to M {rawGreedy.Count} at log π {rawForward:F2}; descent ({rawPasses} passes) to M {rawLocal.Count} " +
                         $"at {rawScore:F2}, E[σ² | S] {raw.Sigma2(rawLocal):F4}: {Show(rawLocal)}");
        List<Trace> rawTraces = Chains(raw, [("M = 0", Array.Empty<BassTerm>()), ("the descent", rawLocal)]);

        // The inputs rescaled to [0, 1], no domain: greedy and its descent; the raw structure's image under the affine map, which
        // spans its columns exactly; and the image's descent.
        var box = new Scorer(BassModel.Create(Config(), Fixture(declared: false)));
        (List<BassTerm> greedy, double forward) = box.Forward([]);
        (List<BassTerm> local, double localScore, int localPasses) = box.Descend(greedy);
        output.WriteLine($"inputs rescaled to [0, 1], no domain: greedy to M {greedy.Count} at log π {forward:F2}; descent ({localPasses} passes) to M {local.Count} at {localScore:F2}, " +
                         $"E[σ² | S] {box.Sigma2(local):F4}: {Show(local)}");
        List<BassTerm> image = Closure(rawLocal);
        double imageScore = box.Score(image);
        output.WriteLine($"  the raw structure's image: M {image.Count}, log π {imageScore:F2}" + (double.IsFinite(imageScore) ? $", E[σ² | S] {box.Sigma2(image):F4}" : ", outside F"));
        (List<BassTerm> rich, double richScore, int richPasses) = box.Descend(image);
        output.WriteLine($"  the image's descent ({richPasses} passes): M {rich.Count}, log π {richScore:F2}, E[σ² | S] {box.Sigma2(rich):F4}: {Show(rich)}");
        output.WriteLine($"  {box.Calls + raw.Calls:N0} scores, {clock.Elapsed.TotalSeconds:F0} s");
        List<Trace> boxTraces = Chains(box, [("M = 0", Array.Empty<BassTerm>()), ("the image's descent", rich)]);

        // Real, not a trap: the image, though it fits better, scores far below the sparse mode; its descent returns to that mode;
        // and chains from M = 0 and from the richest state agree in each box.
        Assert.True(imageScore < localScore - 50.0, $"the image scores {imageScore:F2} against the mode's {localScore:F2}");
        Assert.True(rich.ToHashSet().SetEquals(local), "the image's descent reaches the sparse mode");
        foreach (List<Trace> traces in new[] { rawTraces, boxTraces })
        {
            double[] m = [.. traces.GroupBy(t => t.Start).Select(g => g.Average(t => t.M.Average()))];
            Assert.True(Math.Abs(m[0] - m[1]) < 0.3, $"M {m[0]:F2} from M = 0 against {m[1]:F2} from the richest state");
        }
        Assert.InRange(rawTraces.Average(t => t.Sigma2.Average()), 0.20, 0.26);
        Assert.InRange(boxTraces.Average(t => t.Sigma2.Average()), 0.80, 1.00);
    }

    // Ishigami is symmetric about ½ on the raw inputs: sin z₁ odd, sin² z₂ and z₃⁴ even, so its raw expansion takes odd degrees in
    // x₀ and even in x₁ and x₂. A box centred on ½ keeps that parity, since ψ_α(½ + c(x − ½)) spans the ψ_β of β ≤ α and α's parity;
    // a shifted box breaks it. Two declared boxes of one width, 1.04: centred on ½, and shifted by 0.01.
    [Fact]
    public void Symmetry_CentredAndShiftedBoxes()
    {
        if (!BassQualification.Enabled) return;   // heavy
        foreach ((string name, InputDomain domain) in new[] { ("centred on ½, [−0.02, 1.02]", new InputDomain(-0.02, 1.02)), ("shifted by 0.01, [−0.03, 1.01]", new InputDomain(-0.03, 1.01)) })
        {
            var scorer = new Scorer(BassModel.Create(Config(), Fixture(domain, rescale: false)));
            (List<BassTerm> greedy, _) = scorer.Forward([]);
            (List<BassTerm> local, double score, int passes) = scorer.Descend(greedy);
            output.WriteLine($"{name}: greedy to M {greedy.Count}, descent ({passes} passes) to M {local.Count} at log π {score:F2}, E[σ² | S] {scorer.Sigma2(local):F4}: {Show(local)}");
            List<Trace> traces = Chains(scorer, [("M = 0", Array.Empty<BassTerm>())]);
            double m = traces.Average(t => t.M.Average()), sigma2 = traces.Average(t => t.Sigma2.Average());
            if (Math.Abs(domain.Lower + domain.Upper - 1.0) < 1e-12)
            {
                Assert.True(m > 9.5 && sigma2 < 0.3, $"the centred box reads M {m:F2}, σ² {sigma2:F4}");
            }
            else
            {
                Assert.True(sigma2 > 0.5, $"the shifted box reads M {m:F2}, σ² {sigma2:F4}");
            }
        }
    }
}
