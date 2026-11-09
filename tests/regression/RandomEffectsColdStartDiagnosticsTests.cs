using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>Enumerated exits from the E21 slice 1 cold-start trap, with frozen nuisance coordinates and conditional quadrature.
/// The production kernel is exercised separately from the enumerated row. Runs behind <see cref="BassQualification"/>.</summary>
public sealed class RandomEffectsColdStartDiagnosticsTests(ITestOutputHelper output)
{
    private static readonly BassTerm[] Trap =
    [
        new([new HingeFactor(0, -1, 2), new HingeFactor(1, 1, 0)]),
        new([new HingeFactor(0, -1, 2), new HingeFactor(1, -1, 2)]),
    ];
    private static readonly string TrapKey = BassPriorSamplerTests.Key(Trap);

    private static BassModel Model(double scale, double logTau, bool factors = false, bool grouping = true, bool holdTau = true)
    {
        (BassModel held, BassConfig config) = RandomInterceptTests.Fixture("hinge");
        config = config with
        {
            Coefficients = holdTau ? config.Coefficients with { Scale = TauPrior.Held(Math.Exp(logTau)) } : config.Coefficients,
            RandomIntercept = grouping ? new BassRandomIntercept(scale, ScalePrior: new BassGroupScalePrior()) : null,
            Sampling = config.Sampling with { FactorMoves = factors ? new FactorMoves(1, 1, 0) : null },
        };
        BassData data = grouping ? held.Data : new BassData(held.Data.Predictors, held.Data.Response);
        return BassModel.Create(config, data);
    }

    private static string Key(IEnumerable<BassTerm> terms) => BassPriorSamplerTests.Key(terms);

    private sealed record Edge(string Move, BassTerm[] Terms, string Key, double Q, double ReverseQ, BassRefusal? Refusal,
                               double LogTargetRatio, double LogAcceptance, double Accepted);

    // The existing enumeration builds candidates from the move definitions. Aggregate all paths to an atom and obtain
    // its reverse from the same enumerator. Check the production within-move ratio against that independent reversal.
    private static List<Edge> Row(BassKernelEnumeration.Space space, int from)
    {
        BassKernelEnumeration.State source = space.States[from];
        var moves = new BassStructureMoves(space.Model, source.Terms, source.Loaded.LogTau).Moves;
        var palette = moves.ToDictionary(m => m.Key);
        double total = moves.Sum(m => m.Weight(source.Loaded));
        var row = new List<Edge>();
        foreach (IRjMove<BassState> move in moves.Where(m => m.Weight(source.Loaded) > 0))
        {
            double select = move.Weight(source.Loaded) / total;
            var atoms = BassKernelEnumeration.All(space, move.Key, source, null).Where(e => e.Candidate is not null)
                .GroupBy(e => Key(e.Candidate!));
            foreach (var atom in atoms)
            {
                BassTerm[] terms = atom.First().Candidate!;
                double q = select * atom.Sum(e => Math.Exp(e.LogProbability));
                int index = space.Find(terms);
                BassRefusal? refusal = null;
                if (index < 0)
                {
                    BassDesign strict = space.Model.NewDesign();
                    refusal = strict.Load(terms);
                    if (refusal is null) refusal = space.Model.Target.Refusal(BassStructureFit.Of(strict, terms.Sum(space.Model.Prior.LogPrior)));
                    Assert.NotNull(refusal);
                }
                if (terms.Distinct().Count() != terms.Length)
                {
                    row.Add(new Edge(move.Key, terms, atom.Key, q, 0, refusal, double.NegativeInfinity, double.NegativeInfinity, 0));
                    continue;
                }
                var candidate = new BassState(terms, BassStructureFit.Refused(BassRefusal.Numerical), source.Loaded.LogTau);
                BassKernelEnumeration.State destination = index >= 0 ? space.States[index]
                    : new BassKernelEnumeration.State(terms, [], candidate, BassMixing.Counts(space.Model, candidate),
                                                      double.NegativeInfinity, double.NegativeInfinity);
                IRjMove<BassState> reverseMove = palette[move.ReverseKey];
                double reverseSelect = reverseMove.Weight(destination.Loaded) / moves.Sum(m => m.Weight(destination.Loaded));
                double reverse = reverseSelect * BassKernelEnumeration.All(space, reverseMove.Key, destination, null)
                    .Where(e => e.Candidate is not null && Key(e.Candidate) == Key(source.Terms)).Sum(e => Math.Exp(e.LogProbability));
                Assert.True(reverse > 0, $"No reverse for {move.Key} to {atom.Key}");
                double logRatio = Math.Log(reverse / q);
                foreach (var path in atom)
                    Assert.InRange(Math.Abs(path.LogRatio + Math.Log(reverseSelect / select) - logRatio), 0, 1e-11);
                double targetRatio = destination.LogTarget - source.LogTarget;
                double logAcceptance = targetRatio + logRatio;
                row.Add(new Edge(move.Key, terms, atom.Key, q, reverse, refusal, targetRatio, logAcceptance, q * Math.Exp(Math.Min(0, logAcceptance))));
            }
        }
        Assert.InRange(row.Sum(e => e.Q), 0, 1.0 + 1e-12);
        return row;
    }

    [Theory]
    [InlineData(1.3, -4.0, 1.0, false, true)]
    [InlineData(0.05894, -3.94330, 1.0, false, true)]
    [InlineData(1.11035, -2.62045, 1.0, false, true)]
    [InlineData(1.3, -6.0, 1.0, false, true)]
    [InlineData(1.3, -4.0, 0.4, false, true)]
    [InlineData(1.3, -4.0, 1.0, true, true)]
    [InlineData(1.3, -4.0, 1.0, false, false)]
    public void FrozenExitRow_MatchesTheProductionProposals(double scale, double logTau, double b, bool factors, bool grouping)
    {
        if (!BassQualification.Enabled) return;
        BassModel model = Model(scale, logTau, factors, grouping);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model, b);
        int from = space.Find(Trap);
        Assert.True(from >= 0);
        List<Edge> row = Row(space, from);
        double exit = row.Sum(e => e.Accepted);
        double[] posterior = space.Posterior();
        output.WriteLine($"s={scale:G7}, ln τ={logTau:G7}, b={b}, factor moves={factors}, grouping={grouping}, |F|={space.States.Count}, " +
                         $"P(trap | s, τ, y)={posterior[from]:G8}, P(exit per structure step)={exit:G8}, reciprocal={1 / exit:G8}");
        foreach (var group in row.GroupBy(e => e.Move))
            output.WriteLine($"{group.Key}: proposed mass {group.Sum(e => e.Q):G7}, refused {group.Where(e => e.Refusal is not null).Sum(e => e.Q):G7}, " +
                             $"accepted exit {group.Sum(e => e.Accepted):G8}");
        foreach (Edge edge in row.Where(e => e.Refusal is null).OrderByDescending(e => e.Accepted).Take(8))
            output.WriteLine($"{edge.Move} -> {edge.Key}: q={edge.Q:G7}, reverse={edge.ReverseQ:G7}, Δlog target={edge.LogTargetRatio:F5}, " +
                             $"log A={edge.LogAcceptance:F5}, qα={edge.Accepted:G8}");

        // The reverse birth of either removed interaction on this two-predictor fixture has one forced variable set,
        // six sign/knot choices per factor, and the pseudo-count order probability (w+1)/(2w+1).
        foreach (Edge death in row.Where(e => e.Move == "death"))
        {
            double w = model.Config.Structure.OrderWeight;
            Assert.Equal(0.25, death.Q, 12);
            Assert.Equal((w + 1) / (2 * w + 1) / 36 / 3, death.ReverseQ, 12);
        }
        CheckProduction(space, from, row, b);
        if (scale == 1.3 && logTau == -4) DescribePaths(space, from);
    }

    private void DescribePaths(BassKernelEnumeration.Space space, int from)
    {
        var flows = BassKernelEnumeration.Kernel(space, (move, state) => BassKernelEnumeration.All(space, move, state, null));
        (_, double gap) = BassKernelEnumeration.Check(space, flows);
        var edges = new Dictionary<(int From, int To), double>();
        foreach (var group in flows.Values)
            foreach (var edge in group)
                edges[edge.Key] = edges.TryGetValue(edge.Key, out double old)
                    ? LogAdd(old, edge.Value) : edge.Value;
        var adjacent = edges.GroupBy(e => e.Key.From).ToDictionary(g => g.Key, g => g.ToArray());
        var reached = new HashSet<int> { from };
        var pending = new Queue<int>();
        pending.Enqueue(from);
        while (pending.TryDequeue(out int current))
            foreach (var edge in adjacent.GetValueOrDefault(current) ?? [])
                if (reached.Add(edge.Key.To)) pending.Enqueue(edge.Key.To);
        Assert.Equal(space.States.Count, reached.Count);
        int target = space.Find([new BassTerm([new HingeFactor(0, 1, 0)])]);
        Assert.True(target >= 0);
        double[] cost = [.. space.States.Select(_ => double.PositiveInfinity)];
        int[] previous = [.. space.States.Select(_ => -1)];
        var settled = new HashSet<int>();
        cost[from] = 0;
        for (int i = 0; i < cost.Length; i++)
        {
            int at = Enumerable.Range(0, cost.Length).Where(k => !settled.Contains(k)).MinBy(k => cost[k]);
            if (at == target) break;
            settled.Add(at);
            foreach (var edge in adjacent.GetValueOrDefault(at) ?? [])
            {
                double candidate = cost[at] - edge.Value;
                if (candidate >= cost[edge.Key.To]) continue;
                cost[edge.Key.To] = candidate;
                previous[edge.Key.To] = at;
            }
        }
        var path = new List<int> { target };
        while (path[^1] != from) path.Add(previous[path[^1]]);
        path.Reverse();
        output.WriteLine($"Graph: {reached.Count}/{space.States.Count} states reachable; detailed-balance log gap {gap:E3}; " +
                         $"largest-product path to x0+k0 has {path.Count - 1} steps and log product {-cost[target]:F5}");
        for (int i = 1; i < path.Count; i++)
            output.WriteLine($"path K={Math.Exp(edges[(path[i - 1], path[i])]):G8}: {Key(space.States[path[i - 1]].Terms)} -> {Key(space.States[path[i]].Terms)}");
        var outsidePairSector = Enumerable.Range(0, space.States.Count)
            .Where(i => space.States[i].Terms.Length != 2 || space.States[i].Terms.Any(t => t.Continuous.Count != 2)).ToHashSet();
        var mainEffects = new HashSet<int>
        {
            target, space.Find([new BassTerm([new HingeFactor(0, -1, 2)])]),
        };
        var sector = MeanPassage(space, edges, outsidePairSector);
        var mode = MeanPassage(space, edges, mainEffects);
        output.WriteLine($"Mean first passage from trap at frozen s, τ: leave the (2,2)-order sector " +
                         $"{sector.Times[from]:G9} structure steps; reach a one-term linear x0 representation {mode.Times[from]:G9}; " +
                         $"largest passage-equation residual {Math.Max(sector.Residual, mode.Residual):E3}");
    }

    // Reversibility symmetrizes the killed chain's rate matrix. The diagonal is accumulated from off-diagonal
    // transitions directly, avoiding subtraction of an almost-unit self-loop.
    private static (double[] Times, double Residual) MeanPassage(BassKernelEnumeration.Space space, Dictionary<(int From, int To), double> edges,
                                                                 HashSet<int> targets)
    {
        int[] transient = [.. Enumerable.Range(0, space.States.Count).Where(i => !targets.Contains(i))];
        var index = transient.Select((state, i) => (state, i)).ToDictionary(x => x.state, x => x.i);
        var matrix = new double[transient.Length, transient.Length];
        double top = transient.Max(i => space.States[i].LogPi);
        double[] roots = [.. transient.Select(i => Math.Exp((space.States[i].LogPi - top) / 2))];
        foreach (var edge in edges)
        {
            if (!index.TryGetValue(edge.Key.From, out int i)) continue;
            double p = Math.Exp(edge.Value);
            matrix[i, i] += p;
            if (index.TryGetValue(edge.Key.To, out int j)) matrix[i, j] -= p * roots[i] / roots[j];
        }
        var factor = new DenseCholesky(transient.Length);
        factor.Decompose(matrix);
        Assert.Equal(-1, factor.FailedAt);
        double[] solution = factor.Solve(roots);
        var result = new double[space.States.Count];
        for (int i = 0; i < transient.Length; i++) result[transient[i]] = solution[i] / roots[i];
        double worst = 0;
        foreach (int state in transient)
        {
            double residual = edges.Where(e => e.Key.From == state)
                .Sum(e => Math.Exp(e.Value) * (result[state] - result[e.Key.To])) - 1;
            worst = Math.Max(worst, Math.Abs(residual));
        }
        Assert.InRange(worst, 0, 1e-7);
        return (result, worst);
    }

    private static double LogAdd(double a, double b)
    {
        double top = Math.Max(a, b);
        return top + Math.Log(Math.Exp(a - top) + Math.Exp(b - top));
    }

    private sealed record Integrated(double LogMass, double LogScale, double LogTau, double Exit, double Death, double Prune);

    [Fact]
    public void ConditionalNuisanceLaw_PricesTheTrapsExits()
    {
        if (!BassQualification.Enabled) return;
        foreach (bool factors in new[] { false, true })
        {
            Integrated[] coarse = Integrate(factors, fixedScale: null, step: 0.2, lowerScale: -18, upperScale: 10, lowerTau: -35, upperTau: 4);
            Integrated[] fine = Integrate(factors, fixedScale: null, step: 0.1, lowerScale: -22, upperScale: 12, lowerTau: -40, upperTau: 5);
            Integrated[] held = Integrate(factors, fixedScale: 1.3, step: 0.1, lowerScale: 0, upperScale: 0, lowerTau: -40, upperTau: 5);
            for (int i = 0; i < fine.Length; i++)
            {
                Integrated a = coarse[i], b = fine[i];
                double change = new[] { a.LogMass - b.LogMass, a.LogScale - b.LogScale, a.LogTau - b.LogTau, a.Exit - b.Exit, a.Death - b.Death, a.Prune - b.Prune }
                    .Max(Math.Abs);
                output.WriteLine($"Nuisance conditional on trap, factors={factors}, b={(i == 0 ? 1 : 0.4)}: max quadrature change {change:E3}; " +
                                 $"E[ln s]={b.LogScale:F7}, E[ln τ]={b.LogTau:F7}, exit={b.Exit:G9}, death={b.Death:G9}, prune={b.Prune:G9}; " +
                                 $"at s=1.3 with τ integrated: exit={held[i].Exit:G9}, death={held[i].Death:G9}, prune={held[i].Prune:G9}");
                Assert.InRange(change, 0, 1e-5);
            }
        }
    }

    private static Integrated[] Integrate(bool factors, double? fixedScale, double step, double lowerScale, double upperScale,
                                          double lowerTau, double upperTau)
    {
        BassModel frozen = Model(1.3, -4, factors);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(frozen);
        List<Edge> edges = Row(space, space.Find(Trap));
        BassModel model = Model(1.3, -4, factors, holdTau: false);
        double[] scales = fixedScale is double held ? [Math.Log(held)]
            : [.. Enumerable.Range(0, (int)Math.Round((upperScale - lowerScale) / step) + 1).Select(i => lowerScale + i * step)];
        double[] taus = [.. Enumerable.Range(0, (int)Math.Round((upperTau - lowerTau) / step) + 1).Select(i => lowerTau + i * step)];
        double[] rungs = [1.0, 0.4];
        var mass = new double[rungs.Length];
        var scaleMean = new double[rungs.Length];
        var tauMean = new double[rungs.Length];
        var escape = new double[rungs.Length];
        var deaths = new double[rungs.Length];
        var prunes = new double[rungs.Length];
        BassDesign pivot = model.NewDesign();
        Assert.Null(pivot.Load(Trap));
        BassStructureFit pivotFit = BassStructureFit.Of(pivot, Trap.Sum(model.Prior.LogPrior));
        BassEvaluation pivotValue = model.Target.Evaluate(pivotFit, -4);
        double[] shifts = [.. rungs.Select(b => pivotValue.LogPrior + pivotValue.LogLikelihood(b))];
        for (int v = 0; v < scales.Length; v++)
        {
            WorkingResponse working = model.RandomIntercept!.At(Math.Exp(scales[v])).Whitened(model.BaseReading);
            BassDesign design = model.NewDesign(working);
            Assert.Null(design.Load(Trap));
            BassStructureFit source = BassStructureFit.Of(design, Trap.Sum(model.Prior.LogPrior));
            var fits = new BassStructureFit?[edges.Count];
            for (int k = 0; k < edges.Count; k++)
                if (design.Load(edges[k].Terms) is null)
                    fits[k] = BassStructureFit.Of(design, edges[k].Terms.Sum(model.Prior.LogPrior));
            for (int j = 0; j < taus.Length; j++)
            {
                BassEvaluation current = model.Target.Evaluate(source, taus[j]);
                double quadrature = (fixedScale is not null || v > 0 && v < scales.Length - 1 ? 1 : 0.5)
                                    * (j > 0 && j < taus.Length - 1 ? 1 : 0.5);
                double[] density = [.. rungs.Select(b => current.LogPrior + current.LogLikelihood(b))];
                double[] weight = [.. density.Select((l, r) => Math.Exp(l - shifts[r]) * quadrature)];
                for (int r = 0; r < rungs.Length; r++)
                {
                    mass[r] += weight[r];
                    scaleMean[r] += weight[r] * scales[v];
                    tauMean[r] += weight[r] * taus[j];
                }
                for (int k = 0; k < edges.Count; k++)
                {
                    if (fits[k] is not { } fit) continue;
                    Edge edge = edges[k];
                    BassEvaluation next = model.Target.Evaluate(fit, taus[j]);
                    double ratio = Math.Log(edge.ReverseQ / edge.Q);
                    for (int r = 0; r < rungs.Length; r++)
                    {
                        double acceptance = Math.Exp(Math.Min(0, next.LogPrior + next.LogLikelihood(rungs[r]) - density[r] + ratio));
                        double contribution = weight[r] * edge.Q * acceptance;
                        escape[r] += contribution;
                        if (edge.Move == "death") deaths[r] += contribution;
                        if (edge.Move == "prune") prunes[r] += contribution;
                    }
                }
            }
        }
        double logCell = Math.Log(step) * (fixedScale is null ? 2 : 1);
        return [.. Enumerable.Range(0, rungs.Length).Select(r =>
            new Integrated(shifts[r] + Math.Log(mass[r]) + logCell, scaleMean[r] / mass[r], tauMean[r] / mass[r],
                           escape[r] / mass[r], deaths[r] / mass[r], prunes[r] / mass[r]))];
    }

    private void CheckProduction(BassKernelEnumeration.Space space, int from, List<Edge> row, double b)
    {
        BassModel model = space.Model;
        BassKernelEnumeration.State source = space.States[from];
        var structure = new BassStructureMoves(model, source.Terms, source.Loaded.LogTau);
        var temperature = new Temperature(b, TemperingPath.Likelihood);
        var chain = new ReversibleJumpChain<BassState>(structure.Moves, model.Target, structure.Start,
                                                       new Xoshiro256PlusPlus(SeedTree.Derive(20261129, 2)[0]), temperature);
        BassState start = chain.Current;
        var rng = new Xoshiro256PlusPlus(SeedTree.Derive(20261129, 2)[1]);
        var expected = row.ToDictionary(e => (e.Move, e.Key));
        var counts = new Dictionary<(string Move, string Key), int>();
        var moves = structure.Moves;
        double total = moves.Sum(m => m.Weight(start));
        const int proposals = 40_000, trials = 150_000;
        double worstRatio = 0, worstTarget = 0;
        for (int i = 0; i < proposals; i++)
        {
            double pick = rng.NextDouble() * total, cumulative = 0;
            IRjMove<BassState> move = moves.Last();
            foreach (var option in moves)
            {
                cumulative += option.Weight(start);
                if (pick < cumulative) { move = option; break; }
            }
            if (move.Propose(start, rng) is not { } proposed || Key(proposed.Candidate.Terms) == TrapKey) continue;
            var key = (move.Key, Key(proposed.Candidate.Terms));
            Assert.True(expected.TryGetValue(key, out Edge? edge), $"Unenumerated proposal {key}");
            counts[key] = counts.GetValueOrDefault(key) + 1;
            Assert.Equal(edge!.Refusal, model.Target.Refusal(proposed.Candidate.Fit));
            if (edge.Refusal is not null) continue;
            var reverse = moves.Single(m => m.Key == move.ReverseKey);
            double selectRatio = Math.Log(reverse.Weight(proposed.Candidate) / moves.Sum(m => m.Weight(proposed.Candidate)))
                                 - Math.Log(move.Weight(start) / total);
            double ratioGap = Math.Abs(proposed.LogProposalRatio + selectRatio - Math.Log(edge.ReverseQ / edge.Q));
            double targetGap = Math.Abs(temperature.LogDensity(model.Target.Evaluate(proposed.Candidate)) - chain.CurrentLogDensity - edge.LogTargetRatio);
            worstRatio = Math.Max(worstRatio, ratioGap);
            worstTarget = Math.Max(worstTarget, targetGap);
        }
        foreach (Edge edge in row)
        {
            int n = counts.GetValueOrDefault((edge.Move, edge.Key));
            double error = Math.Sqrt(proposals * edge.Q * (1 - edge.Q));
            Assert.True(Math.Abs(n - proposals * edge.Q) <= 6 * error + 2, $"Proposal frequency {edge.Move}, {edge.Key}: {n}/{proposals}, expected {edge.Q}");
        }
        int exits = 0;
        var destinations = new Dictionary<string, int>();
        for (int i = 0; i < trials; i++)
        {
            chain.SetState(start, model.Target.Evaluate(start));
            chain.Step();
            if (Key(chain.Current.Terms) == TrapKey) continue;
            exits++;
            string key = Key(chain.Current.Terms);
            destinations[key] = destinations.GetValueOrDefault(key) + 1;
        }
        double probability = row.Sum(e => e.Accepted);
        double se = Math.Sqrt(trials * probability * (1 - probability));
        output.WriteLine($"Production frozen trials: {exits}/{trials} exits, expected {trials * probability:G7}; " +
                         $"largest log proposal-ratio gap {worstRatio:E3}, target gap {worstTarget:E3}");
        foreach (var destination in destinations.OrderByDescending(p => p.Value).Take(4))
            output.WriteLine($"accepted {destination.Key}: {destination.Value}");
        Assert.InRange(worstRatio, 0, 1e-10);
        Assert.InRange(worstTarget, 0, 1e-7);
        Assert.True(Math.Abs(exits - trials * probability) <= 6 * se + 2, $"Observed {exits}, expected {trials * probability}");
    }
}
