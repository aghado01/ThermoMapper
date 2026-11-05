using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.Bgp;
using Maths.Regression.Spline;
using Maths.Regression.Spline.Baps;
using Maths.Regression.Spline.Bars;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Evidence;
using Maths.Samplers.FlatHistogram;
using Maths.Samplers.Rjmcmc;
using Maths.Samplers.Smc;

namespace Maths.Regression.Tests;

// ---- The cheap toy: a birth-death chain on a count k, target ∝ Poisson(4) · exp(−0.05 (k − 4)²) ----

internal sealed class CountEvaluation(double logPrior, double ell) : ITargetEvaluation
{
    public double LogPrior { get; } = logPrior;
    public double LogLikelihood(double inverseTemperature) => inverseTemperature * ell;
}

internal sealed class CountTarget : ISamplingTarget<int>
{
    internal static readonly double[] LogFactorial = BuildLogFactorial();
    private static double[] BuildLogFactorial()
    {
        var t = new double[4096];
        for (int k = 1; k < t.Length; k++) t[k] = t[k - 1] + Math.Log(k);
        return t;
    }

    internal static double LogPrior(int k) => k >= 0 && k < LogFactorial.Length ? k * Math.Log(4.0) - LogFactorial[k] : double.NegativeInfinity;
    internal static double Ell(int k) => -0.05 * (k - 4) * (k - 4);
    public bool HasCompleteValues => false;
    public ITargetEvaluation Evaluate(int k) => new CountEvaluation(LogPrior(k), Ell(k));
}

internal sealed class CountBirth(double weight = 0.4) : IRjMove<int>
{
    public string Key => "birth";
    public string ReverseKey => "death";
    public double Weight(int state) => weight;
    public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng) => new Proposal<int>(current + 1, 0.0, 0.0);
}

internal sealed class CountDeath(double weight = 0.4) : IRjMove<int>
{
    public string Key => "death";
    public string ReverseKey => "birth";
    public double Weight(int state) => weight;
    public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng) => current == 0 ? null : new Proposal<int>(current - 1, 0.0, 0.0);
}

internal sealed class CountWalk(double weight = 0.2) : IRjMove<int>
{
    public string Key => "walk";
    public string ReverseKey => "walk";
    public double Weight(int state) => weight;
    public Proposal<int>? Propose(int current, Xoshiro256PlusPlus rng)
    {
        if (current == 0) return new Proposal<int>(1, -Math.Log(2.0), 0.0);
        bool up = rng.NextDouble() < 0.5;
        return new Proposal<int>(up ? current + 1 : current - 1, 0.0, 0.0);
    }
}

// ---- A struct-constrained, static-abstract replica of the same step ----

internal interface IPalette { static abstract int Count { get; } static abstract double Weight(int move, int k); static abstract bool Propose(int move, int k, Xoshiro256PlusPlus rng, out int candidate, out double logRatio); static abstract int Reverse(int move); }
internal interface IDensityOf { static abstract double Prior(int k); static abstract double Ell(int k); }

internal readonly struct CountPalette : IPalette
{
    public static int Count => 3;
    public static double Weight(int move, int k) => move switch { 0 => 0.4, 1 => 0.4, _ => 0.2 };
    public static int Reverse(int move) => move switch { 0 => 1, 1 => 0, _ => 2 };
    public static bool Propose(int move, int k, Xoshiro256PlusPlus rng, out int candidate, out double logRatio)
    {
        logRatio = 0.0;
        switch (move)
        {
            case 0: candidate = k + 1; return true;
            case 1: candidate = k - 1; return k != 0;
            default:
                if (k == 0) { candidate = 1; logRatio = -Math.Log(2.0); return true; }
                candidate = rng.NextDouble() < 0.5 ? k + 1 : k - 1;
                return true;
        }
    }
}

internal readonly struct CountDensity : IDensityOf
{
    public static double Prior(int k) => CountTarget.LogPrior(k);
    public static double Ell(int k) => CountTarget.Ell(k);
}

internal sealed class GenericCountChain<TP, TD> where TP : struct, IPalette where TD : struct, IDensityOf
{
    private readonly Xoshiro256PlusPlus _rng;
    private double _currentLogDensity;
    public int Current { get; private set; }
    public long Attempts { get; private set; }
    public long Accepted { get; private set; }

    public GenericCountChain(int start, Xoshiro256PlusPlus rng)
    {
        Current = start;
        _rng = rng;
        _currentLogDensity = TD.Prior(start) + 1.0 * TD.Ell(start);
    }

    public int Step()
    {
        double totalCur = 0.0;
        for (int i = 0; i < TP.Count; i++) { double w = TP.Weight(i, Current); if (w > 0.0) totalCur += w; }
        if (!(totalCur > 0.0)) return Current;
        int idx = TP.Count - 1;
        {
            double u = _rng.NextDouble() * totalCur, cumulative = 0.0;
            for (int i = 0; i < TP.Count; i++)
            {
                double w = TP.Weight(i, Current);
                if (w <= 0.0) continue;
                cumulative += w;
                if (u < cumulative) { idx = i; break; }
            }
        }
        if (!TP.Propose(idx, Current, _rng, out int candidate, out double logProposal)) return Current;
        Attempts++;
        double totalCand = 0.0;
        for (int i = 0; i < TP.Count; i++) { double w = TP.Weight(i, candidate); if (w > 0.0) totalCand += w; }
        double wReverse = TP.Weight(TP.Reverse(idx), candidate);
        if (!(wReverse > 0.0) || !(totalCand > 0.0)) return Current;
        double logSelection = Math.Log(wReverse / totalCand) - Math.Log(TP.Weight(idx, Current) / totalCur);
        double prior = TD.Prior(candidate);
        double candidateLogDensity = double.IsNegativeInfinity(prior) ? double.NegativeInfinity : prior + 1.0 * TD.Ell(candidate);
        double logAccept = candidateLogDensity > double.NegativeInfinity
            ? (candidateLogDensity - _currentLogDensity) + logProposal + 0.0 + logSelection
            : double.NegativeInfinity;
        if (logAccept >= 0.0 || _rng.NextDouble() < Math.Exp(logAccept))
        {
            Current = candidate;
            _currentLogDensity = candidateLogDensity;
            Accepted++;
        }
        return Current;
    }
}

/// <summary>
/// <see cref="ReversibleJumpChain{TState}.Step"/> copied verbatim (the interface everywhere), with two switches that remove
/// what the read of the code flags: the string-keyed reverse lookup, and the weights computed afresh at the current state.
/// </summary>
internal sealed class StepReplica(IReadOnlyList<IRjMove<int>> moves, ISamplingTarget<int> target, int start, Xoshiro256PlusPlus rng, bool dictionary, bool cacheWeights)
{
    private readonly Dictionary<string, IRjMove<int>> _byKey = moves.ToDictionary(m => m.Key);
    private readonly int[] _reverse = moves.Select(m => moves.ToList().FindIndex(r => r.Key == m.ReverseKey)).ToArray();
    private readonly double[] _weights = new double[moves.Count];
    private ITargetEvaluation _currentEvaluation = target.Evaluate(start);
    private double _currentLogDensity = Temperature.Posterior.LogDensity(target.Evaluate(start));
    public int Current { get; private set; } = start;
    public long Accepted { get; private set; }

    public int Step()
    {
        double totalCur = TotalWeight(Current, fill: cacheWeights);
        if (!(totalCur > 0.0)) return Current;
        int idx = SelectMove(Current, totalCur);
        IRjMove<int> move = moves[idx];
        if (move.Propose(Current, rng) is not Proposal<int> p) return Current;
        double totalCand = TotalWeight(p.Candidate, fill: false);
        double wReverse = (dictionary ? _byKey[move.ReverseKey] : moves[_reverse[idx]]).Weight(p.Candidate);
        if (!(wReverse > 0.0) || !(totalCand > 0.0)) return Current;
        double wCurrent = cacheWeights ? _weights[idx] : move.Weight(Current);
        double logSelection = Math.Log(wReverse / totalCand) - Math.Log(wCurrent / totalCur);
        ITargetEvaluation candidate = target.Evaluate(p.Candidate);
        double candidateLogDensity = Temperature.Posterior.LogDensity(candidate);
        double logAccept = candidateLogDensity > double.NegativeInfinity
            ? (candidateLogDensity - _currentLogDensity) + p.LogProposalRatio + p.LogJacobian + logSelection
            : double.NegativeInfinity;
        if (logAccept >= 0.0 || rng.NextDouble() < Math.Exp(logAccept))
        {
            Current = p.Candidate;
            _currentEvaluation = candidate;
            _currentLogDensity = candidateLogDensity;
            Accepted++;
        }
        return Current;
    }

    private double TotalWeight(int state, bool fill)
    {
        double total = 0.0;
        for (int i = 0; i < moves.Count; i++)
        {
            double w = moves[i].Weight(state);
            if (fill) _weights[i] = w;
            if (w > 0.0) total += w;
        }
        return total;
    }

    private int SelectMove(int state, double total)
    {
        if (moves.Count == 1) return 0;
        double u = rng.NextDouble() * total, cumulative = 0.0;
        for (int i = 0; i < moves.Count; i++)
        {
            double w = cacheWeights ? _weights[i] : moves[i].Weight(state);
            if (w <= 0.0) continue;
            cumulative += w;
            if (u < cumulative) return i;
        }
        return moves.Count - 1;
    }
}

// ---- Densities for the block samplers ----

internal interface IScalarDensity { double Eval(double x); }
internal readonly struct GaussianDensity : IScalarDensity { public double Eval(double x) => -0.5 * x * x; }
internal readonly struct RemlDensity(PenalizedSpline spline) : IScalarDensity { public double Eval(double u) => spline.RemlLogEvidence(Math.Exp(u)); }
internal readonly struct GpDensity(GpRegression gp) : IScalarDensity { public double Eval(double u) => gp.Fit(Math.Exp(u)).LogMarginal; }

internal static class GenericBlocks
{
    // SliceSampler.Step, copied, over a struct density: the same draws in the same order.
    internal static int Slice<TD>(ref double x, ref double logTarget, double width, int maxSteps, TD density, Xoshiro256PlusPlus rng) where TD : struct, IScalarDensity
    {
        double level = logTarget + Math.Log(1.0 - rng.NextDouble());
        double left = x - width * rng.NextDouble(), right = left + width;
        int toLeft = (int)(maxSteps * rng.NextDouble()), toRight = maxSteps - 1 - toLeft;
        int evaluations = 0;
        while (toLeft > 0) { evaluations++; if (!(density.Eval(left) > level)) break; left -= width; toLeft--; }
        while (toRight > 0) { evaluations++; if (!(density.Eval(right) > level)) break; right += width; toRight--; }
        while (true)
        {
            double candidate = left + rng.NextDouble() * (right - left);
            double value = density.Eval(candidate);
            evaluations++;
            if (value > level) { x = candidate; logTarget = value; return evaluations; }
            if (candidate < x) left = candidate; else right = candidate;
            if (!(right - left > 1e-14 * (1.0 + Math.Abs(x)))) return evaluations;
        }
    }

    internal static bool RandomWalk<TD>(ref double x, ref double logTarget, double step, TD density, Xoshiro256PlusPlus rng) where TD : struct, IScalarDensity
    {
        double proposal = x + step * Normal.SampleStandard(rng);
        double proposalTarget = density.Eval(proposal);
        if (!(Math.Log(1.0 - rng.NextDouble()) < proposalTarget - logTarget)) return false;
        x = proposal;
        logTarget = proposalTarget;
        return true;
    }
}

// ---- A tempered Gaussian target for SMC, cheap or with an expensive likelihood ----

internal sealed class GaussEvaluation(double logPrior, double ell) : ITargetEvaluation
{
    public double LogPrior { get; } = logPrior;
    public double LogLikelihood(double b) => b * ell;
    public double LogLikelihoodSlope(double b) => ell;
}

internal sealed class GaussTarget(int terms) : ISamplingTarget<double[]>
{
    private readonly double[] _centres = Enumerable.Range(0, Math.Max(terms, 1)).Select(j => -3.0 + 6.0 * j / Math.Max(terms, 1)).ToArray();
    public bool HasCompleteValues => false;
    public ITargetEvaluation Evaluate(double[] x)
    {
        double prior = 0, ell = 0;
        foreach (double v in x) { prior -= 0.5 * v * v; double z = (v - 1.0) / 0.3; ell -= 0.5 * z * z; }
        if (terms > 0)
        {
            double s = 0;
            for (int j = 0; j < terms; j++) { double d = x[0] - _centres[j]; s += Math.Exp(-d * d); }
            ell += 0.001 * s;
        }
        return new GaussEvaluation(prior, ell);
    }
}

internal sealed class GaussParticles(int dim) : ISmcParticles<double[]>
{
    public double[] Draw(Xoshiro256PlusPlus rng) { var x = new double[dim]; for (int i = 0; i < dim; i++) x[i] = Normal.SampleStandard(rng); return x; }
    public double[] Fork(double[] state) => (double[])state.Clone();
}

internal sealed class IidChainRun : IChainRun<double>, IChain<double>
{
    private readonly Xoshiro256PlusPlus _rng;
    public IidChainRun(int seed) => _rng = new Xoshiro256PlusPlus(seed);
    public IChain<double> Chain => this;
    public long Accepted => 0;
    public long Attempts => 0;
    public double Step() => Normal.SampleStandard(_rng);
    public void Burn(int steps) { for (int i = 0; i < steps; i++) Step(); }
    public void Accumulate(in double draw, Span<double> functionals) => functionals[0] = draw;
}

internal sealed class IidModel : IEnsembleModel<double>
{
    public int FunctionalDim => 1;
    public int MonitoredDim => 1;
    public IChainRun<double> StartChain(int chainSeed, int readoutSeed) => new IidChainRun(chainSeed);
}

internal static class CostSamplers
{
    private static double Sink { get => CostKernels.Sink; set => CostKernels.Sink = value; }
    private static object? ObjectSink { set => CostKernels.ObjectSink = value; }

    private static IReadOnlyList<IRjMove<int>> CountMoves() => new IRjMove<int>[] { new CountBirth(), new CountDeath(), new CountWalk() };

    // ================================================================ sampler:rj

    internal static void RjBlock()
    {
        const string block = "sampler:rj";
        var target = new CountTarget();
        const int check = 20_000;

        // Trajectories first: every replica must walk the production chain's path.
        var production = new ReversibleJumpChain<int>(CountMoves(), target, 4, new Xoshiro256PlusPlus(11));
        var r0 = new StepReplica(CountMoves(), target, 4, new Xoshiro256PlusPlus(11), dictionary: true, cacheWeights: false);
        var r1 = new StepReplica(CountMoves(), target, 4, new Xoshiro256PlusPlus(11), dictionary: false, cacheWeights: true);
        var r2 = new GenericCountChain<CountPalette, CountDensity>(4, new Xoshiro256PlusPlus(11));
        long hash0 = 0, hash1 = 0, hash2 = 0, hashP = 0;
        for (int i = 0; i < check; i++)
        {
            hashP = hashP * 31 + production.Step();
            hash0 = hash0 * 31 + r0.Step();
            hash1 = hash1 * 31 + r1.Step();
            hash2 = hash2 * 31 + r2.Step();
        }
        if (!(hashP == hash0 && hashP == hash1 && hashP == hash2 && production.Accepted == r0.Accepted && production.Accepted == r1.Accepted && production.Accepted == r2.Accepted))
            throw new InvalidOperationException("A step replica left the production chain's path.");
        CostProfile.Row(block, "cheap", "(check)", "trajectories_identical_steps", check);

        var chain = new ReversibleJumpChain<int>(CountMoves(), target, 4, new Xoshiro256PlusPlus(12));
        var rep0 = new StepReplica(CountMoves(), target, 4, new Xoshiro256PlusPlus(12), true, false);
        var rep1 = new StepReplica(CountMoves(), target, 4, new Xoshiro256PlusPlus(12), false, true);
        var gen = new GenericCountChain<CountPalette, CountDensity>(4, new Xoshiro256PlusPlus(12));
        IReadOnlyList<IRjMove<int>> palette = CountMoves();
        var byKey = palette.ToDictionary(m => m.Key);
        var rng = new Xoshiro256PlusPlus(5);
        var parts = new List<CostProfile.Part>
        {
            new("step_production", c => { for (int i = 0; i < c; i++) chain.Step(); }, 1, "step"),
            new("step_replica_verbatim", c => { for (int i = 0; i < c; i++) rep0.Step(); }, 1, "step"),
            new("step_replica_no_dictionary_cached_weights", c => { for (int i = 0; i < c; i++) rep1.Step(); }, 1, "step"),
            new("step_replica_struct_generic_no_allocation", c => { for (int i = 0; i < c; i++) gen.Step(); }, 1, "step"),
            new("total_weight_3_moves", c => { for (int i = 0; i < c; i++) { double t = 0; for (int m = 0; m < palette.Count; m++) { double w = palette[m].Weight(i); if (w > 0) t += w; } Sink = t; } }, 3, "weight call"),
            new("reverse_key_dictionary_lookup", c => { for (int i = 0; i < c; i++) Sink = byKey[palette[i % 3].ReverseKey].Weight(i); }, 1, "lookup"),
            new("target_evaluate_allocating", c => { for (int i = 0; i < c; i++) ObjectSink = target.Evaluate(i & 15); }, 1, "evaluation"),
            new("temperature_log_density", c => { var ev = target.Evaluate(5); for (int i = 0; i < c; i++) Sink = Temperature.Posterior.LogDensity(ev); }, 1, "call"),
            new("move_propose", c => { for (int i = 0; i < c; i++) ObjectSink = palette[2].Propose(i & 15, rng); }, 1, "proposal"),
        };
        CostProfile.Point(block, "cheap target (count)", parts);
        CostProfile.Row(block, "cheap target (count)", "(count)", "acceptance_rate", chain.Accepted / (double)Math.Max(1, chain.Attempts));
        CostProfile.Row(block, "cheap target (count)", "(count)", "attempts_per_step_after_run", chain.Attempts > 0 ? 1.0 : 0.0);

        // The expensive target: BARS's evaluation, with the same chain code. The step's overhead is the cheap target's, whatever
        // the evaluation costs; the share is read from the two.
        foreach (int n in new[] { 100, 500 })
        {
            var fx = new CostEngines.BarsFixture(n);
            ReversibleJumpChain<KnotConfig> bars = fx.Chain(3);
            for (int i = 0; i < 1000; i++) bars.Step();
            CostProfile.Point(block, $"expensive target (BARS n={n})", new List<CostProfile.Part>
            {
                new("step_production", c => { for (int i = 0; i < c; i++) bars.Step(); }, 1, "step"),
                new("target_evaluate", c => { KnotConfig cfg = bars.Current; for (int i = 0; i < c; i++) ObjectSink = fx.Target.Evaluate(cfg); }, 1, "evaluation"),
            });
        }
    }

    // ================================================================ sampler:blocks

    internal static void BlocksBlock()
    {
        const string block = "sampler:blocks";
        var cases = new List<(string Name, Func<double, double> F, IScalarDensity D, double Start)>();
        cases.Add(("cheap Gaussian", x => -0.5 * x * x, new GaussianDensity(), 0.3));
        {
            var rng = new Xoshiro256PlusPlus(91);
            int n = 1000, nu = 104;
            double[] xs = Enumerable.Range(0, n).Select(i => (i + 0.5) / n).ToArray();
            double[] ys = xs.Select(t => Math.Sin(2 * Math.PI * t) + 0.15 * (rng.NextDouble() - 0.5)).ToArray();
            double[,] z = new SplineBasis(3).Design(CostEngines.BarsFixture.Equispaced(nu - 4), xs);
            var ps = new PenalizedSpline(z, ys, new DifferencePenalty(2));
            cases.Add(("expensive: BAPS REML evidence, nu=104", u => ps.RemlLogEvidence(Math.Exp(u)), new RemlDensity(ps), 1.0));
        }
        {
            (double[,] x, double[] y) = (new double[200, 3], new double[200]);
            var rng = new Xoshiro256PlusPlus(92);
            for (int i = 0; i < 200; i++) { for (int c = 0; c < 3; c++) x[i, c] = rng.NextDouble(); y[i] = Math.Sin(6 * x[i, 0]) + 0.1 * rng.NextDouble(); }
            var gp = new GpRegression(x, y, 0.01, new SquaredExponentialKernel());
            cases.Add(("expensive: GP evidence, n=200", u => gp.Fit(Math.Exp(u)).LogMarginal, new GpDensity(gp), -2.0));
        }
        foreach ((string name, Func<double, double> f, IScalarDensity d, double start) in cases)
        {
            double x1 = start, lt1 = f(start), x2 = start, lt2 = f(start), x3 = start, lt3 = f(start), x4 = start, lt4 = f(start);
            var rngA = new Xoshiro256PlusPlus(5); var rngB = new Xoshiro256PlusPlus(5); var rngC = new Xoshiro256PlusPlus(6); var rngD = new Xoshiro256PlusPlus(6);
            long evalsA = 0, evalsB = 0, steps = 0;
            // Path identity: the generic replicas draw what production draws.
            double xa = start, la = f(start), xb = start, lb = f(start);
            for (int i = 0; i < 200; i++)
            {
                evalsA += SliceSampler.Step(ref xa, ref la, 1.0, 10, f, rngA);
                evalsB += d switch
                {
                    GaussianDensity g => GenericBlocks.Slice(ref xb, ref lb, 1.0, 10, g, rngB),
                    RemlDensity r => GenericBlocks.Slice(ref xb, ref lb, 1.0, 10, r, rngB),
                    GpDensity p => GenericBlocks.Slice(ref xb, ref lb, 1.0, 10, p, rngB),
                    _ => 0,
                };
                steps++;
            }
            if (BitConverter.DoubleToInt64Bits(xa) != BitConverter.DoubleToInt64Bits(xb) || evalsA != evalsB)
                throw new InvalidOperationException($"The generic slice replica left the production path ({name}).");
            CostProfile.Row(block, name, "(count)", "evaluations_per_slice_step", evalsA / (double)steps);
            var parts = new List<CostProfile.Part>
            {
                new("slice_delegate_production", c => { for (int i = 0; i < c; i++) SliceSampler.Step(ref x1, ref lt1, 1.0, 10, f, rngC); }, 1, "step"),
                new("slice_struct_generic", c => { for (int i = 0; i < c; i++) _ = d switch { GaussianDensity g => GenericBlocks.Slice(ref x2, ref lt2, 1.0, 10, g, rngD), RemlDensity r => GenericBlocks.Slice(ref x2, ref lt2, 1.0, 10, r, rngD), GpDensity p => GenericBlocks.Slice(ref x2, ref lt2, 1.0, 10, p, rngD), _ => 0 }; }, 1, "step"),
                new("random_walk_delegate_production", c => { for (int i = 0; i < c; i++) RandomWalkMetropolis.Step(ref x3, ref lt3, 0.5, f, rngA); }, 1, "step"),
                new("random_walk_struct_generic", c => { for (int i = 0; i < c; i++) _ = d switch { GaussianDensity g => GenericBlocks.RandomWalk(ref x4, ref lt4, 0.5, g, rngB), RemlDensity r => GenericBlocks.RandomWalk(ref x4, ref lt4, 0.5, r, rngB), GpDensity p => GenericBlocks.RandomWalk(ref x4, ref lt4, 0.5, p, rngB), _ => false }; }, 1, "step"),
                new("density_evaluation_alone", c => { for (int i = 0; i < c; i++) Sink = f(start + 1e-3 * (i & 7)); }, 1, "evaluation"),
            };
            CostProfile.Point(block, name, parts);
        }
    }

    // ================================================================ sampler:pt

    internal static void TemperingBlock()
    {
        const string block = "sampler:pt";
        foreach (int levels in new[] { 5, 10, 20 })
        {
            var cheap = new ParallelTempering<int>(_ => CountMoves(), new CountTarget(), 4, ParallelTempering<int>.GeometricLadder(levels, 0.01),
                                                   TemperingPath.Likelihood, 5);
            for (int i = 0; i < 200; i++) cheap.Step();
            long swaps0 = cheap.SwapAttempts;
            var parts = new List<CostProfile.Part>
            {
                new("scan_whole", c => { for (int i = 0; i < c; i++) cheap.Step(); }, levels, "rung step"),
                new("one_replica_step", c => { for (int i = 0; i < c; i++) cheap.Replicas[i % levels].Step(); }, 1, "rung step"),
            };
            CostProfile.Point(block, $"cheap target, {levels} rungs", parts);
            CostProfile.Row(block, $"cheap target, {levels} rungs", "(count)", "swap_attempts_per_scan", (cheap.SwapAttempts - swaps0) / (double)Math.Max(1, cheap.Scans - 200));
        }
        foreach (int levels in new[] { 5, 10 })
        {
            var fx = new CostEngines.BarsFixture(100);
            var ladder = new ParallelTempering<KnotConfig>(_ => fx.Moves(), fx.Target, CostEngines.BarsFixture.Equispaced(4),
                                                           ParallelTempering<KnotConfig>.GeometricLadder(levels, 0.01), TemperingPath.Likelihood, 7);
            for (int i = 0; i < 300; i++) ladder.Step();
            CostProfile.Point(block, $"BARS n=100, {levels} rungs", new List<CostProfile.Part>
            {
                new("scan_whole", c => { for (int i = 0; i < c; i++) ladder.Step(); }, levels, "rung step"),
                new("one_replica_step", c => { for (int i = 0; i < c; i++) ladder.Replicas[i % levels].Step(); }, 1, "rung step"),
            });
        }
    }

    // ================================================================ sampler:smc

    internal static void SmcBlock()
    {
        const string block = "sampler:smc";
        foreach ((int particles, int dim, int terms) in new[] { (1000, 2, 0), (1000, 10, 0), (4000, 2, 0), (4000, 10, 0), (1000, 10, 2000), (4000, 10, 2000) })
        {
            var target = new GaussTarget(terms);
            var source = new GaussParticles(dim);
            var kernel = new RandomWalkKernel();
            var parts = new List<CostProfile.Part>();
            foreach ((int mutation, double resampleBelow) in new[] { (0, 0.0), (0, 1.0), (1, 0.0), (1, 1.0), (2, 0.0) })
            {
                var options = new SmcOptions { Particles = particles, MutationSteps = mutation, ResampleBelow = resampleBelow, MaxStages = 1_000_000, MaxDegreeOfParallelism = 1 };
                SmcRun<double[]> run = SmcRun<double[]>.Start(new SmcStage<double[]>(target, new Temperature(0.05, TemperingPath.Likelihood)), source, kernel, options, 17);
                for (int s = 1; s <= 15; s++) run.Advance(new SmcStage<double[]>(target, new Temperature(0.05 + 0.95 * s / 15.0, TemperingPath.Likelihood)));
                var stage = new SmcStage<double[]>(target, new Temperature(1.0, TemperingPath.Likelihood));
                int m = mutation;
                double r = resampleBelow;
                parts.Add(new($"advance_mutation{m}_resample{(r >= 1.0 ? "always" : "never")}", c => { for (int i = 0; i < c; i++) run.Advance(stage); }, particles, "particle"));
            }
            string name = $"N={particles} d={dim} {(terms == 0 ? "cheap" : $"expensive ({terms} terms)")}";
            CostProfile.Point(block, name, parts);

            // The pieces of a stage, timed alone beside the advances.
            var population = Enumerable.Range(0, particles).Select(i => source.Draw(new Xoshiro256PlusPlus(i + 1))).ToArray();
            var evaluations = population.Select(target.Evaluate).ToArray();
            var stage1 = new SmcStage<double[]>(target, new Temperature(1.0, TemperingPath.Likelihood));
            var weights = new double[particles];
            Array.Fill(weights, 1.0 / particles);
            var ancestors = new int[particles];
            var rng = new Xoshiro256PlusPlus(23);
            CostProfile.Point(block, name + " pieces", new List<CostProfile.Part>
            {
                new("evaluate_one_particle", c => { for (int i = 0; i < c; i++) ObjectSink = target.Evaluate(population[i % particles]); }, 1, "evaluation"),
                new("kernel_move_one_particle", c => { for (int i = 0; i < c; i++) ObjectSink = kernel.Move(population[i % particles], evaluations[i % particles], stage1, rng); }, 1, "move"),
                new("fork_one_particle", c => { for (int i = 0; i < c; i++) ObjectSink = source.Fork(population[i % particles]); }, 1, "fork"),
                new("kernel_prepare", c => { for (int i = 0; i < c; i++) kernel.Prepare(stage1, population, weights); }, (double)particles * dim, "element"),
                new("ancestors_multinomial", c => { for (int i = 0; i < c; i++) Resampling.Ancestors(ResamplingScheme.Multinomial, weights, rng, ancestors); }, particles, "particle"),
            });
        }

        // The fan-out: one mutation sweep at 1, 4 and 8 degrees of parallelism.
        {
            var target = new GaussTarget(2000);
            var source = new GaussParticles(10);
            var kernel = new RandomWalkKernel();
            var parts = new List<CostProfile.Part>();
            foreach (int dop in new[] { 1, 4, 8 })
            {
                var options = new SmcOptions { Particles = 4000, MutationSteps = 1, ResampleBelow = 1.0, MaxStages = 1_000_000, MaxDegreeOfParallelism = dop };
                SmcRun<double[]> run = SmcRun<double[]>.Start(new SmcStage<double[]>(target, new Temperature(0.05, TemperingPath.Likelihood)), source, kernel, options, 19);
                for (int s = 1; s <= 10; s++) run.Advance(new SmcStage<double[]>(target, new Temperature(0.05 + 0.95 * s / 10.0, TemperingPath.Likelihood)));
                var stage = new SmcStage<double[]>(target, new Temperature(1.0, TemperingPath.Likelihood));
                parts.Add(new($"advance_dop{dop}", c => { for (int i = 0; i < c; i++) run.Advance(stage); }, 4000, "particle"));
            }
            CostProfile.Point(block, "fan-out N=4000 d=10 expensive", parts);
        }
    }

    /// <summary>BASS's SMC kernel at the cost brief's base point: a move, the chain it builds, the account and its lock.</summary>
    internal static void BassSmcBlock()
    {
        const string block = "sampler:smc-bass";
        (BassModel model, _, _) = BassCost.Burned(BassCost.Base);
        var source = new BassParticles(model);
        var kernel = new BassSmcKernel(model);
        var stage = new SmcStage<BassParticle>(new BassParticleTarget(model.Target), new Temperature(0.6, TemperingPath.Likelihood));
        var rng = new Xoshiro256PlusPlus(41);
        BassParticle p = source.Draw(rng);
        ITargetEvaluation ev = stage.Target.Evaluate(p);
        for (int i = 0; i < 100; i++) ev = kernel.Move(p, ev, stage, rng).Evaluation;
        var chain = new ReversibleJumpChain<BassState>(p.Moves.Moves, model.Target, p.State, rng, stage.Temperature);
        chain.Step();
        IReadOnlyList<MoveStat> stats = chain.MoveStats();
        var gate = new Lock();
        var account = new BassSmcAccount(model.Prior);
        CostProfile.Point(block, "base point, one thread", new List<CostProfile.Part>
        {
            new("kernel_move", c => { for (int i = 0; i < c; i++) ev = kernel.Move(p, ev, stage, rng).Evaluation; }, 1, "move"),
            new("reversible_jump_chain_construction", c => { for (int i = 0; i < c; i++) ObjectSink = new ReversibleJumpChain<BassState>(p.Moves.Moves, model.Target, p.State, rng, stage.Temperature); }, 1, "chain"),
            new("chain_step", c => { for (int i = 0; i < c; i++) chain.Step(); }, 1, "step"),
            new("account_add", c => { for (int i = 0; i < c; i++) account.Add(stats, p.Moves); }, 1, "add"),
            new("account_add_under_uncontended_lock", c => { for (int i = 0; i < c; i++) lock (gate) account.Add(stats, p.Moves); }, 1, "add"),
            new("fork", c => { for (int i = 0; i < c; i++) ObjectSink = source.Fork(p); }, 1, "fork"),
        });

        // The lock under contention: T threads each adding the same stats; ns per add, wall over all adds.
        foreach (int threads in new[] { 1, 2, 4, 8 })
        {
            const int adds = 200_000;
            var times = new List<double>();
            for (int rep = 0; rep < 5; rep++)
            {
                var barrier = new Barrier(threads + 1);
                var workers = Enumerable.Range(0, threads).Select(_ => new Thread(() =>
                {
                    barrier.SignalAndWait();
                    for (int i = 0; i < adds; i++) lock (gate) account.Add(stats, p.Moves);
                    barrier.SignalAndWait();
                })).ToArray();
                foreach (Thread w in workers) w.Start();
                barrier.SignalAndWait();
                long t0 = Stopwatch.GetTimestamp();
                barrier.SignalAndWait();
                times.Add((Stopwatch.GetTimestamp() - t0) * 1e9 / Stopwatch.Frequency / ((double)adds * threads));
                foreach (Thread w in workers) w.Join();
            }
            CostProfile.Row(block, "account lock contention", $"threads={threads}", "ns_per_add_wall_over_all_adds", CostProfile.Median(times));
        }
    }

    // ================================================================ sampler:samc, wham, ladder, ensemble

    internal static void SamcBlock()
    {
        const string block = "sampler:samc";
        foreach (int nodes in new[] { 50, 200, 1000, 5000 })
            foreach (int chains in new[] { 4, 16 })
            {
                MulticanonicalWeights start = MulticanonicalWeights.Flat(-50.0, 0.0, 1.0, nodes);
                var learning = new SamcLearning(start, 1.0, 100);
                var rng = new Xoshiro256PlusPlus(3);
                var levels = new double[chains];
                for (int i = 0; i < chains; i++) levels[i] = -50.0 * rng.NextDouble();
                var target = new CountTarget();
                ITargetEvaluation inner = target.Evaluate(7);
                MulticanonicalWeights weights = learning.Weights;
                var mEval = new MulticanonicalEvaluation(inner, weights);
                var mTarget = new MulticanonicalTarget<int>(target, start);
                var plain = new ReversibleJumpChain<int>(CountMoves(), target, 4, new Xoshiro256PlusPlus(8));
                var multi = new ReversibleJumpChain<int>(CountMoves(), mTarget, 4, new Xoshiro256PlusPlus(8));
                var parts = new List<CostProfile.Part>
                {
                    new("learning_update", c => { for (int i = 0; i < c; i++) weights = learning.Update(levels); }, chains + nodes, "visit + node"),
                    new("weights_log_weight_lookup", c => { for (int i = 0; i < c; i++) Sink = weights.LogWeight(-50.0 * (i & 1023) / 1024.0); }, 1, "lookup"),
                    new("weights_try_locate", c => { for (int i = 0; i < c; i++) Sink = weights.TryLocate(-50.0 * (i & 1023) / 1024.0, out int seg, out double frac) ? seg + frac : 0; }, 1, "lookup"),
                    new("evaluation_reweighed", c => { for (int i = 0; i < c; i++) ObjectSink = mEval.Reweighed(i % 2 == 0 ? weights : start); }, 1, "reweigh"),
                    new("plain_chain_step", c => { for (int i = 0; i < c; i++) plain.Step(); }, 1, "step"),
                    new("multicanonical_chain_step", c => { for (int i = 0; i < c; i++) multi.Step(); }, 1, "step"),
                };
                CostProfile.Point(block, $"nodes={nodes} chains={chains}", parts);
            }
    }

    private static TemperingPilot SyntheticPilot(int rungs, int draws)
    {
        double[] b = Enumerable.Range(0, rungs).Select(i => Math.Pow(0.05, i / (double)(rungs - 1))).ToArray();
        b[0] = 1.0;
        var rng = new Xoshiro256PlusPlus(5 + rungs);
        double[][] levels = b.Select(bi => Enumerable.Range(0, draws).Select(_ => -20.0 / Math.Sqrt(bi) + 3.0 / Math.Sqrt(bi) * Normal.SampleStandard(rng)).ToArray()).ToArray();
        return new TemperingPilot(b, levels);
    }

    internal static void WhamBlock()
    {
        const string block = "sampler:wham";
        foreach ((int rungs, int draws) in new[] { (8, 2000), (16, 2000), (32, 2000), (8, 20_000), (16, 20_000), (32, 20_000), (8, 200_000) })
        {
            if (CostProfile.FreeMemoryGb() < 1.0) { CostProfile.Row(block, $"rungs={rungs} draws={draws}", "", "skipped", "free memory below 1 GB"); continue; }
            TemperingPilot pilot = SyntheticPilot(rungs, draws);
            WhamEstimate estimate = Wham.Estimate(pilot);
            CostProfile.Point(block, $"rungs={rungs} draws={draws}", new List<CostProfile.Part>
            {
                new("wham_estimate", c => { for (int i = 0; i < c; i++) ObjectSink = Wham.Estimate(pilot); }, (double)rungs * draws, "draw"),
            });
            GC.KeepAlive(estimate);
        }
    }

    internal static void LadderBlock()
    {
        const string block = "sampler:ladder";
        foreach ((int rungs, int draws) in new[] { (10, 2000), (30, 2000), (60, 2000), (10, 20_000), (30, 20_000), (10, 200_000) })
        {
            double[] b = Enumerable.Range(0, rungs).Select(i => (double)i / (rungs - 1)).ToArray();
            var rng = new Xoshiro256PlusPlus(9 + rungs);
            var ratios = new List<RungDraws?>();
            var slopes = new List<RungDraws>();
            for (int i = 0; i < rungs; i++)
            {
                double[][] chains = Enumerable.Range(0, 4).Select(_ => Enumerable.Range(0, draws / 4).Select(_ => -10.0 + 2.0 * Normal.SampleStandard(rng)).ToArray()).ToArray();
                double step = i + 1 < rungs ? b[i + 1] - b[i] : 0.0;
                ratios.Add(i + 1 < rungs ? RungDraws.FromChains(chains.Select(c => c.Select(v => step * v).ToArray()).ToArray()) : null);
                slopes.Add(RungDraws.FromChains(chains));
            }
            CostProfile.Point(block, $"rungs={rungs} draws={draws}", new List<CostProfile.Part>
            {
                new("stepping_stone", c => { for (int i = 0; i < c; i++) ObjectSink = LadderEvidence.SteppingStone(b, ratios); }, (double)rungs * draws, "draw"),
                new("path_sampling", c => { for (int i = 0; i < c; i++) ObjectSink = LadderEvidence.PathSampling(b, slopes); }, (double)rungs * draws, "draw"),
            });
        }
    }

    internal static void EnsembleBlock()
    {
        const string block = "sampler:ensemble";
        var rng = new Xoshiro256PlusPlus(4);
        foreach ((int draws, int functionals) in new[] { (1000, 1), (10_000, 1), (100_000, 1), (10_000, 10) })
        {
            double[][] chains = Enumerable.Range(0, 4).Select(_ => Enumerable.Range(0, draws).Select(_ => Normal.SampleStandard(rng)).ToArray()).ToArray();
            CostProfile.Point(block, $"diagnose 4 chains x {draws} draws x {functionals} functionals", new List<CostProfile.Part>
            {
                new("chain_diagnostics_diagnose", c => { for (int i = 0; i < c; i++) for (int f = 0; f < functionals; f++) ObjectSink = ChainDiagnostics.Diagnose(chains); }, 4.0 * draws * functionals, "draw"),
                new("chain_diagnostics_mean_ess", c => { for (int i = 0; i < c; i++) for (int f = 0; f < functionals; f++) ObjectSink = ChainDiagnostics.MeanEss(chains); }, 4.0 * draws * functionals, "draw"),
            });
        }
        var ensemble = new ChainEnsemble();
        var model = new IidModel();
        var iid = new IidChainRun(1);
        CostProfile.Point(block, "ensemble run 4 chains x 20000 samples", new List<CostProfile.Part>
        {
            new("ensemble_run_serial", c => { for (int i = 0; i < c; i++) ObjectSink = ensemble.Run(model, 4, 3, 200, 20_000, parallel: new System.Threading.Tasks.ParallelOptions { MaxDegreeOfParallelism = 1 }); }, 80_000, "draw"),
            new("ensemble_run_parallel", c => { for (int i = 0; i < c; i++) ObjectSink = ensemble.Run(model, 4, 3, 200, 20_000); }, 80_000, "draw"),
            new("iid_step_alone", c => { for (int i = 0; i < c; i++) Sink = iid.Step(); }, 1, "draw"),
        });
    }
}
