using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Blocks;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The efficiency pass's lever 5, the split before any change (<c>efficiency-pass-brief.md</c>): what τ's slice step and the composition are
/// made of at a burned state. One τ evaluation is a state (<c>WithLogTau</c>), an evaluation object (<c>Evaluate</c>: the prior's terms,
/// the refusal), and the marginal's arithmetic (<c>Marginal</c> then <c>LogMarginal</c>: the logs and the log-gamma); the step adds the slice's own
/// stepping and a closure per evaluation. Each part is timed on its own, in batches of calls round-robin with the others', the median of
/// nine, and set against the whole block. The composition is split the same way: the coefficient block's construction, the conditional at
/// (τ, 1), and the draw. Behind <c>BASS_QUALIFICATION</c> = 1; rows go to <c>BASS_COST_OUT</c>; times wait for a quiet machine.
/// </summary>
public sealed class BassTauSplitTests
{
    private readonly ITestOutputHelper _out;
    public BassTauSplitTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Split_AtTheBurnedState()
    {
        if (!BassQualification.Enabled) return;
        var rows = new BassCost.Rows("tau-split", "split");
        foreach ((string q, string v) in BassCost.Machine()) rows.Add("machine", q, v);
        var load = new BassCost.Load();
        foreach (BassCost.Point x in new[] { BassCost.Base, BassCost.Base with { N = 200 }, BassCost.Base with { Order = 5 }, BassCost.Base with { P = 200 } })
        {
            (BassModel model, BassChain chain, _) = BassCost.Burned(x);
            BassState state = chain.Current;
            BassTarget target = model.Target;
            BassDesign design = chain.Structure.DesignAt(state);
            BassChainSnapshot fixedChain = chain.Snapshot();
            var counting = new Counting(model.Target);
            BassChain moved = BassChain.Restore(model, fixedChain, counting);
            var tau = new BassTauBlock();
            var rng = new Xoshiro256PlusPlus(seed: 913);
            var table = new BassTermTable();
            moved.Compose(rng, table);
            long updates = 0;
            double sink = 0.0;
            int k = 0;
            double u0 = state.LogTau;
            Temperature temperature = Temperature.Posterior;
            BassState[] states = [.. Enumerable.Range(0, 64).Select(i => state.WithLogTau(u0 + 0.02 * (i - 32)))];
            BassEvaluation[] evaluations = [.. states.Select(target.Evaluate)];
            GPriorConditional conditional = target.CoefficientBlock(design)(state.Tau, 1.0);
            var draw = new double[conditional.Columns];

            var timers = new List<(string Name, Func<int, double> Batch)>
            {
                ("tau block", Timer(200, () => { tau.Update(moved.Moves, rng); updates++; return 0.0; })),
                ("with log tau", Timer(2000, () => states[k++ & 63 ^ 1].WithLogTau(u0 + 1e-3 * k).LogTau)),
                ("evaluate", Timer(2000, () => target.Evaluate(states[k++ & 63]).Tau)),
                ("evaluate and density", Timer(2000, () => temperature.LogDensity(target.Evaluate(states[k++ & 63])))),
                ("marginal only", Timer(2000, () => target.Marginal(states[k++ & 63].Fit, states[k & 63].Tau, 1.0).LogMarginal(evaluations[k & 63].Sigma2))),
                ("slice, constant density", Timer(2000, () =>
                {
                    double u = u0, d = 0.0;
                    SliceSampler.Step(ref u, ref d, BassTauBlock.Width, BassTauBlock.MaxSteps, v => -0.5 * (v - u0) * (v - u0), rng);
                    return u;
                })),
                ("compose", Timer(500, () => moved.Compose(rng, table).Tau)),
                ("coefficient block", Timer(2000, () => target.CoefficientBlock(design).GetHashCode())),
                ("conditional at (tau, 1)", Timer(2000, () => target.CoefficientBlock(design)(state.Tau, 1.0).Columns)),
            };
            Dictionary<string, double[]> times = BassCost.Interleave(timers);
            double evaluationsPerStep = (double)counting.Evaluations / Math.Max(1, updates);
            string key = x.Key;
            rows.Add(key, "M", state.Count);
            rows.Add(key, "τ evaluations per step", evaluationsPerStep);
            foreach ((string name, double[] values) in times)
            {
                rows.Add(key, $"{name} µs", BassCost.Median(values));
                rows.Add(key, $"{name} µs IQR", BassCost.Iqr(values));
            }
            double oneEvaluation = BassCost.Median(times["evaluate and density"]) + BassCost.Median(times["with log tau"]);
            rows.Add(key, "evaluations' share of the τ block (with + evaluate + density, per step)", oneEvaluation * evaluationsPerStep / BassCost.Median(times["tau block"]));
            rows.Add(key, "marginal's share of one evaluation (with + evaluate + density)", BassCost.Median(times["marginal only"]) / oneEvaluation);
            _out.WriteLine($"{key}: M {state.Count}, {evaluationsPerStep:F1} evaluations per step; τ block {BassCost.Median(times["tau block"]):F2} µs = " +
                           $"{evaluationsPerStep:F1} × (with {BassCost.Median(times["with log tau"]) * 1000:F0} ns + evaluate/density {BassCost.Median(times["evaluate and density"]) * 1000:F0} ns, " +
                           $"of which the marginal {BassCost.Median(times["marginal only"]) * 1000:F0} ns) + the slice {BassCost.Median(times["slice, constant density"]):F2} µs; " +
                           $"compose {BassCost.Median(times["compose"]):F2} µs (block {BassCost.Median(times["coefficient block"]):F2}, conditional {BassCost.Median(times["conditional at (tau, 1)"]):F2})");
            Assert.True(double.IsFinite(sink));
            rows.Flush();
        }
        (string[] end, double cpu, double seconds, double others) = load.End();
        rows.Add("load", "other test hosts at the start", string.Join(" ", load.AtStart));
        rows.Add("load", "other test hosts at the end", string.Join(" ", end));
        rows.Add("load", "other processes' cores busy", others);
        rows.Add("load", "flagged", load.AtStart.Length > 0 || end.Length > 0 ? "yes" : "no");
        _out.WriteLine($"load: others {others:F2} cores; rows {rows.Flush()}");
    }

    private static Func<int, double> Timer(int calls, Func<double> call)
    {
        double sink = 0.0;
        return [MethodImpl(MethodImplOptions.AggressiveOptimization)] (int _) =>
        {
            long t0 = System.Diagnostics.Stopwatch.GetTimestamp();
            for (int c = 0; c < calls; c++) sink += call();
            return (System.Diagnostics.Stopwatch.GetTimestamp() - t0) * 1e6 / System.Diagnostics.Stopwatch.Frequency / calls + sink * 0.0;
        };
    }

    private sealed class Counting(BassTarget target) : ISamplingTarget<BassState>
    {
        public long Evaluations { get; private set; }

        public ITargetEvaluation Evaluate(BassState state)
        {
            Evaluations++;
            return target.Evaluate(state);
        }
    }
}
