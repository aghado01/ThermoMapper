using System;
using System.Collections.Generic;
using Maths.Regression.Spline.Bass;
using Maths.Regression.Spline.Bass.Sensitivity;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

public sealed class MultiOutputReadoutCostTests(ITestOutputHelper output)
{
    [Fact]
    public void JointReadoutsAtFourHundredDraws()
    {
        if (!BassQualification.Enabled) return;
        // The first complete pass warms shared readout/JIT paths before any form is reported.
        for (int pass = 0; pass < 2; pass++)
        foreach (string form in new[] { "g", "ridge", "modified" })
        {
            var (fit, data) = MultiOutputReadoutTests.Fixed(form, count: 400);
            var cases = MultiOutputReadoutTests.Cases(); object? sink = null;
            var timers = new List<(string Name, Func<int, double> Batch)>();
            void Add(string name, Action action) => timers.Add((name, BassCost.Batch(BassCost.Meter.Time, action, batch: 10)));
            Add("joint draws", () => sink = BassMultiOutput.Predict(fit, data, cases));
            Add("matrix-t components", () => sink = BassMultiOutput.Conditionals(fit, data, cases));
            Add("covariance summary", () => sink = BassMultiOutput.Covariance(fit, data));
            Add("projected Sobol", () => sink = BassMultiOutput.Sensitivity(fit, data, new double[] { 1, -0.4 }));
            Add("co-active exact", () => sink = BassMultiOutput.CoActive(fit, data, new CoActiveSpec { Output = new double[] { 1, 0 }, OtherOutput = new double[] { 0, 1 } }));
            var readings = BassCost.Interleave(timers, repeats: 5);
            if (pass == 0) continue;
            foreach (var timer in timers) output.WriteLine($"{form}, N=40, M=3, K=2, draws=400, cases=2, {timer.Name}: {BassCost.Median(readings[timer.Name]):F2} μs, IQR {BassCost.Iqr(readings[timer.Name]):F2}");
            output.WriteLine($"clock={(BassCost.Cycles ? "thread cycles" : "wall")}, indicative; other test hosts={string.Join(", ", BassCost.OtherTestHosts())}");
            GC.KeepAlive(sink);
        }
    }
}
