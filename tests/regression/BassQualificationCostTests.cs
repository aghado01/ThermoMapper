using System;
using System.Diagnostics;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// S10's cost baseline: an untempered run's wall time per sweep by the observations n, the inputs p and the terms M, one
/// thread. A sweep's move costs O(N·M + M²) for the bordered update and τ's slice O(1) per step given the fit, and each
/// retained draw's composition O(N·M + M³), so the time should grow about linearly in n and in M at these sizes.
/// <list type="bullet">
///   <item><b>Data:</b> 10 sin(πx₁x₂) + 10x₂² + 3 sin(6πx₁) + noise sd 0.05 over p uniform inputs (the rest inert),
///   n ∈ {200, 1,000, 5,000}, p ∈ {2, 10}: a surface that asks for many terms.</item>
///   <item><b>M:</b> capped at M_max ∈ {5, 20, 50} under the default count prior, and read back as the draws' mean,
///   which is the cap where the data ask for more. A count prior pushed toward the cap was tried first: its exact
///   starts drew collections of about M_max terms, which on 200 points almost never meet the support rule, and the
///   prior sampler refused, as it should.</item>
///   <item><b>Runs:</b> two chains of 500 burn-in and 500 draws, thin 1, at one degree of parallelism, after one
///   untimed warm-up fit.</item>
/// </list>
/// It reads, it does not judge: a qualification run (<see cref="BassQualification.Enabled"/>), whose table the fitting
/// checklist records. Run it alone; other load on the machine moves the times.
/// </summary>
public sealed class BassQualificationCostTests
{
    private readonly ITestOutputHelper _out;
    public BassQualificationCostTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void TimePerSweep_ByObservationsInputsAndTerms()
    {
        if (!BassQualification.Enabled) return;   // heavy
        var run = new BassRunOptions { MaxDegreeOfParallelism = 1 };
        Bass.Fit(Config(5, 100), Data(200, 2), run);   // warm-up: JIT and caches
        _out.WriteLine("n | p | M_max | mean M | µs per sweep");
        foreach (int n in new[] { 200, 1000, 5000 })
            foreach (int p in new[] { 2, 10 })
                foreach (int cap in new[] { 5, 20, 50 })
                {
                    BassData data = Data(n, p);
                    BassConfig config = Config(cap, 500);
                    var clock = Stopwatch.StartNew();
                    BassFit fit = Bass.Fit(config, data, run);
                    double micro = clock.Elapsed.TotalMilliseconds * 1000.0 / (2 * 1000);
                    double meanM = fit.Draws.SelectMany(c => c).Average(d => d.Structure.TermIds.Length);
                    Assert.True(double.IsFinite(micro));
                    _out.WriteLine($"{n} | {p} | {cap} | {meanM:F1} | {micro:F0}");
                }
    }

    private static BassConfig Config(int cap, int length) => new()
    {
        Structure = new BassStructureSpec { MaxBasis = cap },
        Sampling = new BassSampling { Chains = 2, MasterSeed = 71, BurnIn = length, Samples = length },
    };

    private static BassData Data(int n, int p)
    {
        var rng = new Xoshiro256PlusPlus(seed: 73);
        var x = new double[p][];
        for (int v = 0; v < p; v++) x[v] = new double[n];
        var y = new double[n];
        for (int i = 0; i < n; i++)
        {
            for (int v = 0; v < p; v++) x[v][i] = rng.NextDouble();
            y[i] = 10 * Math.Sin(Math.PI * x[0][i] * x[1][i]) + 10 * x[1][i] * x[1][i] + 3 * Math.Sin(6 * Math.PI * x[0][i])
                   + 0.05 * BassQualification.Gaussian(rng);
        }
        return new BassData([.. Enumerable.Range(0, p).Select(v => (BassPredictor)new ContinuousPredictor($"x{v + 1}", x[v]))], new ScalarResponse(y));
    }
}
