using System;
using System.Diagnostics;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// F-M2's Phase 1 smoke reading (the local-move options brief): exploration, not the pass, so it recommends nothing. S10's
/// three fixtures at S10's lengths, each with its S10 master seed and four chains, untempered: the options off (S10's
/// chains), the local steps, the fixed-M family, and all of them; and Friedman tempered, self-sizing, without and with the
/// split rule.
/// <list type="bullet">
/// <item>Per move kind: the acceptance of the attempts and of the admitted proposals, and the null and refused shares.</item>
/// <item>Each chain's band of M over its sampling phase, and M's R̂ and bulk ESS as context.</item>
/// <item>The time per sweep, the fit run on one thread; tempered, per scan and per replica sweep.</item>
/// <item>Tempered: each chain's settled count, its steep pairs (rejection above 0.9), the ideal and observed round-trip
/// rates, and its rungs in [0.10, 0.40] with the summed rejection of the pairs there, beside Phase 0's Δb = 0.01 profile
/// (Σr 4.7–5.2 over that window).</item>
/// </list>
/// It asserts only its own soundness. A qualification run (<see cref="BassQualification.Enabled"/>).
/// </summary>
public sealed class BassLocalMoveSmokeTests
{
    private readonly ITestOutputHelper _out;
    public BassLocalMoveSmokeTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(BassMixing.Friedman)]
    [InlineData(BassMixing.Mixed)]
    [InlineData(BassMixing.Functional)]
    public void Untempered(string name)
    {
        if (!BassQualification.Enabled) return;   // heavy
        foreach (string options in new[] { "none", "steps", "family", "all" })
        {
            BassMixing.Fixture f = BassMixing.Get(name);
            Read($"{name}, {options}", f.Config with { Sampling = BassLocalStepTests.OptionsFor(options)(f.Config.Sampling) }, f.Data);
        }
    }

    [Fact]
    public void FriedmanTempered()
    {
        if (!BassQualification.Enabled) return;   // heavy
        foreach (SplitRule? split in new[] { null, new SplitRule() })
        {
            BassMixing.Fixture f = BassMixing.Get(BassMixing.Friedman, new TemperingSpec { Split = split });
            Read($"{BassMixing.Friedman}, tempered{(split is null ? "" : ", 4b")}", f.Config, f.Data);
        }
    }

    private void Read(string label, BassConfig config, BassData data)
    {
        var watch = Stopwatch.StartNew();
        BassFit fit = Bass.Fit(config, data, new BassRunOptions { MaxDegreeOfParallelism = 1 });
        watch.Stop();
        BassSampling s = config.Sampling;
        BassDiagnostics d = fit.Diagnostics;
        double sweeps = s.Chains * (double)(s.BurnIn + s.Samples * s.Thin);
        double perSweep = watch.Elapsed.TotalMicroseconds / sweeps;
        string time = d.Tempering is { } t
            ? $"{perSweep:F0} µs per scan, {perSweep / t.Levels.Average():F1} µs per replica sweep"
            : $"{perSweep:F1} µs per sweep";
        _out.WriteLine($"{label}: {watch.Elapsed.TotalSeconds:F1} s, {time}; M's R̂ {d.Count.RHat}, bulk ESS {d.Count.BulkEss}");
        for (int c = 0; c < fit.Draws.Count; c++)
        {
            int[] m = [.. fit.Draws[c].Select(draw => draw.Structure.TermIds.Length).Order()];
            Assert.Equal(s.Samples, m.Length);
            _out.WriteLine($"  chain {c}: M in [{m[0]}, {m[^1]}], median {m[m.Length / 2]}");
        }
        foreach (BassKindCount k in d.MoveKinds.Entries.Where(k => k.Attempts > 0))
        {
            long admitted = k.Attempts - k.Null - k.Refused;
            Assert.True(admitted >= 0 && k.Rejected >= 0);
            _out.WriteLine($"  {k.Kind}: {k.Attempts} attempts, α {k.AcceptanceRate:F4}, α admitted {(admitted > 0 ? (double)k.Accepted / admitted : double.NaN):F4}, " +
                           $"null {(double)k.Null / k.Attempts:F3}, refused {(double)k.Refused / k.Attempts:F3}");
        }
        if (d.Tempering is not { } r) return;
        _out.WriteLine($"  ladder: barrier {r.Barrier:F2}, max rejection {r.MaxPairRejection:F3}, ideal round trips {r.IdealRoundTripRate:F4}/scan, " +
                       $"observed {r.RoundTripsPerScan:F4}/scan, at max levels {r.AtMaxLevels}");
        for (int c = 0; c < r.Ladders.Length; c++)
        {
            double[] b = r.Ladders[c], rejection = r.PairRejection[c];
            int steep = rejection.Count(x => x > 0.9), inWindow = b.Count(x => x >= 0.10 && x <= 0.40);
            double windowSum = Enumerable.Range(0, rejection.Length).Where(i => b[i] <= 0.40 && b[i + 1] >= 0.10).Sum(i => rejection[i]);
            _out.WriteLine($"  chain {c}: {r.Levels[c]} rungs, {steep} steep pairs, {inWindow} rungs in [0.10, 0.40] with Σr {windowSum:F2} there; " +
                           $"rejections [{string.Join(", ", rejection.Select(x => x.ToString("F2")))}]");
        }
    }
}
