using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Maths.Rng;
using Repo.TestHarness;
using Xunit.Abstractions;

namespace Maths.Samplers.Tests;

/// <summary>One of a test's comparative assertions, which holds when <see cref="Lower"/> &lt; <see cref="Upper"/>.</summary>
internal readonly record struct Comparison(string Assertion, double Lower, double Upper)
{
    public bool Holds => Lower < Upper;
}

/// <summary>
/// A test's comparative assertions as a function of its seeds: <see cref="Measure"/>(dataSeed, chainSeeds) reruns the
/// test's computation and returns each comparison. <see cref="DataSeed"/> and <see cref="ChainSeeds"/> are the test's
/// own: the seed of its simulated data (null if it has none) and those of its samplers (empty for an exact computation).
/// </summary>
internal sealed record SweptClaim(string Test, int? DataSeed, int[] ChainSeeds, Func<int, int[], Comparison[]> Measure);

/// <summary>
/// The multi-seed sweep of comparative assertions (step 3.5). A test asserting "A beats B" at one seed passes by its
/// margin there; the sweep measures that margin across seeds on two axes:
/// <list type="bullet">
/// <item>the chain seeds alone, on the test's own data: whether the assertion is stable against the sampler's draws, so
/// that the test checks the method rather than one trajectory;</item>
/// <item>the data and chain seeds together: whether the claim holds for the data-generating process, as the test's
/// wording says, or only for its dataset.</item>
/// </list>
/// Each axis's seeds are derived by <see cref="SeedTree"/> from the test's own seeds, and each claim also runs once at
/// those, as a check that it reproduces its test. Guarded by the <c>SEED_SWEEP</c> environment variable, the number of
/// seeds per axis, so the default suite stays fast. The per-run table and the summary are written under
/// <c>artifacts/test-runs/SeedSweep/</c>.
/// </summary>
internal static class SeedSweep
{
    public static int? Seeds => int.TryParse(Environment.GetEnvironmentVariable("SEED_SWEEP"), out int n) && n > 0 ? n : null;

    public static void Run(string suite, IReadOnlyList<SweptClaim> claims, int seeds, ITestOutputHelper output)
    {
        var jobs = new List<(SweptClaim Claim, string Axis, int Seed, int DataSeed, int[] ChainSeeds)>();
        foreach (SweptClaim claim in claims)
        {
            jobs.Add((claim, "own", 0, claim.DataSeed ?? 0, claim.ChainSeeds));
            int m = claim.ChainSeeds.Length;
            // Two blocks of chain seeds, one per axis, from the test's first chain seed (or its data seed).
            int[] chain = m == 0 ? [] : SeedTree.Derive(claim.ChainSeeds[0], 2 * seeds * m);
            int[] data = claim.DataSeed is int d ? SeedTree.Derive(d, seeds) : [];
            for (int s = 0; s < seeds; s++)
            {
                if (m > 0) jobs.Add((claim, "chain", s + 1, claim.DataSeed ?? 0, chain[(s * m)..((s + 1) * m)]));
                if (claim.DataSeed is not null) jobs.Add((claim, "data", s + 1, data[s], chain[((seeds + s) * m)..((seeds + s + 1) * m)]));
            }
        }

        var results = new Comparison[jobs.Count][];
        Parallel.For(0, jobs.Count, j => results[j] = jobs[j].Claim.Measure(jobs[j].DataSeed, jobs[j].ChainSeeds));

        static string F(double v) => v.ToString("R", CultureInfo.InvariantCulture);
        var table = new StringBuilder("test\taxis\tseed\tdata_seed\tchain_seeds\tassertion\tlower\tupper\tholds\n");
        var groups = new Dictionary<(string Test, string Assertion, string Axis), List<Comparison>>();
        var order = new List<(string, string, string)>();
        for (int j = 0; j < jobs.Count; j++)
            foreach (Comparison c in results[j])
            {
                (SweptClaim claim, string axis, int seed, int dataSeed, int[] chainSeeds) = jobs[j];
                table.Append(claim.Test).Append('\t').Append(axis).Append('\t').Append(seed).Append('\t').Append(dataSeed).Append('\t')
                     .Append(string.Join(",", chainSeeds)).Append('\t').Append(c.Assertion).Append('\t').Append(F(c.Lower)).Append('\t')
                     .Append(F(c.Upper)).Append('\t').Append(c.Holds).Append('\n');
                var key = (claim.Test, c.Assertion, axis);
                if (!groups.TryGetValue(key, out List<Comparison>? list))
                {
                    groups[key] = list = [];
                    order.Add(key);
                }
                list.Add(c);
            }

        // Per assertion and axis: how often it holds, its worst margin, and the median and worst ratio of its sides
        // where both are positive (an MSE or an ESS).
        var summary = new StringBuilder("test\tassertion\taxis\tholds\truns\tworst_margin\tmedian_ratio\tworst_ratio\n");
        foreach ((string test, string assertion, string axis) in order)
        {
            List<Comparison> list = groups[(test, assertion, axis)];
            int holds = list.Count(c => c.Holds);
            double worstMargin = list.Min(c => c.Upper - c.Lower);
            double[]? ratios = list.All(c => c.Lower > 0.0 && c.Upper > 0.0) ? list.Select(c => c.Lower / c.Upper).Order().ToArray() : null;
            summary.Append(test).Append('\t').Append(assertion).Append('\t').Append(axis).Append('\t').Append(holds).Append('\t')
                   .Append(list.Count).Append('\t').Append(F(worstMargin)).Append('\t')
                   .Append(ratios is null ? "" : F(ratios[ratios.Length / 2])).Append('\t')
                   .Append(ratios is null ? "" : F(ratios[^1])).Append('\n');
            output.WriteLine($"[seed-sweep] {holds,2}/{list.Count,-2} {axis,-5} {test}: {assertion}; worst margin {worstMargin:G3}" +
                             (ratios is null ? "" : $", ratio median {ratios[ratios.Length / 2]:F3}, worst {ratios[^1]:F3}"));
        }

        ArtifactRun run = HarnessArtifacts.Create("test-runs", "SeedSweep", suite,
                                                  new Dictionary<string, object?> { ["Seeds"] = seeds, ["Claims"] = claims.Count });
        output.WriteLine($"[seed-sweep] table {run.WriteText("seed-sweep.tsv", table.ToString())}");
        output.WriteLine($"[seed-sweep] summary {run.WriteText("summary.tsv", summary.ToString())}");
    }
}
