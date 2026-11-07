using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// KHAOS's primitives in khaos 2.1.1 (<c>khaos_primitives_oracle.R</c>; E19's build brief, slice 0), each against a reference computed
/// on both sides:
/// <list type="bullet">
/// <item><b>The standardized shifted Legendre polynomials.</b> <c>ss_legendre_poly</c>'s values equal Bonnet's recurrence here, are
/// orthonormal under the engine's <see cref="GaussLegendre"/> rule and R's, and ψ_J vanishes at the J-point nodes. Its recursion is
/// exponential in the degree, a cost the engine's lanes do not share.</item>
/// <item><b>|𝒜|.</b> <c>A_size</c> equals the sum here in exact integers, and where small enough the multi-indices walked here, R's
/// brute count and the package's own enumerator.</item>
/// <item><b>The coins.</b> <c>make_weights</c>' probabilities lie in (0, 1), rise with the term counts and sum to q₀ only to its three
/// passes' accuracy; the exact root sums to q₀.</item>
/// <item><b>random_partition's law.</b> The law found by enumerating its draws equals the sorted-cuts count here, and its own draws
/// carry it.</item>
/// </list>
/// </summary>
public sealed class KhaosPrimitiveParityTests(ITestOutputHelper output)
{
    private const int Draws = 40_000;

    // kind, then parameters: 1 (jmax); 2 (p, d, q); 3 (q0, η…); 4 (d, q).
    private static readonly double[][] Cases =
    [
        [1, 8], [1, 24],
        [2, 3, 4, 2], [2, 4, 6, 3], [2, 5, 5, 5], [2, 6, 3, 6], [2, 3, 4, 5], [2, 10, 8, 3], [2, 20, 15, 5], [2, 30, 15, 5],
        [3, 2, 3, 0, 1, 0], [3, 3, 0, 0, 0, 0], [3, 1, 5, 0, 0, 2, 1, 0], [3, 3, 5, 0, 0, 2, 1, 0], [3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12],
        [4, 5, 3], [4, 6, 3], [4, 7, 4], [4, 4, 2], [4, 6, 1],
    ];

    private static readonly Lazy<JsonElement?> Fixture = new(() => GbassOracle.Run("khaos_primitives_oracle.R", KhaosOracle.Padded(Cases), Draws.ToString(), "20261107"));

    [Fact]
    public void Legendre_ValuesMatchTheRecurrence_AndAreOrthonormal()
    {
        if (Fixture.Value is not { } r) return;
        int checkedCases = 0;
        foreach (JsonElement c in r.GetProperty("legendre").EnumerateArray())
        {
            int jmax = c.GetProperty("jmax").GetInt32();
            double[] x = GbassOracle.Numbers(c.GetProperty("x"));
            double[][] values = [.. c.GetProperty("values").EnumerateArray().Select(GbassOracle.Numbers)];
            double worst = 0;
            for (int j = 0; j <= jmax; j++)
                for (int i = 0; i < x.Length; i++)
                    worst = Math.Max(worst, Math.Abs(values[j][i] - KhaosOracle.Psi(j, x[i])) / Math.Sqrt(2 * j + 1));
            Assert.True(worst < 1e-12, $"ss_legendre_poly against the recurrence to degree {jmax}: {worst:G3}");

            // The engine's rule mapped to [0, 1] is R's (Golub–Welsch there, Newton here), and the recurrence is orthonormal under it.
            GaussLegendre rule = GaussLegendre.Of(jmax + 1);
            double[] nodes = GbassOracle.Numbers(c.GetProperty("nodes"));
            for (int k = 0; k < rule.Count; k++) Assert.Equal(nodes[k], 0.5 * (rule.Nodes[k] + 1.0), 13);
            double gramError = 0;
            for (int a = 0; a <= jmax; a++)
                for (int b = a; b <= jmax; b++)
                {
                    double sum = 0;
                    for (int k = 0; k < rule.Count; k++)
                    {
                        double u = 0.5 * (rule.Nodes[k] + 1.0);
                        sum += 0.5 * rule.Weights[k] * KhaosOracle.Psi(a, u) * KhaosOracle.Psi(b, u);
                    }
                    gramError = Math.Max(gramError, Math.Abs(sum - (a == b ? 1.0 : 0.0)));
                }
            double rDiag = c.GetProperty("gram_diag_error").GetDouble(), rOff = c.GetProperty("gram_offdiag_max").GetDouble();
            double zeros = c.GetProperty("zeros_max").GetDouble(), seconds = c.GetProperty("seconds").GetDouble();
            output.WriteLine($"ψ to degree {jmax}: |R − recurrence| {worst:G3}; Gram error here {gramError:G3}, in R {Math.Max(rDiag, rOff):G3}; " +
                             $"ψ_{jmax} at its nodes {zeros:G3}; one ss_legendre_poly call at degree {jmax}: {seconds:F2} s");
            Assert.True(gramError < 1e-12 && rDiag < 1e-12 && rOff < 1e-12 && zeros < 1e-12);
            checkedCases++;
        }
        Assert.Equal(Cases.Count(c => c[0] == 1), checkedCases);
    }

    [Fact]
    public void ASize_IsTheCount_OfTheMultiIndices()
    {
        if (Fixture.Value is not { } r) return;
        int checkedCases = 0;
        foreach (JsonElement c in r.GetProperty("a_size").EnumerateArray())
        {
            int p = c.GetProperty("p").GetInt32(), d = c.GetProperty("d").GetInt32(), q = c.GetProperty("q").GetInt32();
            long bySum = KhaosOracle.CountBySum(p, d, q);
            Assert.Equal(bySum, (long)c.GetProperty("code").GetDouble());
            if (bySum <= 60_000_000) Assert.Equal(bySum, KhaosOracle.CountByEnumeration(p, d, q));
            foreach (string k in new[] { "generate", "brute" })
                if (c.GetProperty(k).ValueKind != JsonValueKind.Null) Assert.Equal(bySum, (long)c.GetProperty(k).GetDouble());
            output.WriteLine($"|𝒜({p}, {d}, {q})| = {bySum}");
            checkedCases++;
        }
        Assert.Equal(Cases.Count(c => c[0] == 2), checkedCases);
    }

    [Fact]
    public void MakeWeights_AreCoinProbabilities_SummingToQ0OnlyApproximately()
    {
        if (Fixture.Value is not { } r) return;
        int checkedCases = 0;
        double worstGap = 0;
        foreach (JsonElement c in r.GetProperty("weights").EnumerateArray())
        {
            double q0 = c.GetProperty("q0").GetDouble();
            double[] eta = GbassOracle.Numbers(c.GetProperty("eta")), code = GbassOracle.Numbers(c.GetProperty("code")), exact = GbassOracle.Numbers(c.GetProperty("exact"));
            int p = eta.Length;
            Assert.All(code, w => Assert.InRange(w, 1e-300, 1.0 - 1e-15));
            for (int i = 0; i < p; i++)
                for (int j = 0; j < p; j++)
                    if (eta[i] > eta[j]) Assert.True(code[i] > code[j], $"weights at η {string.Join(",", eta)} do not rise with the counts");
            if (eta.Distinct().Count() == 1) Assert.All(code, w => Assert.Equal(q0 / p, w, 12));
            // The exact root here, by bisection on Σ v^β = q0 over the same v.
            double[] v = [.. eta.Select(e => (e + 2.0) / eta.Sum(u => u + 2.0) * q0 / p)];
            double lo = 1e-12, hi = 1e3;
            for (int k = 0; k < 200; k++) { double mid = 0.5 * (lo + hi); if (v.Sum(u => Math.Pow(u, mid)) > q0) lo = mid; else hi = mid; }
            Assert.Equal(c.GetProperty("beta_exact").GetDouble(), 0.5 * (lo + hi), 9);
            Assert.Equal(q0, exact.Sum(), 10);
            double gap = c.GetProperty("code_sum").GetDouble() - q0;
            worstGap = Math.Max(worstGap, Math.Abs(gap));
            output.WriteLine($"make_weights(η = {string.Join(",", eta)}, q0 = {q0}): Σw − q0 = {gap:G3}");
            checkedCases++;
        }
        Assert.Equal(Cases.Count(c => c[0] == 3), checkedCases);
        Assert.True(worstGap > 0, "make_weights sums to q0 exactly in every case");
    }

    [Fact]
    public void RandomPartition_LawIsTheSortedCutsCount_AndItsDrawsCarryIt()
    {
        if (Fixture.Value is not { } r) return;
        int checkedCases = 0;
        foreach (JsonElement c in r.GetProperty("partition").EnumerateArray())
        {
            int d = c.GetProperty("d").GetInt32(), q = c.GetProperty("q").GetInt32(), n = c.GetProperty("n").GetInt32();
            List<int[]> mine = KhaosOracle.Compositions(d, q);
            int[][] theirs = [.. c.GetProperty("compositions").EnumerateArray().Select(e => GbassOracle.Numbers(e).Select(v => (int)v).ToArray())];
            Assert.Equal(mine.Count, theirs.Length);
            double[] law = GbassOracle.Numbers(c.GetProperty("law")), freq = GbassOracle.Numbers(c.GetProperty("freq"));
            double worstZ = 0;
            for (int k = 0; k < mine.Count; k++)
            {
                Assert.Equal(mine[k], theirs[k]);
                Assert.Equal(KhaosOracle.SortedCutsProbability(mine[k]), law[k], 14);
                if (law[k] < 1) worstZ = Math.Max(worstZ, Math.Abs(freq[k] - law[k]) / Math.Sqrt(law[k] * (1 - law[k]) / n));
            }
            Assert.Equal(1.0, law.Sum(), 12);
            output.WriteLine($"random_partition({d}, {q}): {mine.Count} compositions, law from {law.Min():G4} to {law.Max():G4} " +
                             $"(uniform {1.0 / KhaosOracle.Choose(d - 1, q - 1):G4}); its draws' worst z {worstZ:F2}");
            Assert.True(worstZ < 5.0);
            checkedCases++;
        }
        Assert.Equal(Cases.Count(c => c[0] == 4), checkedCases);
    }
}
