using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Primitives;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The accumulator's outputs at fixed seeds, recorded at <c>6233f93</c> before the per-draw retention landed: the
/// running sums and RNG words as exact bit patterns, and a SHA-256 of every array (the histogram, the spins, the per-edge
/// counts, the landscapes). Retention reads what the cycle already computes, so it must leave every case bit-identical:
/// a change here means the RNG stream or an existing output moved.
/// </summary>
public sealed class SwAccumulatorProbeTests
{
    [Theory]
    [MemberData(nameof(CaseNames))]
    public void Outputs_MatchTheRecord(string name)
    {
        ProbeCase probe = Cases.Single(c => c.Name == name);

        Assert.Equal(Recorded[name], Fingerprint(Run(probe)));
    }

    public static TheoryData<string> CaseNames()
    {
        var data = new TheoryData<string>();
        foreach (ProbeCase probe in Cases) data.Add(probe.Name);
        return data;
    }

    // ── The record (6233f93) ─────────────────────────────────────────────────

    private static readonly Dictionary<string, string> Recorded = new(StringComparer.Ordinal)
    {
        ["triangle/T0.3/none"] = // m̄ = 0.999
            "draws=400 sums=40AC100000000000,4000000000000000,C092B80000000000,40AC100000000000,4078F8E38E38E38E,4078F4F0329161FA rng=D6D78BA709EC0AEA,72BF68C25DE2AB66,300374B953CC375F,1A7DA408F21CF3F8 hist=D58DFABB87CBECF6 spins=636952D3023D8CF5 bond=- agree=- comember=- sizeLandscape=- giantLandscape=-",
        ["triangle/T0.3/full"] = // m̄ = 0.999
            "draws=400 sums=40AC100000000000,4000000000000000,C092B80000000000,40AC100000000000,4078F8E38E38E38E,4078F4F0329161FA rng=D6D78BA709EC0AEA,72BF68C25DE2AB66,300374B953CC375F,1A7DA408F21CF3F8 hist=D58DFABB87CBECF6 spins=636952D3023D8CF5 bond=5DEC186B6E3E39C8 agree=9B0A13BA95CCB792 comember=4D8F29C3F116C982 sizeLandscape=07DE354BB968686B giantLandscape=F8387A635AFBE340",
        ["triangle/T1/none"] = // m̄ = 0.704
            "draws=400 sums=40A2100000000000,4077800000000000,C086400000000000,409BD00000000000,40719C71C71C71AE,406D2F684BDA12DB rng=EC032A6DBE12586C,09989BA04AB4CEC1,31411F94AD263409,6751C8430734FDC8 hist=681C6CB09BE8DDD3 spins=94D8967691A2A042 bond=- agree=- comember=- sizeLandscape=- giantLandscape=-",
        ["triangle/T1/full"] = // m̄ = 0.704
            "draws=400 sums=40A2100000000000,4077800000000000,C086400000000000,409BD00000000000,40719C71C71C71AE,406D2F684BDA12DB rng=EC032A6DBE12586C,09989BA04AB4CEC1,31411F94AD263409,6751C8430734FDC8 hist=681C6CB09BE8DDD3 spins=94D8967691A2A042 bond=364777013AD445E5 agree=46B51F4EAA675E83 comember=80A9A29CE43D334F sizeLandscape=F18F73C3F2629202 giantLandscape=323146E7CB16E8B2",
        ["triangle/T3/none"] = // m̄ = 0.507
            "draws=400 sums=4097280000000000,4084F80000000000,C07A400000000000,4088C00000000000,406955555555553E,4060D097B425ED0B rng=F8C55045F1D27922,23A1DAEC9E7DEB0A,2CB28AC0FBC1CBA4,103DA4FC265E823F hist=6F7427E3DDCB5A49 spins=87FA498592C87CCE bond=- agree=- comember=- sizeLandscape=- giantLandscape=-",
        ["triangle/T3/full"] = // m̄ = 0.507
            "draws=400 sums=4097280000000000,4084F80000000000,C07A400000000000,4088C00000000000,406955555555553E,4060D097B425ED0B rng=F8C55045F1D27922,23A1DAEC9E7DEB0A,2CB28AC0FBC1CBA4,103DA4FC265E823F hist=6F7427E3DDCB5A49 spins=87FA498592C87CCE bond=1DA9C65F7AF95A1F agree=41A5A1D0932B47EB comember=F9F775B708895AFF sizeLandscape=38E462AD00AADD6D giantLandscape=92DF93E2A0F290FB",
        ["random16/T0.1/none"] = // m̄ = 1.000
            "draws=400 sums=40F9000000000000,0000000000000000,C0C6A5FFFFFFFFD9,411484A028F5C265,4079000000000000,4079000000000000 rng=06AEFC41C62BF090,CF8152D54CC021F3,B54418487A37C050,5264985237DE9E82 hist=4D4834430F89B728 spins=89AE7806201B1EEA bond=- agree=- comember=- sizeLandscape=- giantLandscape=-",
        ["random16/T0.1/full"] = // m̄ = 1.000
            "draws=400 sums=40F9000000000000,0000000000000000,C0C6A5FFFFFFFFD9,411484A028F5C265,4079000000000000,4079000000000000 rng=06AEFC41C62BF090,CF8152D54CC021F3,B54418487A37C050,5264985237DE9E82 hist=4D4834430F89B728 spins=89AE7806201B1EEA bond=365517A39F7B07FF agree=641CD1E1952E3308 comember=641CD1E1952E3308 sizeLandscape=B47EED90447E7EBB giantLandscape=594549C100547139",
        ["random16/T0.6/none"] = // m̄ = 0.474
            "draws=400 sums=40E17B8000000000,40B3E30000000000,C0B544EB851EB84A,40F8CB18A8C154C2,4067AF286BCA1AE8,405E843038BB368F rng=D9EAD947DCD7B197,B65A11DB9D029673,CA4E1B882D3C0AD6,C297BE766B06A4BF hist=2F177618377E03D4 spins=EDDCF51D25D2091D bond=- agree=- comember=- sizeLandscape=- giantLandscape=-",
        ["random16/T0.6/full"] = // m̄ = 0.474
            "draws=400 sums=40E17B8000000000,40B3E30000000000,C0B544EB851EB84A,40F8CB18A8C154C2,4067AF286BCA1AE8,405E843038BB368F rng=D9EAD947DCD7B197,B65A11DB9D029673,CA4E1B882D3C0AD6,C297BE766B06A4BF hist=2F177618377E03D4 spins=EDDCF51D25D2091D bond=4B52961E098B5C8B agree=7A425AEDCD6922BB comember=110DAE7CE732EF31 sizeLandscape=7707CFF9E8F468A1 giantLandscape=719203161F185AC6",
        ["random16/T2/none"] = // m̄ = 0.151
            "draws=400 sums=40BC440000000000,40B72B0000000000,C08EBAF5C28F5C2C,40A9D76000000008,404E3286BCA1AF04,4024987722BF77C7 rng=C96358A9EAD39F4E,96C882A2B351BE94,0F8CA3DAAB8DE0EC,A4AA0BC6D146431B hist=B9807E00508ADFC4 spins=56A70631AABBEAB5 bond=- agree=- comember=- sizeLandscape=- giantLandscape=-",
        ["random16/T2/full"] = // m̄ = 0.151
            "draws=400 sums=40BC440000000000,40B72B0000000000,C08EBAF5C28F5C2C,40A9D76000000008,404E3286BCA1AF04,4024987722BF77C7 rng=C96358A9EAD39F4E,96C882A2B351BE94,0F8CA3DAAB8DE0EC,A4AA0BC6D146431B hist=B9807E00508ADFC4 spins=56A70631AABBEAB5 bond=658AACA49D1DB207 agree=5DD36036A5F3ED72 comember=31683C4275BD6C0E sizeLandscape=C6967091B08294E1 giantLandscape=4AE937F9C9212E63",
    };

    // ── The cases ────────────────────────────────────────────────────────────

    internal sealed record ProbeCase(string Name, CsrGraph Graph, int Q, double Temperature, AccumulationSpec Accumulation);

    internal static readonly RunBudget Budget = new(BurnIn: 50, Cycles: 400);
    internal const int Seed = 20261029;

    /// <summary>Both landscapes beside the three per-edge counts: every optional array the v8 accumulator carries.</summary>
    internal static readonly AccumulationSpec Full =
        AccumulationSpec.Currencies with { ClusterSizeLandscape = true, OrderLandscape = true };

    /// <summary>
    /// A triangle (q = 4) and a 16-node random graph (q = 20, SPC's default), each at a cold, a near-transition and a
    /// hot temperature, under <see cref="AccumulationSpec.None"/> and <see cref="Full"/>.
    /// </summary>
    internal static readonly IReadOnlyList<ProbeCase> Cases = BuildCases();

    private static List<ProbeCase> BuildCases()
    {
        var cases = new List<ProbeCase>();
        foreach ((string graphName, CsrGraph graph, int q, double[] temperatures) in new[]
                 {
                     ("triangle", Triangle(), 4, new[] { 0.3, 1.0, 3.0 }),
                     ("random16", Random16(), 20, new[] { 0.1, 0.6, 2.0 }),
                 })
        {
            foreach (double t in temperatures)
            {
                string tName = t.ToString("0.###", CultureInfo.InvariantCulture);
                cases.Add(new ProbeCase($"{graphName}/T{tName}/none", graph, q, t, AccumulationSpec.None));
                cases.Add(new ProbeCase($"{graphName}/T{tName}/full", graph, q, t, Full));
            }
        }
        return cases;
    }

    internal static Accumulator Run(ProbeCase probe, AccumulationSpec? accumulation = null)
        => SwRunner.Run(new SwRunSpec
        {
            Graph        = probe.Graph,
            Temperature  = probe.Temperature,
            Q            = probe.Q,
            Accumulation = accumulation ?? probe.Accumulation,
            Seed         = Seed,
            Budget       = Budget,
        }).Accumulator;

    internal static CsrGraph Triangle()
        => CsrGraph.FromEdges([new(0, 1, 1.0), new(1, 2, 1.0), new(0, 2, 1.0)], 3);

    /// <summary>A ring of 16 nodes plus 18 random chords, couplings in [0.2, 1.4] (drawn once, kept as a literal).</summary>
    internal static CsrGraph Random16() => CsrGraph.FromEdges(
    [
        new(0, 1, 0.79), new(0, 13, 1.27), new(0, 15, 1.33), new(1, 2, 0.25), new(2, 3, 0.82), new(2, 7, 0.96),
        new(2, 8, 0.96), new(3, 4, 1.2), new(3, 5, 0.71), new(3, 6, 0.59), new(3, 9, 0.73), new(3, 11, 0.97),
        new(4, 5, 1.08), new(4, 8, 0.4), new(5, 6, 0.96), new(5, 9, 0.64), new(6, 7, 0.84), new(6, 12, 0.31),
        new(6, 14, 1.39), new(6, 15, 0.51), new(7, 8, 0.77), new(7, 13, 1.13), new(7, 14, 0.72), new(8, 9, 0.97),
        new(8, 10, 1.35), new(9, 10, 0.81), new(9, 13, 0.4), new(10, 11, 1.02), new(11, 12, 0.72), new(11, 14, 1.22),
        new(12, 13, 0.78), new(12, 15, 0.61), new(13, 14, 0.38), new(14, 15, 1.4),
    ], 16);

    // ── The fingerprint ──────────────────────────────────────────────────────

    /// <summary>
    /// Every output the v8 accumulator carries, but its timestamp: the draw count, the six running sums and four RNG
    /// words as hex bit patterns, and the first 16 hex digits of each array's SHA-256 (<c>-</c> for an absent array).
    /// </summary>
    internal static string Fingerprint(Accumulator acc)
        => string.Join(' ',
            $"draws={acc.DrawCount}",
            "sums=" + string.Join(',', new[]
            {
                acc.RunningSumSqClusterSizes, acc.RunningSumSqClusterSizesExcl, acc.RunningSumEnergy,
                acc.RunningSumEnergySq, acc.RunningSumMag, acc.RunningSumMagSq,
            }.Select(d => BitConverter.DoubleToInt64Bits(d).ToString("X16"))),
            "rng=" + string.Join(',', new[] { acc.RngState0, acc.RngState1, acc.RngState2, acc.RngState3 }
                .Select(w => w.ToString("X16"))),
            $"hist={Hash(acc.ClusterSizeHistogram)}",
            $"spins={Hash(acc.Spins)}",
            $"bond={Hash(acc.BondFormedCount)}",
            $"agree={Hash(acc.SpinAgreementCount)}",
            $"comember={Hash(acc.CoMembershipCount)}",
            $"sizeLandscape={Hash(acc.SumClusterSizePerNode)}",
            $"giantLandscape={Hash(acc.SumInGiantClusterPerNode)}");

    private static string Hash<T>(T[]? values) where T : struct
        => values is null
            ? "-"
            : Convert.ToHexString(SHA256.HashData(MemoryMarshal.AsBytes(values.AsSpan())))[..16];
}
