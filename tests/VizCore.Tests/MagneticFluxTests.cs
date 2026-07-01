using System;
using Graphs.Spectral;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// M2/M3 for the magnetic Laplacian at the flux (Aharonov–Bohm holonomy) level —
/// the graph-intrinsic slice that does not yet need the persistence filtration.
/// <para>
/// <b>M2 (reach grader, §4):</b> flux ∝ the backbone span a return bridges, so it
/// separates equal-weight returns that a span-blind measure cannot, and it vanishes
/// at <c>q=0</c>. <b>M3 (return vs revisit, §3):</b> a coherently directed loop
/// encloses flux; an out-and-back excursion cancels to zero.
/// </para>
/// The persistence-coupled separation (the magnetic <em>persistent</em> Laplacian,
/// brief §7) is deferred research and lives downstream of the TDA.Ph filtration.
/// </summary>
public sealed class MagneticFluxTests
{
    [Fact]
    public void Flux_GradesReturnsByBackboneReach()
    {
        const int L = 12;
        const double q = 0.1;
        (int, int)[] backbone = DirectedPath(L); // 0→1→…→11

        // Two returns, identical chord weight (equal persistence proxy), different span:
        // chord 0—2 bridges reach 2; chord 0—9 bridges reach 9.
        var shortReturn = MagneticLaplacianOperator.FromBackboneAndChords(L, backbone, new[] { (0, 2) }, q);
        var longReturn = MagneticLaplacianOperator.FromBackboneAndChords(L, backbone, new[] { (0, 9) }, q);

        double fluxShort = shortReturn.EnclosedFlux(PathWalk(2)); // forward 0→2, chord back to 0
        double fluxLong = longReturn.EnclosedFlux(PathWalk(9));

        // Flux reads the bridged reach directly: Φ = q · span.
        Assert.InRange(Math.Abs(fluxShort - q * 2), 0.0, 1e-12);
        Assert.InRange(Math.Abs(fluxLong - q * 9), 0.0, 1e-12);

        // The flux axis separates equal-weight returns by reach; a span-blind measure can't.
        Assert.True(Math.Abs(fluxLong - fluxShort) > 1e-6,
            $"flux failed to grade reach: short={fluxShort}, long={fluxLong}");

        // The undirected field (q=0) is blind to reach — flux vanishes for both.
        var undirected = MagneticLaplacianOperator.FromBackboneAndChords(L, backbone, new[] { (0, 9) }, 0.0);
        Assert.InRange(Math.Abs(undirected.EnclosedFlux(PathWalk(9))), 0.0, 1e-12);
    }

    [Fact]
    public void Flux_DistinguishesDirectedReturnFromRevisit()
    {
        const double q = 0.1;

        // A coherently directed loop 0→1→2→0 encloses flux 3q (≠ 0).
        var loop = MagneticLaplacianOperator.FromDirectedEdges(3, new[] { (0, 1), (1, 2), (2, 0) }, q);
        double loopFlux = loop.EnclosedFlux(new[] { 0, 1, 2 });
        Assert.InRange(Math.Abs(loopFlux - 3.0 * q), 0.0, 1e-12);
        Assert.True(Math.Abs(loopFlux) > 1e-6, "directed loop should enclose flux");

        // An out-and-back excursion 0→1→0 cancels: the forward and reverse phases negate.
        var excursion = MagneticLaplacianOperator.FromDirectedEdges(2, new[] { (0, 1) }, q);
        double revisitFlux = excursion.EnclosedFlux(new[] { 0, 1 });
        Assert.InRange(Math.Abs(revisitFlux), 0.0, 1e-12);
    }

    private static (int, int)[] DirectedPath(int length)
    {
        var edges = new (int, int)[length - 1];
        for (int i = 0; i < length - 1; i++) edges[i] = (i, i + 1);
        return edges;
    }

    // Closed walk for a return spanning reach Δ: forward 0,1,…,Δ; EnclosedFlux wraps
    // Δ→0 through the (Θ=0) chord, so the flux is purely the backbone span.
    private static int[] PathWalk(int span)
    {
        var walk = new int[span + 1];
        for (int i = 0; i <= span; i++) walk[i] = i;
        return walk;
    }
}
