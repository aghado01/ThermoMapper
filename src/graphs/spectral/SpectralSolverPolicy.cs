// src/graphs/spectral/SpectralSolverPolicy.cs
#nullable enable

namespace Graphs.Spectral;

/// <summary>
/// The single decision point behind <see cref="SolverKind.Auto"/>: maps a
/// bottom-K spectral request <c>(n, k)</c> to a concrete <see cref="SolverKind"/>.
/// Centralising it here means tuning the policy later touches one function, not
/// every call site.
/// </summary>
/// <remarks>
/// <para>
/// <b>Scope.</b> This policy is for the <em>matrix-free graph path</em>
/// (<see cref="Spectral.ComputeBottomK"/>), where LOBPCG's win is that it never
/// materialises the n×n Laplacian (O(n²) memory). It is deliberately <em>not</em>
/// applied to pre-materialised dense inputs — once the matrix exists the memory
/// advantage is already spent.
/// </para>
/// <para>
/// <b>Threshold tuning.</b> Set from the eigen benchmark report
/// (<c>artifacts/reports/eigen-benchmark-report.md</c>, run 20260602): graph-Laplacian
/// bottom-k showed LOBPCG ≈22× faster than dense at n=500 with near-machine-precision
/// agreement, so the graph crossover is well below 500; the harder random-dense case
/// crossed over at ≈256. Since this policy governs the <em>sparse graph</em> path —
/// where LOBPCG wins earlier than on dense matrices — <see cref="IterativeMinNodes"/>
/// is set to the dense crossover (256) as a conservative proxy. The benchmark was a
/// single run on a Balanced power plan (provisional per the report), so the threshold
/// is deliberately not pushed lower than the measured number; revisit with a
/// multi-run graph sweep below n=500 if a lower crossover is wanted.
/// </para>
/// </remarks>
public static class SpectralSolverPolicy
{
    /// <summary>
    /// At/above this node count <see cref="SolverKind.Auto"/> routes to the
    /// matrix-free iterative solver; below it the dense path is faster and more
    /// robust. Anchored to the random-dense crossover (≈256) from the benchmark as
    /// a conservative proxy for the (earlier) sparse-graph crossover.
    /// </summary>
    public const int IterativeMinNodes = 256;

    /// <summary>
    /// Above <see cref="IterativeMinNodes"/>, only route to the iterative solver
    /// when few modes are wanted; LOBPCG degrades as k approaches n. The benchmark
    /// exercised k=8; values beyond ~8 are an (reasonable) extrapolation.
    /// </summary>
    public const int IterativeMaxK = 32;

    /// <summary>
    /// Resolves an <see cref="SolverKind.Auto"/> request to a concrete solver from
    /// the problem shape. Conservative: returns <see cref="SolverKind.Iterative"/>
    /// only for large, low-rank requests; otherwise <see cref="SolverKind.Dense"/>.
    /// </summary>
    public static SolverKind Resolve(int n, int k)
    {
        if (n >= IterativeMinNodes && k > 0 && k <= IterativeMaxK)
            return SolverKind.Iterative;

        return SolverKind.Dense;
    }
}
