using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace Maths.LinAlg;

/// <summary>The kernel a decomposition ran in.</summary>
public enum DecompositionKernel
{
    /// <summary><see cref="Eigen"/>, the reference cyclic Jacobi.</summary>
    ReferenceJacobi,

    /// <summary><see cref="EigenFast"/>, the production cyclic Jacobi.</summary>
    FastJacobi,

    /// <summary><see cref="ThinSvd"/>.</summary>
    ThinSvd,
}

/// <summary>One decomposition: the kernel it ran in and the matrix's shape.</summary>
public readonly record struct Decomposition(DecompositionKernel Kernel, int Rows, int Columns);

/// <summary>
/// The decompositions run inside a scope, recorded at the executing kernels, so a test can check the work a computation performed:
/// that a readout decomposes a matrix once, or that a solve runs no reference decomposition (spectral-solver brief). A scope belongs
/// to the async flow that opened it, so tests running in parallel do not see each other's work, and nothing is recorded outside one.
/// Scopes nest, each recording what runs inside it.
/// </summary>
public sealed class DecompositionLog : IDisposable
{
    private static readonly AsyncLocal<DecompositionLog?> Current = new();
    private readonly DecompositionLog? _outer;
    private readonly List<Decomposition> _entries = [];
    private bool _closed;

    private DecompositionLog(DecompositionLog? outer) => _outer = outer;

    /// <summary>Opens a scope on the current async flow; disposing it closes the scope.</summary>
    public static DecompositionLog Begin()
    {
        var log = new DecompositionLog(Current.Value);
        Current.Value = log;
        return log;
    }

    /// <summary>The decompositions recorded so far, in the order they started.</summary>
    public IReadOnlyList<Decomposition> Entries
    {
        get
        {
            lock (_entries) return [.. _entries];
        }
    }

    /// <summary>The number of decompositions that ran in <paramref name="kernel"/>.</summary>
    public int Count(DecompositionKernel kernel) => Entries.Count(d => d.Kernel == kernel);

    /// <summary>Closes the scope; the enclosing scope, if any, records from here on.</summary>
    public void Dispose()
    {
        if (_closed) return;
        _closed = true;
        Current.Value = _outer;
    }

    internal static void Record(DecompositionKernel kernel, int rows, int columns)
    {
        for (DecompositionLog? log = Current.Value; log is not null; log = log._outer)
            lock (log._entries) log._entries.Add(new Decomposition(kernel, rows, columns));
    }
}
