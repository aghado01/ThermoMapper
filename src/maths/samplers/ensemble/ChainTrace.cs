using System;
using System.Collections.Generic;

namespace Maths.Samplers.Ensemble;

/// <summary>
/// One chain's ordered draws of a fixed set of scalar functionals: what rank-normalized split R̂ and autocorrelation
/// ESS read (<see cref="ChainDiagnostics.Diagnose(IReadOnlyList{ChainTrace}, int)"/>). A trace is not a statistic: the
/// order of its draws carries the information, so it never merges with another chain's; pooling keeps one trace per
/// chain. It is MCMC diagnostic data, distinct from mergeable readout statistics and from restart state, and its ESS
/// is autocorrelation ESS, distinct from a particle population's weight ESS.
/// </summary>
public sealed class ChainTrace
{
    private readonly List<double>[] _draws;

    /// <param name="functionals">How many functionals each draw records.</param>
    /// <param name="capacity">Draws to reserve room for.</param>
    public ChainTrace(int functionals, int capacity = 0)
    {
        if (functionals < 0) throw new ArgumentOutOfRangeException(nameof(functionals));
        _draws = new List<double>[functionals];
        for (int f = 0; f < functionals; f++) _draws[f] = new List<double>(Math.Max(0, capacity));
    }

    /// <summary>The number of functionals each draw records.</summary>
    public int Functionals => _draws.Length;

    /// <summary>The draws recorded so far.</summary>
    public int Count { get; private set; }

    /// <summary>Record one draw: the leading <see cref="Functionals"/> values of <paramref name="values"/>.</summary>
    public void Add(ReadOnlySpan<double> values)
    {
        if (values.Length < _draws.Length)
            throw new ArgumentException($"A draw needs {_draws.Length} values.", nameof(values));
        for (int f = 0; f < _draws.Length; f++) _draws[f].Add(values[f]);
        Count++;
    }

    /// <summary>Record one draw of a single-functional trace.</summary>
    public void Add(double value) => Add([value]);

    /// <summary>Functional <paramref name="functional"/>'s draws, in order.</summary>
    public double[] Draws(int functional) => _draws[functional].ToArray();

    /// <summary>An independent copy.</summary>
    public ChainTrace Copy()
    {
        var c = new ChainTrace(_draws.Length, Count) { Count = Count };
        for (int f = 0; f < _draws.Length; f++) c._draws[f].AddRange(_draws[f]);
        return c;
    }
}

public static partial class ChainDiagnostics
{
    /// <summary>
    /// The rank-normalized diagnostics of functional <paramref name="functional"/> across the chains' traces, which
    /// must all hold the same number of draws (<see cref="Diagnose(double[][])"/>).
    /// </summary>
    public static ConvergenceDiagnostic Diagnose(IReadOnlyList<ChainTrace> traces, int functional)
    {
        ArgumentNullException.ThrowIfNull(traces);
        var draws = new double[traces.Count][];
        for (int c = 0; c < traces.Count; c++)
        {
            ArgumentNullException.ThrowIfNull(traces[c]);
            draws[c] = traces[c].Draws(functional);
        }
        return Diagnose(draws);
    }
}
