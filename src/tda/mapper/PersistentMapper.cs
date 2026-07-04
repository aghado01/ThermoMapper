#nullable enable
using System;
using System.Collections.Generic;
using Maths.Topology;
using TDA.Ph.Nerves;

using TDA.Ph;
namespace TDA.Mapper;

/// <summary>
/// Persistent Mapper over a scalar parameter sweep (e.g., SPC temperature T or
/// epsilon distance scale).
/// <para>
/// Runs <c>Mapper.Build</c> once per parameter value, builds a
/// <see cref="NerveFiltration"/> from the resulting nerves, and computes a
/// persistence barcode by tracking nerve topology across frames. T is the sweep
/// axis — matching SPC's temperature sweep so nerve topology changes correspond
/// directly to SPC phase transitions on the same axis.
/// </para>
/// </summary>
public static class PersistentMapper
{
    /// <summary>
    /// Build a <see cref="NerveFiltration"/> by running Mapper once per parameter value.
    /// <paramref name="buildFrame"/> is called with each value in <paramref name="parameters"/>
    /// (in the order provided) and must return a <see cref="MapperResult"/> for that value.
    /// Parameters must be provided in non-decreasing order.
    /// </summary>
    public static NerveFiltration BuildFiltration(
        IEnumerable<double> parameters,
        Func<double, MapperResult> buildFrame,
        string parameterLabel = "T")
    {
        ArgumentNullException.ThrowIfNull(parameters);
        ArgumentNullException.ThrowIfNull(buildFrame);

        var frames = new List<NerveFiltrationFrame>();
        int idx = 0;
        foreach (double t in parameters)
        {
            MapperResult result = buildFrame(t);
            int nodeCount = result.Nodes.Count;
            var memberIndices = new int[nodeCount][];
            for (int i = 0; i < nodeCount; i++)
                memberIndices[i] = result.Nodes[i].MemberIndices;
            frames.Add(new NerveFiltrationFrame(t, result.Nerve, memberIndices, idx++));
        }
        return new NerveFiltration(frames, parameterLabel);
    }

    /// <summary>
    /// Build the filtration and immediately compute the H0 persistence barcode.
    /// Convenience overload combining <see cref="BuildFiltration"/> and
    /// <see cref="PersistenceBarcode.ComputeH0"/>.
    /// </summary>
    public static Barcode SweepH0(
        IEnumerable<double> parameters,
        Func<double, MapperResult> buildFrame,
        string parameterLabel = "T")
    {
        NerveFiltration filtration = BuildFiltration(parameters, buildFrame, parameterLabel);
        return PersistenceBarcode.ComputeH0(filtration);
    }
}
