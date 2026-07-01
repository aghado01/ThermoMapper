#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Fast persistent cohomology (Z/2Z) with H0 union-find clearing — Ripser-style barcode path.
/// </summary>
public static class PersistentCohomology
{
    public static Barcode Compute(
        SimplicialFiltration filtration,
        int maxDimension = int.MaxValue,
        double cutoff = 0.0,
        bool representatives = false) =>
        Compute((IFiltration)filtration, maxDimension, cutoff, representatives);

    public static Barcode Compute(
        IFiltration filtration,
        int maxDimension = int.MaxValue,
        double cutoff = 0.0,
        bool representatives = false)
    {
        ArgumentNullException.ThrowIfNull(filtration);

        int n = filtration.Count;
        var bars = new List<Bar>();
        (List<Bar> h0Bars, HashSet<int> toReduce, HashSet<int> toSkip) =
            PersistenceClearing.ComputeH0(filtration, maxDimension, cutoff);
        bars.AddRange(h0Bars);

        if (maxDimension <= 0)
            return new Barcode(bars, filtration.Label);

        int maxSimplexDim = 0;
        for (int i = 0; i < n; i++)
            if (filtration.GetDimension(i) > maxSimplexDim)
                maxSimplexDim = filtration.GetDimension(i);

        int dimCap = Math.Min(maxDimension, maxSimplexDim);
        var reducer = new FiltrationCohomologyReducer(filtration, toReduce, toSkip);

        for (int dim = 1; dim <= dimCap; dim++)
        {
            reducer.RunDimension(dim, maxDimension, cutoff, representatives, bars, null, null);
            if (!reducer.AdvanceDimension(dim, dimCap))
                break;
        }

        return new Barcode(bars, filtration.Label);
    }
}
