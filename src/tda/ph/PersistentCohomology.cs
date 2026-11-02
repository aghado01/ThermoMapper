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

    /// <summary>
    /// The persistence pairs as cell indices, for a reading that needs to know which cells a bar joins. These are the
    /// pairs of the filtration's own order, so they equal a column reducer's, zero-length ones included.
    /// </summary>
    /// <returns>Birth cell, death cell or −1 for an essential class, and degree.</returns>
    internal static List<(int Birth, int Death, int Dimension)> ComputePairs(IFiltration filtration, int maxDimension)
    {
        var h0 = new List<(int Birth, int Death)>();
        (_, HashSet<int> toReduce, HashSet<int> toSkip) = PersistenceClearing.ComputeH0(filtration, maxDimension, 0.0, h0);

        var pairs = new List<(int Birth, int Death, int Dimension)>(filtration.Count);
        foreach ((int birth, int death) in h0)
            pairs.Add((birth, death, 0));

        if (maxDimension <= 0)
            return pairs;

        int maxCellDim = 0;
        for (int i = 0; i < filtration.Count; i++)
            maxCellDim = Math.Max(maxCellDim, filtration.GetDimension(i));

        int dimCap = Math.Min(maxDimension, maxCellDim);
        var reducer = new FiltrationCohomologyReducer(filtration, toReduce, toSkip);
        var deaths = new List<int>();
        var essential = new List<int>();

        for (int dim = 1; dim <= dimCap; dim++)
        {
            deaths.Clear();
            essential.Clear();
            reducer.RunDimension(dim, maxDimension, 0.0, false, null, deaths, essential);

            foreach (int death in deaths)
                pairs.Add((reducer.BirthOf(death), death, dim));
            foreach (int birth in essential)
                pairs.Add((birth, -1, dim));

            if (!reducer.AdvanceDimension(dim, dimCap))
                break;
        }

        return pairs;
    }
}
