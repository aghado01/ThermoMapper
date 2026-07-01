using System;

namespace TDA.Mapper.Diagnostics;

public static class MapperWarnings
{
    public static string[] From(
        NerveTopologyReport topology,
        MapperNodeReport nodeStats,
        CoverageReport coverage)
    {
        if (topology.NerveNodeCount == 0)
        {
            return new[]
            {
                "MAPPER produced no nodes — data may be degenerate or filter is invalid",
            };
        }

        var warnings = new string[4];
        int count = 0;

        if (topology.ConnectedComponents > 1)
        {
            warnings[count++] =
                $"Nerve has {topology.ConnectedComponents} disconnected components — possible disjoint manifolds, cover gap, or insufficient cover overlap";
        }

        if (nodeStats.NodeCount > 0 && nodeStats.MinSize < 3)
        {
            warnings[count++] =
                $"Some nodes have <3 members (min: {nodeStats.MinSize}) — increase clusterer K or cover overlap";
        }

        if (coverage.EmptyBinCount > 0)
        {
            warnings[count++] =
                $"{coverage.EmptyBinCount} cover bins are empty — filter range exceeds data range";
        }

        if (topology.LoopCount > 0)
        {
            warnings[count++] =
                $"Nerve contains {topology.LoopCount} loops — data has non-tree topology";
        }

        return count == 0 ? Array.Empty<string>() : warnings[..count];
    }
}
