using System;

namespace TDA.Mapper.Diagnostics;

public readonly record struct CoverageReport(
    int TotalOriginalPoints,
    int PointsCoveredAtLeastOnce,
    int PointsCoveredMultipleTimes,
    int PointsUncovered,
    double CoverageFraction,
    double MeanOverlapMultiplicity,
    int EmptyBinCount);

public static class Coverage
{
    public static CoverageReport From(MapperResult result, int totalOriginalPoints)
    {
        if (totalOriginalPoints <= 0)
        {
            return new CoverageReport(
                TotalOriginalPoints: 0,
                PointsCoveredAtLeastOnce: 0,
                PointsCoveredMultipleTimes: 0,
                PointsUncovered: 0,
                CoverageFraction: 0.0,
                MeanOverlapMultiplicity: 0.0,
                EmptyBinCount: result.EmptyBinCount);
        }

        var coverageCounts = new int[totalOriginalPoints];
        foreach (var node in result.Nodes)
        {
            foreach (int memberIndex in node.MemberIndices)
            {
                if ((uint)memberIndex >= (uint)totalOriginalPoints)
                    throw new ArgumentOutOfRangeException(nameof(totalOriginalPoints),
                        $"Member index {memberIndex} is outside the original point range [0, {totalOriginalPoints}).");

                coverageCounts[memberIndex]++;
            }
        }

        int coveredAtLeastOnce = 0;
        int coveredMultipleTimes = 0;
        long overlapTotal = 0;

        for (int i = 0; i < coverageCounts.Length; i++)
        {
            int count = coverageCounts[i];
            if (count <= 0)
                continue;

            coveredAtLeastOnce++;
            overlapTotal += count;

            if (count > 1)
                coveredMultipleTimes++;
        }

        return new CoverageReport(
            TotalOriginalPoints: totalOriginalPoints,
            PointsCoveredAtLeastOnce: coveredAtLeastOnce,
            PointsCoveredMultipleTimes: coveredMultipleTimes,
            PointsUncovered: totalOriginalPoints - coveredAtLeastOnce,
            CoverageFraction: coveredAtLeastOnce / (double)totalOriginalPoints,
            MeanOverlapMultiplicity: coveredAtLeastOnce > 0 ? overlapTotal / (double)coveredAtLeastOnce : 0.0,
            EmptyBinCount: result.EmptyBinCount);
    }
}
