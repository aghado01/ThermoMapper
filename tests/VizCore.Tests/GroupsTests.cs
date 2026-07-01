using System;
using Clustering.Primitives;
using Xunit;

namespace VizCore.Tests;

public sealed class GroupsTests
{
    [Fact]
    public void Indexer_Row_ExposeRowMajorMembership()
    {
        // 3 points × 2 groups, row-major: [p0g0,p0g1, p1g0,p1g1, p2g0,p2g1]
        var g = new Groups(new[] { 0.9, 0.1, 0.2, 0.8, 0.5, 0.5 }, pointCount: 3, groupCount: 2);

        Assert.Equal(3, g.PointCount);
        Assert.Equal(2, g.GroupCount);
        Assert.Equal(0.8, g[1, 1], 12);
        Assert.Equal(new[] { 0.5, 0.5 }, g.Row(2).ToArray());
    }

    [Fact]
    public void Argmax_AssignsEachPointToItsMaxGroup_TieToFirst()
    {
        var g = new Groups(new[] { 0.9, 0.1, 0.2, 0.8, 0.5, 0.5 }, pointCount: 3, groupCount: 2);

        Assignment a = g.Argmax();

        Assert.Equal(new[] { 0, 1, 0 }, a.Labels);  // p2's tie resolves to group 0
        Assert.Equal(2, a.Count);
        Assert.Equal(1.0, a.Coverage);              // argmax never abstains
    }

    [Fact]
    public void Argmax_WithNoGroups_LeavesEveryPointUnassigned()
    {
        var g = new Groups(Array.Empty<double>(), pointCount: 2, groupCount: 0);

        Assignment a = g.Argmax();

        Assert.Equal(new[] { Assignment.Unassigned, Assignment.Unassigned }, a.Labels);
        Assert.Equal(0.0, a.Coverage);
    }

    [Fact]
    public void Constructor_RejectsLengthMismatch()
    {
        Assert.Throws<ArgumentException>(() => new Groups(new[] { 0.1, 0.2, 0.3 }, pointCount: 2, groupCount: 2));
    }
}
