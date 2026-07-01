using System;
using Clustering.Primitives;
using Xunit;

namespace VizCore.Tests;

public sealed class AssignmentTests
{
    [Fact]
    public void IsAssigned_AssignedCount_Coverage_RespectUnassignedSentinel()
    {
        var a = new Assignment { Labels = new[] { 0, 1, Assignment.Unassigned, 1, 0 }, Count = 2 };

        Assert.Equal(5, a.PointCount);
        Assert.Equal(4, a.AssignedCount);
        Assert.Equal(0.8, a.Coverage, 12);
        Assert.True(a.IsAssigned(0));
        Assert.False(a.IsAssigned(2));
    }

    [Fact]
    public void Assigned_YieldsOnlyAssignedPairs_SkippingUnassigned()
    {
        var a = new Assignment { Labels = new[] { 0, 1, Assignment.Unassigned, 1, 0 }, Count = 2 };

        Assert.Equal(new[] { (0, 0), (1, 1), (3, 1), (4, 0) }, a.Assigned);
    }

    [Fact]
    public void FromLabels_DerivesCount_IgnoringUnassigned()
    {
        var a = Assignment.FromLabels(new[] { 2, 0, Assignment.Unassigned, 1 });

        Assert.Equal(3, a.Count);
        Assert.InRange(Math.Abs(a.Coverage - 0.75), 0.0, 1e-12);
    }

    [Fact]
    public void EmptyAssignment_HasZeroCoverageAndNoAssignedPairs()
    {
        var a = new Assignment { Labels = Array.Empty<int>(), Count = 0 };

        Assert.Equal(0, a.PointCount);
        Assert.Equal(0.0, a.Coverage);
        Assert.Empty(a.Assigned);
    }
}
