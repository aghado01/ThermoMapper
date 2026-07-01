using System;
using Maths.Geometry;
using Xunit;

namespace Maths.Geometry.Tests;

public sealed class ScaledManifoldTests
{
    [Fact]
    public void RescalesDistanceAndNorm_LeavesLogUnchanged()
    {
        var inner = new EuclideanVectorManifold(3);
        var scaled = new ScaledManifold<EuclideanVectorManifold>(inner, 4.0); // √4 = 2

        double[] p = { 1.0, 2.0, 3.0 };
        double[] q = { 4.0, 6.0, 3.0 };
        double[] v = { 1.0, 1.0, 1.0 };

        Assert.Equal(2.0 * inner.Distance(p, q), scaled.Distance(p, q), 12);
        Assert.Equal(2.0 * inner.Norm(p, v), scaled.Norm(p, v), 12);

        var logInner = new double[3];
        var logScaled = new double[3];
        inner.LogMap(p, q, logInner);
        scaled.LogMap(p, q, logScaled);
        for (int i = 0; i < 3; i++)
            Assert.Equal(logInner[i], logScaled[i], 12);   // metric scaling leaves Log unchanged
    }
}
