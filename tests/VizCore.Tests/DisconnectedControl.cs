using System;

namespace VizCore.Tests;

internal static class DisconnectedControl
{
    internal sealed record Fixture(double[][] Points);

    public static Fixture Generate(int pointsPerComponent = 48, double separation = 14.0)
    {
        var points = new double[pointsPerComponent * 2][];

        for (int i = 0; i < pointsPerComponent; i++)
        {
            double x = 0.25 * i;
            points[i] = new[]
            {
                x,
                0.0,
                0.0
            };
        }

        for (int i = 0; i < pointsPerComponent; i++)
        {
            double x = separation + 0.2 * i;
            points[pointsPerComponent + i] = new[]
            {
                x,
                1.5,
                0.0
            };
        }

        return new Fixture(points);
    }
}
