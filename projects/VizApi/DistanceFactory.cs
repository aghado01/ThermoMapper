using System;
using Graphs.Distance;
using Graphs.Distance.Euclidean;
using Graphs.Distance.Geodesic;
using Maths.Distance;
using Maths.Distance.Geodesic;
using Viz;
using VizApi;

internal static class DistanceFactory
{
    public static Func<int, int, double> Create(double[] features, int d, string metricSpec)
    {
        ArgumentNullException.ThrowIfNull(features);
        if (d < 1) throw new ArgumentOutOfRangeException(nameof(d));

        if (TryCreateIntegratedMetric(features, d, metricSpec, out var integrated))
            return integrated;

        if (!Enum.TryParse(metricSpec, true, out VizApiExtraMetric metric))
            throw new NotSupportedException($"Metric '{metricSpec}' is not a registry metric and has no VizApi-only implementation.");

        return metric switch
        {
            VizApiExtraMetric.FisherRaoSimplex => (i, j) => FisherRaoSimplex.Distance(RowOf(features, d, i), RowOf(features, d, j)),
            VizApiExtraMetric.FisherRaoHalfPlane => (i, j) => FisherRaoHalfPlane.Distance(RowOf(features, d, i), RowOf(features, d, j)),
            VizApiExtraMetric.Wasserstein => (i, j) => Wasserstein1.Distance(RowOf(features, d, i), RowOf(features, d, j)),
            VizApiExtraMetric.Jaccard => (i, j) => Jaccard.Distance(RowOf(features, d, i), RowOf(features, d, j)),
            _ => throw new NotSupportedException($"Metric {metric} is not yet wired into graph binding.")
        };
    }

    private static bool TryCreateIntegratedMetric(double[] features, int d, string metricSpec, out Func<int, int, double> distance)
    {
        try
        {
            Graphs.Distance.IDistanceMetric metric = MetricRegistry.Create(metricSpec);
            distance = (i, j) => metric.Distance(features.AsSpan(i * d, d), features.AsSpan(j * d, d));
            return true;
        }
        catch (ArgumentException)
        {
            distance = null!;
            return false;
        }
    }

    private static double[] RowOf(double[] features, int d, int row)
    {
        var result = new double[d];
        Array.Copy(features, row * d, result, 0, d);
        return result;
    }
}
