using Graphs.Distance;

namespace UserRepl.Commands;

/// <summary>
/// Thin CLI-boundary alias for <see cref="MetricRegistry.Create(string)"/>. The
/// metric vocabulary and spec parsing live once in <see cref="MetricRegistry"/>
/// (Graphs.Distance); this preserves the existing call sites (SpcCommand,
/// ExtractCommand, GraphHealthCommand) without re-listing the metric families —
/// the A3 drift between this and HDBSCAN's struct-generic dispatch is gone
/// because both now derive from the registry.
/// </summary>
public static class DistanceMetricFactory
{
    public static IDistanceMetric Create(string spec) => MetricRegistry.Create(spec);
}
