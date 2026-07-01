using System.Collections.Generic;

namespace Clustering.Graphical.SPC.Runtime.Scheduling;

public sealed record SchedulePartitionRollup(
    double Temperature,
    int ReplicaCount,
    int ClusterCount,
    int PooledCycleCount,
    IReadOnlyList<double> Purities,
    IReadOnlyList<string> LevelNames);
