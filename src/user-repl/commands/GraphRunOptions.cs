using System;
using System.Collections.Generic;
using System.Linq;
using Graphs;

namespace UserRepl.Commands;

/// <summary>
/// The dataset + coupling-graph options every graph-based clustering command
/// shares (<c>spc</c>, <c>cumulative-field</c>): one dataset source, one graph
/// compiler configuration, one output root. Each command's options derive from
/// it and add their method-specific flags; <see cref="GraphRun"/> parses and
/// executes the shared part.
/// </summary>
internal abstract class GraphRunOptions
{
    public bool ShowHelp { get; set; }
    public bool ListGenerators { get; set; }
    public string? GeneratorSchema { get; set; }

    // Dataset
    public string DatasetKind { get; set; } = GraphRun.DatasetKinds.FirstOrDefault() ?? "unknown";
    public string? DatasetFile { get; set; }
    public string? LabelColumn { get; set; }
    public bool HasHeader { get; set; } = true;
    public char Delimiter { get; set; } = ',';
    public Dictionary<string, object?> GeneratorParameters { get; } = new(StringComparer.OrdinalIgnoreCase);

    // Output / reproducibility
    public string? BaseDirectory { get; set; }
    public string? RunName { get; set; }
    public int Seed { get; set; } = 42;

    // Graph topology
    public int K { get; set; } = 10;
    public double Epsilon { get; set; } = 0.0;
    public TopologyKind TopologyKind { get; set; } = TopologyKind.Knn;
    public FilterKind FilterKind { get; set; } = FilterKind.OrRule;
    public string? DistanceMetricSpec { get; set; }
    public bool EnsureConnected { get; set; }
    public bool ApplyLmp { get; set; }

    // Coupling kernel
    public KernelType KernelType { get; set; } = KernelType.Gaussian;
    public double Bandwidth { get; set; } = 0.0;
    public string? MixtureSpec { get; set; }
    public string? MixtureBandwidthSpec { get; set; }

    /// <summary>Null = resolve from metric properties / the MAD default;
    /// <c>MeanEdgeDistance</c> selects the BWD replication kernel
    /// (mean-pair <c>a</c> + 1/K̂).</summary>
    public Graphs.Distance.BandwidthStrategy? BandwidthStrategy { get; set; }
}
