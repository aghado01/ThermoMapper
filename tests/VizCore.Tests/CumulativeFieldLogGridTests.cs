using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using Archivory;
using Clustering.Dendrograms;
using Clustering.Graphical.CumulativeField;
using Clustering.Graphical.CumulativeField.Export;
using Graphs;
using Graphs.Coupling;
using Graphs.Observables;
using Graphs.Primitives;
using UserRepl;
using UserRepl.Commands;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// The cumulative-field command and session on a grid in ln T: the <c>lnspace:</c> spec, the thermal dendrogram's axis
/// label, the session running on the log ladder (a Gaussian gap below the double range resolves; the exports gain their
/// ln T columns), and the temperature grid's outputs staying as they were.
/// </summary>
public sealed class CumulativeFieldLogGridTests
{
    private static SpcUserDataset Dataset()
    {
        double[][] points = CumulativeFieldLogTests.LinePoints();
        int[] labels = Enumerable.Range(0, points.Length).Select(i => i / 10).ToArray();
        return new SpcUserDataset(points, labels, 2, null, new Dictionary<string, object?>());
    }

    private static (CumulativeFieldUserRunResult Result, string Directory) Run(
        GraphBuildResult build, TemperatureGrid grid, double? cutTemperature = null, double? cutLogTemperature = null)
    {
        string directory = Path.Combine(Path.GetTempPath(), "cumulative-field-log-grid-" + Guid.NewGuid().ToString("N"));
        var paths = new CumulativeFieldRunPaths(ArtifactScope.Root(directory, "test", "stamp").EnsureDirectory());
        try
        {
            CumulativeFieldUserRunResult result = CumulativeFieldUserSession.Run(
                Dataset(), paths, build.Graph, EdgeWeightKind.Coupling, Field.Mean, SymmetrizationRule.Mutual,
                grid, theta: 0.5, cutTemperature: cutTemperature, cutLogTemperature: cutLogTemperature);
            return (result, directory);
        }
        catch
        {
            Directory.Delete(directory, recursive: true);
            throw;
        }
    }

    private static string[] Header(string csvPath) => File.ReadLines(csvPath).First().Split(',');

    private static double[] Column(string csvPath, string name)
    {
        string[] rows = File.ReadAllLines(csvPath);
        int index = Array.IndexOf(rows[0].Split(','), name);
        Assert.True(index >= 0, $"{name} is not a column of {csvPath}");
        return rows.Skip(1).Select(r => double.Parse(r.Split(',')[index], System.Globalization.CultureInfo.InvariantCulture)).ToArray();
    }

    [Fact]
    public void Lnspace_ParsesIntoALogGrid_AndTheTemperatureParserRefusesIt()
    {
        TemperatureGrid grid = TemperatureGridSpec.ParseGrid("lnspace:-1300,0,14");

        Assert.True(grid.IsLog);
        Assert.Equal(14, grid.Values.Length);
        Assert.Equal(-1300.0, grid.Values[0]);
        Assert.Equal(0.0, grid.Values[^1]);
        Assert.Equal(100.0, grid.Values[1] - grid.Values[0], 9);
        Assert.Equal(0.0, grid.Temperature(0));                       // exp(-1300) underflows; the log is the record
        Assert.Equal(-1300.0, grid.LogTemperature(0));

        TemperatureGrid plain = TemperatureGridSpec.ParseGrid("logspace:0.01,1,5");
        Assert.False(plain.IsLog);
        Assert.Equal(Math.Log(0.01), plain.LogTemperature(0), 12);

        Assert.Single(TemperatureGridSpec.ParseGrid("lnspace:-5,-5,1").Values);
        Assert.Throws<ArgumentException>(() => TemperatureGridSpec.Parse("lnspace:-1300,0,14"));
        Assert.Throws<ArgumentException>(() => TemperatureGridSpec.ParseGrid("lnspace:0,-5,4"));
        Assert.Throws<ArgumentException>(() => TemperatureGridSpec.ParseGrid("lnspace:-5,Infinity,4"));
        Assert.Throws<ArgumentException>(() => TemperatureGridSpec.ParseGrid("lnspace:-5,0"));
    }

    [Fact]
    public void TheThermalDendrogram_CarriesTheAxisItWasGiven()
    {
        // Path 0-1-2. Edge (0,1) holds through the whole grid, edge (1,2) through its first two points.
        CsrGraph graph = CsrGraph.FromEdges(new[] { new Edge(0, 1, 1.0), new Edge(1, 2, 1.0) }, 3);
        double[] grid = { -1300.0, -1000.0, -700.0 };
        var columns = new double[3][];
        for (int t = 0; t < 3; t++)
        {
            columns[t] = new double[graph.Targets.Length];
            columns[t][graph.FindSlot(0, 1)] = 1.0;
            columns[t][graph.FindSlot(1, 2)] = t < 2 ? 1.0 : 0.0;
        }

        Dendrogram log = ThermalDendrogram.FromEdgeCurves(graph, grid, columns, 0.5, ThermalDendrogram.LogTemperatureAxis);
        Dendrogram plain = ThermalDendrogram.FromEdgeCurves(graph, new[] { 0.1, 0.2, 0.3 }, columns, 0.5);

        Assert.Equal(ThermalDendrogram.LogTemperatureAxis, log.CostAxis);
        Assert.Equal(ThermalDendrogram.TemperatureAxis, plain.CostAxis);
        Assert.Equal(new[] { -700.0, -1000.0 }, log.Merges.Select(m => m.Distance));   // hottest decoupling first
    }

    [Fact]
    public void TheSession_OnALogGrid_ResolvesAGaussianGapBelowTheDoubleRange()
    {
        GraphBuildResult build = CumulativeFieldLogTests.TwoLines(new Gaussian());
        double gap = MergeHeights.Compute(build.Graph)[^1].LogCoupling;
        Assert.True(gap < -1000.0);
        double[] values = Enumerable.Range(0, 60).Select(i => gap - 10.0 + i * ((3.0 - (gap - 10.0)) / 59.0)).ToArray();

        (CumulativeFieldUserRunResult result, string directory) = Run(build, new TemperatureGrid(values, IsLog: true));
        try
        {
            // The two lines are two clusters: a stage over a thousand e-folds wide, which no temperature grid reaches.
            Assert.Equal(2, result.ClusterCount);
            Assert.Equal(0, result.Abstained);

            Assert.Equal(new[] { "log_temperature", "temperature", "cluster_count" }, Header(result.SweepCsvPath));
            Assert.Equal(values, Column(result.SweepCsvPath, "log_temperature"));
            double[] counts = Column(result.SweepCsvPath, "cluster_count");
            Assert.Equal(1.0, counts[0]);
            Assert.True(counts[^1] > 2);
            Assert.Equal(0.0, Column(result.SweepCsvPath, "temperature")[0]);

            // The edge snapshot is at the coldest point, and its log_weight holds the bridge the linear weight cannot.
            Assert.Equal(values[0], result.RepresentativeLogTemperature);
            Assert.Contains("log_weight", Header(result.EdgeCurrenciesCsvPath!));
            Assert.InRange(Column(result.EdgeCurrenciesCsvPath!, "log_weight").Min(), gap - 1e-6, gap + 1e-6);

            using JsonDocument summary = JsonDocument.Parse(File.ReadAllText(result.SummaryJsonPath));
            Assert.Equal(ThermalDendrogram.LogTemperatureAxis, summary.RootElement.GetProperty("temperatureAxis").GetString());
            Assert.Equal(values[0], summary.RootElement.GetProperty("representativeLogTemperature").GetDouble());
            Assert.Equal(60, summary.RootElement.GetProperty("temperatureCount").GetInt32());
        }
        finally { Directory.Delete(directory, recursive: true); }
    }

    [Fact]
    public void TheSession_OnATemperatureGrid_WritesTheFilesItAlwaysDid()
    {
        GraphBuildResult build = CumulativeFieldLogTests.TwoLines(new Cauchy());
        double[] values = TemperatureGridSpec.Parse("logspace:0.0001,1,24");

        (CumulativeFieldUserRunResult result, string directory) = Run(build, new TemperatureGrid(values, IsLog: false));
        try
        {
            Assert.Equal(new[] { "temperature", "cluster_count" }, Header(result.SweepCsvPath));
            Assert.DoesNotContain("log_weight", Header(result.EdgeCurrenciesCsvPath!));
            Assert.Equal(values[0], result.RepresentativeTemperature);

            using JsonDocument summary = JsonDocument.Parse(File.ReadAllText(result.SummaryJsonPath));
            Assert.Equal(ThermalDendrogram.TemperatureAxis, summary.RootElement.GetProperty("temperatureAxis").GetString());
        }
        finally { Directory.Delete(directory, recursive: true); }
    }

    [Fact]
    public void ALogGrid_AndTheTemperatureGridItIsTheLogOf_CountTheSameClustersAtEveryPoint()
    {
        GraphBuildResult build = CumulativeFieldLogTests.TwoLines(new Cauchy());
        double[] temperatures = TemperatureGridSpec.Parse("logspace:0.0001,1,24");
        double[] logs = temperatures.Select(t => Math.Log(t)).ToArray();

        (CumulativeFieldUserRunResult plain, string plainDirectory) = Run(build, new TemperatureGrid(temperatures, IsLog: false));
        (CumulativeFieldUserRunResult log, string logDirectory) = Run(build, new TemperatureGrid(logs, IsLog: true));
        try
        {
            Assert.Equal(Column(plain.SweepCsvPath, "cluster_count"), Column(log.SweepCsvPath, "cluster_count"));
        }
        finally
        {
            Directory.Delete(plainDirectory, recursive: true);
            Directory.Delete(logDirectory, recursive: true);
        }
    }

    [Fact]
    public void TheEdgeSnapshot_IsPickedInTheGridsOwnCoordinate()
    {
        GraphBuildResult build = CumulativeFieldLogTests.TwoLines(new Cauchy());
        var log = new TemperatureGrid(new[] { -9.0, -6.0, -3.0, 0.0 }, IsLog: true);
        var plain = new TemperatureGrid(new[] { 0.001, 0.01, 0.1, 1.0 }, IsLog: false);

        var directories = new List<string>();
        try
        {
            (CumulativeFieldUserRunResult a, string da) = Run(build, log, cutLogTemperature: -5.0);
            (CumulativeFieldUserRunResult b, string db) = Run(build, log, cutTemperature: Math.Exp(-2.9));
            (CumulativeFieldUserRunResult c, string dc) = Run(build, plain, cutLogTemperature: Math.Log(0.012));
            directories.AddRange(new[] { da, db, dc });

            Assert.Equal(-6.0, a.RepresentativeLogTemperature);
            Assert.Equal(-3.0, b.RepresentativeLogTemperature);
            Assert.Equal(0.01, c.RepresentativeTemperature);
        }
        finally { foreach (string d in directories) Directory.Delete(d, recursive: true); }

        Assert.Throws<ArgumentException>(() => Run(build, log, cutTemperature: 0.1, cutLogTemperature: -2.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Run(build, log, cutTemperature: 0.0));
    }

    [Fact]
    public void TheCommand_ResolvesAutoAsATemperatureGrid_AndAutoLnAsTheSameBracketInLogs()
    {
        CsrGraph graph = CumulativeFieldLogTests.TwoLines(new Cauchy()).Graph;
        var options = new CumulativeFieldCommand.Options();

        options.TemperaturesSpec = "auto";
        (TemperatureGrid auto, string autoSpec, string autoBasis) = CumulativeFieldCommand.ResolveTemperatureGrid(options, graph);
        options.TemperaturesSpec = "auto-ln";
        (TemperatureGrid autoLn, string autoLnSpec, string autoLnBasis) = CumulativeFieldCommand.ResolveTemperatureGrid(options, graph);

        Assert.False(auto.IsLog);
        Assert.StartsWith("logspace:", autoSpec);
        Assert.Equal("own-ladder", autoBasis);

        Assert.True(autoLn.IsLog);
        Assert.StartsWith("lnspace:", autoLnSpec);
        Assert.Equal("own-ladder", autoLnBasis);
        Assert.Equal(48, autoLn.Values.Length);
        Assert.Equal(auto.Values.Select(t => Math.Log(t)), autoLn.Values, new RelativeComparer(1e-9));

        // The resolved spec records the grid exactly: re-parsing it reproduces it.
        Assert.Equal(autoLn.Values, TemperatureGridSpec.ParseGrid(autoLnSpec).Values);

        options.TemperaturesSpec = "lnspace:-40,-4,10";
        (TemperatureGrid explicitLog, _, string explicitBasis) = CumulativeFieldCommand.ResolveTemperatureGrid(options, graph);
        Assert.True(explicitLog.IsLog);
        Assert.Equal("explicit", explicitBasis);
    }

    [Fact]
    public void TheCommand_TakesTheCutAsTOrAsLnT_ButNotBoth()
    {
        Assert.Equal(-3.5, CumulativeFieldCommand.ParseArguments(new[] { "--cut-log-temperature", "-3.5" }).CutLogTemperature);
        Assert.Equal(0.1, CumulativeFieldCommand.ParseArguments(new[] { "--cut-temperature", "0.1" }).CutTemperature);
        Assert.Throws<ArgumentException>(() => CumulativeFieldCommand.ParseArguments(
            new[] { "--cut-temperature", "0.1", "--cut-log-temperature", "-2" }));
    }

    private sealed class RelativeComparer : IEqualityComparer<double>
    {
        private readonly double _tolerance;
        public RelativeComparer(double tolerance) => _tolerance = tolerance;
        public bool Equals(double x, double y) => Math.Abs(x - y) <= _tolerance * Math.Max(1.0, Math.Abs(x));
        public int GetHashCode(double value) => 0;
    }
}
