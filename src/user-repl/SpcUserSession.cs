using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text.Json;
using System.Xml.Linq;
using Clustering.Evaluation.External;
using UserRepl.Commands;
using Clustering.Graphical.SPC;
using Clustering.Evaluation.Internal;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Scheduling;
using Clustering.Graphical.SPC.Export;
using Clustering.Graphical.SPC.Partitions;
using Clustering.Graphical.SPC.Partitions.Hierarchical;
using Clustering.Graphical.SPC.Partitions.Strategies;
using Clustering.Graphical.SPC.Profiling.Signals;
using Graphs;
using Graphs.Distance;
using Graphs.Primitives;
using Synthetic;

namespace UserRepl;

public static class SpcUserSession
{
    [RequiresUnreferencedCode(
    "Synthetic generator discovery scans the Synthetic assembly for sealed static classes with a public static Generate method. " +
    "The Synthetic assembly is preserved via TrimmerRootAssembly in the publish csproj.")]
    public static IReadOnlyList<SpcGeneratorInfo> ListAvailableSyntheticGenerators() => SyntheticGeneratorCatalog.Generators;

    [RequiresUnreferencedCode(
    "Synthetic generator discovery scans the Synthetic assembly for sealed static classes with a public static Generate method. " +
    "The Synthetic assembly is preserved via TrimmerRootAssembly in the publish csproj.")]
    public static IReadOnlyList<string> ListAvailableSyntheticGeneratorNames() => SyntheticGeneratorCatalog.Generators.Select(g => g.GeneratorName).ToArray();

    [RequiresUnreferencedCode(
    "Synthetic generator invocation uses descriptors discovered from the Synthetic assembly. " +
    "The Synthetic assembly is preserved via TrimmerRootAssembly in the publish csproj.")]
    public static SpcUserDataset GenerateDataset(string generatorName, IDictionary<string, object?>? parameters = null)
    {
        var dataset = SyntheticGeneratorCatalog.Invoke(generatorName, parameters);
        return SpcUserDataset.FromSyntheticDataset(dataset, parameters == null ? null : new Dictionary<string, object?>(parameters, StringComparer.OrdinalIgnoreCase));
    }

    public static SpcUserDataset FromFeatures(
        double[][] features,
        int[] labels,
        int clusterCount = 0,
        int[][]? labelsByLevel = null,
        IReadOnlyDictionary<string, object?>? metadata = null)
        => SpcUserDataset.FromFeatures(features, labels, clusterCount, labelsByLevel, metadata);

    public static SpcUserDataset FromCsv(
        string path,
        string? labelColumn = null,
        bool hasHeader = true,
        char separator = ',')
        => SpcUserDataset.FromCsv(path, labelColumn, hasHeader, separator);

    public static SpcUserRunResult Run(
        SpcUserDataset dataset,
        GraphCompilerConfig graphConfig,
        IDistanceMetric? metric,
        SpcRunPaths paths,
        IEnumerable<IExternalClusterEvaluator>? externalEvaluators = null,
        IEnumerable<IGraphPartitionEvaluator>? spcEvaluators = null,
        int[]? referenceLabels = null,
        IPartitionStrategy? partitionStrategy = null,
        ISignalAnalyzer? analyzer = null,
        ISweepStrategy? sweepStrategy = null,
        IHierarchicalPartitionStrategy? hierarchicalStrategy = null,
        CsrGraph? prebuiltGraph = null)
    {
        if (dataset is null)
            throw new ArgumentNullException(nameof(dataset));
        if (paths is null)
            throw new ArgumentNullException(nameof(paths));

        // The scope is the single owner of "where"; we just materialize the
        // directories we write into. No run-directory creation here.
        string runDirectory = paths.Scope.EnsureDirectory().Dir;
        var csv = paths.Csv.EnsureDirectory();
        string sweepPath          = csv.File(SpcCsvWriter.SweepFileName);
        string partitionPath      = csv.File(SpcCsvWriter.PartitionFileName);
        string criteriaPath       = csv.File(SpcCsvWriter.CriteriaFileName);
        string sessionPath        = csv.File(SpcCsvWriter.SessionFileName);
        string replicaTracesPath  = csv.File(SpcCsvWriter.ReplicaTracesFileName);
        string equilibriumEdgesPath = csv.File(SpcCsvWriter.EquilibriumEdgesFileName);

        // When the caller has already built the graph (typically to print a
        // pre-run health recommendation via Graphs.Diagnostics.GraphHealth),
        // hand it straight to the CsrGraph overload of SpcClusteringSession.Run.
        // Otherwise the legacy features→graph path inside the session does
        // the build itself.
        SpcSessionResult result = prebuiltGraph is CsrGraph graph
            ? SpcClusteringSession.Run(
                graph,
                partitionStrategy: partitionStrategy,
                analyzer: analyzer,
                sweepStrategy: sweepStrategy,
                hierarchicalStrategy: hierarchicalStrategy,
                externalEvaluators: externalEvaluators,
                spcEvaluators: spcEvaluators,
                referenceLabels: referenceLabels)
            : SpcClusteringSession.Run(
                dataset.Features,
                graphConfig: graphConfig,
                partitionStrategy: partitionStrategy,
                analyzer: analyzer,
                sweepStrategy: sweepStrategy,
                hierarchicalStrategy: hierarchicalStrategy,
                externalEvaluators: externalEvaluators,
                spcEvaluators: spcEvaluators,
                referenceLabels: referenceLabels);

        SpcCsvWriter.WriteSweepProfile(result.Profile, sweepPath);
        SpcCsvWriter.WritePartition(result.Partition, partitionPath, features: dataset.Features, trueLabels: dataset.Labels);
        SpcCsvWriter.WriteCriteria(result.ProfileCriteria, criteriaPath);
        SpcCsvWriter.WriteSessionSummary(result, sessionPath);
        SpcCsvWriter.WriteReplicaTraces(result.SweepRuns, replicaTracesPath);

        // Chosen-T edge currencies — returns null on a trivial graph (G[].Length == 0).
        string? equilibriumEdgesWritten = SpcCsvWriter.WriteEquilibriumEdges(
            result.Graph, result.ChosenAffinities, result.ChosenAlignments, equilibriumEdgesPath,
            coMembership: result.ChosenCoMembership);

        // Persist the full cross-T partition hierarchy as a sidecar when a
        // hierarchical strategy ran. Flat partition stays the primary cut
        // (consumed by evaluators + plotters); hierarchy is supplemental.
        string? partitionHierarchyPath = null;
        if (result.Hierarchy is not null)
        {
            partitionHierarchyPath = paths.Hierarchy;
            UserReplJson.Writer.WriteDocumentToFile(result.Hierarchy, partitionHierarchyPath);
        }

        string summaryPath = paths.Summary;
        var summaryPayload = new SpcRunSummary(
            Dataset: dataset.Metadata,
            Graph: graphConfig,
            Analyzer: (analyzer ?? new ChiPeakSignalAnalyzer()).GetType().Name,
            PartitionStrategy: (partitionStrategy ?? new ThresholdSpinAgreement { Theta = 0.5 }).GetType().Name,
            ReferenceLabels: referenceLabels,
            Partition: new SpcSummaryPartition(
                ClusterCount: result.Partition.Count,
                EvaluatorScores: result.EvaluatorScores),
            Run: new SpcSummaryRun(
                RunDirectory: runDirectory,
                SweepCsv: sweepPath,
                PartitionCsv: partitionPath,
                CriteriaCsv: criteriaPath,
                SessionCsv: sessionPath,
                ReplicaTracesCsv: replicaTracesPath,
                FinalEdgesCsv: equilibriumEdgesWritten));

        File.WriteAllText(summaryPath, JsonSerializer.Serialize(summaryPayload, UserReplJson.IndentedOptions));

        return new SpcUserRunResult(
            SessionResult: result,
            Dataset: dataset,
            RunDirectory: runDirectory,
            SweepCsvPath: sweepPath,
            PartitionCsvPath: partitionPath,
            CriteriaCsvPath: criteriaPath,
            SessionCsvPath: sessionPath,
            SummaryJsonPath: summaryPath,
            ReplicaTracesCsvPath: replicaTracesPath,
            EquilibriumEdgesCsvPath: equilibriumEdgesWritten,
            PartitionHierarchyJsonPath: partitionHierarchyPath);
    }
}

internal static class SyntheticGeneratorCatalog
{
    private static readonly Lazy<IReadOnlyList<SyntheticGeneratorDescriptor>> s_generators = new(Discover);

    public static IReadOnlyList<SpcGeneratorInfo> Generators =>
        s_generators.Value.Select(descriptor => new SpcGeneratorInfo(
            descriptor.GeneratorName,
            descriptor.TypeName,
            descriptor.DocInfo?.Summary,
            descriptor.GenerateMethod.GetParameters()
                .Select(parameter => new SpcGeneratorParameter(
                    parameter.Name ?? string.Empty,
                    GetDisplayType(parameter.ParameterType),
                    parameter.HasDefaultValue,
                    parameter.HasDefaultValue ? parameter.DefaultValue : null,
                    descriptor.DocInfo?.ParameterDescriptions.TryGetValue(parameter.Name ?? string.Empty, out var description) == true
                        ? description
                        : null))
                .ToArray()))
        .ToArray();

    [RequiresUnreferencedCode(
        "Synthetic generator invocation uses descriptors discovered from the Synthetic assembly. " +
        "The Synthetic assembly is preserved via TrimmerRootAssembly in the publish csproj.")]
    public static SyntheticDataset Invoke(string generatorName, IDictionary<string, object?>? parameters)
    {
        if (generatorName is null)
            throw new ArgumentNullException(nameof(generatorName));

        var descriptor = s_generators.Value.FirstOrDefault(descriptor =>
            string.Equals(descriptor.GeneratorName, generatorName, StringComparison.OrdinalIgnoreCase) ||
            string.Equals(descriptor.TypeName, generatorName, StringComparison.OrdinalIgnoreCase));

        if (descriptor is null)
            throw new ArgumentException($"Unknown synthetic generator '{generatorName}'.", nameof(generatorName));

        var parameterInfos = descriptor.GenerateMethod.GetParameters();
        var args = new object?[parameterInfos.Length];
        var unknownKeys = new HashSet<string>(StringComparer.OrdinalIgnoreCase);

        if (parameters != null)
        {
            foreach (var key in parameters.Keys)
                unknownKeys.Add(key);
        }

        for (int i = 0; i < parameterInfos.Length; i++)
        {
            var parameter = parameterInfos[i];
            if (parameters != null && parameters.TryGetValue(parameter.Name ?? string.Empty, out var value))
            {
                args[i] = ConvertParameterValue(value, parameter.ParameterType);
                unknownKeys.Remove(parameter.Name ?? string.Empty);
            }
            else if (parameter.HasDefaultValue)
            {
                args[i] = parameter.DefaultValue;
            }
            else
            {
                throw new ArgumentException($"Missing required parameter '{parameter.Name}' for generator '{generatorName}'.");
            }
        }

        if (unknownKeys.Count > 0)
            throw new ArgumentException($"Unknown parameter(s) for generator '{generatorName}': {string.Join(", ", unknownKeys)}.");

        return (SyntheticDataset)descriptor.GenerateMethod.Invoke(null, args)!;
    }

    private static string GetDisplayType(Type type)
    {
        if (type.IsArray)
            return $"{GetDisplayType(type.GetElementType()!)}[]";
        return type.Name;
    }

    private static object? ConvertParameterValue(object? value, Type targetType)
    {
        if (value is JsonElement jsonElement)
            return ConvertJsonElement(jsonElement, targetType);

        if (value is null)
        {
            if (!targetType.IsValueType || Nullable.GetUnderlyingType(targetType) != null)
                return null;
            throw new ArgumentException($"Cannot assign null to parameter type {targetType.Name}.");
        }

        var effectiveType = Nullable.GetUnderlyingType(targetType) ?? targetType;

        if (effectiveType.IsAssignableFrom(value.GetType()))
            return value;

        if (effectiveType.IsArray)
        {
            var elementType = effectiveType.GetElementType()!;
            string[] values;

            if (value is string valueString)
            {
                values = valueString.Split(',', StringSplitOptions.RemoveEmptyEntries)
                    .Select(entry => entry.Trim())
                    .ToArray();
            }
            else if (value is Array inputArray)
            {
                values = new string[inputArray.Length];
                for (int i = 0; i < inputArray.Length; i++)
                    values[i] = inputArray.GetValue(i)?.ToString() ?? string.Empty;
            }
            else
            {
                throw new ArgumentException($"Cannot convert parameter value of type {value.GetType().Name} to array type {effectiveType.Name}.");
            }

            var resultArray = Array.CreateInstance(elementType, values.Length);
            for (int i = 0; i < values.Length; i++)
                resultArray.SetValue(ConvertParameterValue(values[i], elementType), i);
            return resultArray;
        }

        try
        {
            return Convert.ChangeType(value, effectiveType, CultureInfo.InvariantCulture);
        }
        catch (Exception ex)
        {
            throw new ArgumentException($"Cannot convert parameter value '{value}' to type {effectiveType.Name}.", ex);
        }
    }

    private static object? ConvertJsonElement(JsonElement value, Type targetType)
    {
        switch (value.ValueKind)
        {
            case JsonValueKind.Undefined:
            case JsonValueKind.Null:
                return ConvertParameterValue(null, targetType);

            case JsonValueKind.Array:
            {
                var effectiveType = Nullable.GetUnderlyingType(targetType) ?? targetType;
                if (!effectiveType.IsArray)
                    throw new ArgumentException(
                        $"Cannot convert JSON array parameter to non-array type {effectiveType.Name}.");

                var elementType = effectiveType.GetElementType()!;
                JsonElement[] items = value.EnumerateArray().ToArray();
                Array result = Array.CreateInstance(elementType, items.Length);
                for (int i = 0; i < items.Length; i++)
                    result.SetValue(ConvertParameterValue(items[i], elementType), i);
                return result;
            }

            case JsonValueKind.String:
                return ConvertParameterValue(value.GetString(), targetType);

            case JsonValueKind.Number:
            case JsonValueKind.True:
            case JsonValueKind.False:
                return ConvertParameterValue(value.GetRawText(), targetType);

            default:
                throw new ArgumentException(
                    $"Cannot convert JSON {value.ValueKind} parameter to type {(Nullable.GetUnderlyingType(targetType) ?? targetType).Name}.");
        }
    }

    [RequiresUnreferencedCode(
        "Synthetic generator discovery scans the Synthetic assembly for sealed static classes with a public static Generate method. " +
        "The Synthetic assembly is preserved via TrimmerRootAssembly in the publish csproj.")]
    private static IReadOnlyList<SyntheticGeneratorDescriptor> Discover()
    {
        var syntheticAssembly = typeof(SyntheticDataset).Assembly;
        var descriptors = new List<SyntheticGeneratorDescriptor>();

        foreach (var type in syntheticAssembly.GetTypes())
        {
            if (!type.IsClass || !type.IsAbstract || !type.IsSealed)
                continue;

            if (type.GetMethod("Generate", BindingFlags.Public | BindingFlags.Static) is not MethodInfo method)
                continue;

            if (!typeof(SyntheticDataset).IsAssignableFrom(method.ReturnType))
                continue;

            if (string.IsNullOrWhiteSpace(type.Name))
                continue;

            descriptors.Add(new SyntheticGeneratorDescriptor(
                GeneratorName: type.Name,
                TypeName: type.FullName ?? type.Name,
                GenerateMethod: method,
                DocInfo: TryGetMethodDoc(method)));
        }

        return descriptors.OrderBy(d => d.GeneratorName, StringComparer.OrdinalIgnoreCase).ToArray();
    }

    private static MethodDocInfo? TryGetMethodDoc(MethodInfo method)
    {
        string? sourceRoot = FindSyntheticSourceRoot();
        if (sourceRoot is null)
            return null;

        string fileName = method.DeclaringType?.Name + ".cs";
        if (fileName is null)
            return null;

        string[] files = Directory.GetFiles(sourceRoot, fileName, SearchOption.AllDirectories);
        if (files.Length == 0)
            return null;

        string filePath = files[0];
        string[] lines = File.ReadAllLines(filePath);

        for (int i = 0; i < lines.Length; i++)
        {
            if (lines[i].Contains("public static SyntheticDataset Generate(") ||
                lines[i].Contains("public static SyntheticDataset Generate ("))
            {
                var commentLines = new List<string>();
                for (int j = i - 1; j >= 0; j--)
                {
                    string trimmed = lines[j].TrimStart();
                    if (trimmed.StartsWith("///"))
                    {
                        commentLines.Insert(0, trimmed.Substring(3).Trim());
                        continue;
                    }

                    if (string.IsNullOrWhiteSpace(trimmed))
                        continue;

                    break;
                }

                if (commentLines.Count == 0)
                    return null;

                string xml = "<root>" + string.Join("\n", commentLines) + "</root>";
                try
                {
                    var doc = XDocument.Parse(xml);
                    string? summary = doc.Root?.Element("summary")?.Value.Trim();
                    var paramDescriptions = doc.Root?.Elements("param")
                        .Where(e => e.Attribute("name") != null)
                        .ToDictionary(
                            e => e.Attribute("name")!.Value,
                            e => e.Value.Trim(),
                            StringComparer.OrdinalIgnoreCase)
                        ?? new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);

                    return new MethodDocInfo(summary, paramDescriptions);
                }
                catch
                {
                    return null;
                }
            }
        }

        return null;
    }

    private static string? FindSyntheticSourceRoot()
    {
        var current = new DirectoryInfo(AppContext.BaseDirectory);
        while (current != null)
        {
            string candidate = Path.Combine(current.FullName, "src", "synthetic");
            if (Directory.Exists(candidate))
                return candidate;
            current = current.Parent;
        }

        return null;
    }

    private sealed record MethodDocInfo(
        string? Summary,
        IReadOnlyDictionary<string, string> ParameterDescriptions);

    private sealed record SyntheticGeneratorDescriptor(
        string GeneratorName,
        string TypeName,
        MethodInfo GenerateMethod,
        MethodDocInfo? DocInfo);
}
