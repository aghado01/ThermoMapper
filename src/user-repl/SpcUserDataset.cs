using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Clustering.Graphical.SPC.Export;
using Synthetic;
using Archivory.Tabular;

namespace UserRepl;

public sealed record SpcUserDataset(
    double[][] Features,
    int[] Labels,
    int ClusterCount,
    int[][]? LabelsByLevel,
    IReadOnlyDictionary<string, object?> Metadata)
{
    public static SpcUserDataset FromSyntheticDataset(SyntheticDataset dataset, IReadOnlyDictionary<string, object?>? extraMetadata = null)
    {
        if (dataset is null)
            throw new ArgumentNullException(nameof(dataset));

        var metadata = new Dictionary<string, object?>(StringComparer.OrdinalIgnoreCase);

        if (dataset.Metadata != null)
        {
            metadata["Generator"] = dataset.Metadata.GeneratorName;
            metadata["TopologyTag"] = dataset.Metadata.TopologyTag;
            metadata["GeometryClass"] = dataset.Metadata.GeometryClass;
            metadata["HierarchyTag"] = dataset.Metadata.HierarchyTag;
            metadata["GTNumClusters"] = dataset.Metadata.GTNumClusters;
            metadata["AmbientDimensionality"] = dataset.Metadata.AmbientDimensionality;
            if (!string.IsNullOrWhiteSpace(dataset.Metadata.LiteratureReference))
                metadata["LiteratureReference"] = dataset.Metadata.LiteratureReference;
            if (!string.IsNullOrWhiteSpace(dataset.Metadata.SuggestedMetric))
                metadata["SuggestedMetric"] = dataset.Metadata.SuggestedMetric;
            if (!string.IsNullOrWhiteSpace(dataset.Metadata.FutureMetric))
                metadata["FutureMetric"] = dataset.Metadata.FutureMetric;
        }

        if (extraMetadata is not null)
        {
            foreach (var kvp in extraMetadata)
                metadata[kvp.Key] = kvp.Value;
        }

        return new SpcUserDataset(
            dataset.Features,
            dataset.Labels,
            dataset.ClusterCount,
            dataset.LabelsByLevel,
            metadata);
    }

    public static SpcUserDataset FromFeatures(
        double[][] features,
        int[] labels,
        int clusterCount = 0,
        int[][]? labelsByLevel = null,
        IReadOnlyDictionary<string, object?>? metadata = null)
    {
        if (features is null)
            throw new ArgumentNullException(nameof(features));
        if (labels is null)
            throw new ArgumentNullException(nameof(labels));
        if (features.Length == 0)
            throw new ArgumentException("Features cannot be empty.", nameof(features));
        if (labels.Length != features.Length)
            throw new ArgumentException("Labels length must match feature row count.", nameof(labels));

        int dimension = features[0].Length;
        for (int i = 0; i < features.Length; i++)
        {
            if (features[i] is null)
                throw new ArgumentException($"Feature row {i} is null.", nameof(features));
            if (features[i].Length != dimension)
                throw new ArgumentException("All feature rows must have the same dimensionality.", nameof(features));
        }

        var metadataValues = new Dictionary<string, object?>(StringComparer.OrdinalIgnoreCase);
        if (metadata is not null)
        {
            foreach (var kvp in metadata)
                metadataValues[kvp.Key] = kvp.Value;
        }

        int inferredClusterCount = clusterCount > 0
            ? clusterCount
            : labels.Distinct().Count();

        if (inferredClusterCount <= 0)
            inferredClusterCount = labels.Distinct().Count();

        return new SpcUserDataset(
            features,
            labels,
            inferredClusterCount,
            labelsByLevel,
            metadataValues);
    }

    public TabularProjection ToTabular(string tableName = "spc_dataset")
        => SpcTabularProjections.CreateDatasetProjection(Features, Labels, tableName);

    public static SpcUserDataset FromCsv(
        string path,
        string? labelColumn = null,
        bool hasHeader = true,
        char separator = ',')
    {
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path must be provided.", nameof(path));
        if (!File.Exists(path))
            throw new FileNotFoundException("CSV file not found.", path);

        string[] rawLines = File.ReadAllLines(path);
        var lines = rawLines
            .Select(line => line.Trim())
            .Where(line => !string.IsNullOrWhiteSpace(line) && !line.TrimStart().StartsWith("#", StringComparison.Ordinal))
            .ToArray();

        if (lines.Length == 0)
            throw new ArgumentException("CSV file is empty or contains only comments/blank lines.", nameof(path));

        string[]? headers = null;
        int firstDataLine = 0;
        if (hasHeader)
        {
            headers = ParseCsvLine(lines[0], separator);
            if (headers.Length == 0)
                throw new ArgumentException("CSV header row is empty.", nameof(path));
            firstDataLine = 1;
            if (lines.Length == 1)
                throw new ArgumentException("CSV file contains header only and no data rows.", nameof(path));
        }

        var rows = new List<string[]>(lines.Length - firstDataLine);
        for (int i = firstDataLine; i < lines.Length; i++)
        {
            var fields = ParseCsvLine(lines[i], separator);
            if (fields.Length == 0)
                continue;
            rows.Add(fields);
        }

        if (rows.Count == 0)
            throw new ArgumentException("CSV file contains no valid data rows.", nameof(path));

        int fieldCount = rows[0].Length;
        if (rows.Any(row => row.Length != fieldCount))
            throw new ArgumentException("CSV rows have inconsistent field counts.", nameof(path));

        int labelIndex = ResolveLabelIndex(labelColumn, headers, fieldCount);

        var featureRows = new List<double[]>(rows.Count);
        var labelValues = new List<string>(rows.Count);

        for (int rowIndex = 0; rowIndex < rows.Count; rowIndex++)
        {
            string[] row = rows[rowIndex];
            if (row.Length != fieldCount)
                throw new ArgumentException($"CSV row {rowIndex + firstDataLine + 1} has an unexpected number of fields.", nameof(path));

            if (labelIndex >= 0)
            {
                labelValues.Add(row[labelIndex].Trim());
            }

            var featureRow = new double[fieldCount - (labelIndex >= 0 ? 1 : 0)];
            int featureIndex = 0;
            for (int columnIndex = 0; columnIndex < fieldCount; columnIndex++)
            {
                if (columnIndex == labelIndex)
                    continue;

                string valueText = row[columnIndex].Trim();
                if (!double.TryParse(valueText, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out double value))
                    throw new FormatException($"Unable to parse CSV value '{valueText}' as a number at row {rowIndex + firstDataLine + 1}, column {columnIndex + 1}.");

                featureRow[featureIndex++] = value;
            }

            featureRows.Add(featureRow);
        }

        int[] labels = ConvertLabelValues(labelValues);
        var metadataValues = new Dictionary<string, object?>(StringComparer.OrdinalIgnoreCase)
        {
            ["Source"] = "csv",
            ["CsvPath"] = Path.GetFullPath(path),
            ["HasHeader"] = hasHeader,
            ["Separator"] = separator.ToString(),
        };

        if (labelIndex >= 0)
        {
            metadataValues["LabelColumn"] = labelColumn ?? (headers is not null ? headers[^1] : (fieldCount - 1).ToString());
            if (labelValues.Count > 0 && labels.Length > 0)
                metadataValues["LabelCardinality"] = labels.Distinct().Count();
        }

        return FromFeatures(featureRows.ToArray(), labels, clusterCount: labels.Distinct().Count(), metadata: metadataValues);
    }

    private static int ResolveLabelIndex(string? labelColumn, string[]? headers, int fieldCount)
    {
        if (!string.IsNullOrWhiteSpace(labelColumn))
        {
            if (int.TryParse(labelColumn, NumberStyles.Integer, CultureInfo.InvariantCulture, out int index))
            {
                if (index < 0 || index >= fieldCount)
                    throw new ArgumentException($"Label column index {index} is outside the valid range [0, {fieldCount - 1}].", nameof(labelColumn));
                return index;
            }

            if (headers is null)
                throw new ArgumentException("Label column name cannot be resolved without a CSV header.", nameof(labelColumn));

            int foundIndex = Array.FindIndex(headers, header => string.Equals(header, labelColumn, StringComparison.OrdinalIgnoreCase));
            if (foundIndex < 0)
                throw new ArgumentException($"Label column '{labelColumn}' was not found in the CSV header.", nameof(labelColumn));

            return foundIndex;
        }

        return fieldCount - 1;
    }

    private static int[] ConvertLabelValues(IReadOnlyList<string> rawLabels)
    {
        if (rawLabels is null)
            throw new ArgumentNullException(nameof(rawLabels));

        if (rawLabels.Count == 0)
            return Array.Empty<int>();

        var values = new int[rawLabels.Count];
        bool allIntegers = true;

        for (int i = 0; i < rawLabels.Count; i++)
        {
            if (!int.TryParse(rawLabels[i], NumberStyles.Integer, CultureInfo.InvariantCulture, out int parsed))
            {
                allIntegers = false;
                break;
            }

            values[i] = parsed;
        }

        if (allIntegers)
            return values;

        var map = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);
        int nextLabel = 0;
        for (int i = 0; i < rawLabels.Count; i++)
        {
            string key = rawLabels[i];
            if (!map.TryGetValue(key, out int code))
            {
                code = nextLabel++;
                map[key] = code;
            }
            values[i] = code;
        }

        return values;
    }

    private static string[] ParseCsvLine(string line, char separator)
    {
        if (line is null)
            throw new ArgumentNullException(nameof(line));

        var fields = new List<string>();
        var builder = new System.Text.StringBuilder();
        bool inQuotes = false;

        for (int i = 0; i < line.Length; i++)
        {
            char c = line[i];
            if (inQuotes)
            {
                if (c == '"')
                {
                    if (i + 1 < line.Length && line[i + 1] == '"')
                    {
                        builder.Append('"');
                        i++;
                    }
                    else
                    {
                        inQuotes = false;
                    }
                }
                else
                {
                    builder.Append(c);
                }
            }
            else if (c == '"')
            {
                inQuotes = true;
            }
            else if (c == separator)
            {
                fields.Add(builder.ToString());
                builder.Clear();
            }
            else
            {
                builder.Append(c);
            }
        }

        fields.Add(builder.ToString());
        return fields.ToArray();
    }
}
