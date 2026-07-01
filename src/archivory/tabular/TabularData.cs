using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Archivory.Tabular;

public sealed class TabularData
{
    public TabularData(IReadOnlyList<string> columnNames, IReadOnlyList<IReadOnlyList<object?>> rows)
    {
        ColumnNames = columnNames ?? throw new ArgumentNullException(nameof(columnNames));
        Rows = rows ?? throw new ArgumentNullException(nameof(rows));

        if (ColumnNames.Any(string.IsNullOrWhiteSpace))
            throw new ArgumentException("Column names cannot be null, empty, or whitespace.", nameof(columnNames));

        if (Rows.Any(row => row.Count != ColumnNames.Count))
            throw new ArgumentException("Each row must have the same number of values as the column count.", nameof(rows));
    }

    public IReadOnlyList<string> ColumnNames { get; }
    public IReadOnlyList<IReadOnlyList<object?>> Rows { get; }

    public int ColumnCount => ColumnNames.Count;
    public int RowCount => Rows.Count;

    public void WriteCsv(TextWriter writer, char delimiter = ',', string? lineEnding = null)
    {
        if (writer is null)
            throw new ArgumentNullException(nameof(writer));

        string lineSeparator = lineEnding ?? Environment.NewLine;
        writer.Write(string.Join(delimiter, ColumnNames.Select(Escape)));
        writer.Write(lineSeparator);

        foreach (var row in Rows)
        {
            writer.Write(string.Join(delimiter, row.Select(Escape)));
            writer.Write(lineSeparator);
        }

        writer.Flush();

        string Escape(object? value)
        {
            string text = value?.ToString() ?? string.Empty;
            if (text.IndexOfAny(new[] { '"', '\r', '\n', delimiter }) >= 0)
                return '"' + text.Replace("\"", "\"\"") + '"';

            return text;
        }
    }

    public string ToCsv(char delimiter = ',', string? lineEnding = null)
    {
        using var writer = new StringWriter();
        WriteCsv(writer, delimiter, lineEnding);
        return writer.ToString();
    }
}
