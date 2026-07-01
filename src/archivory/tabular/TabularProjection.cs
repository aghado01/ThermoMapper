using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Archivory.Tabular;

public interface ITabularProjection
{
    string TableName { get; }

    IReadOnlyList<string> Columns { get; }

    IEnumerable<IReadOnlyList<object?>> Rows { get; }
}

public sealed class TabularProjection : ITabularProjection
{
    public TabularProjection(
        string tableName,
        IReadOnlyList<string> columns,
        IEnumerable<IReadOnlyList<object?>> rows)
    {
        if (string.IsNullOrWhiteSpace(tableName))
            throw new ArgumentException("Table name cannot be null, empty, or whitespace.", nameof(tableName));

        Columns = columns ?? throw new ArgumentNullException(nameof(columns));
        if (Columns.Count == 0)
            throw new ArgumentException("At least one column must be provided.", nameof(columns));

        if (Columns.Any(string.IsNullOrWhiteSpace))
            throw new ArgumentException("Column names cannot be null, empty, or whitespace.", nameof(columns));

        if (Columns.Count != Columns.Distinct(StringComparer.Ordinal).Count())
            throw new ArgumentException("Column names must be unique.", nameof(columns));

        if (rows is null)
            throw new ArgumentNullException(nameof(rows));

        var rowArray = rows.ToArray();
        if (rowArray.Any(row => row.Count != Columns.Count))
            throw new ArgumentException("Each row must have the same number of values as the column count.", nameof(rows));

        TableName = tableName;
        Rows = Array.AsReadOnly(rowArray);
    }

    public string TableName { get; }

    public IReadOnlyList<string> Columns { get; }

    public IEnumerable<IReadOnlyList<object?>> Rows { get; }
}

public static class TabularProjectionExtensions
{
    public static TabularData ToTabularData(this ITabularProjection projection)
    {
        if (projection is null)
            throw new ArgumentNullException(nameof(projection));

        return new TabularData(projection.Columns, projection.Rows.ToArray());
    }

    public static void WriteCsv(this ITabularProjection projection, TextWriter writer, char delimiter = ',', string? lineEnding = null)
    {
        if (projection is null)
            throw new ArgumentNullException(nameof(projection));

        projection.ToTabularData().WriteCsv(writer, delimiter, lineEnding);
    }

    public static string ToCsv(this ITabularProjection projection, char delimiter = ',', string? lineEnding = null)
    {
        if (projection is null)
            throw new ArgumentNullException(nameof(projection));

        return projection.ToTabularData().ToCsv(delimiter, lineEnding);
    }

    /// <summary>
    /// Atomically write the projection to <paramref name="path"/>. Writes to
    /// <c>{path}.tmp</c> first, then renames; a partial write never leaves a
    /// half-written file at the canonical path. Existing files at
    /// <paramref name="path"/> are overwritten on success.
    /// </summary>
    public static void WriteToFile(this ITabularProjection projection, string path, char delimiter = ',', string? lineEnding = null)
    {
        if (projection is null)
            throw new ArgumentNullException(nameof(projection));
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path cannot be null or whitespace.", nameof(path));

        string tmpPath = path + ".tmp";
        using (var writer = new StreamWriter(tmpPath, append: false))
        {
            projection.WriteCsv(writer, delimiter, lineEnding);
        }
        File.Move(tmpPath, path, overwrite: true);
    }
}

public static class TabularProjectionFactory
{
    public static TabularProjection Create<T>(
        string tableName,
        IEnumerable<string> columns,
        IEnumerable<T> items,
        params Func<T, object?>[] selectors)
        => Create(tableName, columns, items, (IEnumerable<Func<T, object?>>)selectors);

    public static TabularProjection Create<T>(
        string tableName,
        IEnumerable<string> columns,
        IEnumerable<T> items,
        IEnumerable<Func<T, object?>> selectors)
    {
        if (string.IsNullOrWhiteSpace(tableName))
            throw new ArgumentException("Table name cannot be null, empty, or whitespace.", nameof(tableName));

        if (columns is null)
            throw new ArgumentNullException(nameof(columns));

        if (selectors is null)
            throw new ArgumentNullException(nameof(selectors));

        var columnList = columns.ToArray();
        if (columnList.Length == 0)
            throw new ArgumentException("At least one column must be provided.", nameof(columns));

        if (columnList.Any(string.IsNullOrWhiteSpace))
            throw new ArgumentException("Column names cannot be null, empty, or whitespace.", nameof(columns));

        if (columnList.Length != columnList.Distinct(StringComparer.Ordinal).Count())
            throw new ArgumentException("Column names must be unique.", nameof(columns));

        var selectorList = selectors.ToArray();
        if (selectorList.Length != columnList.Length)
            throw new ArgumentException("The number of selectors must match the number of columns.", nameof(selectors));

        if (items is null)
            throw new ArgumentNullException(nameof(items));

        var rows = items.Select(item => BuildRow(item, selectorList)).ToArray();
        return new TabularProjection(tableName, Array.AsReadOnly(columnList), rows);
    }

    private static IReadOnlyList<object?> BuildRow<T>(T item, Func<T, object?>[] selectors)
    {
        var row = new object?[selectors.Length];
        for (int i = 0; i < selectors.Length; i++)
            row[i] = selectors[i](item);

        return Array.AsReadOnly(row);
    }

    /// <summary>
    /// Build a projection by iterating <paramref name="rowCount"/> indices
    /// and applying per-index column selectors. Suited to parallel-array
    /// shapes where each "row" corresponds to an index into multiple
    /// columnar arrays (e.g. a sweep profile with <c>Temperatures[i]</c>,
    /// <c>Susceptibility[i]</c>, <c>SpecificHeat[i]</c>).
    /// </summary>
    public static IndexedProjectionBuilder CreateIndexed(string tableName, int rowCount)
    {
        if (string.IsNullOrWhiteSpace(tableName))
            throw new ArgumentException("Table name cannot be null, empty, or whitespace.", nameof(tableName));
        if (rowCount < 0)
            throw new ArgumentOutOfRangeException(nameof(rowCount), "Row count cannot be negative.");

        return new IndexedProjectionBuilder(tableName, rowCount);
    }

    /// <summary>
    /// Build a single-row projection by assembling named columns fluently.
    /// Suited to scalar summaries (criteria, session totals, etc.) where
    /// the source is a single instance rather than an iterable.
    /// </summary>
    public static ScalarProjectionBuilder CreateScalar(string tableName)
    {
        if (string.IsNullOrWhiteSpace(tableName))
            throw new ArgumentException("Table name cannot be null, empty, or whitespace.", nameof(tableName));

        return new ScalarProjectionBuilder(tableName);
    }
}

public static class EnumerableTabularProjectionExtensions
{
    public static TabularProjection ProjectToTabular<T>(
        this IEnumerable<T> items,
        string tableName,
        IEnumerable<string> columns,
        params Func<T, object?>[] selectors)
    {
        return TabularProjectionFactory.Create(tableName, columns, items, selectors);
    }
}

/// <summary>
/// Fluent builder for a per-index <see cref="TabularProjection"/>. Use when
/// "row" means "an index into multiple parallel data arrays" rather than
/// "an item in a collection."
/// </summary>
/// <remarks>
/// <para>Composes three column families:</para>
/// <list type="bullet">
///   <item><see cref="Column"/> — always-on named column with a per-index selector.</item>
///   <item><see cref="ColumnIf"/> — conditional column included only when the
///     predicate is true (e.g. an optional aggregate that may or may not
///     be present on the source).</item>
///   <item><see cref="ColumnsFromDictionary{TValue}"/> — flatten a dictionary
///     into one column per key. Keys are sorted by the supplied comparer
///     (ordinal by default) so the column order is deterministic across runs.</item>
/// </list>
/// </remarks>
public sealed class IndexedProjectionBuilder
{
    private readonly string _tableName;
    private readonly int _rowCount;
    private readonly List<(string Name, Func<int, object?> Selector)> _columns = new();

    internal IndexedProjectionBuilder(string tableName, int rowCount)
    {
        _tableName = tableName;
        _rowCount = rowCount;
    }

    public IndexedProjectionBuilder Column(string name, Func<int, object?> selector)
    {
        if (string.IsNullOrWhiteSpace(name))
            throw new ArgumentException("Column name cannot be null, empty, or whitespace.", nameof(name));
        if (selector is null)
            throw new ArgumentNullException(nameof(selector));
        if (_columns.Any(c => string.Equals(c.Name, name, StringComparison.Ordinal)))
            throw new ArgumentException($"Duplicate column name '{name}'.", nameof(name));

        _columns.Add((name, selector));
        return this;
    }

    public IndexedProjectionBuilder ColumnIf(bool include, string name, Func<int, object?> selector)
        => include ? Column(name, selector) : this;

    /// <summary>
    /// Bulk-add named columns. Equivalent to calling <see cref="Column"/>
    /// for each entry but more convenient when the column list is itself
    /// produced from a loop (e.g. one column per feature dimension).
    /// </summary>
    public IndexedProjectionBuilder Columns(IEnumerable<(string Name, Func<int, object?> Selector)> columnSpecs)
    {
        if (columnSpecs is null)
            throw new ArgumentNullException(nameof(columnSpecs));

        foreach (var (name, selector) in columnSpecs)
            Column(name, selector);

        return this;
    }

    /// <summary>
    /// Append one column per key in <paramref name="dictionary"/>. The
    /// <paramref name="valueAccessor"/> takes the dictionary value plus the
    /// current row index and returns the cell value (useful when the
    /// dictionary's values are themselves per-index sequences). Keys are
    /// sorted by <paramref name="keyOrder"/> (ordinal by default) for
    /// deterministic column ordering.
    /// </summary>
    public IndexedProjectionBuilder ColumnsFromDictionary<TValue>(
        IReadOnlyDictionary<string, TValue> dictionary,
        Func<TValue, int, object?> valueAccessor,
        IComparer<string>? keyOrder = null)
    {
        if (dictionary is null)
            throw new ArgumentNullException(nameof(dictionary));
        if (valueAccessor is null)
            throw new ArgumentNullException(nameof(valueAccessor));

        var ordered = dictionary.Keys.OrderBy(k => k, keyOrder ?? StringComparer.Ordinal);
        foreach (var key in ordered)
        {
            TValue value = dictionary[key];
            Column(key, i => valueAccessor(value, i));
        }
        return this;
    }

    public TabularProjection Build()
    {
        if (_columns.Count == 0)
            throw new InvalidOperationException("At least one column must be defined before Build().");

        var columnNames = _columns.Select(c => c.Name).ToArray();
        var selectors = _columns.Select(c => c.Selector).ToArray();

        var rows = new IReadOnlyList<object?>[_rowCount];
        for (int i = 0; i < _rowCount; i++)
        {
            var row = new object?[columnNames.Length];
            for (int c = 0; c < columnNames.Length; c++)
                row[c] = selectors[c](i);
            rows[i] = Array.AsReadOnly(row);
        }

        return new TabularProjection(_tableName, Array.AsReadOnly(columnNames), rows);
    }
}

/// <summary>
/// Fluent builder for a single-row <see cref="TabularProjection"/>. Use for
/// scalar/summary outputs where each column carries one value from a single
/// source instance.
/// </summary>
public sealed class ScalarProjectionBuilder
{
    private readonly string _tableName;
    private readonly List<(string Name, object? Value)> _columns = new();

    internal ScalarProjectionBuilder(string tableName)
    {
        _tableName = tableName;
    }

    public ScalarProjectionBuilder Column(string name, object? value)
    {
        if (string.IsNullOrWhiteSpace(name))
            throw new ArgumentException("Column name cannot be null, empty, or whitespace.", nameof(name));
        if (_columns.Any(c => string.Equals(c.Name, name, StringComparison.Ordinal)))
            throw new ArgumentException($"Duplicate column name '{name}'.", nameof(name));

        _columns.Add((name, value));
        return this;
    }

    public ScalarProjectionBuilder ColumnIf(bool include, string name, object? value)
        => include ? Column(name, value) : this;

    /// <summary>
    /// Append one column per key in <paramref name="dictionary"/>. Keys are
    /// sorted by <paramref name="keyOrder"/> (ordinal by default) for
    /// deterministic column ordering.
    /// </summary>
    public ScalarProjectionBuilder ColumnsFromDictionary<TValue>(
        IReadOnlyDictionary<string, TValue> dictionary,
        IComparer<string>? keyOrder = null)
    {
        if (dictionary is null)
            throw new ArgumentNullException(nameof(dictionary));

        var ordered = dictionary.Keys.OrderBy(k => k, keyOrder ?? StringComparer.Ordinal);
        foreach (var key in ordered)
            Column(key, dictionary[key]);

        return this;
    }

    public TabularProjection Build()
    {
        if (_columns.Count == 0)
            throw new InvalidOperationException("At least one column must be defined before Build().");

        var columnNames = _columns.Select(c => c.Name).ToArray();
        var row = _columns.Select(c => c.Value).ToArray();

        return new TabularProjection(_tableName, Array.AsReadOnly(columnNames), new[] { Array.AsReadOnly(row) });
    }
}
