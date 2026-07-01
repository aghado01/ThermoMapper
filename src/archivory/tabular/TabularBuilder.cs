using System;
using System.Collections.Generic;
using System.Linq;

namespace Archivory.Tabular;

public abstract class TabularBuilder<TSelf>
    where TSelf : TabularBuilder<TSelf>
{
    private readonly List<string> _columnNames = new();
    private readonly List<object?[]> _rows = new();

    protected TabularBuilder()
    {
        Self = (TSelf)this;
    }

    protected TSelf Self { get; }

    public IReadOnlyList<string> ColumnNames => _columnNames;
    public IReadOnlyList<IReadOnlyList<object?>> Rows => _rows.Select(row => Array.AsReadOnly(row)).ToArray();

    public bool HasColumns => _columnNames.Count > 0;
    public int ColumnCount => _columnNames.Count;
    public int RowCount => _rows.Count;

    public TSelf WithColumns(params string[] names)
        => WithColumns((IEnumerable<string>)names);

    public TSelf WithColumns(IEnumerable<string> names)
    {
        if (names is null)
            throw new ArgumentNullException(nameof(names));

        var candidates = names.ToArray();
        ValidateColumns(candidates);
        if (_columnNames.Count > 0)
            throw new InvalidOperationException("Columns have already been defined.");

        _columnNames.AddRange(candidates);
        OnColumnsDefined(_columnNames);
        return Self;
    }

    public TSelf AddColumn(string name)
    {
        if (string.IsNullOrWhiteSpace(name))
            throw new ArgumentException("Column name cannot be null, empty, or whitespace.", nameof(name));

        if (_rows.Count > 0)
            throw new InvalidOperationException("Cannot add columns after rows have been added.");

        if (_columnNames.Contains(name, StringComparer.Ordinal))
            throw new ArgumentException("Duplicate column name.", nameof(name));

        _columnNames.Add(name);
        OnColumnsDefined(_columnNames);
        return Self;
    }

    public TSelf AddRow(params object?[] values)
        => AddRow((IEnumerable<object?>)values);

    public TSelf AddRow(IEnumerable<object?> values)
    {
        if (values is null)
            throw new ArgumentNullException(nameof(values));

        var row = values.ToArray();
        EnsureColumns(row.Length);
        _rows.Add(row);
        OnRowAdded(row);
        return Self;
    }

    public TSelf AddRow<T>(T item, params Func<T, object?>[] selectors)
        => AddRow(item, (IEnumerable<Func<T, object?>>)selectors);

    public TSelf AddRow<T>(T item, IEnumerable<Func<T, object?>> selectors)
    {
        if (selectors is null)
            throw new ArgumentNullException(nameof(selectors));

        return AddRow(selectors.Select(selector => selector(item)));
    }

    public TSelf AddRows<T>(IEnumerable<T> items, IEnumerable<Func<T, object?>> selectors)
    {
        if (items is null)
            throw new ArgumentNullException(nameof(items));

        if (selectors is null)
            throw new ArgumentNullException(nameof(selectors));

        foreach (T item in items)
            AddRow(item, selectors);

        return Self;
    }

    public TSelf Clear()
    {
        _columnNames.Clear();
        _rows.Clear();
        return Self;
    }

    public TabularData Build()
        => new TabularData(_columnNames.AsReadOnly(), _rows.Select(row => Array.AsReadOnly(row)).ToArray());

    protected virtual void OnColumnsDefined(IReadOnlyList<string> columnNames)
    {
    }

    protected virtual void OnRowAdded(object?[] row)
    {
    }

    private void EnsureColumns(int valueCount)
    {
        if (_columnNames.Count == 0)
            throw new InvalidOperationException("Define columns before adding rows.");

        if (valueCount != _columnNames.Count)
            throw new ArgumentException(
                $"Row has {valueCount} values but table has {_columnNames.Count} columns.",
                nameof(valueCount));
    }

    private static void ValidateColumns(string[] names)
    {
        if (names.Length == 0)
            throw new ArgumentException("At least one column must be defined.", nameof(names));

        if (names.Any(string.IsNullOrWhiteSpace))
            throw new ArgumentException("Column names cannot be null, empty, or whitespace.", nameof(names));

        if (names.Length != names.Distinct(StringComparer.Ordinal).Count())
            throw new ArgumentException("Column names must be unique.", nameof(names));
    }
}

public sealed class TabularBuilder : TabularBuilder<TabularBuilder>
{
}
