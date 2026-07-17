using System;

namespace Viz.Scene;

/// <summary>Stable identity of a visual layer within a panel scene.</summary>
public readonly record struct LayerId
{
    public LayerId(string value)
    {
        if (string.IsNullOrWhiteSpace(value))
            throw new ArgumentException("A layer identifier cannot be empty.", nameof(value));

        Value = value;
    }

    public string Value { get; }

    public override string ToString() => Value;
}
