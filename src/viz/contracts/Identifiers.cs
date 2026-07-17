using System;

namespace Viz.Contracts;

internal static class IdentifierGuard
{
    internal static string Require(string value, string parameterName)
    {
        if (string.IsNullOrWhiteSpace(value))
            throw new ArgumentException("A visualization identifier cannot be empty.", parameterName);

        return value;
    }
}

/// <summary>An extensible semantic identifier such as <c>graph.csr</c>.</summary>
public readonly record struct SemanticId
{
    public SemanticId(string value) => Value = IdentifierGuard.Require(value, nameof(value));

    public string Value { get; }

    public override string ToString() => Value;
}

public readonly record struct StudyId
{
    public StudyId(string value) => Value = IdentifierGuard.Require(value, nameof(value));

    public string Value { get; }

    public override string ToString() => Value;
}

public readonly record struct EntitySetId
{
    public EntitySetId(string value) => Value = IdentifierGuard.Require(value, nameof(value));

    public string Value { get; }

    public override string ToString() => Value;
}

public readonly record struct ArtifactId
{
    public ArtifactId(string value) => Value = IdentifierGuard.Require(value, nameof(value));

    public string Value { get; }

    public override string ToString() => Value;
}

public readonly record struct CoordinateSpaceId
{
    public CoordinateSpaceId(string value) => Value = IdentifierGuard.Require(value, nameof(value));

    public string Value { get; }

    public override string ToString() => Value;
}

public readonly record struct PanelId
{
    public PanelId(string value) => Value = IdentifierGuard.Require(value, nameof(value));

    public string Value { get; }

    public override string ToString() => Value;
}

public readonly record struct RunId
{
    public RunId(string value) => Value = IdentifierGuard.Require(value, nameof(value));

    public string Value { get; }

    public override string ToString() => Value;
}

/// <summary>A stable reference to one entity within an entity set.</summary>
public readonly record struct EntityReference
{
    public EntityReference(EntitySetId entitySet, long index)
    {
        if (index < 0)
            throw new ArgumentOutOfRangeException(nameof(index), "An entity index cannot be negative.");

        EntitySet = entitySet;
        Index = index;
    }

    public EntitySetId EntitySet { get; }

    public long Index { get; }
}

/// <summary>Version of the durable Viz contract, independent of assembly version.</summary>
public readonly record struct ContractVersion
{
    public ContractVersion(int major, int minor)
    {
        if (major < 0) throw new ArgumentOutOfRangeException(nameof(major));
        if (minor < 0) throw new ArgumentOutOfRangeException(nameof(minor));

        Major = major;
        Minor = minor;
    }

    public int Major { get; }

    public int Minor { get; }

    public override string ToString() => $"{Major}.{Minor}";
}
