#nullable enable
using System;
using System.Collections.Generic;

namespace Hashish;

/// <summary>
/// Minimal dispatch surface for pairwise measures that expose both a distance
/// and a similarity view over the same input domain.
/// </summary>
public interface IMeasure<T>
{
    double Distance(T a, T b);
    double Similarity(T a, T b);
}

/// <summary>Dispatch adapter for <see cref="Levenshtein"/> over strings.</summary>
public readonly struct LevenshteinMeasure : IMeasure<string?>
{
    public double Distance(string? a, string? b)
        => Levenshtein.Distance(a, b);

    public double Similarity(string? a, string? b)
        => Levenshtein.Similarity(a, b);
}

/// <summary>Dispatch adapter for <see cref="CosineVectors"/> over dense vectors.</summary>
public readonly struct CosineVectorMeasure : IMeasure<double[]>
{
    public double Distance(double[] a, double[] b)
        => CosineVectors.Distance(a, b);

    public double Similarity(double[] a, double[] b)
        => CosineVectors.Similarity(a, b);
}

/// <summary>Dispatch adapter for exact Jaccard set overlap.</summary>
public readonly struct JaccardMeasure<T> : IMeasure<IEnumerable<T>>
    where T : notnull
{
    private readonly IEqualityComparer<T>? _comparer;

    public JaccardMeasure(IEqualityComparer<T>? comparer = null)
        => _comparer = comparer;

    public double Distance(IEnumerable<T> a, IEnumerable<T> b)
        => JaccardContainment.Distance(a, b, _comparer);

    public double Similarity(IEnumerable<T> a, IEnumerable<T> b)
        => JaccardContainment.Similarity(a, b, _comparer);
}

/// <summary>Dispatch adapter for exact Sorensen-Dice set overlap.</summary>
public readonly struct DiceMeasure<T> : IMeasure<IEnumerable<T>>
    where T : notnull
{
    private readonly IEqualityComparer<T>? _comparer;

    public DiceMeasure(IEqualityComparer<T>? comparer = null)
        => _comparer = comparer;

    public double Distance(IEnumerable<T> a, IEnumerable<T> b)
        => JaccardContainment.DiceDistance(a, b, _comparer);

    public double Similarity(IEnumerable<T> a, IEnumerable<T> b)
        => JaccardContainment.DiceSimilarity(a, b, _comparer);
}
