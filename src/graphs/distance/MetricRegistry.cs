using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Graphs.Distance.Euclidean;
using Graphs.Distance.Geodesic;

namespace Graphs.Distance;

/// <summary>
/// Canonical identity of a distance metric. The single enum every layer agrees
/// on — parsing, dispatch, and the discoverable list all key off this.
/// </summary>
public enum MetricKind
{
    Euclidean,
    Manhattan,
    Minkowski,
    Hamming,
    Poincare,
    Cosine,
}

/// <summary>Parsed metric request: a <see cref="MetricKind"/> plus the only
/// parameter any metric currently carries (the Minkowski exponent; ignored for
/// the others).</summary>
public readonly record struct MetricSpec(MetricKind Kind, double MinkowskiP = 2.0);

/// <summary>
/// Self-describing entry for one metric — what a UI dropdown enumerates so that
/// only properly-integrated metrics are ever selectable, and so it can render a
/// parameter field (and, via <see cref="MetricRegistry.GetProperties(MetricSpec)"/>,
/// validate compatibility like "requires unit-norm").
/// </summary>
public sealed record MetricDescriptor(
    MetricKind            Kind,
    string                Id,
    string                DisplayName,
    IReadOnlyList<string> Aliases,
    bool                  RequiresParameter,
    string?               ParameterName);

/// <summary>
/// CPS/visitor seam that lets a caller act on the concrete struct metric type
/// without the registry knowing what the caller does with it. The struct-generic
/// <c>TMetric</c> (a type parameter of <see cref="Visit{TMetric}"/>) is what
/// preserves JIT inlining for hot paths like HDBSCAN's Prim loop
/// (<c>runner.Run&lt;TMetric&gt;</c>) — an <see cref="IDistanceMetric"/> interface
/// instance could not.
/// </summary>
public interface IMetricVisitor<TResult>
{
    TResult Visit<TMetric>(TMetric metric) where TMetric : struct, IDistanceMetric;
}

/// <summary>
/// The one place the metric vocabulary lives: the kind ↔ spelling list, the spec
/// parser, and the kind → concrete-struct switch. Every consumer derives from
/// this instead of re-listing the families, so they cannot drift —
/// <c>DistanceMetricFactory</c> (interface path) and HDBSCAN's struct-generic
/// dispatch are both thin clients, and a future VizCore dropdown enumerates
/// <see cref="Available"/>.
/// </summary>
public static class MetricRegistry
{
    /// <summary>Discoverable metric set — the authoritative list of what's
    /// integrated. UI dropdowns enumerate this.</summary>
    public static IReadOnlyList<MetricDescriptor> Available { get; } = new MetricDescriptor[]
    {
        new(MetricKind.Euclidean, "euclidean", "Euclidean (L2)",            Array.Empty<string>(), false, null),
        new(MetricKind.Manhattan, "manhattan", "Manhattan (L1)",           new[] { "l1" },        false, null),
        new(MetricKind.Minkowski, "minkowski", "Minkowski (Lᵖ)",           Array.Empty<string>(), true,  "p"),
        new(MetricKind.Hamming,   "hamming",   "Hamming (symbol)",          Array.Empty<string>(), false, null),
        new(MetricKind.Poincare,  "poincare",  "Poincaré ball (geodesic)", Array.Empty<string>(), false, null),
        new(MetricKind.Cosine,    "cosine",    "Cosine (spherical geodesic)", Array.Empty<string>(), false, null),
    };

    /// <summary>
    /// Parse a spec string: <c>euclidean</c> | <c>manhattan</c>/<c>l1</c> |
    /// <c>minkowski:p=N</c> (or <c>minkowski:N</c>) | <c>hamming</c> |
    /// <c>poincare</c> | <c>cosine</c>. Case-insensitive.
    /// </summary>
    public static MetricSpec Parse(string spec)
    {
        if (string.IsNullOrWhiteSpace(spec))
            throw new ArgumentException("Distance metric spec must be non-empty.", nameof(spec));

        string trimmed = spec.Trim();
        int colon = trimmed.IndexOf(':');
        string head = (colon > 0 ? trimmed[..colon] : trimmed).Trim().ToLowerInvariant();
        string? tail = colon > 0 ? trimmed[(colon + 1)..] : null;

        return head switch
        {
            "euclidean"         => new MetricSpec(MetricKind.Euclidean),
            "manhattan" or "l1" => new MetricSpec(MetricKind.Manhattan),
            "minkowski"         => new MetricSpec(MetricKind.Minkowski, ParseMinkowskiExponent(tail)),
            "hamming"           => new MetricSpec(MetricKind.Hamming),
            "poincare"          => new MetricSpec(MetricKind.Poincare),
            "cosine"            => new MetricSpec(MetricKind.Cosine),
            _ => throw new ArgumentException(
                $"Unknown distance metric '{spec}'. Valid: {ValidList()}."),
        };
    }

    /// <summary>The single kind → concrete-struct switch. Hands the struct to
    /// <paramref name="visitor"/> so the caller specializes on the real type.</summary>
    public static TResult Invoke<TResult>(MetricSpec spec, IMetricVisitor<TResult> visitor)
    {
        if (visitor is null) throw new ArgumentNullException(nameof(visitor));
        return spec.Kind switch
        {
            MetricKind.Euclidean => visitor.Visit(default(EuclideanMetric)),
            MetricKind.Manhattan => visitor.Visit(default(ManhattanMetric)),
            MetricKind.Minkowski => visitor.Visit(new MinkowskiMetric(spec.MinkowskiP)),
            MetricKind.Hamming   => visitor.Visit(default(HammingMetric)),
            MetricKind.Poincare  => visitor.Visit(default(PoincareMetric)),
            MetricKind.Cosine    => visitor.Visit(default(SphericalGeodesicMetric)),
            _ => throw new ArgumentOutOfRangeException(nameof(spec), spec.Kind, "Unhandled metric kind."),
        };
    }

    /// <summary>Boxing path: a concrete <see cref="IDistanceMetric"/> instance for
    /// callers that don't need struct specialization (the SPC graph builder, etc.).
    /// Replaces the old <c>DistanceMetricFactory</c> switch.</summary>
    public static IDistanceMetric Create(MetricSpec spec) => Invoke(spec, BoxingVisitor.Instance);

    public static IDistanceMetric Create(string spec) => Create(Parse(spec));

    /// <summary>Metric properties (bandwidth strategy, unit-norm / probability
    /// requirements) for the given spec — derived through the same dispatch, so a
    /// dropdown can validate compatibility without a second switch.</summary>
    public static MetricProperties GetProperties(MetricSpec spec) => Invoke(spec, PropertiesVisitor.Instance);

    public static MetricProperties GetProperties(string spec) => GetProperties(Parse(spec));

    private static double ParseMinkowskiExponent(string? tail)
    {
        // Accept "minkowski:p=1.5" or "minkowski:1.5".
        if (string.IsNullOrWhiteSpace(tail))
            throw new ArgumentException("minkowski requires an exponent, e.g. minkowski:p=1.5.");

        string body = tail.Trim();
        int eq = body.IndexOf('=');
        if (eq >= 0) body = body[(eq + 1)..].Trim();

        if (!double.TryParse(body, NumberStyles.Float, CultureInfo.InvariantCulture, out double p))
            throw new ArgumentException($"Cannot parse Minkowski exponent: '{tail}'.");

        return p;
    }

    private static string ValidList() => string.Join(", ",
        Available.Select(d => d.RequiresParameter ? $"{d.Id}:{d.ParameterName}=N" : d.Id));

    private sealed class BoxingVisitor : IMetricVisitor<IDistanceMetric>
    {
        public static readonly BoxingVisitor Instance = new();
        public IDistanceMetric Visit<TMetric>(TMetric metric) where TMetric : struct, IDistanceMetric => metric;
    }

    private sealed class PropertiesVisitor : IMetricVisitor<MetricProperties>
    {
        public static readonly PropertiesVisitor Instance = new();
        public MetricProperties Visit<TMetric>(TMetric metric) where TMetric : struct, IDistanceMetric => metric.Properties;
    }
}
