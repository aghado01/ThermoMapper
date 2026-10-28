using System;
using Archivory;

namespace Clustering.Graphical.CumulativeField.Export;

/// <summary>
/// Binds a <see cref="CumulativeFieldRunLayout"/> schema to a concrete
/// <see cref="ArtifactScope"/>, resolving the bundle's canonical files and
/// sub-scopes — the layout says <em>what</em>, this says <em>where</em>.
/// </summary>
public sealed class CumulativeFieldRunPaths
{
    private readonly ArtifactScope _scope;
    private readonly CumulativeFieldRunLayout _layout;

    public CumulativeFieldRunPaths(ArtifactScope scope, CumulativeFieldRunLayout? layout = null)
    {
        _scope = scope ?? throw new ArgumentNullException(nameof(scope));
        _layout = layout ?? new CumulativeFieldRunLayout();
    }

    /// <summary>The run's root scope (the <c>{family}/{stamp}</c> directory).</summary>
    public ArtifactScope Scope => _scope;

    // ── Root metadata files ────────────────────────────────────────────────
    public string Manifest => _scope.File(_layout.Manifest);
    public string Summary  => _scope.File(_layout.Summary);
    public string Health   => _scope.File(_layout.Health);

    // ── Bulk sub-scopes ────────────────────────────────────────────────────
    // A child scope, not a path: the writer calls EnsureDirectory() and composes
    // filenames from their owner (Csv.File(CumulativeFieldCsvWriter.SweepFileName)).
    public ArtifactScope Csv => _scope.Child(_layout.Csv);
}
