using System;
using Archivory;

namespace Clustering.Graphical.SPC.Export;

/// <summary>
/// Binds an <see cref="SpcRunLayout"/> schema to a concrete <see cref="ArtifactScope"/>,
/// resolving the bundle's canonical files and sub-scopes. The code half of the
/// schema/code split — <see cref="SpcRunLayout"/> says <em>what</em>, this says
/// <em>where</em>. Supersedes the scattered path literals and side-effecting getters
/// of the former <c>SpcOutputPathHelper</c> with one scope-rooted owner threaded
/// down from the top of a run.
/// </summary>
public sealed class SpcRunPaths
{
    private readonly ArtifactScope _scope;
    private readonly SpcRunLayout _layout;

    public SpcRunPaths(ArtifactScope scope, SpcRunLayout? layout = null)
    {
        _scope = scope ?? throw new ArgumentNullException(nameof(scope));
        _layout = layout ?? new SpcRunLayout();
    }

    /// <summary>The run's root scope (the <c>{family}/{stamp}</c> directory).</summary>
    public ArtifactScope Scope => _scope;

    // ── Root metadata files ────────────────────────────────────────────────
    public string Manifest  => _scope.File(_layout.Manifest);
    public string Summary   => _scope.File(_layout.Summary);
    public string Health    => _scope.File(_layout.Health);
    public string Hierarchy => _scope.File(_layout.Hierarchy);

    // ── Bulk sub-scopes ────────────────────────────────────────────────────
    // Each is a child scope, not a path: the writer calls EnsureDirectory() and
    // composes filenames from their owner (e.g. Csv.File(SpcCsvWriter.SweepFileName),
    // or hands Csv.Dir to SpcCsvWriter.WriteAll). Checkpoints feeds the executor's
    // checkpoint directory; SPCX/SPCE filenames stay with SpcExecutor.
    public ArtifactScope Csv         => _scope.Child(_layout.Csv);
    public ArtifactScope Tabular     => _scope.Child(_layout.Tabular);
    public ArtifactScope Checkpoints => _scope.Child(_layout.Checkpoints);
}
