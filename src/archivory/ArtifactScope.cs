using System;
using System.IO;

namespace Archivory;

/// <summary>
/// An immutable node in a run's artifact tree — the "matryoshka doll". A scope
/// names a directory; <see cref="Child"/> nests a role-named sub-scope and
/// <see cref="File"/> resolves a canonical file beneath it. The same writing
/// code can target a <see cref="Root"/> scope (a standalone run) or a
/// <see cref="Child"/> scope (an arm of a larger run) without knowing which —
/// the scope is the single owner of "where", which is what dissolves the old
/// multiple-owners directory scatter.
/// </summary>
/// <remarks>
/// Path-only by construction: building a scope or a tree of scopes touches no
/// disk. Directories are materialized explicitly by <see cref="EnsureDirectory"/>
/// (or by a writer at write time), so a getter never has a side effect.
/// </remarks>
public sealed class ArtifactScope
{
    private ArtifactScope(string dir) => Dir = dir;

    /// <summary>Absolute directory this scope names.</summary>
    public string Dir { get; }

    /// <summary>
    /// Root scope at <c>{baseDirectory}/{family}/{stamp}</c>. <paramref name="family"/>
    /// is the resolved run identity; <paramref name="stamp"/> the invocation timestamp
    /// (see <see cref="RunStamp"/>). Per the repository convention, callers default
    /// <paramref name="baseDirectory"/> to <c>artifacts/</c>.
    /// </summary>
    public static ArtifactScope Root(string baseDirectory, string family, string stamp)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(baseDirectory);
        ArgumentException.ThrowIfNullOrWhiteSpace(family);
        ArgumentException.ThrowIfNullOrWhiteSpace(stamp);
        return new ArtifactScope(Path.Combine(Path.GetFullPath(baseDirectory), family, stamp));
    }

    /// <summary>
    /// A role-named child scope (<c>{Dir}/{role}</c>). The child carries no stamp of
    /// its own — its position is its identity. <paramref name="role"/> is a single
    /// canonical token (<c>spc</c>, <c>csv</c>, <c>checkpoints</c>); chain for depth
    /// (<c>scope.Child("checkpoints").Child("probes")</c>).
    /// </summary>
    public ArtifactScope Child(string role)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(role);
        return new ArtifactScope(Path.Combine(Dir, role));
    }

    /// <summary>Absolute path to a canonical file directly under this scope. No side effect.</summary>
    public string File(string name)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(name);
        return Path.Combine(Dir, name);
    }

    /// <summary>Create this scope's directory (idempotent). Returns this scope for chaining.</summary>
    public ArtifactScope EnsureDirectory()
    {
        Directory.CreateDirectory(Dir);
        return this;
    }
}
