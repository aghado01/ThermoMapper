using System.Collections.Generic;

namespace RepoAudit;

public enum TypeKind { Class, Interface, Struct, Record, RecordStruct, Enum, Delegate }

public sealed record TypeEntry(
    string Name,
    TypeKind Kind,
    IReadOnlyList<string> Modifiers,
    bool IsPartial,
    IReadOnlyList<string> ContributingFiles,
    IReadOnlyList<string> BaseTypes,
    int GenericArity,
    string ContainingNamespace);

public sealed record FileAnalysis(
    string FilePath,
    IReadOnlyList<string> DeclaredNamespaces,
    IReadOnlyList<string> ImportedNamespaces,
    IReadOnlyList<string> ReferencedTypeNames);

public sealed record ProjectAnalysis(
    string CsprojPath,
    string AssemblyName,
    IReadOnlySet<string> DeclaredNamespaces,
    IReadOnlySet<string> ImportedNamespaces,
    IReadOnlyList<string> ProjectReferences,
    IReadOnlyList<TypeEntry> Types,
    IReadOnlyList<FileAnalysis> FileAnalyses,
    IReadOnlyList<PackageRef> PackageReferences,
    CsprojCompilerSettings CompilerSettings);

internal sealed record SourceFileAnalysis(
    string FilePath,
    IReadOnlyList<string> DeclaredNamespaces,
    IReadOnlyList<string> ImportedNamespaces,
    IReadOnlyList<RawTypeInfo> RawTypes,
    IReadOnlyList<string> ReferencedTypeNames);

internal sealed record RawTypeInfo(
    string Name,
    TypeKind Kind,
    IReadOnlyList<string> Modifiers,
    IReadOnlyList<string> BaseTypes,
    int GenericArity,
    string ContainingNamespace,
    string DeclaredInFile);
