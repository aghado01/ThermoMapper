using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Microsoft.CodeAnalysis;
using Microsoft.CodeAnalysis.CSharp;
using Microsoft.CodeAnalysis.CSharp.Syntax;

namespace RepoAudit;

internal static class SyntaxWalker
{
    public static SourceFileAnalysis Analyze(string filePath)
    {
        string source;
        try { source = File.ReadAllText(filePath); }
        catch (Exception)
        {
            return new SourceFileAnalysis(filePath,
                Array.Empty<string>(), Array.Empty<string>(), Array.Empty<RawTypeInfo>(),
                Array.Empty<string>());
        }

        SyntaxTree tree = CSharpSyntaxTree.ParseText(source, path: filePath);
        CompilationUnitSyntax root = tree.GetCompilationUnitRoot();

        var importedNs = new HashSet<string>(StringComparer.Ordinal);
        foreach (UsingDirectiveSyntax u in root.DescendantNodes().OfType<UsingDirectiveSyntax>())
        {
            if (u.StaticKeyword.IsKind(SyntaxKind.StaticKeyword)) continue;
            if (u.Alias != null) continue;
            string? name = u.Name?.ToString();
            if (name != null) importedNs.Add(name);
        }

        var declaredNs = new HashSet<string>(StringComparer.Ordinal);
        foreach (BaseNamespaceDeclarationSyntax ns in root.DescendantNodes().OfType<BaseNamespaceDeclarationSyntax>())
            declaredNs.Add(ns.Name.ToString());

        var rawTypes = new List<RawTypeInfo>();
        foreach (MemberDeclarationSyntax member in root.DescendantNodes().OfType<MemberDeclarationSyntax>())
        {
            if (member.Parent is not (BaseNamespaceDeclarationSyntax or CompilationUnitSyntax)) continue;

            RawTypeInfo? info = member switch
            {
                ClassDeclarationSyntax c    => Make(c.Identifier.ValueText, TypeKind.Class,
                                                   c.Modifiers, c.BaseList, c.TypeParameterList, c, filePath),
                InterfaceDeclarationSyntax i => Make(i.Identifier.ValueText, TypeKind.Interface,
                                                   i.Modifiers, i.BaseList, i.TypeParameterList, i, filePath),
                StructDeclarationSyntax s   => Make(s.Identifier.ValueText, TypeKind.Struct,
                                                   s.Modifiers, s.BaseList, s.TypeParameterList, s, filePath),
                RecordDeclarationSyntax r   => Make(r.Identifier.ValueText,
                                                   r.ClassOrStructKeyword.IsKind(SyntaxKind.StructKeyword)
                                                       ? TypeKind.RecordStruct : TypeKind.Record,
                                                   r.Modifiers, r.BaseList, r.TypeParameterList, r, filePath),
                EnumDeclarationSyntax e     => new RawTypeInfo(e.Identifier.ValueText, TypeKind.Enum,
                                                   e.Modifiers.Select(m => m.ValueText).ToArray(),
                                                   Array.Empty<string>(), 0,
                                                   ContainingNamespace(e), filePath),
                DelegateDeclarationSyntax d => new RawTypeInfo(d.Identifier.ValueText, TypeKind.Delegate,
                                                   d.Modifiers.Select(m => m.ValueText).ToArray(),
                                                   Array.Empty<string>(),
                                                   d.TypeParameterList?.Parameters.Count ?? 0,
                                                   ContainingNamespace(d), filePath),
                _ => null
            };

            if (info != null) rawTypes.Add(info);
        }

        var referencedTypeNames = new HashSet<string>(StringComparer.Ordinal);
        CollectTypeReferences(root, referencedTypeNames);

        return new SourceFileAnalysis(
            filePath,
            declaredNs.ToArray(),
            importedNs.ToArray(),
            rawTypes.ToArray(),
            referencedTypeNames.OrderBy(n => n).ToArray());
    }

    // ── Type reference collection ─────────────────────────────────────────
    // Walks explicit type-position nodes and XML doc crefs; gathers simple
    // type names as syntactic evidence of usage.  No semantic binding needed.

    private static void CollectTypeReferences(SyntaxNode root, HashSet<string> names)
    {
        foreach (SyntaxNode node in root.DescendantNodes())
        {
            TypeSyntax? t = node switch
            {
                VariableDeclarationSyntax v       => v.Type,
                ParameterSyntax p                 => p.Type,
                MethodDeclarationSyntax m         => m.ReturnType,
                PropertyDeclarationSyntax pr      => pr.Type,
                FieldDeclarationSyntax f          => f.Declaration.Type,
                SimpleBaseTypeSyntax b            => b.Type,
                ObjectCreationExpressionSyntax o  => o.Type,
                CastExpressionSyntax c            => c.Type,
                TypeOfExpressionSyntax tof        => tof.Type,
                TypePatternSyntax tp              => tp.Type,
                _ => null
            };
            if (t != null) ExtractSimpleNames(t, names);
        }

        CollectDocCrefs(root, names);
    }

    // ── XML documentation cref references ────────────────────────────────
    // Visits every <see cref="…"/>, <exception cref="…"/>,
    // <param type="…"/> etc. in single- and multi-line doc comments and
    // extracts the simple type names referenced there.  This catches types
    // that only appear in API documentation (e.g. when callers use var for
    // local variables but document the concrete return type in remarks).

    private static void CollectDocCrefs(SyntaxNode root, HashSet<string> names)
    {
        foreach (SyntaxTrivia trivia in root.DescendantTrivia())
        {
            if (!trivia.HasStructure) continue;
            if (trivia.GetStructure() is not DocumentationCommentTriviaSyntax docComment)
                continue;

            foreach (XmlCrefAttributeSyntax crefAttr in docComment.DescendantNodes()
                         .OfType<XmlCrefAttributeSyntax>())
                ExtractCrefNames(crefAttr.Cref, names);
        }
    }

    private static void ExtractCrefNames(CrefSyntax cref, HashSet<string> names)
    {
        switch (cref)
        {
            case NameMemberCrefSyntax nm:
                // Name covers: Foo, Foo.Bar, Foo<T>, N.Foo etc.
                ExtractSimpleNames(nm.Name, names);
                // Also harvest parameter types: cref="Foo.Method(Bar b, Baz z)"
                if (nm.Parameters != null)
                    foreach (CrefParameterSyntax p in nm.Parameters.Parameters)
                        ExtractSimpleNames(p.Type, names);
                break;

            case QualifiedCrefSyntax qc:
                // Container.Member — recurse into both sides
                ExtractSimpleNames(qc.Container, names);
                ExtractCrefNames(qc.Member, names);
                break;

            case TypeCrefSyntax tc:
                ExtractSimpleNames(tc.Type, names);
                break;

            // IndexerMemberCrefSyntax, OperatorMemberCrefSyntax,
            // ConversionOperatorMemberCrefSyntax — no useful type name to extract
        }
    }

    private static void ExtractSimpleNames(TypeSyntax typeSyntax, HashSet<string> names)
    {
        switch (typeSyntax)
        {
            case PredefinedTypeSyntax:
                break; // int, string, bool, etc. — not project types

            case IdentifierNameSyntax id:
                string text = id.Identifier.ValueText;
                if (text.Length > 0 && text != "var" && text != "dynamic")
                    names.Add(text);
                break;

            case GenericNameSyntax gen:
                names.Add(gen.Identifier.ValueText);
                foreach (TypeSyntax arg in gen.TypeArgumentList.Arguments)
                    ExtractSimpleNames(arg, names);
                break;

            case QualifiedNameSyntax qual:
                // Fully-qualified: take the rightmost simple segment
                ExtractSimpleNames(qual.Right, names);
                break;

            case NullableTypeSyntax nullable:
                ExtractSimpleNames(nullable.ElementType, names);
                break;

            case ArrayTypeSyntax arr:
                ExtractSimpleNames(arr.ElementType, names);
                break;

            case TupleTypeSyntax tuple:
                foreach (TupleElementSyntax el in tuple.Elements)
                    ExtractSimpleNames(el.Type, names);
                break;
        }
    }

    private static RawTypeInfo Make(
        string name, TypeKind kind, SyntaxTokenList modifiers,
        BaseListSyntax? baseList, TypeParameterListSyntax? typeParams,
        MemberDeclarationSyntax node, string filePath) =>
        new(name, kind,
            modifiers.Select(m => m.ValueText).ToArray(),
            baseList?.Types.Select(t => t.Type.ToString()).ToArray() ?? Array.Empty<string>(),
            typeParams?.Parameters.Count ?? 0,
            ContainingNamespace(node),
            filePath);

    private static string ContainingNamespace(SyntaxNode node)
    {
        var parts = new List<string>();
        SyntaxNode? current = node.Parent;
        while (current != null)
        {
            if (current is BaseNamespaceDeclarationSyntax ns)
                parts.Insert(0, ns.Name.ToString());
            current = current.Parent;
        }
        return string.Join(".", parts);
    }
}

public static class ProjectAnalyzer
{
    public static ProjectAnalysis Analyze(ProjectOwnership ownership, CsprojGlobs globs)
    {
        var fileAnalyses = ownership.OwnedFiles
            .Where(f => f.EndsWith(".cs", StringComparison.OrdinalIgnoreCase))
            .Select(SyntaxWalker.Analyze)
            .ToList();

        var declaredNs = new HashSet<string>(StringComparer.Ordinal);
        var importedNs = new HashSet<string>(StringComparer.Ordinal);

        foreach (SourceFileAnalysis fa in fileAnalyses)
        {
            foreach (string ns in fa.DeclaredNamespaces) declaredNs.Add(ns);
            foreach (string ns in fa.ImportedNamespaces)  importedNs.Add(ns);
        }

        var fileAnalysesPublic = fileAnalyses
            .Select(fa => new FileAnalysis(fa.FilePath, fa.DeclaredNamespaces, fa.ImportedNamespaces,
                                           fa.ReferencedTypeNames))
            .ToList();

        return new ProjectAnalysis(
            ownership.CsprojPath,
            ownership.AssemblyName,
            declaredNs,
            importedNs,
            globs.ProjectReferences,
            GroupTypes(fileAnalyses),
            fileAnalysesPublic,
            globs.PackageReferences,
            globs.CompilerSettings);
    }

    private static IReadOnlyList<TypeEntry> GroupTypes(IEnumerable<SourceFileAnalysis> analyses)
    {
        var result = new List<TypeEntry>();

        foreach (var group in analyses
            .SelectMany(fa => fa.RawTypes)
            .GroupBy(t => (t.ContainingNamespace, t.Name, t.Kind)))
        {
            var members = group.ToList();
            bool isPartial = members.Count > 1 ||
                             members.Any(t => t.Modifiers.Contains("partial", StringComparer.Ordinal));

            result.Add(new TypeEntry(
                group.Key.Name,
                group.Key.Kind,
                members.SelectMany(t => t.Modifiers).Distinct(StringComparer.Ordinal)
                       .OrderBy(m => m).ToArray(),
                isPartial,
                members.Select(t => t.DeclaredInFile).Distinct().OrderBy(f => f).ToArray(),
                members.SelectMany(t => t.BaseTypes).Distinct(StringComparer.Ordinal).ToArray(),
                members[0].GenericArity,
                group.Key.ContainingNamespace));
        }

        return result.AsReadOnly();
    }
}
