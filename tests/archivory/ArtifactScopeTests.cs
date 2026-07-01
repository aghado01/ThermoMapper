using System;
using System.IO;
using Archivory;
using Xunit;

namespace Archivory.Tests;

public sealed class ArtifactScopeTests
{
    // ── Pure path composition (no I/O → fully parallel-safe per fact) ──────────

    [Fact]
    public void Root_ComposesBaseFamilyStamp()
    {
        var scope = ArtifactScope.Root("artifacts", "spc_n", "20260526T004514.123Z");
        Assert.True(Path.IsPathRooted(scope.Dir));
        Assert.EndsWith(Path.Combine("spc_n", "20260526T004514.123Z"), scope.Dir);
    }

    [Fact]
    public void Child_AppendsRoleUnderParent()
    {
        var root = ArtifactScope.Root("artifacts", "spc_n", "stamp");
        Assert.Equal(Path.Combine(root.Dir, "spc"), root.Child("spc").Dir);
    }

    [Fact]
    public void Child_ChainsForDepth()
    {
        var root = ArtifactScope.Root("artifacts", "spc_n", "stamp");
        ArtifactScope probes = root.Child("checkpoints").Child("probes");
        Assert.Equal(Path.Combine(root.Dir, "checkpoints", "probes"), probes.Dir);
    }

    [Fact]
    public void File_ResolvesDirectlyUnderScope()
    {
        var scope = ArtifactScope.Root("artifacts", "spc_n", "stamp").Child("spc");
        Assert.Equal(Path.Combine(scope.Dir, "manifest.json"), scope.File("manifest.json"));
    }

    [Fact]
    public void BuildingAScope_TouchesNoDisk()
    {
        // A unique base that is never created; pure path composition must not materialize it.
        string ghostBase = Path.Combine(Path.GetTempPath(), "archivory-ghost-" + Guid.NewGuid().ToString("N"));
        var scope = ArtifactScope.Root(ghostBase, "fam", "stamp").Child("spc");
        _ = scope.File("x.json");
        Assert.False(Directory.Exists(scope.Dir));
    }

    [Theory]
    [InlineData("")]
    [InlineData("   ")]
    public void Root_RejectsBlankSegments(string blank)
    {
        Assert.Throws<ArgumentException>(() => ArtifactScope.Root(blank, "fam", "stamp"));
        Assert.Throws<ArgumentException>(() => ArtifactScope.Root("base", blank, "stamp"));
        Assert.Throws<ArgumentException>(() => ArtifactScope.Root("base", "fam", blank));
    }

    // ── The one I/O fact: isolated under a per-fact unique temp dir so concurrent
    //    harness workers never race the same path. ──────────────────────────────

    [Fact]
    public void EnsureDirectory_CreatesTheDirectory()
    {
        string baseDir = Path.Combine(Path.GetTempPath(), "archivory-tests-" + Guid.NewGuid().ToString("N"));
        try
        {
            ArtifactScope scope = ArtifactScope.Root(baseDir, "fam", "stamp").Child("csv").EnsureDirectory();
            Assert.True(Directory.Exists(scope.Dir));
        }
        finally
        {
            if (Directory.Exists(baseDir))
                Directory.Delete(baseDir, recursive: true);
        }
    }
}
