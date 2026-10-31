using System;
using System.IO;
using System.Linq;
using System.Text.Json;
using StorageAudit;
using Xunit;
using static StorageAudit.Tests.Fixture;

namespace StorageAudit.Tests;

public class ArchiveTests
{
    [Fact]
    public void ArchiveCopiesVerifiesAndWritesAManifestBeforeRemovingTheSource()
    {
        using var f = new Fixture();
        f.Write("artifacts/spc/run/a.json", 1500, (byte)'a');
        f.Write("artifacts/spc/b.json", 700, (byte)'b');
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/spc/", Category.Archivable)));
        string archive = Path.Combine(f.Scratch, "thumb");

        ApplyResult r = Assert.Single(f.Apply(execute: true, archiveTo: archive));

        string destination = Path.Combine(archive, "repo", "artifacts", "spc");
        Assert.Equal("archived", r.Outcome);
        Assert.Equal(destination, r.Destination);
        Assert.False(Directory.Exists(f.Full("artifacts/spc")));
        Assert.Equal(1500, new FileInfo(Path.Combine(destination, "run", "a.json")).Length);
        Assert.Equal(700, new FileInfo(Path.Combine(destination, "b.json")).Length);

        string manifestPath = Assert.Single(Directory.GetFiles(archive, "storage-audit-manifest-*.json"));
        ArchiveManifest manifest = JsonSerializer.Deserialize<ArchiveManifest>(File.ReadAllText(manifestPath), StorageRules.JsonOptions)!;
        ManifestItem item = Assert.Single(manifest.Items);
        Assert.Equal(2, item.Files);
        Assert.Equal(2200, item.Bytes);
        Assert.Empty(ArchiveVerifier.Verify(item.FileList, destination));
    }

    [Fact]
    public void TheVerifierCatchesAMissingFileAndAnAlteredOne()
    {
        using var f = new Fixture();
        string dest = Path.Combine(f.Scratch, "copy");
        string a = f.Write(Path.Combine(dest, "a.bin"), 100, (byte)'a');
        string b = f.Write(Path.Combine(dest, "b.bin"), 100, (byte)'b');
        ManifestFile[] manifest =
        [
            new("a.bin", 100, ArchiveVerifier.HashFile(a)),
            new("b.bin", 100, ArchiveVerifier.HashFile(b)),
        ];
        Assert.Empty(ArchiveVerifier.Verify(manifest, dest));

        f.Write(b, 100, (byte)'c');                     // same size, different bytes
        File.Delete(a);

        var problems = ArchiveVerifier.Verify(manifest, dest);
        Assert.Contains(problems, p => p.StartsWith("missing a.bin", StringComparison.Ordinal));
        Assert.Contains(problems, p => p.StartsWith("hash of b.bin", StringComparison.Ordinal));
    }

    [Fact]
    public void ArchiveRefusesAnExistingDestination_AnArchiveInsideTheRepo_AndExecutingWithoutAnArchive()
    {
        using var f = new Fixture();
        f.Write("artifacts/spc/a.json");
        f.Age(10);
        f.Plan(Rules(Rule("artifacts/spc/", Category.Archivable)));

        string thumb = Path.Combine(f.Scratch, "thumb");
        f.Write(Path.Combine(thumb, "repo", "artifacts", "spc", "older.json"));
        Assert.Contains("already holds", Refused(f.Apply(execute: true, archiveTo: thumb)));
        Assert.Contains("inside the root 'repo'", Refused(f.Apply(execute: true, archiveTo: f.Full("archive"))));
        Assert.Contains("needs --archive-to", Refused(f.Apply(execute: true)));

        Assert.Equal("would-archive", Assert.Single(f.Apply()).Outcome);   // a dry run may leave it out
        Assert.True(File.Exists(f.Full("artifacts/spc/a.json")));
    }
}
