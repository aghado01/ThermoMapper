using System;
using StorageAudit;
using Xunit;

namespace StorageAudit.Tests;

public class PathGlobTests
{
    [Fact]
    public void StarStaysWithinOneSegment()
    {
        var glob = new PathGlob("artifacts/*/");
        Assert.True(glob.IsMatch("artifacts/bin", isDirectory: true));
        Assert.True(glob.IsMatch("artifacts/bin/", isDirectory: true));
        Assert.False(glob.IsMatch("artifacts/bin/Debug", isDirectory: true));
        Assert.False(glob.IsMatch("artifacts", isDirectory: true));
    }

    [Fact]
    public void DoubleStarSpansSegments_IncludingNone()
    {
        var glob = new PathGlob("**/.snapshot/*/");
        Assert.True(glob.IsMatch(".snapshot/20260827", isDirectory: true));
        Assert.True(glob.IsMatch("src/.snapshot/20260827", isDirectory: true));
        Assert.True(glob.IsMatch("a/b/c/.snapshot/x", isDirectory: true));
        Assert.False(glob.IsMatch("src/.snapshot", isDirectory: true));
    }

    [Fact]
    public void TrailingSlashMatchesFoldersOnly_OtherwiseFilesOnly()
    {
        var folder = new PathGlob("artifacts/bin/");
        Assert.True(folder.IsMatch("artifacts/bin", isDirectory: true));
        Assert.False(folder.IsMatch("artifacts/bin", isDirectory: false));

        var files = new PathGlob("artifacts/*.tsv");
        Assert.True(files.IsMatch("artifacts/parity.tsv", isDirectory: false));
        Assert.False(files.IsMatch("artifacts/parity.tsv", isDirectory: true));
        Assert.False(files.IsMatch("artifacts/sub/parity.tsv", isDirectory: false));
    }

    [Fact]
    public void QuestionMarkIsOneCharacter_AndMatchingIgnoresCase()
    {
        var glob = new PathGlob("chip-??/");
        Assert.True(glob.IsMatch("CHIP-FH", isDirectory: true));
        Assert.False(glob.IsMatch("chip-fh3", isDirectory: true));
    }

    [Fact]
    public void RegexCharactersAreLiteral_AndBackslashesAreSeparators()
    {
        var glob = new PathGlob(@"a+b\(c).txt");
        Assert.True(glob.IsMatch("a+b/(c).txt", isDirectory: false));
        Assert.False(glob.IsMatch("aab/(c)xtxt", isDirectory: false));
    }

    [Fact]
    public void EmptyOrRootGlobsAreRejected()
    {
        Assert.Throws<ArgumentException>(() => new PathGlob(""));
        Assert.Throws<ArgumentException>(() => new PathGlob("/"));
    }
}
