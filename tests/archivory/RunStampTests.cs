using System;
using System.IO;
using Archivory;
using Xunit;

namespace Archivory.Tests;

// Pure facts — no I/O, no shared state — so every one is safe under the parallel
// fact harness with no fixture or isolation.
public sealed class RunStampTests
{
    [Fact]
    public void For_FormatsIsoBasicMillisecondUtc()
    {
        var instant = new DateTime(2026, 5, 26, 0, 45, 14, 123, DateTimeKind.Utc);
        Assert.Equal("20260526T004514.123Z", RunStamp.For(instant));
    }

    [Fact]
    public void For_PadsMillisecondsToThreeDigits()
    {
        var instant = new DateTime(2026, 1, 2, 3, 4, 5, 5, DateTimeKind.Utc);
        Assert.Equal("20260102T030405.005Z", RunStamp.For(instant));
    }

    [Fact]
    public void For_ConvertsNonUtcInstantToUtc()
    {
        // Derive a local instant from a known UTC so the assertion is host-TZ-independent:
        // formatting the local form must yield the same stamp as the UTC it came from.
        var utc = new DateTime(2026, 5, 26, 0, 45, 14, 123, DateTimeKind.Utc);
        DateTime asLocal = utc.ToLocalTime();
        Assert.Equal(RunStamp.For(utc), RunStamp.For(asLocal));
    }

    [Fact]
    public void Now_MatchesTheCanonicalShape()
        => Assert.Matches(@"^\d{8}T\d{6}\.\d{3}Z$", RunStamp.Now());

    [Fact]
    public void Now_IsAValidDirectoryName()
        => Assert.Equal(-1, RunStamp.Now().IndexOfAny(Path.GetInvalidFileNameChars()));
}
