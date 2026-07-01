using System;
using Archivory;
using Xunit;

namespace Archivory.Tests;

// Pure resolution truth table — no I/O, no shared state — parallel-safe per fact.
public sealed class RunIdentityTests
{
    [Fact]
    public void ExplicitName_WinsAndIsRecorded()
    {
        var id = RunIdentity.Resolve("my_study", callerStub: "spc");
        Assert.Equal("my_study", id.Family);
        Assert.Equal("explicit", id.Source);
        Assert.Equal("my_study", id.Requested);
    }

    [Theory]
    [InlineData(null)]
    [InlineData("")]
    [InlineData("   ")]
    public void NoRequest_FallsBackToCallerStub(string? requested)
    {
        var id = RunIdentity.Resolve(requested, callerStub: "spc");
        Assert.Equal("spc", id.Family);
        Assert.Equal("auto:caller=spc", id.Source);
        Assert.Null(id.Requested);
    }

    [Fact]
    public void ExplicitName_IsSanitizedToOnePathSegment()
    {
        // Space and '/' are separators on every platform; collapse to a single '_'.
        var id = RunIdentity.Resolve("my run/v2", callerStub: "spc");
        Assert.Equal("my_run_v2", id.Family);
        Assert.DoesNotContain("__", id.Family);
        Assert.Equal("explicit", id.Source);
    }

    [Fact]
    public void BlankCallerStub_Throws()
        => Assert.Throws<ArgumentException>(() => RunIdentity.Resolve(null, callerStub: "  "));
}
