using System;
using Xunit;

namespace TestHarness.Runner.Tests;

/// <summary>
/// The harness selects each case itself, through a per-case runsettings file; a pass-through argument that
/// replaces or filters that selection would make a case run the wrong tests, or none, and still report.
/// </summary>
public sealed class HarnessOptionsTests
{
    [Theory]
    [InlineData("--filter")]
    [InlineData("--filter:Category=Slow")]
    [InlineData("--settings")]
    [InlineData("--settings:other.runsettings")]
    [InlineData("-s")]
    [InlineData("RunConfiguration.TestCaseFilter=Name~X")]
    [InlineData("-p:VSTestTestCaseFilter=Name~X")]
    [InlineData("/p:VSTestSetting=other.runsettings")]
    public void CommonArg_ThatOverridesCaseSelection_IsRejected(string arg)
    {
        Assert.Throws<ArgumentException>(() =>
            HarnessOptions.Parse(["--project", "p.csproj", "--common-arg", arg]));
    }

    [Theory]
    [InlineData("--blame")]
    [InlineData("--logger:trx")]
    [InlineData("-p:ParallelizeTestCollections=false")]
    public void CommonArg_ThatLeavesCaseSelectionAlone_IsAccepted(string arg)
    {
        HarnessOptions opts = HarnessOptions.Parse(["--project", "p.csproj", "--common-arg", arg]);
        Assert.Equal(new[] { arg }, opts.CommonArgs);
    }
}
