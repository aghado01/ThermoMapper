using System.Linq;
using System.Xml.Linq;
using Xunit;

namespace TestHarness.Runner.Tests;

/// <summary>
/// Per-case selection: the runsettings filter for a display name, the suite a case is filed under, and the
/// display names the harness cannot tell apart.
/// </summary>
public sealed class CaseSelectionTests
{
    [Theory]
    [InlineData("Ns.C.M", @"DisplayName=Ns.C.M")]
    [InlineData("Ns.C.M(shape: \"sine\")", @"DisplayName=Ns.C.M\(shape: ""sine""\)")]
    [InlineData("Ns.C.M(s: \"a\\b;c|d&e!f~g=h\")", @"DisplayName=Ns.C.M\(s: ""a\\b;c\|d\&e\!f\~g\=h""\)")]
    public void RunSettings_FilterEscapesTheDisplayName_AndTreatsNoMatchAsError(string displayName, string filter)
    {
        XElement config = FactRunner.CaseRunSettings(displayName).Element("RunConfiguration")!;
        Assert.Equal(filter, config.Element("TestCaseFilter")!.Value);
        Assert.Equal("true", config.Element("TreatNoTestsAsError")!.Value);
    }

    [Theory]
    [InlineData("Ns.Sub.Suite.Method", "Suite")]
    [InlineData("Ns.Suite.Method(x: 0.5)", "Suite")]
    [InlineData("Ns.Suite.Method(s: \"a.b.c\")", "Suite")]
    public void SuiteName_ComesFromTheMethodPath_NotTheArguments(string displayName, string suite)
    {
        Assert.Equal(suite, PathUtil.SuiteNameFromFact(displayName));
    }

    [Fact]
    public void CaseCollisions_GroupNamesThatDifferOnlyInLetterCase()
    {
        string[] facts =
        [
            "Ns.C.M(s: \"abc\")",
            "Ns.C.M(s: \"ABC\")",
            "Ns.C.M(s: \"xyz\")",
            "Ns.C.Other",
            "Ns.C.other",
        ];

        var groups = TestDiscovery.CaseCollisions(facts);

        Assert.Equal(2, groups.Count);
        Assert.Contains(groups, g => g.Count == 2 && g.Contains("Ns.C.M(s: \"abc\")") && g.Contains("Ns.C.M(s: \"ABC\")"));
        Assert.Contains(groups, g => g.Count == 2 && g.Contains("Ns.C.Other") && g.Contains("Ns.C.other"));
    }

    [Fact]
    public void CaseCollisions_NoneAmongDistinctNames()
    {
        Assert.Empty(TestDiscovery.CaseCollisions(["Ns.C.A", "Ns.C.B", "Ns.C.M(s: \"a\")"]));
    }
}
