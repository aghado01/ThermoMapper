using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.Conformal;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The split form with held-out groups as its units (E16's follow-up for new groups, agreed with the GBASS chip for E21's slice 5): singleton
/// groups are the case-level split form bit for bit; on a two-level model the worst-case group score holds a new group's cases all at once
/// and the one-case score holds a single case, each at 1 − α, where calibrating on cases holds the group jointly far less often.
/// </summary>
public sealed class GroupConformalTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    // A one-draw seam: each case's log density under N(0, sd²) at the values of its grid.
    private static CandidateLogDensity Seam(double[][] grids, double sd, int[]? strata = null)
    {
        int columns = grids.Sum(g => g.Length);
        var values = new double[1, columns];
        int c = 0;
        foreach (double[] grid in grids)
            foreach (double v in grid) values[0, c++] = Normal.MeanSd(0.0, sd).LogPdf(v);
        return new CandidateLogDensity(values, grids, [], strata);
    }

    [Fact]
    public void SingletonGroups_AreTheCaseLevelSplitForm()
    {
        var rng = new Xoshiro256PlusPlus(2201);
        double[] y = [.. Enumerable.Range(0, 39).Select(_ => 1.2 * Gaussian(rng))];
        double[] grid = [.. Enumerable.Range(0, 41).Select(g => -4.0 + 0.2 * g)];
        CandidateLogDensity calibration = Seam([.. y.Select(v => new[] { v })], 1.1);
        CandidateLogDensity candidates = Seam([grid, grid], 1.1);
        ConformalPrediction cases = SplitConformal.Density(calibration, candidates);
        int[] own = [.. Enumerable.Range(0, y.Length).Select(i => 100 - 2 * i)];
        foreach (GroupScore score in Enum.GetValues<GroupScore>())
        {
            ConformalPrediction groups = GroupConformal.Density(calibration, own, candidates, score: score);
            for (int t = 0; t < 2; t++) Assert.Equal(cases.Sets[t].PValues, groups.Sets[t].PValues);
        }
    }

    [Fact]
    public void NewGroups_AreCoveredByGroupUnits()
    {
        // y = u_g + ε, u_g ~ N(0, 1), ε ~ N(0, ½²); the fit the true new-group predictive N(0, 1.25). 79 held-out calibration groups of 3 to 8
        // cases, and a new group drawn alike, its size included; 2,000 replications, α 0.1.
        const int groups = 79, replications = 2000;
        var rng = new Xoshiro256PlusPlus(2026_1136);
        double sd = Math.Sqrt(1.25);
        int worstJoint = 0, oneCase = 0, casesJoint = 0, casesSingle = 0;
        for (int r = 0; r < replications; r++)
        {
            var y = new List<double>();
            var label = new List<int>();
            for (int g = 0; g < groups; g++)
            {
                double u = Gaussian(rng);
                int size = 3 + rng.NextInt(6);
                for (int i = 0; i < size; i++) { y.Add(u + 0.5 * Gaussian(rng)); label.Add(g); }
            }
            double un = Gaussian(rng);
            int testCases = 3 + rng.NextInt(6);
            double[] test = [.. Enumerable.Range(0, testCases).Select(_ => un + 0.5 * Gaussian(rng))];
            CandidateLogDensity calibration = Seam([.. y.Select(v => new[] { v })], sd);
            CandidateLogDensity candidates = Seam([.. test.Select(v => new[] { v })], sd);
            ConformalPrediction worst = GroupConformal.Density(calibration, label, candidates, score: GroupScore.WorstCase);
            ConformalPrediction one = GroupConformal.Density(calibration, label, candidates, score: GroupScore.OneCase, seed: r);
            ConformalPrediction byCase = SplitConformal.Density(calibration, candidates);
            if (worst.Sets.All(s => s.Accepts(0))) worstJoint++;
            if (one.Sets[0].Accepts(0)) oneCase++;
            if (byCase.Sets.All(s => s.Accepts(0))) casesJoint++;
            if (byCase.Sets[0].Accepts(0)) casesSingle++;
        }
        double se = Math.Sqrt(0.09 / replications);
        double Rate(int count) => count / (double)replications;
        output.WriteLine($"a new group of 3 to 8: worst-case group units hold it whole {Rate(worstJoint):F4}; one-case group units hold one case {Rate(oneCase):F4}; "
                         + $"case units hold it whole {Rate(casesJoint):F4} and one case {Rate(casesSingle):F4}; se {se:F4}");
        Assert.True(Rate(worstJoint) >= 0.9 - 3.5 * se);
        Assert.True(Rate(oneCase) >= 0.9 - 3.5 * se);
        Assert.True(Rate(casesJoint) < Rate(worstJoint) - 3.0 * se, "case units should hold a new group whole less often");
    }

    [Fact]
    public void OneCase_IsSeeded_AndWhatGroupsDoNotServe_IsRefused()
    {
        var rng = new Xoshiro256PlusPlus(2202);
        double[] y = [.. Enumerable.Range(0, 30).Select(_ => Gaussian(rng))];
        int[] label = [.. Enumerable.Range(0, 30).Select(i => i / 3)];
        double[] grid = [-1.0, 0.0, 1.0];
        CandidateLogDensity calibration = Seam([.. y.Select(v => new[] { v })], 1.0);
        CandidateLogDensity candidates = Seam([grid], 1.0);
        ConformalPrediction a = GroupConformal.Density(calibration, label, candidates, score: GroupScore.OneCase, seed: 4);
        ConformalPrediction b = GroupConformal.Density(calibration, label, candidates, score: GroupScore.OneCase, seed: 4);
        Assert.Equal(a.Sets[0].PValues, b.Sets[0].PValues);
        Assert.Equal(11, a.Sets[0].Units);

        Assert.Throws<ArgumentException>(() => GroupConformal.Density(calibration, label[..29], candidates));
        Assert.Throws<ArgumentOutOfRangeException>(() => GroupConformal.Density(calibration, label, candidates, score: (GroupScore)3));
        Assert.Throws<ArgumentException>(() => GroupConformal.Density(calibration, label, Seam([grid], 1.0, [1])));
        Assert.Throws<ArgumentException>(() => GroupConformal.Density(calibration with { NewGroup = [.. y.Select((_, i) => i != 5)] }, label, candidates));
        Assert.Throws<ArgumentException>(() => GroupConformal.Density(calibration, label, candidates with { NewGroup = [false] }));
        // All new groups: served.
        ConformalPrediction flagged = GroupConformal.Density(calibration with { NewGroup = [.. y.Select(_ => true)] }, label, candidates with { NewGroup = [true] });
        Assert.Equal(GroupConformal.Density(calibration, label, candidates).Sets[0].PValues, flagged.Sets[0].PValues);
    }
}
