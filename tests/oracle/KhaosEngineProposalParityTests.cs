using System;
using System.Linq;
using System.Text.Json;
using Maths.Regression.Spline.Bass;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// The engine's KHAOS proposal (the KHAOS brief's slice 1) against <c>khaos</c> 2.1.1 at slice 0's frozen states
/// (<c>khaos_conditionals_oracle.R</c>): <see cref="LegendreTermProposal.InclusionProbabilities"/> is <c>make_weights</c> as
/// written, and the inclusion's law is the exact forward probability, Σ_{q₀} J(q₀)·P(χ | q₀)/P(valid | q₀), that slice 0 showed
/// <c>khaos</c>'s own draw follows; <c>khaos</c>'s ratios use other forms, which slice 0 confirmed as defects.
/// </summary>
public sealed class KhaosEngineProposalParityTests(ITestOutputHelper output)
{
    [Fact]
    public void Inclusion_IsKhaosMakeWeights_AndItsExactForwardLaw()
    {
        if (KhaosOracle.Conditionals.Value is not { } r) return;
        int checkedCases = 0;
        foreach (JsonElement c in r.GetProperty("birth").EnumerateArray())
        {
            int order = c.GetProperty("order").GetInt32(), degree = c.GetProperty("degree").GetInt32();
            int[] eta = [.. GbassOracle.Numbers(c.GetProperty("eta")).Select(v => (int)v)];
            double[][] W = [.. c.GetProperty("weights").EnumerateArray().Select(GbassOracle.Numbers)];
            int p = eta.Length;
            var prior = new BassStructurePrior(new LegendreTermSpace(p, degree, order), new BassCountPrior(4, 0.5));
            var proposal = new LegendreTermProposal(prior, new LegendreBasis { MaxDegree = degree });
            double worstWeight = 0.0, worstLaw = 0.0;
            var engine = new double[p];
            for (int q0 = 1; q0 <= order; q0++)
            {
                proposal.InclusionProbabilities(eta, q0, engine);
                for (int j = 0; j < p; j++) worstWeight = Math.Max(worstWeight, Math.Abs(engine[j] / W[q0 - 1][j] - 1.0));
            }
            foreach (JsonElement e in c.GetProperty("chis").EnumerateArray())
            {
                bool[] chi = [.. GbassOracle.Numbers(e.GetProperty("chi")).Select(v => v != 0.0)];
                double exact = e.GetProperty("exact").GetDouble();
                worstLaw = Math.Max(worstLaw, Math.Abs(Math.Exp(proposal.LogInclusion(eta, chi)) / exact - 1.0));
            }
            output.WriteLine($"birth at counts {string.Join(",", eta)}, q_max {order}: make_weights' coins within {worstWeight:E2}, the forward law within {worstLaw:E2} (relative)");
            Assert.True(worstWeight < 1e-13, $"make_weights differs by {worstWeight:E2}");
            Assert.True(worstLaw < 1e-12, $"the forward law differs by {worstLaw:E2}");
            checkedCases++;
        }
        Assert.Equal(KhaosOracle.ProposalCases.Count(c => c[0] == 1), checkedCases);
    }
}
