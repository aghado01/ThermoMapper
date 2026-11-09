using System;
using System.Linq;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers;
using Maths.Samplers.Ensemble;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// E21's ridge blocks beyond KHAOS's first parts (<c>issues/regression/bass/ridge-blocks-brief.md</c>), slice by slice.
/// <list type="bullet">
/// <item><b>Slice 0, the whitening on the working response:</b> a held-s block's fits keep their draws bit for bit, pinned by hash on a
/// hinge chain, a Legendre chain and a tempered hinge fit with its reference rung; the hinge prefix tables stand aside under a block, whose
/// columns they would read unwhitened; a design follows a reading at another scale as it follows moved weights, to the fit of a design
/// built there.</item>
/// </list>
/// </summary>
public sealed class RandomEffectsTests(ITestOutputHelper output)
{
    // A fit's hash: the bit-identity tests' (structures, τ, σ², coefficients, the move counts, the ladder), with each draw's group effects.
    private static string Hash(BassFit fit)
    {
        var hash = new RunFingerprint().Add(BassBitIdentityTests.Hash(fit));
        foreach (DrawStore<BassDraw> chain in fit.Draws)
            foreach (BassDraw d in chain) hash.Add(d.GroupEffects ?? []);
        return hash.ToString();
    }

    private static (BassConfig Config, BassData Data) Pinned(string fit)
    {
        var sampling = new BassSampling { Chains = 2, MasterSeed = 2108, BurnIn = 300, Samples = 200 };
        if (fit == "tempered hinge")
        {
            var spec = new TemperingSpec { Levels = 3, HotEnd = 0.2, ReferenceRung = true };
            (BassModel model, BassConfig config) = RandomInterceptTests.Fixture("hinge", spec);
            return (config with { Sampling = sampling with { Tempering = spec } }, model.Data);
        }
        (BassModel m, BassConfig c) = RandomInterceptTests.Fixture(fit);
        return (c with { Sampling = sampling }, m.Data);
    }

    // Slice 0's gate: moving the whitening onto the working response leaves a held-s block's fits as they were.
    [Theory]
    [InlineData("hinge", "b9b95aee57ef459c")]
    [InlineData("legendre", "a1f4528886afa5f7")]
    [InlineData("tempered hinge", "05122582e8f56f75")]
    public void AHeldScaleFit_KeepsItsDraws(string fit, string expected)
    {
        (BassConfig config, BassData data) = Pinned(fit);
        string hash = Hash(Bass.Fit(config, data));
        output.WriteLine($"{fit}: {hash}");
        Assert.Equal(expected, hash);
    }

    // The design follows a reading at another scale (WorkingResponse.Whitening) as it follows moved weights: every column whitened afresh and
    // the Gram factored strictly, so its fit is a design's built on that reading and loaded with the same terms. Back on its first reading
    // it holds its first fit.
    [Fact]
    public void TheDesign_FollowsTheWhitening_ToAnotherScale()
    {
        (BassModel model, _) = RandomInterceptTests.Fixture("hinge");
        WorkingResponse At(double s) => RandomInterceptBlock.Of(model.Data.Grouping!, new BassRandomIntercept(s)).Whitened(new WorkingResponse(model.Data.Response));
        BassDesign moving = model.NewDesign();
        var rng = new Xoshiro256PlusPlus(seed: 41);
        for (int k = 0; k < 200 && moving.Count < 4; k++)
        {
            var counts = new BassTermCounts(model.Prior);
            foreach (BassTerm t in moving.Terms) counts.Add(t);
            if (model.Proposal.DrawBirth(rng, counts) is { } term && moving.ProposeBirth(term).Admitted) moving.Accept();
        }
        Assert.True(moving.Count >= 2, $"the design holds {moving.Count} terms");
        double rss0 = moving.ResidualSquares, fit0 = moving.FittedSquares;
        double[] beta0 = moving.Coefficients.ToArray();

        WorkingResponse other = At(0.4);
        Assert.True(moving.ProposeFollow(other).Admitted);
        moving.Accept();
        Assert.Same(other.Whitening, moving.Working.Whitening);
        BassDesign fresh = model.NewDesign(other);
        foreach (BassTerm t in moving.Terms)
        {
            Assert.True(fresh.ProposeBirth(t).Admitted);
            fresh.Accept();
        }
        double Rel(double a, double b) => Math.Abs(a - b) / Math.Max(Math.Abs(a), Math.Abs(b));
        double worst = Math.Max(Rel(moving.ResidualSquares, fresh.ResidualSquares), Rel(moving.FittedSquares, fresh.FittedSquares));
        double scale = Math.Max(1.0, fresh.Coefficients.ToArray().Max(Math.Abs));
        for (int j = 0; j < fresh.Count; j++) worst = Math.Max(worst, Math.Abs(moving.Coefficients[j] - fresh.Coefficients[j]) / scale);
        Assert.NotEqual(rss0, moving.ResidualSquares);   // the scale moved the fit
        output.WriteLine($"{moving.Count} terms; at s = 0.4 against a design built there, the worst relative gap {worst:E2}; ‖r‖² {rss0:G6} at 1.3, {moving.ResidualSquares:G6} at 0.4");
        Assert.True(worst <= 1e-12, $"{worst:E2}");

        Assert.True(moving.ProposeFollow(model.Working).Admitted);
        moving.Accept();
        double back = Math.Max(Rel(moving.ResidualSquares, rss0), Rel(moving.FittedSquares, fit0));
        for (int j = 0; j < beta0.Length; j++) back = Math.Max(back, Math.Abs(moving.Coefficients[j] - beta0[j]) / scale);
        Assert.True(back <= 1e-12, $"back at s = 1.3: {back:E2}");
        Assert.Throws<ArgumentException>(() => moving.ProposeFollow(new WorkingResponse(model.Data.Response)));   // no whitening: another layout
    }

    // Under a block the design holds whitened columns, V^(−½)xⱼ, and a birth's border is their inner product with the new column's
    // whitened image. The hinge prefix tables sum a held column against the raw hinge, which would give xⱼᵀV^(−½)h for xⱼᵀV⁻¹h, so they
    // stand aside: with the tables asked for, every proposal is the dots' own.
    [Fact]
    public void UnderABlock_TheTablesStandAside()
    {
        var point = BassCost.Base with { N = 600, P = 6, Cap = 15 };
        BassData plain = BassCost.Data(point);
        var data = new BassData(plain.Predictors, plain.Response, RandomInterceptTests.Grouping(point.N, 8, seed: 11));
        BassModel model = BassModel.Create(BassCost.Config(point) with { RandomIntercept = new BassRandomIntercept(1.3) }, data);
        BassDesign off = model.NewDesign(), on = model.NewDesign();
        off.TablesOn = false;
        on.TablesOn = true;
        var rng = new Xoshiro256PlusPlus(seed: 5151);
        int proposals = 0;
        double worst = 0.0;
        for (int step = 0; step < 1500; step++)
        {
            var counts = new BassTermCounts(model.Prior);
            foreach (BassTerm t in off.Terms) counts.Add(t);
            DesignProposal a, b;
            if ((rng.NextDouble() < 0.6 || off.Count == 0) && off.Count < model.MaxBasis)
            {
                if (model.Proposal.DrawBirth(rng, counts) is not { } t) continue;
                a = off.ProposeBirth(t);
                b = on.ProposeBirth(t);
            }
            else if (off.Count > 0)
            {
                int i = (int)(rng.NextUInt64() % (ulong)off.Count);
                a = off.ProposeDeath(i);
                b = on.ProposeDeath(i);
            }
            else continue;
            proposals++;
            Assert.Equal(a.Refusal, b.Refusal);
            if (!a.Admitted) continue;
            worst = Math.Max(worst, Math.Abs(a.ResidualSquares - b.ResidualSquares) / a.ResidualSquares);
            if (rng.NextDouble() < 0.5)
            {
                off.Accept();
                on.Accept();
            }
        }
        output.WriteLine($"{proposals} proposals; the tables read {on.TableLookups} entries; the worst relative ‖r‖² between tables and dots {worst:E2}");
        Assert.Equal(0L, on.TableLookups);
        Assert.Equal(0.0, worst);
    }
}
