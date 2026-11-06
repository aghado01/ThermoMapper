using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.Json;
using Maths.Distributions;
using Maths.Rng;
using Maths.Samplers.Blocks;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Oracle.Tests;

/// <summary>
/// GBASS's full conditionals at a frozen state (<c>gbass_conditionals_oracle.R</c>; E1's build brief, slice 0). The oracle evaluates
/// <c>gbass</c>'s and <c>nwbass</c>'s own expressions as written and holds each against the joint GBASS states, by the spread over eight
/// points of the joint less the conditional's log kernel. Here each is held against the engine's forms:
/// <list type="bullet">
/// <item><b>v under t mixing</b> is <see cref="StudentTWeights"/>' draw: λ = 1/v ~ Gamma((ν + 1)/2, (ν + r²/(wc))/2), from
/// <see cref="ConjugateUpdates.NormalPrecision(double, double, double, double)"/>, number for number;</item>
/// <item><b>τ</b> is <see cref="ConjugateUpdates.NormalVariance(double, double, double, double)"/> over the M + 1 coefficients, and
/// <b>w at β = 0</b> the GIG form of the same update (N-X3) over the N residuals and the M + 1 coefficients;</item>
/// <item><b>β and γ</b> are the normals of their closed forms; γ's is its conditional on γ &gt; 0 alone, since the Wald's normalizer is
/// e^|γ|, and its mass below 0 is the size of the untruncated draw's error (E1's defect);</item>
/// <item><b>w at β ≠ 0</b> is the modified half-normal on w^(−½) when w's prior has a = 0, and misses the joint when a &gt; 0, the term
/// GBASS drops (E1's defect);</item>
/// <item><b>a and a birth's likelihood</b> as written equal dense algebra on the stated joint.</item>
/// </list>
/// </summary>
public sealed class GbassConditionalParityTests(ITestOutputHelper output)
{
    private const int N = 40;
    private static readonly double[] Coefficients = [0.9, 2.1, -1.1];

    // A two-hinge design with the intercept inside the prior, as GBASS holds it, and latent scales and a response drawn once.
    private static (double[] Y, double[] V, double[][] B) State(int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        var y = new double[N]; var v = new double[N]; var b = new double[N][];
        for (int i = 0; i < N; i++)
        {
            double x = rng.NextDouble();
            b[i] = [1.0, Math.Max(x - 0.3, 0) / 0.7, Math.Max(0.6 - x, 0) / 0.6];
            v[i] = Gamma.ShapeRate(1.0, 1.0).Sample(rng);
            y[i] = b[i].Zip(Coefficients, (p, q) => p * q).Sum() + 0.3 * Math.Sqrt(v[i]) * Normal.Standard.Sample(rng);
        }
        return (y, v, b);
    }

    private sealed record Run(JsonElement R, double[] Y, double[] V, double[][] B, Dictionary<string, double> P);

    private static Run? Evaluate(int seed, Dictionary<string, double> p)
    {
        (double[] y, double[] v, double[][] b) = State(seed);
        string[] tail = [$"a={string.Join(";", Coefficients.Select(c => c.ToString("R", CultureInfo.InvariantCulture)))}",
                         .. p.Select(kv => $"{kv.Key}={kv.Value.ToString("R", CultureInfo.InvariantCulture)}")];
        JsonElement? r = GbassOracle.Run("gbass_conditionals_oracle.R", Enumerable.Range(0, N).Select(i => (double[])[y[i], v[i], .. b[i]]), tail);
        return r is { } j ? new Run(j, y, v, b, p) : null;
    }

    private static Dictionary<string, double> Keys(double beta, double wA, double vP, double vA, double vB) => new()
    {
        ["w"] = 0.08, ["beta"] = beta, ["tau"] = 3, ["gamma"] = 0.8, ["scale"] = 1.3, ["w_p"] = -0.1, ["w_a"] = wA, ["w_b"] = 0.1,
        ["v_p"] = vP, ["v_a"] = vA, ["v_b"] = vB, ["a_tau"] = 0.5, ["b_tau"] = 20, ["m_beta"] = 0.2, ["s_beta"] = 2, ["m_gamma"] = 1,
        ["s_gamma"] = 3, ["lower_bound"] = 0.01,
    };

    // A skewed state with general GIG mixing and w's prior free of a; and a t(5) state at β = 0 with a > 0 in w's prior.
    private static readonly Lazy<Run?> Skewed = new(() => Evaluate(20261105, Keys(0.7, 0.0, -1.5, 0.6, 3.0)));
    private static readonly Lazy<Run?> Symmetric = new(() => Evaluate(20261106, Keys(0.0, 0.4, -2.5, 0.0, 5.0)));

    private const double Rounding = 1e-9;   // nats over eight points of a log joint near 10²

    [Fact]
    public void AsWritten_TheConditionalsAreTheJoints_SaveTheTwoDefects()
    {
        if (Skewed.Value is not { } s || Symmetric.Value is not { } t) return;
        foreach ((string name, Run run) in new[] { ("skewed", s), ("symmetric", t) })
        {
            JsonElement r = run.R;
            double Spread(string block, string field = "spread") => GbassOracle.Number(r.GetProperty(block).GetProperty(field));
            output.WriteLine($"{name}: spreads v {Spread("v"):E2}, beta {Spread("beta"):E2}, tau {Spread("tau"):E2}, w_gig {Spread("w_gig"):E2}, " +
                             $"w_mhn {Spread("w_mhn"):E2}, gamma+ {Spread("gamma", "spread_positive"):E2}, gamma- {Spread("gamma", "spread_negative"):E2}");
            Assert.True(Spread("v") < Rounding);
            Assert.True(Spread("beta") < Rounding);
            Assert.True(Spread("tau") < Rounding);
            Assert.True(Spread("w_gig") < Rounding);
            Assert.True(Spread("gamma", "spread_positive") < Rounding);
            // γ's normal is no conditional below 0: the joint's normalizer there is e^{−γN}, not e^{γN}.
            Assert.True(Spread("gamma", "spread_negative") > 1.0);
            JsonElement a = r.GetProperty("a");
            Assert.Equal(GbassOracle.Numbers(a.GetProperty("mean_dense")), GbassOracle.Numbers(a.GetProperty("mean_code")), (x, y) => Math.Abs(x - y) < 1e-10);
            Assert.True(a.GetProperty("cov_max_diff").GetDouble() < 1e-12 * a.GetProperty("cov_scale").GetDouble());
            JsonElement birth = r.GetProperty("birth");
            Assert.Equal(birth.GetProperty("dense").GetDouble(), birth.GetProperty("code").GetDouble(), 9);
        }
        // w at β ≠ 0: exact while w's prior has a = 0 (the skewed state), off the joint once it has a > 0 (the symmetric state's).
        Assert.True(GbassOracle.Number(s.R.GetProperty("w_mhn").GetProperty("spread")) < Rounding);
        Assert.True(GbassOracle.Number(t.R.GetProperty("w_mhn").GetProperty("spread")) > 1e-3);
    }

    [Fact]
    public void TMixing_VIsStudentTWeightsDraw_NumberForNumber()
    {
        if (Symmetric.Value is not { } t) return;
        double w = t.P["w"], c = t.P["scale"], nu = t.P["v_b"];
        Assert.Equal(-0.5 * nu, t.P["v_p"]);
        double[] r = Residuals(t);
        Assert.Equal(r, GbassOracle.Numbers(t.R.GetProperty("r")), (x, y) => Math.Abs(x - y) < 1e-13);
        JsonElement v = t.R.GetProperty("v");
        double p = v.GetProperty("p").GetDouble();
        Assert.Equal(0.0, v.GetProperty("a").GetDouble());
        double[] b = GbassOracle.Numbers(v.GetProperty("b"));
        for (int i = 0; i < N; i++)
        {
            // StudentTWeights' step at σ² = wc: the precision λ = 1/v ~ Gamma((ν + 1)/2, (ν + r²/σ²)/2).
            Gamma lambda = ConjugateUpdates.NormalPrecision(0.5 * nu, 0.5 * nu, 1.0, r[i] * r[i] / (w * c));
            Assert.Equal(-p, lambda.Shape);
            Assert.Equal(1.0, 0.5 * b[i] / lambda.Rate, 14);
        }
    }

    [Fact]
    public void TauAndW_AreTheConjugateUpdates()
    {
        if (Skewed.Value is not { } s || Symmetric.Value is not { } t) return;
        foreach (Run run in new[] { s, t })
        {
            JsonElement r = run.R;
            double w = run.P["w"], tau = run.P["tau"], c = run.P["scale"];
            int m = Coefficients.Length;
            double pen = Coefficients.Sum(a => a * a);
            Assert.Equal(pen, r.GetProperty("pen").GetDouble(), 13);
            InverseGamma tauPost = ConjugateUpdates.NormalVariance(run.P["a_tau"], run.P["b_tau"], m, pen / w);
            Assert.Equal(tauPost.Shape, r.GetProperty("tau").GetProperty("shape").GetDouble());
            Assert.Equal(1.0, r.GetProperty("tau").GetProperty("rate").GetDouble() / tauPost.Scale, 14);
            double[] res = Residuals(run);
            double rssOverV = res.Select((e, i) => e * e / run.V[i]).Sum() / c;
            Assert.Equal(1.0, r.GetProperty("rss_over_v").GetDouble() / rssOverV, 13);
            GigParameters wPost = ConjugateUpdates.NormalVariance(new GigParameters(run.P["w_p"], run.P["w_a"], run.P["w_b"]), N + m, rssOverV + pen / tau);
            JsonElement g = r.GetProperty("w_gig");
            Assert.Equal(wPost.P, g.GetProperty("p").GetDouble(), 13);
            Assert.Equal(wPost.A, g.GetProperty("a").GetDouble());
            Assert.Equal(1.0, g.GetProperty("b").GetDouble() / wPost.B, 13);
            // The modified half-normal on w^(−½): α = N + (M + 1) − 2p_w, β = ½·b's term, γ = βΣr/c, w's prior a dropped.
            JsonElement h = r.GetProperty("w_mhn");
            Assert.Equal(N + m - 2 * run.P["w_p"], h.GetProperty("alpha").GetDouble(), 12);
            Assert.Equal(1.0, h.GetProperty("beta").GetDouble() / (0.5 * (run.P["w_b"] + rssOverV + pen / tau)), 13);
            Assert.Equal(run.P["beta"] * res.Sum() / c, h.GetProperty("gamma").GetDouble(), 12);
        }
    }

    [Fact]
    public void BetaAndGamma_AreTheirNormals_AndGammasDrawCrossesZero()
    {
        if (Skewed.Value is not { } s) return;
        double w = s.P["w"], c = s.P["scale"], sb = s.P["s_beta"], mb = s.P["m_beta"], sg = s.P["s_gamma"], mg = s.P["m_gamma"];
        double[] res = Residuals(s);
        double sumV = s.V.Sum();
        double precision = sumV / c + 1.0 / (sb * sb);
        JsonElement beta = s.R.GetProperty("beta");
        Assert.Equal(1.0, beta.GetProperty("mean").GetDouble() / ((res.Sum() / (c * Math.Sqrt(w)) + mb / (sb * sb)) / precision), 12);
        Assert.Equal(1.0, beta.GetProperty("sd").GetDouble() * Math.Sqrt(precision), 12);
        JsonElement gamma = s.R.GetProperty("gamma");
        double gm = (sg * sg * N + mg) / (sg * sg * sumV + 1), gs = sg / Math.Sqrt(sg * sg * sumV + 1);
        Assert.Equal(1.0, gamma.GetProperty("mean").GetDouble() / gm, 12);
        Assert.Equal(1.0, gamma.GetProperty("sd").GetDouble() / gs, 12);
        Assert.Equal(Normal.MeanSd(gm, gs).Cdf(0.0), gamma.GetProperty("mass_below_zero").GetDouble(), 12);
    }

    private static double[] Residuals(Run run) => [.. run.Y.Select((yi, i) => yi - run.B[i].Zip(Coefficients, (p, q) => p * q).Sum())];
}
