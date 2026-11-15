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
/// The score's mixture model on a training split (E16, the conformal layer's slice 15b; Zhang and Candès's matching algorithm, arXiv
/// 2409.19712 §6): the simplex step meets its optimality conditions on a seeded sweep; on a two-cluster mixture the rules choose J = 2, and the
/// memberships and the clusters' distribution functions are recovered up to labels; posterior conformal prediction with the estimated
/// memberships covers exactly marginally, as a fixed function of the inputs must, and per cluster near nominal.
/// </summary>
public sealed class ScoreMixtureTests(ITestOutputHelper output)
{
    private static double Gaussian(Xoshiro256PlusPlus rng) => Normal.SampleStandard(rng);

    [Fact]
    public void TheSimplexStep_MeetsItsOptimalityConditions_OnASweep()
    {
        // 600 problems: J from 1 to 6 clusters, s from 1 to 9 levels, some profiles repeated, scales from 10⁻³ to 10³. At the returned π the
        // gradient g = Hπ − c is level on the support and no lower off it (the KKT conditions), and no simplex point does better.
        var rng = new Xoshiro256PlusPlus(1515);
        double worst = 0.0;
        for (int trial = 0; trial < 600; trial++)
        {
            int k = 1 + rng.NextInt(6), s = 1 + rng.NextInt(9);
            double scale = Math.Pow(10.0, -3.0 + 6.0 * rng.NextDouble());
            var gamma = new double[k, s];
            for (int a = 0; a < k; a++)
                for (int t = 0; t < s; t++) gamma[a, t] = a > 0 && rng.NextDouble() < 0.2 ? gamma[a - 1, t] : scale * rng.NextDouble();
            double[] tau = [.. Enumerable.Range(0, s).Select(_ => scale * (1.5 * rng.NextDouble() - 0.25))];
            var h = new double[k, k];
            var c = new double[k];
            double trace = 0.0;
            for (int a = 0; a < k; a++)
            {
                for (int b = 0; b < k; b++)
                    for (int t = 0; t < s; t++) h[a, b] += gamma[a, t] * gamma[b, t];
                for (int t = 0; t < s; t++) c[a] += gamma[a, t] * tau[t];
                trace += h[a, a];
            }
            for (int a = 0; a < k; a++) h[a, a] += 1e-12 * trace / k;
            double[] pi = [.. Enumerable.Repeat(1.0 / k, k)];
            ScoreMixture.Project(h, c, pi);
            Assert.All(pi, v => Assert.True(v >= 0.0));
            Assert.Equal(1.0, pi.Sum(), 12);
            double[] g = [.. Enumerable.Range(0, k).Select(a => Enumerable.Range(0, k).Sum(b => h[a, b] * pi[b]) - c[a])];
            double level = g.Min(), size = 1.0 + Math.Max(trace, c.Max(Math.Abs)), tolerance = 1e-8 * size;
            for (int a = 0; a < k; a++)
                if (pi[a] > 1e-9)
                {
                    worst = Math.Max(worst, (g[a] - level) / size);
                    Assert.True(g[a] - level <= tolerance, $"trial {trial}: the gradient on the support departs by {g[a] - level}");
                }
            double Objective(double[] p) => 0.5 * Enumerable.Range(0, k).Sum(a => p[a] * Enumerable.Range(0, k).Sum(b => h[a, b] * p[b])) - Enumerable.Range(0, k).Sum(a => c[a] * p[a]);
            double best = Objective(pi);
            for (int r = 0; r < 200; r++)
            {
                double[] e = [.. Enumerable.Range(0, k).Select(_ => -Math.Log(1.0 - rng.NextDouble()))];
                double sum = e.Sum();
                Assert.True(Objective([.. e.Select(v => v / sum)]) >= best - tolerance * 1e-3 * size);
            }
        }
        output.WriteLine($"the gradient's largest spread on the support, relative: {worst:E2}");
    }

    // ---- The two-cluster mixture of PosteriorConformalTests ----

    private static double Curve(double x) => Math.Sin(2.0 * Math.PI * x);
    private static double Membership(double x) => 1.0 / (1.0 + Math.Exp((x - 0.5) / 0.06));

    // x and a linear spline basis of it, knots at 0.1, …, 0.9: the logistic fits' features.
    private static double[,] Features(IReadOnlyList<double> x)
    {
        var f = new double[x.Count, 10];
        for (int i = 0; i < x.Count; i++)
        {
            f[i, 0] = x[i];
            for (int k = 1; k < 10; k++) f[i, k] = Math.Max(0.0, x[i] - 0.1 * k);
        }
        return f;
    }

    // n draws of the fixture: inputs, responses and absolute residuals from the true curve.
    private static (double[] X, double[] Y, double[] Scores) Draw(int n, Xoshiro256PlusPlus rng)
    {
        double[] x = [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())];
        double[] y = [.. x.Select(v => Curve(v) + (rng.NextDouble() < Membership(v) ? 0.2 : 1.0) * Gaussian(rng))];
        return (x, y, [.. x.Select((v, i) => Math.Abs(y[i] - Curve(v)))]);
    }

    [Fact]
    public void TwoClusters_TheMixtureIsRecovered_UpToLabels()
    {
        var rng = new Xoshiro256PlusPlus(2026_1118);
        (double[] x, _, double[] scores) = Draw(1000, rng);
        var clock = System.Diagnostics.Stopwatch.StartNew();
        ScoreMixture fit = ScoreMixture.Fit(Features(x), scores, new ScoreMixtureOptions { Seed = 3 });
        output.WriteLine($"the fit took {clock.Elapsed.TotalSeconds:F2} s: {fit.Cost}");
        output.WriteLine($"R² by J: {string.Join(", ", fit.RSquared.Select(v => v.ToString("F4")))}; J {fit.Clusters}; m {fit.Precision} "
                         + $"(mean ESS {fit.MeanEffectiveSampleSize:F1}, mean self-weight {fit.MeanSelfWeight:F4}, rule met {fit.PrecisionRuleMet}); {fit.Levels.Count} levels");
        Assert.Equal(2, fit.Clusters);
        Assert.True(fit.PrecisionRuleMet);

        // The memberships along x, matched to the truth's labels.
        double[] grid = [.. Enumerable.Range(1, 99).Select(g => g / 100.0)];
        double[,] pi = fit.Memberships(Features(grid));
        double Error(int label) => grid.Select((v, g) => Math.Abs(pi[g, label] - Membership(v))).Average();
        int narrow = Error(0) <= Error(1) ? 0 : 1;
        // Away from the meeting, within the knots: beyond the last knot the spline's slope rests on a tenth of the inputs, read on its own.
        double mean = Error(narrow), far = grid.Select((v, g) => Math.Abs(v - 0.5) > 0.15 && v >= 0.1 && v <= 0.9 ? Math.Abs(pi[g, narrow] - Membership(v)) : 0.0).Max();
        double edge = grid.Select((v, g) => v < 0.1 || v > 0.9 ? Math.Abs(pi[g, narrow] - Membership(v)) : 0.0).Max();
        // The profiles against the clusters' distribution functions at the levels, 2Φ(ξ/σ) − 1.
        double[,] gamma = fit.Profiles;
        double profile = 0.0;
        for (int t = 0; t < fit.Levels.Count; t++)
        {
            profile = Math.Max(profile, Math.Abs(gamma[narrow, t] - (2.0 * Math.Exp(Normal.StandardLogCdf(fit.Levels[t] / 0.2)) - 1.0)));
            profile = Math.Max(profile, Math.Abs(gamma[1 - narrow, t] - (2.0 * Math.Exp(Normal.StandardLogCdf(fit.Levels[t] / 1.0)) - 1.0)));
        }
        output.WriteLine($"memberships: mean error {mean:F4}, largest away from the meeting within the knots {far:F4}, beyond them {edge:F4}; profiles' largest error {profile:F4}");
        output.WriteLine("along x (x, estimate, truth): " + string.Join("; ", grid.Select((v, g) => (v, g)).Where(p => p.g % 5 == 0 || Math.Abs(p.v - 0.5) < 0.2 && p.g % 2 == 0)
                                                                          .Select(p => $"{p.v:F2} {pi[p.g, narrow]:F3} {Membership(p.v):F3}")));
        // The factorization itself recovers π to 0.01 from the true τ (TheFactorization_RecoversTheTrueMemberships); what remains is the
        // logistic fits' noise in τ at 1,000 points, largest toward the inputs' edges.
        Assert.True(mean < 0.08, $"mean membership error {mean:F4}");
        Assert.True(far < 0.2, $"membership error away from the meeting {far:F4}");
        Assert.True(edge < 0.25, $"membership error beyond the knots {edge:F4}");
        Assert.True(profile < 0.1, $"profile error {profile:F4}");
    }

    [Fact]
    public void AMiddleCluster_TheKnotsRecoverIt_FromTheRawInput_WhereLinearTauCannot()
    {
        // The narrow cluster holds the middle of x, |x − ½| below about 0.2: its membership is not monotone in x, so τ linear in x cannot
        // hold it, where the knots option's splines of x can (τ through richer features).
        static double Middle(double x) => 1.0 / (1.0 + Math.Exp((Math.Abs(x - 0.5) - 0.2) / 0.04));
        var rng = new Xoshiro256PlusPlus(2026_1119);
        double[] x = [.. Enumerable.Range(0, 1000).Select(_ => rng.NextDouble())];
        double[] scores = [.. x.Select(v => Math.Abs((rng.NextDouble() < Middle(v) ? 0.2 : 1.0) * Gaussian(rng)))];
        static double[,] Raw(IReadOnlyList<double> v)
        {
            var f = new double[v.Count, 1];
            for (int i = 0; i < v.Count; i++) f[i, 0] = v[i];
            return f;
        }
        double[] grid = [.. Enumerable.Range(1, 99).Select(g => g / 100.0)];
        double Error(ScoreMixture fit)
        {
            double[,] pi = fit.Memberships(Raw(grid));
            return Enumerable.Range(0, fit.Clusters).Min(k => grid.Select((v, g) => Math.Abs(pi[g, k] - Middle(v))).Average());
        }
        var clock = System.Diagnostics.Stopwatch.StartNew();
        ScoreMixture splines = ScoreMixture.Fit(Raw(x), scores, new ScoreMixtureOptions { Seed = 3, Knots = 6 });
        double seconds = clock.Elapsed.TotalSeconds;
        ScoreMixture linear = ScoreMixture.Fit(Raw(x), scores, new ScoreMixtureOptions { Seed = 3 });
        double splineError = Error(splines), linearError = Error(linear);
        output.WriteLine($"knots 6: J {splines.Clusters}, mean membership error {splineError:F4}, R² by J {string.Join(", ", splines.RSquared.Select(v => v.ToString("F4")))}, "
                         + $"{seconds:F2} s; linear: J {linear.Clusters}, mean membership error {linearError:F4}, R² by J {string.Join(", ", linear.RSquared.Select(v => v.ToString("F4")))}");
        Assert.Equal(2, splines.Clusters);
        Assert.True(splineError < 0.08, $"the splines' mean membership error {splineError:F4}");
        Assert.True(linearError > 0.15, $"the linear fit's mean membership error {linearError:F4}");
        Assert.Throws<ArgumentOutOfRangeException>(() => ScoreMixture.Fit(Raw(x), scores, new ScoreMixtureOptions { Knots = 0 }));
    }

    [Fact]
    public void TheFactorization_RecoversTheTrueMemberships()
    {
        // The true τ_t(x) = π₁(x)F₁(ξ_t) + π₂(x)F₂(ξ_t) at 1,000 inputs and the deciles of the mixture's scores: the factorization at J = 2
        // returns π up to labels and the clusters' distribution functions, and J = 3 gains nothing in R².
        var rng = new Xoshiro256PlusPlus(2026_1120);
        double[] x = [.. Enumerable.Range(0, 1000).Select(_ => rng.NextDouble())];
        double[] levels = [0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.65, 0.9, 1.3];
        double F(double level, double sigma) => 2.0 * Math.Exp(Normal.StandardLogCdf(level / sigma)) - 1.0;
        double[][] tau = [.. x.Select(v => levels.Select(l => Membership(v) * F(l, 0.2) + (1.0 - Membership(v)) * F(l, 1.0)).ToArray())];
        var clock = System.Diagnostics.Stopwatch.StartNew();
        (double[,] gamma, double r2) = ScoreMixture.Factor(tau, 2, 5, new Xoshiro256PlusPlus(1))!.Value;
        (_, double r3) = ScoreMixture.Factor(tau, 3, 5, new Xoshiro256PlusPlus(1))!.Value;
        output.WriteLine($"R² at J = 2: {r2:F8}, at J = 3: {r3:F8} ({clock.Elapsed.TotalSeconds:F2} s)");
        Assert.True(r2 > 1.0 - 1e-6);
        Assert.True(r3 - r2 < 1e-6);

        int narrow = gamma[0, 0] > gamma[1, 0] ? 0 : 1;
        double profile = 0.0, member = 0.0;
        for (int t = 0; t < levels.Length; t++)
        {
            profile = Math.Max(profile, Math.Abs(gamma[narrow, t] - F(levels[t], 0.2)));
            profile = Math.Max(profile, Math.Abs(gamma[1 - narrow, t] - F(levels[t], 1.0)));
        }
        var h = new double[2, 2];
        for (int a = 0; a < 2; a++)
            for (int b = 0; b < 2; b++)
                for (int t = 0; t < levels.Length; t++) h[a, b] += gamma[a, t] * gamma[b, t];
        for (int i = 0; i < x.Length; i++)
        {
            double[] c = [.. Enumerable.Range(0, 2).Select(a => Enumerable.Range(0, levels.Length).Sum(t => gamma[a, t] * tau[i][t]))];
            double[] pi = [0.5, 0.5];
            ScoreMixture.Project(h, c, pi);
            member = Math.Max(member, Math.Abs(pi[narrow] - Membership(x[i])));
        }
        output.WriteLine($"the profiles' largest error {profile:E2}, the memberships' {member:E2}");
        Assert.True(profile < 0.01, $"profile error {profile:E2}");
        Assert.True(member < 0.01, $"membership error {member:E2}");
    }

    [Fact]
    public void TheExtrapolatedAlternation_EndsNoHigher_InFewerPasses()
    {
        // The true two-cluster τ at 1,000 inputs, exact and with noise 0.01 clamped into (0, 1): at J 2 to 4, five restarts each, the
        // extrapolated alternation's residual is not above the plain one's and it takes no more passes.
        foreach (double noise in new[] { 0.0, 0.01 })
        {
            var rng = new Xoshiro256PlusPlus(2026_1160 + 1000);
            double[] x = [.. Enumerable.Range(0, 1000).Select(_ => rng.NextDouble())];
            double[] levels = [0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.65, 0.9, 1.3];
            double F(double level, double sigma) => 2.0 * Math.Exp(Normal.StandardLogCdf(level / sigma)) - 1.0;
            double[][] tau = [.. x.Select(v => levels.Select(l => Math.Clamp(Membership(v) * F(l, 0.2) + (1.0 - Membership(v)) * F(l, 1.0)
                                                                                + noise * Gaussian(rng), 1e-6, 1 - 1e-6)).ToArray())];
            foreach (int j in new[] { 2, 3, 4 })
            {
                var plain = new ScoreMixtureCost();
                var pushed = new ScoreMixtureCost();
                (_, double plainR2) = ScoreMixture.Factor(tau, j, 5, new Xoshiro256PlusPlus(7), plain, extrapolate: false)!.Value;
                (_, double pushedR2) = ScoreMixture.Factor(tau, j, 5, new Xoshiro256PlusPlus(7), pushed)!.Value;
                output.WriteLine($"noise {noise}, J {j}: plain 1 − R² {1 - plainR2:E6} in {plain.Passes} passes; extrapolated {1 - pushedR2:E6} in {pushed.Passes} "
                                 + $"({pushed.PushesKept} of {pushed.Pushes} pushes kept)");
                Assert.True(1 - pushedR2 <= (1 - plainR2) * (1 + 1e-4), $"J {j}, noise {noise}");
                Assert.True(pushed.Passes <= plain.Passes, $"J {j}, noise {noise}: {pushed.Passes} passes against {plain.Passes}");
            }
        }
    }

    [Fact]
    public void EndToEnd_TheEstimatedMembershipsCover()
    {
        // The mixture fitted once on a training split of 400, the calibration split's size, as their D′ and D are: the precision rule reads
        // its effective sample size on the training split. Then 2,000 replications of 399 calibration units and a test case, the randomized
        // weights at the estimated memberships beside the unweighted split.
        const int n = 399, replications = 2000;
        const double alpha = 0.1;
        var rng = new Xoshiro256PlusPlus(2026_1119);
        (double[] tx, _, double[] ts) = Draw(400, rng);
        ScoreMixture fit = ScoreMixture.Fit(Features(tx), ts, new ScoreMixtureOptions { Seed = 5 });
        var covered = new int[2, 2];
        var widths = new double[2, 2];
        var finite = new int[2, 2];
        var tests = new int[2];
        for (int r = 0; r < replications; r++)
        {
            (double[] x, double[] y, _) = Draw(n, rng);
            double xt = rng.NextDouble(), yt = Curve(xt) + (rng.NextDouble() < Membership(xt) ? 0.2 : 1.0) * Gaussian(rng);
            int region = Membership(xt) >= 0.5 ? 0 : 1;
            tests[region]++;
            var centre = new ResidualReadings([.. x.Select(Curve)]);
            var at = new ResidualReadings([Curve(xt)]);
            ConformalInterval[] read =
            [
                SplitConformal.Residuals(centre, y, at, alpha).Intervals[0],
                SplitConformal.Residuals(centre, y, at, alpha, fit.Weights(Features(x), Features([xt]), seed: r).Weights).Intervals[0],
            ];
            for (int f = 0; f < 2; f++)
            {
                if (read[f].Contains(yt)) covered[f, region]++;
                if (!read[f].WholeSpace) { widths[f, region] += read[f].Width; finite[f, region]++; }
            }
        }
        double nominal = 1 - alpha, se = Math.Sqrt(alpha * nominal / replications);
        string[] rules = ["unweighted", "estimated memberships"];
        for (int f = 0; f < 2; f++)
            output.WriteLine($"{rules[f]}: marginal {(covered[f, 0] + covered[f, 1]) / (double)replications:F4}; cluster 1 {covered[f, 0] / (double)tests[0]:F4}, "
                             + $"cluster 2 {covered[f, 1] / (double)tests[1]:F4}; mean width {widths[f, 0] / finite[f, 0]:F3} and {widths[f, 1] / finite[f, 1]:F3}, "
                             + $"whole space {tests[0] - finite[f, 0]} and {tests[1] - finite[f, 1]}");
        output.WriteLine($"J {fit.Clusters}, m {fit.Precision}; {tests[0]} and {tests[1]} tests; se {se:F4}");
        double marginal = (covered[1, 0] + covered[1, 1]) / (double)replications;
        Assert.True(marginal >= nominal - 3.5 * se, $"marginal {marginal:F4}");
        double wide = covered[1, 1] / (double)tests[1], plain = covered[0, 1] / (double)tests[1], seWide = Math.Sqrt(alpha * nominal / tests[1]);
        Assert.True(wide > plain + 3.0 * seWide, $"the wide cluster: {wide:F4} against the unweighted {plain:F4}");
    }

    [Fact]
    public void UnequalSplits_TheRuleReadsTheCalibrationSize()
    {
        // Slice 15's finding: fitted on 1,000 and calibrated on 399, the training's own reading chose m near its ceiling and gave a share of the
        // sets the whole space where the clusters meet. The rule at the calibration's size takes a smaller m. Both m's on the same 2,000
        // replications of 399 calibration units and a test case.
        const int n = 399, replications = 2000;
        const double alpha = 0.1;
        var rng = new Xoshiro256PlusPlus(2026_1137);
        (double[] tx, _, double[] ts) = Draw(1000, rng);
        ScoreMixture fit = ScoreMixture.Fit(Features(tx), ts, new ScoreMixtureOptions { Seed = 5 });
        PrecisionChoice atCalibration = fit.PrecisionAt(n);
        int[] precisions = [fit.Precision, atCalibration.Precision];
        var covered = new int[2, 2];
        var whole = new int[2];
        var widths = new double[2];
        var tests = new int[2];
        for (int r = 0; r < replications; r++)
        {
            (double[] x, double[] y, _) = Draw(n, rng);
            double xt = rng.NextDouble(), yt = Curve(xt) + (rng.NextDouble() < Membership(xt) ? 0.2 : 1.0) * Gaussian(rng);
            int region = Membership(xt) >= 0.5 ? 0 : 1;
            tests[region]++;
            var centre = new ResidualReadings([.. x.Select(Curve)]);
            var at = new ResidualReadings([Curve(xt)]);
            double[,] units = fit.Memberships(Features(x)), cases = fit.Memberships(Features([xt]));
            for (int f = 0; f < 2; f++)
            {
                ShiftWeights weights = PosteriorConformal.Weights(units, cases, precisions[f], seed: r).Weights;
                ConformalInterval read = SplitConformal.Residuals(centre, y, at, alpha, weights).Intervals[0];
                if (read.Contains(yt)) covered[f, region]++;
                if (read.WholeSpace) whole[f]++; else widths[f] += read.Width;
            }
        }
        double nominal = 1 - alpha, se = Math.Sqrt(alpha * nominal / replications);
        string[] rules = ["the training's own reading", "the calibration's size"];
        for (int f = 0; f < 2; f++)
            output.WriteLine($"{rules[f]}: m {precisions[f]}; marginal {(covered[f, 0] + covered[f, 1]) / (double)replications:F4}, cluster 1 "
                             + $"{covered[f, 0] / (double)tests[0]:F4}, cluster 2 {covered[f, 1] / (double)tests[1]:F4}; whole space {whole[f] / (double)replications:P1}; "
                             + $"mean finite width {widths[f] / (replications - whole[f]):F3}");
        output.WriteLine($"at {n} units: mean ESS {atCalibration.MeanEffectiveSampleSize:F1}, mean self-weight {atCalibration.MeanSelfWeight:F4}, rule met {atCalibration.RuleMet}; "
                         + $"the training's: {fit.MeanEffectiveSampleSize:F1} and {fit.MeanSelfWeight:F4}; se {se:F4}");
        Assert.True(atCalibration.Precision < fit.Precision);
        Assert.True(atCalibration.RuleMet);
        Assert.True(whole[1] < whole[0], "the calibration's size should leave fewer whole-space sets");
        Assert.True(whole[1] <= 0.01 * replications, $"whole space {whole[1]} of {replications}");
        double marginal = (covered[1, 0] + covered[1, 1]) / (double)replications;
        Assert.True(marginal >= nominal - 3.5 * se, $"marginal {marginal:F4}");
        foreach (int region in new[] { 0, 1 })
            Assert.True(covered[1, region] / (double)tests[region] >= nominal - 3.5 * Math.Sqrt(alpha * nominal / tests[region]), $"cluster {region + 1}");
    }

    [Fact]
    public void TheRescaledReading_MatchesPoolsOfTheCalibrationsSize_AndEqualSplitsReadAsBefore()
    {
        var rng = new Xoshiro256PlusPlus(2026_1138);
        (double[] tx, _, double[] ts) = Draw(1000, rng);
        ScoreMixture fit = ScoreMixture.Fit(Features(tx), ts, new ScoreMixtureOptions { Seed = 5 });
        double[,] pi = fit.Memberships(Features(tx));
        double[][] rows = [.. Enumerable.Range(0, tx.Length).Select(i => Enumerable.Range(0, fit.Clusters).Select(k => pi[i, k]).ToArray())];

        // The rescaled sums against pools of 399 others drawn at random, each point weighed as the reading weighs it. The rescaling reads the
        // expected sums, and a point's self-weight is convex in its others' sum, so it reads the self-weight a little low where the weights
        // concentrate: closely at the m the rule takes, within 5% at the ceiling.
        foreach ((int m, double tolerance) in new[] { (fit.PrecisionAt(399).Precision, 0.01), (fit.Precision, 0.05) })
        {
            (double ess, double self) = ScoreMixture.Reading(rows, m, 5, 399.0 / 999.0);
            var pick = new Xoshiro256PlusPlus(11);
            double poolEss = 0.0, poolSelf = 0.0;
            for (int i = 0; i < rows.Length; i++)
            {
                int[] draw = PosteriorConformal.Draw(rows[i], m, PosteriorConformal.Stream(5, i));
                int[] others = [.. Enumerable.Range(0, rows.Length).Where(j => j != i)];
                for (int a = 0; a < 399; a++)
                {
                    int b = a + pick.NextInt(others.Length - a);
                    (others[a], others[b]) = (others[b], others[a]);
                }
                double LogWeight(int j) => Enumerable.Range(0, fit.Clusters).Where(c => draw[c] > 0).Sum(c => draw[c] * Math.Log(rows[j][c]));
                double[] logs = [LogWeight(i), .. others.Take(399).Select(LogWeight)];
                double top = logs.Max();
                double[] w = [.. logs.Select(l => Math.Exp(l - top))];
                poolEss += w.Sum() * w.Sum() / w.Sum(v => v * v);
                poolSelf += w[0] / w.Sum();
            }
            poolEss /= rows.Length;
            poolSelf /= rows.Length;
            output.WriteLine($"m {m}: rescaled ESS {ess:F2} and self-weight {self:F5}; pools of 399 drawn: {poolEss:F2} and {poolSelf:F5}");
            Assert.True(Math.Abs(ess / poolEss - 1.0) < tolerance, $"ESS {ess:F2} against {poolEss:F2}");
            Assert.True(Math.Abs(self / poolSelf - 1.0) < tolerance, $"self-weight {self:F5} against {poolSelf:F5}");
            Assert.True(self <= poolSelf, "the rescaled self-weight should read low, by convexity");
        }

        // Equal splits, n − 1 calibration units, read as before, bit for bit; a smaller split takes a smaller m.
        PrecisionChoice own = fit.PrecisionAt(999);
        Assert.Equal((fit.Precision, fit.MeanEffectiveSampleSize, fit.MeanSelfWeight, fit.PrecisionRuleMet), (own.Precision, own.MeanEffectiveSampleSize, own.MeanSelfWeight, own.RuleMet));
        (double[] cx, _, _) = Draw(999, rng);
        PosteriorWeights viaFit = fit.Weights(Features(cx), Features([0.3, 0.7]), seed: 2);
        PosteriorWeights direct = PosteriorConformal.Weights(fit.Memberships(Features(cx)), fit.Memberships(Features([0.3, 0.7])), fit.Precision, seed: 2);
        Assert.Equal(direct.SelfWeights, viaFit.SelfWeights);
        Assert.Equal(direct.EffectiveSampleSizes, viaFit.EffectiveSampleSizes);
        PosteriorWeights smaller = fit.Weights(Features(cx[..199]), Features([0.3, 0.7]), seed: 2);
        PosteriorWeights smallerDirect = PosteriorConformal.Weights(fit.Memberships(Features(cx[..199])), fit.Memberships(Features([0.3, 0.7])),
                                                                    fit.PrecisionAt(199).Precision, seed: 2);
        Assert.Equal(smallerDirect.SelfWeights, smaller.SelfWeights);
        Assert.True(fit.PrecisionAt(199).Precision <= fit.PrecisionAt(399).Precision);
        Assert.True(fit.PrecisionAt(399).Precision <= fit.Precision);
        Assert.Same(fit.PrecisionAt(399), fit.PrecisionAt(399));
        Assert.Throws<ArgumentOutOfRangeException>(() => fit.PrecisionAt(0));

        // A fixed m stands at every size, its readings and the rule's verdict read there.
        ScoreMixture fixedFit = ScoreMixture.Fit(Features(tx), ts, new ScoreMixtureOptions { Seed = 5, Clusters = 2, Precision = 400 });
        PrecisionChoice fixedSmall = fixedFit.PrecisionAt(99);
        Assert.Equal(400, fixedSmall.Precision);
        Assert.False(fixedSmall.RuleMet);
        Assert.True(fixedSmall.MeanSelfWeight > fixedFit.MeanSelfWeight);
    }

    [Fact]
    public void WhatTheMixtureDoesNotServe_IsRefused()
    {
        var rng = new Xoshiro256PlusPlus(17);
        (double[] x, _, double[] scores) = Draw(200, rng);
        double[,] features = Features(x);
        Assert.Throws<ArgumentException>(() => ScoreMixture.Fit(features, scores[..199]));
        Assert.Throws<ArgumentException>(() => ScoreMixture.Fit(features, [.. scores.Select((v, i) => i == 7 ? double.NaN : v)]));
        Assert.Contains("one value", Assert.Throws<ArgumentException>(() => ScoreMixture.Fit(features, [.. scores.Select(_ => 0.3)])).Message);
        foreach (ScoreMixtureOptions bad in new[]
        {
            new ScoreMixtureOptions { Levels = 0 },
            new ScoreMixtureOptions { Clusters = 0 },
            new ScoreMixtureOptions { MaxClusters = 0 },
            new ScoreMixtureOptions { ClusterGain = -0.1 },
            new ScoreMixtureOptions { Precision = 0 },
            new ScoreMixtureOptions { MinPrecision = 0 },
            new ScoreMixtureOptions { MinPrecision = 50, MaxPrecision = 20 },
            new ScoreMixtureOptions { MaxSelfWeight = 0.0 },
            new ScoreMixtureOptions { Restarts = 0 },
        })
            Assert.Throws<ArgumentOutOfRangeException>(() => ScoreMixture.Fit(features, scores, bad));

        // The fit is the seed's.
        ScoreMixture one = ScoreMixture.Fit(features, scores, new ScoreMixtureOptions { Seed = 4 }), two = ScoreMixture.Fit(features, scores, new ScoreMixtureOptions { Seed = 4 });
        Assert.Equal(one.Profiles, two.Profiles);
        Assert.Equal(one.Precision, two.Precision);
        Assert.Throws<ArgumentException>(() => ScoreMixture.Fit(features, scores, new ScoreMixtureOptions { Clusters = 201 }));

        // J and m fixed: the rules stand aside.
        ScoreMixture fixedFit = ScoreMixture.Fit(features, scores, new ScoreMixtureOptions { Clusters = 3, Precision = 17 });
        Assert.Equal(3, fixedFit.Clusters);
        Assert.Equal(17, fixedFit.Precision);
        Assert.Single(fixedFit.RSquared);
        Assert.Throws<ArgumentException>(() => fixedFit.Memberships(new double[3, 4]));
        double[,] pi = fixedFit.Memberships(Features([0.1, 0.5, 0.9]));
        for (int i = 0; i < 3; i++) Assert.Equal(1.0, pi[i, 0] + pi[i, 1] + pi[i, 2], 12);
    }
}
