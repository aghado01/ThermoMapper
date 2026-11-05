using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Maths.Samplers.Ensemble;
using Maths.Samplers.Rjmcmc;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// The efficiency pass's lever 1 (<c>efficiency-pass-brief.md</c>): a one-way hinge's cross products with the held columns from prefix
/// tables. The tables are a pure cache, so a proposal's numbers do not depend on which tables a design has built; they differ from the
/// dots' by rounding, not in bits. The checks:
/// <list type="bullet">
/// <item>the scan against the explicit dot, to 4·K·ε·‖z‖‖u‖ (K the knots), over sizes, distributions, signs, degrees and held columns of
/// every kind, unit and precision weights;</item>
/// <item>a design's proposals with tables on against the same design with them off: admission equal, the account within the tolerance,
/// the committed design within it, over scripted births, deaths and changes;</item>
/// <item>the chain's law: the enumerated fixtures' chains with tables on against exact π, and paired seeds, tables on and off, on a rich
/// surface and a weighted one and on a tempered ladder;</item>
/// <item>purity: a chain restored from a snapshot, with no tables, continues bit for bit.</item>
/// </list>
/// </summary>
public sealed class BassPrefixTableTests
{
    private readonly ITestOutputHelper _out;
    public BassPrefixTableTests(ITestOutputHelper output) => _out = output;

    // ---- The scan ----

    private static double[] DrawX(Xoshiro256PlusPlus rng, int n, string kind) => kind switch
    {
        "mass" => [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble() < 0.99 ? 0.5 : rng.NextDouble())],
        "ties" => [.. Enumerable.Range(0, n).Select(_ => Math.Floor(rng.NextDouble() * 10.0) / 10.0)],
        "ends" => [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble() < 0.5 ? 0.01 * rng.NextDouble() : 1.0 - 0.01 * rng.NextDouble())],
        "integers" => [.. Enumerable.Range(0, n).Select(_ => (double)(rng.NextUInt64() % (ulong)Math.Max(2, n)))],
        _ => [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble())],
    };

    private static double[] DrawColumn(Xoshiro256PlusPlus rng, double[] x, string kind)
    {
        int n = x.Length;
        double[] z = kind switch
        {
            "heavy" => [.. Enumerable.Range(0, n).Select(_ => BassQualification.Gaussian(rng) / Math.Max(0.02, rng.NextDouble()))],
            "collinear" => [.. x.Select(v => v + 1e-6 * BassQualification.Gaussian(rng))],
            _ => [.. Enumerable.Range(0, n).Select(_ => BassQualification.Gaussian(rng))],
        };
        double mean = PairwiseSum.Of(z) / n;
        for (int i = 0; i < n; i++) z[i] -= mean;
        return z;
    }

    [Fact]
    public void Scan_AgainstTheDot_WithinTheStatedBound()
    {
        var rng = new Xoshiro256PlusPlus(seed: 9100);
        double worst = 0.0;
        long checkedSums = 0;
        foreach (int n in new[] { 3, 4, 5, 10, 100, 1000, 5000, 20000 })
            foreach (string xKind in new[] { "uniform", "mass", "ties", "ends", "integers" })
            {
                double[] x = DrawX(rng, n, xKind);
                double[] knots = [.. x.Distinct().Order()];
                HingeOrder order = HingeOrder.Of(x, knots) ?? throw new InvalidOperationException("knots are the distinct values");
                int k = knots.Length;
                foreach (string zKind in new[] { "gaussian", "heavy", "collinear" })
                    foreach (bool weighted in new[] { false, true })
                    {
                        double[] z = DrawColumn(rng, x, zKind);
                        double[]? roots = weighted ? [.. Enumerable.Range(0, n).Select(_ => 0.5 + 1.5 * rng.NextDouble())] : null;
                        double[] g = roots is null ? z : [.. z.Zip(roots, (a, b) => a * b)];
                        double normZ = Math.Sqrt(g.Sum(v => v * v));
                        foreach (int degree in new[] { 1, 2 })
                        {
                            double[] above = HingeScan.Above(z, roots, order, degree), below = HingeScan.Below(z, roots, order, degree);
                            int[] picks = [.. new[] { 0, 1, k / 2, k - 2, k - 1 }.Concat(Enumerable.Range(0, 4).Select(_ => (int)(rng.NextUInt64() % (ulong)k))).Where(p => p >= 0 && p < k).Distinct()];
                            foreach (int knot in picks)
                                foreach (int sign in new[] { 1, -1 })
                                {
                                    double t = knots[knot];
                                    double[] u = [.. x.Select(v => Math.Pow(Math.Max(sign * (v - t), 0.0), degree))];
                                    double dot = PairwiseSum.Dot(g, u), normU = Math.Sqrt(u.Sum(v => v * v));
                                    double table = sign > 0 ? above[knot] : below[knot];
                                    double bound = (4.0 * k + 8.0 * Math.Log2(n + 1)) * 1.1102230246251565e-16 * normZ * normU;
                                    double error = Math.Abs(table - dot);
                                    Assert.True(error <= bound + 1e-300, $"n={n} {xKind} {zKind} weighted={weighted} q={degree} sign={sign} knot={knot}/{k}: |{table:R} − {dot:R}| = {error:E2} > {bound:E2}");
                                    if (bound > 0) worst = Math.Max(worst, error / (normZ * normU) / 1.1102230246251565e-16);
                                    checkedSums++;
                                }
                        }
                    }
            }
        _out.WriteLine($"{checkedSums} sums within the bound; the worst error was {worst:F1} ε of ‖z‖‖u‖");
    }

    [Fact]
    public void Order_RefusesKnotsThatAreNotTheDistinctValues()
    {
        Assert.Null(HingeOrder.Of([0.0, 0.5, 1.0], [0.0, 1.0]));
        Assert.Null(HingeOrder.Of([0.0, 0.5, 1.0], [0.0, 0.25, 0.5, 1.0]));
        Assert.NotNull(HingeOrder.Of([0.0, 0.5, 0.5, 1.0], [0.0, 0.5, 1.0]));
    }

    // ---- A design's proposals, tables on against off ----

    private static BassModel Model(int n, int p, int cap, bool weighted, out BassData data, bool ridge = false, bool categorical = false,
                                   WeightedPriorGram gram = WeightedPriorGram.Weighted)
    {
        var point = categorical ? new BassCost.Point(BassCost.Surface.Categorical, n, 7, cap, BassCost.BaseOrder, L: 8) : BassCost.Base with { N = n, P = p, Cap = cap };
        data = BassCost.Data(point);
        BassConfig config = BassCost.Config(point);
        if (ridge) config = config with { Coefficients = BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma()) };
        if (!weighted) return BassModel.Create(config, data);
        double[] w = BassWeightedCoreTests.Weights(n, seed: 21);
        return BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, weights: w, priorGram: gram)));
    }

    private static double Rel(double a, double b) => a.Equals(b) ? 0.0 : Math.Abs(a - b) / Math.Max(1e-300, Math.Max(Math.Abs(a), Math.Abs(b)));   // equal infinities and NaNs agree

    [Theory]
    [InlineData(false, false, false)]
    [InlineData(true, false, false)]
    [InlineData(false, true, false)]   // the normalized ridge
    [InlineData(false, false, true)]   // categorical terms among the held columns
    [InlineData(true, false, true)]
    public void Design_ProposalsFromTables_AgreeWithTheDots(bool weighted, bool ridge, bool categorical)
    {
        BassModel model = Model(600, 6, 15, weighted, out _, ridge, categorical);
        BassDesign off = model.NewDesign(), on = model.NewDesign();
        off.TablesOn = false;
        on.TablesOn = true;
        var rng = new Xoshiro256PlusPlus(seed: 5150 + (weighted ? 1 : 0));
        int proposals = 0, admitted = 0, accepted = 0;
        double worstRss = 0, worstFit = 0, worstAngle = 0;
        for (int step = 0; step < 4000; step++)
        {
            var counts = new BassTermCounts(model.Prior);
            foreach (BassTerm t in off.Terms) counts.Add(t);
            int m = off.Count;
            double kind = rng.NextDouble();
            DesignProposal a, b;
            if (kind < 0.4 && m < model.MaxBasis)
            {
                if (model.Proposal.DrawBirth(rng, counts) is not { } t) continue;
                a = off.ProposeBirth(t);
                b = on.ProposeBirth(t);
            }
            else if (kind < 0.55 && m > 0)
            {
                int i = (int)(rng.NextUInt64() % (ulong)m);
                a = off.ProposeDeath(i);
                b = on.ProposeDeath(i);
            }
            else if (m > 0)
            {
                int i = (int)(rng.NextUInt64() % (ulong)m);
                if (model.Proposal.ProposeChange(rng, off.Terms[i]).Term is not { } t || t.Equals(off.Terms[i])) continue;
                a = off.ProposeChange(i, t);
                b = on.ProposeChange(i, t);
            }
            else continue;
            proposals++;
            Assert.Equal(a.Refusal, b.Refusal);
            if (a.Admitted)
            {
                admitted++;
                Assert.Equal(a.Columns, b.Columns);
                worstRss = Math.Max(worstRss, Rel(a.ResidualSquares, b.ResidualSquares));
                worstFit = Math.Max(worstFit, Rel(a.FittedSquares, b.FittedSquares));
                worstAngle = Math.Max(worstAngle, Rel(a.MinAngle, b.MinAngle));
                if (!double.IsNaN(a.NewestAngle)) worstAngle = Math.Max(worstAngle, Rel(a.NewestAngle, b.NewestAngle));
                if (rng.NextDouble() < 0.5)
                {
                    off.Accept();
                    on.Accept();
                    accepted++;
                }
            }
        }
        Assert.True(on.TableLookups > 0, "the tables were never read");
        Assert.True(worstRss <= 1e-8 && worstFit <= 1e-8 && worstAngle <= 1e-6, $"worst relative differences: ‖r‖² {worstRss:E2}, Q {worstFit:E2}, angles {worstAngle:E2}");
        Assert.Equal(off.Terms.Select(t => t.ToString()), on.Terms.Select(t => t.ToString()));
        double scale = Math.Max(1.0, off.Coefficients.ToArray().Max(Math.Abs));
        for (int j = 0; j < off.Count; j++) Assert.True(Math.Abs(off.Coefficients[j] - on.Coefficients[j]) <= 1e-8 * scale, $"β̂ {j}");
        var fa = (BorderedCholesky)off.Factor;
        var fb = (BorderedCholesky)on.Factor;
        double worstEntry = 0.0;
        for (int i = 0; i < off.Count; i++)
            for (int j = 0; j < off.Count; j++) worstEntry = Math.Max(worstEntry, Math.Abs(fa.Entry(i, j) - fb.Entry(i, j)) / Math.Sqrt(fa.Entry(i, i) * fa.Entry(j, j)));
        Assert.True(worstEntry <= 1e-11, $"a Gram entry differs by {worstEntry:E2} of the columns' norms");
        _out.WriteLine($"{(weighted ? "weighted" : "unit")}{(ridge ? ", ridge" : "")}{(categorical ? ", categorical" : "")}: {proposals} proposals, {admitted} admitted, {accepted} accepted; {on.TableLookups} table entries from {on.TableBuilds} tables, " +
                       $"{on.TableFallbacks} proposals re-decided by the dots; worst relative ‖r‖² {worstRss:E1}, Q {worstFit:E1}, angles {worstAngle:E1}, Gram entry {worstEntry:E1}");
    }

    // Where tables cannot apply the design keeps the dots, and a proposal is the same bit for bit whether tables are asked for or not: a
    // prior apart from the Gram (each term also needs its prior-side column's products) is one.
    [Fact]
    public void Design_PriorApartFromTheGram_KeepsTheDots()
    {
        BassModel model = Model(600, 6, 12, weighted: true, out _, gram: WeightedPriorGram.Unweighted);
        BassDesign off = model.NewDesign(), on = model.NewDesign();
        off.TablesOn = false;
        on.TablesOn = true;
        var rng = new Xoshiro256PlusPlus(seed: 5400);
        for (int step = 0; step < 800; step++)
        {
            var counts = new BassTermCounts(model.Prior);
            foreach (BassTerm t in off.Terms) counts.Add(t);
            if (off.Count >= model.MaxBasis || model.Proposal.DrawBirth(rng, counts) is not { } born) continue;
            DesignProposal a = off.ProposeBirth(born), b = on.ProposeBirth(born);
            Assert.Equal(BitConverter.DoubleToInt64Bits(a.ResidualSquares), BitConverter.DoubleToInt64Bits(b.ResidualSquares));
            Assert.Equal(BitConverter.DoubleToInt64Bits(a.MinAngle), BitConverter.DoubleToInt64Bits(b.MinAngle));
            if (a.Admitted && rng.NextDouble() < 0.5)
            {
                off.Accept();
                on.Accept();
            }
        }
        Assert.Equal(0, on.TableLookups);
        Assert.True(off.Count >= 3);
    }

    [Fact]
    public void Design_BorderFromTables_IsTheDotsWithinTheBound()
    {
        BassModel model = Model(1000, 8, 20, weighted: false, out _);
        BassDesign design = model.NewDesign();
        design.TablesOn = true;
        var rng = new Xoshiro256PlusPlus(seed: 5300);
        for (int step = 0; step < 600 && design.Count < 20; step++)
        {
            var counts = new BassTermCounts(model.Prior);
            foreach (BassTerm t in design.Terms) counts.Add(t);
            if (model.Proposal.DrawBirth(rng, counts) is { } born && design.ProposeBirth(born).Admitted) design.Accept();
        }
        Assert.True(design.Count >= 10);
        double worst = 0.0;
        int compared = 0;
        for (int trial = 0; trial < 400; trial++)
        {
            var one = new BassTerm([new HingeFactor((int)(rng.NextUInt64() % 8), rng.NextDouble() < 0.5 ? 1 : -1, (int)(rng.NextUInt64() % 1000))]);
            if (design.BorderComparison(one) is not var (tabled, dotted)) continue;
            for (int j = 0; j < tabled.Length; j++)
            {
                worst = Math.Max(worst, Math.Abs(tabled[j] - dotted[j]) / Math.Max(Math.Sqrt(((BorderedCholesky)design.Factor).Entry(j, j)), 1e-300));
                compared++;
            }
        }
        Assert.True(compared > 0);
        _out.WriteLine($"{compared} border entries; the worst difference is {worst:E2} of √G_jj (with ‖u‖ ≤ √n)");
        Assert.True(worst <= 4000 * 1.1102230246251565e-16 * Math.Sqrt(1000));
    }

    // ---- The chain's law ----

    private static BassChain EnumerationChain(BassModel model, int c, bool tables)
    {
        var chain = new BassChain(model, [], Math.Log(0.1), new Xoshiro256PlusPlus(400 + c), new Xoshiro256PlusPlus(500 + c));
        chain.Structure.DesignAt(chain.Current).TablesOn = tables;
        return chain;
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Chain_WithTables_WalksTheEnumeratedPosterior(bool weighted)
    {
        BassData data = BassSweepTests.Data();
        BassConfig config = BassWeightedCoreTests.EnumerationConfig(Sigma2Prior.InverseGamma(2.0, 0.3));
        BassModel model = weighted
            ? BassModel.Create(config, data, new BassPriorReferenceTests.Snapshot(r => new WorkingResponse(r, weights: BassWeightedCoreTests.Weights(r.Observations, seed: 21), priorGram: WeightedPriorGram.Weighted)))
            : BassModel.Create(config, data);
        BassKernelEnumeration.Space space = BassKernelEnumeration.Enumerate(model);
        var designs = new List<BassDesign>();
        BassWeightedCoreTests.CheckVisits(space, space.Posterior(), c =>
        {
            BassChain chain = EnumerationChain(model, c, tables: true);
            designs.Add(chain.Structure.DesignAt(chain.Current));
            return chain;
        }, _out, $"tables on, {(weighted ? "weighted" : "unit")}");
        Assert.True(designs.Sum(d => d.TableLookups) > 0, "the tables were never read");
        _out.WriteLine($"{designs.Sum(d => d.TableLookups)} table entries read from {designs.Sum(d => d.TableBuilds)} tables; {designs.Sum(d => d.TableFallbacks)} proposals re-decided by the dots");
    }

    private static (string[] Structures, double[] LogDensity, double[] LogTau, long Lookups, long Fallbacks) Run(BassModel model, bool tables, int sweeps, int chainIndex = 0)
    {
        BassChain chain = BassMixing.Chain(model, chainIndex);
        BassDesign design = chain.Structure.DesignAt(chain.Current);
        design.TablesOn = tables;
        var structures = new string[sweeps];
        var density = new double[sweeps];
        var tau = new double[sweeps];
        for (int s = 0; s < sweeps; s++)
        {
            chain.Step();
            structures[s] = string.Join("|", chain.Current.Terms.Select(t => t.ToString()));
            density[s] = chain.Moves.CurrentLogDensity;
            tau[s] = chain.Current.LogTau;
        }
        return (structures, density, tau, design.TableLookups, design.TableFallbacks);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void Chain_PairedSeeds_TablesOnAndOff(bool weighted)
    {
        BassModel model = Model(600, 6, 15, weighted, out _);
        const int sweeps = 6000;
        var off = Run(model, tables: false, sweeps);
        var on = Run(model, tables: true, sweeps);
        int firstDifference = -1;
        double worstDensity = 0.0, worstTau = 0.0;
        for (int s = 0; s < sweeps; s++)
        {
            if (off.Structures[s] != on.Structures[s] && firstDifference < 0) firstDifference = s;
            if (firstDifference >= 0) break;
            worstDensity = Math.Max(worstDensity, Rel(off.LogDensity[s], on.LogDensity[s]));
            worstTau = Math.Max(worstTau, Math.Abs(off.LogTau[s] - on.LogTau[s]));
        }
        _out.WriteLine($"{(weighted ? "weighted" : "unit")}: {sweeps} sweeps, structures {(firstDifference < 0 ? "identical throughout" : $"differ first at sweep {firstDifference}")}; " +
                       $"worst relative ln density difference {worstDensity:E2}, worst |Δ ln τ| {worstTau:E2}; M ends at {on.Structures[^1].Split('|').Length}; " +
                       $"{on.Lookups} table entries read, {on.Fallbacks} proposals re-decided by the dots");
        Assert.True(on.Lookups > 0 && off.Lookups == 0, "tables on should read tables and off should not");
        Assert.True(firstDifference < 0, $"the paths part at sweep {firstDifference}");
        Assert.True(worstDensity <= 1e-9 && worstTau <= 1e-8);
    }

    [Fact]
    public void Ladder_PairedSeeds_TablesOnAndOff()
    {
        BassModel model = Model(600, 6, 12, weighted: false, out _);
        (BassState start, int[] seeds) = BassMixing.Start(model, 0);
        const int levels = 4, scans = 500;
        double hot = new BassLadder(model, new TemperingSpec(), start.Terms, start.LogTau, seeds[3]).Tempering.InverseTemperatures.Min();
        double[] betas = [.. Enumerable.Range(0, levels).Select(k => Math.Pow(hot, k / (levels - 1.0)))];
        long lookups = 0;
        string[][] Trace(bool tables)
        {
            var ladder = new BassLadder(model, new TemperingSpec { Levels = levels, Adaptive = false }, betas, start.Terms, start.LogTau, seeds[3]);
            for (int j = 0; j < levels; j++) ladder.Structures[j].DesignAt(ladder.Tempering.Replicas[j].Current).TablesOn = tables;
            var trace = new string[scans][];
            for (int s = 0; s < scans; s++)
            {
                ladder.Step();
                trace[s] = [.. ladder.Tempering.Replicas.Select(r => string.Join("|", r.Current.Terms.Select(t => t.ToString())))];
            }
            if (tables) lookups = Enumerable.Range(0, levels).Sum(j => ladder.Structures[j].DesignAt(ladder.Tempering.Replicas[j].Current).TableLookups);
            return trace;
        }
        string[][] off = Trace(false), on = Trace(true);
        int first = Enumerable.Range(0, scans).FirstOrDefault(s => !off[s].SequenceEqual(on[s]), -1);
        _out.WriteLine($"tempered, {levels} rungs, {scans} scans: replicas {(first < 0 ? "identical throughout" : $"differ first at scan {first}")}; {lookups} table entries read");
        Assert.True(lookups > 0, "the tables were never read");
        Assert.True(first < 0, $"the ladders part at scan {first}");
    }

    // The cost grid's base point, eight chains per path, long: the brief's check (c). Paired seeds make it a stronger reading than two
    // laws compared by their means: the paths either keep one structure sequence or part, and where they part the means are compared.
    [Fact]
    public void Chain_AtTheBasePoint_PairedSeeds_Long()
    {
        if (!BassQualification.Enabled) return;
        var point = BassCost.Base;
        BassModel model = BassModel.Create(BassCost.Config(point, chains: 8), BassCost.Data(point));
        const int sweeps = 20000, chains = 8;
        var parted = new List<string>();
        double worstDensity = 0.0;
        long lookups = 0, fallbacks = 0;
        for (int c = 0; c < chains; c++)
        {
            var off = Run(model, tables: false, sweeps, c);
            var on = Run(model, tables: true, sweeps, c);
            lookups += on.Lookups;
            fallbacks += on.Fallbacks;
            int first = -1;
            for (int s = 0; s < sweeps && first < 0; s++)
            {
                if (off.Structures[s] != on.Structures[s]) first = s;
                else worstDensity = Math.Max(worstDensity, Rel(off.LogDensity[s], on.LogDensity[s]));
            }
            if (first >= 0)
            {
                parted.Add($"chain {c} at sweep {first}");
                double mOff = off.Structures.Skip(first).Average(t => t.Split('|').Length), mOn = on.Structures.Skip(first).Average(t => t.Split('|').Length);
                double dOff = off.LogDensity.Skip(first).Average(), dOn = on.LogDensity.Skip(first).Average();
                _out.WriteLine($"chain {c} parts at sweep {first}: mean M {mOff:F2} / {mOn:F2}, mean ln density {dOff:F2} / {dOn:F2}");
            }
        }
        _out.WriteLine($"base point, {chains} chains × {sweeps} sweeps, tables on against off: {(parted.Count == 0 ? "identical structure sequences throughout" : "parted: " + string.Join("; ", parted))}; " +
                       $"worst relative ln density difference before parting {worstDensity:E2}; {lookups} table entries read, {fallbacks} proposals re-decided by the dots");
        Assert.True(lookups > 0);
        Assert.True(worstDensity <= 1e-9);
    }

    // ---- Purity ----

    [Fact]
    public void Chain_RestoredWithoutTables_ContinuesBitForBit()
    {
        BassModel model = Model(600, 6, 15, weighted: false, out _);
        BassChain live = BassMixing.Chain(model, 0);
        live.Structure.DesignAt(live.Current).TablesOn = true;
        for (int s = 0; s < 500; s++) live.Step();
        BassChainSnapshot snapshot = live.Snapshot();
        BassChain restored = BassChain.Restore(model, snapshot);
        restored.Structure.DesignAt(restored.Current).TablesOn = true;
        long built = live.Structure.DesignAt(live.Current).TableBuilds;
        Assert.True(built > 0, "the live chain never built a table");
        for (int s = 0; s < 800; s++)
        {
            live.Step();
            restored.Step();
            Assert.Equal(BitConverter.DoubleToInt64Bits(live.Moves.CurrentLogDensity), BitConverter.DoubleToInt64Bits(restored.Moves.CurrentLogDensity));
            Assert.Equal(BitConverter.DoubleToInt64Bits(live.Current.LogTau), BitConverter.DoubleToInt64Bits(restored.Current.LogTau));
            Assert.True(live.Current.Terms.SequenceEqual(restored.Current.Terms), $"sweep {s}");
        }
        _out.WriteLine($"restored without tables, 800 sweeps bit for bit; the live chain had built {built} tables by the snapshot");
    }
}
