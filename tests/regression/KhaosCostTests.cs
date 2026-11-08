using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bass;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.Regression.Tests;

/// <summary>
/// KHAOS's cost row (the KHAOS brief's slice 1), behind <see cref="BassQualification"/>, on <see cref="BassCost.Tick"/>'s clock
/// (<c>BASS_COST_CLOCK</c> = cycles for the thread's cycle counter, indicative under load): a Legendre product column against a hinge
/// product column at the cost report's n; and a sweep at the report's base (its rich surface, n 1,000, p 10, M_max 20, order 2) with
/// d_max 15, and at p 20, d_max 15, q_max 5, under KHAOS's raw ridge and the normalized ridge, beside the hinge fit's at the base.
/// Each reading is the median over nine batches of a batch's mean per call, after a warm-up batch and the chain's burn-in.
/// </summary>
public sealed class KhaosCostTests(ITestOutputHelper output)
{
    private static double Median(double[] values)
    {
        double[] sorted = [.. values.Order()];
        return sorted[sorted.Length / 2];
    }

    // µs per call of each of `calls`, round-robin: nine batches each after a warm-up round, the median per call, so a change of load
    // between them reaches every reading alike.
    private static double[] RoundRobin(Action[] calls, int size)
    {
        var batches = new double[calls.Length][];
        for (int c = 0; c < calls.Length; c++) batches[c] = new double[10];
        for (int b = 0; b < 10; b++)
            for (int c = 0; c < calls.Length; c++)
            {
                long start = BassCost.Tick();
                for (int k = 0; k < size; k++) calls[c]();
                batches[c][b] = (BassCost.Tick() - start) / BassCost.TicksPerSecond / size * 1e6;
            }
        return [.. batches.Select(t => Median(t[1..]))];
    }

    // The first form's arithmetic, for the record: each point's chain of α steps, dividing at each, here scalar.
    private static void PerPoint(int degree, ReadOnlySpan<double> x, Span<double> column)
    {
        double root = Math.Sqrt(2.0 * degree + 1.0);
        for (int i = 0; i < column.Length; i++)
        {
            double u = 2.0 * x[i] - 1.0, p0 = 1.0, p1 = u;
            for (int k = 1; k < degree; k++)
            {
                double p2 = ((2 * k + 1) * u * p1 - k * p0) / (k + 1);
                p0 = p1;
                p1 = p2;
            }
            column[i] *= root * p1;
        }
    }

    [Fact]
    public void ColumnsAndSweeps_OnTheClock()
    {
        if (!BassQualification.Enabled) return;   // heavy
        BassData data = BassCost.Rich(BassCost.BaseN, BassCost.BaseP);
        BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
        var basis = new BassBasis(inputs, 1);
        var column = new double[inputs.Cases.Cases];
        var legendre = new BassTerm(legendre: [new(0, 8), new(3, 7)]);
        var hinge = new BassTerm([new HingeFactor(0, 1, 400), new HingeFactor(3, -1, 600)]);
        var psi15 = new BassTerm(legendre: [new(0, 15)]);
        double[] t = RoundRobin([
            () => basis.CaseColumn(legendre, inputs.Cases, column),
            () => basis.CaseColumn(psi15, inputs.Cases, column),
            () => basis.CaseColumn(hinge, inputs.Cases, column),
            () => { column.AsSpan().Fill(1.0); PerPoint(15, inputs.Cases.Continuous(0), column); },
        ], 2000);
        output.WriteLine($"clock {(BassCost.Cycles ? "cycles" : "wall")}; n {inputs.Cases.Cases}, round-robin: Legendre ψ8·ψ7 column {t[0]:F2} µs, ψ15 {t[1]:F2} µs " +
                         $"(scalar per point, dividing: {t[3]:F2}); hinge product of two {t[2]:F2} µs; ψ8·ψ7 over the hinges {t[0] / t[2]:F2}");

        // The sweeps, each chain burned, then read round-robin.
        var names = new List<string>();
        var chains = new List<BassChain>();
        BassCoefficientPrior raw = BassCoefficientPrior.RawRidge(TauPrior.Held(1e-5)), normalized = BassCoefficientPrior.Ridge(TauPrior.UnitInformationGamma());
        var modified = new BassCoefficientPrior(new ModifiedGForm(1.0), TauPrior.Gamma(0.001, 1000));   // khaos's g-prior, g₀² drawn
        foreach ((string name, int p, int order, BassCoefficientPrior coefficients) in new[] {
                     ("base, raw ridge", BassCost.BaseP, BassCost.BaseOrder, raw), ("base, normalized ridge", BassCost.BaseP, BassCost.BaseOrder, normalized),
                     ("base, modified g-prior", BassCost.BaseP, BassCost.BaseOrder, modified),
                     ("p 20, q_max 5, raw ridge", 20, 5, raw), ("p 20, q_max 5, normalized ridge", 20, 5, normalized), ("p 20, q_max 5, modified g-prior", 20, 5, modified) })
        {
            BassData d = BassCost.Rich(BassCost.BaseN, p);
            var config = new BassConfig
            {
                Structure = new BassStructureSpec { MaxInteraction = order, MaxBasis = BassCost.BaseCap, Legendre = new LegendreBasis { MaxDegree = 15 } },
                Coefficients = coefficients,
                Intercept = coefficients.Form is RidgeForm ? GPriorCentre.Mean : GPriorCentre.Zero,
                Sampling = new BassSampling { Chains = 1, MasterSeed = 3 },
            };
            names.Add(name);
            chains.Add(BassMixing.Chain(BassModel.Create(config, d), 0));
        }
        names.Add("base, hinges, g-form");
        chains.Add(BassMixing.Chain(BassModel.Create(new BassConfig { Structure = new BassStructureSpec { MaxInteraction = BassCost.BaseOrder, MaxBasis = BassCost.BaseCap }, Sampling = new BassSampling { Chains = 1, MasterSeed = 3 } }, data), 0));
        foreach (BassChain chain in chains)
            for (int s = 0; s < 4000; s++) chain.Step();
        double[] sweeps = RoundRobin([.. chains.Select(chain => (Action)(() => chain.Step()))], 200);
        for (int c = 0; c < chains.Count; c++)
            output.WriteLine($"{names[c]}: a sweep {sweeps[c]:F1} µs at M {chains[c].Current.Count}; acceptance {string.Join(", ", chains[c].Structure.Kinds.Entries.Where(e => e.Attempts > 0).Select(e => $"{e.Kind} {(double)e.Accepted / e.Attempts:F3}"))}");
    }
}
