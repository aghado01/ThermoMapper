using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC.Runtime.Execution;
using Clustering.Graphical.SPC.Runtime.Core.Sampler;
using Graphs.Observables;
using PottsObs = Graphs.Models.Potts.Observables;

namespace Clustering.Graphical.SPC.Profiling;

public sealed record SweepProfile(
    IReadOnlyList<double> Temperatures,
    IReadOnlyList<double> Susceptibility,
    IReadOnlyList<double> SpecificHeat,
    IReadOnlyList<double> LabelEntropy,
    IReadOnlyList<double>? BondEntropy,
    IReadOnlyDictionary<string, IReadOnlyList<double>> AdditionalChannels)
{
    public static SweepProfile Empty { get; } = new(
        Array.Empty<double>(),
        Array.Empty<double>(),
        Array.Empty<double>(),
        Array.Empty<double>(),
        null,
        new Dictionary<string, IReadOnlyList<double>>());

    public int Count => Temperatures.Count;

    public bool IsEmpty => Count == 0;

    public SweepProfile WithChannel(string name, IReadOnlyList<double> values)
    {
        if (name is null) throw new ArgumentNullException(nameof(name));
        if (values is null) throw new ArgumentNullException(nameof(values));
        if (values.Count != Count)
            throw new ArgumentException("Channel values must match temperature length.", nameof(values));

        var channels = AdditionalChannels.ToDictionary(kvp => kvp.Key, kvp => kvp.Value);
        channels[name] = values;
        return this with { AdditionalChannels = channels };
    }

    public static SweepProfile From(
        IReadOnlyList<SpcRunResult> runs,
        SusceptibilityKind kind = SusceptibilityKind.FkCluster)
    {
        ArgumentNullException.ThrowIfNull(runs);
        if (runs.Count == 0)
            return Empty;

        var frames = runs.Select(run => run.Accumulator).ToArray();

        // Sweep-tier assembly: each model reduction says what a frame means; the
        // assembler groups those per-draw values by T and averages the replicas.
        var (temperatures, chi) = SweepCurves.ByTemperature(
            frames, f => PottsObs.Susceptibility.Reduce(f.RunningSumSqClusterSizes, f.DrawCount, f.Spins.Length));
        // FK reduced: the same reduction over the giant-excluded sum (peaks vanish once the
        // giant percolates — the cleanest SP-transition detector). Linear, so reduce-then-average is exact.
        var (_, chiReduced) = SweepCurves.ByTemperature(
            frames, f => PottsObs.Susceptibility.Reduce(f.RunningSumSqClusterSizesExcl, f.DrawCount, f.Spins.Length));
        var (_, meanEnergy) = SweepCurves.ByTemperature(
            frames, f => PottsObs.MeanEnergy.Reduce(f.RunningSumEnergy, f.DrawCount));
        var (_, meanEnergySq) = SweepCurves.ByTemperature(
            frames, f => f.DrawCount > 0 ? f.RunningSumEnergySq / f.DrawCount : 0.0);
        var (_, meanMag) = SweepCurves.ByTemperature(
            frames, f => PottsObs.Magnetization.Reduce(f.RunningSumMag, f.RunningSumMagSq, f.DrawCount).Mean);
        var (_, magSecondMoment) = SweepCurves.ByTemperature(
            frames, f => PottsObs.Magnetization.Reduce(f.RunningSumMag, f.RunningSumMagSq, f.DrawCount).SecondMoment);

        // Variance is computed from the *averaged* moments (not averaged per
        // frame) — the order-parameter fluctuation of the equilibrium ensemble.
        var magVariance = new double[meanMag.Length];
        for (int i = 0; i < magVariance.Length; i++)
            magVariance[i] = Math.Max(0.0, magSecondMoment[i] - meanMag[i] * meanMag[i]);

        // Specific heat C_v = (⟨E²⟩−⟨E⟩²)/T², reduced from the *pooled* energy moments
        // (nonlinear — pool-then-reduce, matching the magnetization channels; not
        // per-frame-then-averaged, which drops the between-replica variance for Replicas>1).
        var cv = new double[temperatures.Length];
        for (int i = 0; i < cv.Length; i++)
            cv[i] = temperatures[i] > 0.0
                ? Math.Max(0.0, meanEnergySq[i] - meanEnergy[i] * meanEnergy[i]) / (temperatures[i] * temperatures[i])
                : 0.0;

        // Magnetization susceptibility χ_m = (N/T)(⟨m²⟩−⟨m⟩²), the literal BWD/Domany paper χ.
        // Reduced from the *pooled* moments above (nonlinear — never per-frame-then-averaged).
        int siteCount = frames[0].Spins.Length;
        var chiMag = new double[temperatures.Length];
        for (int i = 0; i < chiMag.Length; i++)
            chiMag[i] = PottsObs.MagnetizationSusceptibility.Reduce(
                meanMag[i], magSecondMoment[i], siteCount, temperatures[i]);

        // The configured kind drives peak-finding; all three ride along as channels below for
        // free comparison (FK vs magnetization disagree only in edge cases — now testable).
        IReadOnlyList<double> primaryChi = kind switch
        {
            SusceptibilityKind.FkReduced             => chiReduced,
            SusceptibilityKind.Magnetization         => chiMag,
            SusceptibilityKind.MagnetizationVariance => magVariance,
            _                                        => chi,
        };

        var labelEntropy = LabelEntropyCurve.From(frames);

        var channels = new Dictionary<string, IReadOnlyList<double>>(StringComparer.Ordinal)
        {
            ["MeanEnergy"] = meanEnergy,
            ["MeanMagnetization"] = meanMag,
            ["MagnetizationSecondMoment"] = magSecondMoment,
            ["MagnetizationVariance"] = magVariance,
            ["SusceptibilityFkCluster"] = chi,
            ["SusceptibilityFkReduced"] = chiReduced,
            ["SusceptibilityMagnetization"] = chiMag,
        };

        // Bond entropy: the dispersion of the bond-survival field. Mint the
        // model-agnostic currency from the SW accumulator by pooling the raw
        // bond counts across replicas per temperature, then compute entropy.
        IReadOnlyList<double>? bondEntropyCurve = null;
        if (frames.Length > 0 && frames[0].BondFormedCount is not null)
        {
            var bondEntropyArr = new double[temperatures.Length];
            var groupedFrames = frames.GroupBy(f => f.Temperature).OrderBy(g => g.Key).ToList();
            int tIndex = 0;
            foreach (var group in groupedFrames)
            {
                int csrLength = group.First().BondFormedCount!.Length;
                var pooledBonds = new int[csrLength];
                int pooledDraws = 0;
                
                foreach (var f in group)
                {
                    pooledDraws += f.DrawCount;
                    for (int i = 0; i < csrLength; i++)
                        pooledBonds[i] += f.BondFormedCount![i];
                }
                
                if (pooledDraws > 0)
                {
                    double invDraws = 1.0 / pooledDraws;
                    var pooledRates = new double[csrLength];
                    for (int i = 0; i < csrLength; i++)
                        pooledRates[i] = pooledBonds[i] * invDraws;
                    bondEntropyArr[tIndex] = AffinityEntropy.EntropyNats(pooledRates);
                }
                else
                {
                    bondEntropyArr[tIndex] = 0.0;
                }
                tIndex++;
            }
            bondEntropyCurve = bondEntropyArr;
        }

        return new SweepProfile(
            temperatures,
            primaryChi,
            cv,
            labelEntropy.Entropy,
            bondEntropyCurve,
            channels);
    }

}
