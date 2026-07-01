using System;
using System.Collections.Generic;
using System.Linq;
using Clustering.Graphical.SPC.Profiling;

namespace Clustering.Graphical.SPC.Profiling.Signals;

public sealed class ChiPeakSignalAnalyzer : ISignalAnalyzer
{
    public ProfileCriteria Analyze(SweepProfile profile)
    {
        if (profile is null)
            throw new ArgumentNullException(nameof(profile));

        if (profile.Count == 0)
            throw new ArgumentException("SweepProfile must contain at least one point.", nameof(profile));

        SpPlateauResult plateau = SpcProfileAnalysis.SpPlateau(profile);
        var band = SpcProfileAnalysis.ComputeHalfMaximumBand(profile);
        double stability = SpcProfileAnalysis.ComputeStability(profile);
        double peakChi = profile.Susceptibility.Count > 0
            ? profile.Susceptibility.Max()
            : 0.0;
        double bondEntropyVal = double.NaN;
        if (profile.BondEntropy != null && profile.Temperatures.Count > 0)
        {
            int bestIdx = 0;
            double bestDist = double.MaxValue;
            for (int i = 0; i < profile.Temperatures.Count; i++)
            {
                double dist = Math.Abs(profile.Temperatures[i] - plateau.TClus);
                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestIdx = i;
                }
            }
            bondEntropyVal = profile.BondEntropy[bestIdx];
        }

        return new ProfileCriteria(
            AnchorTemperature: plateau.TFs,
            AnchorBand: band,
            RefinedTemperature: plateau.TClus,
            CorroborationScore: stability,
            Enrichments: new Dictionary<string, double>
            {
                ["Stability"]   = stability,
                ["BandWidth"]   = band.Hi - band.Lo,
                ["PeakChi"]     = peakChi,
                ["BondEntropy"] = bondEntropyVal,
                ["TFs"]         = plateau.TFs,
                ["TPs"]         = plateau.TPs,
                ["TClus"]       = plateau.TClus,
                ["CliffFound"]  = plateau.CliffFound ? 1.0 : 0.0,
            });
    }
}
