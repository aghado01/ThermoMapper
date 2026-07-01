using System.Collections.Generic;

namespace Clustering.Graphical.SPC.Profiling;

public sealed record ProfileCriteria(
    double AnchorTemperature,
    (double Lo, double Hi) AnchorBand,
    double RefinedTemperature,
    double CorroborationScore,
    IReadOnlyDictionary<string, double> Enrichments);
