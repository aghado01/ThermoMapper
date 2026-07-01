using Clustering.Graphical.SPC.Profiling;

namespace Clustering.Graphical.SPC.Profiling.Signals;

public interface ISignalAnalyzer
{
    ProfileCriteria Analyze(SweepProfile profile);
}
