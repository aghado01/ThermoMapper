using System;
using Maths.Distance;

namespace Graphs.Distance;

/// <summary>
/// Graph-construction distance contract: extends the maths primitive with
/// <see cref="MetricProperties"/> for bandwidth dispatch and input validation.
/// </summary>
public interface IDistanceMetric : Maths.Distance.IDistanceMetric
{
    new MetricProperties Properties { get; }
}
