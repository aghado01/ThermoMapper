#nullable enable

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Classifies the endpoint of a persistence interval, specifically required
/// for Zigzag Persistent Homology where intervals can have open or closed ends
/// corresponding to add or delete operations.
/// </summary>
public enum IntervalEnd
{
    /// <summary>
    /// The interval endpoint is closed (inclusive).
    /// Birth closed: an add operation.
    /// Death closed: a backward/delete operation (per FastZigzag Def 2).
    /// Default for standard Persistent Homology intervals.
    /// </summary>
    Closed,

    /// <summary>
    /// The interval endpoint is open (exclusive).
    /// </summary>
    Open
}
