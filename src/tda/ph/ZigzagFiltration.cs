#nullable enable
using System;
using System.Collections;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Direction of a zigzag operation: either adding or deleting a cell.
/// </summary>
public enum ZigzagDirection
{
    Add,
    Delete
}

/// <summary>
/// Represents a single step in a zigzag filtration sequence.
/// </summary>
public readonly struct ZigzagStep
{
    /// <summary>
    /// Permanent global cell identifier. Uniquely identifies the cell across the entire zigzag sequence.
    /// </summary>
    public int GlobalCellId { get; }

    /// <summary>
    /// Whether this step adds or deletes the cell.
    /// </summary>
    public ZigzagDirection Direction { get; }

    /// <summary>
    /// For <see cref="ZigzagDirection.Add"/>, the boundary of the cell as indices into cells ALREADY present
    /// in the complex immediately prior to this addition. 
    /// For <see cref="ZigzagDirection.Delete"/>, this must be null.
    /// </summary>
    public int[]? BoundaryAtAdd { get; }

    public ZigzagStep(int globalCellId, ZigzagDirection direction, int[]? boundaryAtAdd)
    {
        if (direction == ZigzagDirection.Delete && boundaryAtAdd != null)
        {
            throw new ArgumentException("Delete steps must not carry a boundary payload.", nameof(boundaryAtAdd));
        }

        if (direction == ZigzagDirection.Add && boundaryAtAdd == null)
        {
            throw new ArgumentException("Add steps must carry a boundary payload.", nameof(boundaryAtAdd));
        }

        GlobalCellId = globalCellId;
        Direction = direction;
        BoundaryAtAdd = boundaryAtAdd;
    }
}

/// <summary>
/// A sequence of cell additions and deletions representing a zigzag filtration.
/// Convention: Starts and ends empty.
/// </summary>
public sealed class ZigzagFiltration : IReadOnlyList<ZigzagStep>
{
    private readonly List<ZigzagStep> _steps;

    public ZigzagFiltration()
    {
        _steps = new List<ZigzagStep>();
    }

    public ZigzagFiltration(IEnumerable<ZigzagStep> steps)
    {
        _steps = new List<ZigzagStep>(steps);
    }

    public int Count => _steps.Count;

    public ZigzagStep this[int index] => _steps[index];

    public void Add(ZigzagStep step)
    {
        _steps.Add(step);
    }

    public void Add(int globalCellId, int[] boundaryAtAdd)
    {
        _steps.Add(new ZigzagStep(globalCellId, ZigzagDirection.Add, boundaryAtAdd));
    }

    public void Delete(int globalCellId)
    {
        _steps.Add(new ZigzagStep(globalCellId, ZigzagDirection.Delete, null));
    }

    public IEnumerator<ZigzagStep> GetEnumerator() => _steps.GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}
