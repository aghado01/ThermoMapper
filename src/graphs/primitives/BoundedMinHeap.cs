using System;

namespace Graphs.Primitives;

/// <summary>
/// Fixed-capacity max-heap of the K smallest distances seen so far.
/// Used by directed KNN and core-distance computation.
/// </summary>
public sealed class BoundedMinHeap
{
    private readonly Neighbor[] _heap;
    private int _count;
    private readonly int _capacity;

    public BoundedMinHeap(int k)
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        _capacity = k;
        _heap = new Neighbor[k];
    }

    public void TryAdd(int index, double distance)
    {
        if (_count < _capacity)
        {
            _heap[_count] = new Neighbor { Index = index, Distance = distance };
            _count++;
            if (_count == _capacity) BuildHeap();
        }
        else if (distance < _heap[0].Distance)
        {
            _heap[0] = new Neighbor { Index = index, Distance = distance };
            SiftDown(0);
        }
    }

    public Neighbor[] GetSorted()
    {
        var result = new Neighbor[_count];
        Array.Copy(_heap, result, _count);
        Array.Sort(result, (a, b) => a.Distance.CompareTo(b.Distance));
        return result;
    }

    private void BuildHeap()
    {
        for (int i = _count / 2 - 1; i >= 0; i--)
            SiftDown(i);
    }

    private void SiftDown(int i)
    {
        while (true)
        {
            int largest = i;
            int left = 2 * i + 1;
            int right = 2 * i + 2;

            if (left < _count && _heap[left].Distance > _heap[largest].Distance) largest = left;
            if (right < _count && _heap[right].Distance > _heap[largest].Distance) largest = right;
            if (largest == i) break;

            var tmp = _heap[i];
            _heap[i] = _heap[largest];
            _heap[largest] = tmp;
            i = largest;
        }
    }
}
