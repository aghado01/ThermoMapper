using System;

namespace Maths.Regression.Conformal;

/// <summary>
/// An order-statistic tree: an AVL tree whose nodes count their subtrees (Cormen, Leiserson, Rivest and Stein, ch. 14). It inserts
/// a key and ranks it among those held before, or ranks any key without inserting it, as (Below, Tied), each in O(log n) and one
/// walk down. Equal keys share one node with its multiplicity. The conformal test martingale's p-values (E35) read it, each new
/// score ranked among every score so far, where a sorted array built once serves the split form's fixed calibration. Keys order by
/// their <see cref="IComparable{T}"/>; a caller holding keys that do not order, such as NaN, refuses them first. Each node's fields
/// sit together in one array of structs, which doubles as it fills, so a walk reads one cache line per level. Not thread-safe.
/// </summary>
internal sealed class OrderStatisticTree<T> where T : IComparable<T>
{
    private struct Node
    {
        public T Key;
        public int Count, Size, Height, Left, Right;
    }

    private Node[] _node = new Node[16];
    private int _nodes, _root = -1;

    /// <summary>The keys held, each counted with its multiplicity.</summary>
    public int Count => Size(_root);

    /// <summary>The distinct keys held.</summary>
    public int Distinct => _nodes;

    /// <summary>The tree's height, at most about 1.44·log₂(Distinct + 2).</summary>
    public int Height => HeightOf(_root);

    /// <summary>Inserts <paramref name="key"/>, and gives the keys held before it strictly below it and equal to it, O(log n).</summary>
    public (int Below, int Tied) Add(T key)
    {
        int below = 0, tied = 0;
        _root = Insert(_root, key, ref below, ref tied);
        return (below, tied);
    }

    /// <summary>The keys held strictly below <paramref name="key"/>, and those equal to it, O(log n).</summary>
    public (int Below, int Tied) Rank(T key)
    {
        int below = 0, node = _root;
        while (node >= 0)
        {
            ref Node at = ref _node[node];
            int c = key.CompareTo(at.Key);
            if (c < 0) node = at.Left;
            else if (c > 0)
            {
                below += Size(at.Left) + at.Count;
                node = at.Right;
            }
            else return (below + Size(at.Left), at.Count);
        }
        return (below, 0);
    }

    private int Size(int node) => node < 0 ? 0 : _node[node].Size;

    private int HeightOf(int node) => node < 0 ? 0 : _node[node].Height;

    private int Insert(int node, T key, ref int below, ref int tied)
    {
        if (node < 0) return Leaf(key);
        int c = key.CompareTo(_node[node].Key);
        if (c == 0)
        {
            below += Size(_node[node].Left);
            tied = _node[node].Count++;
            _node[node].Size++;
            return node;
        }
        // The child is written back after the call, through the array as it then stands: an insert below may grow it.
        if (c < 0)
        {
            int left = Insert(_node[node].Left, key, ref below, ref tied);
            _node[node].Left = left;
        }
        else
        {
            below += Size(_node[node].Left) + _node[node].Count;
            int right = Insert(_node[node].Right, key, ref below, ref tied);
            _node[node].Right = right;
        }
        return Balance(node);
    }

    private int Leaf(T key)
    {
        if (_nodes == _node.Length) Array.Resize(ref _node, 2 * _node.Length);
        int node = _nodes++;
        _node[node] = new Node { Key = key, Count = 1, Size = 1, Height = 1, Left = -1, Right = -1 };
        return node;
    }

    private void Update(int node)
    {
        ref Node at = ref _node[node];
        at.Size = Size(at.Left) + Size(at.Right) + at.Count;
        at.Height = 1 + Math.Max(HeightOf(at.Left), HeightOf(at.Right));
    }

    private int Balance(int node)
    {
        Update(node);
        int left = _node[node].Left, right = _node[node].Right;
        int lean = HeightOf(left) - HeightOf(right);
        if (lean > 1)
        {
            if (HeightOf(_node[left].Left) < HeightOf(_node[left].Right)) _node[node].Left = RotateLeft(left);
            return RotateRight(node);
        }
        if (lean < -1)
        {
            if (HeightOf(_node[right].Right) < HeightOf(_node[right].Left)) _node[node].Right = RotateRight(right);
            return RotateLeft(node);
        }
        return node;
    }

    private int RotateRight(int node)
    {
        int left = _node[node].Left;
        _node[node].Left = _node[left].Right;
        _node[left].Right = node;
        Update(node);
        Update(left);
        return left;
    }

    private int RotateLeft(int node)
    {
        int right = _node[node].Right;
        _node[node].Right = _node[right].Left;
        _node[right].Left = node;
        Update(node);
        Update(right);
        return right;
    }
}
