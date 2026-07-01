namespace Graphs.Primitives
{
    /// <summary>
    /// Path-compressed union-find with union-by-size.
    /// Used for graph connectivity diagnostics in GraphBuilder.Validate.
    /// </summary>
    public sealed class UnionFind
    {
        private readonly int[] _parent;
        private readonly int[] _size;

        public UnionFind(int count)
        {
            _parent = new int[count];
            _size = new int[count];
            for (int i = 0; i < count; i++)
            {
                _parent[i] = i;
                _size[i] = 1;
            }
        }

        public int Find(int item)
        {
            int root = item;
            while (root != _parent[root])
                root = _parent[root];

            // Path compression
            while (item != root)
            {
                int next = _parent[item];
                _parent[item] = root;
                item = next;
            }

            return root;
        }

        public void Union(int a, int b)
        {
            int ra = Find(a);
            int rb = Find(b);
            if (ra == rb) return;

            if (_size[ra] < _size[rb])
            {
                _parent[ra] = rb;
                _size[rb] += _size[ra];
            }
            else
            {
                _parent[rb] = ra;
                _size[ra] += _size[rb];
            }
        }

        /// <summary>
        /// Returns the size of the component whose root is <paramref name="root"/>.
        /// Caller must pass a root (e.g. the return value of <see cref="Find"/>);
        /// no Find is performed here.
        /// </summary>
        public int Size(int root) => _size[root];

        /// <summary>
        /// Relabels the current root of <paramref name="member"/>'s component as
        /// <paramref name="newId"/> and stamps the component size. Used by HDBSCAN's
        /// Kruskal pass to assign fresh internal-node ids (in [N, 2N-2]) after each
        /// merge. The structure must have been allocated with capacity for
        /// <paramref name="newId"/>.
        /// </summary>
        public void Reroot(int member, int newId, int newSize)
        {
            int r = Find(member);
            _parent[newId] = newId;
            _parent[r] = newId;
            _size[newId] = newSize;
        }

        /// <summary>
        /// Resets all nodes to singleton components without reallocation.
        /// Used by PottsModel to reuse the structure across Swendsen-Wang sweeps.
        /// </summary>
        public void Reset()
        {
            int count = _parent.Length;
            for (int i = 0; i < count; i++)
            {
                _parent[i] = i;
                _size[i] = 1;
            }
        }

        /// <summary>
        /// Returns the canonical root label for each node.
        /// </summary>
        public int[] GetLabels()
        {
            int count = _parent.Length;
            var labels = new int[count];
            for (int i = 0; i < count; i++)
                labels[i] = Find(i);
            return labels;
        }

        /// <summary>
        /// Walks the parent array once and writes the size of each component
        /// (one entry per UF root) into <paramref name="sizes"/>. Returns the
        /// number of roots written. <paramref name="sizes"/> must have capacity
        /// at least N.
        /// </summary>
        /// <remarks>
        /// O(N), zero allocation. Sizes are written in the order roots are
        /// encountered while scanning indices <c>0..N-1</c>; this ordering
        /// carries no semantic meaning — callers should treat the output as a
        /// multiset. Relies on <c>_size[root]</c> being maintained at union
        /// time, which is why this UF carries union-by-size in the first place.
        /// </remarks>
        public int WriteRootSizesTo(System.Span<int> sizes)
        {
            int n = _parent.Length;
            if (sizes.Length < n)
                throw new System.ArgumentException("Destination span must have capacity at least N.", nameof(sizes));

            int[] parent = _parent;
            int[] size = _size;

            int rootCount = 0;
            for (int i = 0; i < n; i++)
            {
                if (parent[i] == i)
                    sizes[rootCount++] = size[i];
            }
            return rootCount;
        }
    }
}
