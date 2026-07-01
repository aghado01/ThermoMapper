#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Naive, obviously-correct oracle for computing Zigzag Persistent Homology.
/// Uses generalized-rank inclusion-exclusion to decompose the A_n module into intervals.
/// Designed for correctness and validation of faster algorithms, not for performance.
/// </summary>
public static class ZigzagBarcodeNaive
{
    public static Barcode Compute(ZigzagFiltration f, int maxDimension = int.MaxValue)
    {
        int numSteps = f.Count;
        if (numSteps == 0) return new Barcode(Array.Empty<Bar>(), "Zigzag");

        // Determine universe of cells
        int maxCellId = -1;
        foreach (var step in f)
        {
            if (step.GlobalCellId > maxCellId) maxCellId = step.GlobalCellId;
        }
        int N = maxCellId + 1;

        // Cell metadata (Dimension, Boundary)
        int[] dimOf = new int[N];
        int[][] bndOf = new int[N][];
        for (int i = 0; i < N; i++) bndOf[i] = Array.Empty<int>();

        // Materialize complexes K_0 ... K_m
        // K_0 is empty. K_{i} is after step i-1.
        var K = new List<HashSet<int>>();
        K.Add(new HashSet<int>());
        
        var currentK = new HashSet<int>();
        for (int i = 0; i < numSteps; i++)
        {
            var step = f[i];
            if (step.Direction == ZigzagDirection.Add)
            {
                currentK.Add(step.GlobalCellId);
                // We assume dimension is derived from boundary length or we need to pass dimension.
                // Wait, ZigzagStep doesn't include Dimension!
                // Let's deduce dimension from boundary. If boundary is empty, dim 0.
                // Else, dim is boundary_cell's dim + 1.
                int dim = 0;
                if (step.BoundaryAtAdd != null && step.BoundaryAtAdd.Length > 0)
                {
                    dim = dimOf[step.BoundaryAtAdd[0]] + 1;
                }
                dimOf[step.GlobalCellId] = dim;
                bndOf[step.GlobalCellId] = step.BoundaryAtAdd ?? Array.Empty<int>();
            }
            else
            {
                currentK.Remove(step.GlobalCellId);
            }
            K.Add(new HashSet<int>(currentK));
        }

        var isForward = new bool[numSteps];
        for (int i = 0; i < numSteps; i++)
            isForward[i] = f[i].Direction == ZigzagDirection.Add;

        var bars = new List<Bar>();
        for (int p = 0; p <= maxDimension; p++)
        {
            // Homology bases H_p(K_i) for every complex.
            var V = new List<List<bool[]>>();
            for (int i = 0; i <= numSteps; i++)
                V.Add(ComputeHomologyBasis(K[i], p, dimOf, bndOf, N));

            // Induced maps per step (here the maps are inclusions, so cycles carry over directly).
            var M = new List<bool[,]>();
            for (int i = 0; i < numSteps; i++)
                M.Add(isForward[i]
                    ? ComputeInducedMapMatrix(V[i], V[i + 1], K[i + 1], p + 1, dimOf, bndOf, N)
                    : ComputeInducedMapMatrix(V[i + 1], V[i], K[i], p + 1, dimOf, bndOf, N));

            bars.AddRange(DecomposeDimension(numSteps, V, M, isForward, p));
        }

        return new Barcode(bars, "Zigzag Step");
    }

    /// <summary>
    /// Decompose one dimension's zigzag module (bases V, induced maps M, arrow directions isForward)
    /// into 4-type interval bars via generalized-rank inclusion-exclusion. Shared by the naive oracle
    /// and the general simplicial-zigzag oracle (<see cref="ZigzagMapBarcode"/>).
    /// </summary>
    internal static List<Bar> DecomposeDimension(int numSteps, List<List<bool[]>> V, List<bool[,]> M, bool[] isForward, int p)
    {
        var bars = new List<Bar>();

        int[,] GR = new int[numSteps + 1, numSteps + 1];
        for (int b = 0; b <= numSteps; b++)
            for (int d = b; d <= numSteps; d++)
                GR[b, d] = ComputeGeneralizedRank(b, d, V, M, isForward);

        for (int b = 0; b <= numSteps; b++)
            for (int d = b; d <= numSteps; d++)
            {
                int mult = GR[b, d]
                         - (b > 0 ? GR[b - 1, d] : 0)
                         - (d < numSteps ? GR[b, d + 1] : 0)
                         + (b > 0 && d < numSteps ? GR[b - 1, d + 1] : 0);

                for (int m = 0; m < mult; m++)
                {
                    // birth closed iff the arrow into K_b is forward; death closed iff the arrow out
                    // of K_d is backward. Report by step: birth = creating step b-1, death = step d.
                    IntervalEnd bEnd = (b > 0 && isForward[b - 1]) ? IntervalEnd.Closed : IntervalEnd.Open;
                    IntervalEnd dEnd = (d < numSteps && !isForward[d]) ? IntervalEnd.Closed : IntervalEnd.Open;
                    bars.Add(new Bar(b - 1, d, p, null, null, null, bEnd, dEnd));
                }
            }

        return bars;
    }

    internal static List<bool[]> ComputeHomologyBasis(HashSet<int> K, int p, int[] dimOf, int[][] bndOf, int N)
    {
        var cellsP = K.Where(c => dimOf[c] == p).ToList();
        var cellsPMinus1 = K.Where(c => dimOf[c] == p - 1).ToList();
        var cellsPPlus1 = K.Where(c => dimOf[c] == p + 1).ToList();

        // Build boundary matrix d_p: rows = (p-1)-cells, cols = p-cells
        bool[,] dp = new bool[cellsPMinus1.Count, cellsP.Count];
        for (int j = 0; j < cellsP.Count; j++)
        {
            int cell = cellsP[j];
            foreach (int face in bndOf[cell])
            {
                int r = cellsPMinus1.IndexOf(face);
                if (r >= 0) dp[r, j] = true;
            }
        }

        // Z = ker d_p
        var Z_local = Z2LinearAlgebra.Nullspace(dp);
        var Z = new List<bool[]>();
        foreach (var z in Z_local)
        {
            bool[] vec = new bool[N];
            for (int j = 0; j < cellsP.Count; j++)
            {
                if (z[j]) vec[cellsP[j]] = true;
            }
            Z.Add(vec);
        }

        if (p == 0)
        {
            // For H0, ker d_0 is everything
            Z.Clear();
            foreach (int cell in cellsP)
            {
                bool[] vec = new bool[N];
                vec[cell] = true;
                Z.Add(vec);
            }
        }

        // B = im d_{p+1}
        var B = new List<bool[]>();
        foreach (int cell in cellsPPlus1)
        {
            bool[] vec = new bool[N];
            foreach (int face in bndOf[cell])
            {
                if (K.Contains(face)) vec[face] = true;
            }
            B.Add(vec);
        }

        // Basis for H_p = Z / B
        var reducedZ = Z2LinearAlgebra.ReduceModuloSpan(Z, B);
        var basis = new List<bool[]>();
        
        // We need to form a linearly independent subset of the original Z that spans Z/B.
        // We can do this by keeping a running span of B + selected elements.
        var currentSpan = new List<bool[]>(B);
        for (int i = 0; i < Z.Count; i++)
        {
            // Reduce Z[i] against currentSpan
            var check = Z2LinearAlgebra.ReduceModuloSpan(new List<bool[]> { Z[i] }, currentSpan);
            bool isZero = true;
            for (int k = 0; k < N; k++) if (check[0][k]) { isZero = false; break; }

            if (!isZero)
            {
                basis.Add(Z[i]);
                currentSpan.Add(Z[i]);
            }
        }

        return basis;
    }

    internal static bool[,] ComputeInducedMapMatrix(List<bool[]> sourceBasis, List<bool[]> targetBasis, HashSet<int> targetK, int pPlus1, int[] dimOf, int[][] bndOf, int N)
    {
        int cols = sourceBasis.Count;
        int rows = targetBasis.Count;
        bool[,] M = new bool[rows, cols];

        if (cols == 0 || rows == 0) return M;

        var cellsPPlus1 = targetK.Where(c => dimOf[c] == pPlus1).ToList();
        var B = new List<bool[]>();
        foreach (int cell in cellsPPlus1)
        {
            bool[] vec = new bool[N];
            foreach (int face in bndOf[cell])
            {
                if (targetK.Contains(face)) vec[face] = true;
            }
            B.Add(vec);
        }

        // To express a vector v in terms of targetBasis modulo B:
        // We form a matrix columns: [targetBasis | B], solve [targetBasis | B] * x = v
        int targetBasisSize = targetBasis.Count;
        int BSize = B.Count;
        int totalCols = targetBasisSize + BSize;

        // Find pivot rows for [targetBasis | B] to solve efficiently
        bool[,] sys = new bool[N, totalCols];
        for (int j = 0; j < targetBasisSize; j++)
        {
            for (int i = 0; i < N; i++) sys[i, j] = targetBasis[j][i];
        }
        for (int j = 0; j < BSize; j++)
        {
            for (int i = 0; i < N; i++) sys[i, targetBasisSize + j] = B[j][i];
        }

        // Row reduction of sys keeping track of operations to apply to rhs
        int[] pivotCols = new int[N];
        for (int i = 0; i < N; i++) pivotCols[i] = -1;

        int r = 0;
        for (int c = 0; c < totalCols && r < N; c++)
        {
            int pivot = -1;
            for (int i = r; i < N; i++)
            {
                if (sys[i, c]) { pivot = i; break; }
            }
            if (pivot == -1) continue;

            if (pivot != r)
            {
                for (int j = c; j < totalCols; j++) (sys[r, j], sys[pivot, j]) = (sys[pivot, j], sys[r, j]);
            }

            pivotCols[r] = c;

            for (int i = 0; i < N; i++)
            {
                if (i != r && sys[i, c])
                {
                    for (int j = c; j < totalCols; j++) sys[i, j] ^= sys[r, j];
                }
            }
            r++;
        }

        // Now for each source vector v, we reduce it using the same operations
        for (int j = 0; j < cols; j++)
        {
            bool[] v = (bool[])sourceBasis[j].Clone();
            
            // Wait, we need to apply the SAME row swaps to v that we applied to sys.
            // Since we didn't store row operations, let's just augment sys with all source vectors!
            // Actually, much easier: augment [targetBasis | B | sourceBasis] and reduce fully.
        }

        // Let's do the augmented matrix approach directly.
        int augCols = targetBasisSize + BSize + cols;
        bool[,] aug = new bool[N, augCols];
        for (int j = 0; j < targetBasisSize; j++)
            for (int i = 0; i < N; i++) aug[i, j] = targetBasis[j][i];
        for (int j = 0; j < BSize; j++)
            for (int i = 0; i < N; i++) aug[i, targetBasisSize + j] = B[j][i];
        for (int j = 0; j < cols; j++)
            for (int i = 0; i < N; i++) aug[i, targetBasisSize + BSize + j] = sourceBasis[j][i];

        r = 0;
        int[] pCols = new int[N];
        for (int i = 0; i < N; i++) pCols[i] = -1;

        for (int c = 0; c < targetBasisSize + BSize && r < N; c++)
        {
            int pivot = -1;
            for (int i = r; i < N; i++)
            {
                if (aug[i, c]) { pivot = i; break; }
            }
            if (pivot == -1) continue;

            if (pivot != r)
            {
                for (int j = c; j < augCols; j++) (aug[r, j], aug[pivot, j]) = (aug[pivot, j], aug[r, j]);
            }

            pCols[r] = c;

            for (int i = 0; i < N; i++)
            {
                if (i != r && aug[i, c])
                {
                    for (int j = c; j < augCols; j++) aug[i, j] ^= aug[r, j];
                }
            }
            r++;
        }

        // Now the equations are solved. For each source vector (column in augmented part),
        // its expression in terms of the pivot columns is directly readable.
        // Specifically, if pCols[i] = c, then the c-th variable has value aug[i, rhs_col].
        for (int j = 0; j < cols; j++)
        {
            int rhsCol = targetBasisSize + BSize + j;
            for (int i = 0; i < r; i++)
            {
                int c = pCols[i];
                if (c < targetBasisSize)
                {
                    M[c, j] = aug[i, rhsCol];
                }
            }
        }

        return M;
    }

    internal static int ComputeGeneralizedRank(int b, int d, List<List<bool[]>> V, List<bool[,]> M, bool[] isForward)
    {
        // Block dimensions and offsets in the direct sum D = ⊕_{k=b}^d V_k.
        int len = d - b + 1;
        int[] n = new int[len];
        int[] off = new int[len];
        int totalVars = 0;
        for (int t = 0; t < len; t++)
        {
            n[t] = V[b + t].Count;
            off[t] = totalVars;
            totalVars += n[t];
        }
        if (totalVars == 0) return 0;

        // Limit L = compatible tuples = nullspace of the arrow-compatibility equations.
        int totalEqs = 0;
        for (int k = b; k < d; k++)
            totalEqs += (isForward[k]) ? n[k + 1 - b] : n[k - b];

        List<bool[]> L;
        if (totalEqs == 0)
        {
            // No arrows in range: the limit is the whole space.
            L = new List<bool[]>(totalVars);
            for (int i = 0; i < totalVars; i++) { var e = new bool[totalVars]; e[i] = true; L.Add(e); }
        }
        else
        {
            bool[,] E = new bool[totalEqs, totalVars];
            int eqIdx = 0;
            for (int k = b; k < d; k++)
            {
                int t = k - b;
                int nk = n[t], nk1 = n[t + 1];
                bool[,] Mk = M[k];
                if (isForward[k])
                {
                    // Forward V_k -> V_{k+1}:  M_k v_k - v_{k+1} = 0
                    for (int row = 0; row < nk1; row++)
                    {
                        for (int col = 0; col < nk; col++) E[eqIdx + row, off[t] + col] = Mk[row, col];
                        E[eqIdx + row, off[t + 1] + row] = true;
                    }
                    eqIdx += nk1;
                }
                else
                {
                    // Backward V_{k+1} -> V_k:  M_k v_{k+1} - v_k = 0
                    for (int row = 0; row < nk; row++)
                    {
                        for (int col = 0; col < nk1; col++) E[eqIdx + row, off[t + 1] + col] = Mk[row, col];
                        E[eqIdx + row, off[t] + row] = true;
                    }
                    eqIdx += nk;
                }
            }
            L = Z2LinearAlgebra.Nullspace(E);
        }

        // Colimit relations R in D: each arrow identifies a source generator with its image.
        var R = new List<bool[]>();
        for (int k = b; k < d; k++)
        {
            int t = k - b;
            int nk = n[t], nk1 = n[t + 1];
            bool[,] Mk = M[k];
            if (isForward[k])
            {
                for (int j = 0; j < nk; j++)
                {
                    var r = new bool[totalVars];
                    r[off[t] + j] = true;
                    for (int row = 0; row < nk1; row++) if (Mk[row, j]) r[off[t + 1] + row] = true;
                    R.Add(r);
                }
            }
            else
            {
                for (int j = 0; j < nk1; j++)
                {
                    var r = new bool[totalVars];
                    r[off[t + 1] + j] = true;
                    for (int row = 0; row < nk; row++) if (Mk[row, j]) r[off[t] + row] = true;
                    R.Add(r);
                }
            }
        }

        // Generalized rank = rank of the canonical map lim -> colim. That map equals
        // (project a section to its first block V_b) then include into the colimit D/R.
        // rank(image in D/R) = rank([images | R]) - rank(R).
        int nb = n[0];
        var combined = new List<bool[]>(L.Count + R.Count);
        foreach (var lv in L)
        {
            var img = new bool[totalVars];
            for (int i = 0; i < nb; i++) img[i] = lv[i];
            combined.Add(img);
        }
        combined.AddRange(R);

        return RankOfRows(combined, totalVars) - RankOfRows(R, totalVars);
    }

    internal static int RankOfRows(List<bool[]> rows, int width)
    {
        if (rows.Count == 0) return 0;
        bool[,] m = new bool[rows.Count, width];
        for (int i = 0; i < rows.Count; i++)
            for (int j = 0; j < width; j++) m[i, j] = rows[i][j];
        return Z2LinearAlgebra.Rank(m);
    }
}
