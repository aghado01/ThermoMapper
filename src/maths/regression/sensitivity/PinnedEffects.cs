using System;
using System.Collections.Generic;
using System.Linq;

namespace Maths.Regression.Sensitivity;

/// <summary>A pinned block: the inputs pinned, ascending, and their values; none for the base block.</summary>
public sealed record PinnedBlock(int[] Inputs, double[] Values);

/// <summary>
/// The pinned design for one effect over p inputs (the sampled-sensitivity brief, decision 5): the subset u, ascending, and a grid per
/// input of u. The effect at a cell g of the grids' product is the ANOVA component f_u(g) = Σ_{v⊆u} (−1)^{|u|−|v|} E[f | x_v = g_v],
/// v = ∅ giving E[f], so one input's is the main effect E[f | x_v = g] − E[f] and two inputs' the two-way effect. Each conditional mean
/// averages f over the base rows with x_v pinned; the blocks are the base itself and every nonempty v ⊆ u at every cell of v's grids.
/// Cells run over the grids with the last input fastest.
/// </summary>
public sealed class PinnedEffectDesign
{
    private readonly int[][] _terms, _signs;

    /// <summary>The design for subset <paramref name="subset"/> of <paramref name="inputs"/> inputs, with one nonempty finite grid per
    /// input of the subset in its order (ascending inputs); refused above <paramref name="budget"/> blocks.</summary>
    public PinnedEffectDesign(int inputs, int[] subset, IReadOnlyList<double[]> grids, int budget = 1_000_000)
    {
        ArgumentNullException.ThrowIfNull(subset);
        ArgumentNullException.ThrowIfNull(grids);
        if (inputs < 1) throw new ArgumentOutOfRangeException(nameof(inputs), "A design needs an input.");
        if (subset.Length == 0) throw new ArgumentException("An effect names an input.", nameof(subset));
        if (subset.Length > 16) throw new ArgumentException("An effect names at most 16 inputs.", nameof(subset));
        for (int k = 0; k < subset.Length; k++)
        {
            if (subset[k] < 0 || subset[k] >= inputs) throw new ArgumentException($"Input {subset[k]} is outside the {inputs} inputs.", nameof(subset));
            if (k > 0 && subset[k] <= subset[k - 1]) throw new ArgumentException("The subset's inputs ascend, without repeats.", nameof(subset));
        }
        if (grids.Count != subset.Length) throw new ArgumentException($"The subset has {subset.Length} inputs and {grids.Count} grids.", nameof(grids));
        for (int k = 0; k < grids.Count; k++)
        {
            if (grids[k] is null || grids[k].Length == 0) throw new ArgumentException($"Grid {k} is empty.", nameof(grids));
            if (grids[k].Any(x => !double.IsFinite(x))) throw new ArgumentException($"Grid {k} holds a value that is not finite.", nameof(grids));
        }
        int q = subset.Length;
        int[] sizes = [.. grids.Select(g => g.Length)];
        long count = 1;
        for (int mask = 1; mask < 1 << q; mask++) count += Product(sizes, mask);
        if (count > budget) throw new ArgumentException($"The effect needs {count} blocks, above the budget of {budget}.", nameof(grids));

        // Blocks: the base, then each nonempty mask over u's positions with each cell of its grids, last position fastest.
        var blocks = new List<PinnedBlock> { new([], []) };
        var first = new int[1 << q];
        for (int mask = 1; mask < 1 << q; mask++)
        {
            first[mask] = blocks.Count;
            int[] positions = [.. Enumerable.Range(0, q).Where(t => (mask >> t & 1) != 0)];
            int cells = (int)Product(sizes, mask);
            for (int c = 0; c < cells; c++)
            {
                var values = new double[positions.Length];
                int rest = c;
                for (int t = positions.Length - 1; t >= 0; t--)
                {
                    values[t] = grids[positions[t]][rest % sizes[positions[t]]];
                    rest /= sizes[positions[t]];
                }
                blocks.Add(new PinnedBlock([.. positions.Select(t => subset[t])], values));
            }
        }

        // Each cell's alternating sum: for each mask, the block at the cell's indices on the mask's positions.
        int full = (int)Product(sizes, (1 << q) - 1);
        _terms = new int[full][];
        _signs = new int[full][];
        var index = new int[q];
        for (int c = 0; c < full; c++)
        {
            int rest = c;
            for (int t = q - 1; t >= 0; t--)
            {
                index[t] = rest % sizes[t];
                rest /= sizes[t];
            }
            _terms[c] = new int[1 << q];
            _signs[c] = new int[1 << q];
            for (int mask = 0; mask < 1 << q; mask++)
            {
                int cell = 0;
                for (int t = 0; t < q; t++)
                    if ((mask >> t & 1) != 0) cell = cell * sizes[t] + index[t];
                _terms[c][mask] = mask == 0 ? 0 : first[mask] + cell;
                _signs[c][mask] = (q - System.Numerics.BitOperations.PopCount((uint)mask)) % 2 == 0 ? 1 : -1;
            }
        }
        Inputs = inputs;
        Subset = [.. subset];
        Grids = Array.AsReadOnly(grids.Select(g => (double[])g.Clone()).ToArray());
        Blocks = blocks.AsReadOnly();
        Cells = full;
    }

    /// <summary>The inputs.</summary>
    public int Inputs { get; }

    /// <summary>The effect's inputs, ascending.</summary>
    public int[] Subset { get; }

    /// <summary>Each input's grid, in the subset's order.</summary>
    public IReadOnlyList<double[]> Grids { get; }

    /// <summary>The blocks: the base first, then the pinned ones.</summary>
    public IReadOnlyList<PinnedBlock> Blocks { get; }

    /// <summary>The grids' cells.</summary>
    public int Cells { get; }

    /// <summary>Block <paramref name="block"/>'s points: <paramref name="basePoints"/> (N × p) with its inputs pinned.</summary>
    public double[,] Points(double[,] basePoints, int block)
    {
        CheckBase(basePoints);
        var points = (double[,])basePoints.Clone();
        PinnedBlock pinned = Blocks[block];
        for (int r = 0; r < points.GetLength(0); r++)
            for (int t = 0; t < pinned.Inputs.Length; t++) points[r, pinned.Inputs[t]] = pinned.Values[t];
        return points;
    }

    internal int CheckBase(double[,] basePoints)
    {
        ArgumentNullException.ThrowIfNull(basePoints);
        if (basePoints.GetLength(1) != Inputs) throw new ArgumentException($"The base points have one column per input, {Inputs}.", nameof(basePoints));
        return basePoints.GetLength(0);
    }

    internal int[] Terms(int cell) => _terms[cell];

    internal int[] Signs(int cell) => _signs[cell];

    private static long Product(int[] sizes, int mask)
    {
        long product = 1;
        for (int t = 0; t < sizes.Length; t++)
            if ((mask >> t & 1) != 0) product *= sizes[t];
        return product;
    }
}

/// <summary>An effect's values per cell, each with its Monte Carlo standard error.</summary>
public sealed record PinnedEffectEstimate(double[] Values, double[] Errors);

/// <summary>
/// An effect of any function by pinning (the sampled-sensitivity brief, decision 5): it accumulates the function's values on the base
/// rows and on each pinned block, and reduces them, row by row, to each cell's ANOVA component. Per row the alternating sum
/// Σ_{v⊆u} (−1)^{|u|−|v|} f(g_v, x_{−v}) is one value with expectation f_u(g) under a product law, so each cell's estimate is a mean of
/// N independent terms, unbiased, its error their sample deviation over √N. The sum's signs add to zero, so the estimate is
/// shift-invariant, and the base rows' variation outside u cancels where f is additive across u.
/// </summary>
public static class PinnedEffects
{
    /// <summary>The effect from the values at <paramref name="design"/>'s blocks, one array of N ≥ 2 values per block in its
    /// order.</summary>
    public static PinnedEffectEstimate Estimate(PinnedEffectDesign design, IReadOnlyList<double[]> blockValues)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(blockValues);
        if (blockValues.Count != design.Blocks.Count)
            throw new ArgumentException($"The design has {design.Blocks.Count} blocks and {blockValues.Count} were given.", nameof(blockValues));
        if (blockValues[0] is null) throw new ArgumentException("Block 0 is null.", nameof(blockValues));
        int n = blockValues[0].Length;
        if (n < 2) throw new ArgumentException("An effect needs two base rows.", nameof(blockValues));
        for (int k = 0; k < blockValues.Count; k++)
        {
            if (blockValues[k] is null || blockValues[k].Length != n) throw new ArgumentException($"Block {k} does not hold {n} values.", nameof(blockValues));
            for (int j = 0; j < n; j++)
                if (!double.IsFinite(blockValues[k][j])) throw new ArgumentException($"The function's value on block {k} at row {j} is not finite.", nameof(blockValues));
        }
        var values = new double[design.Cells];
        var errors = new double[design.Cells];
        var e = new double[n];
        for (int c = 0; c < design.Cells; c++)
        {
            int[] terms = design.Terms(c), signs = design.Signs(c);
            for (int j = 0; j < n; j++)
            {
                double s = 0.0;
                for (int t = 0; t < terms.Length; t++) s += signs[t] * blockValues[terms[t]][j];
                e[j] = s;
            }
            values[c] = PickFreeze.Mean(e);
            errors[c] = PickFreeze.Deviation(e) / Math.Sqrt(n);
        }
        return new PinnedEffectEstimate(values, errors);
    }

    /// <summary>The effect of <paramref name="f"/> at <paramref name="design"/>'s blocks over <paramref name="basePoints"/> (N × p,
    /// rows drawn independently from the product law).</summary>
    public static PinnedEffectEstimate Estimate(PinnedEffectDesign design, Func<ReadOnlySpan<double>, double> f, double[,] basePoints)
    {
        ArgumentNullException.ThrowIfNull(design);
        ArgumentNullException.ThrowIfNull(f);
        int n = design.CheckBase(basePoints), p = design.Inputs;
        var row = new double[p];
        var blockValues = new double[design.Blocks.Count][];
        for (int k = 0; k < blockValues.Length; k++)
        {
            PinnedBlock pinned = design.Blocks[k];
            var values = new double[n];
            for (int j = 0; j < n; j++)
            {
                for (int v = 0; v < p; v++) row[v] = basePoints[j, v];
                for (int t = 0; t < pinned.Inputs.Length; t++) row[pinned.Inputs[t]] = pinned.Values[t];
                values[j] = f(row);
            }
            blockValues[k] = values;
        }
        return Estimate(design, blockValues);
    }
}
