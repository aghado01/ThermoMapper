#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// A <see cref="FastZigzag"/> barcode with a compatible cycle for each bar at every level it covers: consecutive levels'
/// cycles are homologous in the larger of the two complexes, and at each level the alive bars' classes are a basis of
/// homology, so this is the interval decomposition itself. Dey and Hou's set of representatives, the object Dey, Hou and
/// Morozov (arXiv:2410.20565) compute with wires and bundles; here it is read from the coned run's <c>V</c>.
/// <para>A bar keeps its death column's chain <c>V = ω·S + W</c>, and a cycle is a view of that chain at a level, so a bar
/// takes one chain and not one cycle per level. Born by an addition, <c>z_ℓ = S + ∂(W's cells added at step ≥ ℓ)</c>. Born by
/// a deletion, <c>z_ℓ = ∂(S's cells deleted at step &lt; ℓ)</c>. The proof is in <c>docs/TDA/ZigzagPersistence.md</c>.</para>
/// </summary>
public sealed class ZigzagCycles
{
    readonly record struct BarChain(int Dimension, double Death, bool BornByAddition, int[] Added, int[] Coned);

    readonly Dictionary<int, BarChain> _chains;
    readonly IReadOnlyList<int[]> _boundary;
    readonly IReadOnlyList<int> _addStep;
    readonly IReadOnlyList<int> _deleteStep;
    readonly IReadOnlyList<int> _cell;

    /// <summary>The bars. Each carries its <see cref="Bar.Cycle"/>, this type's cycle at the bar's first level.</summary>
    public Barcode Barcode { get; }

    internal ZigzagCycles(
        Barcode barcode,
        IReadOnlyList<int[]> boundary,
        IReadOnlyList<int> addStep,
        IReadOnlyList<int> deleteStep,
        IReadOnlyList<int> cell,
        Dictionary<int, (bool BornByAddition, int[] Added, int[] Coned)> chainsByBirthStep)
    {
        Barcode = barcode;
        _boundary = boundary;
        _addStep = addStep;
        _deleteStep = deleteStep;
        _cell = cell;

        _chains = new Dictionary<int, BarChain>(barcode.Bars.Count);
        foreach (Bar bar in barcode.Bars)
            if (chainsByBirthStep.TryGetValue((int)bar.Birth, out var chain))
                _chains[(int)bar.Birth] = new BarChain(bar.Dimension, bar.Death, chain.BornByAddition, chain.Added, chain.Coned);
    }

    /// <summary>
    /// The bar's cycle at a level, as sorted global cell ids, all of them cells of the complex there. Level ℓ is the complex
    /// after ℓ steps, and a bar covers levels <c>Birth + 1</c> through <c>Death</c>.
    /// </summary>
    /// <param name="bar">A bar of <see cref="Barcode"/>. A step begins at most one bar, so its birth step names it.</param>
    /// <param name="level">A level the bar covers.</param>
    public int[] CycleAt(in Bar bar, int level)
    {
        if (!_chains.TryGetValue((int)bar.Birth, out BarChain chain) || chain.Dimension != bar.Dimension || chain.Death != bar.Death)
            throw new ArgumentException($"[{bar.Birth}, {bar.Death}] in degree {bar.Dimension} is not a bar of this barcode.", nameof(bar));
        if (level < bar.Birth + 1 || level > bar.Death)
            throw new ArgumentOutOfRangeException(nameof(level), $"Level {level} is outside the bar's levels {bar.Birth + 1} through {bar.Death}.");

        var incarnations = new HashSet<int>();
        void Toggle(int incarnation)
        {
            if (!incarnations.Remove(incarnation)) incarnations.Add(incarnation);
        }

        if (chain.BornByAddition)
        {
            foreach (int coned in chain.Coned)
                Toggle(coned);
            foreach (int added in chain.Added)
                if (_addStep[added] >= level)
                    foreach (int face in _boundary[added])
                        Toggle(face);
        }
        else
        {
            foreach (int coned in chain.Coned)
                if (_deleteStep[coned] < level)
                    foreach (int face in _boundary[coned])
                        Toggle(face);
        }

        // Each is alive at this level, so at most one incarnation of a cell is here and the ids are distinct.
        var cycle = new int[incarnations.Count];
        int next = 0;
        foreach (int incarnation in incarnations)
            cycle[next++] = _cell[incarnation];
        Array.Sort(cycle);
        return cycle;
    }
}
