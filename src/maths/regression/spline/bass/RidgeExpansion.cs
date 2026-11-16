using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A ridge basis's expansion (E36's slice 3): each structure term, a directional term (θ, t₀), into its block of K column terms
/// b₁, …, b_K over θ's projection, with the knots t₁, …, t_{K+1} the projected training cases above t₀ give
/// (<see cref="DirectionalTermSpace.RidgeKnots"/>), carried by every column (<see cref="RidgeColumn"/>). The design expands the structure it
/// holds and reads it back from the columns, so the chain, the prior and the proposals see structure terms and the fit sees columns. A
/// term's categorical factors (slice 5) gate every column of its block, and its knots are read over its cells; a term of categorical factors
/// alone is one column, itself, so blocks vary in width and the read-back walks them term by term. A term's block is computed once and kept
/// beside it, a pure cache. Immutable otherwise, so one serves every chain of a fit.
/// </summary>
public sealed class RidgeExpansion
{
    private readonly DirectionalTermSpace _space;
    private readonly ConditionalWeakTable<BassTerm, StrongBox<BassTerm[]?>> _blocks = new();

    /// <summary>The expansion into <paramref name="columns"/> = K ≥ 2 columns over <paramref name="space"/>'s training cases.</summary>
    public RidgeExpansion(DirectionalTermSpace space, int columns, bool mixedWidths = false)
    {
        ArgumentNullException.ThrowIfNull(space);
        if (columns < 2) throw new ArgumentOutOfRangeException(nameof(columns), "A ridge block has K ≥ 2 columns; K = 1 is the directional hinge alone.");
        _space = space;
        Columns = columns;
        MixedWidths = mixedWidths;
    }

    /// <summary>K, a block's columns.</summary>
    public int Columns { get; }

    /// <summary>Whether each directional structure term declares its own width.</summary>
    public bool MixedWidths { get; }

    /// <summary>The block of structure term <paramref name="unit"/>: a directional term (θ, t₀), with any categorical factors gating it, as K
    /// column terms in order ℓ = 1, …, K, each holding those factors; a term of categorical factors alone as itself. Null where the knots leave
    /// a column undefined, a term outside the support.</summary>
    public IReadOnlyList<BassTerm>? Expand(BassTerm unit)
    {
        ArgumentNullException.ThrowIfNull(unit);
        if (unit.Ridge is not null) throw new ArgumentException("A ridge basis expands structure terms, not columns.", nameof(unit));
        if (unit.Directional is null)
        {
            if (unit.Categorical.Count == 0) throw new ArgumentException("A ridge basis expands directional terms (θ, t₀) and categorical ones.", nameof(unit));
            return [unit];
        }
        if (MixedWidths && unit.Directional.Columns == 1) return [unit];
        return _blocks.GetValue(unit, u => new StrongBox<BassTerm[]?>(Block(u))).Value;
    }

    private BassTerm[]? Block(BassTerm unit)
    {
        DirectionalHinge hinge = unit.Directional!;
        if (_space.RidgeKnots(hinge, unit.Categorical, Columns) is not { } knots) return null;
        var block = new BassTerm[Columns];
        for (int l = 1; l <= Columns; l++) block[l - 1] = new BassTerm(categorical: unit.Categorical, ridge: new RidgeColumn(hinge.Direction, knots, l));
        return block;
    }

    /// <summary>The structure terms of <paramref name="columns"/>, whole blocks in order, read back term by term: a block of K from its first
    /// column, (θ, t₀) with its categorical factors; a column without a ridge, a term of categorical factors alone, as itself.</summary>
    public IReadOnlyList<BassTerm> Units(IReadOnlyList<BassTerm> columns)
    {
        ArgumentNullException.ThrowIfNull(columns);
        var units = new List<BassTerm>(columns.Count / Columns + 1);
        for (int j = 0; j < columns.Count;)
        {
            BassTerm column = columns[j];
            if (column.Ridge is not { } first)
            {
                if (column.Directional is not null && (!MixedWidths || column.Directional.Columns != 1)) throw new ArgumentException("A ridge basis's columns are ridge columns or categorical terms.", nameof(columns));
                units.Add(column);
                j++;
                continue;
            }
            if (first.Column != 1 || first.Columns != Columns || j + Columns > columns.Count)
                throw new ArgumentException("The columns are not whole blocks in order.", nameof(columns));
            for (int l = 1; l < Columns; l++)
                if (columns[j + l].Ridge is not { } next || next.Column != l + 1 || next.Columns != Columns
                    || !System.Linq.Enumerable.SequenceEqual(next.Knots, first.Knots)
                    || !System.Linq.Enumerable.SequenceEqual(columns[j + l].Categorical, column.Categorical)
                    || !ReferenceEquals(next.Direction, first.Direction) && !next.Direction.Equals(first.Direction))
                    throw new ArgumentException("The columns are not whole blocks in order.", nameof(columns));
            units.Add(new BassTerm(categorical: column.Categorical, directional: new DirectionalHinge(first.Direction, first.Knots[0], MixedWidths ? Columns : 1)));
            j += Columns;
        }
        return units;
    }

    /// <summary>The columns structure term <paramref name="unit"/> expands into: K for a direction, 1 for categorical factors alone.</summary>
    public int Width(BassTerm unit) => unit.Directional is null ? 1 : MixedWidths ? unit.Directional.Columns : Columns;
}
