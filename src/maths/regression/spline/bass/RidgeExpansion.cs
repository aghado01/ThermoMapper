using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A ridge basis's expansion (E36's slice 3): each structure term, a directional term (θ, t₀), into its block of K column terms
/// b₁, …, b_K over θ's projection, with the knots t₁, …, t_{K+1} the projected training cases above t₀ give
/// (<see cref="DirectionalTermSpace.RidgeKnots"/>), carried by every column (<see cref="RidgeColumn"/>). The design expands the structure it
/// holds and reads it back from the columns, so the chain, the prior and the proposals see structure terms and the fit sees columns. A
/// term's block is computed once and kept beside it, a pure cache. Immutable otherwise, so one serves every chain of a fit.
/// </summary>
public sealed class RidgeExpansion
{
    private readonly DirectionalTermSpace _space;
    private readonly ConditionalWeakTable<BassTerm, StrongBox<BassTerm[]?>> _blocks = new();

    /// <summary>The expansion into <paramref name="columns"/> = K ≥ 2 columns over <paramref name="space"/>'s training cases.</summary>
    public RidgeExpansion(DirectionalTermSpace space, int columns)
    {
        ArgumentNullException.ThrowIfNull(space);
        if (columns < 2) throw new ArgumentOutOfRangeException(nameof(columns), "A ridge block has K ≥ 2 columns; K = 1 is the directional hinge alone.");
        _space = space;
        Columns = columns;
    }

    /// <summary>K, a block's columns.</summary>
    public int Columns { get; }

    /// <summary>The block of structure term <paramref name="unit"/>, a directional term, as K column terms in order ℓ = 1, …, K; null where its
    /// knots leave a column undefined, a term outside the support.</summary>
    public IReadOnlyList<BassTerm>? Expand(BassTerm unit)
    {
        ArgumentNullException.ThrowIfNull(unit);
        if (unit.Directional is not { } hinge) throw new ArgumentException("A ridge basis expands directional terms (θ, t₀).", nameof(unit));
        return _blocks.GetValue(unit, u => new StrongBox<BassTerm[]?>(Block(hinge))).Value;
    }

    private BassTerm[]? Block(DirectionalHinge hinge)
    {
        if (_space.RidgeKnots(hinge, Columns) is not { } knots) return null;
        var block = new BassTerm[Columns];
        for (int l = 1; l <= Columns; l++) block[l - 1] = new BassTerm(ridge: new RidgeColumn(hinge.Direction, knots, l));
        return block;
    }

    /// <summary>The structure terms of <paramref name="columns"/>, blocks of K in order, read back from each block's first column: (θ, t₀).</summary>
    public IReadOnlyList<BassTerm> Units(IReadOnlyList<BassTerm> columns)
    {
        ArgumentNullException.ThrowIfNull(columns);
        if (columns.Count % Columns != 0) throw new ArgumentException($"{columns.Count} columns are not blocks of {Columns}.", nameof(columns));
        var units = new BassTerm[columns.Count / Columns];
        for (int b = 0; b < units.Length; b++)
        {
            RidgeColumn first = columns[b * Columns].Ridge ?? throw new ArgumentException("A ridge basis's columns are ridge columns.", nameof(columns));
            if (first.Column != 1 || first.Columns != Columns) throw new ArgumentException("The columns are not whole blocks in order.", nameof(columns));
            units[b] = new BassTerm(directional: new DirectionalHinge(first.Direction, first.Knots[0]));
        }
        return units;
    }
}
