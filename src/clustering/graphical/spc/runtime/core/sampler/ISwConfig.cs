namespace Clustering.Graphical.SPC.Runtime.Core.Sampler;

/// <summary>
/// Compile-time configuration for the Swendsen–Wang sampler. Each
/// <c>static abstract</c> member folds to a per-specialization constant, letting the JIT
/// dead-code-eliminate the accumulation branch whose guard is <c>false</c> — a branch-free
/// hot loop per specialization.
/// </summary>
/// <remarks>
/// <para><b>Currency gates, independently monomorphized.</b> The two per-edge accumulations are the
/// ones that justify a JIT gate — each carries <c>O(|E|)</c> memory and a scatter-write per bond in the
/// hot loop. They are now gated <i>independently</i> (<see cref="Affinities"/> ⟂ <see cref="Alignments"/>):
/// the bond-survival counts that reduce to the <c>Affinities</c> currency and the spin-agreement counts
/// that reduce to the <c>Alignments</c> currency are distinct fields, selectable on their own. Four
/// specializations result (<see cref="NoCurrencies"/> / <see cref="AffinitiesOnly"/> /
/// <see cref="AlignmentsOnly"/> / <see cref="BothCurrencies"/>); the declarative
/// <see cref="AccumulationSpec"/> selects one at the cold run-setup boundary.</para>
///
/// <para><b>Everything else is unconditional or runtime-gated.</b> The scalar moments (FK χ, specific
/// heat, magnetization, cluster-size histogram) are free byproducts of the union-find pass and always
/// populated. The per-node landscapes are an <c>O(N)</c> post-pass — cheap enough to live behind a
/// <i>runtime</i> flag on the sampler, not a generic parameter, avoiding further specialization growth.
/// Only an inner-loop scatter-write earns a compile-time bit here.</para>
/// </remarks>
public interface ISwConfig
{
    /// <summary>
    /// When <c>true</c>, accumulate the per-edge bond-survival counts (<c>BondFormedCount</c>) that
    /// reduce to the <c>Affinities</c> currency. <c>O(|E|)</c> memory + a hot-loop scatter-write.
    /// </summary>
    static abstract bool Affinities { get; }

    /// <summary>
    /// When <c>true</c>, accumulate the per-edge spin-agreement counts (<c>SpinAgreementCount</c>) that
    /// reduce to the <c>Alignments</c> currency. <c>O(|E|)</c> memory + a hot-loop scatter-write.
    /// </summary>
    static abstract bool Alignments { get; }
}

/// <summary>Scalar moments + cluster-size histogram only — no per-edge currency accumulation. The lightweight default.</summary>
public readonly struct NoCurrencies : ISwConfig
{
    public static bool Affinities => false;
    public static bool Alignments => false;
}

/// <summary>Adds the per-edge bond-survival counts (the <c>Affinities</c> currency precursor) only.</summary>
public readonly struct AffinitiesOnly : ISwConfig
{
    public static bool Affinities => true;
    public static bool Alignments => false;
}

/// <summary>Adds the per-edge spin-agreement counts (the <c>Alignments</c> currency precursor) only.</summary>
public readonly struct AlignmentsOnly : ISwConfig
{
    public static bool Affinities => false;
    public static bool Alignments => true;
}

/// <summary>Adds both per-edge currency precursors (<c>Affinities</c> + <c>Alignments</c>).</summary>
public readonly struct BothCurrencies : ISwConfig
{
    public static bool Affinities => true;
    public static bool Alignments => true;
}
