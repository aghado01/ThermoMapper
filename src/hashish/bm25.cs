#nullable enable
using System.Collections.Frozen;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Hashish;

/// <summary>
/// Computes BM25-style document-frequency statistics from a text corpus.
/// Thin shim over <see cref="InverseDocumentFrequency.Compute"/> retained for
/// the legacy <see cref="SimHash"/> entry point that expects the (avg length,
/// frozen IDF map) tuple shape. New callers should consume <see cref="IdfModel"/>
/// directly — it carries DF counts, average length, and weights in one object.
/// </summary>
public static class Bm25Stats
{
    /// <summary>
    /// Computes average document length and smoothed IDF weights from a corpus.
    /// Output feeds directly into <see cref="SimHash"/> as an immutable IDF map.
    /// IDF formula: log((N + 1) / (df + 1)) + 1  (Lucene/BM25+ smoothed variant).
    /// </summary>
    /// <param name="documents">Input corpus.</param>
    /// <param name="ignoreCase">Fold tokens to lowercase before counting. Default true.</param>
    /// <returns>
    /// A tuple of <c>AvgDocLength</c> (mean token count) and a frozen IDF map.
    /// </returns>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static (double AvgDocLength, FrozenDictionary<string, double> IdfMap)
        Compute(IReadOnlyList<string> documents, bool ignoreCase = true)
    {
        // Preserve historical Bm25Stats semantics: no Unicode compatibility folding,
        // since pre-consolidation Bm25Stats matched the raw input directly. Callers
        // who want FormKC normalization should use InverseDocumentFrequency.Compute.
        IdfModel model = InverseDocumentFrequency.Compute(
            documents,
            formula: IdfFormula.Smooth,
            ignoreCase: ignoreCase,
            normalizeCompatibility: false,
            minTokenLength: 1);

        return (model.AverageDocumentLength, model.Weights);
    }
}
