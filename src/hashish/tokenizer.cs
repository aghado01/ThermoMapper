#nullable enable
using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;
using System.Text.RegularExpressions;

namespace Hashish;

/// <summary>
/// Lightweight text normalization and word tokenization for Hashish primitives.
/// Keeps preprocessing explicit so hashing, sketching, and vectorization paths can
/// share the same token stream without depending on SPC metric code.
/// TODO: this will need to be extended with more sophisticated tokenization and normalization options
/// </summary>
public static class TokenizerPreprocessing
{
    /// <summary>
    /// Single shared word-token regex: one compiled instance for every Hashish primitive
    /// that splits on Unicode word characters. Consumers that want to avoid Match allocation
    /// can call <c>EnumerateMatches</c> on the normalized text directly.
    /// </summary>
    internal static readonly Regex WordRegex = new(
        @"\w+",
        RegexOptions.Compiled | RegexOptions.CultureInvariant | RegexOptions.NonBacktracking
    );

    /// <summary>
    /// Normalizes text with optional Unicode compatibility normalization and case folding.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static string Normalize(
        string text,
        bool ignoreCase = true,
        bool normalizeCompatibility = true,
        bool trim = true)
    {
        ArgumentNullException.ThrowIfNull(text);

        string normalized = normalizeCompatibility
            ? text.Normalize(NormalizationForm.FormKC)
            : text.Normalize();

        if (trim)
            normalized = normalized.Trim();

        return ignoreCase ? normalized.ToLowerInvariant() : normalized;
    }

    /// <summary>
    /// Extracts word tokens from text after normalization. Empty input returns an empty array.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static string[] TokenizeWords(
        string text,
        bool ignoreCase = true,
        bool normalizeCompatibility = true,
        int minTokenLength = 1)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentOutOfRangeException.ThrowIfLessThan(minTokenLength, 1, nameof(minTokenLength));

        if (string.IsNullOrWhiteSpace(text))
            return Array.Empty<string>();

        string normalized = Normalize(text, ignoreCase, normalizeCompatibility);
        var matches = WordRegex.Matches(normalized);
        if (matches.Count == 0)
            return Array.Empty<string>();

        var tokens = new List<string>(matches.Count);
        foreach (Match match in matches)
        {
            if (match.Length >= minTokenLength)
                tokens.Add(match.Value);
        }

        return tokens.Count == 0 ? Array.Empty<string>() : tokens.ToArray();
    }
}
