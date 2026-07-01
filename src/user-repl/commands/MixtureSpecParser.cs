using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Graphs.Coupling;

namespace UserRepl.Commands;

/// <summary>
/// Parses kernel-mixture flag values (<c>--mixture</c>,
/// <c>--mixture-bandwidth</c>) into the
/// <see cref="MixtureWeights"/> / <see cref="MixtureBandwidth"/> records
/// the graph builder consumes.
/// </summary>
/// <remarks>
/// <para><b>Format.</b> Comma-separated <c>name=value</c> pairs. Names
/// are case-insensitive; recognized names are <c>gauss</c>/<c>gaussian</c>,
/// <c>cauchy</c>, <c>laplace</c>/<c>laplacian</c>. Missing components
/// default to zero (weight) or NaN (bandwidth — see
/// <see cref="ParseBandwidth"/>).</para>
///
/// <para><b>Examples.</b> <c>"gauss=0.5,cauchy=0.3,laplace=0.2"</c> →
/// <see cref="MixtureWeights"/>(0.5, 0.3, 0.2). The graph builder
/// normalizes weights internally, so callers don't have to.</para>
/// </remarks>
public static class MixtureSpecParser
{
    public static MixtureWeights ParseWeights(string spec)
    {
        var dict = ParseKeyValuePairs(spec);
        return new MixtureWeights(
            Gaussian:  Get(dict, "gauss", "gaussian", "g"),
            Cauchy:    Get(dict, "cauchy", "c"),
            Laplacian: Get(dict, "laplace", "laplacian", "l"));
    }

    /// <summary>
    /// Parse a per-kernel bandwidth spec. Components not listed default
    /// to <c>0.0</c> — the graph builder treats zero as "auto-estimate
    /// this component from nearest-neighbor distance statistics."
    /// </summary>
    public static MixtureBandwidth ParseBandwidth(string spec)
    {
        var dict = ParseKeyValuePairs(spec);
        return new MixtureBandwidth(
            Gaussian:  Get(dict, "gauss", "gaussian", "g"),
            Cauchy:    Get(dict, "cauchy", "c"),
            Laplacian: Get(dict, "laplace", "laplacian", "l"));
    }

    private static Dictionary<string, double> ParseKeyValuePairs(string spec)
    {
        if (string.IsNullOrWhiteSpace(spec))
            throw new ArgumentException("Mixture spec must be non-empty.", nameof(spec));

        var dict = new Dictionary<string, double>(StringComparer.OrdinalIgnoreCase);
        foreach (string raw in spec.Split(',', StringSplitOptions.RemoveEmptyEntries))
        {
            string token = raw.Trim();
            int eq = token.IndexOf('=');
            if (eq <= 0)
                throw new ArgumentException(
                    $"Mixture token '{token}' must be in the form name=value.");

            string key = token.Substring(0, eq).Trim();
            string val = token.Substring(eq + 1).Trim();

            if (!double.TryParse(val, NumberStyles.Float, CultureInfo.InvariantCulture, out double v))
                throw new ArgumentException(
                    $"Cannot parse mixture value '{val}' for component '{key}' as double.");

            dict[key] = v;
        }
        return dict;
    }

    private static double Get(Dictionary<string, double> dict, params string[] aliases)
    {
        foreach (string a in aliases)
            if (dict.TryGetValue(a, out double v))
                return v;
        return 0.0;
    }
}
