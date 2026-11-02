using System;
using System.Globalization;
using System.Linq;

namespace UserRepl.Commands;

/// <summary>
/// Parses the <c>--temperatures</c> argument value into a concrete
/// <c>double[]</c> grid. Supported forms:
/// <list type="bullet">
///   <item><c>linspace:Tmin,Tmax,N</c> — evenly-spaced grid of N points</item>
///   <item><c>logspace:Tmin,Tmax,N</c> — log-spaced grid of N points (both bounds must be positive)</item>
///   <item><c>lnspace:lnTmin,lnTmax,N</c> — N points evenly spaced in <c>ln T</c> (bounds may be negative);
///         a grid of <c>ln T</c> reaches temperatures below the double range, and only
///         <see cref="ParseGrid"/> returns it — <see cref="Parse"/> is for commands that take a T</item>
///   <item>A bare comma-separated list of doubles, e.g. <c>0.01,0.05,0.1,0.5</c></item>
/// </list>
/// </summary>
/// <remarks>
/// <para>Endpoints are inclusive: <c>linspace:0,1,3</c> yields
/// <c>[0, 0.5, 1]</c>. Single-point grids (<c>linspace:0.1,0.1,1</c>)
/// are allowed; the strategy will treat them as a one-temperature
/// sweep.</para>
///
/// <para>Validation is strict — malformed input throws
/// <see cref="ArgumentException"/> with a message that points at the
/// offending token. Callers should let the CLI's top-level error
/// handler print the message.</para>
/// </remarks>
public static class TemperatureGridSpec
{
    /// <summary>
    /// <see cref="Parse"/>, or the <c>lnspace:</c> grid of <c>ln T</c>; the result says which it holds.
    /// </summary>
    public static TemperatureGrid ParseGrid(string spec)
    {
        if (string.IsNullOrWhiteSpace(spec))
            throw new ArgumentException("Temperature spec must be non-empty.", nameof(spec));

        string trimmed = spec.Trim();
        if (trimmed.StartsWith(LogSpacingPrefix, StringComparison.OrdinalIgnoreCase))
            return new TemperatureGrid(ParseSpacing(trimmed.Substring(LogSpacingPrefix.Length), "lnspace"), IsLog: true);
        return new TemperatureGrid(Parse(trimmed), IsLog: false);
    }

    private const string LogSpacingPrefix = "lnspace:";

    public static double[] Parse(string spec)
    {
        if (string.IsNullOrWhiteSpace(spec))
            throw new ArgumentException("Temperature spec must be non-empty.", nameof(spec));

        string trimmed = spec.Trim();
        int colonIdx = trimmed.IndexOf(':');
        if (colonIdx > 0)
        {
            string mode = trimmed.Substring(0, colonIdx).Trim().ToLowerInvariant();
            string body = trimmed.Substring(colonIdx + 1);
            return mode switch
            {
                "linspace" => ParseSpacing(body, mode),
                "logspace" => ParseSpacing(body, mode),
                "lnspace" => throw new ArgumentException(
                    "A 'lnspace' grid holds ln T, which this command does not take; it takes a temperature grid " +
                    "('linspace:Tmin,Tmax,N', 'logspace:Tmin,Tmax,N' or an explicit list)."),
                _ => throw new ArgumentException(
                    $"Unknown temperature spec mode '{mode}'. " +
                    "Use 'auto', 'linspace:Tmin,Tmax,N', 'logspace:Tmin,Tmax,N', or an explicit comma-separated list."),
            };
        }

        return ParseExplicit(trimmed);
    }

    /// <param name="mode"><c>linspace</c>, <c>logspace</c> or <c>lnspace</c>; the last evenly spaces ln T, so its bounds are named ln T.</param>
    private static double[] ParseSpacing(string body, string mode)
    {
        bool logScale = mode == "logspace";
        string lo = mode == "lnspace" ? "lnTmin" : "Tmin";
        string hi = mode == "lnspace" ? "lnTmax" : "Tmax";

        var parts = body.Split(',', StringSplitOptions.RemoveEmptyEntries)
            .Select(p => p.Trim()).ToArray();
        if (parts.Length != 3)
            throw new ArgumentException($"{mode} requires three values: {lo},{hi},N.");

        double tMin = ParseDouble(parts[0], lo);
        double tMax = ParseDouble(parts[1], hi);
        int n       = ParseInt(parts[2], "N");

        if (n < 1) throw new ArgumentException($"N ({n}) must be at least 1.");
        if (tMax < tMin) throw new ArgumentException($"{hi} ({tMax}) must be >= {lo} ({tMin}).");

        if (n == 1) return new[] { tMin };

        var grid = new double[n];
        if (logScale)
        {
            if (tMin <= 0.0 || tMax <= 0.0)
                throw new ArgumentException(
                    $"logspace requires positive bounds; got Tmin={tMin}, Tmax={tMax}.");

            double logMin = Math.Log(tMin);
            double logMax = Math.Log(tMax);
            double step = (logMax - logMin) / (n - 1);
            for (int i = 0; i < n; i++)
                grid[i] = Math.Exp(logMin + i * step);
        }
        else
        {
            double step = (tMax - tMin) / (n - 1);
            for (int i = 0; i < n; i++)
                grid[i] = tMin + i * step;
        }

        return grid;
    }

    private static double[] ParseExplicit(string body)
    {
        var parts = body.Split(',', StringSplitOptions.RemoveEmptyEntries);
        if (parts.Length == 0)
            throw new ArgumentException("Explicit temperature list must contain at least one value.");

        var values = new double[parts.Length];
        for (int i = 0; i < parts.Length; i++)
            values[i] = ParseDouble(parts[i].Trim(), $"value[{i}]");
        return values;
    }

    private static double ParseDouble(string text, string fieldName)
    {
        if (!double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out double v) || !double.IsFinite(v))
            throw new ArgumentException($"Cannot parse {fieldName} as a finite double: '{text}'.");
        return v;
    }

    private static int ParseInt(string text, string fieldName)
    {
        if (!int.TryParse(text, NumberStyles.Integer, CultureInfo.InvariantCulture, out int v))
            throw new ArgumentException($"Cannot parse {fieldName} as integer: '{text}'.");
        return v;
    }
}

/// <summary>
/// A parsed temperature grid and what its values are: temperatures T, or their logarithms ln T
/// (<see cref="TemperatureGridSpec.ParseGrid"/>'s <c>lnspace:</c>, which reaches below the double range).
/// </summary>
public sealed record TemperatureGrid(double[] Values, bool IsLog)
{
    /// <summary>ln T of grid point <paramref name="i"/>.</summary>
    public double LogTemperature(int i) => IsLog ? Values[i] : Math.Log(Values[i]);

    /// <summary>T of grid point <paramref name="i"/>; 0 once a log grid's ln T falls below about −745.</summary>
    public double Temperature(int i) => IsLog ? Math.Exp(Values[i]) : Values[i];
}
