using System;
using System.Buffers.Binary;
using System.Collections.Generic;
using System.Security.Cryptography;

namespace Maths.Regression;

/// <summary>
/// The evidence seam (E16, CBMA's evidence from the engines' fits): a model's log marginal likelihood ln p(y) of the response as given, in
/// its own units with any internal transform's Jacobian taken in, the estimate's Monte Carlo error (0 where exact), the engine's own flags
/// on it, what the value depends on beyond the data (<see cref="Priors"/>: an engine whose priors read the training response says so
/// there), and the key of what it was read on (<see cref="KeyOf"/>: the response as the engine read it, and its weights). An engine builds it from its fit, as it builds
/// <see cref="CandidateLogDensity"/>; the conformal layer's model averaging reads it (<see cref="Conformal.ModelAveraging"/>), which
/// compares evidences only on one response.
/// </summary>
public sealed record ModelEvidence
{
    /// <param name="logEvidence">ln p(y), finite.</param>
    /// <param name="monteCarloError">Its Monte Carlo error, finite and nonnegative; 0 where exact.</param>
    /// <param name="responseKey">The key of the response it was read on (<see cref="KeyOf"/>).</param>
    /// <param name="flags">The engine's flags on the estimate, each a short reason; none by default.</param>
    /// <param name="priors">What the value depends on beyond the data, as the engine describes it; empty by default.</param>
    public ModelEvidence(double logEvidence, double monteCarloError, string responseKey, IReadOnlyList<string>? flags = null, string? priors = null)
    {
        if (!double.IsFinite(logEvidence)) throw new ArgumentOutOfRangeException(nameof(logEvidence), "A log evidence is finite.");
        if (!(monteCarloError >= 0.0 && double.IsFinite(monteCarloError)))
            throw new ArgumentOutOfRangeException(nameof(monteCarloError), "A Monte Carlo error is finite and nonnegative.");
        ArgumentException.ThrowIfNullOrEmpty(responseKey);
        LogEvidence = logEvidence;
        MonteCarloError = monteCarloError;
        ResponseKey = responseKey;
        Flags = flags is null ? [] : [.. flags];
        Priors = priors ?? "";
    }

    /// <summary>ln p(y).</summary>
    public double LogEvidence { get; }

    /// <summary>Its Monte Carlo error; 0 where exact.</summary>
    public double MonteCarloError { get; }

    /// <summary>The key of the response it was read on.</summary>
    public string ResponseKey { get; }

    /// <summary>The engine's flags on the estimate.</summary>
    public IReadOnlyList<string> Flags { get; }

    /// <summary>What the value depends on beyond the data, as the engine describes it.</summary>
    public string Priors { get; }

    /// <summary>The key of what an evidence was read on: SHA-256, as lowercase hex, of the response's values as the engine read them, in
    /// order, then, where the fit weighs its observations, a marker and the weights, each value's bits little-endian, then
    /// <paramref name="reading"/>, how the weights are read (precision or replication, say), as UTF-8. Two evidences share a key only where
    /// they read the same values, weighed alike; a curve's or several outputs' values are keyed as the engine lays them out.</summary>
    public static string KeyOf(IReadOnlyList<double> response, IReadOnlyList<double>? weights = null, string? reading = null)
    {
        ArgumentNullException.ThrowIfNull(response);
        byte[] text = reading is null ? [] : System.Text.Encoding.UTF8.GetBytes(reading);
        int count = response.Count + (weights is null ? 0 : weights.Count + 1);
        var bytes = new byte[8 * count + (reading is null ? 0 : 8 + text.Length)];
        for (int i = 0; i < response.Count; i++) BinaryPrimitives.WriteInt64LittleEndian(bytes.AsSpan(8 * i), BitConverter.DoubleToInt64Bits(response[i]));
        if (weights is not null)
        {
            BinaryPrimitives.WriteInt64LittleEndian(bytes.AsSpan(8 * response.Count), long.MinValue + weights.Count);
            for (int i = 0; i < weights.Count; i++)
                BinaryPrimitives.WriteInt64LittleEndian(bytes.AsSpan(8 * (response.Count + 1 + i)), BitConverter.DoubleToInt64Bits(weights[i]));
        }
        if (reading is not null)
        {
            BinaryPrimitives.WriteInt64LittleEndian(bytes.AsSpan(8 * count), long.MaxValue - text.Length);
            text.CopyTo(bytes.AsSpan(8 * count + 8));
        }
        return Convert.ToHexString(SHA256.HashData(bytes)).ToLowerInvariant();
    }
}
