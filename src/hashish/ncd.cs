#nullable enable
using System;
using System.IO;
using System.IO.Compression;
using System.Text;

namespace Hashish;

public enum CompressionCodec { Brotli, Deflate, Gzip, ZLib }

/// <summary>
/// Standalone Normalized Compression Distance primitive for pairwise byte/text comparison.
/// This is intentionally independent from SpcCore's NCD graph-builder path.
/// </summary>
public static class NormalizedCompressionDistance
{
    public static double Compute(string first, string second, CompressionCodec codec = CompressionCodec.Brotli, Encoding? encoding = null)
    {
        ArgumentNullException.ThrowIfNull(first);
        ArgumentNullException.ThrowIfNull(second);
        encoding ??= Encoding.UTF8;
        return Compute(encoding.GetBytes(first), encoding.GetBytes(second), codec);
    }

    public static double Compute(byte[] first, byte[] second, CompressionCodec codec = CompressionCodec.Brotli)
    {
        ArgumentNullException.ThrowIfNull(first);
        ArgumentNullException.ThrowIfNull(second);

        int ca = CompressedSize(first, codec);
        int cb = CompressedSize(second, codec);
        int cab = CompressedSizeJoint(first, second, codec);
        return Ratio(ca, cb, cab);
    }

    public static int CompressedSize(byte[] data, CompressionCodec codec = CompressionCodec.Brotli)
    {
        ArgumentNullException.ThrowIfNull(data);

        using var output = new MemoryStream();
        using (Stream stream = CreateCompressor(output, codec))
            stream.Write(data, 0, data.Length);

        return (int)output.Length;
    }

    public static int CompressedSizeJoint(byte[] first, byte[] second, CompressionCodec codec = CompressionCodec.Brotli)
    {
        ArgumentNullException.ThrowIfNull(first);
        ArgumentNullException.ThrowIfNull(second);

        using var output = new MemoryStream();
        using (Stream stream = CreateCompressor(output, codec))
        {
            stream.Write(first, 0, first.Length);
            stream.Write(second, 0, second.Length);
        }

        return (int)output.Length;
    }

    public static double Ratio(int compressedFirst, int compressedSecond, int compressedJoint)
    {
        int min = Math.Min(compressedFirst, compressedSecond);
        int max = Math.Max(compressedFirst, compressedSecond);
        return max == 0 ? 0.0 : (double)(compressedJoint - min) / max;
    }

    private static Stream CreateCompressor(Stream output, CompressionCodec codec)
    {
        return codec switch
        {
            CompressionCodec.Brotli => new BrotliStream(output, CompressionLevel.Optimal, leaveOpen: true),
            CompressionCodec.Deflate => new DeflateStream(output, CompressionLevel.Optimal, leaveOpen: true),
            CompressionCodec.Gzip => new GZipStream(output, CompressionLevel.Optimal, leaveOpen: true),
            CompressionCodec.ZLib => new ZLibStream(output, CompressionLevel.Optimal, leaveOpen: true),
            _ => throw new NotSupportedException($"Compression codec {codec} is not supported.")
        };
    }
}
