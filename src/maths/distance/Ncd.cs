using System;
using System.IO;
using System.IO.Compression;

namespace Maths.Distance
{
    /// <summary>
    /// Compressor selection for NCD. All options are System.IO.Compression
    /// streaming compressors — no external dependencies.
    /// Brotli is the default; designed for text and ships with a static
    /// English/HTML/CSS dictionary that helps short-input compression.
    /// Deflate / GZip / ZLib are the calibration-clean alternatives
    /// (no static dictionary), differing only in stream framing overhead.
    /// </summary>
    public enum NcdCompressor { Brotli, Deflate, Gzip }

    public static class NormalizedCompressionDistance
    {
        /// <summary>
        /// Pairwise NCD between two byte payloads.
        /// Public entry point for non-graph consumers (e.g. per-segment text
        /// scoring in downstream tools). Caller does not need to construct
        /// an SpcBatchRequest.
        /// Other metrics in this library keep their pairwise distance functions
        /// private; NCD exposes its primitive because direct pairwise scoring
        /// is its primary external use case.
        /// </summary>
        public static double Distance(
            byte[] a, byte[] b, NcdCompressor compressor = NcdCompressor.Brotli)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            int ca = CompressedSize(a, compressor);
            int cb = CompressedSize(b, compressor);
            int cab = CompressedSizeJoint(a, b, compressor);
            return NcdRatio(ca, cb, cab);
        }

        private static double NcdRatio(int ca, int cb, int cab)
        {
            int min = Math.Min(ca, cb);
            int max = Math.Max(ca, cb);
            if (max == 0) return 0.0;
            return (double)(cab - min) / max;
        }

        private static int CompressedSize(byte[] data, NcdCompressor compressor)
        {
            using var output = new MemoryStream();
            using (Stream stream = CreateCompressor(output, compressor))
            {
                stream.Write(data, 0, data.Length);
            }
            return (int)output.Length;
        }

        // Joint compression via two sequential Write calls — avoids
        // allocating a concatenated byte[]. All four BCL streaming
        // compressors process successive Writes as one continuous input
        // (no per-write flush by default), so the result is identical to
        // C(a·b).
        private static int CompressedSizeJoint(byte[] a, byte[] b, NcdCompressor compressor)
        {
            using var output = new MemoryStream();
            using (Stream stream = CreateCompressor(output, compressor))
            {
                stream.Write(a, 0, a.Length);
                stream.Write(b, 0, b.Length);
            }
            return (int)output.Length;
        }

        private static Stream CreateCompressor(Stream output, NcdCompressor compressor)
        {
            return compressor switch
            {
                NcdCompressor.Brotli => new BrotliStream(output, CompressionLevel.Optimal, leaveOpen: true),
                NcdCompressor.Deflate => new DeflateStream(output, CompressionLevel.Optimal, leaveOpen: true),
                NcdCompressor.Gzip => new GZipStream(output, CompressionLevel.Optimal, leaveOpen: true),
                _ => throw new NotSupportedException($"NcdCompressor {compressor} not implemented.")
            };
        }
    }
}
