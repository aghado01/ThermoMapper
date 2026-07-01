using System;
using System.Buffers;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;

namespace Hashish;

/// <summary>
/// MinHash signature for Jaccard similarity estimation via shingling.
/// Stateful instance; construct once with a fixed <c>numHashes</c> and <c>shingleSize</c>,
/// then call <see cref="Compute(string)"/> per document.
/// </summary>
public sealed class MinHash
{
    // FNV-1a constants — 32-bit variant for signature slots.
    private const uint Fnv32OffsetBasis = 2166136261u;
    private const uint Fnv32Prime = 16777619u;

    // Stackalloc threshold for UTF-8 token byte buffer.
    private const int Utf8StackLimit = 256;

    private readonly int _numHashes;
    private readonly int _shingleSize;

    public int SignatureLength => _numHashes;
    public int ShingleSize => _shingleSize;

    /// <param name="numHashes">Number of hash functions (signature length). Must be &gt; 0.</param>
    /// <param name="shingleSize">Character n-gram width. Must be &gt; 0.</param>
    public MinHash(int numHashes = 128, int shingleSize = 3)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(numHashes, 0, nameof(numHashes));
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(shingleSize, 0, nameof(shingleSize));
        _numHashes = numHashes;
        _shingleSize = shingleSize;
    }

    /// <summary>
    /// Computes a MinHash signature for <paramref name="content"/>.
    /// Returns an array of <c>numHashes</c> minimum hash values (all <c>uint.MaxValue</c> if input is too short).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public uint[] Compute(string content)
    {
        var signature = new uint[_numHashes];
        signature.AsSpan().Fill(uint.MaxValue);

        if (string.IsNullOrEmpty(content) || content.Length < _shingleSize)
            return signature;

        var shingles = BuildShingles(content.AsSpan(), _shingleSize);
        if (shingles.Count == 0) return signature;

        foreach (string shingle in shingles)
        {
            for (int i = 0; i < _numHashes; i++)
            {
                uint h = HashWithSeed(shingle.AsSpan(), (uint)i);
                if (h < signature[i]) signature[i] = h;
            }
        }

        return signature;
    }

    /// <summary>
    /// Creates a band/row LSH index whose expected signature length matches this MinHash instance.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public MinHashLshIndex CreateLshIndex(int bands = 32, int rowsPerBand = 4)
        => new MinHashLshIndex(this, bands, rowsPerBand);

    /// <summary>
    /// Estimated Jaccard similarity in [0, 1] between two signatures of equal length.
    /// </summary>
    public static double EstimateJaccard(ReadOnlySpan<uint> sig1, ReadOnlySpan<uint> sig2)
    {
        if (sig1.Length != sig2.Length)
            throw new ArgumentException("Signature lengths must match.");

        int matches = 0;
        for (int i = 0; i < sig1.Length; i++)
            if (sig1[i] == sig2[i]) matches++;

        return (double)matches / sig1.Length;
    }

    /// <summary>Array overload.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double EstimateJaccard(uint[] sig1, uint[] sig2)
        => EstimateJaccard(sig1.AsSpan(), sig2.AsSpan());

    // Build a deduplicated shingle set using AsSpan slices to avoid Substring allocations.
    // Strings are only created for HashSet insertion.
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static HashSet<string> BuildShingles(ReadOnlySpan<char> text, int size)
    {
        var shingles = new HashSet<string>(
            capacity: Math.Max(text.Length - size + 1, 0),
            comparer: StringComparer.Ordinal
        );

        int limit = text.Length - size;
        for (int i = 0; i <= limit; i++)
            shingles.Add(new string(text.Slice(i, size)));

        return shingles;
    }

    // Seeded FNV-1a over UTF-8 bytes.
    // For tokens <= Utf8StackLimit bytes: encode into a stackalloc buffer (zero heap).
    // Seed offsets the initial basis so each hash function is independent.
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static uint HashWithSeed(ReadOnlySpan<char> token, uint seed)
    {
        int maxBytes = Encoding.UTF8.GetMaxByteCount(token.Length);
        uint hash = Fnv32OffsetBasis + seed;

        if (maxBytes <= Utf8StackLimit)
        {
            Span<byte> buf = stackalloc byte[Utf8StackLimit];
            int written = Encoding.UTF8.GetBytes(token, buf);
            ReadOnlySpan<byte> bytes = buf[..written];
            for (int i = 0; i < bytes.Length; i++)
            {
                hash ^= bytes[i];
                hash *= Fnv32Prime;
            }
        }
        else
        {
            byte[] rented = ArrayPool<byte>.Shared.Rent(maxBytes);
            try
            {
                int written = Encoding.UTF8.GetBytes(token, rented);
                ReadOnlySpan<byte> bytes = rented.AsSpan(0, written);
                for (int i = 0; i < bytes.Length; i++)
                {
                    hash ^= bytes[i];
                    hash *= Fnv32Prime;
                }
            }
            finally { ArrayPool<byte>.Shared.Return(rented); }
        }

        return hash;
    }

    // Static convenience overloads mirror the PS API surface.

    /// <summary>Static overload: compute with default parameters.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static uint[] Compute(string content, int numHashes = 128, int shingleSize = 3)
        => new MinHash(numHashes, shingleSize).Compute(content);
}

/// <summary>
/// Band/row locality-sensitive index for MinHash signatures.
/// Returns candidate ids that share at least one exact band with the query.
/// </summary>
public sealed class MinHashLshIndex
{
    private readonly int _bands;
    private readonly int _rowsPerBand;
    private readonly int _signatureLength;
    private readonly Dictionary<BandKey, List<int>> _buckets = new();

    public MinHashLshIndex(MinHash minHash, int bands = 32, int rowsPerBand = 4)
        : this(minHash?.SignatureLength ?? throw new ArgumentNullException(nameof(minHash)), bands, rowsPerBand)
    {
    }

    public MinHashLshIndex(int signatureLength, int bands = 32, int rowsPerBand = 4)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(signatureLength, 0, nameof(signatureLength));
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(bands, 0, nameof(bands));
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(rowsPerBand, 0, nameof(rowsPerBand));

        if (checked(bands * rowsPerBand) != signatureLength)
            throw new ArgumentException("Bands x rowsPerBand must equal the MinHash signature length.");

        _bands = bands;
        _rowsPerBand = rowsPerBand;
        _signatureLength = signatureLength;
    }

    public int Bands => _bands;
    public int RowsPerBand => _rowsPerBand;
    public int SignatureLength => _signatureLength;
    public int BucketCount => _buckets.Count;

    /// <summary>Adds a document id and MinHash signature to the index.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public void Add(int id, ReadOnlySpan<uint> signature)
    {
        ValidateSignature(signature);

        for (int band = 0; band < _bands; band++)
        {
            BandKey key = BuildKey(signature, band);
            if (!_buckets.TryGetValue(key, out List<int>? ids))
            {
                ids = new List<int>();
                _buckets.Add(key, ids);
            }

            ids.Add(id);
        }
    }

    /// <summary>Adds a batch of signatures where array index is used as document id.</summary>
    public void AddRange(IReadOnlyList<uint[]> signatures)
    {
        ArgumentNullException.ThrowIfNull(signatures);
        for (int i = 0; i < signatures.Count; i++)
            Add(i, signatures[i]);
    }

    /// <summary>Returns candidate ids sharing at least one band with <paramref name="signature"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public int[] GetCandidateIds(ReadOnlySpan<uint> signature)
    {
        ValidateSignature(signature);

        var candidates = new HashSet<int>();
        for (int band = 0; band < _bands; band++)
        {
            if (_buckets.TryGetValue(BuildKey(signature, band), out List<int>? ids))
            {
                for (int i = 0; i < ids.Count; i++)
                    candidates.Add(ids[i]);
            }
        }

        if (candidates.Count == 0)
            return Array.Empty<int>();

        var result = new int[candidates.Count];
        candidates.CopyTo(result);
        Array.Sort(result);
        return result;
    }

    public void Clear() => _buckets.Clear();

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void ValidateSignature(ReadOnlySpan<uint> signature)
    {
        if (signature.Length != _signatureLength)
            throw new ArgumentException($"Signature length must be {_signatureLength} for {_bands} bands x {_rowsPerBand} rows.");
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private BandKey BuildKey(ReadOnlySpan<uint> signature, int band)
    {
        int offset = band * _rowsPerBand;
        ulong hash = SeededHash.Fnv1a(signature.Slice(offset, _rowsPerBand), SeededHash.Seed(band));
        return new BandKey(band, hash);
    }

    private readonly struct BandKey : IEquatable<BandKey>
    {
        private readonly int _band;
        private readonly ulong _hash;

        public BandKey(int band, ulong hash)
        {
            _band = band;
            _hash = hash;
        }

        public bool Equals(BandKey other) => _band == other._band && _hash == other._hash;
        public override bool Equals(object? obj) => obj is BandKey other && Equals(other);
        public override int GetHashCode() => HashCode.Combine(_band, _hash);
    }
}
