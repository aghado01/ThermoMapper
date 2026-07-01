using System;
using System.IO;

namespace Archivory;

/// <summary>
/// Generic stream-based binary serializer contract. Implementations own
/// a self-describing on-disk format (typically a magic + version header)
/// and define how a value of <typeparamref name="T"/> is written to /
/// read from a stream.
/// </summary>
/// <remarks>
/// <para>The first primitive in the <c>Archivory</c> project — the
/// long-term home for the repository's persistence layer. The current
/// shape is intentionally minimal (stream-write + stream-read +
/// canonical file extension); a fuller blittable-archive primitive
/// (manifest + compression + seek-able partial reads) is planned to
/// land here alongside this interface, with concrete consumers
/// migrating onto it as it matures.</para>
/// <para>The <c>Tabular</c> project is conceptually a specialization
/// within <c>Archivory</c>'s broader mission — CSV is a lossy,
/// human-readable projection — and is expected to fold in here over
/// time. For now the two coexist as siblings while <c>Archivory</c>'s
/// surface stabilizes.</para>
/// </remarks>
/// <typeparam name="T">Domain type being serialized.</typeparam>
public interface IBinarySerializer<T>
{
    /// <summary>Canonical file extension this serializer expects/produces.</summary>
    string DefaultFileExtension { get; }

    void WriteTo(T value, Stream stream);

    T ReadFrom(Stream stream);
}

/// <summary>
/// Base implementation of <see cref="IBinarySerializer{T}"/> that
/// provides the file-write and file-read instance helpers — including
/// atomic write-via-temp-rename. Concrete subclasses implement the
/// per-format <see cref="WriteTo"/> / <see cref="ReadFrom"/> stream
/// methods and the <see cref="DefaultFileExtension"/> property.
/// </summary>
public abstract class BinarySerializerBase<T> : IBinarySerializer<T>
{
    public abstract string DefaultFileExtension { get; }

    public abstract void WriteTo(T value, Stream stream);

    public abstract T ReadFrom(Stream stream);

    /// <summary>
    /// Atomically write <paramref name="value"/> to
    /// <paramref name="path"/>. Writes to <c>{path}.tmp</c> first, then
    /// renames; a partial write never leaves a half-written file at the
    /// canonical path. The rename fails (throws) if the destination
    /// already exists — checkpoint files are treated as write-once.
    /// </summary>
    public void WriteToFile(T value, string path)
    {
        if (value is null)
            throw new ArgumentNullException(nameof(value));

        // Checkpoints are write-once (overwrite: false). Shares the atomic
        // temp-rename primitive with every other Archivory artifact writer.
        ArtifactFile.WriteAtomic(path, stream => WriteTo(value, stream), overwrite: false);
    }

    public T ReadFromFile(string path)
    {
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path must be provided.", nameof(path));
        if (!File.Exists(path))
            throw new FileNotFoundException("Serialized file not found.", path);

        using var fs = File.OpenRead(path);
        return ReadFrom(fs);
    }
}
