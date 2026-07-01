using System;
using System.IO;

namespace Archivory;

/// <summary>
/// Shared durable-write primitive for every Archivory artifact writer
/// (binary <see cref="BinarySerializerBase{T}"/>, tabular, and JSON/JSONL).
/// Writes to <c>{path}.tmp</c>, flushes, then atomically renames into place
/// so a partial or interrupted write never leaves a half-written file at the
/// canonical path.
/// </summary>
/// <remarks>
/// One temp-rename implementation, not one per format. <c>overwrite</c>
/// is <see langword="false"/> for write-once checkpoints (a re-run must not clobber
/// an existing artifact) and <see langword="true"/> for regenerable artifacts
/// (manifests, summaries) that a re-run is expected to replace.
/// </remarks>
public static class ArtifactFile
{
    public static void WriteAtomic(string path, Action<Stream> write, bool overwrite = true)
    {
        if (string.IsNullOrWhiteSpace(path))
            throw new ArgumentException("Path must be provided.", nameof(path));
        if (write is null)
            throw new ArgumentNullException(nameof(write));

        string tmpPath = path + ".tmp";
        using (var fs = File.Create(tmpPath))
        {
            write(fs);
            fs.Flush();
        }

        File.Move(tmpPath, path, overwrite);
    }
}
