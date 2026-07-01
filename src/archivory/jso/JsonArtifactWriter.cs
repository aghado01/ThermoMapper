using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json;
using System.Text.Json.Nodes;
using System.Text.Json.Serialization.Metadata;

namespace Archivory.Jso;

/// <summary>
/// Generic, schema-agnostic JSON / JSONL artifact writer. Knows nothing about
/// any specific artifact type — a new config or payload type needs zero changes
/// here; it is just another type <c>T</c> resolved through the
/// supplied type-info resolver.
/// </summary>
/// <remarks>
/// <para>Patterns borrowed from the <c>jso-jackson</c> PowerShell engine:
/// the write path never hand-rolls JSON strings or escaping (always
/// <see cref="JsonSerializer"/> / <see cref="Utf8JsonWriter"/>); JSONL is a
/// first-class <i>streaming</i> mode (one compact record per line) rather than
/// "JSON in a loop"; and <see cref="WriteNode"/> is a format-preserving
/// passthrough with no POCO round-trip.</para>
///
/// <para>Construct once per (resolver) and reuse. The instance holds an
/// indented profile for single documents and a compact profile for JSONL
/// records, both derived from the same <see cref="JsonArtifactConventions"/>.</para>
/// </remarks>
public sealed class JsonArtifactWriter
{
    private readonly JsonSerializerOptions _document;   // indented — one structured object
    private readonly JsonSerializerOptions _record;     // compact  — one JSONL line

    /// <summary>Build from the canonical conventions over a domain resolver.</summary>
    public JsonArtifactWriter(IJsonTypeInfoResolver? typeInfoResolver = null)
    {
        _document = JsonArtifactConventions.Create(typeInfoResolver, indented: true);
        _record   = JsonArtifactConventions.Create(typeInfoResolver, indented: false);
    }

    /// <summary>Build from an explicit document-options bundle. The compact
    /// JSONL profile is derived by turning indentation off; everything else
    /// (resolver, naming, number handling) is preserved.</summary>
    public JsonArtifactWriter(JsonSerializerOptions documentOptions)
    {
        _document = documentOptions ?? throw new ArgumentNullException(nameof(documentOptions));
        _record   = new JsonSerializerOptions(documentOptions) { WriteIndented = false };
    }

    // ── Single structured document (manifest, summary, config) ──────────────

    public void WriteDocument<T>(T value, Stream stream) =>
        JsonSerializer.Serialize(stream, value, _document);

    /// <summary>Atomic write — never leaves a half-written file at <paramref name="path"/>.</summary>
    public void WriteDocumentToFile<T>(T value, string path, bool overwrite = true) =>
        ArtifactFile.WriteAtomic(path, stream => WriteDocument(value, stream), overwrite);

    // ── JSONL stream — one compact record per line ──────────────────────────

    /// <summary>Append one compact record + newline to an open stream. Use for
    /// open-ended JSONL sinks (per-record results, <c>errors.jsonl</c>, frame
    /// streams) where records arrive over time.</summary>
    public void AppendRecord<T>(T value, Stream stream)
    {
        JsonSerializer.Serialize(stream, value, _record);
        stream.WriteByte((byte)'\n');
    }

    /// <summary>Stream a sequence to a JSONL file atomically. Returns the path
    /// and the number of records written.</summary>
    public JsonlWriteResult WriteRecords<T>(IEnumerable<T> records, string path, bool overwrite = true)
    {
        if (records is null) throw new ArgumentNullException(nameof(records));

        long count = 0;
        ArtifactFile.WriteAtomic(path, stream =>
        {
            foreach (T record in records)
            {
                AppendRecord(record, stream);
                count++;
            }
        }, overwrite);

        return new JsonlWriteResult(path, count);
    }

    // ── Format-preserving passthrough (no POCO round-trip) ──────────────────

    /// <summary>Write an existing <see cref="JsonNode"/> verbatim — property
    /// names and structure are preserved exactly; only indentation is controlled.
    /// For transforming/relaying existing JSON without a typed model.</summary>
    public void WriteNode(JsonNode node, Stream stream, bool indented = true)
    {
        if (node is null) throw new ArgumentNullException(nameof(node));
        using var writer = new Utf8JsonWriter(stream, new JsonWriterOptions { Indented = indented });
        node.WriteTo(writer);
        writer.Flush();
    }
}

/// <summary>Outcome of a <see cref="JsonArtifactWriter.WriteRecords{T}"/> call.</summary>
public sealed record JsonlWriteResult(string OutputPath, long RecordCount);
