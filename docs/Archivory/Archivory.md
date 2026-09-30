# Archivory

> Archivory is the persistence layer for run artifacts. It covers:
>
> - a deterministic run-directory tree;
> - one atomic write protocol;
> - write-once binary checkpoints;
> - JSON and JSONL conventions;
> - tabular projections that write CSV.
>
> **Projects:** `Archivory`, `Archivory.Jso`, `Archivory.Tabular`
>
> **Maturity:** core. The binary serializer is deliberately minimal.

## Overview

The three projects sit at the bottom of the [architecture](./Architecture.md): `Archivory` and `Archivory.Tabular` are in layer 0 with no references, and `Archivory.Jso` (layer 1) references `Archivory`. The SPC stack uses them for checkpoints and CSV projections. The out-of-scope REPL uses them for run directories and manifests.

A run writes into `{base}/{family}/{stamp}/…`. The file names below are illustrative:

```text
artifacts/                     ← baseDirectory
└── spc/                       ← family (sanitized RunIdentity)
    └── 20260923T101530.123Z/  ← RunStamp (UTC, milliseconds)
        ├── manifest.json      ← JsonArtifactWriter.WriteDocumentToFile
        ├── csv/levels.csv     ← TabularProjectionExtensions.WriteToFile
        └── ckpt/state.bin     ← BinarySerializerBase.WriteToFile (write-once)
```

## Key concepts

- **Scope.** An `ArtifactScope` is an immutable, path-only node of the run tree. Creating one never touches the disk. `EnsureDirectory()` is its only side effect.
- **Run identity.** `RunIdentity.Resolve(requested, callerStub)` returns a `RunIdentity(Family, Source, Requested)` record. A non-blank `requested` wins; otherwise the (required) `callerStub` names the family.
  - `Family` is the sanitized folder name: every run of whitespace or invalid filename characters becomes one `_`, leading and trailing `_` are trimmed, an empty result becomes `"run"`, and case is preserved.
  - `Source` is `explicit` or `auto:caller={sanitized stub}`.
  - `Requested` keeps the raw string, or `null` when no explicit name was given.
- **Stamp.** `RunStamp` formats UTC time as `yyyyMMdd'T'HHmmss'.'fff'Z'` with the invariant culture. Sub-millisecond digits are truncated, not rounded. An unspecified-kind `DateTime` is treated as local time.
- **Atomic replace.** A reader of a path sees either the old file or the complete new one, never a partial write.

## Design

Each concern has a single owner:

| Concern | Owner | Notes |
|---|---|---|
| *Where* to write | `ArtifactScope` (`Root`, `Child`, `File`) | `Dir` is absolute, resolved against the working directory at `Root`. Writers accept a path from `scope.File(name)`. |
| *How* to write safely | `ArtifactFile.WriteAtomic(path, Action<Stream>, overwrite = true)` | Format writers pass a stream callback. |
| JSON conventions | `JsonArtifactConventions.Create(typeInfoResolver, indented, snakeCase, allowNamedFloatingPointLiterals)` | Schema-agnostic. Callers supply their own source-generated resolver. |

**Write policies differ by format.**

- Binary checkpoints are **write-once**: `overwrite: false` is fixed.
- JSON documents and JSONL default to overwriting.
- CSV `WriteToFile` always overwrites; there is no write-once option.

**Two tabular models.**

- `TabularBuilder<TSelf>` is a fluent builder that uses the curiously recurring template pattern (CRTP). Subclasses extend it through the protected hooks `OnColumnsDefined` and `OnRowAdded`. It produces a `TabularData`, which writes CSV.
- `ITabularProjection`/`TabularProjection` add a table name. They are built by `TabularProjectionFactory.Create<T>(tableName, columns, items, selectors)` (or the `items.ProjectToTabular(tableName, columns, selectors)` extension), by `CreateIndexed(tableName, rowCount)` for parallel arrays, or by `CreateScalar(tableName)` for one-row summaries. A projection converts to `TabularData` for output.

`Archivory.Tabular` does not reference `Archivory`, so its `WriteToFile` carries its own copy of the temp-file-and-rename logic. That copy hard-codes `overwrite: true`.

## API at a glance

| Area | Entry points |
|---|---|
| Run tree | `ArtifactScope.Root(baseDirectory, family, stamp)`, `.Child(role)`, `.File(name)`, `.EnsureDirectory()`, `.Dir`; `RunIdentity.Resolve(requested, callerStub)`; `RunStamp.Pattern`, `.Now()`, `.For(DateTime)` |
| Atomic write | `ArtifactFile.WriteAtomic(path, write, overwrite = true)` |
| Binary | `IBinarySerializer<T>` {`DefaultFileExtension`, `WriteTo`, `ReadFrom`}; `BinarySerializerBase<T>` {`WriteToFile`, `ReadFromFile`} |
| JSON | `JsonArtifactConventions.Create(…)`; `JsonArtifactWriter(IJsonTypeInfoResolver?)` or `JsonArtifactWriter(JsonSerializerOptions)`; `WriteDocument`, `WriteDocumentToFile`, `AppendRecord`, `WriteRecords` → `JsonlWriteResult(OutputPath, RecordCount)`, `WriteNode` |
| Tabular | `TabularBuilder<TSelf>` (abstract CRTP base) and the sealed `TabularBuilder`, `TabularData` {`WriteCsv`, `ToCsv`}; `ITabularProjection`, `TabularProjection`, `TabularProjectionExtensions` {`ToTabularData`, `WriteCsv`, `ToCsv`, `WriteToFile`}; `TabularProjectionFactory` {`Create<T>`, `CreateIndexed`, `CreateScalar`}; `EnumerableTabularProjectionExtensions.ProjectToTabular`; `IndexedProjectionBuilder`, `ScalarProjectionBuilder` |

## Algorithms

### Atomic write protocol

1. `File.Create(path + ".tmp")`
2. Invoke the write callback.
3. `Flush()`, then dispose.
4. `File.Move(tmp, path, overwrite)`. With `overwrite: false`, an existing target makes the move throw. This is how checkpoints stay write-once.

The protocol guarantees atomic *replacement*, not durability. It calls `Flush()`, not `Flush(true)`, so there is no fsync and a power loss can still lose the new contents. It also:

- does not create parent directories (call `EnsureDirectory` first);
- leaves the `.tmp` file behind if the callback throws or the target already exists;
- uses one fixed `.tmp` name per path, so concurrent writers to the same path collide.

### JSON and JSONL

**`JsonArtifactConventions.Create`.** The defaults are:

- the default resolver when none is given;
- indented output;
- snake_case property names (`SnakeCaseLower`);
- named floating-point literals allowed, so NaN and ∞ round-trip.

**`JsonArtifactWriter` keeps two profiles, one for documents and one for records.** How they are derived depends on the constructor:

- **Resolver constructor:** two independent `Create` calls, one indented and one not. Snake_case and named floats are always on, so opting out requires the options constructor.
- **Options constructor:** the caller's options are used unchanged for documents, and a copy with `WriteIndented = false` is used for records.

JSONL records are written one per line, and each line ends in a single `\n` on every operating system. The document and bulk-record file writers go through the atomic protocol.

Everything else is inherited from System.Text.Json:

- enums serialize as numbers;
- dictionary keys are not snake-cased;
- nulls are written;
- non-ASCII text is escaped.

### CSV

The output is a header row followed by one line per row.

- A cell is quoted iff its text contains `"`, CR, LF or the delimiter. Embedded quotes are doubled, and null becomes an empty cell.
- The delimiter defaults to `,` and the line ending to `Environment.NewLine`.
- `WriteToFile` writes UTF-8 without a BOM.
- Cell text comes from `ToString()` under the **current culture**. On a comma-decimal locale, `0.5` becomes `"0,5"`. It is quoted because it contains the delimiter, so the file stays parseable, but its contents depend on the locale.

`ColumnsFromDictionary` (indexed and scalar builders) sorts keys with `StringComparer.Ordinal` by default, overridable through `keyOrder`, so column order is deterministic.

## Configuration

| Setting | Default | Meaning |
|---|---|---|
| `ArtifactFile.WriteAtomic` `overwrite` | true | false means write-once |
| `BinarySerializerBase.WriteToFile` | overwrite false (fixed) | Checkpoints are write-once |
| `RunStamp.Pattern` | `yyyyMMdd'T'HHmmss'.'fff'Z'` | UTC, milliseconds |
| `JsonArtifactConventions.Create` | default resolver, indented, snake_case, named floats | JSON profile |
| `JsonArtifactWriter.WriteNode` `indented` | true | — |
| CSV delimiter / line ending | `,` / `Environment.NewLine` | — |

## Status and caveats

- **Atomic, not durable.** Without an fsync, a crash right after a write can lose it. Leftover `.tmp` files after a failure are not cleaned up.
- **Scope paths are not sandboxed.** `Child` and `File` pass their arguments to `Path.Combine`, so a rooted path or a `..` component can escape the scope. `RunIdentity` normalizes separators but does not reject `..`, and the set of invalid characters depends on the operating system.
- **Validation depends on the tabular entry point.** The `TabularData` constructor checks only non-blank names and row width. `TabularProjection`, the factory and the indexed and scalar builders also require at least one column and unique names. `TabularBuilder` checks uniqueness as columns are added (and ≥ 1 column in `WithColumns`), but `Build()` on a builder with no columns returns an empty table. `TabularBuilder.Build()` hands out a *live* read-only view of its column list, so clearing the builder later changes the built table.
- **Culture-dependent CSV.** Numbers are formatted in the current culture. Set an invariant culture on the thread for portable files.
- **Deliberately minimal binary serialization.** `IBinarySerializer<T>` is a stream contract. The source plans a fuller archive primitive and a merge of `Tabular` into `Archivory`.

## Source map

| Role | Files (`ThermoMapper/src/archivory/`) |
|---|---|
| Run tree | `ArtifactScope.cs`, `RunIdentity.cs`, `RunStamp.cs` |
| Write primitive and binary | `ArtifactFile.cs`, `BinarySerialization.cs` |
| JSON | `jso/JsonArtifactConventions.cs`, `jso/JsonArtifactWriter.cs` |
| Tabular | `tabular/TabularBuilder.cs`, `tabular/TabularData.cs`, `tabular/TabularProjection.cs` |

**See also:** [SpcAnalysis](./SpcAnalysis.md) · [SyntheticDatasets](./SyntheticDatasets.md) · [DesignPrinciples](./DesignPrinciples.md)
