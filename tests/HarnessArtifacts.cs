using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text.Json;
using System.Text;

namespace Repo.TestHarness;

internal sealed class ArtifactRun
{
    private static readonly JsonSerializerOptions JsonOptions = new()
    {
        WriteIndented = true,
    };

    public ArtifactRun(string runDirectory, string runName, string manifestPath)
    {
        RunDirectory = runDirectory;
        RunName = runName;
        ManifestPath = manifestPath;
    }

    public string RunDirectory { get; }

    public string RunName { get; }

    public string ManifestPath { get; }

    public string WriteJson<T>(string fileName, T value)
    {
        string path = Path.Combine(RunDirectory, fileName);
        File.WriteAllText(path, JsonSerializer.Serialize(value, JsonOptions));
        return path;
    }

    public string WriteText(string fileName, string content)
    {
        string path = Path.Combine(RunDirectory, fileName);
        File.WriteAllText(path, content);
        return path;
    }

    public string WriteRunJson<T>(string artifactName, T value)
    {
        return WriteJson(HarnessArtifacts.BuildArtifactFileName(RunName, artifactName, "json"), value);
    }

    public string WriteRunText(string artifactName, string content)
    {
        return WriteText(HarnessArtifacts.BuildArtifactFileName(RunName, artifactName, "txt"), content);
    }
}

internal static class HarnessArtifacts
{
    private const string SharedRunDirectoryEnvVar = "PWSHSPC_SHARED_RUN_DIRECTORY";
    private const string SharedManifestPathEnvVar = "PWSHSPC_SHARED_MANIFEST_PATH";

    public static ArtifactRun Create(
        string runKind,
        string suiteName,
        string runName,
        IReadOnlyDictionary<string, object?>? metadata = null)
    {
        string repositoryRoot = FindRepositoryRoot();
        string suiteDirectory = Path.Combine(
            repositoryRoot,
            "artifacts",
            SanitizePathSegment(runKind),
            SanitizePathSegment(suiteName));
        Directory.CreateDirectory(suiteDirectory);

        string? sharedRunDirectory = Environment.GetEnvironmentVariable(SharedRunDirectoryEnvVar);
        if (!string.IsNullOrWhiteSpace(sharedRunDirectory))
        {
            string sharedDirectory = Path.GetFullPath(sharedRunDirectory);
            string? sharedManifestPath = Environment.GetEnvironmentVariable(SharedManifestPathEnvVar);
            bool writeManifest = string.IsNullOrWhiteSpace(sharedManifestPath);
            string manifestPath = writeManifest
                ? Path.Combine(sharedDirectory, BuildArtifactFileName(runName, "manifest", "json"))
                : Path.GetFullPath(sharedManifestPath!);

            return CreateCore(runKind, suiteName, runName, sharedDirectory, repositoryRoot, metadata, manifestPath, writeManifest, sharedRunDirectory, sharedManifestPath);
        }

        string runDirectory = Path.Combine(
            suiteDirectory,
            $"{DateTime.Now:yyyyMMdd_HHmmssfff}__{SanitizePathSegment(runName)}__{Guid.NewGuid():N}");

        string perTestManifestPath = Path.Combine(runDirectory, BuildArtifactFileName(runName, "manifest", "json"));
        return CreateCore(runKind, suiteName, runName, runDirectory, repositoryRoot, metadata, perTestManifestPath, writeManifest: true, sharedRunDirectory: null, sharedManifestPath: null);
    }

    public static ArtifactRun Attach(
        string runKind,
        string suiteName,
        string runName,
        string runDirectory,
        IReadOnlyDictionary<string, object?>? metadata = null)
    {
        string repositoryRoot = FindRepositoryRoot();
        string fullRunDirectory = Path.GetFullPath(runDirectory);
        string manifestPath = Path.Combine(fullRunDirectory, BuildArtifactFileName(runName, "manifest", "json"));
        return CreateCore(runKind, suiteName, runName, fullRunDirectory, repositoryRoot, metadata, manifestPath, writeManifest: true, sharedRunDirectory: null, sharedManifestPath: null);
    }

    private static ArtifactRun CreateCore(
        string runKind,
        string suiteName,
        string runName,
        string runDirectory,
        string repositoryRoot,
        IReadOnlyDictionary<string, object?>? metadata,
        string manifestPath,
        bool writeManifest,
        string? sharedRunDirectory,
        string? sharedManifestPath)
    {
        Directory.CreateDirectory(runDirectory);

        DateTime createdAt = DateTime.Now;
        string runStamp = createdAt.ToString("yyyyMMdd_HHmmssfff", CultureInfo.InvariantCulture);
        Guid runId = Guid.NewGuid();

        var metadataCopy = metadata is null
            ? new Dictionary<string, object?>()
            : new Dictionary<string, object?>(metadata);

        if (!metadataCopy.ContainsKey("RunId"))
            metadataCopy["RunId"] = runId;
        if (!metadataCopy.ContainsKey("RunStamp"))
            metadataCopy["RunStamp"] = runStamp;
        if (!metadataCopy.ContainsKey("SharedRunDirectory") && !string.IsNullOrWhiteSpace(sharedRunDirectory))
            metadataCopy["SharedRunDirectory"] = Path.GetFullPath(sharedRunDirectory);
        if (!metadataCopy.ContainsKey("SharedManifestPath") && !string.IsNullOrWhiteSpace(sharedManifestPath))
            metadataCopy["SharedManifestPath"] = Path.GetFullPath(sharedManifestPath!);

        var manifest = new HarnessRunManifest(
            RunId: runId,
            RunStamp: runStamp,
            RunKind: runKind,
            SuiteName: suiteName,
            RunName: runName,
            CreatedAtLocal: createdAt,
            MachineName: Environment.MachineName,
            ProcessId: Environment.ProcessId,
            FrameworkDescription: System.Runtime.InteropServices.RuntimeInformation.FrameworkDescription,
            RepositoryRoot: repositoryRoot,
            RunDirectory: runDirectory,
            SharedRunDirectory: string.IsNullOrWhiteSpace(sharedRunDirectory) ? null : Path.GetFullPath(sharedRunDirectory),
            SharedManifestPath: string.IsNullOrWhiteSpace(sharedManifestPath) ? null : Path.GetFullPath(sharedManifestPath!),
            Metadata: metadataCopy);

        if (writeManifest)
        {
            File.WriteAllText(
                manifestPath,
                JsonSerializer.Serialize(
                    manifest,
                    new JsonSerializerOptions
                    {
                        WriteIndented = true,
                    }));
        }

        return new ArtifactRun(runDirectory, runName, manifestPath);
    }

    public static string GetSuiteArtifactDirectory(
        string runKind,
        string suiteName,
        string runStamp)
    {
        string repositoryRoot = FindRepositoryRoot();
        return Path.Combine(
            repositoryRoot,
            "artifacts",
            SanitizePathSegment(runKind),
            SanitizePathSegment(suiteName),
            SanitizePathSegment(runStamp));
    }

    private static string FindRepositoryRoot()
    {
        DirectoryInfo? current = new DirectoryInfo(AppContext.BaseDirectory);
        while (current is not null)
        {
            if (File.Exists(Path.Combine(current.FullName, "Directory.Build.props")) &&
                File.Exists(Path.Combine(current.FullName, "changelog.md")))
            {
                return current.FullName;
            }

            current = current.Parent;
        }

        throw new DirectoryNotFoundException("Could not locate the repository root from the harness output directory.");
    }

    internal static string BuildArtifactFileName(string runName, string artifactName, string extension)
    {
        string runStem = SanitizePathSegment(runName);
        string artifactStem = SanitizePathSegment(artifactName);
        string normalizedExtension = extension.Trim().TrimStart('.');
        return $"{runStem}.{artifactStem}.{normalizedExtension}";
    }

    private static string SanitizePathSegment(string value)
    {
        char[] invalid = Path.GetInvalidFileNameChars();
        var builder = new StringBuilder(value.Length);
        bool lastWasSeparator = false;
        for (int i = 0; i < value.Length; i++)
        {
            char c = value[i];
            bool isInvalid = Array.IndexOf(invalid, c) >= 0;
            bool isSeparator = isInvalid || char.IsWhiteSpace(c) || c == '.';
            if (isSeparator)
            {
                if (!lastWasSeparator)
                {
                    builder.Append('-');
                    lastWasSeparator = true;
                }

                continue;
            }

            builder.Append(c);
            lastWasSeparator = false;
        }

        string sanitized = builder.ToString().Trim('-');
        return string.IsNullOrWhiteSpace(sanitized) ? "artifact" : sanitized;
    }

    private sealed record HarnessRunManifest(
        Guid RunId,
        string RunStamp,
        string RunKind,
        string SuiteName,
        string RunName,
        DateTime CreatedAtLocal,
        string MachineName,
        int ProcessId,
        string FrameworkDescription,
        string RepositoryRoot,
        string RunDirectory,
        string? SharedRunDirectory,
        string? SharedManifestPath,
        Dictionary<string, object?> Metadata);
}

/// <summary>
/// Marks a test class as a named harness fixture that the runner can discover and fan out.
/// The runner derives the test filter from the class name: FullyQualifiedName~ClassName.
/// </summary>
[AttributeUsage(AttributeTargets.Class, Inherited = false)]
public sealed class HarnessFixtureAttribute : Attribute
{
    public HarnessFixtureAttribute(string description)
    {
        Description = description;
    }

    public string Description { get; }

    /// <summary>Override the build configuration. Default: Debug.</summary>
    public string? Configuration { get; init; }

    /// <summary>Extra MSBuild/dotnet-test arguments passed through to the child process.</summary>
    public string[]? CommonArguments { get; init; }
}
