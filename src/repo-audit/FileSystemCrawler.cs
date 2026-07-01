using System;
using System.Collections.Generic;
using System.IO;
using System.Security;

namespace RepoAudit;

public sealed record CrawlFile(string AbsolutePath, long SizeBytes);

public sealed record CrawlNode(
    string NodePath,
    string AbsolutePath,
    int NodeDepth,
    List<CrawlFile> Files);

public sealed record SkippedEntry(string Path, string Reason, string? Detail = null);

public sealed record CrawlerResult(
    string RootPath,
    Dictionary<string, CrawlNode> Graph,
    int DirectoryCount,
    int FileCount,
    IReadOnlyList<SkippedEntry> Skipped);

public sealed class FileSystemCrawler
{
    private readonly string _rootPath;
    private readonly Dictionary<string, CrawlNode> _graph;
    private readonly List<SkippedEntry> _skipped;
    private int _dirCount;
    private int _fileCount;
    private bool _hasRun;

    public FileSystemCrawler(string rootPath)
    {
        string resolved = Path.GetFullPath(rootPath).Replace('\\', '/');
        if (!resolved.EndsWith('/')) resolved += '/';

        if (!Directory.Exists(resolved))
            throw new DirectoryNotFoundException($"Root path not found: {rootPath}");

        _rootPath = resolved;
        _graph = new Dictionary<string, CrawlNode>(StringComparer.Ordinal);
        _skipped = new List<SkippedEntry>();
    }

    public CrawlerResult Invoke()
    {
        if (_hasRun)
            throw new InvalidOperationException("Crawler already invoked. Create a new instance.");

        _hasRun = true;

        var queue = new Queue<(string Path, int Depth)>();
        string rootDirPath = _rootPath.TrimEnd('/');

        _graph[""] = new CrawlNode("", _rootPath, 0, new List<CrawlFile>());
        _dirCount = 1;
        queue.Enqueue((rootDirPath, 0));

        while (queue.Count > 0)
        {
            var (dir, depth) = queue.Dequeue();
            string nodePath = ToNodePath(dir);

            try
            {
                foreach (string entryRaw in Directory.EnumerateFileSystemEntries(dir))
                {
                    string entry = entryRaw.Replace('\\', '/');

                    FileAttributes attrs;
                    try { attrs = File.GetAttributes(entry); }
                    catch (Exception ex)
                    {
                        _skipped.Add(new SkippedEntry(entry, "AttributeReadFailed", ex.GetType().Name));
                        continue;
                    }

                    if (attrs.HasFlag(FileAttributes.Directory))
                    {
                        if (attrs.HasFlag(FileAttributes.ReparsePoint))
                        {
                            _skipped.Add(new SkippedEntry(entry, "ReparsePoint"));
                            continue;
                        }

                        string childNodePath = ToNodePath(entry);
                        _graph[childNodePath] = new CrawlNode(childNodePath, entry + "/", depth + 1, new List<CrawlFile>());
                        _dirCount++;
                        queue.Enqueue((entry, depth + 1));
                    }
                    else
                    {
                        long len;
                        try { len = new FileInfo(entry).Length; }
                        catch (Exception ex)
                        {
                            _skipped.Add(new SkippedEntry(entry, "FileSizeReadFailed", ex.GetType().Name));
                            continue;
                        }

                        _graph[nodePath].Files.Add(new CrawlFile(entry, len));
                        _fileCount++;
                    }
                }
            }
            catch (UnauthorizedAccessException) { _skipped.Add(new SkippedEntry(dir, "AccessDenied")); }
            catch (PathTooLongException)         { _skipped.Add(new SkippedEntry(dir, "PathTooLong")); }
            catch (IOException ex)               { _skipped.Add(new SkippedEntry(dir, "IOException", ex.Message)); }
            catch (SecurityException)            { _skipped.Add(new SkippedEntry(dir, "SecurityException")); }
            catch (Exception ex)                 { _skipped.Add(new SkippedEntry(dir, ex.GetType().Name, ex.Message)); }
        }

        return new CrawlerResult(_rootPath, _graph, _dirCount, _fileCount, _skipped.AsReadOnly());
    }

    private string ToNodePath(string fwdSlashPath)
    {
        string rootTrimmed = _rootPath.TrimEnd('/');
        if (fwdSlashPath == rootTrimmed) return "";
        string rel = Path.GetRelativePath(rootTrimmed, fwdSlashPath).Replace('\\', '/');
        if (!rel.EndsWith('/')) rel += '/';
        return rel;
    }
}
