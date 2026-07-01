using System;
using System.IO;

namespace UserRepl;

/// <summary>
/// Resolver for the project's centralized canonical-dataset store
/// (<c>&lt;project-root&gt;/datasets/</c>). The single home for the real-world
/// benchmark assets (Iris, Landsat, ISOLET, the SPC_N reference impl) plus
/// their provenance and prep scripts — so tests, an interactive REPL, and demo
/// staging all reach them the same way instead of through bin-depth-pinned
/// <c>../../../../../</c> relative paths.
/// </summary>
/// <remarks>
/// The root is located by walking up from the running assembly's base directory
/// until a <c>datasets/</c> folder holding the canonical <c>iris.csv</c> anchor
/// is found — robust to where the exe/test bin sits. Synthetic canonical sets
/// (Bwd1995Toy, BlattHierarchy, EyeTorusToy, …) are generated in code under
/// <c>src/synthetic/</c> and are NOT file assets; see <c>datasets/README.md</c>.
/// </remarks>
public static class Datasets
{
    private const string AnchorFile = "iris.csv";
    private static readonly Lazy<string> RootLazy = new(LocateRoot);

    /// <summary>Absolute path to the centralized <c>datasets/</c> directory.</summary>
    public static string Root => RootLazy.Value;

    /// <summary>Absolute path to a dataset asset by its name relative to the
    /// datasets root, e.g. <c>Datasets.Path("iris.csv")</c> or
    /// <c>Datasets.Path("reference/spc_n/data")</c>.</summary>
    public static string Path(string relativeName)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(relativeName);
        return System.IO.Path.Combine(Root, relativeName);
    }

    private static string LocateRoot()
    {
        for (DirectoryInfo? dir = new(AppContext.BaseDirectory); dir is not null; dir = dir.Parent)
        {
            string candidate = System.IO.Path.Combine(dir.FullName, "datasets");
            if (File.Exists(System.IO.Path.Combine(candidate, AnchorFile)))
                return candidate;
        }

        throw new DirectoryNotFoundException(
            $"Could not locate the canonical datasets/ root (no datasets/{AnchorFile} found walking up " +
            $"from '{AppContext.BaseDirectory}').");
    }
}
