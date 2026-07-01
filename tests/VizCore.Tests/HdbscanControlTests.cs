using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Linq;
using Clustering.Graphical.HdbScan;
using Maths.LinAlg;
using UserRepl;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// HDBSCAN as an <b>independent density-based control</b> for the SPC parity
/// benchmarks. HDBSCAN shares NONE of the SPC/Potts machinery (no temperature
/// sweep, no co-membership currency), so it triangulates: where SPC/lineage
/// recovers a gold-standard structure, an independent method should corroborate
/// it; where the data defeats SPC (ISOLET), an independent method should fail
/// too — which is what tells us the wall is the data, not the pipeline.
/// <para>Also the seed of HDBSCAN's own fact battery (it currently has none
/// akin to the SPC parity suite).</para>
/// <para><b>Opt-out:</b> tagged <c>Category=Control</c> — exclude with
/// <c>dotnet test --filter "Category!=Control"</c>, run alone with
/// <c>--filter "Category=Control"</c>.</para>
/// Full experiment + the SPC-vs-HDBSCAN reading lives in
/// <c>.discussion/issues/spc-parity/isolet-pca-wall.md</c>.
/// </summary>
[Trait("Category", "Control")]
public sealed class HdbscanControlTests
{
    /// <summary>
    /// One HDBSCAN run per benchmark, asserting the control's three claims:
    /// (1) it corroborates the gold standard where SPC works (toy stripes, the
    /// Iris setosa split, the Landsat land-cover), and (2) it CONFIRMS THE ISOLET
    /// WALL — an independent method also fails to separate the 26 letters
    /// unsupervised (≤8 of 26), so parking ISOLET is a data limit, not an SPC
    /// defect. Oracle = the true labels (validation-independent), never SPC's own
    /// output. A diagnostic table is also dumped for the growing battery.
    /// </summary>
    [Fact]
    public void Hdbscan_DensityControl_CorroboratesGoldStandard_AndConfirmsIsoletWall()
    {
        var iris = SpcUserSession.FromCsv(Datasets.Path("iris.csv"), null, false, ',');
        var landsat = SpcUserSession.FromCsv(Datasets.Path("landsat.csv"), null, false, ',');
        var toy = SpcUserDataset.FromSyntheticDataset(Synthetic.Euclidean.Bwd1995Toy.Generate(seed: 42), null);
        var (iso617, isoLabels) = LoadIsoletGz(Datasets.Path("isolet.csv.gz"));
        PcaResult pca = Pca.Compute(iso617, numComponents: 50, center: true, whiten: false);
        Func<double[], double[]> proj = Pca.MakeProjector(pca);
        double[][] isoPca50 = iso617.Select(proj).ToArray();

        // toy/iris use raw features (well-scaled); Landsat needs z-norm because
        // HDBSCAN-euclidean is scale-sensitive and has no built-in normalization
        // (SPC's mean-edge bandwidth absorbs raw scale; HDBSCAN does not — see the
        // enhancement notes). ISOLET features are pre-scaled to [-1,1].
        var toyM = Measure("toy",          toy.Features,         toy.Labels,     mcs: 50, mp: 5);
        var irisM = Measure("iris",        iris.Features,        iris.Labels,    mcs: 10, mp: 5);
        var landsatM = Measure("landsat",  ZNorm(landsat.Features), landsat.Labels, mcs: 50, mp: 5);
        var isoM = Measure("isolet-pca50", isoPca50,             isoLabels,      mcs: 50, mp: 5);

        var dump = new[] { "bench\tnClasses\tbigClusters\tbigPurity\tcovered\tnoise\tclassesCaptured" }
            .Concat(new[] { toyM, irisM, landsatM, isoM }.Select(m =>
                $"{m.Name}\t{m.NClasses}\t{m.BigClusters}\t{m.BigPurity:F3}\t{m.Covered}\t{m.Noise}\t{m.ClassesCaptured}"));
        string outPath = Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "../../../../../artifacts/hdbscan-control.tsv"));
        Directory.CreateDirectory(Path.GetDirectoryName(outPath)!);
        File.WriteAllLines(outPath, dump);

        // ── Corroboration where SPC/lineage works (independent of SPC) ──
        Assert.True(toyM.BigClusters >= 3 && toyM.BigPurity >= 0.80,
            $"HDBSCAN should recover the 3 toy stripes at high purity. {toyM}");
        Assert.True(irisM.BigClusters >= 2 && irisM.BigPurity >= 0.60,
            $"HDBSCAN should recover the Iris setosa split (the versicolor/virginica overlap caps purity). {irisM}");
        Assert.True(landsatM.ClassesCaptured >= 3 && landsatM.BigPurity >= 0.55,
            $"HDBSCAN (z-normed) should recover ≥3 land-cover types — independent baseline (SPC-lineage does better). {landsatM}");

        // ── The ISOLET wall, confirmed by an independent method ──
        // The falsifiable control: an unsupervised density method ALSO cannot
        // separate the 26 spoken letters (it captures far fewer). If a future
        // change makes HDBSCAN separate many letters unsupervised, this fires —
        // a signal the wall moved, worth re-opening the parity.
        Assert.Equal(26, isoM.NClasses);
        Assert.True(isoM.ClassesCaptured <= 8,
            $"Control: unsupervised HDBSCAN must NOT separate the ISOLET letters (the wall). {isoM}");
    }

    private readonly record struct BenchMetrics(
        string Name, int NClasses, int Clusters, int BigClusters, double BigPurity, int Covered, int Noise, int ClassesCaptured)
    {
        public override string ToString() =>
            $"{Name}: classes={NClasses}, clusters={Clusters}, bigClusters={BigClusters}, purity={BigPurity:F3}, " +
            $"covered={Covered}, noise={Noise}, captured={ClassesCaptured}";
    }

    /// <summary>
    /// Leaf selection should recover MORE Landsat land-cover structure than EOM —
    /// the diagnosis (EOM under-segments multi-class data) made concrete on the
    /// shared selector axis. Validation-independent: purity/capture are scored
    /// against the TRUE land-cover labels, never SPC's output. Landsat-only so it
    /// stays cheap (no ISOLET PCA).
    /// </summary>
    [Fact]
    public void Hdbscan_LeafSelection_RecoversMoreLandsatStructure_ThanEom()
    {
        var landsat = SpcUserSession.FromCsv(Datasets.Path("landsat.csv"), null, false, ',');
        double[][] z = ZNorm(landsat.Features);

        var eom  = Measure("landsat-eom",  z, landsat.Labels, mcs: 50, mp: 5, method: ClusterSelectionMethod.Eom);
        var leaf = Measure("landsat-leaf", z, landsat.Labels, mcs: 50, mp: 5, method: ClusterSelectionMethod.Leaf);

        string outPath = Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "../../../../../artifacts/hdbscan-leaf-vs-eom.tsv"));
        Directory.CreateDirectory(Path.GetDirectoryName(outPath)!);
        File.WriteAllLines(outPath, new[]
        {
            "selector\tclusters\tbigClusters\tbigPurity\tcovered\tnoise\tclassesCaptured",
            $"eom\t{eom.Clusters}\t{eom.BigClusters}\t{eom.BigPurity:F3}\t{eom.Covered}\t{eom.Noise}\t{eom.ClassesCaptured}",
            $"leaf\t{leaf.Clusters}\t{leaf.BigClusters}\t{leaf.BigPurity:F3}\t{leaf.Covered}\t{leaf.Noise}\t{leaf.ClassesCaptured}",
        });

        // Leaf is structurally never coarser than EOM, and on under-segmented
        // Landsat it recovers strictly finer structure.
        Assert.True(leaf.Clusters > eom.Clusters,
            $"leaf should recover more clusters than EOM on Landsat. eom={eom}, leaf={leaf}");
        // Recovering more clusters must not cost land-cover coverage.
        Assert.True(leaf.ClassesCaptured >= eom.ClassesCaptured,
            $"leaf should capture at least as many land-cover types as EOM. eom={eom}, leaf={leaf}");
    }

    private static BenchMetrics Measure(
        string name, double[][] features, int[] trueLabels, int mcs, int mp,
        ClusterSelectionMethod method = ClusterSelectionMethod.Eom)
    {
        var res = HdbscanSession.Run(features, new HdbscanSettings
        {
            MinPts = mp,
            MinClusterSize = mcs,
            AllowSingleCluster = false,
            ClusterSelectionMethod = method,
            Metric = "euclidean"
        });

        var big = res.Labels
            .Select((l, i) => (l, i))
            .Where(x => x.l >= 0)
            .GroupBy(x => x.l)
            .Select(g =>
            {
                int size = g.Count();
                var maj = g.GroupBy(x => trueLabels[x.i]).OrderByDescending(s => s.Count()).First();
                return (Size: size, Major: maj.Key, MajCount: maj.Count());
            })
            .Where(c => c.Size >= mcs)
            .ToList();

        int covered = big.Sum(c => c.Size);
        double purity = covered == 0 ? 0.0 : (double)big.Sum(c => c.MajCount) / covered;
        int captured = big.Select(c => c.Major).Distinct().Count();
        int noise = res.Labels.Count(l => l < 0);
        return new BenchMetrics(name, trueLabels.Distinct().Count(), res.ClusterCount, big.Count, purity, covered, noise, captured);
    }

    /// <summary>Per-feature z-score — HDBSCAN-euclidean is scale-sensitive and has
    /// no built-in normalization, so raw multi-variance features can degenerate.</summary>
    private static double[][] ZNorm(double[][] x)
    {
        int n = x.Length, d = x[0].Length;
        var mean = new double[d];
        var sd = new double[d];
        foreach (var r in x) for (int j = 0; j < d; j++) mean[j] += r[j];
        for (int j = 0; j < d; j++) mean[j] /= n;
        foreach (var r in x) for (int j = 0; j < d; j++) { double e = r[j] - mean[j]; sd[j] += e * e; }
        for (int j = 0; j < d; j++) sd[j] = Math.Sqrt(sd[j] / n);
        return x.Select(r =>
        {
            var o = new double[d];
            for (int j = 0; j < d; j++) o[j] = sd[j] > 1e-12 ? (r[j] - mean[j]) / sd[j] : 0.0;
            return o;
        }).ToArray();
    }

    /// <summary>Stream the gzipped 617-dim ISOLET CSV: 617 features + an integer
    /// letter label (1..26) per row.</summary>
    private static (double[][] Features, int[] Labels) LoadIsoletGz(string gzPath)
    {
        using var fs = File.OpenRead(gzPath);
        using var gz = new GZipStream(fs, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var feats = new List<double[]>(8000);
        var labels = new List<int>(8000);
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (string.IsNullOrWhiteSpace(line)) continue;
            string[] parts = line.Split(',');
            int d = parts.Length - 1;
            var row = new double[d];
            for (int j = 0; j < d; j++)
                row[j] = double.Parse(parts[j], CultureInfo.InvariantCulture);
            feats.Add(row);
            labels.Add((int)Math.Round(double.Parse(parts[d], CultureInfo.InvariantCulture)));
        }
        return (feats.ToArray(), labels.ToArray());
    }
}
