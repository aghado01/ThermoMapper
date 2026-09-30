using System;
using System.Collections.Generic;
using System.IO;
using Repo.TestHarness;
using Synthetic;
using Synthetic.Euclidean;
using Viz;
using Viz.Adapters.Synthetic;
using Viz.Renderers;

internal static class FixtureGallery
{
    // Point-cloud portraits: no graph construction or clustering is run.
    public static void Render(string? outputDirectory)
    {
        var run = outputDirectory is null
            ? HarnessArtifacts.Create("smoke-runs", "VizCoreSmoke", "Gallery")
            : HarnessArtifacts.Attach("smoke-runs", "VizCoreSmoke", "Gallery", outputDirectory);
        var scenes = new List<object>();
        var eye = EyeTorusToy.Generate();

        Write("eye-torus", "EyeTorusToy — seed 42", eye);
        Write("eye-torus-xy", "EyeTorusToy — XY projection of the same samples",
            EyeTorusToy.FlattenToPlane(eye));
        Write("eye-torus-depth", "EyeTorusToy — full tube depth, no background",
            EyeTorusToy.Generate(new EyeTorusToy.EyeTorusToyConfig
            {
                ZThickness = 1.0,
                BackgroundDensityRatio = 0.0,
                Seed = 42,
            }));
        Write("crescent-ellipsoid", "Crescent + ellipsoid — seed 42",
            CrescentAndEllipsoid.Generate(seed: 42));
        Write("mobius-ellipsoid", "Möbius + ellipsoid — seed 42",
            MobiusAndEllipsoid.Generate(placement: MobiusPlacement.CenterCrossOrtho,
                dimensions: 3, seed: 42));

        string summary = run.WriteJson("summary.json", new { Scenes = scenes });
        Console.WriteLine($"RunRoot\t{run.RunDirectory}");
        Console.WriteLine($"Manifest\t{run.ManifestPath}");
        Console.WriteLine($"Summary\t{summary}");

        void Write(string stem, string title, SyntheticDataset dataset)
        {
            var scene = SceneBuilder.Build(new SyntheticDatasetAdapter().Adapt(dataset),
                new SceneDescriptor
                {
                    Title = title,
                    Hints = new SceneRenderHints { ShowSpineOverlays = false, ShowTangentBases = false },
                });
            string path = Path.Combine(run.RunDirectory, stem + ".html");
            using (var stream = File.Create(path))
                new ThreeJsHtmlRenderTarget().Render(scene, stream);
            using (var stream = File.Create(Path.Combine(run.RunDirectory, stem + ".json")))
                new JsonExportRenderTarget().Render(scene, stream);
            scenes.Add(new
            {
                Title = title,
                Html = Path.GetFileName(path),
                PointCount = scene.Points.N,
                Dimension = scene.Points.D,
                dataset.Parameters,
            });
        }
    }
}
