using System;
using System.Collections.Generic;
using System.Text.Json;
using Synthetic.Euclidean;
using Viz;
using Viz.Adapters.Synthetic;

namespace VizApi;

// Initial composition for the editor; the scientific generator's defaults stay separate.
internal sealed record EyeTorusRequest(
    int CentralPoints = 4000, double CentralMajorR = 1.7, double CentralMinorR = 0.35,
    double[]? CentralOffset = null,
    int UpperPoints = 2500, double UpperMajorR = 3.1, double UpperMinorR = 0.22,
    double UpperAngle = 90, double UpperSpan = 150, double[]? UpperOffset = null,
    int LowerPoints = 3500, double LowerMajorR = 3.8, double LowerMinorR = 0.4,
    double LowerAngle = 270, double LowerSpan = 130, double[]? LowerOffset = null,
    int PupilPoints = 2000, double PupilDilation = 0.75, double[]? PupilOffset = null,
    double ZThickness = 1, double HalfArcTaper = 0.65,
    string CrossSection = "Shell", double StructureNoiseSigma = 0.025,
    double DensityGradientStrength = 0.4, double BackgroundDensityRatio = 0,
    double SpectralExponent = 1, double BackgroundContrast = 1, int Seed = 42);

internal static class EyeTorusScene
{
    public static ScenePackage Build(EyeTorusRequest request)
    {
        Validate(request);
        var config = new EyeTorusToy.EyeTorusToyConfig
        {
            CentralPoints = request.CentralPoints,
            CentralMajorR = request.CentralMajorR,
            CentralMinorR = request.CentralMinorR,
            CentralOffset = request.CentralOffset ?? [0, 0, 0],
            UpperPoints = request.UpperPoints,
            UpperMajorR = request.UpperMajorR,
            UpperMinorR = request.UpperMinorR,
            UpperOffset = request.UpperOffset ?? [0, -0.5, 0],
            UpperArcStart = (request.UpperAngle - request.UpperSpan / 2) * Math.PI / 180,
            UpperArcEnd = (request.UpperAngle + request.UpperSpan / 2) * Math.PI / 180,
            LowerPoints = request.LowerPoints,
            LowerMajorR = request.LowerMajorR,
            LowerMinorR = request.LowerMinorR,
            LowerOffset = request.LowerOffset ?? [0, 1, 0],
            LowerArcStart = (request.LowerAngle - request.LowerSpan / 2) * Math.PI / 180,
            LowerArcEnd = (request.LowerAngle + request.LowerSpan / 2) * Math.PI / 180,
            PupilPoints = request.PupilPoints,
            PupilDilation = request.PupilDilation,
            PupilOffset = request.PupilOffset ?? [0, 0, 0],
            ZThickness = request.ZThickness,
            HalfArcTaper = request.HalfArcTaper,
            CrossSection = Enum.Parse<CrossSectionShape>(request.CrossSection),
            StructureNoiseSigma = request.StructureNoiseSigma,
            DensityGradientStrength = request.DensityGradientStrength,
            BackgroundDensityRatio = request.BackgroundDensityRatio,
            SpectralExponent = request.SpectralExponent,
            BackgroundContrast = request.BackgroundContrast,
            Seed = request.Seed,
        };
        var dataset = EyeTorusToy.Generate(config);
        var adapted = new SyntheticDatasetAdapter().Adapt(dataset);
        var parameters = new Dictionary<string, object>(dataset.Parameters);
        var json = JsonSerializer.SerializeToElement(request, new JsonSerializerOptions(JsonSerializerDefaults.Web));
        foreach (var property in json.EnumerateObject()) parameters[property.Name] = property.Value.Clone();
        parameters["centralOffset"] = config.CentralOffset;
        parameters["upperOffset"] = config.UpperOffset;
        parameters["lowerOffset"] = config.LowerOffset;
        parameters["pupilOffset"] = config.PupilOffset;
        int backgroundLabel = (int)dataset.Parameters["backgroundLabel"];
        var names = new Dictionary<int, string> { [0] = "Iris", [1] = "Upper lid", [2] = "Lower lid" };
        var colors = new Dictionary<int, string> { [0] = "#00d2df", [1] = "#ffab00", [2] = "#8195ad" };
        if (request.PupilPoints > 0) { names[3] = "Pupil"; colors[3] = "#a78bfa"; }
        names[backgroundLabel] = "Background";
        colors[backgroundLabel] = "#48505f";
        parameters["componentNames"] = names;
        parameters["componentColors"] = colors;
        var view = new VizDataset(adapted.Points, adapted.LabelLayers, adapted.NodeSignalLayers,
            adapted.EdgeLayers, adapted.GaussianLayers, adapted.TemporalSequences, adapted.SpineLayers,
            generatorParams: parameters, generatorParamSchema: SchemaCatalog.ForGenerator("EyeTorusToy"));
        return SceneBuilder.Build(view, new SceneDescriptor
        {
            Title = "EyeTorusToy — live editor",
            Hints = new SceneRenderHints { ShowSpineOverlays = false, ShowTangentBases = false },
        });
    }

    private static void Validate(EyeTorusRequest r)
    {
        foreach (int count in new[] { r.CentralPoints, r.UpperPoints, r.LowerPoints, r.PupilPoints })
            if (count < 0 || count > 20000) throw new ArgumentException("Each component needs 0–20,000 points.");
        if (r.CentralPoints + r.UpperPoints + r.LowerPoints + r.PupilPoints == 0)
            throw new ArgumentException("Give at least one component a nonzero point count.");
        foreach (var (name, major, minor) in new[] {
            ("Iris", r.CentralMajorR, r.CentralMinorR),
            ("Upper lid", r.UpperMajorR, r.UpperMinorR),
            ("Lower lid", r.LowerMajorR, r.LowerMinorR) })
        {
            Range(name + " ring radius", major, 0.1, 8);
            Range(name + " tube radius", minor, 0.01, 2);
            if (minor >= major) throw new ArgumentException(name + ": tube radius must be smaller than ring radius.");
        }
        foreach (var offset in new[] { r.CentralOffset, r.UpperOffset, r.LowerOffset, r.PupilOffset })
        {
            if (offset is null) continue;
            if (offset.Length != 3) throw new ArgumentException("Component positions need X, Y, and Z.");
            foreach (double value in offset) Range("Position", value, -8, 8);
        }
        Range("Upper angle", r.UpperAngle, -360, 360);
        Range("Lower angle", r.LowerAngle, -360, 360);
        Range("Upper span", r.UpperSpan, 1, 360);
        Range("Lower span", r.LowerSpan, 1, 360);
        Range("Pupil dilation", r.PupilDilation, 0, 0.98);
        Range("Tube depth", r.ZThickness, 0, 2);
        Range("Taper", r.HalfArcTaper, 0, 1);
        Range("Jitter", r.StructureNoiseSigma, 0, 0.5);
        Range("Density gradient", r.DensityGradientStrength, 0, 1);
        Range("Background ratio", r.BackgroundDensityRatio, 0, 0.5);
        Range("Noise spectrum", r.SpectralExponent, 0, 3);
        Range("Background contrast", r.BackgroundContrast, 0, 2);
        if (!Enum.TryParse<CrossSectionShape>(r.CrossSection, out var shape) || !Enum.IsDefined(shape))
            throw new ArgumentException("Cross-section must be Shell, Solid, or Ribbon.");
    }

    private static void Range(string name, double value, double min, double max)
    {
        if (!double.IsFinite(value) || value < min || value > max)
            throw new ArgumentException($"{name} must be between {min} and {max}.");
    }
}
