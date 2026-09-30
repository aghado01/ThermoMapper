namespace Viz;

internal static class EyeTorusSchema
{
    public static readonly GeneratorParamSchema Value = new()
    {
        Sections =
        [
            new() { Label = "Iris · full torus", Params = [
                Count("centralPoints", "Iris points"),
                Float("centralMajorR", "Iris ring radius", 0.1, 8),
                Float("centralMinorR", "Iris tube radius", 0.01, 2),
                Position("centralOffset", "Iris") ] },
            new() { Label = "Upper lid · gold arc", Params = [
                Count("upperPoints", "Upper points"),
                Float("upperMajorR", "Upper ring radius", 0.1, 8),
                Float("upperMinorR", "Upper tube radius", 0.01, 2),
                Float("upperAngle", "Upper angle °", -360, 360, 1),
                Float("upperSpan", "Upper arc span °", 1, 360, 1),
                Position("upperOffset", "Upper") ] },
            new() { Label = "Lower lid · large arc", Params = [
                Count("lowerPoints", "Lower points"),
                Float("lowerMajorR", "Lower ring radius", 0.1, 8),
                Float("lowerMinorR", "Lower tube radius", 0.01, 2),
                Float("lowerAngle", "Lower angle °", -360, 360, 1),
                Float("lowerSpan", "Lower arc span °", 1, 360, 1),
                Position("lowerOffset", "Lower") ] },
            new() { Label = "Pupil · solid sphere", Params = [
                Count("pupilPoints", "Pupil points"),
                Float("pupilDilation", "Pupil / iris bore", 0, 0.98),
                Position("pupilOffset", "Pupil") ] },
            new() { Label = "Sampling", Params = [
                Float("zThickness", "Tube depth", 0, 2),
                Float("halfArcTaper", "Lid taper", 0, 1),
                new() { Key = "crossSection", Label = "Cross-section", Type = "enum", EnumValues = ["Shell", "Solid", "Ribbon"] },
                Float("structureNoiseSigma", "Point jitter", 0, 0.5, 0.005),
                Float("densityGradientStrength", "Inner-edge bias", 0, 1),
                new() { Key = "seed", Label = "Seed", Type = "int", Min = 0, Max = 99999 } ] },
            new() { Label = "Background", Params = [
                Float("backgroundDensityRatio", "Background ratio", 0, 0.5),
                Float("spectralExponent", "Noise spectrum", 0, 3),
                Float("backgroundContrast", "Density contrast", 0, 2) ] },
        ],
    };

    private static ParamSpec Count(string key, string label) =>
        new() { Key = key, Label = label, Type = "int", Min = 0, Max = 20000 };
    private static ParamSpec Float(string key, string label, double min, double max, double step = 0.01) =>
        new() { Key = key, Label = label, Type = "float", Min = min, Max = max, Step = step };
    private static ParamSpec Position(string key, string label) =>
        new() { Key = key, Label = label, Type = "vec3", Min = -8, Max = 8, Step = 0.05, VecLabels = ["X", "Y", "Z"] };
}
