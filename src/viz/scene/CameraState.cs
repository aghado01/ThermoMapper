namespace Viz.Scene;

public readonly record struct Vector3(double X, double Y, double Z);

public enum CameraProjection
{
    Perspective,
    Orthographic,
}

/// <summary>Durable camera state; geometry-model transforms remain coordinate artifacts.</summary>
public sealed record CameraState(
    Vector3 Position,
    Vector3 Target,
    Vector3 Up,
    CameraProjection Projection,
    double FieldOfViewOrScale,
    double Near,
    double Far);
