namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A BASS run's account: the draws per chain it used and its moves' acceptance over the sampling phase.
/// </summary>
public sealed record BassDiagnostics(int SamplesUsed, double AcceptanceRate);
