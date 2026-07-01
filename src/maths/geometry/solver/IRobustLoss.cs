// ============================================================================
// Losses/IRobustLoss.cs
// ============================================================================
namespace Maths.Geometry.Solver;

public interface IRobustLoss
{
    /// <summary>L² returns true (single Karcher iteration suffices). Iterative losses return false.</summary>
    static abstract bool IsClosedForm { get; }

    /// <summary>L¹ returns true (loss is non-smooth at coincident data points). Smooth losses (L², Huber, Tukey on the bulk) return false.</summary>
    static abstract bool IsSingularAtZero { get; }

    /// <summary>IRLS weight w(r) = ψ(r)/r. Used for both the location iteration and the scatter accumulation.</summary>
    static abstract double Weight(double r);
}

public readonly struct L2Loss : IRobustLoss
{
    public static bool IsClosedForm => true;
    public static bool IsSingularAtZero => false;
    public static double Weight(double r) => 1.0;
}

public readonly struct L1Loss : IRobustLoss
{
    public static bool IsClosedForm => false;
    public static bool IsSingularAtZero => true;
    public static double Weight(double r) => 1.0 / r;   // caller responsible for r > 0
}
