using System;

namespace Maths.Geometry.DimReduction
{
    /// <summary>
    /// Declarative configuration for <see cref="SubspaceAnnealer.Compute"/>: the proposal mixture,
    /// the acceptance-adaptive step controller, and the cooling schedule. A pure data record — the
    /// annealer consumes it and owns all behavior. Defaults encode the two-plane-primary design
    /// with a small isotropic share as ergodicity insurance.
    /// </summary>
    public sealed record SubspaceAnnealerOptions
    {
        /// <summary>
        /// Fraction of proposals drawn as isotropic full-rank horizontal tangents instead of
        /// two-plane Givens rotations. 0 anneals with two-plane moves only; 1 recovers the
        /// isotropic-only proposal kind. Each retained column's Givens angle and the isotropic
        /// kind adapt their own step scale within the shared bounds.
        /// </summary>
        public double IsotropicFraction { get; init; } = 0.1;

        /// <summary>
        /// Acceptance rate the step controller steers toward — the classic random-walk
        /// Metropolis operating point. Higher acceptance grows the step, lower shrinks it.
        /// </summary>
        public double TargetAcceptance { get; init; } = 0.25;

        /// <summary>Geodesic step scale at iteration 0, clamped into [<see cref="StepFloor"/>,
        /// <see cref="StepCeiling"/>] before the anneal starts.</summary>
        public double InitialStep { get; init; } = 0.1;

        /// <summary>Lower bound on the adaptive step scale, keeping proposals above pure
        /// diffusion when acceptance runs low.</summary>
        public double StepFloor { get; init; } = 1e-3;

        /// <summary>Upper bound on the adaptive step scale. The default π/2 is a quarter turn —
        /// the principal-angle injectivity radius, past which a two-plane move walks back toward
        /// the start of its geodesic circle.</summary>
        public double StepCeiling { get; init; } = Math.PI / 2.0;

        /// <summary>Metropolis temperature at iteration 0.</summary>
        public double InitialTemperature { get; init; } = 1.0;

        /// <summary>Geometric cooling factor per iteration. Cooling governs the Metropolis
        /// temperature only; the step scale is governed by acceptance.</summary>
        public double CoolingRate { get; init; } = 0.99;
    }
}
