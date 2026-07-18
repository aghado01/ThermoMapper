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

        /// <summary>
        /// Reject records the anneal cannot run on. Lives on the DTO so every consumer shares one
        /// gate: <see cref="SubspaceAnnealer.Compute"/> calls it on entry, and drivers that do
        /// expensive setup before annealing (SPRED builds its reference barcode first) call it up
        /// front to fail before that work. Checks are written so NaN fails them — the relational
        /// forms NaN slips through were the hole.
        /// </summary>
        public void Validate()
        {
            if (IsotropicFraction is not (>= 0.0 and <= 1.0))
                throw new ArgumentOutOfRangeException(nameof(IsotropicFraction),
                    "IsotropicFraction is a mixture probability: it must lie in [0, 1].");
            if (TargetAcceptance is not (> 0.0 and < 1.0))
                throw new ArgumentOutOfRangeException(nameof(TargetAcceptance),
                    "TargetAcceptance must lie strictly inside (0, 1) — the step controller's " +
                    "zero-drift update law is degenerate at either endpoint.");
            if (!double.IsFinite(StepFloor) || StepFloor <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(StepFloor),
                    "StepFloor is a geodesic length: it must be finite and positive.");
            if (!double.IsFinite(StepCeiling) || StepCeiling < StepFloor)
                throw new ArgumentOutOfRangeException(nameof(StepCeiling),
                    "StepCeiling must be finite and at least StepFloor — step scales clamp into " +
                    "[floor, ceiling], and an unbounded ceiling lets the controller wrap geodesic circles.");
            if (!double.IsFinite(InitialStep) || InitialStep <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(InitialStep),
                    "InitialStep is a geodesic length: it must be finite and positive.");
            if (!double.IsFinite(InitialTemperature) || InitialTemperature <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(InitialTemperature),
                    "InitialTemperature must be finite and positive — an infinite temperature accepts " +
                    "every proposal and a NaN one silently freezes Metropolis into greedy descent.");
            if (CoolingRate is not (> 0.0 and <= 1.0))
                throw new ArgumentOutOfRangeException(nameof(CoolingRate),
                    "CoolingRate must lie in (0, 1].");
        }
    }
}
