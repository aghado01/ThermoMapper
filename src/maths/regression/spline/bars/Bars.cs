using System;

using Maths.Regression.Spline;

namespace Maths.Regression.Spline.Bars;

/// <summary>
/// Declarative configuration for a BARS run — the "state of the universe" the strict core consumes
/// (<see cref="Bars.Run"/>). Components are supplied as instances (model / prior / kernel / robustness) with
/// sensible defaults; the ensemble and stopping knobs follow. A fluent/chained surface, if wanted, belongs at
/// the REPL boundary, not here.
/// </summary>
public sealed record BarsConfig
{
    /// <summary>Spline degree (3 = cubic) — used when <see cref="Basis"/> is null.</summary>
    public int Degree { get; init; } = 3;

    /// <summary>Carrier realization; null = cubic <see cref="SplineBasis"/>(<see cref="Degree"/>). Use a <see cref="StepBasis"/> for step-function fits.</summary>
    public IBasis? Basis { get; init; }

    /// <summary>Observation model — the likelihood (default exact weighted-Normal).</summary>
    public IObservationModel Model { get; init; } = new WeightedNormalModel();

    /// <summary>Prior on the knot count.</summary>
    public IComplexityPrior Prior { get; init; } = new PoissonPrior(5.0);

    /// <summary>Proposal locality kernel.</summary>
    public IKnotKernel Kernel { get; init; } = new LocalBetaKernel(50.0);

    /// <summary>Optional scale-mixture robustness (null = non-robust least squares).</summary>
    public IWeightResampler? Robustness { get; init; }

    public int Chains { get; init; } = 4;
    public int MasterSeed { get; init; }
    public int BurnIn { get; init; } = 2000;
    public int MaxSamples { get; init; } = 2000;
    public int StartDispersion { get; init; } = 4;
    public int WeightEvery { get; init; } = 25;

    /// <summary>Stop once R̂ ≤ this on all reported functionals (0 = run the full <see cref="MaxSamples"/>).</summary>
    public double RHatTarget { get; init; }

    /// <summary>Samples per round between R̂ checks (≤0 = a single round).</summary>
    public int BatchSize { get; init; }

    /// <summary>Relative prominence (fraction of the fit's range) for the significant-peak count.</summary>
    public double PeakProminence { get; init; } = 0.1;

    /// <summary>
    /// Drop fraction defining each peak's span: the span runs to where the curve falls this fraction of the peak's
    /// prominence below the apex (½ = the FWHM analogue). A <c>DOMAIN-PREMISE</c> — the consumer's notion of "the
    /// span around the peak", not BARS's.
    /// </summary>
    public double SpanDropFraction { get; init; } = 0.5;
}

/// <summary>Strict-core entry point: build the ensemble from a <see cref="BarsConfig"/> and run it.</summary>
public static class Bars
{
    public static BarsResult Run(BarsConfig config, double[] x, double[] y, double[] grid)
    {
        ArgumentNullException.ThrowIfNull(config);
        IBasis basis = config.Basis ?? new SplineBasis(config.Degree);
        var ensemble = new BarsEnsemble(basis, config.Model, config.Prior, config.Kernel, config.Robustness);
        return ensemble.Run(x, y, grid,
            config.Chains, config.MasterSeed, config.BurnIn, config.MaxSamples,
            config.StartDispersion, config.WeightEvery, config.RHatTarget, config.BatchSize, config.PeakProminence,
            config.SpanDropFraction);
    }
}
