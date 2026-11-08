using System;
using System.Collections.Generic;
using Maths.Regression.GPrior;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// One chain's state (D04): the structure, its terms in the design's column order with their fit, and u = ln τ, the
/// coordinate τ's slice step moves (ln of a held τ), with, under a hierarchical σ² prior, v = ln r, the coordinate the rate's step
/// moves. Immutable: a τ update shares the structure, and a state's structure identity, which the chain's design tracks, is its
/// terms' array.
/// </summary>
public sealed class BassState
{
    private readonly BassTerm[] _terms;

    internal BassState(BassTerm[] terms, BassStructureFit fit, double logTau, double logRate = double.NaN)
    {
        _terms = terms;
        Fit = fit;
        LogTau = logTau;
        LogRate = logRate;
    }

    /// <summary>The terms, in the design's column order.</summary>
    public IReadOnlyList<BassTerm> Terms => _terms;

    /// <summary>M, the terms.</summary>
    public int Count => _terms.Length;

    /// <summary>The structure's fit in the centred basis, or its refusal.</summary>
    public BassStructureFit Fit { get; }

    /// <summary>u = ln τ.</summary>
    public double LogTau { get; }

    /// <summary>τ = eᵘ.</summary>
    public double Tau => Math.Exp(LogTau);

    /// <summary>v = ln r, σ²'s rate under a hierarchical prior; NaN under every other, which has no such coordinate.</summary>
    public double LogRate { get; }

    /// <summary>v for a chain started from this prior draw: a draw far below the doubles starts at their end, since from there it would
    /// take the rate's step long to climb, and burn-in forgets a start. A prior rung and SMC's particles keep the draw itself.</summary>
    internal double StartLogRate => Math.Max(LogRate, ResolvedSigma2Rate.LogSmallestRate);

    /// <summary>This structure at another τ.</summary>
    public BassState WithLogTau(double logTau) => new(_terms, Fit, logTau, LogRate);

    /// <summary>This structure at another rate.</summary>
    public BassState WithLogRate(double logRate) => new(_terms, Fit, LogTau, logRate);

    internal BassTerm[] TermArray => _terms;

    /// <summary>Whether <paramref name="other"/> holds this state's structure, the same terms' array.</summary>
    internal bool SharesStructure(BassState other) => ReferenceEquals(_terms, other._terms);
}

/// <summary>A retained structure: sorted ids into the fit's <see cref="BassTermTable"/>, an unordered collection
/// (D01).</summary>
public readonly record struct BassStructure(int[] TermIds);

/// <summary>
/// One retained draw, its generative state (D04): the structure, τ, σ², the coefficients and λ, and σ²'s rate under a
/// hierarchical prior. The coefficients are in
/// the uncentred basis, R's: the intercept, then one per term in <see cref="BassStructure.TermIds"/>' order, so the
/// draw's function is f = a₀ + Σⱼ aⱼBⱼ, with Bⱼ the term's column (its case side times its output side for a direct
/// functional response). <see cref="Working"/> is the working response the draw was composed on, its offset and weights being part of
/// that state, shared with the chain that held it; null when it is the response itself, or when <see cref="Latents"/> stand for it.
/// <see cref="Latents"/> are the observation model's latents, kept alone by a model that rebuilds its working response from them
/// (<see cref="IBassObservationModel.Rebuild"/>): the generative state, at N doubles where the working response holds up to 4N.
/// <see cref="GroupEffects"/> are a random-intercept fit's u over its observed levels (E21; <see cref="RandomInterceptBlock"/>), drawn
/// with the rest, so a case of level g reads f + u_g; null for a fit without the block.
/// </summary>
public sealed record BassDraw(BassStructure Structure, double Tau, double Sigma2, double[] Coefficients, double Lambda, double? Sigma2Rate = null,
                              WorkingResponse? Working = null, ObservationLatents? Latents = null, double[]? GroupEffects = null)
{
    /// <summary>The working response this draw reads on <paramref name="model"/>: its own when it was composed on that response, or rebuilt
    /// from its latents; the model's otherwise, since a draw read on other data has no offset or weights there.</summary>
    internal WorkingResponse WorkingOn(BassModel model)
    {
        if (Working is not { } own)
            return Latents is { } latents && model.Observation.Rebuild(model.Working, latents) is { } rebuilt ? rebuilt : model.Working;
        BassResponse response = model.Data.Response;
        return ReferenceEquals(own.Response, response) || own.Response.Values.SequenceEqual(response.Values) ? own : model.Working;
    }

    /// <summary>The latents this draw reads on <paramref name="model"/>, its own or its working response's, without a rebuild.</summary>
    internal ObservationLatents? LatentsOn(BassModel model) => Latents ?? WorkingOn(model).Latents;
}
