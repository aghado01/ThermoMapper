using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// One chain's state (D04): the structure, its terms in the design's column order with their fit, and u = ln τ, the
/// coordinate τ's slice step moves (ln of a held τ). Immutable: a τ update shares the structure, and a
/// state's structure identity, which the chain's design tracks, is its terms' array.
/// </summary>
public sealed class BassState
{
    private readonly BassTerm[] _terms;

    internal BassState(BassTerm[] terms, BassStructureFit fit, double logTau)
    {
        _terms = terms;
        Fit = fit;
        LogTau = logTau;
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

    /// <summary>This structure at another τ.</summary>
    public BassState WithLogTau(double logTau) => new(_terms, Fit, logTau);

    internal BassTerm[] TermArray => _terms;

    /// <summary>Whether <paramref name="other"/> holds this state's structure, the same terms' array.</summary>
    internal bool SharesStructure(BassState other) => ReferenceEquals(_terms, other._terms);
}

/// <summary>A retained structure: sorted ids into the fit's <see cref="BassTermTable"/>, an unordered collection
/// (D01).</summary>
public readonly record struct BassStructure(int[] TermIds);

/// <summary>
/// One retained draw, its generative state (D04): the structure, τ, σ², the coefficients and λ. The coefficients are in
/// the uncentred basis, R's: the intercept, then one per term in <see cref="BassStructure.TermIds"/>' order, so the
/// draw's function is f = a₀ + Σⱼ aⱼBⱼ, with Bⱼ the term's column (its case side times its output side for a direct
/// functional response).
/// </summary>
public sealed record BassDraw(BassStructure Structure, double Tau, double Sigma2, double[] Coefficients, double Lambda);
