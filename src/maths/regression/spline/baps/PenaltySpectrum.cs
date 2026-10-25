using System;
using System.Collections.Generic;

namespace Maths.Regression.Spline.Baps;

/// <summary>
/// The joint spectrum of an anisotropic penalty's components, <c>P(λ) = Σⱼ λⱼPⱼ</c>, which share one eigenbasis:
/// component j's eigenvalue <c>γⱼₗ</c> on mode l. P(λ)'s non-zero eigenvalues are <c>Sₗ(λ) = Σⱼ λⱼγⱼₗ</c> over the modes
/// where some <c>γⱼₗ &gt; 0</c>, so its log pseudo-determinant and that log's derivatives in ln λ cost O(ν) at any λ
/// (Bach &amp; Klein 2025, Corollary 3.2 and Proposition 3.4). The modes where every γ is 0 span the penalty's null space.
/// </summary>
public sealed class PenaltySpectrum
{
    private readonly double[][] _gamma;   // [component][non-null mode], in mode order

    /// <param name="gamma">Each component's eigenvalue on each mode, <c>gamma[j][l]</c>, every row over the same modes.</param>
    public PenaltySpectrum(double[][] gamma)
    {
        ArgumentNullException.ThrowIfNull(gamma);
        if (gamma.Length == 0) throw new ArgumentException("A spectrum needs at least one component.", nameof(gamma));
        int modes = gamma[0]?.Length ?? throw new ArgumentException("A component's eigenvalues are null.", nameof(gamma));
        foreach (double[] row in gamma)
        {
            if (row is null || row.Length != modes)
                throw new ArgumentException("Every component needs an eigenvalue on every mode.", nameof(gamma));
            foreach (double g in row)
                if (!(g >= 0.0) || double.IsInfinity(g))
                    throw new ArgumentException("Penalty eigenvalues must be finite and non-negative.", nameof(gamma));
        }

        var active = new List<int>(modes);
        for (int l = 0; l < modes; l++)
        {
            bool penalized = false;
            foreach (double[] row in gamma) penalized |= row[l] > 0.0;
            if (penalized) active.Add(l);
        }
        _gamma = new double[gamma.Length][];
        for (int j = 0; j < gamma.Length; j++)
        {
            _gamma[j] = new double[active.Count];
            for (int a = 0; a < active.Count; a++) _gamma[j][a] = gamma[j][active[a]];
        }
        Modes = modes;
        Rank = active.Count;
    }

    /// <summary>The number of components.</summary>
    public int Components => _gamma.Length;

    /// <summary>The number of modes, ν.</summary>
    public int Modes { get; }

    /// <summary>The number of non-null modes: P(λ)'s rank at any λ with every λⱼ &gt; 0.</summary>
    public int Rank { get; }

    /// <summary>The null space's dimension, ν − rank.</summary>
    public int Nullity => Modes - Rank;

    /// <summary><c>ln Det P(λ) = Σₗ ln Σⱼ λⱼγⱼₗ</c> over the non-null modes.</summary>
    public double LogPseudoDeterminant(ReadOnlySpan<double> lambdas)
    {
        CheckLambdas(lambdas);
        double sum = 0.0;
        for (int l = 0; l < Rank; l++)
        {
            double s = 0.0;
            for (int j = 0; j < _gamma.Length; j++) s += lambdas[j] * _gamma[j][l];
            sum += Math.Log(s);
        }
        return sum;
    }

    /// <summary>
    /// <c>ln Det P(λ)</c> with its gradient and Hessian in <c>uⱼ = ln λⱼ</c>. With <c>pⱼₗ = λⱼγⱼₗ/Sₗ</c>, the gradient is
    /// <c>Σₗ pⱼₗ</c> and the Hessian <c>Σₗ (δⱼₖpⱼₗ − pⱼₗpₖₗ)</c>, written row-major into <paramref name="hessian"/>. Each mode's
    /// term is a multinomial covariance, so the Hessian is positive semi-definite: the log pseudo-determinant is convex in
    /// ln λ.
    /// </summary>
    public double LogPseudoDeterminant(ReadOnlySpan<double> lambdas, Span<double> gradient, Span<double> hessian)
    {
        CheckLambdas(lambdas);
        int k = _gamma.Length;
        if (gradient.Length != k || hessian.Length != k * k)
            throw new ArgumentException($"The gradient needs {k} entries and the Hessian {k * k}.", nameof(gradient));
        gradient.Clear();
        hessian.Clear();
        Span<double> p = k <= 16 ? stackalloc double[k] : new double[k];
        double sum = 0.0;
        for (int l = 0; l < Rank; l++)
        {
            double s = 0.0;
            for (int j = 0; j < k; j++) { p[j] = lambdas[j] * _gamma[j][l]; s += p[j]; }
            sum += Math.Log(s);
            for (int j = 0; j < k; j++) p[j] /= s;
            for (int j = 0; j < k; j++)
            {
                gradient[j] += p[j];
                hessian[j * k + j] += p[j];
                for (int i = 0; i < k; i++) hessian[j * k + i] -= p[j] * p[i];
            }
        }
        return sum;
    }

    private void CheckLambdas(ReadOnlySpan<double> lambdas)
    {
        if (lambdas.Length != _gamma.Length)
            throw new ArgumentException($"Need one smoothing per component ({_gamma.Length}).", nameof(lambdas));
        foreach (double lambda in lambdas)
            if (!(lambda > 0.0) || double.IsInfinity(lambda))
                throw new ArgumentOutOfRangeException(nameof(lambdas), "Smoothing parameters must be positive and finite.");
    }
}
