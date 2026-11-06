using System;
using Maths.Rng;

namespace Maths.Distributions;

/// <summary>
/// The modified half-normal MHN(α, β, γ), density ∝ x^(α−1)·exp(−βx² + γx) on x &gt; 0, α &gt; 0, β &gt; 0 (Sun, Kong and Pal 2021): the
/// law of w^(−½) for a normal variance-mean mixture's global scale w when the mean carries √w·β (E1). Three exact draws:
/// <list type="bullet">
/// <item><b>α ≥ 1:</b> log-concave in x, so the tangent hat of <see cref="LogConcaveEnvelope"/> about the mode, the log ratio
/// (α − 1)(ln(1 + t) − t) − βδ², t = δ/x*, with the mode's identity cancelling the linear terms exactly.</item>
/// <item><b>α &lt; 1, γ ≤ 0:</b> log-concave in y = ln x, where the log density is αy − βe^(2y) + γe^y; the hat there, its ratio
/// βz*²(2δ − (e^(2δ) − 1)) + γz*((e^δ − 1) − δ), z* = e^(y*).</item>
/// <item><b>α &lt; 1, γ &gt; 0:</b> not log-concave in either. Expanding e^(γx) gives a mixture: k with weight ∝ (γ/√β)ᵏΓ((α + k)/2)/k!,
/// then x² ~ Gamma((α + k)/2, rate β). The weights are summed in logs until a term falls e⁻⁴⁰ below the largest past the
/// peak, the rest below the rounding of 1.</item>
/// </list>
/// The paper's own samplers are rejection schemes per branch; these are exact by other routes, checked against quadrature on
/// every branch. No normalizing constant lives here. A readonly value type.
/// </summary>
public readonly struct ModifiedHalfNormal
{
    private const int MaxSeriesTerms = 1_000_000;

    private ModifiedHalfNormal(double alpha, double beta, double gamma)
    {
        Alpha = alpha;
        Beta = beta;
        Gamma = gamma;
    }

    /// <summary>MHN(<paramref name="alpha"/>, <paramref name="beta"/>, <paramref name="gamma"/>): α and β positive and finite, γ finite.</summary>
    public static ModifiedHalfNormal Of(double alpha, double beta, double gamma)
    {
        if (!(alpha > 0.0 && double.IsFinite(alpha))) throw new ArgumentOutOfRangeException(nameof(alpha), "α must be positive and finite.");
        if (!(beta > 0.0 && double.IsFinite(beta))) throw new ArgumentOutOfRangeException(nameof(beta), "β must be positive and finite.");
        if (!double.IsFinite(gamma)) throw new ArgumentOutOfRangeException(nameof(gamma), "γ must be finite.");
        return new ModifiedHalfNormal(alpha, beta, gamma);
    }

    /// <summary>α, the power: x^(α−1).</summary>
    public double Alpha { get; }

    /// <summary>β, the quadratic's coefficient.</summary>
    public double Beta { get; }

    /// <summary>γ, the linear term's coefficient.</summary>
    public double Gamma { get; }

    /// <summary>ln of the unnormalized density, (α − 1)·ln x − βx² + γx; −∞ for x ≤ 0.</summary>
    public double LogKernel(double x) => x > 0.0 ? (Alpha - 1.0) * Math.Log(x) - Beta * x * x + Gamma * x : double.NegativeInfinity;

    /// <summary>The density's mode: the positive root of 2βx² − γx − (α − 1) = 0 for α ≥ 1 (0 when α = 1 and γ ≤ 0); 0 below α = 1,
    /// where the density is unbounded there.</summary>
    public double Mode => Alpha < 1.0 ? 0.0 : PositiveRoot(Alpha - 1.0, Beta, Gamma);

    // The positive root of 2βx² − γx − c = 0, c ≥ 0, without cancellation: (γ + √(γ² + 8βc))/(4β), or 2c/(√(γ² + 8βc) − γ) for γ < 0.
    private static double PositiveRoot(double c, double beta, double gamma)
    {
        double root = Math.Sqrt(gamma * gamma + 8.0 * beta * c);
        return gamma >= 0.0 ? (gamma + root) / (4.0 * beta) : 2.0 * c / (root - gamma);
    }

    /// <summary>One draw, by the branch the parameters fall in.</summary>
    public double Sample(Xoshiro256PlusPlus rng)
    {
        ArgumentNullException.ThrowIfNull(rng);
        if (Alpha >= 1.0)
        {
            double mode = PositiveRoot(Alpha - 1.0, Beta, Gamma);
            if (mode == 0.0)
            {
                // α = 1, γ ≤ 0: the density falls from 0, h(δ) = −βδ² + γδ.
                var edge = new EdgeRatio(Beta, Gamma);
                double reach = LogConcaveEnvelope.Reach(-Gamma, 2.0 * Beta);
                return Math.Max(LogConcaveEnvelope.Draw(rng, edge, 0.0, double.PositiveInfinity, reach, reach), 0.0);
            }
            var h = new LinearRatio(Alpha - 1.0, Beta, mode);
            double curvature = (Alpha - 1.0) / (mode * mode) + 2.0 * Beta, step = LogConcaveEnvelope.Reach(0.0, curvature);
            // The left reach stays inside the span: at δ = −x* the log density is −∞, and a design point there has no tangent.
            double delta = LogConcaveEnvelope.Draw(rng, h, mode, double.PositiveInfinity, Math.Min(step, 0.5 * mode), step);
            return Math.Max(mode + delta, double.Epsilon);
        }
        if (Gamma <= 0.0)
        {
            double z = PositiveRoot(Alpha, Beta, Gamma);   // e^(y*): 2βz² − γz − α = 0
            var h = new LogRatio(Beta * z * z, Gamma * z);
            double reach = LogConcaveEnvelope.Reach(0.0, 4.0 * Beta * z * z - Gamma * z);
            double delta = LogConcaveEnvelope.Draw(rng, h, double.PositiveInfinity, double.PositiveInfinity, reach, reach);
            return z * Math.Exp(delta);
        }
        return SeriesDraw(rng);
    }

    // α < 1, γ > 0: k from its weights, then x = √g with g ~ Gamma((α + k)/2, rate β).
    private double SeriesDraw(Xoshiro256PlusPlus rng)
    {
        double logRatio = Math.Log(Gamma) - 0.5 * Math.Log(Beta);
        var logs = new System.Collections.Generic.List<double>();
        double best = double.NegativeInfinity;
        for (int k = 0; k < MaxSeriesTerms; k++)
        {
            double term = k * logRatio + SpecialFunctions.LogGamma(0.5 * (Alpha + k)) - SpecialFunctions.LogGamma(k + 1.0);
            logs.Add(term);
            if (term > best) best = term;
            else if (term < best - 40.0) break;
            if (k == MaxSeriesTerms - 1)
                throw new NotSupportedException($"MHN({Alpha}, {Beta}, {Gamma}) below α = 1 needs more than {MaxSeriesTerms} terms of its series.");
        }
        double total = 0.0;
        foreach (double t in logs) total += Math.Exp(t - best);
        double pick = rng.NextDouble() * total, run = 0.0;
        int chosen = logs.Count - 1;
        for (int k = 0; k < logs.Count; k++)
        {
            run += Math.Exp(logs[k] - best);
            if (pick < run)
            {
                chosen = k;
                break;
            }
        }
        return Math.Sqrt(Distributions.Gamma.ShapeRate(0.5 * (Alpha + chosen), Beta).Sample(rng));
    }

    // About the mode x* of α ≥ 1: (α − 1)(ln(1 + t) − t) − βδ², t = δ/x*, the linear terms cancelled by 2βx* − γ = (α − 1)/x*.
    internal readonly struct LinearRatio(double power, double beta, double mode) : IConcaveLogRatio
    {
        public double Value(double delta) => power * LogOnePlusLess(delta / mode) - beta * delta * delta;

        public double Slope(double delta)
        {
            double t = delta / mode;
            return -power * t / ((1.0 + t) * mode) - 2.0 * beta * delta;
        }
    }

    // α = 1 with γ ≤ 0, about 0: −βδ² + γδ.
    internal readonly struct EdgeRatio(double beta, double gamma) : IConcaveLogRatio
    {
        public double Value(double delta) => delta * (gamma - beta * delta);

        public double Slope(double delta) => gamma - 2.0 * beta * delta;
    }

    // About y* = ln z* for α < 1, γ ≤ 0: βz*²(2δ − (e^(2δ) − 1)) + γz*((e^δ − 1) − δ), with q = βz*², l = γz*.
    internal readonly struct LogRatio(double q, double l) : IConcaveLogRatio
    {
        public double Value(double delta) => q * -ExpLessLinear(2.0 * delta) + l * ExpLessLinear(delta);

        public double Slope(double delta) => q * 2.0 * -SpecialFunctions.Expm1(2.0 * delta) + l * SpecialFunctions.Expm1(delta);
    }

    // ln(1 + t) − t, by its series where the two would cancel.
    private static double LogOnePlusLess(double t)
    {
        if (Math.Abs(t) > 0.1) return SpecialFunctions.Log1p(t) - t;
        double t2 = t * t;
        return t2 * (-0.5 + t * (1.0 / 3.0 + t * (-0.25 + t * (0.2 + t * (-1.0 / 6.0 + t * (1.0 / 7.0 + t * (-0.125 + t * (1.0 / 9.0 - t / 10.0))))))));
    }

    // eˣ − 1 − x, by its series where the two would cancel.
    private static double ExpLessLinear(double x)
    {
        if (Math.Abs(x) > 0.1) return SpecialFunctions.Expm1(x) - x;
        double x2 = x * x;
        return x2 * (0.5 + x * (1.0 / 6.0 + x * (1.0 / 24.0 + x * (1.0 / 120.0 + x * (1.0 / 720.0 + x * (1.0 / 5040.0 + x / 40320.0))))));
    }
}
