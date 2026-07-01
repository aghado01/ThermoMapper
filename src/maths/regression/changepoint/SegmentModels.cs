using System;

namespace Maths.Regression.Changepoint;

/// <summary>
/// The per-segment model behind the exact change-point DP (Hutter 2005, "any noise / level prior"): for a segment
/// covering data points a..b it supplies the single-segment evidence (the level integrated out) and the posterior
/// level mean/variance (for the regression curve). The DP itself is model-agnostic — it only reads these tables —
/// so swapping a Gaussian likelihood for a heavy-tailed (outlier-robust) one is a model swap, not a DP change.
/// </summary>
public interface ISegmentModel
{
    /// <summary>Number of data points n.</summary>
    int Count { get; }

    /// <summary>
    /// Fill the upper triangle (1 ≤ a ≤ b ≤ n): <paramref name="logEvidence"/>[a,b] = log D(a,b) (level integrated
    /// out), <paramref name="levelMean"/>[a,b] = E[μ | segment], <paramref name="levelVar"/>[a,b] = Var[μ | segment].
    /// Tables are (n+1)×(n+1).
    /// </summary>
    void Fill(double[,] logEvidence, double[,] levelMean, double[,] levelVar);
}

/// <summary>
/// Gaussian noise + Gaussian segment-level prior — the conjugate, closed-form model. The single-segment evidence
/// is the marginal <c>N(y_{a..b}; ν1, σ²I + ρ²11ᵀ)</c>, and the posterior level is the standard conjugate
/// Gaussian; both are O(1) per segment from prefix sums.
/// </summary>
public sealed class GaussianSegmentModel : ISegmentModel
{
    private readonly double[] _y;
    private readonly double _sigma2;
    private readonly double _nu;
    private readonly double _rho2;

    public GaussianSegmentModel(double[] y, double sigma2, double nu, double rho2)
    {
        ArgumentNullException.ThrowIfNull(y);
        if (!(sigma2 > 0.0)) throw new ArgumentOutOfRangeException(nameof(sigma2));
        if (!(rho2 > 0.0)) throw new ArgumentOutOfRangeException(nameof(rho2));
        _y = y; _sigma2 = sigma2; _nu = nu; _rho2 = rho2;
    }

    public int Count => _y.Length;

    public void Fill(double[,] logEvidence, double[,] levelMean, double[,] levelVar)
    {
        int n = _y.Length;
        var rPre = new double[n + 1];
        var qPre = new double[n + 1];
        var yPre = new double[n + 1];
        for (int i = 1; i <= n; i++)
        {
            double ri = _y[i - 1] - _nu;
            rPre[i] = rPre[i - 1] + ri;
            qPre[i] = qPre[i - 1] + ri * ri;
            yPre[i] = yPre[i - 1] + _y[i - 1];
        }

        double logSigma2 = Math.Log(_sigma2);
        double half2Pi = 0.5 * Math.Log(2.0 * Math.PI);
        double invRho2 = 1.0 / _rho2;

        for (int a = 1; a <= n; a++)
            for (int b = a; b <= n; b++)
            {
                int len = b - a + 1;
                double sr = rPre[b] - rPre[a - 1];
                double sq = qPre[b] - qPre[a - 1];
                double t = _sigma2 + len * _rho2;
                logEvidence[a, b] = -len * (half2Pi + 0.5 * logSigma2)
                                    - 0.5 * (Math.Log(t) - logSigma2)
                                    - 0.5 / _sigma2 * (sq - _rho2 / t * sr * sr);
                double v = 1.0 / (invRho2 + len / _sigma2);
                levelVar[a, b] = v;
                levelMean[a, b] = (_nu * invRho2 + (yPre[b] - yPre[a - 1]) / _sigma2) * v;
            }
    }
}

/// <summary>
/// Cauchy (heavy-tailed) noise + Gaussian level prior — the outlier-robust model (Hutter: "Cauchy can handle
/// outliers"). A single gross outlier contributes only boundedly to a Cauchy likelihood, so it neither inflates
/// the segment level nor spawns a spurious one-point segment. No conjugacy, so the single-segment evidence and
/// posterior level moments are computed by 1-D quadrature over the level μ on a shared grid (<c>scale</c> is the
/// Cauchy noise width, the robust analogue of σ) — the per-point log-likelihoods are prefix-summed across the
/// grid, so each segment integral is O(gridPoints).
/// </summary>
public sealed class CauchySegmentModel : ISegmentModel
{
    private readonly double[] _y;
    private readonly double _scale;
    private readonly double _nu;
    private readonly double _rho2;
    private readonly int _grid;

    public CauchySegmentModel(double[] y, double scale, double nu, double rho2, int gridPoints = 512)
    {
        ArgumentNullException.ThrowIfNull(y);
        if (!(scale > 0.0)) throw new ArgumentOutOfRangeException(nameof(scale));
        if (!(rho2 > 0.0)) throw new ArgumentOutOfRangeException(nameof(rho2));
        if (gridPoints < 8) throw new ArgumentOutOfRangeException(nameof(gridPoints));
        _y = y; _scale = scale; _nu = nu; _rho2 = rho2; _grid = gridPoints;
    }

    public int Count => _y.Length;

    public void Fill(double[,] logEvidence, double[,] levelMean, double[,] levelVar)
    {
        int n = _y.Length, g = _grid;
        double span = 6.0 * Math.Sqrt(_rho2);                              // the level lives within the prior
        double lo = _nu - span, hi = _nu + span;
        double dmu = (hi - lo) / (g - 1);
        double logScalePi = Math.Log(Math.PI * _scale);
        double logPriorNorm = -0.5 * Math.Log(2.0 * Math.PI * _rho2);

        var mu = new double[g];
        var logPrior = new double[g];
        var cum = new double[g][];                                          // cum[h][i] = Σ_{j≤i} log Cauchy(y_j|μ_h)
        for (int h = 0; h < g; h++)
        {
            double m = lo + h * dmu;
            mu[h] = m;
            logPrior[h] = logPriorNorm - 0.5 * (m - _nu) * (m - _nu) / _rho2;
            var c = new double[n + 1];
            for (int i = 1; i <= n; i++)
            {
                double z = (_y[i - 1] - m) / _scale;
                c[i] = c[i - 1] + (-logScalePi - Math.Log(1.0 + z * z));
            }
            cum[h] = c;
        }

        for (int a = 1; a <= n; a++)
            for (int b = a; b <= n; b++)
            {
                double maxL = double.NegativeInfinity;
                for (int h = 0; h < g; h++)
                {
                    double l = cum[h][b] - cum[h][a - 1] + logPrior[h];
                    if (l > maxL) maxL = l;
                }
                double s0 = 0.0, s1 = 0.0, s2 = 0.0;
                for (int h = 0; h < g; h++)
                {
                    double l = cum[h][b] - cum[h][a - 1] + logPrior[h];
                    double e = Math.Exp(l - maxL);
                    s0 += e; s1 += mu[h] * e; s2 += mu[h] * mu[h] * e;
                }
                logEvidence[a, b] = maxL + Math.Log(s0 * dmu);             // log ∫ exp(·) dμ
                double mean = s1 / s0;
                levelMean[a, b] = mean;
                levelVar[a, b] = Math.Max(0.0, s2 / s0 - mean * mean);
            }
    }
}
