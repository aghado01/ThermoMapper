using System;
using System.Collections.Generic;

namespace Maths.Regression.Conformal;

/// <summary>
/// A betting function of a conformal test martingale (E35, M-2): the next step's bet fₙ on [0, 1], chosen from the p-values before
/// it. A member gives ln fₙ(p) through <see cref="LogBet"/> before pₙ is seen, and <see cref="Observe"/> then records pₙ, which only
/// later bets read. Each bet is ≥ 0 and integrates to at most 1 on [0, 1], exactly 1 where <see cref="IntegratesToOne"/> holds. That
/// makes the capital Sₙ = Π fᵢ(pᵢ) a supermartingale under uniform p-values, a martingale on equality (<c>lintegral_capital</c>,
/// <c>lintegral_capital_le</c>).
/// </summary>
public interface IBettingFunction
{
    /// <summary>ln fₙ(p) for the step about to be taken, a function of p₁, …, pₙ₋₁ only.</summary>
    double LogBet(double p);

    /// <summary>Records pₙ once the bet on it is placed.</summary>
    void Observe(double p);

    /// <summary>Whether every bet integrates to exactly 1, so the capital is a martingale and its value at a stopping time an
    /// e-value with equality; otherwise at most 1.</summary>
    bool IntegratesToOne { get; }
}

/// <summary>The power bet f(p) = εp^(ε−1), ε in (0, 1) (Vovk, Nouretdinov and Gammerman 2003): it reads no history and bets on small
/// p-values, more sharply as ε falls. It integrates to 1.</summary>
public sealed class PowerBet : IBettingFunction
{
    private readonly double _logEpsilon;

    /// <param name="epsilon">ε, in (0, 1).</param>
    public PowerBet(double epsilon)
    {
        if (!(epsilon > 0.0 && epsilon < 1.0)) throw new ArgumentOutOfRangeException(nameof(epsilon), "ε lies in (0, 1).");
        Epsilon = epsilon;
        _logEpsilon = Math.Log(epsilon);
    }

    /// <summary>ε.</summary>
    public double Epsilon { get; }

    /// <inheritdoc/>
    public double LogBet(double p) => _logEpsilon + (Epsilon - 1.0) * Math.Log(p);

    /// <inheritdoc/>
    public void Observe(double p) { }

    /// <inheritdoc/>
    public bool IntegratesToOne => true;
}

/// <summary>
/// The power bets mixed uniformly over ε, Vovk, Nouretdinov and Gammerman's simple mixture: f(p) = ∫₀¹ εp^(ε−1) dε = (eˣ − 1 − x)/x²
/// with x = −ln p, ½ at p = 1. It reads no history and integrates to 1. Below x = 1 it is summed as Σ xʲ/(j + 2)!, whose terms are all
/// positive.
/// </summary>
public sealed class PowerMixtureBet : IBettingFunction
{
    /// <inheritdoc/>
    public double LogBet(double p)
    {
        double x = -Math.Log(p);
        if (x >= 1.0) return Math.Log(Math.Exp(x) - 1.0 - x) - 2.0 * Math.Log(x);   // e^x − 1 − x ≥ e − 2 here: no cancellation
        double term = 0.5, sum = 0.5;   // xʲ/(j + 2)!, from j = 0
        for (int j = 1; j <= 20; j++)
        {
            term *= x / (j + 2);
            sum += term;
        }
        return Math.Log(sum);
    }

    /// <inheritdoc/>
    public void Observe(double p) { }

    /// <inheritdoc/>
    public bool IntegratesToOne => true;
}

/// <summary>
/// Simple Jumper (Vovk, Petej, Nouretdinov, Ahlberg, Carlsson and Gammerman 2021, Algorithm 1): three bets f_ε(p) = 1 + ε(p − ½), ε in
/// {−1, 0, 1}, whose capitals pass a share J of their sum among themselves before each step, as a Markov chain jumping at rate J. The
/// step's bet is their mixture, f(p) = 1 + (w₁ − w₋₁)(p − ½), with the weights w_ε from the past. It integrates to 1 and lies in [½, 3/2].
/// </summary>
public sealed class SimpleJumperBet : IBettingFunction
{
    private readonly double _jump;
    private double _down = 1.0 / 3.0, _flat = 1.0 / 3.0, _up = 1.0 / 3.0;

    /// <param name="jump">J, the share passed at each step, in [0, 1]; Vovk et al.'s 0.01 by default.</param>
    public SimpleJumperBet(double jump = 0.01)
    {
        if (!(jump >= 0.0 && jump <= 1.0)) throw new ArgumentOutOfRangeException(nameof(jump), "J lies in [0, 1].");
        _jump = jump;
    }

    /// <inheritdoc/>
    public double LogBet(double p) => Math.Log(1.0 + (_up - _down) * (p - 0.5));

    /// <inheritdoc/>
    public void Observe(double p)
    {
        double f = 1.0 + (_up - _down) * (p - 0.5);
        _down *= (1.0 - (p - 0.5)) / f;
        _up *= (1.0 + (p - 0.5)) / f;
        _flat = 1.0 - _down - _up;
        _down = (1.0 - _jump) * _down + _jump / 3.0;   // the next step's jump
        _flat = (1.0 - _jump) * _flat + _jump / 3.0;
        _up = (1.0 - _jump) * _up + _jump / 3.0;
    }

    /// <inheritdoc/>
    public bool IntegratesToOne => true;
}

/// <summary>
/// The plug-in bet (Fedorova, Gammerman, Nouretdinov and Vovk 2012): a kernel density of the p-values so far, reflected at 0 and 1, as
/// the next step's bet, which asymptotically does at least as well as any power bet.
/// <list type="bullet">
/// <item><b>The density.</b> The p-values are counted in equal bins, and a Gaussian kernel of the bins, with both reflections, is held at
/// the bins' centres, linear between them and flat beyond the outer ones. Its bandwidth is Silverman's rule, 0.9·min(sd, IQR/1.34)·n^(−1/5),
/// from the binned values, at least a bin. It is refitted when the count reaches the warm-up and each time it doubles; between refits
/// each new value adds its bin's kernel row, O(bins).</item>
/// <item><b>The bet.</b> The density over its integral, exact for the piecewise-linear form, scaled down by 10⁻¹² so rounding cannot lift
/// the integral above 1 (Ville needs only ∫ f ≤ 1), then mixed with the uniform at weight <c>floor</c>, so a bet never reaches 0. Before
/// the warm-up it bets 1.</item>
/// </list>
/// </summary>
public sealed class PlugInBet : IBettingFunction
{
    private const double Shrink = 1.0 - 1e-12;
    private readonly int _bins, _warmUp;
    private readonly double _floor;
    private readonly int[] _counts;
    private readonly double[] _density;
    private double[,] _rows = new double[0, 0];
    private double[] _rowIntegrals = [];
    private double _integral;
    private int _count, _nextRefit;

    /// <param name="warmUp">The p-values seen before the first kernel bet; until then each bet is 1.</param>
    /// <param name="bins">The bins the p-values are counted in, and the density's grid.</param>
    /// <param name="floor">The uniform's weight in the bet, in [0, 1): 0 is the paper's pure kernel bet.</param>
    public PlugInBet(int warmUp = 10, int bins = 100, double floor = 0.0)
    {
        if (warmUp < 2) throw new ArgumentOutOfRangeException(nameof(warmUp), "The warm-up needs at least two p-values for a spread.");
        if (bins < 2) throw new ArgumentOutOfRangeException(nameof(bins));
        if (!(floor >= 0.0 && floor < 1.0)) throw new ArgumentOutOfRangeException(nameof(floor), "The uniform's weight lies in [0, 1).");
        _warmUp = warmUp;
        _bins = bins;
        _floor = floor;
        _counts = new int[bins];
        _density = new double[bins];
        _nextRefit = warmUp;
    }

    /// <summary>The kernel's bandwidth at the last refit; NaN before the warm-up.</summary>
    public double Bandwidth { get; private set; } = double.NaN;

    /// <inheritdoc/>
    public double LogBet(double p)
    {
        if (_count < _warmUp) return 0.0;
        double u = p * _bins - 0.5, d;
        if (u <= 0.0) d = _density[0];
        else if (u >= _bins - 1) d = _density[_bins - 1];
        else
        {
            int g = (int)u;
            double t = u - g;
            d = (1.0 - t) * _density[g] + t * _density[g + 1];
        }
        return Math.Log(_floor + (1.0 - _floor) * Shrink * d / _integral);
    }

    /// <inheritdoc/>
    public void Observe(double p)
    {
        int bin = Math.Clamp((int)(p * _bins), 0, _bins - 1);
        _counts[bin]++;
        _count++;
        if (_count == _nextRefit)
        {
            Refit();
            _nextRefit = 2 * _count;
        }
        else if (_count > _warmUp)
        {
            for (int g = 0; g < _bins; g++) _density[g] += _rows[bin, g];
            _integral += _rowIntegrals[bin];
        }
    }

    /// <inheritdoc/>
    public bool IntegratesToOne => false;

    // Silverman's bandwidth from the bins, each bin's kernel row at the grid, and the density and its integral from the counts.
    private void Refit()
    {
        double width = 1.0 / _bins, mean = 0.0, square = 0.0;
        for (int b = 0; b < _bins; b++)
        {
            double c = (b + 0.5) * width;
            mean += _counts[b] * c;
            square += _counts[b] * c * c;
        }
        mean /= _count;
        double sd = Math.Sqrt(Math.Max(0.0, square / _count - mean * mean) * _count / (_count - 1.0));
        double iqr = Quantile(0.75) - Quantile(0.25), spread = iqr > 0.0 ? Math.Min(sd, iqr / 1.34) : sd;
        double h = Math.Max(width, 0.9 * spread * Math.Pow(_count, -0.2));
        Bandwidth = h;
        _rows = new double[_bins, _bins];
        _rowIntegrals = new double[_bins];
        for (int b = 0; b < _bins; b++)
        {
            double c = (b + 0.5) * width;
            for (int g = 0; g < _bins; g++)
            {
                double x = (g + 0.5) * width;
                _rows[b, g] = Gauss((x - c) / h) + Gauss((x + c) / h) + Gauss((x - 2.0 + c) / h);
            }
            _rowIntegrals[b] = Integral(b);
        }
        Array.Clear(_density);
        _integral = 0.0;
        for (int b = 0; b < _bins; b++)
        {
            if (_counts[b] == 0) continue;
            for (int g = 0; g < _bins; g++) _density[g] += _counts[b] * _rows[b, g];
            _integral += _counts[b] * _rowIntegrals[b];
        }
    }

    private static double Gauss(double z) => Math.Exp(-0.5 * z * z);

    // ∫₀¹ of row b's piecewise-linear form: flat to the first centre and from the last, trapezoids between.
    private double Integral(int b)
    {
        double width = 1.0 / _bins, sum = 0.5 * width * (_rows[b, 0] + _rows[b, _bins - 1]);
        for (int g = 0; g + 1 < _bins; g++) sum += 0.5 * width * (_rows[b, g] + _rows[b, g + 1]);
        return sum;
    }

    // The binned values' quantile, by the bins' centres.
    private double Quantile(double q)
    {
        double target = q * _count;
        int seen = 0;
        for (int b = 0; b < _bins; b++)
        {
            seen += _counts[b];
            if (seen >= target) return (b + 0.5) / _bins;
        }
        return 1.0 - 0.5 / _bins;
    }
}

/// <summary>One step of the capital: the bet's log at the p-value, ln fₙ(pₙ), which shows the observation's pull on the martingale, and
/// the log capital after it, ln Sₙ.</summary>
public readonly record struct CapitalStep(int Count, double P, double LogBet, double LogCapital);

/// <summary>
/// The capital of a betting function over a sequence of p-values (E35, M-2): Sₙ = Π fᵢ(pᵢ) with S₀ = 1, held in logs. Each bet is placed
/// before its p-value is observed. Under i.i.d. uniform p-values the capital is a martingale if each bet integrates to 1, a supermartingale
/// if at most 1 (<c>lintegral_capital</c>), and its value at a stopping time that reads only the past is an e-value
/// (<c>lintegral_stoppedCapital</c>). Not thread-safe.
/// </summary>
public sealed class BettingCapital
{
    private readonly IBettingFunction _bet;

    public BettingCapital(IBettingFunction bet) => _bet = bet ?? throw new ArgumentNullException(nameof(bet));

    /// <summary>The p-values bet on so far.</summary>
    public int Count { get; private set; }

    /// <summary>ln Sₙ, 0 before the first p-value.</summary>
    public double LogCapital { get; private set; }

    /// <summary>Bets on <paramref name="p"/>, in (0, 1], then lets the betting function observe it.</summary>
    public CapitalStep Add(double p)
    {
        if (!(p > 0.0 && p <= 1.0)) throw new ArgumentOutOfRangeException(nameof(p), "A p-value lies in (0, 1].");
        double logBet = _bet.LogBet(p);
        if (double.IsNaN(logBet) || double.IsPositiveInfinity(logBet))
            throw new InvalidOperationException($"The bet at p = {p:R} is {Math.Exp(logBet)}; a bet is finite and at least 0.");
        _bet.Observe(p);
        Count++;
        LogCapital += logBet;
        return new CapitalStep(Count, p, logBet, LogCapital);
    }
}

/// <summary>The procedure that reads a conformal test martingale's capital for an alarm (Vovk, <i>Statistical Science</i> 2021).</summary>
public enum AlarmProcedure
{
    /// <summary>Ville's: Sₙ/S_τ ≥ c, the capital since the last alarm τ (0 at the start). Under exchangeability the probability of any
    /// false alarm in a run is at most 1/c, proved (<c>ville</c>).</summary>
    Ville,

    /// <summary>CUSUM: max over τ ≤ i &lt; n of Sₙ/Sᵢ ≥ c, by Page's recursion γₙ = fₙ·max(γₙ₋₁, 1). It reads past a capital lost before a
    /// change. The mean time to a false alarm is at least c, and, repeated, the long-run frequency of false alarms at most 1/c: posed
    /// (Vovk 2021, Propositions 4.2 and 4.4), not proved.</summary>
    Cusum,

    /// <summary>Shiryaev–Roberts: the sum over τ ≤ i &lt; n of Sₙ/Sᵢ ≥ c, by Rₙ = fₙ·(Rₙ₋₁ + 1). Its guarantees are CUSUM's, posed, not
    /// proved; CUSUM alarms no more often.</summary>
    ShiryaevRoberts,
}

/// <summary>What an alarm guarantees under exchangeability, in its procedure's terms, and whether it is proved (M-2) or posed.</summary>
public sealed record AlarmGuarantee(string Statement, bool Proved);

/// <summary>
/// An alarm on a conformal test martingale's capital: its procedure (<see cref="AlarmProcedure"/>) and threshold c &gt; 1, its statistic in
/// logs, and the arrivals at which it fired. Each statistic reads only the capital so far. After each alarm it restarts from the capital
/// then, so each alarm tests the stream since the last; its guarantee is per run. A capital lost to a zero bet stays lost for Ville's,
/// and CUSUM and Shiryaev–Roberts read past it from the next step. Not thread-safe.
/// </summary>
public sealed class MartingaleAlarm
{
    private readonly double _logThreshold;
    private readonly List<int> _times = [];
    private int _count;

    /// <param name="procedure">The procedure.</param>
    /// <param name="threshold">c, above 1.</param>
    public MartingaleAlarm(AlarmProcedure procedure, double threshold)
    {
        if (!(threshold > 1.0 && double.IsFinite(threshold))) throw new ArgumentOutOfRangeException(nameof(threshold), "The threshold c is finite and above 1.");
        if (!Enum.IsDefined(procedure)) throw new ArgumentOutOfRangeException(nameof(procedure));
        Procedure = procedure;
        Threshold = threshold;
        _logThreshold = Math.Log(threshold);
        LogStatistic = Start;
    }

    /// <summary>The procedure.</summary>
    public AlarmProcedure Procedure { get; }

    /// <summary>c.</summary>
    public double Threshold { get; }

    /// <summary>The statistic's log after the latest step: ln(Sₙ/S_τ), ln γₙ or ln Rₙ.</summary>
    public double LogStatistic { get; private set; }

    /// <summary>The arrivals at which it fired, counted from 1.</summary>
    public IReadOnlyList<int> AlarmTimes => _times;

    /// <summary>Its guarantee under exchangeability, in its procedure's terms.</summary>
    public AlarmGuarantee Guarantee => Procedure == AlarmProcedure.Ville
        ? new AlarmGuarantee($"The probability of any false alarm in a run since the last alarm is at most 1/c = {1.0 / Threshold:G4} (Ville's inequality).", true)
        : new AlarmGuarantee($"The mean time to a false alarm is at least c = {Threshold:G4}, and repeated, the long-run frequency of false alarms is at most 1/c "
                             + $"= {1.0 / Threshold:G4} (Vovk 2021, Propositions 4.2 and 4.4): posed, not proved.", false);

    // Ville's run starts at S/S = 1; CUSUM's and Shiryaev–Roberts's at γ₀ = R₀ = 0.
    private double Start => Procedure == AlarmProcedure.Ville ? 0.0 : double.NegativeInfinity;

    /// <summary>Advances by one step's ln fₙ(pₙ); true when the alarm fires at this step, after which it restarts.</summary>
    public bool Step(double logBet)
    {
        if (double.IsNaN(logBet) || double.IsPositiveInfinity(logBet)) throw new ArgumentException("A bet's log is finite or −∞.", nameof(logBet));
        _count++;
        double previous = LogStatistic;
        LogStatistic = Procedure switch
        {
            AlarmProcedure.Ville => previous + logBet,
            AlarmProcedure.Cusum => logBet + Math.Max(previous, 0.0),
            _ => logBet + SoftPlus(previous),
        };
        if (!(LogStatistic >= _logThreshold)) return false;
        _times.Add(_count);
        LogStatistic = Start;
        return true;
    }

    // ln(1 + eˣ), with ln(1 + 0) = 0 at x = −∞.
    private static double SoftPlus(double x) => x > 0.0 ? x + Math.Log(1.0 + Math.Exp(-x)) : Math.Log(1.0 + Math.Exp(x));
}

/// <summary>One arrival of a conformal test martingale: its smoothed p-value, the bet's log at it, the log capital after it, and the alarms
/// that fired at it (empty when none did).</summary>
public readonly record struct MartingaleStep(SmoothedPValue PValue, double LogBet, double LogCapital, IReadOnlyList<AlarmProcedure>? Fired = null);

/// <summary>
/// A conformal test martingale (E35; Vovk, Nouretdinov and Gammerman 2003; Vovk et al. 2021): each arrival's conformity score is ranked
/// into its smoothed p-value (<see cref="SmoothedPValues"/>, M-1), and a betting function turns the p-values into capital
/// (<see cref="BettingCapital"/>, M-2). Under exchangeability of the scored units the capital is a (super)martingale with expectation
/// at most 1, so its value at any stopping time that reads only the past is an e-value (<see cref="EValue"/>), and growth is evidence
/// against exchangeability. Alarms read the capital (<see cref="MartingaleAlarm"/>): Ville's, proved, and CUSUM and Shiryaev–Roberts, whose
/// guarantees are posed. Not thread-safe: one per stream.
/// </summary>
public sealed class ConformalTestMartingale
{
    private static readonly AlarmProcedure[] None = [];
    private readonly SmoothedPValues _pValues;
    private readonly BettingCapital _capital;
    private readonly List<double> _path = [];
    private readonly MartingaleAlarm[] _alarms;

    /// <param name="bet">The betting function, fresh for this stream.</param>
    /// <param name="seed">The tie-breaking draws' seed.</param>
    /// <param name="alarms">The alarms that read the capital, each fresh for this stream.</param>
    public ConformalTestMartingale(IBettingFunction bet, int seed, params MartingaleAlarm[] alarms)
    {
        _capital = new BettingCapital(bet);
        _pValues = new SmoothedPValues(seed);
        ArgumentNullException.ThrowIfNull(alarms);
        _alarms = [.. alarms];
        foreach (MartingaleAlarm alarm in _alarms) ArgumentNullException.ThrowIfNull(alarm, nameof(alarms));
    }

    /// <summary>The alarms, with their statistics, times and guarantees.</summary>
    public IReadOnlyList<MartingaleAlarm> Alarms => _alarms;

    /// <summary>The arrivals so far.</summary>
    public int Count => _capital.Count;

    /// <summary>ln Sₙ.</summary>
    public double LogCapital => _capital.LogCapital;

    /// <summary>Sₙ, the e-value at the current time when the time was fixed in advance or chosen by a rule that reads only the past.</summary>
    public double EValue => Math.Exp(_capital.LogCapital);

    /// <summary>ln S₁, …, ln Sₙ.</summary>
    public IReadOnlyList<double> LogPath => _path;

    /// <summary>Ranks the arrival's <paramref name="score"/> (higher conforming better), bets on its p-value and gives the step.</summary>
    public MartingaleStep Add(double score)
    {
        SmoothedPValue p = _pValues.Add(score);
        CapitalStep step = _capital.Add(p.P);
        _path.Add(step.LogCapital);
        List<AlarmProcedure>? fired = null;
        foreach (MartingaleAlarm alarm in _alarms)
            if (alarm.Step(step.LogBet)) (fired ??= []).Add(alarm.Procedure);
        return new MartingaleStep(p, step.LogBet, step.LogCapital, fired is null ? None : fired);
    }
}
