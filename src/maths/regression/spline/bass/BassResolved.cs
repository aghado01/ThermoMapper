using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Maths.Distributions;
using Maths.LinAlg;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Samplers;

namespace Maths.Regression.Spline.Bass;

/// <summary>One auto-pick a fit made at entry: the setting, the value it took, and why (CONTRIBUTING's automation log).</summary>
public sealed record BassAutoPick(string Setting, string Value, string Reason);

/// <summary>
/// What a config resolves to on its data, at entry: the counts, the caps and thresholds left to the data, and the
/// priors read from the response. With <see cref="AutoPicks"/>, the log of every auto-pick and its reason, it is the
/// fit's account of what it decided beyond the config. The data-dependent refusals are made here too: a response the
/// intercept convention cannot fit (D03), and a support threshold no term could meet.
/// </summary>
public sealed record BassResolved
{
    /// <summary>n_x, the cases.</summary>
    public required int Cases { get; init; }

    /// <summary>n_r, the output locations; 1 for a scalar response.</summary>
    public required int Locations { get; init; }

    /// <summary>N = n_x · n_r, the observations τ's unit information counts (D03).</summary>
    public required int Observations { get; init; }

    /// <summary>The cap on M.</summary>
    public required int MaxBasis { get; init; }

    /// <summary>The caps on a term's continuous, categorical and functional factors.</summary>
    public required int MaxInteraction { get; init; }

    /// <summary>The cap on a term's categorical factors.</summary>
    public required int MaxCategoricalInteraction { get; init; }

    /// <summary>The cap on a term's functional factors; 0 for a scalar response.</summary>
    public required int MaxFunctionalInteraction { get; init; }

    /// <summary>The support threshold on a term's design part.</summary>
    public required double SupportMinimum { get; init; }

    /// <summary>The support threshold on a term's functional part; null for a scalar response.</summary>
    public required double? SupportMinimumFunctional { get; init; }

    /// <summary>The level the coefficient prior is centred on: the response's mean, or zero under the zero centre.</summary>
    public required double CentreLevel { get; init; }

    /// <summary>ȳ, the response's mean over all N observations, whatever the centre; the zero centre's level energy
    /// Nȳ² reads it.</summary>
    public required double ResponseMean { get; init; }

    /// <summary>s², the response's mean square about <see cref="CentreLevel"/>, which the scaled σ² prior reads.</summary>
    public required double ResponseMeanSquare { get; init; }

    /// <summary>σ²'s prior as resolved; null for the Jeffreys prior.</summary>
    public required InverseGamma? Sigma2 { get; init; }

    /// <summary>τ's prior, Gamma(a, rate a·N), when τ is sampled; null when it is fixed.</summary>
    public required Gamma? TauPrior { get; init; }

    /// <summary>τ's value when it is fixed; null when it is sampled.</summary>
    public required double? FixedTau { get; init; }

    /// <summary>The training points whose fitted values the diagnostics monitor: the config's, a direct response's
    /// missing locations at the middle one, or the defaults.</summary>
    public required IReadOnlyList<BassMonitoredPoint> MonitoredPoints { get; init; }

    /// <summary>Every auto-pick, with its reason.</summary>
    public required IReadOnlyList<BassAutoPick> AutoPicks { get; init; }

    /// <summary>
    /// <paramref name="config"/> resolved on <paramref name="data"/>, after <see cref="BassConfig.Validate"/>. Refused:
    /// a zero response under the zero centre, and a constant one under level centring or the flat intercept, which
    /// leave no spread to fit (D03); a support threshold above its part's size, which no term could meet; and a ladder
    /// with a rung where its tempered form is undefined (D05), its hot end or the reference rung's b = 0.
    /// </summary>
    public static BassResolved Resolve(BassConfig config, BassData data)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(data);
        config.Validate();
        var log = new List<BassAutoPick>();
        BassResponse response = data.Response;
        int cases = response.Cases, locations = response.Locations, n = response.Observations;
        bool functional = response is FunctionalResponse;
        log.Add(new("Observations", Text(n), functional ? $"N = n_x · n_r = {cases} · {locations}" : "N = n_x, a scalar response"));

        // The response the intercept convention can fit, its level and its mean square about it.
        ReadOnlySpan<double> y = response.Values;
        bool constant = true, zero = y[0] == 0.0;
        for (int i = 1; i < y.Length && constant; i++) constant = y[i] == y[0];
        double mean = PairwiseSum.Of(y) / n, level;
        if (config.Intercept == GPriorCentre.Zero)
        {
            if (constant && zero) throw new ArgumentException("The response is zero, which leaves σ no posterior (D03).", nameof(data));
            if (constant && !config.Sigma2.IsProper)
                throw new ArgumentException("Under the Jeffreys σ² prior a constant response is fitted exactly by the intercept, which leaves the posterior improper (D03).", nameof(data));
            level = 0.0;
            log.Add(new("CentreLevel", "0", "the zero centre"));
        }
        else
        {
            if (constant) throw new ArgumentException("The response is constant, which the intercept absorbs whole, leaving no spread to fit (D03).", nameof(data));
            level = mean;
            log.Add(new("CentreLevel", Text(level), config.Intercept == GPriorCentre.Flat ? "the response's mean, the flat intercept's estimate" : "the response's mean, level centring"));
        }
        var squares = new double[n];
        for (int i = 0; i < n; i++) squares[i] = (y[i] - level) * (y[i] - level);
        double meanSquare = PairwiseSum.Of(squares) / n;
        InverseGamma? sigma2 = config.Sigma2.For(meanSquare);
        if (config.Sigma2.ReadsResponse)
            log.Add(new("ResponseMeanSquare", Text(meanSquare), "s², the response's mean square about the centre level, which the scaled σ² prior reads"));

        // A ladder's rungs where its tempered form exists (D05): the integrated form loses b = 0 under the improper options,
        // and under the flat intercept every b up to (1 − 2c)/N.
        if (config.Sampling.Tempering is { Path: TemperingPath.Likelihood } tempering)
        {
            string improper = config.Intercept == GPriorCentre.Flat ? "the flat intercept" : "the Jeffreys σ² prior";
            if (!BassTarget.Defines(tempering.Form, config.Intercept, sigma2, n, tempering.HotEnd))
                throw new ArgumentException(
                    $"Sampling.Tempering.HotEnd is {Text(tempering.HotEnd)}, where the integrated form is undefined under {improper}, which needs " +
                    $"b > (1 − 2c)/N = {Text((1.0 - 2.0 * (sigma2?.Shape ?? 0.0)) / n)}; raise it, or take the linear form (D05).", "config");
            if (tempering.ReferenceRung && !BassTarget.Defines(tempering.Form, config.Intercept, sigma2, n, 0.0))
                throw new ArgumentException(
                    $"Sampling.Tempering.ReferenceRung needs b = 0, where the integrated form is undefined under {improper}; the linear form reaches it (D05).", "config");
        }

        // The caps, by the predictors of each kind.
        int continuous = 0, categorical = 0;
        foreach (BassPredictor p in data.Predictors)
        {
            if (p is ContinuousPredictor) continuous++;
            else categorical++;
            if (p is CategoricalPredictor factor)
                for (int c = 0; c < factor.Levels.Count; c++)
                    if (!factor.Observed[c])
                        log.Add(new($"Levels of '{factor.Name}'", $"'{factor.Levels[c]}' left out", "no training case takes it"));
        }
        int outputs = response is FunctionalResponse f ? f.OutputVariables.Count : 0;
        BassStructureSpec s = config.Structure;
        int maxInteraction = Cap("MaxInteraction", s.MaxInteraction, continuous, "continuous predictors", log);
        int maxCategorical = Cap("MaxCategoricalInteraction", s.MaxCategoricalInteraction, categorical, "categorical predictors", log);
        int maxFunctional = Cap("MaxFunctionalInteraction", s.MaxFunctionalInteraction, outputs, "output variables", log);

        int bound = n - 1;   // the intercept counts toward the dimension bound
        int maxBasis;
        if (s.MaxBasis is { } declared)
        {
            maxBasis = Math.Min(declared, bound);
            if (declared > bound) log.Add(new("MaxBasis", Text(maxBasis), $"the declared {declared} exceeds N − 1, the bound with the intercept"));
        }
        else
        {
            maxBasis = Math.Min(1000, bound);
            log.Add(new("MaxBasis", Text(maxBasis), "min(1000, N − 1): R's 1,000, with the intercept counting toward the dimension bound"));
        }

        SupportRule rule = s.Support;
        double support = Threshold(rule.Minimum, cases, "Structure.Support.Minimum", "SupportMinimum", "n_x", log);
        double? supportFunctional = functional
            ? Threshold(rule.MinimumFunctional, locations, "Structure.Support.MinimumFunctional", "SupportMinimumFunctional", "n_r", log)
            : null;

        IReadOnlyList<BassMonitoredPoint> points = ResolvePoints(config.Sampling.MonitoredPoints, cases, functional ? locations : null, log);

        Gamma? tauPrior = null;
        double? fixedTau = null;
        switch (config.Coefficients)
        {
            case SampledG g:
                tauPrior = Gamma.ShapeRate(g.TauShape, g.TauShape * n);
                log.Add(new("TauPrior", $"Gamma({Text(g.TauShape)}, rate {Text(g.TauShape * n)})", "the unit-information family, rate a·N, mean 1/N"));
                break;
            case FixedG { Tau: { } tau }:
                fixedTau = tau;
                break;
            case FixedG:
                fixedTau = 1.0 / n;
                log.Add(new("FixedTau", Text(1.0 / n), "unit information, 1/N"));
                break;
        }

        return new BassResolved
        {
            Cases = cases, Locations = locations, Observations = n, MaxBasis = maxBasis,
            MaxInteraction = maxInteraction, MaxCategoricalInteraction = maxCategorical, MaxFunctionalInteraction = maxFunctional,
            SupportMinimum = support, SupportMinimumFunctional = supportFunctional,
            MonitoredPoints = points, ResponseMean = mean, CentreLevel = level, ResponseMeanSquare = meanSquare, Sigma2 = sigma2, TauPrior = tauPrior, FixedTau = fixedTau,
            AutoPicks = log.AsReadOnly(),
        };
    }

    // The declared points, checked against the data, or the defaults; a direct response's missing location is the middle one.
    private static IReadOnlyList<BassMonitoredPoint> ResolvePoints(IReadOnlyList<BassMonitoredPoint>? declared, int cases, int? locations,
                                                                     List<BassAutoPick> log)
    {
        if (declared is null)
        {
            IReadOnlyList<BassMonitoredPoint> defaults = BassDiagnostics.DefaultPoints(cases, locations);
            log.Add(new("MonitoredPoints", string.Join("; ", defaults), "training cases at the quartile positions of the data's order" +
                        (locations is null ? "" : ", at the middle output location")));
            return defaults;
        }
        var points = new BassMonitoredPoint[declared.Count];
        for (int k = 0; k < points.Length; k++)
        {
            BassMonitoredPoint p = declared[k];
            if (p.Case >= cases) throw new ArgumentException($"Sampling.MonitoredPoints holds case {p.Case}; the data have {cases}.", "config");
            if (locations is not int m)
            {
                if (p.Location is not null) throw new ArgumentException("Sampling.MonitoredPoints gives a location, but the response is scalar.", "config");
                points[k] = p;
                continue;
            }
            if (p.Location >= m) throw new ArgumentException($"Sampling.MonitoredPoints holds location {p.Location}; the response has {m}.", "config");
            points[k] = p with { Location = p.Location ?? (m - 1) / 2 };
        }
        if (points.Distinct().Count() != points.Length)
            throw new ArgumentException("Sampling.MonitoredPoints names one point twice once the middle location fills a missing one.", "config");
        return points;
    }

    // A cap at most the count of its variables, as R caps maxInt by pdes, pcat and pfunc.
    private static int Cap(string setting, int declared, int available, string what, List<BassAutoPick> log)
    {
        if (declared <= available) return declared;
        log.Add(new(setting, Text(available), $"the declared {declared} exceeds the {available} {what}"));
        return available;
    }

    // R's default min(20, 0.1·n) for a part of size n, or the declared threshold, which must be one some term can meet.
    private static double Threshold(double? declared, int size, string field, string setting, string sizeName, List<BassAutoPick> log)
    {
        if (declared is { } d)
        {
            if (d > size) throw new ArgumentException($"{field} is {Text(d)}, above the {size} a term's part can reach, so no term could meet it.", "config");
            return d;
        }
        double resolved = Math.Min(20.0, 0.1 * size);
        log.Add(new(setting, Text(resolved), $"min(20, 0.1·{sizeName}), R's default"));
        return resolved;
    }

    private static string Text(double value) => value.ToString("R", CultureInfo.InvariantCulture);
    private static string Text(int value) => value.ToString(CultureInfo.InvariantCulture);
}
