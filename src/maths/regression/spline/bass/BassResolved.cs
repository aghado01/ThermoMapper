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

    /// <summary>τ's prior on the coefficient form's unit-information precision τ₀ (1/N for the g-form, 1 for the ridge):
    /// sampled, or held at its declared value or τ₀.</summary>
    public required ResolvedTauPrior Scale { get; init; }

    /// <summary>The training points whose fitted values the diagnostics monitor: the config's, a direct response's
    /// missing locations at the middle one, or the defaults.</summary>
    public required IReadOnlyList<BassMonitoredPoint> MonitoredPoints { get; init; }

    /// <summary>The groups an R̂ stop waits on, by name (<see cref="BassSampling.StoppingScope"/>): every group, or the
    /// predictive group alone. The verdict reports every group either way.</summary>
    public required IReadOnlyList<string> StoppingGroups { get; init; }

    /// <summary>
    /// A tempered run's hot end, its lowest positive b: declared, or auto-picked (<see cref="AutoHotEnd"/>); null when the
    /// run is untempered.
    /// </summary>
    public required double? HotEnd { get; init; }

    /// <summary>
    /// The local steps' windows w = max(1, ⌈ω·n_v⌉) over each hinge variable's n_v knots, continuous and functional
    /// (<see cref="LocalSteps.Window"/>); null when the steps are off.
    /// </summary>
    public required BassStepWindows? StepWindows { get; init; }

    /// <summary>The fixed-M factor moves' shares of the change's slot; null when the family is off.</summary>
    public required BassFactorShares? FactorShares { get; init; }

    /// <summary>A tempered run's split threshold (<see cref="SplitRule"/>); null when the rule is off or the run untempered.</summary>
    public required double? SplitThreshold { get; init; }

    /// <summary>Every auto-pick, with its reason.</summary>
    public required IReadOnlyList<BassAutoPick> AutoPicks { get; init; }

    /// <summary>
    /// <paramref name="config"/> resolved on <paramref name="data"/>, after <see cref="BassConfig.Validate"/>. Refused:
    /// a zero response under the zero centre, and a constant one under level centring or the flat intercept, which
    /// leave no spread to fit (D03); a support threshold above its part's size, which no term could meet; a ladder
    /// with a rung where its tempered form is undefined (D05), its hot end, declared or auto-picked, or the reference
    /// rung's b = 0, and an auto-picked hot end outside (0, 1), which tiny data give; and SMC on the integrated form under
    /// the improper options, which has no b = 0 to start from.
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

        BassCoefficientPrior coefficients = config.Coefficients;
        ResolvedTauPrior scale = coefficients.Scale.For(coefficients.UnitObservations(n));
        string unit = coefficients.Form is RidgeForm ? "1, the normalized ridge's" : "1/N, the g-form's";
        switch (coefficients.Scale.Kind)
        {
            case TauPriorKind.UnitInformationGamma:
                log.Add(new("TauPrior", $"Gamma({Text(scale.Gamma!.Value.Shape)}, rate {Text(scale.Gamma.Value.Rate)})", $"the unit-information Gamma, rate a/τ₀, mean τ₀ = {unit}"));
                break;
            case TauPriorKind.HyperG or TauPriorKind.HyperGOverN:
                log.Add(new("TauPrior", scale.ToString(), $"τ₀ = {unit}" + (scale.IsProper ? "" : "; improper, a ≤ 2")));
                break;
            case TauPriorKind.Held when coefficients.Scale.Tau is null:
                log.Add(new("HeldTau", Text(scale.HeldTau!.Value), $"unit information, τ₀ = {unit}"));
                break;
        }

        double? hotEnd = null;
        if (config.Sampling.Tempering is { } tempering)
        {
            if (tempering.Levels is null)
                log.Add(new("Tempering.Levels", $"{tempering.StartLevels()}, self-sizing within {tempering.MaxLevels}",
                            "Levels null: the ladder sizes itself during burn-in by its barrier, N ≈ 2Λ pairs"));
            hotEnd = tempering.HotEnd;
            if (hotEnd is null)
            {
                double? halving = HalvingTau(coefficients.Form, scale, n);
                double picked = AutoHotEnd(tempering, config.Intercept, n, halving);
                if (!(picked > 0.0 && picked < 1.0))
                    throw new ArgumentException($"Sampling.Tempering.HotEnd auto-picks {Text(picked)} on N = {n} observations, outside (0, 1); set it.", "config");
                hotEnd = picked;
                log.Add(new("Tempering.HotEnd", Text(picked), AutoHotEndReason(tempering, config.Intercept, halving, picked)));
            }

            // A ladder's rungs where its tempered form exists (D05): the integrated form loses b = 0 under the improper
            // options, and under the flat intercept every b up to (1 − 2c)/N.
            if (tempering.Path == TemperingPath.Likelihood)
            {
                string improper = config.Intercept == GPriorCentre.Flat ? "the flat intercept" : "the Jeffreys σ² prior";
                bool setHotEnd = tempering.HotEnd is not null;
                if (!BassTarget.Defines(tempering.Form, config.Intercept, sigma2, n, hotEnd.Value))
                    throw new ArgumentException(
                        $"Sampling.Tempering.HotEnd {(setHotEnd ? "is" : "auto-picks")} {Text(hotEnd.Value)}, where the integrated form is undefined under " +
                        $"{improper}, which needs b > (1 − 2c)/N = {Text((1.0 - 2.0 * (sigma2?.Shape ?? 0.0)) / n)}; " +
                        $"{(setHotEnd ? "raise it" : "set it above")}, or take the linear form (D05).", "config");
                if (tempering.ReferenceRung && !BassTarget.Defines(tempering.Form, config.Intercept, sigma2, n, 0.0))
                    throw new ArgumentException(
                        $"Sampling.Tempering.ReferenceRung needs b = 0, where the integrated form is undefined under {improper}; the linear form reaches it (D05).", "config");
            }
        }

        if (config.Sampling.Smc is { } smc)
        {
            // SMC starts from the prior at b = 0, where the integrated form exists under the proper priors alone (D05).
            if (!BassTarget.Defines(smc.Form, config.Intercept, sigma2, n, 0.0))
                throw new ArgumentException(
                    $"Sampling.Smc.Form is integrated, which has no b = 0 member under {(config.Intercept == GPriorCentre.Flat ? "the flat intercept" : "the Jeffreys σ² prior")}, " +
                    "where SMC starts from the prior; the linear form reaches it (D05).", "config");
            log.Add(new("Sampling.Chains", "not read", $"SMC runs Sampling.Smc.Replicates = {smc.Replicates} independent populations"));
            log.Add(new("Sampling.BurnIn, Samples, Thin, BatchSize, MinEssPerChain, StoppingScope", "not read",
                        $"SMC's draws are each replicate's final population of {smc.Particles} particles, and it ends at b = 1, with no R̂ stop"));
        }

        BassStepWindows? windows = null;
        if (config.Sampling.LocalSteps is { } steps)
        {
            // Each hinge variable's knots, as the prior counts them: its distinct scaled values.
            BassInputs inputs = BassInputs.Prepare(data, BassInputTransform.Fit(data));
            int[] Windows(int count, Func<int, int> knots)
            {
                var w = new int[count];
                for (int v = 0; v < count; v++) w[v] = steps.Window(knots(v));
                return w;
            }
            windows = new BassStepWindows(Windows(inputs.Cases.ContinuousCount, v => inputs.ContinuousKnots(v).Length),
                                          Windows(inputs.Outputs?.Count ?? 0, v => inputs.FunctionalKnots(v).Length));
            log.Add(new("LocalSteps windows", windows.ToString(), $"max(1, ⌈ω·n_v⌉) over each hinge variable's n_v knots, ω = {Text(steps.WindowFraction)}"));
        }
        BassFactorShares? shares = config.Sampling.FactorMoves is { } family ? BassFactorShares.Of(family) : null;

        return new BassResolved
        {
            Cases = cases, Locations = locations, Observations = n, MaxBasis = maxBasis,
            MaxInteraction = maxInteraction, MaxCategoricalInteraction = maxCategorical, MaxFunctionalInteraction = maxFunctional,
            SupportMinimum = support, SupportMinimumFunctional = supportFunctional,
            MonitoredPoints = points, ResponseMean = mean, CentreLevel = level, ResponseMeanSquare = meanSquare, Sigma2 = sigma2, Scale = scale,
            StoppingGroups = BassDiagnostics.StoppingGroups(config.Sampling.StoppingScope),
            HotEnd = hotEnd, StepWindows = windows, FactorShares = shares, SplitThreshold = config.Sampling.Tempering?.Split?.Threshold,
            AutoPicks = log.AsReadOnly(),
        };
    }

    /// <summary>
    /// The hot end a ladder auto-picks on <paramref name="observations"/> observations (decided 2026-10-28): 1/N, the b at
    /// which τ₀ makes a unit-information column's shrinkage, τ/(τ + b) under the g-form and τ/(τ + bN) under the ridge,
    /// one half, so the tempered likelihood weighs one observation, as the unit-information prior does; 2/N on the
    /// likelihood path's integrated form under the flat intercept, whose pole at (1 − 2c)/N it clears with ν = Nb − 1 = 1;
    /// and under a held τ the smaller of these and the b at which the held τ halves it (<paramref name="halvingTau"/>).
    /// </summary>
    internal static double AutoHotEnd(TemperingSpec tempering, GPriorCentre intercept, int observations, double? halvingTau)
    {
        double unit = (Pole(tempering, intercept) ? 2.0 : 1.0) / observations;
        return halvingTau is double tau ? Math.Min(tau, unit) : unit;
    }

    /// <summary>
    /// Under a held τ, the b at which it halves a unit-information column's shrinkage: τ under the g-form, τ/N under the
    /// ridge, whose columns carry N observations' information each; null when τ is sampled.
    /// </summary>
    internal static double? HalvingTau(CoefficientForm form, ResolvedTauPrior scale, int observations)
        => scale.HeldTau is double tau ? (form is RidgeForm ? tau / observations : tau) : null;

    private static bool Pole(TemperingSpec tempering, GPriorCentre intercept)
        => tempering.Path == TemperingPath.Likelihood && tempering.Form == TemperedForm.Integrated && intercept == GPriorCentre.Flat;

    private static string AutoHotEndReason(TemperingSpec tempering, GPriorCentre intercept, double? halvingTau, double picked)
    {
        string unit = Pole(tempering, intercept)
            ? "2/N, past the integrated form's pole under the flat intercept at (1 − 2c)/N, with ν = Nb − 1 = 1"
            : "1/N, where τ₀ halves a unit-information column's shrinkage: the tempered likelihood weighs one observation, as the unit-information prior does";
        if (halvingTau is not double tau) return unit;
        return picked == tau ? $"the b at which the held τ halves that shrinkage, at most {unit}" : $"{unit}, below the b at which the held τ halves that shrinkage";
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

/// <summary>The local steps' window w per hinge variable (<see cref="LocalSteps.Window"/>): the continuous predictors', then
/// the output variables'.</summary>
public sealed record BassStepWindows(int[] Continuous, int[] Functional)
{
    /// <inheritdoc/>
    public override string ToString()
        => $"continuous [{string.Join(", ", Continuous)}]" + (Functional.Length > 0 ? $"; functional [{string.Join(", ", Functional)}]" : "");
}

/// <summary>
/// The fixed-M factor moves' shares of the change's slot (<see cref="FactorMoves"/>), each its weight over their sum:
/// a factor's redraw or local step, a prune, a grow (half the prune-or-grow weight each) and a variable change.
/// </summary>
public sealed record BassFactorShares(double Change, double Prune, double Grow, double VariableChange)
{
    /// <summary>The shares <paramref name="family"/>'s weights give.</summary>
    public static BassFactorShares Of(FactorMoves family)
    {
        ArgumentNullException.ThrowIfNull(family);
        double total = family.Change + family.PruneOrGrow + family.VariableChange;
        return new(family.Change / total, 0.5 * family.PruneOrGrow / total, 0.5 * family.PruneOrGrow / total, family.VariableChange / total);
    }
}
