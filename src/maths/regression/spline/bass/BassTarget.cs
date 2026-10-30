using System;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Samplers;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A structure's account for the collapsed target: M, the sum of its terms' ln π (D01), and its least-squares fit in
/// the centred basis, ‖r‖² and Q_c, with, under the normalized ridge, its columns' ridge terms (a spectrum, or a bordered
/// candidate of one); or the design's refusal, which puts it outside F.
/// </summary>
public readonly record struct BassStructureFit(int Terms, double LogTermPrior, double ResidualSquares, double FittedSquares,
                                               BassRefusal? Refusal = null, RidgeColumns? Ridge = null)
{
    /// <summary>The design's current structure, whose terms' log prior is <paramref name="logTermPrior"/>.</summary>
    public static BassStructureFit Of(BassDesign design, double logTermPrior)
    {
        ArgumentNullException.ThrowIfNull(design);
        return new(design.Count, logTermPrior, design.ResidualSquares, design.FittedSquares, null, design.IsRidge ? design.Spectrum : null);
    }

    /// <summary>A proposed structure, admitted or refused, whose terms' log prior is <paramref name="logTermPrior"/>.</summary>
    public static BassStructureFit Of(DesignProposal proposal, double logTermPrior) => proposal.Refusal is { } reason
        ? Refused(reason)
        : new(proposal.Columns, logTermPrior, proposal.ResidualSquares, proposal.FittedSquares, null, proposal.Ridge);

    /// <summary>A structure outside F for <paramref name="reason"/>.</summary>
    public static BassStructureFit Refused(BassRefusal reason) => new(-1, double.NaN, double.NaN, double.NaN, reason);
}

/// <summary>
/// BASS's collapsed target (D03, D04, D05; stage 4's S5): the posterior over (structure, τ), with the coefficients and
/// σ² integrated out.
/// <list type="bullet">
/// <item><b>The prior,</b> untempered on the likelihood path: D02's count weight and D01's term measure
/// (<see cref="BassStructurePrior"/>), and τ's prior from the shared family (<see cref="TauPrior"/>) on the coefficient
/// form's τ₀, as a density on u = ln τ. That is the coordinate τ's slice step moves (D04), so the powered path powers τ
/// in log τ (D05's recorded choice). A held τ is no coordinate of the state.</item>
/// <item><b>The likelihood,</b> p(y | structure, τ) from the shared block of the coefficient form: the g-form's on the
/// design's centred statistics (<see cref="CollapsedGPrior.CentredMarginalAt"/>), or the normalized ridge's on the
/// structure's ridge terms (<see cref="CollapsedGPrior.CentredRidgeMarginalAt"/>), under the intercept's convention (level
/// centring, Liang's flat intercept or R's zero centre) and σ²'s prior (scaled to the response, Jeffreys, or a declared
/// inverse gamma). On the likelihood path at b it takes D05's integrated form p_b, the coefficients and σ² integrated
/// against the tempered likelihood, or the linear form b·ln p₁.</item>
/// <item><b>The support,</b> F: the design's admission (D06) and, under the Jeffreys σ² prior, no exact fit,
/// ‖r‖ ≤ 10³·ε·‖y − level‖ (D03's propriety rule).</item>
/// </list>
/// Under the proper defaults the values are complete, so a ladder's evidence reads them. Under Jeffreys or the flat
/// intercept they hold up to that prior's constant, which cancels from every within-rung ratio and every swap.
/// Immutable: one serves every chain.
/// </summary>
public sealed class BassTarget : ISamplingTarget<BassState>
{
    private const double ExactFitMultiple = 1e3;
    private const double MachineEpsilon = 2.220446049250313e-16;   // 2⁻⁵², not double.Epsilon, the least subnormal

    private readonly double _exactFitSquares;   // ‖r‖² at or below which Jeffreys refuses; −∞ under a proper prior

    /// <summary>The target <paramref name="config"/> sets on its data (<paramref name="resolved"/>), over
    /// <paramref name="prior"/>'s structures, tempered in <paramref name="form"/>.</summary>
    public BassTarget(BassConfig config, BassResolved resolved, BassStructurePrior prior, TemperedForm form = TemperedForm.Integrated)
    {
        ArgumentNullException.ThrowIfNull(config);
        ArgumentNullException.ThrowIfNull(resolved);
        ArgumentNullException.ThrowIfNull(prior);
        if (form is not (TemperedForm.Integrated or TemperedForm.Linear)) throw new ArgumentOutOfRangeException(nameof(form));
        if (config.Intercept is not (GPriorCentre.Mean or GPriorCentre.Flat or GPriorCentre.Zero))
            throw new ArgumentException("BASS's intercept is level-centred, flat or at zero (D03).", nameof(config));
        if (resolved.Scale is null) throw new ArgumentException("The resolved config holds no τ prior.", nameof(resolved));
        Prior = prior;
        Form = form;
        CoefficientForm = config.Coefficients.Form;
        Intercept = config.Intercept;
        Sigma2 = resolved.Sigma2;
        Scale = resolved.Scale;
        Observations = resolved.Observations;
        ResponseMean = resolved.ResponseMean;
        double floor = ExactFitMultiple * MachineEpsilon;
        _exactFitSquares = Sigma2 is null ? floor * floor * Observations * resolved.ResponseMeanSquare : double.NegativeInfinity;
    }

    /// <summary>The prior on structures.</summary>
    public BassStructurePrior Prior { get; }

    /// <summary>The tempered form on the likelihood path.</summary>
    public TemperedForm Form { get; }

    /// <summary>The coefficient prior's covariance form: the g-form or the normalized ridge.</summary>
    public CoefficientForm CoefficientForm { get; }

    /// <summary>The intercept's convention.</summary>
    public GPriorCentre Intercept { get; }

    /// <summary>σ²'s prior as resolved; null is the Jeffreys 1/σ².</summary>
    public InverseGamma? Sigma2 { get; }

    /// <summary>τ's prior on the coefficient form's τ₀, sampled or held.</summary>
    public ResolvedTauPrior Scale { get; }

    /// <summary>The held τ; null when τ is sampled.</summary>
    public double? HeldTau => Scale.HeldTau;

    /// <summary>Whether τ is a coordinate of the state: sampled, not held.</summary>
    public bool SamplesTau => !Scale.IsHeld;

    /// <summary>N, the observations.</summary>
    public int Observations { get; }

    /// <summary>ȳ, which the zero centre's level energy reads.</summary>
    public double ResponseMean { get; }

    /// <summary>
    /// The log density of u = ln τ under τ's prior, its Jacobian included (<see cref="ResolvedTauPrior.LogDensity"/>); a
    /// sampled τ only. Where eᵘ leaves the doubles, τ is outside the support: −∞.
    /// </summary>
    public double LogTauDensity(double logTau) => Scale.LogDensity(logTau);

    /// <summary>Why <paramref name="fit"/> lies outside F, or null: the design's refusal, or an exact fit under the
    /// Jeffreys σ² prior.</summary>
    public BassRefusal? Refusal(BassStructureFit fit)
        => fit.Refusal ?? (fit.ResidualSquares <= _exactFitSquares ? BassRefusal.ExactFit : null);

    /// <summary>The structure's untempered log prior, D02's count weight plus its terms' (on F).</summary>
    public double LogStructurePrior(BassStructureFit fit)
        => Refusal(fit) is null ? Prior.LogCollectionWeight(fit.Terms) + fit.LogTermPrior : double.NegativeInfinity;

    /// <summary>The coefficient block's marginal for <paramref name="fit"/> at (<paramref name="tau"/>,
    /// <paramref name="inverseTemperature"/>), the integrated form's; σ²'s posterior reads it too.</summary>
    public GPriorMarginal Marginal(BassStructureFit fit, double tau, double inverseTemperature = 1.0)
    {
        if (Refusal(fit) is { } reason) throw new ArgumentException($"The structure lies outside F ({reason}).", nameof(fit));
        if (CoefficientForm is not RidgeForm)
            return CollapsedGPrior.CentredMarginalAt(fit.Terms, fit.ResidualSquares, fit.FittedSquares, Observations, ResponseMean,
                                                     Intercept, tau, inverseTemperature);
        RidgeColumns ridge = fit.Ridge ?? throw new ArgumentException("A structure under the ridge carries its columns' ridge terms.", nameof(fit));
        if (ridge.Count != fit.Terms) throw new ArgumentException("The ridge terms are another structure's.", nameof(fit));
        return CollapsedGPrior.CentredRidgeMarginalAt(ridge, fit.ResidualSquares, Observations, ResponseMean, Intercept, tau, inverseTemperature);
    }

    /// <summary>
    /// The coefficient block of <paramref name="design"/>'s current structure, (τ, b) to its conditional: the g-form's
    /// <see cref="CollapsedGPrior"/> on the design's factor and fit, or the ridge's <see cref="CollapsedRidgePrior"/> on its
    /// spectrum. The retained draws' composition, the collapsed pointwise content and the conditional means read it.
    /// </summary>
    public Func<double, double, GPriorConditional> CoefficientBlock(BassDesign design)
    {
        ArgumentNullException.ThrowIfNull(design);
        if (CoefficientForm is RidgeForm)
        {
            CollapsedRidgePrior ridge = CollapsedRidgePrior.FromCentred(design.Spectrum, design.ResidualSquares, design.ResponseMean, Intercept);
            return ridge.At;
        }
        CollapsedGPrior block = CollapsedGPrior.FromCentredLeastSquares(design.Factor, design.Coefficients.ToArray(), design.ResidualSquares,
                                                                         design.FittedSquares, design.Observations, design.ResponseMean, Intercept);
        return block.At;
    }

    /// <summary>
    /// Whether the integrated form is defined at <paramref name="inverseTemperature"/>: everywhere on [0, 1] under the
    /// proper defaults; b &gt; 0 under Jeffreys; and under the flat intercept, whose α integral spends one observation,
    /// g₁ + (Nb − 1)/2 &gt; 0. The linear form is defined on all of [0, 1].
    /// </summary>
    public bool Defines(double inverseTemperature) => Defines(Form, Intercept, Sigma2, Observations, inverseTemperature);

    /// <summary>
    /// <see cref="Defines(double)"/> for a target in <paramref name="form"/> under <paramref name="intercept"/> and σ²'s
    /// resolved prior <paramref name="sigma2"/> (null: Jeffreys) on <paramref name="observations"/> observations, before
    /// the target is built: what <see cref="BassResolved"/> checks a ladder's rungs by.
    /// </summary>
    public static bool Defines(TemperedForm form, GPriorCentre intercept, InverseGamma? sigma2, int observations, double inverseTemperature)
    {
        double b = inverseTemperature;
        if (!(b >= 0.0 && b <= 1.0)) return false;
        if (form == TemperedForm.Linear) return true;
        double shape = sigma2?.Shape ?? 0.0;
        if (intercept == GPriorCentre.Flat) return b > 0.0 && shape + 0.5 * (observations * b - 1) > 0.0;
        return b > 0.0 || sigma2 is not null;
    }

    /// <summary>A state with τ sampled: <paramref name="fit"/> with τ = e^<paramref name="logTau"/>.</summary>
    public BassEvaluation Evaluate(BassStructureFit fit, double logTau)
    {
        if (!SamplesTau) throw new InvalidOperationException("τ is held; evaluate the structure alone.");
        double structure = LogStructurePrior(fit);
        double prior = double.IsNegativeInfinity(structure) ? structure : structure + LogTauDensity(logTau);
        return new BassEvaluation(this, fit, Math.Exp(logTau), prior, Refusal(fit));
    }

    /// <summary>A state with τ held: <paramref name="fit"/> at the held τ.</summary>
    public BassEvaluation Evaluate(BassStructureFit fit)
    {
        double tau = HeldTau ?? throw new InvalidOperationException("τ is sampled; evaluate the structure with it.");
        return new BassEvaluation(this, fit, tau, LogStructurePrior(fit), Refusal(fit));
    }

    /// <summary>A chain's state: its structure's fit with its τ when sampled, or with the held one.</summary>
    public BassEvaluation Evaluate(BassState state)
    {
        ArgumentNullException.ThrowIfNull(state);
        return SamplesTau ? Evaluate(state.Fit, state.LogTau) : Evaluate(state.Fit);
    }

    ITargetEvaluation ISamplingTarget<BassState>.Evaluate(BassState state) => Evaluate(state);
}

/// <summary>
/// One state's evaluation under a <see cref="BassTarget"/>: its untempered log prior, and its log likelihood at any
/// inverse temperature in the target's form, each an O(1) query on the structure's statistics.
/// </summary>
public sealed class BassEvaluation : ITargetEvaluation
{
    private readonly BassTarget _target;
    private double? _untempered;

    internal BassEvaluation(BassTarget target, BassStructureFit fit, double tau, double logPrior, BassRefusal? refusal)
    {
        _target = target;
        Fit = fit;
        Tau = tau;
        LogPrior = logPrior;
        Refusal = refusal;
    }

    /// <summary>The structure's statistics.</summary>
    public BassStructureFit Fit { get; }

    /// <summary>τ, the state's or the held one.</summary>
    public double Tau { get; }

    /// <summary>Why the state lies outside F, or null.</summary>
    public BassRefusal? Refusal { get; }

    /// <inheritdoc/>
    public double LogPrior { get; }

    /// <summary>Whether the state lies off the support: refused, or its τ outside the prior's.</summary>
    public bool OffSupport => double.IsNegativeInfinity(LogPrior);

    /// <summary>ln p₁(y | structure, τ), the untempered collapsed marginal; −∞ off the support.</summary>
    public double LogMarginal => _untempered ??= OffSupport
        ? double.NegativeInfinity
        : _target.Marginal(Fit, Tau).LogMarginal(_target.Sigma2);

    /// <summary>The g-prior block's marginal at <paramref name="inverseTemperature"/>, for σ²'s posterior.</summary>
    public GPriorMarginal Marginal(double inverseTemperature = 1.0) => _target.Marginal(Fit, Tau, inverseTemperature);

    /// <inheritdoc/>
    public double LogLikelihood(double inverseTemperature)
    {
        double b = Check(inverseTemperature);
        if (OffSupport) return double.NegativeInfinity;
        if (_target.Form == TemperedForm.Linear) return b * LogMarginal;
        return b == 1.0 ? LogMarginal : _target.Marginal(Fit, Tau, b).LogMarginal(_target.Sigma2);
    }

    /// <inheritdoc/>
    public double LogLikelihoodSlope(double inverseTemperature)
    {
        double b = Check(inverseTemperature);
        if (OffSupport) throw new InvalidOperationException($"The state lies off the support ({Refusal?.ToString() ?? "τ"}); it has no likelihood.");
        return _target.Form == TemperedForm.Linear
            ? LogMarginal
            : _target.Marginal(Fit, Tau, b).LogMarginalSlope(_target.Sigma2);
    }

    private double Check(double b)
    {
        if (!(b >= 0.0 && b <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(b), "The inverse temperature must lie in [0, 1].");
        if (!_target.Defines(b))
            throw new ArgumentOutOfRangeException(nameof(b),
                $"The integrated form is undefined at b = {b} under {(_target.Intercept == GPriorCentre.Flat ? "the flat intercept" : "the Jeffreys σ² prior")}; the linear form reaches it.");
        return b;
    }
}
