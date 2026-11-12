using System;
using System.Collections.Generic;
using System.Linq;
using Maths.Distributions;
using Maths.Regression.GPrior;
using Maths.Regression.Spline.Bars;
using Maths.Samplers;
using Maths.Samplers.Rjmcmc;

namespace Maths.Regression.Spline.Bass;

/// <summary>
/// A structure's account for the collapsed target: M, the sum of its terms' ln π (D01), and its least-squares fit in
/// the centred basis, ‖r‖² and Q_c, with, under a prior the g-form's closed form does not serve, its columns' terms: the
/// normalized ridge's (a spectrum, or a bordered candidate of one) or a <see cref="PriorPencil"/>; or the design's refusal,
/// which puts it outside F. It also names the working response it was fitted on, whose ȳ_w and weighting the target reads with these
/// sums, so a state carries them wherever it goes: through a swap, a particle's fork, or a flat-histogram chain. Null reads the
/// target's own.
/// </summary>
public readonly record struct BassStructureFit(int Terms, double LogTermPrior, double ResidualSquares, double FittedSquares,
                                               BassRefusal? Refusal = null, RidgeColumns? Ridge = null, WorkingResponse? Working = null,
                                               OutputForms? Outputs = null)
{
    /// <summary>The design's current structure, whose terms' log prior is <paramref name="logTermPrior"/>. Several outputs (E34) carry
    /// their K × K forms in <see cref="Outputs"/>, the sums of squares then NaN.</summary>
    public static BassStructureFit Of(BassDesign design, double logTermPrior)
    {
        ArgumentNullException.ThrowIfNull(design);
        return new(design.Count, logTermPrior, design.ResidualSquares, design.FittedSquares, null, design.PriorTerms, design.Working,
                   design.OutputForms);
    }

    /// <summary>A proposed structure, admitted or refused, whose terms' log prior is <paramref name="logTermPrior"/>, fitted on
    /// <paramref name="working"/>.</summary>
    public static BassStructureFit Of(DesignProposal proposal, double logTermPrior, WorkingResponse? working = null) => proposal.Refusal is { } reason
        ? Refused(reason)
        : new(proposal.Columns, logTermPrior, proposal.ResidualSquares, proposal.FittedSquares, null, proposal.Ridge, working, proposal.Outputs);

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
/// design's centred statistics (<see cref="CollapsedGPrior.CentredMarginalAt"/>), the normalized ridge's on the
/// structure's ridge terms (<see cref="CollapsedGPrior.CentredRidgeMarginalAt"/>), or under weights with the prior free of
/// them a <see cref="PriorPencil"/>'s (<see cref="CollapsedGPrior.CentredPencilMarginalAt"/>), under the intercept's convention
/// (level centring, Liang's flat intercept or R's zero centre) and σ²'s prior (scaled to the response, Jeffreys, or a declared
/// inverse gamma). On the likelihood path at b it takes D05's integrated form p_b, the coefficients and σ² integrated
/// against the tempered likelihood, or the linear form b·ln p₁. With σ² held (<see cref="GivenSigma2"/>) it is the marginal
/// given σ² from the same S and determinant term, the coefficients alone integrated. Under a frozen read
/// (<see cref="BassPriorReference"/>) it holds α's prior centre and the coefficient prior's constant apart from the working response.</item>
/// <item><b>σ²'s rate,</b> under a hierarchical prior (<see cref="Sigma2Rate"/>): a coordinate of the state as τ is, its density on
/// v = ln r in the prior, the marginal given r. With σ² given, the prior also carries IG(σ²; a, r). One target serves every chain,
/// rung and particle, each state holding its own r.</item>
/// <item><b>A working response that moves</b> (<see cref="BassStructureMoves.Follow"/>) is read off the state too: its fit names the
/// working response it was fitted on, for ȳ_w, the weighting and the exact-fit floor. So the target stays one target, and what the
/// priors read stays frozen in it.</item>
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
    private const double FloorSquared = ExactFitMultiple * MachineEpsilon * (ExactFitMultiple * MachineEpsilon);

    private readonly double _exactFitSquares;   // ‖r‖² at or below which Jeffreys refuses; −∞ under a proper prior
    private readonly WorkingResponse? _resolvedOn;   // the working response `resolved` read, whose floor is the one above

    /// <summary>The target <paramref name="config"/> sets on its data (<paramref name="resolved"/>, read from
    /// <paramref name="working"/>), over <paramref name="prior"/>'s structures, tempered in <paramref name="form"/>.</summary>
    public BassTarget(BassConfig config, BassResolved resolved, BassStructurePrior prior, TemperedForm form = TemperedForm.Integrated, WorkingResponse? working = null)
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
        GivenSigma2 = resolved.GivenSigma2;
        Scale = resolved.Scale;
        Observations = resolved.Observations;
        ResponseMean = resolved.ResponseMean;
        Weighting = resolved.Weighting;
        PriorLevel = resolved.PriorLevel;
        CoefficientConstant = resolved.CoefficientConstant;
        Sigma2Rate = resolved.Sigma2Rate;
        Outputs = resolved.Outputs;
        OutputSkewSd = resolved.Outputs is not null && config.Noise?.SkewPrior is { IsHeld: false } skew ? skew.Sd : null;
        _resolvedOn = working;
        _exactFitSquares = Sigma2 is null && GivenSigma2 is null && Outputs is null ? FloorSquared * Weighting.Mass * resolved.WeightedMeanSquare : double.NegativeInfinity;
    }

    // `source` given σ² = `sigma2`: every other part shared.
    private BassTarget(BassTarget source, double sigma2)
    {
        Prior = source.Prior;
        Form = source.Form;
        CoefficientForm = source.CoefficientForm;
        Intercept = source.Intercept;
        Sigma2 = source.Sigma2;
        GivenSigma2 = sigma2;
        Scale = source.Scale;
        Observations = source.Observations;
        ResponseMean = source.ResponseMean;
        Weighting = source.Weighting;
        PriorLevel = source.PriorLevel;
        CoefficientConstant = source.CoefficientConstant;
        Sigma2Rate = source.Sigma2Rate;
        Outputs = source.Outputs;
        OutputSkewSd = source.OutputSkewSd;
        _resolvedOn = source._resolvedOn;
        _exactFitSquares = double.NegativeInfinity;
    }

    /// <summary>This target with its structure step given σ² = <paramref name="sigma2"/>, the form for a σ² drawn apart from the
    /// structure (D04's generalized-hyperbolic row) or supplied from outside: a block that moves σ² hands the chain this
    /// target (<see cref="ReversibleJumpChain{TState}.ReplaceTarget"/>). <see cref="Sigma2"/> stays the prior such a block
    /// draws against.</summary>
    public BassTarget Given(double sigma2)
    {
        if (!(sigma2 > 0.0) || double.IsPositiveInfinity(sigma2)) throw new ArgumentOutOfRangeException(nameof(sigma2), "σ² must be positive and finite.");
        if (Outputs is not null) throw new InvalidOperationException("Several outputs integrate Σ; a Σ given arrives with E34's noise (slice 4).");
        return new BassTarget(this, sigma2);
    }

    /// <summary>Whether the values are complete: under proper σ² (or Σ) and intercept priors, or with σ² given and the intercept
    /// proper. τ's prior is a separate condition, for an evidence, since a ratio across b does not read it.</summary>
    public bool HasCompleteValues => (Sigma2 is not null || GivenSigma2 is not null || Outputs is not null) && Intercept != GPriorCentre.Flat;

    /// <summary>Several outputs' resolution (E34): their means, levels and Σ's prior, which the marginal reads; null for one response.</summary>
    public BassResolvedOutputs? Outputs { get; }

    /// <summary>
    /// Several outputs' marginal for <paramref name="fit"/> at (<paramref name="tau"/>, <paramref name="inverseTemperature"/>): the g-form's
    /// K × K closed form on the structure's forms (<see cref="OutputMarginal.CentredAt"/>), under the intercept's convention; Σ's posterior
    /// reads it too (E34).
    /// </summary>
    public OutputMarginal OutputsMarginal(BassStructureFit fit, double tau, double inverseTemperature = 1.0)
    {
        if (Outputs is not { } outputs) throw new InvalidOperationException("The target holds one response; its marginal is Marginal's.");
        if (Refusal(fit) is { } reason) throw new ArgumentException($"The structure lies outside F ({reason}).", nameof(fit));
        OutputForms forms = fit.Outputs ?? throw new ArgumentException("A structure of several outputs carries their forms.", nameof(fit));
        // The reading the forms are on (E34's slice 4): its weights' mass, Σ ln wᵢ and count, a random-effect block's −ln|V| among them, and
        // its outputs' means; at the response itself, the resolution's means and unit weights.
        WorkingResponse? reading = fit.Working;
        if (forms.Prior is { } pencil)
            return pencil.At(forms.Residual, Observations, tau / CoefficientConstant, inverseTemperature, reading is null || reading.IsResponseItself ? null : reading.Weighting);
        if (reading is null || reading.IsResponseItself)
            return OutputMarginal.CentredAt(fit.Terms, forms.Residual, forms.Fitted, Observations, outputs.Means, Intercept, tau / CoefficientConstant, inverseTemperature);
        // Under a noise the level is frozen apart from v (the resolution's raw mean), so level centring reads the gap ȳ_w − ℓ.
        return OutputMarginal.CentredAt(fit.Terms, forms.Residual, forms.Fitted, Observations, reading.OutputMeans.ToArray(), Intercept, tau / CoefficientConstant,
                                        inverseTemperature, Intercept == GPriorCentre.Mean && reading.Latents is not null ? outputs.Levels : null, reading.Weighting);
    }

    /// <summary>The prior on structures.</summary>
    public BassStructurePrior Prior { get; }

    /// <summary>Several outputs' drawn skew's s, γ | Σ ~ N(0, s²Σ) (E34's decision 6); null where the skew is not drawn or the response is one.</summary>
    internal double? OutputSkewSd { get; init; }

    /// <summary>
    /// Σ's prior with pseudo-rows, latents whose prior is N(0, c·Σ) per row (E34's slice 4: the skew γ, one row at c = s², and a grouping's
    /// effects U_k, q_k rows at c = s_k²): IW(ν₀ + Σn, Ψ₀ + Σ X/c), X each one's K × K cross product, and their joint density with Σ
    /// integrated, the matrix t
    /// <para>ln p = −(K/2)·Σ n·ln(πc) + ln Γ_K((ν₀ + Σn)/2) − ln Γ_K(ν₀/2) + (ν₀/2)ln|Ψ₀| − ((ν₀ + Σn)/2)ln|Ψ₀ + Σ X/c|.</para>
    /// </summary>
    internal static (InverseWishart Prior, double LogDensity) WithPseudoRows(InverseWishart prior, IEnumerable<(double[,] Cross, int Count, double Scale2)> rows)
    {
        double[,] scale = prior.Scale;
        int k = scale.GetLength(0), count = 0;
        double logScales = 0.0;
        foreach ((double[,] cross, int n, double c) in rows)
        {
            for (int a = 0; a < k; a++)
                for (int d = 0; d < k; d++) scale[a, d] += cross[a, d] / c;
            count += n;
            logScales += n * Math.Log(Math.PI * c);
        }
        double nu = prior.DegreesOfFreedom;
        var updated = InverseWishart.Of(nu + count, scale);
        double logDensity = -0.5 * k * logScales + SpecialFunctions.LogMultivariateGamma(k, 0.5 * (nu + count))
                            - SpecialFunctions.LogMultivariateGamma(k, 0.5 * nu) + 0.5 * nu * prior.LogDeterminantScale - 0.5 * (nu + count) * updated.LogDeterminantScale;
        return (updated, logDensity);
    }

    // Several outputs' pseudo-rows at `working`: γ's where the skew is drawn, each grouping's U where random intercepts are a latent.
    internal IEnumerable<(double[,] Cross, int Count, double Scale2)> PseudoRows(WorkingResponse? working)
    {
        if (OutputSkewSd is double sd && working?.Latents?.OutputSkew is { } gamma)
        {
            int k = gamma.Length;
            var cross = new double[k, k];
            for (int a = 0; a < k; a++)
                for (int d = 0; d < k; d++) cross[a, d] = gamma[a] * gamma[d];
            yield return (cross, 1, sd * sd);
        }
        if (working?.Effects is { OutputEffects: not null } effects)
            foreach ((double[,] Cross, int Count, double Scale2) row in effects.OutputPseudoRows()) yield return row;
    }

    /// <summary>The tempered form on the likelihood path.</summary>
    public TemperedForm Form { get; }

    /// <summary>The coefficient prior's covariance form: the g-form or the normalized ridge.</summary>
    public CoefficientForm CoefficientForm { get; }

    /// <summary>The intercept's convention.</summary>
    public GPriorCentre Intercept { get; }

    /// <summary>σ²'s prior as resolved; null is the Jeffreys 1/σ², or a σ² held by the config. Under a hierarchical prior, the
    /// inverse gamma at its rate's mean; a state's own is <see cref="Sigma2At"/>.</summary>
    public InverseGamma? Sigma2 { get; }

    /// <summary>σ²'s rate under a hierarchical prior; null under every other, whose state has no rate.</summary>
    public ResolvedSigma2Rate? Sigma2Rate { get; }

    /// <summary>σ²'s prior at <paramref name="state"/>: its rate's inverse gamma under a hierarchical prior, <see cref="Sigma2"/>
    /// otherwise.</summary>
    public InverseGamma? Sigma2At(BassState state) => Sigma2Rate is { } rate ? rate.PriorAt(state.LogRate) : Sigma2;

    /// <summary>α's prior centre under a frozen read; null is the convention's own, ȳ under level centring.</summary>
    public double? PriorLevel { get; }

    /// <summary>κ in the coefficient prior's covariance σ²κ(τG)⁻¹ (<see cref="BassPriorReference.CoefficientConstant"/>); 1 is the
    /// Gaussian's.</summary>
    public double CoefficientConstant { get; }

    /// <summary>The σ² the structure step's marginal is given; null when σ² is integrated out. A state may hold its own
    /// (<see cref="GivenAt"/>).</summary>
    public double? GivenSigma2 { get; }

    /// <summary>The σ² <paramref name="fit"/>'s marginal is given: the state's own when its working response's latents hold one, drawn
    /// apart from the structure (E1's decision 1, w on the record), else <see cref="GivenSigma2"/>; null when σ² is integrated out.</summary>
    public double? GivenAt(BassStructureFit fit) => fit.Working?.Latents?.Sigma2 ?? GivenSigma2;

    // ln p(σ²) of a σ² the state holds, the state's prior's part beyond the structure and τ: the inverse gamma's, Jeffreys' −ln σ², and
    // nothing for a held σ², a point mass. A hierarchical prior's reads the rate (BassTarget.Evaluate(BassState)).
    private double LogHeldPrior(BassStructureFit fit)
    {
        if (fit.Working?.Latents?.Sigma2 is not double sigma2 || GivenSigma2 is not null) return 0.0;
        return Sigma2 is { } prior ? prior.LogPdf(sigma2) : -Math.Log(sigma2);
    }

    /// <summary>The count, the mass and Σ ln wᵢ the marginal reads beyond the design's weighted statistics, on the model's working
    /// response; a fit reads the one it was fitted on. The count does not move.</summary>
    public GPriorWeights Weighting { get; }

    /// <summary>τ's prior on the coefficient form's τ₀, sampled or held.</summary>
    public ResolvedTauPrior Scale { get; }

    /// <summary>The held τ; null when τ is sampled.</summary>
    public double? HeldTau => Scale.HeldTau;

    /// <summary>Whether τ is a coordinate of the state: sampled, not held.</summary>
    public bool SamplesTau => !Scale.IsHeld;

    /// <summary>N, the observations.</summary>
    public int Observations { get; }

    /// <summary>ȳ under W, which the zero centre's level energy reads, on the model's working response; a fit reads the one it was
    /// fitted on.</summary>
    public double ResponseMean { get; }

    /// <summary>
    /// The log density of u = ln τ under τ's prior, its Jacobian included (<see cref="ResolvedTauPrior.LogDensity"/>); a
    /// sampled τ only. Where eᵘ leaves the doubles, τ is outside the support: −∞.
    /// </summary>
    public double LogTauDensity(double logTau) => Scale.LogDensity(logTau);

    // ln s is a coordinate of the joint target, like ln τ. Replicas need the same prior as the scale step for powered swaps.
    internal double WithGroupScalePrior(double logPrior, BassStructureFit fit)
        => (fit.Working ?? _resolvedOn)?.Whitening is { } block ? logPrior + block.LogScalePrior : logPrior;

    /// <summary>Why <paramref name="fit"/> lies outside F, or null: the design's refusal, or an exact fit under the
    /// Jeffreys σ² prior.</summary>
    public BassRefusal? Refusal(BassStructureFit fit)
        => fit.Refusal ?? (fit.ResidualSquares <= ExactFitSquares(fit.Working) ? BassRefusal.ExactFit : null);

    // The floor on `working`: the round-off of its own squares about the level, s² standing in under a frozen read that states one.
    private double ExactFitSquares(WorkingResponse? working)
    {
        if (working is null || ReferenceEquals(working, _resolvedOn) || double.IsNegativeInfinity(_exactFitSquares)) return _exactFitSquares;
        double mass = working.Weighting.Mass;
        if (working.Reference?.MeanSquare is double frozen) return FloorSquared * mass * frozen;
        double gap = Intercept == GPriorCentre.Zero ? working.ResponseMean : PriorLevel is double level ? working.ResponseMean - level : 0.0;
        return FloorSquared * (working.ResponseSquares + mass * gap * gap);
    }

    /// <summary>The structure's untempered log prior, D02's count weight plus its terms' (on F).</summary>
    public double LogStructurePrior(BassStructureFit fit)
        => Refusal(fit) is null ? Prior.LogCollectionWeight(fit.Terms) + fit.LogTermPrior : double.NegativeInfinity;

    /// <summary>The coefficient block's marginal for <paramref name="fit"/> at (<paramref name="tau"/>,
    /// <paramref name="inverseTemperature"/>), the integrated form's; σ²'s posterior reads it too.</summary>
    public GPriorMarginal Marginal(BassStructureFit fit, double tau, double inverseTemperature = 1.0)
    {
        if (Outputs is not null) throw new InvalidOperationException("Several outputs' marginal is a K × K block's: read OutputsMarginal.");
        if (Refusal(fit) is { } reason) throw new ArgumentException($"The structure lies outside F ({reason}).", nameof(fit));
        double precision = tau / CoefficientConstant;
        GPriorWeights weighting = fit.Working?.Weighting ?? Weighting;
        if (fit.Ridge is PriorPencil pencil)
        {
            if (pencil.Count != fit.Terms) throw new ArgumentException("The pencil is another structure's.", nameof(fit));
            return CollapsedGPrior.CentredPencilMarginalAt(pencil, fit.ResidualSquares, Observations, weighting, Intercept, precision, inverseTemperature);
        }
        double mean = fit.Working?.ResponseMean ?? ResponseMean;
        if (CoefficientForm is not RidgeForm)
            return CollapsedGPrior.CentredMarginalAt(fit.Terms, fit.ResidualSquares, fit.FittedSquares, Observations, mean,
                                                     Intercept, precision, inverseTemperature, weighting, PriorLevel);
        RidgeColumns ridge = fit.Ridge ?? throw new ArgumentException("A structure under the ridge carries its columns' ridge terms.", nameof(fit));
        if (ridge.Count != fit.Terms) throw new ArgumentException("The ridge terms are another structure's.", nameof(fit));
        return CollapsedGPrior.CentredRidgeMarginalAt(ridge, fit.ResidualSquares, Observations, mean, Intercept, precision, inverseTemperature, weighting, PriorLevel);
    }

    // The marginal's terms that τ and the structure do not move, held once per σ² prior and count (lever 5 of the efficiency pass).
    private readonly Sigma2PriorTerms _sigma2Terms = new();

    // ln p_b of a structure's marginal: σ² integrated against its prior at this state, or given.
    internal double LogValue(in GPriorMarginal marginal, InverseGamma? sigma2, double? given) => given is double g ? marginal.LogMarginalGiven(g) : marginal.LogMarginal(sigma2, _sigma2Terms);

    // ∂/∂b of LogValue.
    internal double LogSlope(in GPriorMarginal marginal, InverseGamma? sigma2, double? given) => given is double g ? marginal.LogMarginalGivenSlope(g) : marginal.LogMarginalSlope(sigma2);

    /// <summary>
    /// The coefficient block of <paramref name="design"/>'s current structure, (τ, b) to its conditional: the g-form's
    /// <see cref="CollapsedGPrior"/> on the design's factor and fit, the ridge's <see cref="CollapsedRidgePrior"/> on its
    /// spectrum, or <see cref="CollapsedPencilPrior"/> on its pencil. The retained draws' composition, the collapsed pointwise
    /// content and the conditional means read it.
    /// </summary>
    public Func<double, double, GPriorConditional> CoefficientBlock(BassDesign design)
    {
        ArgumentNullException.ThrowIfNull(design);
        if (design.Pencil is { } pencil)
        {
            CollapsedPencilPrior apart = CollapsedPencilPrior.FromCentred(pencil, design.ResidualSquares, design.Observations, design.Weighting, design.ResponseMean,
                                                                          design.PriorLevel, Intercept);
            return Precision(apart.At);
        }
        if (CoefficientForm is RidgeForm)
        {
            CollapsedRidgePrior ridge = CollapsedRidgePrior.FromCentred(design.Spectrum, design.ResidualSquares, design.ResponseMean, Intercept, design.Weighting, PriorLevel);
            return Precision(ridge.At);
        }
        CollapsedGPrior block = CollapsedGPrior.FromCentredLeastSquares(design.Factor, design.Coefficients.ToArray(), design.ResidualSquares,
                                                                         design.FittedSquares, design.Observations, design.ResponseMean, Intercept, design.Weighting, PriorLevel);
        return Precision(block.At);
    }

    // A block's (τ, b) with τ read through the coefficient prior's constant.
    private Func<double, double, GPriorConditional> Precision(Func<double, double, GPriorConditional> at)
        => CoefficientConstant == 1.0 ? at : (tau, b) => at(tau / CoefficientConstant, b);

    /// <summary>
    /// Whether the integrated form is defined at <paramref name="inverseTemperature"/>: everywhere on [0, 1] under the
    /// proper defaults; b &gt; 0 under Jeffreys; and under the flat intercept, whose α integral spends one observation,
    /// g₁ + (Nb − 1)/2 &gt; 0, or b &gt; 0 alone with σ² given. The linear form is defined on all of [0, 1].
    /// </summary>
    public bool Defines(double inverseTemperature)
        => Defines(Form, Intercept, Sigma2, Weighting.Count, inverseTemperature, GivenSigma2 is not null, Outputs?.Prior);

    /// <summary>
    /// <see cref="Defines(double)"/> for a target in <paramref name="form"/> under <paramref name="intercept"/> and σ²'s
    /// resolved prior <paramref name="sigma2"/> (null: Jeffreys, unless <paramref name="givenSigma2"/>) on
    /// <paramref name="observations"/> observations, before the target is built: what <see cref="BassResolved"/> checks a
    /// ladder's rungs by.
    /// </summary>
    public static bool Defines(TemperedForm form, GPriorCentre intercept, InverseGamma? sigma2, double observations, double inverseTemperature,
                               bool givenSigma2 = false, InverseWishart? outputsPrior = null)
    {
        double b = inverseTemperature;
        if (!(b >= 0.0 && b <= 1.0)) return false;
        if (form == TemperedForm.Linear) return true;
        if (outputsPrior is { } iw)
            return (intercept != GPriorCentre.Flat || b > 0) &&
                   iw.DegreesOfFreedom + observations * b - (intercept == GPriorCentre.Flat ? 1 : 0) > iw.Dimension - 1;
        if (givenSigma2) return intercept != GPriorCentre.Flat || b > 0.0;
        double shape = sigma2?.Shape ?? 0.0;
        if (intercept == GPriorCentre.Flat) return b > 0.0 && shape + 0.5 * (observations * b - 1) > 0.0;
        return b > 0.0 || sigma2 is not null;
    }

    /// <summary>A state with τ sampled: <paramref name="fit"/> with τ = e^<paramref name="logTau"/>. Under a hierarchical σ² prior the
    /// rate is part of the state, so evaluate that.</summary>
    public BassEvaluation Evaluate(BassStructureFit fit, double logTau)
    {
        if (!SamplesTau) throw new InvalidOperationException("τ is held; evaluate the structure alone.");
        RefuseSampledRate();
        double structure = LogStructurePrior(fit);
        double prior = double.IsNegativeInfinity(structure) ? structure : structure + LogTauDensity(logTau) + LogHeldPrior(fit);
        return new BassEvaluation(this, fit, Math.Exp(logTau), prior, Refusal(fit));
    }

    /// <summary>A state with τ held: <paramref name="fit"/> at the held τ. Under a hierarchical σ² prior the rate is part of the state,
    /// so evaluate that.</summary>
    public BassEvaluation Evaluate(BassStructureFit fit)
    {
        double tau = HeldTau ?? throw new InvalidOperationException("τ is sampled; evaluate the structure with it.");
        RefuseSampledRate();
        double structure = LogStructurePrior(fit);
        return new BassEvaluation(this, fit, tau, double.IsNegativeInfinity(structure) ? structure : structure + LogHeldPrior(fit), Refusal(fit));
    }

    /// <summary>A chain's state: its structure's fit with its τ when sampled, or with the held one, and its rate under a hierarchical
    /// σ² prior.</summary>
    public BassEvaluation Evaluate(BassState state)
    {
        ArgumentNullException.ThrowIfNull(state);
        if (Outputs?.Rates is { } outputRates)
        {
            double logPrior = LogStructurePrior(state.Fit) + outputRates.LogDensity(state.LogOutputRates);
            if (SamplesTau) logPrior += LogTauDensity(state.LogTau);
            var priorAtState = double.IsNegativeInfinity(logPrior) ? null : outputRates.PriorAt(state.LogOutputRates);
            return new BassEvaluation(this, state.Fit, SamplesTau ? state.Tau : HeldTau!.Value, logPrior, Refusal(state.Fit),
                                      outputPrior: priorAtState);
        }
        if (Sigma2Rate is not { } rate) return SamplesTau ? Evaluate(state.Fit, state.LogTau) : Evaluate(state.Fit);
        BassStructureFit fit = state.Fit;
        double prior = LogStructurePrior(fit);
        if (!double.IsNegativeInfinity(prior) && SamplesTau) prior += LogTauDensity(state.LogTau);
        if (!double.IsNegativeInfinity(prior)) prior += rate.LogDensity(state.LogRate);
        if (!double.IsNegativeInfinity(prior) && GivenAt(fit) is double given) prior += rate.LogConditional(given, state.LogRate);
        return new BassEvaluation(this, fit, SamplesTau ? Math.Exp(state.LogTau) : HeldTau!.Value, prior, Refusal(fit), state.LogRate);
    }

    private void RefuseSampledRate()
    {
        if (Sigma2Rate is not null || Outputs?.Rates is not null)
            throw new InvalidOperationException("The covariance prior's rates are sampled; evaluate the state, which holds them.");
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
    private readonly double _logRate;   // σ²'s rate under a hierarchical prior
    private readonly InverseWishart? _outputPrior;
    private double? _untempered;

    internal BassEvaluation(BassTarget target, BassStructureFit fit, double tau, double logPrior, BassRefusal? refusal, double logRate = double.NaN,
                            InverseWishart? outputPrior = null)
    {
        _target = target;
        _logRate = logRate;
        _outputPrior = outputPrior ?? target.Outputs?.Prior;
        // Several outputs' latents with a Σ-scaled prior (E34's slice 4): the drawn skew γ | Σ ~ N(0, s²Σ) and random intercepts' U_k | Σ ~
        // N(0, s_k²·I ⊗ Σ) enter Σ's posterior as pseudo-rows, and their matrix t under Σ's prior joins the log prior, so the target is
        // p(structure, τ, γ, U, Y | v, s).
        if (target.Outputs is not null && fit.Working is { Latents: not null } && _outputPrior is { } basis && !double.IsNegativeInfinity(logPrior))
        {
            var rows = target.PseudoRows(fit.Working).ToList();
            if (rows.Count > 0)
            {
                (_outputPrior, double logRows) = BassTarget.WithPseudoRows(basis, rows);
                logPrior += logRows;
            }
        }
        Fit = fit;
        Tau = tau;
        LogPrior = target.WithGroupScalePrior(logPrior, fit);
        Refusal = refusal;
    }

    /// <summary>The target this evaluation read, which a replica's own target follows when a swap hands it this state.</summary>
    internal BassTarget Target => _target;

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
        : LogValueAt(1.0);

    /// <summary>σ²'s prior at this state (<see cref="BassTarget.Sigma2At"/>); the rate's mean's, off the support.</summary>
    public InverseGamma? Sigma2 => _target.Sigma2Rate is { } rate && !OffSupport ? rate.PriorAt(_logRate) : _target.Sigma2;

    /// <summary>The g-prior block's marginal at <paramref name="inverseTemperature"/>, for σ²'s posterior.</summary>
    public GPriorMarginal Marginal(double inverseTemperature = 1.0) => _target.Marginal(Fit, Tau, inverseTemperature);

    internal InverseWishart? OutputPrior => _outputPrior;

    /// <summary>Several outputs' K × K block's marginal at <paramref name="inverseTemperature"/>, for Σ's posterior (E34).</summary>
    public OutputMarginal OutputsMarginal(double inverseTemperature = 1.0) => _target.OutputsMarginal(Fit, Tau, inverseTemperature);

    /// <inheritdoc/>
    public double LogLikelihood(double inverseTemperature)
    {
        double b = Check(inverseTemperature);
        if (OffSupport) return double.NegativeInfinity;
        if (_target.Form == TemperedForm.Linear) return b * LogMarginal;
        return b == 1.0 ? LogMarginal : LogValueAt(b);
    }

    /// <inheritdoc/>
    public double LogLikelihoodSlope(double inverseTemperature)
    {
        double b = Check(inverseTemperature);
        if (OffSupport) throw new InvalidOperationException($"The state lies off the support ({Refusal?.ToString() ?? "τ"}); it has no likelihood.");
        if (_target.Outputs is { } outputs)
            return _target.Form == TemperedForm.Linear ? LogMarginal : OutputsMarginal(b).LogLikelihoodSlope(_outputPrior);
        return _target.Form == TemperedForm.Linear
            ? LogMarginal
            : _target.LogSlope(_target.Marginal(Fit, Tau, b), Sigma2, _target.GivenAt(Fit));
    }

    // ln p_b at this state, with what a hierarchical prior's reading of r leaves out below the doubles. At b = 0 the marginal reads no r.
    private double LogValueAt(double b)
    {
        if (_target.Outputs is not null) return _target.OutputsMarginal(Fit, Tau, b).LogMarginal(_outputPrior);
        double? given = _target.GivenAt(Fit);
        double value = _target.LogValue(_target.Marginal(Fit, Tau, b), Sigma2, given);
        return _target.Sigma2Rate is { } rate && b > 0.0 && given is null ? value + rate.LogRemainder(_logRate) : value;
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
