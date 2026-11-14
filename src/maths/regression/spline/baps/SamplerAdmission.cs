using System;
using System.Collections.Generic;
using Maths.LinAlg;

namespace Maths.Regression.Spline.Baps;

/// <summary>
/// How often a BAPS sampler estimates its system's condition (shared-code defects brief, slices 6 and 7). Every sampler step factors
/// a system A(c) = Σₖ cₖSₖ: fixed positive semidefinite terms Sₖ and positive coefficients cₖ, drawn or fixed at 1, such as
/// ZᵀZ + λP. <see cref="BapsAdmission"/>'s rule (positive definite, κ(E) ≤ 10¹⁴) is the check, and this setting says
/// where its estimate runs. Every mode checks positive definiteness at every step, and a run every mode admits draws the
/// same chain in each.
/// </summary>
public enum ConditionCheck
{
    /// <summary>The estimate at every factorization, <see cref="BapsAdmission"/>'s rule at each step: the reference the other two
    /// modes are checked against.</summary>
    EveryFactorization,

    /// <summary>
    /// One exact reading per run at the trace-balance reference point c_ref, then a scalar bound at each step, and the
    /// estimate at any step the bound does not clear. With ρₖ = cₖ / c_ref,ₖ every Rayleigh quotient of A(c) lies within
    /// [min ρ, max ρ] times A(c_ref)'s, which carries to the gate's reading as
    /// κ₁(E(c)) ≤ m·((1 + √ν)/2)·(max ρ / min ρ)·‖E(c_ref)⁻¹‖₁, with m the most nonzeros in a column of A, proved in Lean as
    /// <c>cond1_unitScale_le_sharp</c>. A step whose bound, doubled for rounding, is at most 10¹⁴ is admitted unestimated,
    /// and the estimate would have admitted it. So it decides every step as <see cref="EveryFactorization"/> does, and differs
    /// in cost alone. The default.
    /// </summary>
    Screened,

    /// <summary>
    /// One exact reading per run at the reference point, which refuses a design that does not identify the penalty's null
    /// space before the first draw, and positive definiteness alone at each step. The opt-out for a model already known to
    /// be sound: a step whose drawn coefficients lie far from the reference is admitted while it stays positive definite.
    /// </summary>
    ReferenceOnly,
}

/// <summary>
/// What a run's condition checks saw at one of its sampler's systems, over every chain and every step, burn-in included.
/// The fields reduce what the steps saw and leave the draws untouched. They ride in each chain's state, so a resumed run
/// reports what the uninterrupted one does.
/// </summary>
/// <param name="System">The system, as its refusal names it.</param>
/// <param name="Check">The run's <see cref="ConditionCheck"/>.</param>
/// <param name="ReferenceCondition">κ₁(E) at the reference point, read exactly; +∞ when the reference is not positive
/// definite, and NaN under <see cref="ConditionCheck.EveryFactorization"/>, which takes no reading there.</param>
/// <param name="WidestRatio">The largest max ρ / min ρ any step reached: how far the coefficients went from the
/// reference point, in every mode.</param>
/// <param name="Factorizations">The factorizations checked.</param>
/// <param name="Estimates">The condition estimates that ran: one per factorization under
/// <see cref="ConditionCheck.EveryFactorization"/>, one per step the screen did not clear under
/// <see cref="ConditionCheck.Screened"/>, and none under <see cref="ConditionCheck.ReferenceOnly"/>.</param>
/// <param name="LargestCondition">The largest estimate seen; NaN when none ran.</param>
public sealed record ConditionRecord(string System, ConditionCheck Check, double ReferenceCondition, double WidestRatio,
                                     long Factorizations, long Estimates, double LargestCondition);

/// <summary>One chain's condition records at one system, carried in its state.</summary>
internal sealed class ConditionTally
{
    public double WidestRatio = 1.0;
    public long Factorizations, Estimates;
    public double LargestCondition;   // read only when Estimates > 0

    public ConditionTally Copy() => (ConditionTally)MemberwiseClone();

    public static ConditionTally[] Copy(ConditionTally[] tallies)
    {
        var copy = new ConditionTally[tallies.Length];
        for (int i = 0; i < tallies.Length; i++) copy[i] = tallies[i].Copy();
        return copy;
    }

    public static ConditionTally[] Fresh(int systems)
    {
        var tallies = new ConditionTally[systems];
        for (int i = 0; i < systems; i++) tallies[i] = new ConditionTally();
        return tallies;
    }
}

/// <summary>
/// One sampler system's condition check for one run (<see cref="ConditionCheck"/>): the reference point's coefficients,
/// its exact reading where the mode takes one, and the screen. Immutable once built, so the chains share it; each chain's
/// records go to its own <see cref="ConditionTally"/>.
/// </summary>
/// <remarks>
/// The screen's derivation, for A(c) = Σₖ cₖSₖ, ρₖ = cₖ / c_ref,ₖ, s = max ρ / min ρ, and E(c) A(c)'s Jacobi scaling:
/// <list type="bullet">
/// <item>xᵀA(c)x lies within [min ρ, max ρ]·xᵀA(c_ref)x, and diag A(c) within [min ρ, max ρ]·diag A(c_ref). So
/// λ_min(E(c)) ≥ λ_min(E(c_ref)) / s: the two Jacobi scalings differ by factors the ratio already covers.</item>
/// <item>E(c) has a unit diagonal and entries of at most 1, so ‖E(c)‖₁ ≤ m.</item>
/// <item>‖E(c)⁻¹‖₁ ≤ ((1 + √ν)/2)·‖E(c)⁻¹‖₂ = ((1 + √ν)/2) / λ_min(E(c)) for a positive definite matrix, a column step that
/// is attained, and ‖E(c_ref)⁻¹‖₂ ≤ ‖E(c_ref)⁻¹‖₁ for a symmetric one.</item>
/// <item>The estimate errs low, so it never reads above κ₁(E(c)) ≤ m·((1 + √ν)/2)·s·‖E(c_ref)⁻¹‖₁, up to rounding, which the
/// factor 2 covers. Its own shortfall has no bound, so the reference is read exactly, from ν solves.</item>
/// </list>
/// The chain is proved in Lean (<c>lean/Enthymemata/PenalizedAdmission.lean</c>, <c>cond1_unitScale_le_sharp</c>) in exact
/// arithmetic. In the brief's form S·κ_ref·s, S = m·(1 + √ν) / ‖E(c_ref)‖₁: it depends on ν and the pattern.
/// </remarks>
internal sealed class SamplerAdmission
{
    /// <summary>The screen's margin for rounding, in the computed estimate and the computed system.</summary>
    private const double RoundingMargin = 2.0;

    /// <summary>The text <see cref="ConditionCheck.EveryFactorization"/>'s refusal has always given.</summary>
    internal const string Refused = "not positive definite, or its equilibrated condition past 1e14";

    private readonly double[]? _referenceWeights;
    private readonly double _screen;   // 2·m·((1 + √ν)/2)·‖E(c_ref)⁻¹‖₁; +∞ where the screen clears no step

    private SamplerAdmission(ConditionCheck check, string system, double[]? referenceWeights, double referenceCondition,
                             double screen)
    {
        Check = check;
        System = system;
        _referenceWeights = referenceWeights;
        ReferenceCondition = referenceCondition;
        _screen = screen;
    }

    public ConditionCheck Check { get; }

    public string System { get; }

    /// <summary>κ₁(E(c_ref)), read exactly; +∞ for a reference that is not positive definite, NaN where not read.</summary>
    public double ReferenceCondition { get; }

    /// <summary>S·κ_ref = 2·m·((1 + √ν)/2)·‖E(c_ref)⁻¹‖₁: a step at ratio s clears the screen when this times s is at most
    /// 10¹⁴.</summary>
    public double ScreenScale => _screen;

    /// <summary>c_ref, one per term, where the system's ratios are read against coefficients.</summary>
    public ReadOnlySpan<double> ReferenceWeights => _referenceWeights;

    /// <summary>
    /// A banded system's check. <paramref name="referenceWeights"/> are c_ref, one per term, aligned with the weights a step
    /// passes to <see cref="Spread(ReadOnlySpan{double})"/>; <paramref name="reference"/> builds A(c_ref) in lower-band
    /// storage at <paramref name="bandwidth"/>; <paramref name="pattern"/> is m. Under
    /// <see cref="ConditionCheck.ReferenceOnly"/> a reference that is not admissible throws here, before the first draw.
    /// </summary>
    public static SamplerAdmission Banded(ConditionCheck check, string system, double[]? referenceWeights, int dimension,
                                          int bandwidth, int pattern, Func<double[,]> reference)
    {
        if (check == ConditionCheck.EveryFactorization)
            return new SamplerAdmission(check, system, referenceWeights, double.NaN, double.PositiveInfinity);
        double[,] band = reference();
        var factor = new BandCholesky(dimension, bandwidth, BandFactorization.Ldlt);
        factor.DecomposeBanded(band);
        var diagonal = new double[dimension];
        for (int j = 0; j < dimension; j++) diagonal[j] = band[0, j];
        return Read(check, system, referenceWeights, factor, diagonal, BandNorm(band, dimension, bandwidth), pattern);
    }

    /// <summary>A dense system's check, as <see cref="Banded"/>; m is the order.</summary>
    public static SamplerAdmission Dense(ConditionCheck check, string system, int dimension, Func<double[,]> reference)
    {
        if (check == ConditionCheck.EveryFactorization)
            return new SamplerAdmission(check, system, null, double.NaN, double.PositiveInfinity);
        double[,] a = reference();
        var factor = new DenseCholesky(dimension);
        factor.Decompose(a);
        var diagonal = new double[dimension];
        for (int j = 0; j < dimension; j++) diagonal[j] = a[j, j];
        return Read(check, system, null, factor, diagonal, DenseNorm(a, dimension), dimension);
    }

    // The exact reading at the reference: ‖E‖₁ from A, ‖E⁻¹‖₁ from E⁻¹'s columns D^(½)A⁻¹D^(½)eⱼ, one solve each.
    private static SamplerAdmission Read(ConditionCheck check, string system, double[]? referenceWeights, ISpdFactor factor,
                                         double[] diagonal, double norm, int pattern)
    {
        int n = diagonal.Length;
        double inverseNorm = double.PositiveInfinity;
        if (factor.Status.PositiveDefinite && double.IsFinite(norm))
        {
            var root = new double[n];
            for (int j = 0; j < n; j++) root[j] = Math.Sqrt(diagonal[j]);
            var v = new double[n];
            inverseNorm = 0.0;
            for (int j = 0; j < n; j++)
            {
                Array.Clear(v);
                v[j] = root[j];
                factor.SolveInPlace(v);
                double column = 0.0;
                for (int i = 0; i < n; i++) column += Math.Abs(v[i]) * root[i];
                inverseNorm = Math.Max(inverseNorm, column);
            }
            if (!double.IsFinite(inverseNorm)) inverseNorm = double.PositiveInfinity;
        }
        double condition = double.IsFinite(inverseNorm) ? norm * inverseNorm : double.PositiveInfinity;
        if (check == ConditionCheck.ReferenceOnly && !(condition <= BapsAdmission.MaxCondition))
            throw new InvalidOperationException(
                $"{system} is not admissible at the trace-balance reference point: {(factor.Status.PositiveDefinite ? $"κ₁(E) = {condition:G3} past 1e14" : "not positive definite")} ({factor.Status}). The design leaves the penalty's null space unidentified, or nearly.");
        double screen = RoundingMargin * pattern * (0.5 * (1.0 + Math.Sqrt(n))) * inverseNorm;
        return new SamplerAdmission(check, system, referenceWeights, condition,
                                    double.IsFinite(screen) ? screen : double.PositiveInfinity);
    }

    /// <summary>max ρ / min ρ for a step's coefficients, aligned with the reference's; +∞ when a coefficient is 0, NaN
    /// when one is NaN.</summary>
    public double Spread(ReadOnlySpan<double> weights)
    {
        double[] reference = _referenceWeights!;
        double lo = double.PositiveInfinity, hi = 0.0;
        for (int k = 0; k < weights.Length; k++)
        {
            double r = weights[k] / reference[k];
            lo = Math.Min(lo, r);
            hi = Math.Max(hi, r);
        }
        return hi / lo;
    }

    /// <summary>Whether a step at <paramref name="spread"/> runs the estimate: always under
    /// <see cref="ConditionCheck.EveryFactorization"/>, never under <see cref="ConditionCheck.ReferenceOnly"/>, and where
    /// the screen does not clear it under <see cref="ConditionCheck.Screened"/>.</summary>
    public bool Estimates(double spread) => Check switch
    {
        ConditionCheck.EveryFactorization => true,
        ConditionCheck.Screened => !(_screen * spread <= BapsAdmission.MaxCondition),
        _ => false,
    };

    /// <summary>A step's check of its factor at <paramref name="spread"/>, recorded in <paramref name="tally"/>: null
    /// when admitted, else why it was refused.</summary>
    public string? Admit(ISpdFactor factor, double spread, ConditionTally tally)
    {
        tally.Factorizations++;
        tally.WidestRatio = Math.Max(tally.WidestRatio, spread);
        if (!factor.Status.PositiveDefinite)
            return Check == ConditionCheck.EveryFactorization ? Refused : "not positive definite";
        if (!Estimates(spread)) return null;
        double condition = factor.EquilibratedCondition();
        tally.LargestCondition = tally.Estimates++ == 0 ? condition : Math.Max(tally.LargestCondition, condition);
        if (condition <= BapsAdmission.MaxCondition) return null;
        return Check == ConditionCheck.EveryFactorization
            ? Refused
            : $"its equilibrated condition {condition:G3} past 1e14, where the reference point's is {ReferenceCondition:G3} and the coefficients' ratio to it {spread:G3}";
    }

    /// <summary>The run's record at this system, reduced over the chains' tallies.</summary>
    public ConditionRecord Record(IEnumerable<ConditionTally> tallies)
    {
        double widest = 1.0, largest = double.NaN;
        long factorizations = 0, estimates = 0;
        foreach (ConditionTally t in tallies)
        {
            widest = Math.Max(widest, t.WidestRatio);
            factorizations += t.Factorizations;
            if (t.Estimates > 0) largest = estimates == 0 ? t.LargestCondition : Math.Max(largest, t.LargestCondition);
            estimates += t.Estimates;
        }
        return new ConditionRecord(System, Check, ReferenceCondition, widest, factorizations, estimates, largest);
    }

    /// <summary>The trace balance of term k against the leading term: tr S_lead / tr Sₖ, or 1 where a trace is not
    /// positive.</summary>
    public static double Balance(double leading, double term) =>
        leading > 0.0 && term > 0.0 && double.IsFinite(leading / term) ? leading / term : 1.0;

    /// <summary>The trace of a lower-band matrix.</summary>
    public static double Trace(double[,] band, int dimension)
    {
        double t = 0.0;
        for (int j = 0; j < dimension; j++) t += band[0, j];
        return t;
    }

    /// <summary>m: the most entries in a column of A that any of <paramref name="bands"/> (lower-band storage, each
    /// with its own bandwidth) holds nonzero, the diagonal always counted; <paramref name="fullBands"/> count every entry
    /// within their bandwidth.</summary>
    public static int Pattern(int dimension, IEnumerable<double[,]> bands, IEnumerable<int>? fullBands = null)
    {
        int width = 0;
        var lists = new List<double[,]>(bands);
        foreach (double[,] b in lists) width = Math.Max(width, b.GetLength(0) - 1);
        var full = new List<int>(fullBands ?? []);
        foreach (int w in full) width = Math.Max(width, w);
        width = Math.Min(width, dimension - 1);
        var nonzero = new bool[width + 1, dimension];
        for (int j = 0; j < dimension; j++) nonzero[0, j] = true;
        foreach (double[,] b in lists)
            for (int d = 1; d < b.GetLength(0) && d <= width; d++)
                for (int j = 0; j + d < dimension; j++) nonzero[d, j] |= b[d, j] != 0.0;
        foreach (int w in full)
            for (int d = 1; d <= Math.Min(w, width); d++)
                for (int j = 0; j + d < dimension; j++) nonzero[d, j] = true;
        int most = 1;
        for (int j = 0; j < dimension; j++)
        {
            int count = 1;
            for (int d = 1; d <= width; d++)
            {
                if (j + d < dimension && nonzero[d, j]) count++;
                if (j - d >= 0 && nonzero[d, j - d]) count++;
            }
            most = Math.Max(most, count);
        }
        return most;
    }

    // ‖E‖₁ of a lower-band matrix, E its Jacobi scaling; NaN unless every diagonal entry is finite and positive.
    private static double BandNorm(double[,] band, int n, int p)
    {
        var root = new double[n];
        var sums = new double[n];
        for (int j = 0; j < n; j++)
        {
            double d = band[0, j];
            if (!(d > 0.0 && double.IsFinite(d))) return double.NaN;
            root[j] = Math.Sqrt(d);
            sums[j] = 1.0;
        }
        for (int j = 0; j < n; j++)
            for (int d = 1; d <= Math.Min(p, n - 1 - j); d++)
            {
                double e = Math.Abs(band[d, j]) / (root[j] * root[j + d]);
                sums[j] += e;
                sums[j + d] += e;
            }
        double norm = 0.0;
        foreach (double s in sums) norm = Math.Max(norm, s);
        return norm;
    }

    // ‖E‖₁ of a dense symmetric matrix.
    private static double DenseNorm(double[,] a, int n)
    {
        var root = new double[n];
        for (int j = 0; j < n; j++)
        {
            double d = a[j, j];
            if (!(d > 0.0 && double.IsFinite(d))) return double.NaN;
            root[j] = Math.Sqrt(d);
        }
        double norm = 0.0;
        for (int j = 0; j < n; j++)
        {
            double s = 0.0;
            for (int i = 0; i < n; i++) s += Math.Abs(a[i, j]) / (root[i] * root[j]);
            norm = Math.Max(norm, s);
        }
        return norm;
    }
}
