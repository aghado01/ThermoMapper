using System;
using System.Collections.Generic;
using Maths.LinAlg;

namespace Maths.Regression.Spline.Bass;

/// <summary>Which truncation residuals a decomposition's predictive draws add (D07).</summary>
public enum TruncationResidual
{
    /// <summary>
    /// D07's default: each case's residual against the rank-r basis built without it, its rank capped at the other
    /// cases' numerical rank, so it stands for a new case's error outside the basis.
    /// </summary>
    LeaveOneCaseOut,

    /// <summary>R's residual: each case against the basis it helped build, which understates a new case's error and is
    /// zero at full retention.</summary>
    InSample,

    /// <summary>Residuals supplied with a basis, in its standardized units.</summary>
    Supplied,
}

/// <summary>Where a reconstruction map can be evaluated (O-X5).</summary>
public enum ReconstructionDomain
{
    /// <summary>Only at the training output locations, as a PCA or supplied basis is.</summary>
    TrainingLocations,

    /// <summary>At new output coordinates too, as an analytic basis would be (an extension; no map in the port has it).</summary>
    AnyLocation,
}

/// <summary>
/// A response representation's reconstruction map (O-X5): coefficients w and a truncation residual e to an output
/// vector y. The linear map (<see cref="LinearReconstruction"/>) is the port's only member; elastic emulation's
/// composition joins as another.
/// </summary>
public abstract class ReconstructionMap
{
    private protected ReconstructionMap() { }

    /// <summary>The output locations m.</summary>
    public abstract int Locations { get; }

    /// <summary>The coefficients r.</summary>
    public abstract int Coefficients { get; }

    /// <summary>Where the map can be evaluated.</summary>
    public abstract ReconstructionDomain Domain { get; }

    /// <summary>y from <paramref name="w"/> and the residual <paramref name="e"/> (empty for none), into <paramref name="y"/>.</summary>
    public abstract void Reconstruct(ReadOnlySpan<double> w, ReadOnlySpan<double> e, Span<double> y);
}

/// <summary>
/// D07's linear map, y = y_m + y_s ∘ (K w + e): the offsets y_m, the scales y_s and the basis K (m × r), which may be
/// nonorthogonal. Evaluated at the training locations only. Immutable: the arrays are copied in.
/// </summary>
public sealed class LinearReconstruction : ReconstructionMap
{
    private readonly double[,] _basis;
    private readonly double[] _offsets, _scales;

    /// <summary>The map of <paramref name="basis"/> (m × r), <paramref name="offsets"/> and positive
    /// <paramref name="scales"/> (m each).</summary>
    public LinearReconstruction(double[,] basis, double[] offsets, double[] scales)
    {
        ArgumentNullException.ThrowIfNull(basis);
        ArgumentNullException.ThrowIfNull(offsets);
        ArgumentNullException.ThrowIfNull(scales);
        int m = basis.GetLength(0), r = basis.GetLength(1);
        if (m < 1 || r < 1) throw new ArgumentException("A basis needs a location and a component.", nameof(basis));
        foreach (double v in basis)
            if (!double.IsFinite(v)) throw new ArgumentException("The basis has a non-finite entry.", nameof(basis));
        if (offsets.Length != m) throw new ArgumentException($"The offsets need one value per location ({m}).", nameof(offsets));
        if (scales.Length != m) throw new ArgumentException($"The scales need one value per location ({m}).", nameof(scales));
        foreach (double v in offsets)
            if (!double.IsFinite(v)) throw new ArgumentException("The offsets must be finite.", nameof(offsets));
        foreach (double v in scales)
            if (!(v > 0.0) || double.IsPositiveInfinity(v)) throw new ArgumentException("The scales must be positive and finite.", nameof(scales));
        _basis = (double[,])basis.Clone();
        _offsets = (double[])offsets.Clone();
        _scales = (double[])scales.Clone();
    }

    /// <inheritdoc/>
    public override int Locations => _offsets.Length;

    /// <inheritdoc/>
    public override int Coefficients => _basis.GetLength(1);

    /// <inheritdoc/>
    public override ReconstructionDomain Domain => ReconstructionDomain.TrainingLocations;

    /// <summary>y_m.</summary>
    public ReadOnlySpan<double> Offsets => _offsets;

    /// <summary>y_s.</summary>
    public ReadOnlySpan<double> Scales => _scales;

    /// <summary>K at <paramref name="location"/> and <paramref name="component"/>.</summary>
    public double Basis(int location, int component) => _basis[location, component];

    /// <summary>A copy of K.</summary>
    public double[,] BasisMatrix() => (double[,])_basis.Clone();

    internal double[,] BasisArray => _basis;

    /// <inheritdoc/>
    public override void Reconstruct(ReadOnlySpan<double> w, ReadOnlySpan<double> e, Span<double> y)
    {
        int m = Locations, r = Coefficients;
        if (w.Length != r) throw new ArgumentException($"The map takes {r} coefficients.", nameof(w));
        if (!e.IsEmpty && e.Length != m) throw new ArgumentException($"A residual has one value per location ({m}).", nameof(e));
        if (y.Length != m) throw new ArgumentException($"The output has one value per location ({m}).", nameof(y));
        for (int l = 0; l < m; l++)
        {
            double s = 0.0;
            for (int k = 0; k < r; k++) s += _basis[l, k] * w[k];
            if (!e.IsEmpty) s += e[l];
            y[l] = _offsets[l] + _scales[l] * s;
        }
    }

    /// <summary>A covariance over the locations in standardized units to physical ones: diag(y_s)·C·diag(y_s) (D07).</summary>
    public double[,] Physical(double[,] standardized)
    {
        ArgumentNullException.ThrowIfNull(standardized);
        int m = Locations;
        if (standardized.GetLength(0) != m || standardized.GetLength(1) != m)
            throw new ArgumentException($"The covariance must be {m} × {m}.", nameof(standardized));
        var physical = new double[m, m];
        for (int l = 0; l < m; l++)
            for (int j = 0; j < m; j++) physical[l, j] = _scales[l] * standardized[l, j] * _scales[j];
        return physical;
    }
}

/// <summary>
/// The slot for basis-estimation uncertainty (O-X3): no member yet, since the baseline holds the basis fixed (D07's
/// "not in the baseline"); a model of it, such as Fast Bayesian FPCA's, joins here.
/// </summary>
public abstract record BasisUncertainty;

/// <summary>
/// A decomposition's truncation residuals (D07), each case's e_i = y_{c,i} − K wᵢ, observed minus reconstructed, in
/// standardized units, one row per case. <see cref="Values"/> are the ones predictive draws add, of the declared
/// <see cref="Kind"/>; a PCA decomposition keeps both sets beside them, and their RMS ratio, leave-one-case-out over in
/// sample, is the basis-stability diagnostic: near 1 when the basis is settled, large as r approaches n − 1.
/// </summary>
public sealed class TruncationResiduals
{
    internal TruncationResiduals(TruncationResidual kind, double[,] values, double[,]? inSample, double[,]? leaveOneOut,
                                 int[]? leaveOneOutRanks, double tolerance)
    {
        Kind = kind;
        Values = values;
        InSample = inSample;
        LeaveOneOut = leaveOneOut;
        LeaveOneOutRanks = leaveOneOutRanks;
        Rms = RootMeanSquare(values);
        InSampleRms = inSample is null ? null : RootMeanSquare(inSample);
        LeaveOneOutRms = leaveOneOut is null ? null : RootMeanSquare(leaveOneOut);
        if (inSample is not null && leaveOneOut is not null && MaxAbs(inSample) > tolerance)
            StabilityRatio = LeaveOneOutRms / InSampleRms;
    }

    /// <summary>Which residuals <see cref="Values"/> holds.</summary>
    public TruncationResidual Kind { get; }

    /// <summary>The residuals predictive draws add: n × m, case i's in row i, in standardized units.</summary>
    public double[,] Values { get; }

    /// <summary>The in-sample residuals; null when the residuals were supplied.</summary>
    public double[,]? InSample { get; }

    /// <summary>The leave-one-case-out residuals; null for a supplied basis, which has no construction to leave a case out of.</summary>
    public double[,]? LeaveOneOut { get; }

    /// <summary>Each case's leave-one-out rank, min(r, the other cases' numerical rank).</summary>
    public int[]? LeaveOneOutRanks { get; }

    /// <summary>The RMS of <see cref="Values"/>.</summary>
    public double Rms { get; }

    /// <summary>The in-sample residuals' RMS.</summary>
    public double? InSampleRms { get; }

    /// <summary>The leave-one-out residuals' RMS.</summary>
    public double? LeaveOneOutRms { get; }

    /// <summary>Leave-one-out RMS over in-sample RMS; null when either is absent or the in-sample residual is zero to the
    /// rank tolerance, as at full retention.</summary>
    public double? StabilityRatio { get; }

    /// <summary>The cases n.</summary>
    public int Cases => Values.GetLength(0);

    /// <summary>The residuals' mean ē, the m-vector a uniform draw over the cases averages to.</summary>
    public double[] Mean()
    {
        int n = Cases, m = Values.GetLength(1);
        var mean = new double[m];
        var column = new double[n];
        for (int l = 0; l < m; l++)
        {
            for (int i = 0; i < n; i++) column[i] = Values[i, l];
            mean[l] = PairwiseSum.Of(column) / n;
        }
        return mean;
    }

    /// <summary>The covariance of e drawn uniformly over the cases, (1/n)·Σᵢ (eᵢ − ē)(eᵢ − ē)ᵀ, m × m in standardized units.</summary>
    public double[,] Covariance()
    {
        int n = Cases, m = Values.GetLength(1);
        double[] mean = Mean();
        var centred = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++) centred[i, l] = Values[i, l] - mean[l];
        var covariance = new double[m, m];
        for (int l = 0; l < m; l++)
            for (int j = l; j < m; j++)
            {
                double s = 0.0;
                for (int i = 0; i < n; i++) s += centred[i, l] * centred[i, j];
                covariance[l, j] = covariance[j, l] = s / n;
            }
        return covariance;
    }

    private static double RootMeanSquare(double[,] e)
    {
        double s = 0.0;
        foreach (double v in e) s += v * v;
        return Math.Sqrt(s / e.Length);
    }

    private static double MaxAbs(double[,] e)
    {
        double max = 0.0;
        foreach (double v in e) max = Math.Max(max, Math.Abs(v));
        return max;
    }
}

/// <summary>
/// How a PCA decomposition is built (D07; R's <c>bassPCA</c>): the centre y_m (null: the column means), the scale y_s
/// (null: 1; <see cref="ScaleByColumnSd"/> takes the columns' standard deviations, a constant column's 1), the rank (the
/// smallest whose cumulative variance share reaches <see cref="VarianceShare"/>, or <see cref="Rank"/>, each capped at the
/// numerical rank) and the truncation residuals predictive draws add.
/// </summary>
public sealed record BassPcaSpec
{
    /// <summary>The variance share the rank must reach, in (0, 1]; 1 keeps the numerical rank. R's default, 99%.</summary>
    public double VarianceShare { get; init; } = 0.99;

    /// <summary>An explicit rank, overriding <see cref="VarianceShare"/>; capped at the numerical rank.</summary>
    public int? Rank { get; init; }

    /// <summary>y_m, one value per location; null takes the column means.</summary>
    public double[]? Centre { get; init; }

    /// <summary>Whether y_s is each column's standard deviation (n − 1 denominator), 1 for a constant column.</summary>
    public bool ScaleByColumnSd { get; init; }

    /// <summary>y_s, one positive value per location; null takes 1, or the standard deviations.</summary>
    public double[]? Scale { get; init; }

    /// <summary>The residuals predictive draws add: leave-one-case-out by default (D07), or <see cref="TruncationResidual.InSample"/>.</summary>
    public TruncationResidual Residual { get; init; } = TruncationResidual.LeaveOneCaseOut;
}

/// <summary>
/// A PCA decomposition's account: every singular value of the standardized responses, the numerical rank ρ and the
/// tolerance that judged it, max(n, m)·ε·d₁, the variance share asked for and reached at r, whether ρ capped the rank,
/// and the constant output columns.
/// </summary>
public sealed record BassPcaAccount(double[] SingularValues, int NumericalRank, double RankTolerance, double? VarianceShare,
                                    double ShareReached, bool RankCapped, int[] ConstantColumns);

/// <summary>
/// A response decomposition (D07; R's <c>bassPCA</c> and <c>bassBasis</c>): the representation of n cases' m-vectors as
/// coefficient responses W (r × n) under a reconstruction map, y = y_m + y_s ∘ (K w + e), with the truncation residuals
/// e. Each row of W is one component's coefficient response, fitted as a scalar response on the shared predictors.
/// <list type="bullet">
/// <item><b>PCA</b> (<see cref="FromPca"/>): the thin SVD of the standardized responses, Y_cᵀ = U D Vᵀ, with K = U_r D_r and
/// W = V_rᵀ, R's convention, whose rows have unit norm; the rank as <see cref="BassPcaSpec"/> asks, capped at the numerical
/// rank judged at max(n, m)·ε·d₁; each component's sign fixed so that its largest entry of K is positive, a convention
/// for comparisons only. Where singular values repeat only the subspace is defined, so compare reconstructions.</item>
/// <item><b>A supplied basis</b> (<see cref="FromBasis"/>): K, nonorthogonal allowed, with prepared coefficient responses,
/// offsets and scales. The coefficients are not solved for here.</item>
/// <item><b>Truncation residuals</b> (<see cref="TruncationResiduals"/>): leave-one-case-out by default, computed in case
/// coordinates from the one SVD at O(n⁴), independent of m, each case's rank capped at the other cases' numerical rank
/// judged at the same tolerance; the in-sample residual is the option. The scales stay the full data's, so every
/// residual is in the decomposition's standardized units.</item>
/// </list>
/// <see cref="Uncertainty"/> is the slot for basis-estimation uncertainty, null in the baseline (O-X3). Immutable.
/// </summary>
public sealed class BassDecomposition
{
    private readonly double[,] _coefficients;

    private BassDecomposition(LinearReconstruction map, double[,] coefficients, TruncationResiduals truncation, BassPcaAccount? pca)
    {
        Map = map;
        _coefficients = coefficients;
        Truncation = truncation;
        Pca = pca;
    }

    /// <summary>The reconstruction map.</summary>
    public LinearReconstruction Map { get; }

    /// <summary>The components r.</summary>
    public int Rank => Map.Coefficients;

    /// <summary>The cases n.</summary>
    public int Cases => _coefficients.GetLength(1);

    /// <summary>The output locations m.</summary>
    public int Locations => Map.Locations;

    /// <summary>W at component <paramref name="component"/> and case <paramref name="case"/>.</summary>
    public double Coefficient(int component, int @case) => _coefficients[component, @case];

    /// <summary>Component <paramref name="component"/>'s coefficient response, a copy of W's row.</summary>
    public double[] CoefficientResponse(int component)
    {
        var w = new double[Cases];
        for (int i = 0; i < w.Length; i++) w[i] = _coefficients[component, i];
        return w;
    }

    /// <summary>The truncation residuals.</summary>
    public TruncationResiduals Truncation { get; }

    /// <summary>Basis-estimation uncertainty: null, the basis held fixed (the baseline).</summary>
    public BasisUncertainty? Uncertainty => null;

    /// <summary>The PCA account; null for a supplied basis.</summary>
    public BassPcaAccount? Pca { get; }

    /// <summary>The truncation covariance in physical units, diag(y_s)·Cov(e)·diag(y_s), e drawn uniformly over the cases.</summary>
    public double[,] TruncationCovariance() => Map.Physical(Truncation.Covariance());

    /// <summary>
    /// The PCA decomposition of <paramref name="y"/> (n × m, a case per row) under <paramref name="spec"/>. Refused: fewer
    /// than two cases, a non-finite value, a centre or scale of the wrong length, a scale not positive and finite, both
    /// scale options at once, a share outside (0, 1], a rank below 1, and standardized responses that are zero to the
    /// rank tolerance, which leave no component to fit.
    /// </summary>
    public static BassDecomposition FromPca(double[,] y, BassPcaSpec? spec = null)
    {
        ArgumentNullException.ThrowIfNull(y);
        spec ??= new BassPcaSpec();
        int n = y.GetLength(0), m = y.GetLength(1);
        if (n < 2) throw new ArgumentException("A decomposition needs at least two cases.", nameof(y));
        if (m < 1) throw new ArgumentException("A decomposition needs an output location.", nameof(y));
        foreach (double v in y)
            if (!double.IsFinite(v)) throw new ArgumentException("The responses have a non-finite value.", nameof(y));
        if (!(spec.VarianceShare > 0.0 && spec.VarianceShare <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(spec), "The variance share lies in (0, 1].");
        if (spec.Rank is < 1) throw new ArgumentOutOfRangeException(nameof(spec), "An explicit rank is at least 1.");
        if (spec.ScaleByColumnSd && spec.Scale is not null)
            throw new ArgumentException("Give the scales or ask for the columns' standard deviations, not both.", nameof(spec));
        if (spec.Residual == TruncationResidual.Supplied)
            throw new ArgumentException("A PCA decomposition builds its own residuals: leave-one-case-out or in sample.", nameof(spec));

        var constant = new List<int>();
        for (int l = 0; l < m; l++)
        {
            double min = y[0, l], max = y[0, l];
            for (int i = 1; i < n; i++)
            {
                min = Math.Min(min, y[i, l]);
                max = Math.Max(max, y[i, l]);
            }
            if (min == max) constant.Add(l);
        }
        double[] offsets = spec.Centre is { } centre ? Checked(centre, m, "centre", positive: false) : ColumnMeans(y);
        double[] scales = spec.Scale is { } scale ? Checked(scale, m, "scale", positive: true)
                        : spec.ScaleByColumnSd ? ColumnSds(y, offsets: ColumnMeans(y), constant) : Ones(m);

        var yc = new double[m, n];   // Y_cᵀ: a location per row, a case per column
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++) yc[l, i] = (y[i, l] - offsets[l]) / scales[l];
        ThinSvdResult svd = ThinSvd.Decompose(yc);
        int rho = svd.NumericalRank;
        if (rho == 0) throw new ArgumentException("The standardized responses are zero to working precision: no component to fit.", nameof(y));
        double[] d = svd.SingularValues;

        double total = 0.0;
        foreach (double v in d) total += v * v;
        int r;
        if (spec.Rank is { } explicitRank) r = Math.Min(explicitRank, rho);
        else
        {
            // The smallest r whose share reaches the target (≥, not R's strict >, which finds none at 100%), capped at ρ;
            // full retention is ρ itself, whatever rounding does to the last shares.
            r = rho;
            double cumulative = 0.0;
            for (int k = 0; k < rho && spec.VarianceShare < 1.0; k++)
            {
                cumulative += d[k] * d[k];
                if (cumulative / total >= spec.VarianceShare)
                {
                    r = k + 1;
                    break;
                }
            }
        }
        double reached = 0.0;
        for (int k = 0; k < r; k++) reached += d[k] * d[k];
        bool capped = spec.Rank is { } asked ? asked > rho : spec.VarianceShare < 1.0 && r == rho && reached / total < spec.VarianceShare;

        // The sign convention, for comparisons only: each component's largest entry of K is positive.
        double[][] left = svd.LeftVectors, right = svd.RightVectors;
        for (int k = 0; k < rho; k++)
        {
            int largest = 0;
            for (int l = 1; l < m; l++)
                if (Math.Abs(left[k][l]) > Math.Abs(left[k][largest])) largest = l;
            if (left[k][largest] >= 0.0) continue;
            for (int l = 0; l < m; l++) left[k][l] = -left[k][l];
            for (int i = 0; i < n; i++) right[k][i] = -right[k][i];
        }
        var basis = new double[m, r];
        var coefficients = new double[r, n];
        for (int k = 0; k < r; k++)
        {
            for (int l = 0; l < m; l++) basis[l, k] = left[k][l] * d[k];
            for (int i = 0; i < n; i++) coefficients[k, i] = right[k][i];
        }

        double[,] inSample = InSampleResiduals(yc, basis, coefficients);
        double[,] leaveOneOut = LeaveOneOutResiduals(svd, r, centredByMean: spec.Centre is null, out int[] ranks);
        double[,] chosen = spec.Residual == TruncationResidual.InSample ? inSample : leaveOneOut;
        var truncation = new TruncationResiduals(spec.Residual, chosen, inSample, leaveOneOut, ranks, svd.RankTolerance);
        var account = new BassPcaAccount((double[])d.Clone(), rho, svd.RankTolerance, spec.Rank is null ? spec.VarianceShare : null,
                                         reached / total, capped, [.. constant]);
        return new BassDecomposition(new LinearReconstruction(basis, offsets, scales), coefficients, truncation, account);
    }

    /// <summary>
    /// A supplied basis (R's <c>bassBasis</c>): K = <paramref name="basis"/> (m × r, nonorthogonal allowed), the prepared
    /// coefficient responses W = <paramref name="coefficients"/> (r × n), <paramref name="offsets"/> and
    /// <paramref name="scales"/>, with the truncation residuals from the training <paramref name="responses"/> (n × m),
    /// e_i = (y_i − y_m)/y_s − K wᵢ, which a basis fixed in advance did not see, or <paramref name="residuals"/> (n × m,
    /// standardized) as supplied: exactly one of the two.
    /// </summary>
    public static BassDecomposition FromBasis(double[,] basis, double[,] coefficients, double[] offsets, double[] scales,
                                              double[,]? responses = null, double[,]? residuals = null)
    {
        ArgumentNullException.ThrowIfNull(coefficients);
        var map = new LinearReconstruction(basis, offsets, scales);
        int m = map.Locations, r = map.Coefficients, n = coefficients.GetLength(1);
        if (coefficients.GetLength(0) != r) throw new ArgumentException($"The coefficients need one row per component ({r}).", nameof(coefficients));
        if (n < 2) throw new ArgumentException("A decomposition needs at least two cases.", nameof(coefficients));
        foreach (double w in coefficients)
            if (!double.IsFinite(w)) throw new ArgumentException("The coefficients must be finite.", nameof(coefficients));
        if ((responses is null) == (residuals is null))
            throw new ArgumentException("Give the training responses or the residuals, exactly one of them.", nameof(responses));
        double[,] given = responses ?? residuals!;
        if (given.GetLength(0) != n || given.GetLength(1) != m)
            throw new ArgumentException($"The {(responses is null ? "residuals" : "responses")} must be {n} × {m}, a case per row.",
                                        responses is null ? nameof(residuals) : nameof(responses));
        foreach (double v in given)
            if (!double.IsFinite(v)) throw new ArgumentException("The responses and residuals must be finite.", nameof(responses));

        TruncationResiduals truncation;
        if (responses is not null)
        {
            var yc = new double[m, n];
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++) yc[l, i] = (responses[i, l] - map.Offsets[l]) / map.Scales[l];
            double[,] inSample = InSampleResiduals(yc, map.BasisArray, coefficients);
            truncation = new TruncationResiduals(TruncationResidual.InSample, inSample, inSample, null, null, 0.0);
        }
        else truncation = new TruncationResiduals(TruncationResidual.Supplied, (double[,])residuals!.Clone(), null, null, null, 0.0);
        return new BassDecomposition(map, (double[,])coefficients.Clone(), truncation, null);
    }

    // e_i = y_{c,i} − K wᵢ, in output space, so the training reconstruction holds to rounding.
    private static double[,] InSampleResiduals(double[,] yc, double[,] basis, double[,] coefficients)
    {
        int m = yc.GetLength(0), n = yc.GetLength(1), r = basis.GetLength(1);
        var e = new double[n, m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++)
            {
                double s = 0.0;
                for (int k = 0; k < r; k++) s += basis[l, k] * coefficients[k, i];
                e[i, l] = yc[l, i] - s;
            }
        return e;
    }

    // The leave-one-case-out residuals in case coordinates. Y_cᵀ = U_ρ C with C = D_ρ V_ρᵀ, and every quantity lies in
    // span(U_ρ). Under the mean centre, case i centred by the others' mean is n/(n−1)·cᵢ and case j is cⱼ + cᵢ/(n−1);
    // under a fixed centre they are cᵢ and cⱼ. Each refit is the SVD of the others' coordinates (ρ × (n − 1)), its rank
    // min(r, the count above the full data's tolerance), so a zero matrix (n = 2) has rank 0.
    private static double[,] LeaveOneOutResiduals(ThinSvdResult svd, int r, bool centredByMean, out int[] ranks)
    {
        int rho = svd.NumericalRank, m = svd.LeftVectors[0].Length, n = svd.RightVectors[0].Length;
        double tolerance = svd.RankTolerance;
        var c = new double[rho, n];
        for (int k = 0; k < rho; k++)
            for (int i = 0; i < n; i++) c[k, i] = svd.SingularValues[k] * svd.RightVectors[k][i];

        var e = new double[n, m];
        ranks = new int[n];
        var z = new double[rho, n - 1];
        var x = new double[rho];
        for (int i = 0; i < n; i++)
        {
            double shift = centredByMean ? 1.0 / (n - 1) : 0.0;
            for (int k = 0; k < rho; k++)
            {
                for (int j = 0, col = 0; j < n; j++)
                    if (j != i) z[k, col++] = c[k, j] + c[k, i] * shift;
                x[k] = c[k, i] * (centredByMean ? (double)n / (n - 1) : 1.0);
            }
            ThinSvdResult others = ThinSvd.Decompose(z);
            int below = 0;
            while (below < others.SingularValues.Length && others.SingularValues[below] > tolerance) below++;
            int rank = Math.Min(r, below);
            ranks[i] = rank;
            for (int q = 0; q < rank; q++)
            {
                double[] direction = others.LeftVectors[q];
                double projection = 0.0;
                for (int k = 0; k < rho; k++) projection += direction[k] * x[k];
                for (int k = 0; k < rho; k++) x[k] -= projection * direction[k];
            }
            for (int l = 0; l < m; l++)
            {
                double s = 0.0;
                for (int k = 0; k < rho; k++) s += svd.LeftVectors[k][l] * x[k];
                e[i, l] = s;
            }
        }
        return e;
    }

    private static double[] ColumnMeans(double[,] y)
    {
        int n = y.GetLength(0), m = y.GetLength(1);
        var means = new double[m];
        var column = new double[n];
        for (int l = 0; l < m; l++)
        {
            for (int i = 0; i < n; i++) column[i] = y[i, l];
            means[l] = PairwiseSum.Of(column) / n;
        }
        return means;
    }

    // R's sd, the n − 1 denominator; a constant column gets 1 (D07), never a division by zero.
    private static double[] ColumnSds(double[,] y, double[] offsets, List<int> constant)
    {
        int n = y.GetLength(0), m = y.GetLength(1);
        var sds = new double[m];
        var squares = new double[n];
        for (int l = 0; l < m; l++)
        {
            if (constant.Contains(l))
            {
                sds[l] = 1.0;
                continue;
            }
            for (int i = 0; i < n; i++) squares[i] = (y[i, l] - offsets[l]) * (y[i, l] - offsets[l]);
            sds[l] = Math.Sqrt(PairwiseSum.Of(squares) / (n - 1));
        }
        return sds;
    }

    private static double[] Ones(int m)
    {
        var ones = new double[m];
        Array.Fill(ones, 1.0);
        return ones;
    }

    private static double[] Checked(double[] values, int m, string what, bool positive)
    {
        if (values.Length != m) throw new ArgumentException($"The {what} needs one value per location ({m}).", "spec");
        foreach (double v in values)
            if (!double.IsFinite(v) || (positive && !(v > 0.0)))
                throw new ArgumentException($"The {what} must be finite{(positive ? " and positive" : "")}.", "spec");
        return (double[])values.Clone();
    }
}
