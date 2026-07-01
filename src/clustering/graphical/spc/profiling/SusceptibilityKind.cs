namespace Clustering.Graphical.SPC.Profiling;

/// <summary>
/// Selects which susceptibility channel drives chosen-T peak-finding in a
/// <see cref="SweepProfile"/>. All three are always assembled as channels (free comparison —
/// FK and magnetization estimators disagree only in edge cases, now testable programmatically);
/// this picks the primary <see cref="SweepProfile.Susceptibility"/> the analyzer's peak detector reads.
/// </summary>
public enum SusceptibilityKind
{
    /// <summary>FK cluster susceptibility ⟨Σ|c|²⟩/N — SW-native, lower variance. The default.</summary>
    FkCluster,

    /// <summary>
    /// Giant-excluded FK susceptibility ⟨Σ|c|² − |c_max|²⟩/N — the cleanest SP-transition
    /// detector (peaks vanish once the giant cluster percolates).
    /// </summary>
    FkReduced,

    /// <summary>
    /// Magnetization susceptibility (N/T)(⟨m²⟩ − ⟨m⟩²) — the literal BWD/Domany paper χ.
    /// Caution as a landmark driver: the N/T factor amplifies the residual
    /// fluctuations of weakly-coupled (background) spins into a cold-edge
    /// maximum on inhomogeneous data.
    /// </summary>
    Magnetization,

    /// <summary>
    /// Magnetization variance ⟨m²⟩ − ⟨m⟩² = χT/N — the curve the papers
    /// actually PLOT for landmark reading (WBD1998 Fig 8). Same peak /
    /// plateau / cliff semantics as χ with no 1/T cold-edge amplification —
    /// the replication landmark driver.
    /// </summary>
    MagnetizationVariance,
}
