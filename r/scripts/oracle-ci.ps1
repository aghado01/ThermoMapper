<#
.SYNOPSIS
oracle-ci.ps1 — health & reproducibility gate for the R validation oracle.

.DESCRIPTION
Turns R on (the portable toolchain at $PORTABLE_ROOT/rlang), checks the renv
project library is in sync with renv.lock, and smoke-runs a base-R oracle so a
broken bridge fails loudly. The actual parity assertions live in the C# tests,
which reach the same toolchain.

  -Restore   reproduce the pinned library first (renv::restore).

Mirrors lean/scripts/meta-ci.ps1: the in-tree gate for an opt-in toolchain.
#>
param([switch]$Restore)
$ErrorActionPreference = 'Stop'
$rRoot = Split-Path -Parent $PSScriptRoot

# 1. Turn R on. This first honors live/user-level R env vars, then falls back to
# the PDenv layout used by the portable toolchain.
. (Join-Path $PSScriptRoot 'r-session.ps1')
$toolchain = Initialize-ROracleSession -SetAliases

if (-not $toolchain -or -not (Test-Path -LiteralPath $toolchain.RscriptExe -PathType Leaf)) {
    Write-Host "oracle-ci: SKIP - R toolchain not found (set R_HOME, put Rscript on PATH, or install PDenv R under ~/PDenv/rlang)." -ForegroundColor Yellow
    exit 0
}

Push-Location -LiteralPath $rRoot
try {
    Invoke-RscriptChecked -TimeoutSeconds 30 -ArgumentList @('--vanilla', '-e', "cat(R.version.string, '\n')") | Out-Null

    if ($Restore) {
        Invoke-RscriptChecked -TimeoutSeconds 900 -ArgumentList @('-e', 'renv::restore(prompt = FALSE)') | Out-Null
    }
    Invoke-RscriptChecked -TimeoutSeconds 45 -ArgumentList @('-e', 'renv::status()') | Out-Null

    # Smoke: a base-R oracle round-trips on a tiny fixture.
    $fixture = Join-Path $env:TEMP 'oracle_smoke.csv'
    @('1,2', '2,1', '3,5', '5,3', '4,4') | Set-Content -LiteralPath $fixture
    $out = Join-Path $env:TEMP 'oracle_smoke.json'
    Remove-Item -LiteralPath $out -ErrorAction SilentlyContinue
    Invoke-RscriptChecked -TimeoutSeconds 120 -ArgumentList @('oracles/pca_oracle.R', $fixture, $out, '2') | Out-Null

    if (-not (Test-Path -LiteralPath $out)) {
        Write-Host 'oracle-ci: FAIL - pca_oracle smoke produced no output.' -ForegroundColor Red
        exit 1
    }
    Write-Host 'oracle-ci: PASS - R on, renv checked, oracle smoke green.' -ForegroundColor Green
}
finally { Pop-Location }
exit 0
