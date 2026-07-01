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

# 1. Turn R on (single entry point; R is opt-in, not in the default env).
. "$env:PORTABLE_ROOT/UserGithub/PowerShellCore/ps.core.bootstrap/helpers/env-Rlang.ps1"

if (-not (Test-Path -LiteralPath $env:R_HOME)) {
    Write-Host "oracle-ci: SKIP - R toolchain not found at $env:R_HOME (provision with scripts/bootstrap.R)." -ForegroundColor Yellow
    exit 0
}

Push-Location -LiteralPath $rRoot
try {
    Invoke-RVersion

    if ($Restore) { Invoke-RenvRestore }
    Invoke-RenvStatus   # reproducibility: warns if the library drifts from renv.lock

    # Smoke: a base-R oracle round-trips on a tiny fixture.
    $fixture = Join-Path $env:TEMP 'oracle_smoke.csv'
    @('1,2', '2,1', '3,5', '5,3', '4,4') | Set-Content -LiteralPath $fixture
    $out = Join-Path $env:TEMP 'oracle_smoke.json'
    Remove-Item -LiteralPath $out -ErrorAction SilentlyContinue
    Invoke-Rscript 'oracles/pca_oracle.R' $fixture $out 2

    if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $out)) {
        Write-Host 'oracle-ci: FAIL - pca_oracle smoke produced no output.' -ForegroundColor Red
        exit 1
    }
    Write-Host 'oracle-ci: PASS - R on, renv checked, oracle smoke green.' -ForegroundColor Green
}
finally { Pop-Location }
exit 0
