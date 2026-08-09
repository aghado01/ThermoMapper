<#
.SYNOPSIS
meta-ci.ps1 — taxonomy gate for the Lean rigor harness.

.DESCRIPTION
Active path: Protolemmata/ (markdown briefs) -> Enthymemata/ (statements
compile, proofs may still apologize) -> Lemmas/ (verified reusable results),
with optional promotion to Theorems/ for consequential verified deliverables.
Archeion/ is the non-building side exit for superseded material.

Invariants enforced:
  1. `lake build` is green — every active formal stage must compile; an enthymema owes
     proofs, never statements.
  2. Each aggregate module imports every Lean file in its active stage; no file
     can silently evade the build.
  3. Active Lean source uses scoped Mathlib modules; the `import Mathlib`
     umbrella is forbidden.
  4. The Lemmas and Theorems stages never apologize: no `sorry` token in
     their Lean files and no `import Enthymemata.*` (which could launder a
     sorried dependency). Lemmas also cannot import upward from Theorems.
  5. Enthymema ledger, per file:
       unstated         — no declarations yet (a prelemma in costume)
       apologizing(n)   — n sorries outstanding
       PROOF-CLOSED     — declarations present, zero sorries: review the
                          statement and dependencies for promotion to Lemmas/.

Exit 1 on any invariant violation or red build. With -Validate, ledger
notices (unstated / proof-closed) also fail — taxonomy drift is an error
in strict mode. -NoBuild skips the compile gate (for use after lean-action
has already built, e.g. in CI).
#>
param(
    [switch]$Validate,
    [switch]$NoBuild
)

$ErrorActionPreference = 'Stop'
$leanRoot = Split-Path -Parent $PSScriptRoot

if (-not $NoBuild) {
    $portableRootPath = "$env:PORTABLE_ROOT".Trim().Replace('/', '\')
    $lakeExecutable = $null

    if (-not [string]::IsNullOrWhiteSpace($portableRootPath)) {
        $portableElanHome = Join-Path $portableRootPath 'elan'
        $portableLake = Join-Path $portableElanHome 'bin\lake.exe'
        if (Test-Path -LiteralPath $portableLake -PathType Leaf) {
            $lakeExecutable = $portableLake
            if ([string]::IsNullOrWhiteSpace($env:ELAN_HOME)) {
                $env:ELAN_HOME = $portableElanHome
            }
        }
    }

    if ($null -eq $lakeExecutable) {
        $lakeCommand = Get-Command lake -CommandType Application -ErrorAction SilentlyContinue |
            Select-Object -First 1
        if ($null -ne $lakeCommand) {
            $lakeExecutable = $lakeCommand.Source
        }
    }

    if ($null -eq $lakeExecutable) {
        throw 'Lake is not available on PATH or under PORTABLE_ROOT\elan\bin.'
    }

    Push-Location -LiteralPath $leanRoot
    try {
        & $lakeExecutable build
        if ($LASTEXITCODE -ne 0) {
            Write-Host 'meta-ci: FAIL — build is red.' -ForegroundColor Red
            exit 1
        }
    }
    finally { Pop-Location }
}

$violations = @()
$notices    = @()

# --- Invariant 2: aggregate modules cover every active Lean file ------------
foreach ($tier in @('Lemmas', 'Theorems', 'Enthymemata')) {
    $tierRoot = Join-Path $leanRoot $tier
    $aggregatePath = Join-Path $leanRoot "$tier.lean"
    $aggregateImports = @(
        Get-Content -LiteralPath $aggregatePath | ForEach-Object {
            if ($_ -match '^\s*import\s+(\S+)') { $Matches[1] }
        }
    )

    foreach ($f in Get-ChildItem -LiteralPath $tierRoot -Filter *.lean -Recurse) {
        $relative = [IO.Path]::GetRelativePath($tierRoot, $f.FullName)
        $moduleSuffix = ($relative -replace '\.lean$', '') -replace '[\\/]', '.'
        $moduleName = "$tier.$moduleSuffix"
        if ($moduleName -notin $aggregateImports) {
            $violations += "$tier aggregate omits active module: $moduleName"
        }
    }
}

# --- Invariant 3: active source uses scoped Mathlib imports -----------------
foreach ($tier in @('Lemmas', 'Theorems', 'Enthymemata')) {
    $activeFiles = @(
        Get-ChildItem -LiteralPath (Join-Path $leanRoot $tier) -Filter *.lean -Recurse
        Get-Item -LiteralPath (Join-Path $leanRoot "$tier.lean")
    )
    foreach ($f in $activeFiles) {
        foreach ($m in @(Select-String -LiteralPath $f.FullName -Pattern '^\s*import\s+Mathlib(?:\s|$)')) {
            $violations += "$tier uses the Mathlib umbrella: $($f.Name):$($m.LineNumber)"
        }
    }
}

# --- Invariant 4: verified stages don't apologize ---------------------------
foreach ($tier in @('Lemmas', 'Theorems')) {
    $verifiedFiles = @(Get-ChildItem -LiteralPath (Join-Path $leanRoot $tier) -Filter *.lean -Recurse) +
                     @(Get-Item -LiteralPath (Join-Path $leanRoot "$tier.lean"))
    foreach ($f in $verifiedFiles) {
        foreach ($m in @(Select-String -LiteralPath $f.FullName -Pattern '\bsorry\b')) {
            $violations += "$tier tier apologizes: $($f.Name):$($m.LineNumber)"
        }
        foreach ($m in @(Select-String -LiteralPath $f.FullName -Pattern '^\s*import\s+Enthymemata')) {
            $violations += "$tier tier imports Enthymemata: $($f.Name):$($m.LineNumber)"
        }
        if ($tier -eq 'Lemmas') {
            foreach ($m in @(Select-String -LiteralPath $f.FullName -Pattern '^\s*import\s+Theorems')) {
                $violations += "Lemmas tier imports upward from Theorems: $($f.Name):$($m.LineNumber)"
            }
        }
    }
}

# --- Invariant 5: enthymema ledger ------------------------------------------
$declPattern = '^\s*(noncomputable\s+)?(private\s+)?(theorem|lemma|def|abbrev|instance|structure|inductive|axiom|opaque)\b'
$ledger = foreach ($f in Get-ChildItem -LiteralPath (Join-Path $leanRoot 'Enthymemata') -Filter *.lean -Recurse) {
    $sorries = @(Select-String -LiteralPath $f.FullName -Pattern '\bsorry\b').Count
    $decls   = @(Select-String -LiteralPath $f.FullName -Pattern $declPattern).Count
    $state   = if ($decls -eq 0) { 'unstated' }
               elseif ($sorries -eq 0) { 'PROOF-CLOSED' }
               else { "apologizing($sorries)" }
    [pscustomobject]@{ File = $f.Name; Decls = $decls; Sorries = $sorries; State = $state }
}

$ledger | Format-Table -AutoSize | Out-String -Width 120 | Write-Host

foreach ($row in $ledger) {
    if ($row.State -eq 'PROOF-CLOSED') {
        $notices += "$($row.File) is proof-closed — review it for promotion to Lemmas/."
    }
    elseif ($row.State -eq 'unstated') {
        $notices += "$($row.File) has no statements yet — still a protolemma in costume."
    }
}

foreach ($n in $notices)    { Write-Host "notice: $n" -ForegroundColor Yellow }
foreach ($v in $violations) { Write-Host "VIOLATION: $v" -ForegroundColor Red }

if ($violations.Count -gt 0) {
    Write-Host 'meta-ci: FAIL' -ForegroundColor Red
    exit 1
}
if ($Validate -and $notices.Count -gt 0) {
    Write-Host 'meta-ci: FAIL (strict — ledger notices are errors under -Validate)' -ForegroundColor Red
    exit 1
}
Write-Host 'meta-ci: PASS — lemmas and theorems prove; enthymemata may apologize.' -ForegroundColor Green
exit 0
