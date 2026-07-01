<#
.SYNOPSIS
meta-ci.ps1 — taxonomy gate for the Lean rigor harness.

.DESCRIPTION
Tiers: prelemmas/ (markdown briefs) -> Enthymemes/ (statements compile,
proofs still apologize) -> Lemmas/ (no apologies).

Invariants enforced:
  1. `lake build` is green — BOTH tiers must compile; an enthymeme owes
     proofs, never statements.
  2. The Lemmas tier never apologizes: no `sorry` token anywhere under
     Lemmas/ (even in prose — lemmas don't say the word), and no
     `import Enthymemes.*` (an import would let a "proved" lemma lean on a
     sorried declaration without the token ever appearing in its file).
  3. Enthymeme ledger, per file:
       unstated         — no declarations yet (a prelemma in costume)
       apologizing(n)   — n sorries outstanding
       PROMOTION-READY  — declarations present, zero sorries: move it to
                          Lemmas/ (file move; declaration names are stable —
                          tier is location, not identity)

Exit 1 on any invariant-2 violation or red build. With -Validate, ledger
notices (unstated / promotion-ready) also fail — taxonomy drift is an error
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
    Push-Location -LiteralPath $leanRoot
    try {
        lake build
        if ($LASTEXITCODE -ne 0) {
            Write-Host 'meta-ci: FAIL — build is red.' -ForegroundColor Red
            exit 1
        }
    }
    finally { Pop-Location }
}

$violations = @()
$notices    = @()

# --- Invariant 2: lemmas don't apologize -----------------------------------
$lemmaFiles = @(Get-ChildItem -LiteralPath (Join-Path $leanRoot 'Lemmas') -Filter *.lean -Recurse) +
              @(Get-Item -LiteralPath (Join-Path $leanRoot 'Lemmas.lean'))
foreach ($f in $lemmaFiles) {
    foreach ($m in @(Select-String -LiteralPath $f.FullName -Pattern '\bsorry\b')) {
        $violations += "Lemmas tier apologizes: $($f.Name):$($m.LineNumber)"
    }
    foreach ($m in @(Select-String -LiteralPath $f.FullName -Pattern '^\s*import\s+Enthymemes')) {
        $violations += "Lemmas tier imports Enthymemes: $($f.Name):$($m.LineNumber)"
    }
}

# --- Invariant 3: enthymeme ledger ------------------------------------------
$declPattern = '^\s*(noncomputable\s+)?(private\s+)?(theorem|lemma|def|abbrev|instance|structure|inductive)\b'
$ledger = foreach ($f in Get-ChildItem -LiteralPath (Join-Path $leanRoot 'Enthymemes') -Filter *.lean -Recurse) {
    $sorries = @(Select-String -LiteralPath $f.FullName -Pattern '\bsorry\b').Count
    $decls   = @(Select-String -LiteralPath $f.FullName -Pattern $declPattern).Count
    $state   = if ($decls -eq 0) { 'unstated' }
               elseif ($sorries -eq 0) { 'PROMOTION-READY' }
               else { "apologizing($sorries)" }
    [pscustomobject]@{ File = $f.Name; Decls = $decls; Sorries = $sorries; State = $state }
}

$ledger | Format-Table -AutoSize | Out-String -Width 120 | Write-Host

foreach ($row in $ledger) {
    if ($row.State -eq 'PROMOTION-READY') {
        $notices += "$($row.File) stopped apologizing — promote it to Lemmas/."
    }
    elseif ($row.State -eq 'unstated') {
        $notices += "$($row.File) has no statements yet — still a prelemma in costume."
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
Write-Host 'meta-ci: PASS — lemmas prove, enthymemes apologize.' -ForegroundColor Green
exit 0
