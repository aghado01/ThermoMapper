# Verify a Lean slice from a clean archive of <Commit>'s lean/: the checks script, the module build,
# #print axioms, and meta-ci -NoBuild.
#
#   pwsh -NoProfile -File verify-lean.ps1 -Commit 4383916 -AxiomsFile axioms.lean
#
# The axioms file is `import <Module>` followed by one `#print axioms <name>` per new result.
# The snapshot's lean/.lake/packages is a junction to the real packages, so Mathlib is not rebuilt.
# The junction is removed with `cmd /c rmdir` BEFORE the snapshot is deleted, and the script checks
# that the real packages survived: deleting through the junction would take them with it.
# Receipts land in artifacts/verify-lean-<Commit>/receipts/ next to a copy of the axioms file.
# -NoRun prepares and cleans up only; -Keep leaves the snapshot (the junction is removed regardless).
param(
  [Parameter(Mandatory)][string]$Commit,
  [string]$Module = 'Enthymemata.ConformalValidity',
  [string]$Checks = 'Protolemmata\Conformal\conformal-checks.py',
  [string]$AxiomsFile,
  [string]$Repo = 'D:\aghado01\ThermoMapper',
  [string]$Python = 'C:/Users/azrie/PDenv/PyenvPython/pyenv-win/versions/3.12.6/python.exe',
  [switch]$NoRun,
  [switch]$Keep
)
$ErrorActionPreference = 'Continue'
[System.Diagnostics.Process]::GetCurrentProcess().PriorityClass = 'BelowNormal'
Set-Location -LiteralPath $Repo

$snap = Join-Path $Repo "artifacts\verify-lean-$Commit"
$lean = Join-Path $snap 'lean'
$out = Join-Path $snap 'receipts'
$realPackages = Join-Path $Repo 'lean\.lake\packages'
$junction = Join-Path $lean '.lake\packages'
$packagesBefore = (Get-ChildItem -LiteralPath $realPackages -Force).Count
if (Test-Path -LiteralPath $lean) { throw "$lean exists: remove its junction with cmd /c rmdir first, then the folder" }
New-Item -ItemType Directory -Force -Path $snap, $out | Out-Null

$tarball = Join-Path $snap 'lean.tar'
git archive --output $tarball $Commit lean
& (Join-Path $env:SystemRoot 'System32\tar.exe') -xf $tarball -C $snap
Remove-Item -LiteralPath $tarball
New-Item -ItemType Directory -Force -Path (Join-Path $lean '.lake') | Out-Null
cmd /c mklink /J "$junction" "$realPackages" | Out-Null
$summary = Join-Path $out 'summary.txt'
"commit $Commit; packages through the junction $((Get-ChildItem -LiteralPath $junction -Force).Count)" | Out-File $summary

try {
  if (-not $NoRun) {
    $sw = [Diagnostics.Stopwatch]::StartNew()
    & $Python -I -X utf8 (Join-Path $lean $Checks) *> (Join-Path $out 'checks.txt')
    "checks exit $LASTEXITCODE in $([int]$sw.Elapsed.TotalSeconds) s" | Out-File -Append $summary
    Push-Location -LiteralPath $lean
    try {
      $sw.Restart()
      & lake build $Module *> (Join-Path $out 'build.txt')
      "build exit $LASTEXITCODE in $([int]$sw.Elapsed.TotalSeconds) s" | Out-File -Append $summary
      if ($AxiomsFile) {
        Copy-Item -LiteralPath $AxiomsFile -Destination (Join-Path $snap 'axioms.lean')
        & lake env lean (Join-Path $snap 'axioms.lean') *> (Join-Path $out 'axioms.txt')
        "axioms exit $LASTEXITCODE" | Out-File -Append $summary
      }
    } finally { Pop-Location }
    & pwsh -NoProfile -File (Join-Path $lean 'scripts\meta-ci.ps1') -NoBuild *> (Join-Path $out 'meta-ci.txt')
    "meta-ci exit $LASTEXITCODE" | Out-File -Append $summary
  }
} finally {
  cmd /c rmdir "$junction"
  if (Test-Path -LiteralPath $junction) { throw "the junction is still at ${junction}; snapshot kept" }
  $packagesAfter = (Get-ChildItem -LiteralPath $realPackages -Force).Count
  if ($packagesAfter -ne $packagesBefore) { throw "real packages changed: $packagesBefore before, $packagesAfter after" }
  if (-not $Keep) { Remove-Item -LiteralPath $lean -Recurse -Force }
  Get-Content -LiteralPath $summary
}
