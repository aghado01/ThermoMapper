# Verify a C# slice from a clean archive of HEAD: the Commit selection over the slice's changed files.
#
#   git diff --name-only <last verified commit>..HEAD -- src tests projects > changed.txt
#   pwsh -NoProfile -File verify-csharp.ps1 -Name e8-s7d-<head> -ChangedFile changed.txt
#
# Receipts land in artifacts/verify-<Name>/ (commit.log, commit.trx, counters.txt). The snapshot and
# build in artifacts/verify-main/ are removed afterwards unless -Keep. -NoRun builds the snapshot only.
# Runs below normal priority on two processors: the machine is shared.
param(
  [Parameter(Mandatory)][string]$Name,
  [string]$ChangedFile,
  [string]$Repo = 'D:\aghado01\ThermoMapper',
  [switch]$NoRun,
  [switch]$Keep
)
$ErrorActionPreference = 'Continue'
[System.Diagnostics.Process]::GetCurrentProcess().PriorityClass = 'BelowNormal'
$env:MSBUILDDISABLENODEREUSE = '1'
$env:DOTNET_PROCESSOR_COUNT = '2'
Set-Location -LiteralPath $Repo

$work = Join-Path $Repo 'artifacts\verify-main'
$snap = Join-Path $work 'snap'
$out = Join-Path $Repo "artifacts\verify-$Name"
if (Test-Path -LiteralPath $work) { throw "$work exists: another verification is running or was left behind" }
New-Item -ItemType Directory -Force -Path $snap, $out | Out-Null

# The snapshot leaves out r/: the oracle tests find the real r/ by walking up from it.
$tarball = Join-Path $work 'head.tar'
git archive --output $tarball HEAD
& (Join-Path $env:SystemRoot 'System32\tar.exe') -xf $tarball -C $snap --exclude=r
Remove-Item -LiteralPath $tarball
"HEAD $(git rev-parse --short HEAD); snapshot files $((Get-ChildItem -LiteralPath $snap -Recurse -File).Count)" |
  Tee-Object -FilePath (Join-Path $out 'snapshot.txt')

if (-not $NoRun) {
  if (-not $ChangedFile) { throw '-ChangedFile is required unless -NoRun' }
  $changed = Get-Content -LiteralPath $ChangedFile | Where-Object { $_ -and -not $_.StartsWith('#') }
  & .\scripts\test-select.ps1 -SourceRoot $snap -Scope Commit -Changed $changed -ArtifactsRoot "$work\bin\" *> (Join-Path $out 'commit.log')
  "exit $LASTEXITCODE" | Out-File -Append (Join-Path $out 'commit.log')
  $trx = Get-ChildItem -LiteralPath $work -Recurse -Filter *.trx | Sort-Object LastWriteTime | Select-Object -Last 1
  if ($trx) {
    Copy-Item -LiteralPath $trx.FullName -Destination (Join-Path $out 'commit.trx')
    $c = ([xml](Get-Content -LiteralPath $trx.FullName -Raw)).TestRun.ResultSummary.Counters
    "total $($c.total) executed $($c.executed) passed $($c.passed) failed $($c.failed) notExecuted $($c.notExecuted)" |
      Tee-Object -FilePath (Join-Path $out 'counters.txt')
  } else { 'no TRX produced: read commit.log' | Tee-Object -FilePath (Join-Path $out 'counters.txt') }
}

if (-not $Keep) { Remove-Item -LiteralPath $work -Recurse -Force }
