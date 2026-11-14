# Run named qualification tests, the ones gated behind BASS_QUALIFICATION, from a clean archive of HEAD. test-select
# takes no filter, so a slice's gated tests run here, beside verify-csharp.ps1's Commit selection.
#
#   pwsh -NoProfile -File scripts/verify-qualification.ps1 -Name e21-s6-qual `
#     -Filter 'FullyQualifiedName~RandomEffectsNoiseQualificationTests.Geweke_TheSweepLeavesItsJointInvariant'
#
# -Filter is dotnet test's filter: join several tests with |, and narrow a theory's rows with &DisplayName~<text>.
# The snapshot is verify-csharp.ps1 -NoRun -Keep's. Receipts land in artifacts/verify-<Name>/: qual.log, qual.trx,
# counters.txt and outputs.txt, each test's outcome with its written output, which is where measurements report.
# The snapshot and build in artifacts/verify-main/ are removed afterwards unless -Keep.
# Runs below normal priority on two processors. VSTEST_CONNECTION_TIMEOUT is 600 s: under load the test host can take
# longer than the default to connect before a long theory.
param(
  [Parameter(Mandatory)][string]$Name,
  [Parameter(Mandatory)][string]$Filter,
  [string]$Project = 'tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj',
  [string]$RepoRoot = (Split-Path $PSScriptRoot),
  [switch]$Keep
)
$ErrorActionPreference = 'Continue'
[System.Diagnostics.Process]::GetCurrentProcess().PriorityClass = 'BelowNormal'
$env:MSBUILDDISABLENODEREUSE = '1'
$env:DOTNET_PROCESSOR_COUNT = '2'
$env:BASS_QUALIFICATION = '1'
$env:VSTEST_CONNECTION_TIMEOUT = '600'

$work = Join-Path $RepoRoot 'artifacts\verify-main'
$out = Join-Path $RepoRoot "artifacts\verify-$Name"
& pwsh -NoProfile -File (Join-Path $PSScriptRoot 'verify-csharp.ps1') -Name $Name -NoRun -Keep -RepoRoot $RepoRoot
if (-not (Test-Path -LiteralPath (Join-Path $work 'snap'))) { throw "no snapshot in ${work}; read verify-csharp.ps1's error above" }

try {
  $projectPath = Join-Path (Join-Path $work 'snap') $Project
  $results = Join-Path $out 'results'
  $sw = [Diagnostics.Stopwatch]::StartNew()
  & dotnet test $projectPath -c Release "-p:ArtifactsRoot=$work\bin\" --disable-build-servers --filter $Filter `
    --logger 'trx;LogFileName=qual.trx' --results-directory $results *> (Join-Path $out 'qual.log')
  "exit $LASTEXITCODE in $([int]$sw.Elapsed.TotalSeconds) s" | Tee-Object -Append -FilePath (Join-Path $out 'qual.log')

  $trx = Get-ChildItem -LiteralPath $results -Recurse -Filter *.trx -ErrorAction SilentlyContinue | Sort-Object LastWriteTime | Select-Object -Last 1
  if ($trx) {
    Copy-Item -LiteralPath $trx.FullName -Destination (Join-Path $out 'qual.trx')
    [xml]$x = Get-Content -LiteralPath $trx.FullName -Raw
    $c = $x.TestRun.ResultSummary.Counters
    "total $($c.total) executed $($c.executed) passed $($c.passed) failed $($c.failed) notExecuted $($c.notExecuted)" |
      Tee-Object -FilePath (Join-Path $out 'counters.txt')
    $lines = foreach ($r in $x.TestRun.Results.UnitTestResult) {
      "$($r.outcome) $($r.duration) $($r.testName)"
      if ($r.Output.StdOut) { $r.Output.StdOut -split "`n" | ForEach-Object { "    $_".TrimEnd() } }
      if ($r.Output.ErrorInfo.Message) { "    ERROR: $($r.Output.ErrorInfo.Message)" }
    }
    $lines | Set-Content -LiteralPath (Join-Path $out 'outputs.txt')
  }
  else { 'no TRX produced: read qual.log' | Tee-Object -FilePath (Join-Path $out 'counters.txt') }
}
finally {
  if (-not $Keep) { Remove-Item -LiteralPath $work -Recurse -Force -ErrorAction SilentlyContinue }
}
