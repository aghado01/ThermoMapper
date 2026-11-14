# Verify a C# slice from a clean archive of HEAD: the Commit selection over the slice's changed files.
#
#   git diff --name-only <last verified commit>..HEAD -- src tests projects > changed.txt
#   pwsh -NoProfile -File scripts/verify-csharp.ps1 -Name e8-s7d-<head> -ChangedFile changed.txt
#
# Receipts land in artifacts/verify-<Name>/ (commit.log, commit.trx, counters.txt). The snapshot and
# build in artifacts/verify-main/ are removed afterwards unless -Keep. -NoRun stops after extracting
# the snapshot; pair it with -Keep to look at the snapshot. A relative -ChangedFile is read from the
# directory the script is called from.
# The selection runs in the regression test project unless -Project names others
# (tests/<Name>.Tests/<Name>.Tests.csproj, comma-separated or as a list), each in turn on the one snapshot:
#   pwsh -NoProfile -File scripts/verify-csharp.ps1 -Name e8-s7e -ChangedFile changed.txt `
#     -Project tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj,tests/Maths.Geometry.Tests/Maths.Geometry.Tests.csproj
# Name every test project the changed sources reach: `test-select.ps1 -Project <p> -Changed <files> -ListOnly`
# shows each one's selection. With several, the logs and TRX files carry the project's short name and
# counters.txt has one line per project.
# Runs below normal priority on two processors: the machine is shared.
param(
    [Parameter(Mandatory)][string]$Name,
    [string]$ChangedFile,
    [string[]]$Project,
    [string]$RepoRoot = (Split-Path $PSScriptRoot),
    [switch]$NoRun,
    [switch]$Keep
)
$ErrorActionPreference = 'Continue'
# The arguments are checked before anything is written, so a bad call leaves no artifacts/verify-main.
if (-not $NoRun) {
    if (-not $ChangedFile) { throw '-ChangedFile is required unless -NoRun' }
    $ChangedFile = (Resolve-Path -LiteralPath $ChangedFile -ErrorAction Stop).Path
}
[System.Diagnostics.Process]::GetCurrentProcess().PriorityClass = 'BelowNormal'
$env:MSBUILDDISABLENODEREUSE = '1'
$env:DOTNET_PROCESSOR_COUNT = '2'
Set-Location -LiteralPath $RepoRoot

$work = Join-Path $RepoRoot 'artifacts\verify-main'
$snap = Join-Path $work 'snap'
$out = Join-Path $RepoRoot "artifacts\verify-$Name"
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
    $changed = Get-Content -LiteralPath $ChangedFile | Where-Object { $_ -and -not $_.StartsWith('#') }
    # pwsh -File hands a list over as one comma-joined string, so split it here.
    $projects = @($Project | ForEach-Object { $_ -split ',' } | ForEach-Object { $_.Trim() } | Where-Object { $_ })
    if ($projects.Count -eq 0) { $projects = @('') }   # test-select's default, the regression project
    $several = $projects.Count -gt 1
    $counters = Join-Path $out 'counters.txt'
    Remove-Item -LiteralPath $counters -ErrorAction SilentlyContinue
    foreach ($p in $projects) {
        $key = if ($p) { ((Split-Path (Split-Path $p) -Leaf) -replace '\.Tests$', '' -replace '^Maths\.', '').ToLowerInvariant() } else { 'regression' }
        $suffix = if ($several) { "-$key" } else { '' }
        $label = if ($several) { "{0,-14}" -f $key } else { '' }
        $log = Join-Path $out "commit$suffix.log"
        $projectArgs = if ($p) { @{ Project = $p } } else { @{} }
        $started = Get-Date
        & .\scripts\test-select.ps1 -SourceRoot $snap -Scope Commit -Changed $changed -ArtifactsRoot "$work\bin\" @projectArgs *> $log
        "exit $LASTEXITCODE" | Out-File -Append $log
        # Only this run's TRX: a run that produced none must not report the previous project's.
        $trx = Get-ChildItem -LiteralPath $work -Recurse -Filter *.trx | Where-Object LastWriteTime -GE $started |
            Sort-Object LastWriteTime | Select-Object -Last 1
        if ($trx) {
            Copy-Item -LiteralPath $trx.FullName -Destination (Join-Path $out "commit$suffix.trx")
            $c = ([xml](Get-Content -LiteralPath $trx.FullName -Raw)).TestRun.ResultSummary.Counters
            "${label}total $($c.total) executed $($c.executed) passed $($c.passed) failed $($c.failed) notExecuted $($c.notExecuted)" |
                Tee-Object -Append -FilePath $counters
        }
        else { "${label}no TRX produced: read commit$suffix.log" | Tee-Object -Append -FilePath $counters }
    }
}

if (-not $Keep) { Remove-Item -LiteralPath $work -Recurse -Force }
