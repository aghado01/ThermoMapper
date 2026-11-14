# Verify a slice whose files lie outside the regression assembly's compile closure: from a clean archive of HEAD,
# find every test project whose closure holds a changed file, and run each such project whole (scope Slice).
# The regression assembly is left to scripts/verify-csharp.ps1's Commit selection.
#
#   pwsh -NoProfile -File verify-projects.ps1 -Commit 37eb63c -Tag defects-s2-projects
#
# Receipts land in artifacts/verify-<Tag>-<head>/: changed.txt, closure.txt, counters.txt, one log and trx per project.
# Uses artifacts/verify-main as the standard script does, and removes it at the end.
param(
    [Parameter(Mandatory)][string]$Commit,
    [Parameter(Mandatory)][string]$Tag,
    # Project names to run whole. Given, the closure sweep is skipped (closure.txt from an earlier sweep stands).
    [string[]]$Only,
    [string]$RepoRoot = 'D:\aghado01\ThermoMapper'
)
$ErrorActionPreference = 'Continue'
[System.Diagnostics.Process]::GetCurrentProcess().PriorityClass = 'BelowNormal'
$env:MSBUILDDISABLENODEREUSE = '1'
$env:DOTNET_PROCESSOR_COUNT = '2'
Set-Location -LiteralPath $RepoRoot

$head = (git rev-parse --short HEAD)
$name = "$Tag-$head"
$work = Join-Path $RepoRoot 'artifacts\verify-main'
$snap = Join-Path $work 'snap'
$out = Join-Path $RepoRoot "artifacts\verify-$name"

& pwsh -NoProfile -File (Join-Path $RepoRoot 'scripts\verify-csharp.ps1') -Name $name -NoRun -Keep
if (-not (Test-Path -LiteralPath $snap)) { throw "no snapshot at $snap" }

try {
    $changed = @(git diff --name-only "$Commit~1" $Commit -- src tests projects)
    $changed | Set-Content -LiteralPath (Join-Path $out 'changed.txt')

    $projects = Get-ChildItem -LiteralPath (Join-Path $snap 'tests') -Recurse -Filter *.Tests.csproj |
        ForEach-Object { [IO.Path]::GetRelativePath($snap, $_.FullName).Replace('\', '/') } | Sort-Object
    $reach = @()
    if ($Only) {
        # pwsh -File hands a comma list over as one string.
        $Only = @($Only | ForEach-Object { $_ -split ',' } | ForEach-Object { $_.Trim() } | Where-Object { $_ })
        $reach =@($projects | Where-Object { $Only -contains [IO.Path]::GetFileNameWithoutExtension($_) })
        $projects = @()
    }
    foreach ($p in $projects) {
        try {
            $list = & (Join-Path $RepoRoot 'scripts\test-select.ps1') -SourceRoot $snap -Project $p -Scope Commit -Changed $changed -ListOnly *>&1 | Out-String
        }
        catch { $list = "error: $_" }
        $m = [regex]::Match($list, '(\d+) in the compile closure')
        $sel = [regex]::Match($list, 'selected classes: (\d+) of (\d+)')
        $inClosure = if ($m.Success) { [int]$m.Groups[1].Value } else { -1 }
        $note = if ($m.Success) { "closure $inClosure; selected $($sel.Groups[1].Value) of $($sel.Groups[2].Value)" } else { ($list.Trim() -split "`n")[0] }
        "$p : $note" | Tee-Object -Append -FilePath (Join-Path $out 'closure.txt')
        if ($inClosure -gt 0) { $reach += $p }
    }

    foreach ($p in $reach) {
        $tagName = [IO.Path]::GetFileNameWithoutExtension($p)
        if ($tagName -eq 'Maths.Regression.Tests') { "$tagName left to verify-csharp.ps1's Commit selection" | Tee-Object -Append -FilePath (Join-Path $out 'counters.txt'); continue }
        $started = Get-Date
        $log = Join-Path $out "$tagName.log"
        try {
            & (Join-Path $RepoRoot 'scripts\test-select.ps1') -SourceRoot $snap -Project $p -Scope Slice -Changed $changed -ArtifactsRoot "$work\bin\" *> $log
            "exit $LASTEXITCODE" | Out-File -Append $log
        }
        catch { "error: $_" | Out-File -Append $log }
        $trx = Get-ChildItem -LiteralPath $work -Recurse -Filter *.trx | Where-Object { $_.LastWriteTime -ge $started } |
            Sort-Object LastWriteTime | Select-Object -Last 1
        if ($trx) {
            Copy-Item -LiteralPath $trx.FullName -Destination (Join-Path $out "$tagName.trx")
            $c = ([xml](Get-Content -LiteralPath $trx.FullName -Raw)).TestRun.ResultSummary.Counters
            "$tagName total $($c.total) executed $($c.executed) passed $($c.passed) failed $($c.failed) notExecuted $($c.notExecuted)" |
                Tee-Object -Append -FilePath (Join-Path $out 'counters.txt')
        }
        else { "$tagName no TRX produced: read $tagName.log" | Tee-Object -Append -FilePath (Join-Path $out 'counters.txt') }
    }
}
finally {
    Remove-Item -LiteralPath $work -Recurse -Force
    "removed $work; out $out"
}
