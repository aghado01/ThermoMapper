# Verify a slice whose files reach test projects other than the regression one: from a clean archive of HEAD,
# sweep every test project for the classes the changed files select, then run
#   - each project named in -Whole whole (scope Slice), and
#   - every other project with a selected class by its Commit selection (selected classes and the pinned tier).
# The regression project is left to scripts/verify-csharp.ps1's Commit selection.
#
#   pwsh -NoProfile -File verify-projects.ps1 -ChangedFile changed.txt -Tag defects-s4-projects -Whole Maths.LinAlg.Tests,TDA.Mapper.Tests
#
# Receipts land in artifacts/verify-<Tag>-<head>/: changed.txt, closure.txt, counters.txt, one log and trx per project.
# Uses artifacts/verify-main as the standard script does, and removes it at the end.
param(
    [Parameter(Mandatory)][string]$ChangedFile,
    [Parameter(Mandatory)][string]$Tag,
    [string[]]$Whole,
    [string[]]$Skip,
    [string]$RepoRoot = 'D:\aghado01\ThermoMapper'
)
$ErrorActionPreference = 'Continue'
$ChangedFile = (Resolve-Path -LiteralPath $ChangedFile -ErrorAction Stop).Path
# pwsh -File hands a comma list over as one string.
function Split-Names($names) { @($names | ForEach-Object { $_ -split ',' } | ForEach-Object { $_.Trim() } | Where-Object { $_ }) }
$Whole = Split-Names $Whole
$Skip = Split-Names $Skip
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
    $changed = @(Get-Content -LiteralPath $ChangedFile | Where-Object { $_ -and -not $_.StartsWith('#') })
    $changed | Set-Content -LiteralPath (Join-Path $out 'changed.txt')
    $select = Join-Path $RepoRoot 'scripts\test-select.ps1'

    $projects = Get-ChildItem -LiteralPath (Join-Path $snap 'tests') -Recurse -Filter *.Tests.csproj |
        ForEach-Object { [IO.Path]::GetRelativePath($snap, $_.FullName).Replace('\', '/') } | Sort-Object
    $runs = @()
    foreach ($p in $projects) {
        $tagName = [IO.Path]::GetFileNameWithoutExtension($p)
        try { $list = & $select -SourceRoot $snap -Project $p -Scope Commit -Changed $changed -ListOnly *>&1 | Out-String }
        catch { $list = "error: $_" }
        $m = [regex]::Match($list, '(\d+) in the compile closure')
        $sel = [regex]::Match($list, 'selected classes: (\d+) of (\d+)')
        $selected = if ($sel.Success) { [int]$sel.Groups[1].Value } else { -1 }
        $note = if ($m.Success) { "closure $($m.Groups[1].Value); selected $selected of $($sel.Groups[2].Value)" } else { ($list.Trim() -split "`n")[0] }
        $scope = if ($tagName -eq 'Maths.Regression.Tests') { 'left to verify-csharp.ps1' }
                 elseif ($Skip -contains $tagName) { 'skipped by -Skip' }
                 elseif ($Whole -contains $tagName) { 'Slice' }
                 elseif ($selected -gt 0) { 'Commit' }
                 else { 'not run' }
        "$p : $note; $scope" | Tee-Object -Append -FilePath (Join-Path $out 'closure.txt')
        if ($scope -eq 'Slice' -or $scope -eq 'Commit') { $runs += [pscustomobject]@{ Project = $p; Name = $tagName; Scope = $scope } }
    }

    # Whole projects first, then the selections.
    foreach ($r in ($runs | Sort-Object @{ Expression = { $_.Scope -ne 'Slice' } }, Name)) {
        $started = Get-Date
        $log = Join-Path $out "$($r.Name).log"
        try {
            & $select -SourceRoot $snap -Project $r.Project -Scope $r.Scope -Changed $changed -ArtifactsRoot "$work\bin\" *> $log
            "exit $LASTEXITCODE" | Out-File -Append $log
        }
        catch { "error: $_" | Out-File -Append $log }
        $trx = Get-ChildItem -LiteralPath $work -Recurse -Filter *.trx | Where-Object { $_.LastWriteTime -ge $started } |
            Sort-Object LastWriteTime | Select-Object -Last 1
        if ($trx) {
            Copy-Item -LiteralPath $trx.FullName -Destination (Join-Path $out "$($r.Name).trx")
            $x = [xml](Get-Content -LiteralPath $trx.FullName -Raw)
            $c = $x.TestRun.ResultSummary.Counters
            $skipped = @($x.TestRun.Results.UnitTestResult | Where-Object { $_.outcome -eq 'NotExecuted' }).Count
            "$($r.Name) ($($r.Scope)) total $($c.total) passed $($c.passed) failed $($c.failed) skipped $skipped" |
                Tee-Object -Append -FilePath (Join-Path $out 'counters.txt')
        }
        else { "$($r.Name) ($($r.Scope)) no TRX produced: read $($r.Name).log" | Tee-Object -Append -FilePath (Join-Path $out 'counters.txt') }
    }
}
finally {
    Remove-Item -LiteralPath $work -Recurse -Force
    "removed $work; out $out"
}
