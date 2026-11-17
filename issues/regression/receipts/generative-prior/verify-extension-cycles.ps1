param(
    [string]$Revision = '0fe9045',
    [Parameter(Mandatory)][ValidatePattern('^[a-zA-Z0-9_-]+$')][string]$RunName,
    [string]$Analysis = 'artifacts/generative-prior/g1-extensions-primary/analysis/analysis.json'
)
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $true
$repo = (git -C $PSScriptRoot rev-parse --show-toplevel).Trim()
& (Join-Path $PSScriptRoot 'replay.ps1') -Revision $Revision -RunName $RunName -Experiment extensions -Cells d2-mag-k10 -TrialIndices 6 -FullIterations
$run = Join-Path $repo "artifacts/generative-prior/$RunName"
$pin = (git -C $repo rev-parse --verify "$Revision^{commit}").Trim()
$settingsPath = Join-Path $run 'source/issues/regression/receipts/generative-prior/settings/g1-extensions.json'
$config = Get-Content -LiteralPath $settingsPath -Raw | ConvertFrom-Json -AsHashtable
$audit = Get-Content -LiteralPath (Join-Path $repo $Analysis) -Raw | ConvertFrom-Json -AsHashtable
$assembly = Join-Path $run 'build/bin/GenerativePrior/Release/net10.0/GenerativePrior.dll'
$settingsDirectory = Join-Path $run 'control-settings'
$null = New-Item -ItemType Directory -Path $settingsDirectory
$env:DOTNET_PROCESSOR_COUNT = '1'
[Diagnostics.Process]::GetCurrentProcess().PriorityClass = 'BelowNormal'
$plan = @()
foreach ($key in ($audit.cycle_representatives.Keys | Sort-Object { [int]$_ })) {
    $case = $audit.cycle_representatives[$key]
    if ([int]$key -eq 96) {
        $plan += @{ Period = 96; Cell = $case.cell; Sign = $case.sign; Repetition = $case.repetition
            Output = 'results/d2-mag-k10-on/trials.jsonl' }
        continue
    }
    $cell = $config.ExtendedCells | Where-Object { $_.Name -eq $case.cell }
    $solver = $config.Solver | ConvertTo-Json -Depth 10 | ConvertFrom-Json -AsHashtable
    $solver.CodeDimension = $cell.K
    $solver.MasterSeed = $cell.MasterSeed
    $solver.IntermediateWidth = if ($cell.Depth -eq 3) { 400 } else { $null }
    $solver.Magnitude = $cell.Magnitude
    $solver.SignCheck = $case.sign
    $solver.FastForward = $false
    $solver.Step = [Math]::Pow(2, $cell.Depth) / ($cell.Depth * $cell.Depth)
    $setting = Join-Path $settingsDirectory "period-$key.json"
    $solver | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $setting
    $result = Join-Path $run "results/period-$key"
    dotnet $assembly $setting $result $pin ([string]$case.repetition)
    $plan += @{ Period = [int]$key; Cell = $case.cell; Sign = $case.sign; Repetition = $case.repetition
        Output = "results/period-$key/trials.jsonl"; Settings = "control-settings/period-$key.json" }
}
$plan | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Join-Path $run 'control-plan.json')
