param(
    [Parameter(Mandatory)][string]$Revision,
    [Parameter(Mandatory)][ValidateSet('g3','g4','g2')][string]$Experiment,
    [Parameter(Mandatory)][ValidatePattern('^[a-zA-Z0-9_-]+$')][string]$RunName
)
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $true
$repo = (git -C $PSScriptRoot rev-parse --show-toplevel).Trim()
$pin = (git -C $repo rev-parse --verify "$Revision^{commit}").Trim()
$receipts = 'issues/regression/receipts/generative-prior'
$run = Join-Path $repo "artifacts/generative-prior/$RunName"
if (Test-Path -LiteralPath $run) { throw "Use a fresh run label. Directory exists at $run" }
$source = Join-Path $run 'source'
$build = Join-Path $run 'build'
$null = New-Item -ItemType Directory -Path $source
$archive = Join-Path $run 'source.tar'
git -C $repo archive --format=tar "--output=$archive" $pin -- Directory.Build.props projects/Maths.Rng projects/Maths.Distributions projects/Maths.LinAlg projects/Maths.Samplers projects/Maths.Regression src/maths/rng src/maths/distributions src/maths/linalg src/maths/samplers src/maths/regression "$receipts/probe" "$receipts/lab" "$receipts/settings" issues/regression/receipts/latent-recovery/probe/Synthesizer.cs
tar -xf $archive -C $source
Remove-Item -LiteralPath $archive
$env:DOTNET_PROCESSOR_COUNT = '1'
$env:DOTNET_CLI_TELEMETRY_OPTOUT = '1'
[Diagnostics.Process]::GetCurrentProcess().PriorityClass = 'BelowNormal'
$probe = Join-Path $source $receipts
dotnet build (Join-Path $probe 'lab/GenerativePriorLab.csproj') -c Release "-p:ArtifactsRoot=$build/" --disable-build-servers -m:1 -p:UseSharedCompilation=false
$settingsName = switch ($Experiment) {
    g3 { 'g3-weights.json' }
    g4 { 'g4-audio.json' }
    g2 { 'g2-calibration.json' }
}
$settings = Join-Path $probe "settings/$settingsName"
dotnet (Join-Path $build 'bin/GenerativePriorLab/Release/net10.0/GenerativePriorLab.dll') $Experiment $settings (Join-Path $run 'results') $pin
