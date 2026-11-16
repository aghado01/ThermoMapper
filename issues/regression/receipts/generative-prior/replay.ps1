# From any checkout containing this launcher. Each label gets a fresh snapshot and result directory.
param(
    [Parameter(Mandatory)][string]$Revision,
    [Parameter(Mandatory)][ValidatePattern('^[a-zA-Z0-9_-]+$')][string]$RunName
)
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $true
$repo = (git -C $PSScriptRoot rev-parse --show-toplevel).Trim()
$pin = (git -C $repo rev-parse --verify "$Revision^{commit}").Trim()
$receipts = 'issues/regression/receipts/generative-prior'
$run = Join-Path $repo "artifacts/generative-prior/$RunName"
if (Test-Path -LiteralPath $run) { throw "Run directory exists: $run. Use a fresh label." }
$source = Join-Path $run 'source'
$build = Join-Path $run 'build'
$results = Join-Path $run 'results'
$null = New-Item -ItemType Directory -Path $source
$archive = Join-Path $run 'source.tar'
git -C $repo archive --format=tar "--output=$archive" $pin -- Directory.Build.props projects/Maths.Rng projects/Maths.Distributions projects/Maths.LinAlg src/maths/rng src/maths/distributions src/maths/linalg "$receipts/probe" "$receipts/settings"
tar -xf $archive -C $source
Remove-Item -LiteralPath $archive

$env:DOTNET_PROCESSOR_COUNT = '1'
$env:DOTNET_CLI_TELEMETRY_OPTOUT = '1'
[Diagnostics.Process]::GetCurrentProcess().PriorityClass = 'BelowNormal'
$probe = Join-Path $source $receipts
dotnet build (Join-Path $probe 'probe/GenerativePrior.csproj') -c Release "-p:ArtifactsRoot=$build/" --disable-build-servers -m:1 -p:UseSharedCompilation=false
$assembly = Join-Path $build 'bin/GenerativePrior/Release/net10.0/GenerativePrior.dll'
dotnet $assembly (Join-Path $probe 'settings/g1-k10.json') $results $pin
