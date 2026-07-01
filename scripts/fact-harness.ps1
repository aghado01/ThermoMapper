<#
.SYNOPSIS
    Runs a test project through the parallel fact harness (projects/TestHarness.Runner).

.DESCRIPTION
    Convenience wrapper that runs each [Fact] of a test project as its own
    `dotnet test --no-build` process, fanned out across workers. Defaults to the
    VizCore.Tests suite in Release.

    Release is the default configuration on purpose: the harness is a performance
    runner, and unoptimized Debug builds push heavy compute facts (e.g. the 1500-point
    dense Fiedler eigensolve in CrescentFixtureGoldenDrafts) past their fixture timeouts.

.EXAMPLE
    scripts/fact-harness.ps1
    scripts/fact-harness.ps1 -Project tests/Synthetic.Tests/Synthetic.Tests.csproj
    scripts/fact-harness.ps1 -ListOnly
    scripts/fact-harness.ps1 -MaxWorkers 8
#>
[CmdletBinding()]
param(
    [string] $RepoRoot      = (Split-Path $PSScriptRoot),
    [string] $Project       = 'tests/VizCore.Tests/VizCore.Tests.csproj',
    [ValidateSet('Release', 'Debug')]
    [string] $Configuration = 'Release',
    [int]    $MaxWorkers    = 0,
    [switch] $ListOnly
)

$ErrorActionPreference = 'Stop'

$runner = Join-Path $RepoRoot 'projects/TestHarness.Runner/TestHarness.Runner.csproj'
if (-not (Test-Path -LiteralPath $runner)) {
    Write-Error "TestHarness.Runner not found at $runner."
    exit 1
}

$runnerArgs = @('--project', $Project, '--configuration', $Configuration)
if ($MaxWorkers -gt 0) { $runnerArgs += @('--max-workers', "$MaxWorkers") }
if ($ListOnly)         { $runnerArgs += '--list-only' }

dotnet run --project $runner -c Release -- @runnerArgs
exit $LASTEXITCODE
