[CmdletBinding()]
param(
    [Parameter(Mandatory, Position = 0)]
    [ValidateSet('profile', 'plan', 'apply')]
    [string] $Command,
    [Parameter(Position = 1)]
    [string] $PlanPath     = '',
    [string] $RepoRoot     = (Split-Path $PSScriptRoot),
    [string] $Rules        = '',
    [string] $ArchiveTo    = '',
    [string] $ArtifactsRoot = (Join-Path (Split-Path $PSScriptRoot) 'artifacts'),
    [switch] $Execute,
    [switch] $NoDisplay,
    [switch] $Rebuild
)

$ErrorActionPreference = 'Stop'

$project = Join-Path $RepoRoot 'projects\StorageAudit\StorageAudit.csproj'
$binDir  = Join-Path $ArtifactsRoot 'bin\StorageAudit\Debug\net10.0'
$exePath = Join-Path $binDir 'StorageAudit.exe'

if ($Rebuild -or -not (Test-Path -LiteralPath $exePath)) {
    Write-Host 'Building storage-audit...'
    dotnet build $project -v minimal --nologo "-p:ArtifactsRoot=$(Join-Path $ArtifactsRoot '')"
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}

if (-not (Test-Path -LiteralPath $exePath)) {
    Write-Error "StorageAudit.exe not found at $exePath after build."
    exit 1
}

$runArgs = @($Command)
$outDir  = ''
if ($Command -eq 'apply') {
    if (-not $PlanPath) {
        Write-Host 'apply needs the reviewed plan file: storage-audit.ps1 apply <plan.json> [-Execute] [-ArchiveTo <path>]'
        exit 2
    }
    $runArgs += (Resolve-Path -LiteralPath $PlanPath).Path
    $runArgs += '--repo'; $runArgs += $RepoRoot
    if ($Execute)   { $runArgs += '--execute' }
    if ($ArchiveTo) { $runArgs += '--archive-to'; $runArgs += $ArchiveTo }
}
else {
    $outDir = Join-Path $RepoRoot ("artifacts\storage-audit\runs\" + (Get-Date -Format 'yyyyMMdd_HHmmss'))
    $runArgs += '--repo'; $runArgs += $RepoRoot
    $runArgs += '--out';  $runArgs += $outDir
}
if ($Rules) { $runArgs += '--rules'; $runArgs += $Rules }

# Run from a copy under %TEMP%, so apply can delete artifacts\bin without deleting itself.
$shadow = Join-Path ([IO.Path]::GetTempPath()) "storage-audit-shadow\$PID"
New-Item -ItemType Directory -Path $shadow -Force | Out-Null
Get-ChildItem -LiteralPath $binDir | Copy-Item -Destination $shadow -Recurse -Force

try {
    & (Join-Path $shadow 'StorageAudit.exe') @runArgs
    $exitCode = $LASTEXITCODE
}
finally {
    Remove-Item -LiteralPath $shadow -Recurse -Force -ErrorAction SilentlyContinue
}

if ($exitCode -eq 0 -and $outDir -and -not $NoDisplay) {
    $report = Join-Path $outDir ($(if ($Command -eq 'plan') { 'plan.md' } else { 'profile.md' }))
    if (Test-Path -LiteralPath $report) {
        Write-Host ''
        try   { Show-Markdown -LiteralPath $report }
        catch { Get-Content -LiteralPath $report | Write-Host }
    }
}

exit $exitCode
