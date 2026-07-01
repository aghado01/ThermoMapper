[CmdletBinding()]
param(
    [string] $RepoRoot      = (Split-Path $PSScriptRoot),
    [string] $OutputBaseDir = (Join-Path (Split-Path $PSScriptRoot) 'artifacts/repo-audit'),
    [switch] $Validate,
    [switch] $NoGit,
    [string] $Impact        = '',
    [switch] $NoDisplay,
    [switch] $Rebuild
)

$ErrorActionPreference = 'Stop'

$project = Join-Path $RepoRoot 'projects\RepoAudit\RepoAudit.csproj'
$exePath = Join-Path $RepoRoot 'artifacts\bin\RepoAudit\Debug\net10.0\RepoAudit.exe'

if ($Rebuild -or -not (Test-Path -LiteralPath $exePath)) {
    Write-Host 'Building repo-audit...'
    dotnet build $project -v minimal --nologo
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}

if (-not (Test-Path -LiteralPath $exePath)) {
    Write-Error "RepoAudit.exe not found at $exePath after build."
    exit 1
}

$stamp      = Get-Date -Format 'yyyyMMdd_HHmmss'
$attemptDir = Join-Path $OutputBaseDir $stamp
New-Item -ItemType Directory -Path $attemptDir -Force | Out-Null

$logPath = Join-Path $attemptDir 'run.log'

$runArgs = @($RepoRoot, '--attempt-dir', $attemptDir)
if ($Validate) { $runArgs += '--validate' }
if ($NoGit)    { $runArgs += '--no-git' }
if ($Impact)   { $runArgs += '--impact'; $runArgs += $Impact }

& $exePath @runArgs 2>&1 | Tee-Object -FilePath $logPath
$exitCode = $LASTEXITCODE

if ($exitCode -ne 0) {
    Write-Host ''
    Write-Host "repo-audit exited with code $exitCode. Log: $logPath"
    exit $exitCode
}

if (-not $Validate -and -not $Impact -and -not $NoDisplay) {
    $healthPath = Join-Path $attemptDir 'project-health.md'
    if (Test-Path -LiteralPath $healthPath) {
        Write-Host ''
        try   { Show-Markdown -LiteralPath $healthPath }
        catch { Get-Content -LiteralPath $healthPath | Write-Host }
    }
}

exit $exitCode
