. (Join-Path $PSScriptRoot 'scripts/r-session.ps1')
Initialize-ROracleSession -SetAliases | Out-Null
Set-Location $PSScriptRoot
Invoke-Rscript scripts/bootstrap.R
