# Project-local R session helpers for the validation oracles.
#
# This mirrors the personal PDenv aliases, but resolves the R toolchain here so
# repo scripts do not require a pre-sourced shell profile.

$script:RExe = $null
$script:RscriptExe = $null
$script:ROracleToolchain = $null

function Get-ROracleEnvValue {
    param([Parameter(Mandatory)][string]$Name)

    $value = [Environment]::GetEnvironmentVariable($Name, 'Process')
    if (-not [string]::IsNullOrWhiteSpace($value)) {
        return [Environment]::ExpandEnvironmentVariables($value)
    }

    $value = [Environment]::GetEnvironmentVariable($Name, 'User')
    if (-not [string]::IsNullOrWhiteSpace($value)) {
        return [Environment]::ExpandEnvironmentVariables($value)
    }

    return $null
}

function Get-ROracleRscriptFromHome {
    param([Parameter(Mandatory)][string]$RHome)

    foreach ($rel in @('bin\x64\Rscript.exe', 'bin\Rscript.exe')) {
        $candidate = Join-Path $RHome $rel
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
            return $candidate
        }
    }

    return $null
}

function Get-ROracleRExeFromHome {
    param([Parameter(Mandatory)][string]$RHome)

    foreach ($rel in @('bin\x64\R.exe', 'bin\R.exe')) {
        $candidate = Join-Path $RHome $rel
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
            return $candidate
        }
    }

    return $null
}

function Get-ROracleRHomeFromRscript {
    param([Parameter(Mandatory)][string]$Rscript)

    $dir = Split-Path -Parent $Rscript
    if ((Split-Path -Leaf $dir) -ieq 'x64') {
        $bin = Split-Path -Parent $dir
        if ((Split-Path -Leaf $bin) -ieq 'bin') {
            return Split-Path -Parent $bin
        }
    }

    if ((Split-Path -Leaf $dir) -ieq 'bin') {
        return Split-Path -Parent $dir
    }

    return $null
}

function Get-ROracleProfileRoots {
    $roots = @()

    foreach ($name in @('USERPROFILE', 'HOME')) {
        $value = Get-ROracleEnvValue $name
        if ($value) { $roots += $value }
    }

    $profile = [Environment]::GetFolderPath('UserProfile')
    if ($profile) { $roots += $profile }

    foreach ($name in @('CLAUDE_CONFIG_DIR', 'CODEX_HOME')) {
        $value = Get-ROracleEnvValue $name
        if ($value) { $roots += (Split-Path -Parent $value) }
    }

    $roots | Where-Object { $_ } | Select-Object -Unique
}

function Get-ROraclePdenvRoots {
    $roots = @()

    $portableRoot = Get-ROracleEnvValue 'PORTABLE_ROOT'
    if ($portableRoot) { $roots += $portableRoot }

    foreach ($name in @('CLAUDE_CODE_SHELL', 'CLAUDE_CODE_GIT_BASH_PATH')) {
        $value = Get-ROracleEnvValue $name
        if (-not $value) { continue }

        $normalized = $value -replace '/', '\'
        $match = [regex]::Match($normalized, '^(?<root>.*?\\PDenv)(\\|$)')
        if ($match.Success) {
            $roots += $match.Groups['root'].Value
        }
    }

    foreach ($profile in Get-ROracleProfileRoots) {
        $roots += (Join-Path $profile 'PDenv')
    }

    $roots | Where-Object { $_ } | Select-Object -Unique
}

function Get-ROracleCandidateHomes {
    $rHome = Get-ROracleEnvValue 'R_HOME'
    if ($rHome) { $rHome }

    $rlangRoots = @()
    foreach ($root in Get-ROraclePdenvRoots) {
        $rlangRoots += (Join-Path $root 'rlang')
    }

    foreach ($root in ($rlangRoots | Select-Object -Unique)) {
        if (Test-Path -LiteralPath $root -PathType Container) {
            Get-ChildItem -LiteralPath $root -Directory -Filter 'R-*' |
                Sort-Object Name -Descending |
                ForEach-Object { $_.FullName }
        }
    }

    $path = @(
        [Environment]::GetEnvironmentVariable('PATH', 'Process'),
        [Environment]::GetEnvironmentVariable('PATH', 'User')
    ) -join [IO.Path]::PathSeparator

    foreach ($entry in ($path -split [IO.Path]::PathSeparator)) {
        if ([string]::IsNullOrWhiteSpace($entry)) { continue }
        $candidate = Join-Path ([Environment]::ExpandEnvironmentVariables($entry)) 'Rscript.exe'
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
            $rHomeFromPath = Get-ROracleRHomeFromRscript $candidate
            if ($rHomeFromPath) { $rHomeFromPath }
        }
    }
}

function Get-ROracleToolchain {
    foreach ($rHome in (Get-ROracleCandidateHomes | Select-Object -Unique)) {
        $rscript = Get-ROracleRscriptFromHome $rHome
        if (-not $rscript) { continue }

        $rlangRoot = Split-Path -Parent $rHome
        $rLibs = Get-ROracleEnvValue 'R_LIBS'
        if (-not $rLibs -and (Test-Path -LiteralPath (Join-Path $rlangRoot 'library') -PathType Container)) {
            $rLibs = Join-Path $rlangRoot 'library'
        }

        $renvRoot = Get-ROracleEnvValue 'RENV_PATHS_ROOT'
        if (-not $renvRoot -and (Test-Path -LiteralPath (Join-Path $rlangRoot 'renv') -PathType Container)) {
            $renvRoot = Join-Path $rlangRoot 'renv'
        }

        $renvCache = Get-ROracleEnvValue 'RENV_PATHS_CACHE'
        if (-not $renvCache -and $renvRoot) {
            $renvCache = Join-Path $renvRoot 'cache'
        }

        return [pscustomobject]@{
            RHome      = $rHome
            RExe       = Get-ROracleRExeFromHome $rHome
            RscriptExe = $rscript
            RLibs      = $rLibs
            RenvRoot   = $renvRoot
            RenvCache  = $renvCache
        }
    }

    return $null
}

function Initialize-ROracleSession {
    param([switch]$SetAliases)

    $toolchain = Get-ROracleToolchain
    if (-not $toolchain) {
        return $null
    }

    $env:R_HOME = $toolchain.RHome
    if ($toolchain.RLibs) { $env:R_LIBS = $toolchain.RLibs }
    if ($toolchain.RenvRoot) { $env:RENV_PATHS_ROOT = $toolchain.RenvRoot }
    if ($toolchain.RenvCache) { $env:RENV_PATHS_CACHE = $toolchain.RenvCache }
    if (-not (Get-ROracleEnvValue 'RENV_CONFIG_SANDBOX_ENABLED')) {
        $env:RENV_CONFIG_SANDBOX_ENABLED = 'FALSE'
    }

    $rBin = Split-Path -Parent $toolchain.RscriptExe
    $pathEntries = $env:PATH -split [IO.Path]::PathSeparator
    if ($pathEntries -notcontains $rBin) {
        $env:PATH = $rBin + [IO.Path]::PathSeparator + $env:PATH
    }

    $script:ROracleToolchain = $toolchain
    $script:RExe = $toolchain.RExe
    $script:RscriptExe = $toolchain.RscriptExe

    if ($SetAliases) {
        Set-ROracleAliases
    }

    return $toolchain
}

function Assert-ROracleSession {
    if (-not $script:ROracleToolchain) {
        Initialize-ROracleSession | Out-Null
    }

    if (-not $script:RscriptExe) {
        throw 'Rscript.exe was not found. Set R_HOME, put Rscript on PATH, or install PDenv R under ~/PDenv/rlang.'
    }
}

function Invoke-Rscript {
    Assert-ROracleSession
    & $script:RscriptExe @args
}

function Invoke-RscriptTimed {
    param(
        [string[]]$ArgumentList = @(),
        [int]$TimeoutSeconds = 120
    )

    Assert-ROracleSession

    $psi = [Diagnostics.ProcessStartInfo]::new()
    $psi.FileName = $script:RscriptExe
    $psi.WorkingDirectory = (Get-Location).Path
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    foreach ($arg in $ArgumentList) {
        [void]$psi.ArgumentList.Add($arg)
    }

    $proc = [Diagnostics.Process]::Start($psi)
    $stdoutTask = $proc.StandardOutput.ReadToEndAsync()
    $stderrTask = $proc.StandardError.ReadToEndAsync()
    $waitMs = [Math]::Min([int64]$TimeoutSeconds * 1000, [int]::MaxValue)
    if (-not $proc.WaitForExit([int]$waitMs)) {
        try { $proc.Kill($true) } catch { }
        try { $proc.WaitForExit() } catch { }
        return [pscustomobject]@{
            ExitCode = $null
            TimedOut = $true
            Stdout   = $stdoutTask.GetAwaiter().GetResult()
            Stderr   = $stderrTask.GetAwaiter().GetResult()
        }
    }

    [pscustomobject]@{
        ExitCode = $proc.ExitCode
        TimedOut = $false
        Stdout   = $stdoutTask.GetAwaiter().GetResult()
        Stderr   = $stderrTask.GetAwaiter().GetResult()
    }
}

function Invoke-RscriptChecked {
    param(
        [string[]]$ArgumentList = @(),
        [int]$TimeoutSeconds = 120
    )

    $result = Invoke-RscriptTimed -ArgumentList $ArgumentList -TimeoutSeconds $TimeoutSeconds
    if ($result.Stdout) { Write-Host -NoNewline $result.Stdout }
    if ($result.Stderr) { Write-Host -NoNewline $result.Stderr }

    if ($result.TimedOut) {
        throw "Rscript timed out after $TimeoutSeconds seconds: $($ArgumentList -join ' ')"
    }

    if ($result.ExitCode -ne 0) {
        throw "Rscript failed with exit code $($result.ExitCode): $($ArgumentList -join ' ')"
    }

    return $result
}

function Invoke-R {
    Assert-ROracleSession
    if (-not $script:RExe) {
        throw 'R.exe was not found next to the resolved Rscript.exe.'
    }
    & $script:RExe @args
}

function Invoke-RVersion {
    Invoke-Rscript -e "cat(R.version.string, '\n')"
}

function Invoke-RenvStatus {
    Invoke-Rscript -e "renv::status()"
}

function Invoke-RenvRestore {
    Invoke-Rscript -e "renv::restore(prompt = FALSE)"
}

function Invoke-RenvSnapshot {
    Invoke-Rscript -e "renv::snapshot(prompt = FALSE)"
}

function Get-RAliases {
    @{
        Rscript = 'Invoke-Rscript'
        rs      = 'Invoke-Rscript'
        rx      = 'Invoke-R'
        rver    = 'Invoke-RVersion'
        rns     = 'Invoke-RenvStatus'
        rnr     = 'Invoke-RenvRestore'
        rnsnap  = 'Invoke-RenvSnapshot'
    }
}

function Set-ROracleAliases {
    foreach ($alias in (Get-RAliases).GetEnumerator()) {
        Set-Alias -Name $alias.Key -Value $alias.Value -Scope Global -Force -ErrorAction SilentlyContinue
    }
}
