param(
    [Parameter(Mandatory)][string]$Revision,
    [Parameter(Mandatory)][string]$RunName,
    [Parameter(Mandatory)][string]$Python,
    [Parameter(Mandatory)][string]$Rows
)
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $true
$repo = (git -C $PSScriptRoot rev-parse --show-toplevel).Trim()
$pin = (git -C $repo rev-parse --verify "$Revision^{commit}").Trim()
$receipts = 'issues/regression/receipts/generative-prior'
$run = Join-Path $repo "artifacts/generative-prior/$RunName"
if (Test-Path -LiteralPath $run) { throw "Use a fresh run label: $run" }
$source = Join-Path $run 'source'
$null = New-Item -ItemType Directory -Path $source
$archive = Join-Path $run 'source.tar'
git -C $repo archive --format=tar "--output=$archive" $pin -- "$receipts/g5-noise.py" "$receipts/settings/g5-noise.json"
tar -xf $archive -C $source
Remove-Item -LiteralPath $archive
$pythonPath = (Resolve-Path -LiteralPath $Python).Path
$rowsPath = (Resolve-Path -LiteralPath $Rows).Path
[Diagnostics.Process]::GetCurrentProcess().PriorityClass = 'BelowNormal'
$script = Join-Path $source "$receipts/g5-noise.py"
# Verify inherited Windows priority inside the analysis process before any bootstrap.
$launch = 'import ctypes, runpy, sys; kernel = ctypes.windll.kernel32; kernel.GetCurrentProcess.restype = ctypes.c_void_p; kernel.GetPriorityClass.argtypes = [ctypes.c_void_p]; assert kernel.GetPriorityClass(kernel.GetCurrentProcess()) == 0x4000; sys.argv = sys.argv[1:]; runpy.run_path(sys.argv[0], run_name="__main__")'
& $pythonPath -c $launch $script --rows $rowsPath --output (Join-Path $run 'results') *> (Join-Path $run 'console.log')
@{ revision = $pin; priority = 'BelowNormal, checked inside Python'; workers = 1; rows = $rowsPath } |
    ConvertTo-Json | Set-Content -LiteralPath (Join-Path $run 'execution.json') -Encoding utf8
Write-Output "Completed saved-row replay at $run/results"
