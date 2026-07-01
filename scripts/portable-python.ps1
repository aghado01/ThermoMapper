
$venv = "$env:PORTABLE_ROOT/UserGithub/PowerShellCore/ps.core.pwshspc/.venv"

& $py -m venv $venv
& "$venv/Scripts/Activate.ps1"

$py_venv = "$venv/Scripts/python.exe"
