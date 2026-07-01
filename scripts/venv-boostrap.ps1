# bootstrap local python from portable python

$py = "$env:PORTABLE_PYTHON/python.exe"
# $venv = "$env:PORTABLE_ROOT/UserGithub/PowerShellCore/ps.core.pwshspc/.venv"
$venv = "C:\Users\azrie\PDenv\UserGithub\PowerShellCore\ps.core.pwshspc\notebook\mvp\.venv"

& $py -m venv $venv
& "$venv/Scripts/Activate.ps1"

$py_venv = "$venv/Scripts/python.exe"

& $py_venv -m pip install --upgrade pip

& $py_venv -m pip install -r requirements.txt
