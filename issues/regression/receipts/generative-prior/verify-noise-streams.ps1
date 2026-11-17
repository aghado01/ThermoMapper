param(
    [Parameter(Mandatory)][string]$RngAssembly,
    [Parameter(Mandatory)][string]$Python,
    [Parameter(Mandatory)][string]$Output
)
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $true
$assemblyPath = (Resolve-Path -LiteralPath $RngAssembly).Path
$null = [Reflection.Assembly]::LoadFrom($assemblyPath)
$null = New-Item -ItemType Directory -Path $Output -Force
$reference = foreach ($seed in @(181204176, -2025086531, 181204183)) {
    $rng = [Maths.Rng.Xoshiro256PlusPlus]::new([int]$seed)
    $words = @(for ($i = 0; $i -lt 10; $i++) { $rng.NextUInt64().ToString([Globalization.CultureInfo]::InvariantCulture) })
    $rng = [Maths.Rng.Xoshiro256PlusPlus]::new([int]$seed)
    $bounded = @(for ($i = 0; $i -lt 10; $i++) { $rng.NextInt(30) })
    @{ seed = $seed; children = [Maths.Rng.SeedTree]::Derive($seed, 10); words = $words; bounded = $bounded }
}
$referencePath = Join-Path $Output 'csharp-reference.json'
@{ assembly = $assemblyPath; assembly_sha256 = (Get-FileHash -LiteralPath $assemblyPath -Algorithm SHA256).Hash.ToLowerInvariant(); cases = @($reference) } |
    ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $referencePath -Encoding utf8
& $Python (Join-Path $PSScriptRoot 'verify-noise-streams.py') --reference $referencePath --output (Join-Path $Output 'verification.json')
