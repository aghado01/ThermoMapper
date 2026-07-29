[CmdletBinding()]
param (
    [Parameter(Mandatory = $false)]
    [string]$InputDir = ".\reference-manual\.lake\build\doc", # Verso HTML output path

    [Parameter(Mandatory = $false)]
    [string]$OutputDir = ".\codex-scientiae\lean-manual",    # Your corpus destination

    [Parameter(Mandatory = $false)]
    [string]$PandocFormat = "gfm"                            # GitHub Flavored Markdown
)

# Ensure output root exists
if (-not (Test-Path $OutputDir)) {
    New-Item -ItemType Directory -Path $OutputDir -Force | Out-Null
}

# Resolve absolute paths to avoid context loss inside the parallel runspace
$resolvedInDir = (Resolve-Path $InputDir).Path
$resolvedOutDir = (Resolve-Path $OutputDir).Path

$HtmlFiles = Get-ChildItem -Path $resolvedInDir -Filter "*.html" -Recurse

if ($HtmlFiles.Count -eq 0) {
    Write-Warning "No HTML files found in $resolvedInDir. Ensure Verso has run."
    exit
}

Write-Host "Discovered $($HtmlFiles.Count) HTML documents. Initializing parallel Pandoc conversion..." -ForegroundColor Cyan

# Fire up the parallel execution pool
$HtmlFiles | ForEach-Object -Parallel {
    $file = $_
    $inRoot = $using:resolvedInDir
    $outRoot = $using:resolvedOutDir
    $format = $using:PandocFormat

    # Reconstruct the exact directory tree in the output folder
    $relativePath = $file.DirectoryName.Substring($inRoot.Length).TrimStart('\', '/')
    $targetDir = Join-Path $outRoot $relativePath

    if (-not (Test-Path $targetDir)) {
        New-Item -ItemType Directory -Path $targetDir -Force | Out-Null
    }

    $outFilePath = Join-Path $targetDir ($file.BaseName + ".md")

    try {
        # --wrap=none is critical for LLM corpus building; it stops pandoc from 
        # inserting hard line breaks at 80 characters, keeping paragraphs contiguous.
        pandoc.exe -f html -t $format --wrap=none $file.FullName -o $outFilePath
        Write-Host "[OK] Converted: $($file.Name)" -ForegroundColor DarkGray
    }
    catch {
        Write-Error "[FAIL] Pandoc choked on $($file.Name): $_"
    }

} -ThrottleLimit ([Environment]::ProcessorCount)

Write-Host "Corpus successfully generated at $resolvedOutDir." -ForegroundColor Green