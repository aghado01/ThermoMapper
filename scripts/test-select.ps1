<#
.SYNOPSIS
    Runs the tests a change can affect: the test classes a name graph over the project's sources selects, the slow tier
    left out, and the pinned tier added at a commit.

.DESCRIPTION
    Selection. The changed files come from git (the working tree and untracked files against -Base, or the changes one
    -Commit made) or from -Changed. The graph covers the test project's compile closure: the Compile items of its
    csproj, of every project it references, transitively, and of the Directory.Build files above them, with each SDK
    project's default items (every .cs file under its folder, outside bin, obj and hidden folders) unless
    EnableDefaultCompileItems or EnableDefaultItems turns them off. A file is
    affected when it changed or when it names a top-level type, a delegate or an extension method declared in an
    affected file; comments and strings are not read, and a name read as a member (after a dot, or called) is not a
    type's unless its dotted chain starts at a namespace of the closure. The selected classes are the test classes
    declared in affected files. A change to a csproj in the closure or to a Directory.Build file selects every class;
    a changed file outside the closure selects nothing.

    The graph reads names, not symbols. It over-selects where two types share a name, and through hubs: a file that
    names a type is affected whether or not the path it runs reaches the change. It under-selects where code reaches a
    type without naming it on the way, by reflection or by a nested type named bare from a derived class; the pinned
    tier at each commit and the whole assembly at a slice's end are the net for that.

    Scopes, the rule in CONTRIBUTING.md's testing standards:
      Edit    the selected classes, the slow tier left out
      Commit  the selected classes, the slow tier left out, and the pinned tier ([Trait("Tier", "Pinned")])
      Slice   the whole assembly
    The slow tier is the list beside the project, slow-tests.txt, which -WriteSlowList writes from a profile's trx.

    A run builds and tests in Release into -ArtifactsRoot, at below-normal priority, with xUnit capped at -Threads
    and MSBuild's build servers off, and leaves a trx beside its runsettings under <ArtifactsRoot>test-select\.

.EXAMPLE
    scripts/test-select.ps1 -ListOnly
    scripts/test-select.ps1 -Scope Commit -Changed src/maths/regression/spline/bass/BassNoise.cs -ArtifactsRoot D:\aghado01\ThermoMapper\artifacts\chip-x\
    scripts/test-select.ps1 -Commit 37137f2 -ListOnly -Explain
    scripts/test-select.ps1 -Scope Slice -ArtifactsRoot D:\aghado01\ThermoMapper\artifacts\chip-x\
    scripts/test-select.ps1 -WriteSlowList artifacts\chip-x\profile.trx -SlowSeconds 10
#>
[CmdletBinding()]
param(
    [string]   $Project       = 'tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj',
    [ValidateSet('Edit', 'Commit', 'Slice')]
    [string]   $Scope         = 'Edit',
    [string]   $Base          = 'HEAD',
    [string]   $Commit,
    [string[]] $Changed,
    [string]   $SourceRoot,
    [string]   $ArtifactsRoot,
    [int]      $Threads       = 6,
    [switch]   $IncludeSlow,
    [switch]   $NoBuild,
    [switch]   $ListOnly,
    [switch]   $Explain,
    [string]   $WriteSlowList,
    [double]   $SlowSeconds   = 10
)

$ErrorActionPreference = 'Stop'
$GitRoot = Split-Path $PSScriptRoot
if (-not $SourceRoot) { $SourceRoot = $GitRoot }
$SourceRoot = (Resolve-Path -LiteralPath $SourceRoot).Path
$projectPath = (Resolve-Path -LiteralPath (Join-Path $SourceRoot $Project)).Path
$slowListPath = Join-Path (Split-Path $projectPath) 'slow-tests.txt'

function Get-FullPath([string] $path) { [IO.Path]::GetFullPath($path) }
function Get-RepoRelative([string] $full) { [IO.Path]::GetRelativePath($SourceRoot, $full).Replace('\', '/') }

# ── The slow list, written from a profile ────────────────────────────────────────────────────────────────────────────
if ($WriteSlowList) {
    [xml] $trx = Get-Content -LiteralPath $WriteSlowList -Raw
    $ns = @{ t = 'http://microsoft.com/schemas/VisualStudio/TeamTest/2010' }
    $methodOf = @{}
    foreach ($node in (Select-Xml -Xml $trx -XPath '//t:UnitTest' -Namespace $ns)) {
        $m = $node.Node.TestMethod
        $methodOf[$node.Node.id] = "$($m.className).$($m.name)"
    }
    $seconds = @{}; $cases = @{}
    foreach ($node in (Select-Xml -Xml $trx -XPath '//t:UnitTestResult' -Namespace $ns)) {
        $key = $methodOf[$node.Node.testId]
        $seconds[$key] = [double]$seconds[$key] + [TimeSpan]::Parse($node.Node.duration).TotalSeconds
        $cases[$key] = [int]$cases[$key] + 1
    }
    $total = ($seconds.Values | Measure-Object -Sum).Sum
    $slow = $seconds.GetEnumerator() | Where-Object Value -ge $SlowSeconds | Sort-Object Value -Descending
    $slowTotal = ($slow | Measure-Object -Property Value -Sum).Sum
    $times = (Select-Xml -Xml $trx -XPath '//t:Times' -Namespace $ns).Node
    $lines = @(
        "# The slow tier of $(Split-Path $Project -Leaf): each test method whose cases summed to at least $SlowSeconds s in one"
        "# profile of the whole assembly, $(Split-Path $WriteSlowList -Leaf), run $($times.start)."
        "# $($slow.Count) of $($seconds.Count) methods, $([math]::Round($slowTotal)) s of $([math]::Round($total)) s summed over the cases."
        "# Written by scripts/test-select.ps1 -WriteSlowList; rewrite it from a new profile rather than by hand."
        "# Edit and Commit runs leave these out; a slice's end runs them with the rest."
    )
    foreach ($e in $slow) { $lines += ('{0}  # {1:0.0} s over {2} case(s)' -f $e.Key, $e.Value, $cases[$e.Key]) }
    Set-Content -LiteralPath $slowListPath -Value $lines -Encoding utf8NoBOM
    Write-Host "slow tier: $($slow.Count) methods, $([math]::Round($slowTotal)) of $([math]::Round($total)) s -> $slowListPath"
    exit 0
}

# ── The compile closure ──────────────────────────────────────────────────────────────────────────────────────────────
# Each Include pattern is kept as well as expanded, so that a deleted file can be placed in or out of the closure.
$compilePatterns = [Collections.Generic.List[object]]::new()

function Expand-Compile([string] $pattern, [string] $dir, [switch] $Record) {
    $p = $pattern.Replace('/', '\')
    if (-not [IO.Path]::IsPathRooted($p)) { $p = Join-Path $dir $p }
    $leaf = Split-Path $p -Leaf
    if ($p.Contains('**')) {
        $root = Get-FullPath ($p.Substring(0, $p.IndexOf('**')))
        if ($Record) { $compilePatterns.Add(@{ Dir = $root; Recurse = $true; Leaf = $leaf }) }
        if (Test-Path -LiteralPath $root) { Get-ChildItem -LiteralPath $root -Recurse -File -Filter $leaf | ForEach-Object FullName }
    }
    else {
        $parent = Get-FullPath (Split-Path $p)
        if ($Record) { $compilePatterns.Add(@{ Dir = $parent; Recurse = $false; Leaf = $leaf }) }
        if ($leaf.Contains('*') -or $leaf.Contains('?')) {
            if (Test-Path -LiteralPath $parent) { Get-ChildItem -LiteralPath $parent -File -Filter $leaf | ForEach-Object FullName }
        }
        elseif (Test-Path -LiteralPath $p) { Get-FullPath $p }
    }
}

function Test-InCompilePatterns([string] $full) {
    foreach ($c in $compilePatterns) {
        $under = if ($c.Recurse) { $full.StartsWith($c.Dir.TrimEnd('\') + '\', [StringComparison]::OrdinalIgnoreCase) }
                 else { [string]::Equals((Split-Path $full), $c.Dir.TrimEnd('\'), [StringComparison]::OrdinalIgnoreCase) }
        if ($under -and (Split-Path $full -Leaf) -like $c.Leaf) { return $true }
    }
    return $false
}

# The last value an MSBuild file gives a boolean property: $true, $false, or $null where it sets none.
function Get-BooleanProperty([xml] $xml, [string] $name) {
    $node = @($xml.SelectNodes("//PropertyGroup/$name")) | Select-Object -Last 1
    if (-not $node) { return $null }
    return -not ($node.InnerText.Trim() -ieq 'false')
}

# An SDK project compiles every .cs file under its folder unless EnableDefaultCompileItems or EnableDefaultItems is
# false, in the project or else in the nearest Directory.Build.props above it, which MSBuild reads first.
function Test-DefaultCompileItems([string] $path, [xml] $xml) {
    if (-not $xml.Project.Sdk) { return $false }
    $sources = @($xml)
    $d = Split-Path $path
    while ($d -and $d.StartsWith($SourceRoot, [StringComparison]::OrdinalIgnoreCase)) {
        $props = Join-Path $d 'Directory.Build.props'
        if (Test-Path -LiteralPath $props) { $sources += [xml](Get-Content -LiteralPath $props -Raw); break }
        $d = Split-Path $d
    }
    foreach ($source in $sources) {
        if ((Get-BooleanProperty $source 'EnableDefaultItems') -eq $false) { return $false }
        $compile = Get-BooleanProperty $source 'EnableDefaultCompileItems'
        if ($null -ne $compile) { return $compile }
    }
    return $true
}

# The SDK's default compile items: **/*.cs under the project's folder, less bin, obj and hidden folders, as its
# default excludes leave them out.
function Get-DefaultCompileItems([string] $dir) {
    $compilePatterns.Add(@{ Dir = (Get-FullPath $dir); Recurse = $true; Leaf = '*.cs' })
    foreach ($f in (Get-ChildItem -LiteralPath $dir -Recurse -File -Filter '*.cs')) {
        $relative = [IO.Path]::GetRelativePath($dir, $f.FullName)
        $folders = @(Split-Path $relative) -split '[\\/]' | Where-Object { $_ }
        if ($folders | Where-Object { $_ -ieq 'bin' -or $_ -ieq 'obj' -or $_.StartsWith('.') }) { continue }
        $f.FullName
    }
}

function Read-ProjectFile([string] $path, [hashtable] $closure) {
    $dir = Split-Path $path
    [xml] $xml = Get-Content -LiteralPath $path -Raw
    $removes = [Collections.Generic.HashSet[string]]::new([StringComparer]::OrdinalIgnoreCase)
    foreach ($c in $xml.SelectNodes('//Compile[@Remove]')) {
        foreach ($part in $c.Remove.Split(';')) { Expand-Compile $part.Replace('$(MSBuildThisFileDirectory)', "$dir\") $dir | ForEach-Object { [void]$removes.Add($_) } }
    }
    $includes = [Collections.Generic.List[string]]::new()
    if ($path.EndsWith('.csproj', [StringComparison]::OrdinalIgnoreCase) -and (Test-DefaultCompileItems $path $xml)) {
        foreach ($f in (Get-DefaultCompileItems $dir)) { $includes.Add($f) }
    }
    foreach ($c in $xml.SelectNodes('//Compile[@Include]')) {
        foreach ($part in $c.Include.Split(';')) {
            foreach ($f in (Expand-Compile $part.Replace('$(MSBuildThisFileDirectory)', "$dir\") $dir -Record)) { $includes.Add($f) }
        }
    }
    foreach ($f in $includes) { if (-not $removes.Contains($f)) { [void]$closure.Files.Add($f) } }
    return $xml
}

function Get-Closure([string] $csproj) {
    $closure = @{
        Files      = [Collections.Generic.HashSet[string]]::new([StringComparer]::OrdinalIgnoreCase)
        BuildFiles = [Collections.Generic.HashSet[string]]::new([StringComparer]::OrdinalIgnoreCase)
    }
    $pending = [Collections.Generic.Queue[string]]::new(); $pending.Enqueue($csproj)
    while ($pending.Count -gt 0) {
        $path = $pending.Dequeue()
        if (-not $closure.BuildFiles.Add($path)) { continue }
        $xml = Read-ProjectFile $path $closure
        foreach ($r in $xml.SelectNodes('//ProjectReference[@Include]')) {
            $pending.Enqueue((Get-FullPath (Join-Path (Split-Path $path) $r.Include.Replace('/', '\'))))
        }
        # MSBuild imports the nearest Directory.Build.props and .targets above each project.
        foreach ($name in 'Directory.Build.props', 'Directory.Build.targets') {
            $d = Split-Path $path
            while ($d -and $d.StartsWith($SourceRoot, [StringComparison]::OrdinalIgnoreCase)) {
                $candidate = Join-Path $d $name
                if (Test-Path -LiteralPath $candidate) {
                    if ($closure.BuildFiles.Add($candidate)) { [void](Read-ProjectFile $candidate $closure) }
                    break
                }
                $d = Split-Path $d
            }
        }
    }
    return $closure
}

# ── Names: what each file declares and reads (scripts/test-select.cs) ────────────────────────────────────────────────
if (-not ('SourceGraph' -as [type])) { Add-Type -Path (Join-Path $PSScriptRoot 'test-select.cs') }

$closure = Get-Closure $projectPath
$facts = [Collections.Generic.Dictionary[string, SourceFacts]]::new([StringComparer]::OrdinalIgnoreCase)
foreach ($file in $closure.Files) { $facts[$file] = [SourceGraph]::Analyze([IO.File]::ReadAllText($file)) }
$testClasses = @{}
foreach ($e in $facts.GetEnumerator()) { if ($e.Value.TestClasses.Count -gt 0) { $testClasses[$e.Key] = @($e.Value.TestClasses) } }
$allClasses = @($testClasses.Values | ForEach-Object { $_ } | Sort-Object -Unique)

# ── The change ───────────────────────────────────────────────────────────────────────────────────────────────────────
function Invoke-Git { $out = & git -C $GitRoot @args; if ($LASTEXITCODE -ne 0) { throw "git $args failed" }; $out }

$before = $null
if ($Changed) {
    # Commas separate paths too, since pwsh -File passes a list as one string.
    $paths = @($Changed | ForEach-Object { $_ -split ',' } | ForEach-Object { $_.Trim() })
}
elseif ($Commit) {
    $before = "$Commit~1"
    $paths = Invoke-Git diff --name-only --no-renames $before $Commit
}
else {
    $before = $Base
    $paths = @(Invoke-Git diff --name-only --no-renames $Base) + @(Invoke-Git ls-files --others --exclude-standard)
}
$paths = @($paths | Where-Object { $_ } | ForEach-Object { $_.Replace('\', '/') } | Sort-Object -Unique)

$changedFiles = [Collections.Generic.List[string]]::new()
$deletedNames = [Collections.Generic.List[string]]::new()
$why = [Collections.Generic.Dictionary[string, string]]::new([StringComparer]::OrdinalIgnoreCase)
$everything = $false; $inClosure = 0; $outside = @()
foreach ($rel in $paths) {
    $full = Get-FullPath (Join-Path $SourceRoot $rel)
    if ($closure.BuildFiles.Contains($full)) { $everything = $true; $inClosure++; continue }
    $wasCompiled = $closure.Files.Contains($full)
    $wasDeleted = -not $wasCompiled -and -not (Test-Path -LiteralPath $full) -and (Test-InCompilePatterns $full)
    if (-not $wasCompiled -and -not $wasDeleted) { $outside += $rel; continue }
    $inClosure++
    if ($wasCompiled) { $changedFiles.Add($full) }
    elseif ($before) {
        # A deleted file's names still reach the files that named them. A file that stays is read as it is now: a name
        # renamed in it leaves its old name's readers failing to build, which the build reports.
        $old = & git -C $GitRoot show "${before}:$rel" 2>$null
        if ($LASTEXITCODE -eq 0) { $deletedNames.AddRange([SourceGraph]::Analyze($old -join "`n").Declared) }
    }
}
$affected = [SourceGraph]::Affected($facts, $changedFiles, $deletedNames, $why)

$selected = if ($everything) { $allClasses } else {
    @($testClasses.Keys | Where-Object { $affected.Contains($_) } | ForEach-Object { $testClasses[$_] } | Sort-Object -Unique)
}

# ── The filter ───────────────────────────────────────────────────────────────────────────────────────────────────────
$slowMethods = @()
if (Test-Path -LiteralPath $slowListPath) {
    $slowMethods = @(Get-Content -LiteralPath $slowListPath | Where-Object { $_ -and -not $_.StartsWith('#') } | ForEach-Object { ($_ -split '\s+')[0] })
}
$noSlow = ($slowMethods | ForEach-Object { "FullyQualifiedName!=$_" }) -join '&'
$classPart = ($selected | ForEach-Object { "FullyQualifiedName~$_." }) -join '|'
$all = $selected.Count -eq $allClasses.Count

$filter = switch ($Scope) {
    'Slice' { '' }
    default {
        $parts = @()
        if (-not $all) { $parts += "($classPart)" }
        if (-not $IncludeSlow -and $noSlow) { $parts += "($noSlow)" }
        # $null selects no test; '' selects every test, which is the routine filter when every class is selected and no
        # slow method is left out.
        $routine = if ($selected.Count -gt 0) { $parts -join '&' } else { $null }
        if ($Scope -eq 'Commit') {
            if ($null -eq $routine) { 'Tier=Pinned' } elseif ($routine -eq '') { '' } else { "($routine)|(Tier=Pinned)" }
        }
        else { $routine }
    }
}

$source = if ($Changed) { 'the paths given' } elseif ($Commit) { "commit $Commit" } else { "the working tree and untracked files against $Base" }
Write-Host "changes: $($paths.Count) path(s) from $source; $inClosure in the compile closure, $($outside.Count) outside it"
if ($everything) { Write-Host 'a build file changed: every class is selected' }
Write-Host "affected files: $($affected.Count) of $($closure.Files.Count); selected classes: $($selected.Count) of $($allClasses.Count)"
Write-Host "scope ${Scope}: $(if ($Scope -eq 'Slice') { 'the whole assembly' } else { "selected classes$(if (-not $IncludeSlow) { ", without the slow tier's $($slowMethods.Count) methods" })$(if ($Scope -eq 'Commit') { ', with the pinned tier' })" })"

if ($Explain) {
    foreach ($file in ($testClasses.Keys | Where-Object { $affected.Contains($_) } | Sort-Object)) {
        $chain = @(); $f = $file
        for ($i = 0; $i -lt 50 -and $why.ContainsKey($f) -and $why[$f] -ne 'changed'; $i++) {
            $n = $why[$f]; $chain += $n; $f = $why[":$n"]
        }
        Write-Host ("  {0} <- {1}{2}" -f (Get-RepoRelative $file), ($chain -join ' <- '), $(if ($f -and $f -ne $file) { " <- $(Get-RepoRelative $f)" } else { '' }))
    }
}

# An Edit run needs a selected class; a Commit run needs a change in the closure, its pinned tier being the net for
# what the graph cannot see.
$nothing = ($Scope -eq 'Edit' -and $selected.Count -eq 0) -or ($Scope -eq 'Commit' -and $inClosure -eq 0)
if ($ListOnly -or $nothing) {
    if ($ListOnly) { $selected | ForEach-Object { Write-Host "  $_" }; Write-Host "filter: $filter" }
    else { Write-Host 'nothing to run' }
    exit 0
}

# ── The run ──────────────────────────────────────────────────────────────────────────────────────────────────────────
if (-not $ArtifactsRoot) { throw '-ArtifactsRoot is required to run: builds and results go there, as the briefs set it.' }
$ArtifactsRoot = (Get-FullPath $ArtifactsRoot).TrimEnd('\') + '\'
$stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
$runDir = Join-Path $ArtifactsRoot "test-select\$stamp-$($Scope.ToLowerInvariant())"
New-Item -ItemType Directory -Force -Path $runDir | Out-Null
$settings = Join-Path $runDir 'filter.runsettings'
$filterXml = if ($filter) { "<TestCaseFilter>$([Security.SecurityElement]::Escape($filter))</TestCaseFilter>" } else { '' }
Set-Content -LiteralPath $settings -Encoding utf8NoBOM -Value @"
<RunSettings>
  <RunConfiguration>$filterXml</RunConfiguration>
  <xUnit><MaxParallelThreads>$Threads</MaxParallelThreads></xUnit>
</RunSettings>
"@

$psi = [Diagnostics.ProcessStartInfo]::new('dotnet')
foreach ($a in @('test', $projectPath, '-c', 'Release', "-p:ArtifactsRoot=$ArtifactsRoot", '--disable-build-servers',
                 '--settings', $settings, '--logger', "trx;LogFileName=$Scope.trx", '--results-directory', $runDir)) {
    $psi.ArgumentList.Add($a)
}
if ($NoBuild) { $psi.ArgumentList.Add('--no-build') }
$psi.UseShellExecute = $false
$psi.Environment['MSBUILDDISABLENODEREUSE'] = '1'

# The child inherits the priority class it is created under: lower ours for the start, then restore it.
$self = Get-Process -Id $PID
$prior = $self.PriorityClass
$self.PriorityClass = 'BelowNormal'
try { $proc = [Diagnostics.Process]::Start($psi) } finally { $self.PriorityClass = $prior }
$clock = [Diagnostics.Stopwatch]::StartNew()
$proc.WaitForExit()
Write-Host ("test-select: exit {0}, {1:0} s, results in {2}" -f $proc.ExitCode, $clock.Elapsed.TotalSeconds, $runDir)
exit $proc.ExitCode
