## RepoAudit - Static & semantic analysis

`scripts/repo-audit.ps1`

Entry point that builds `src/repo-audit` (project: `projects/RepoAudit/RepoAudit.csproj`) on first run, then runs it. Default run writes `artifacts/project-health.md` and renders it inline. Flags:

- `-Validate` — strict mode; fail on warnings.
- `-NoGit` — skip git-history-based checks.
- `-Impact <path>` — focused analysis around a single path.
- `-NoDisplay` — suppress the auto-render of the health report.
- `-Rebuild` — force-rebuild RepoAudit before running. **Required after editing `src/repo-audit/*`:** the built exe is cached (only rebuilt on `-Rebuild` or if missing), so changes to the analyzer itself won't take effect without it.

Static analysis is still growing — it's a fast way to catch compile-time and integration problems before a full `dotnet build`.

**Usage**: repo-audit is a detector, but not entirely trustworthy. Its a valuable tool. Utilize it, but verify apparent claims against the evidence in the code it flags by inspection, and by targeted `dotnet build` spot-checks. When repo-audit inaccuracies are detected, escalate to user for potential bug fixes and/or enhancements to the tool, in order to make it more trustworthy going foward.

## test-harness - Parallel test runner

`src/test-harness/`

C# xUnit parallel fact runner; project at `projects/TestHarness.Runner/`. Discovers facts via `dotnet test --list-tests`, runs them concurrently up to `--max-workers`, and drops per-suite manifests + a `summary.json` into `artifacts/test-runs/<suite>/<stamp>/`. Invoke via:

```
dotnet run --project projects/TestHarness.Runner --
    --project <test.csproj> [--fixture <Class> | --filter <expr>]
    [--max-workers N] [--list-only]
```

`scripts/Invoke-ParallelFactBattery.ps1` is the earlier PowerShell prototype this C# runner replaced — no callers, and its `TestBatteries.psd1` catalog is missing. Don't use it.

## Shell environment and dotnet

Tools (`dotnet`, `git`, `python`, `npm`) are **not** on the system PATH. They are injected by the bootstrap profile (`PowerShell.InheritedEnv.ps1`). Always use a terminal where `$env:SHARED_ENV_LOADED -eq $true`.

- Use `$env:PORTABLE_ROOT` (not `$env:USERPROFILE`) as the path anchor.
- `netstat`, `taskkill`, and `Get-NetTCPConnection` are unavailable — use `Stop-Process -Name dotnet -Force -ErrorAction SilentlyContinue` to kill a running server.
- Use absolute `--project` paths when cwd may not be the repo root.
- **Never restart VizApi with `--no-build` after editing C# source** — always build first, then start.

```powershell
Stop-Process -Name dotnet -Force -ErrorAction SilentlyContinue
dotnet build projects/VizApi/VizApi.csproj -c Release
dotnet run --project projects/VizApi/VizApi.csproj --no-build -c Release
```

The bootstrap also registers portable-aware dotnet aliases (`dbd`, `drn`, `dtst`, `dw`, `dnew`, `dadd`, `dnv`) that call the portable exe directly — prefer these over bare `dotnet` in shell sessions.
