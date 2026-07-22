---
name: feedback_portable_env_invocation
description: "pdenv tools reach the shell via the USER REGISTRY PATH now (shifted off the bootstrap/SHARED_ENV mechanism); SHARED_ENV_LOADED intentionally unset & not returning; a pdenv tool missing from a tool-shell probe could be stale-PATH OR not-installed — don't conclude from the tool shell"
metadata:
  node_type: memory
  type: feedback
  originSessionId: d1c8433d-0272-46f7-8f5f-b886757661a4
---

**Current mechanism (2026-06-24).** pdenv tools reach the shell via the **Windows user registry env
PATH** — Azriel added the `PORTABLE_ROOT` bins to the user-level PATH in the registry. This *replaced*
the bootstrap-profile / `SHARED_ENV` approach; the shift was forced by the recent Anthropic
app-store/walled-install update to Claude Code. `PowerShell.InheritedEnv.ps1` is dead — never source it.

**SHARED_ENV_LOADED is intentionally unset and is NOT coming back** as the mechanism (this supersedes
the earlier "mid-migration, will return" note). The registry PATH is the direction. Do **not** gate on
`$env:SHARED_ENV_LOADED` or assume the `dbd`/`drn`/`dtst` aliases exist. The project `CLAUDE.md` still
says "use a terminal where `SHARED_ENV_LOADED -eq $true`" / "prefer `dbd`/`drn`" — that guidance is
**stale**; ignore it.

**A tool-shell probe can't tell "not on PATH" from "not installed."** Two reasons a pdenv binary may
not resolve in the Claude tool shell: (a) **stale PATH snapshot** — a registry PATH change only reaches
processes started *after* it, so if the Claude app launched before the change its tool shells won't see
newly-added bins; or (b) the tool simply **isn't installed** in pdenv yet. Observed 2026-06-24: a
`Get-Command llama-server` probe found nothing — and Azriel first said llama.cpp was first-class in
pdenv, then was unsure he'd added it. So a single tool-shell probe is **not** evidence either way.
`$env:PORTABLE_ROOT` (`C:/Users/azrie/PDenv`) IS inherited; per-bin PATH entries may not be.

**How to apply.**
- Do NOT probe the tool shell to "confirm" a pdenv tool exists or is missing — it's ambiguous. Don't
  crawl `PORTABLE_ROOT` hunting for binaries either; ask Azriel if install status matters.
- If a tool command genuinely needs a pdenv binary and it won't resolve, use the absolute path under
  `$env:PORTABLE_ROOT`, or hand the command to Azriel's interactive shell.
- Bare `dotnet` has historically resolved in the tool shell; newer pdenv additions may not.

**llama.cpp — local embedder, now ON the user-registry pdenv PATH (added 2026-06-24).** The
local-embedding provider for Smart Connections (Obsidian "Custom Local (OpenAI Format)" →
`llama-server --embedding`) and the SPCX semantic pass — one binary, not LLamaSharp's bundled native
copy. Resolves in Azriel's interactive shell; the Claude **tool shell may still not see it until the
app restarts** (stale-snapshot caveat above) — so don't probe to "confirm" it. See the memory-graph
brief (vault-root `README.md`, mirrored by [[project_discussion_superrepo]]).

Related: [[feedback_design_before_build]].
