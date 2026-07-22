---
name: LiteralPath default convention
description: Always use -LiteralPath over -Path for filesystem cmdlets in this codebase
type: feedback
---

Default to `-LiteralPath` instead of `-Path` for all filesystem cmdlets: `Get-Content`, `Test-Path`, `Get-Item`, `Move-Item`, `Copy-Item`, `Remove-Item`, etc.

**Why:** `-Path` treats `[`, `]`, `*`, `?` as wildcard characters. Files with bracketed names (user-named files, imported content, discussion docs) silently fail or match wrong paths. Discovered when a file named `[unsure of correct filename].md` broke `Get-Content` in the Colonel benchmark runner.

**How to apply:** Any time filesystem cmdlets operate on paths derived from the filesystem (e.g. from `Get-ChildItem`, user input, crawler output) — use `-LiteralPath`. Only use `-Path` when wildcard expansion is explicitly intended. This applies to the colonel, crawler, ignore modules, and any test/benchmark scripts.
