---
name: naming conventions — prefer functional over project-prefixed names
description: User prefers functional, descriptive names (repo-audit, file-ownership) over project-scoped prefixes (spcx-audit, ps-core-x).
type: feedback
originSessionId: b3770f94-98fe-4d6a-a07c-e783b9e7a2f6
---
Use functional names that describe what a tool/command/namespace does, not which project it belongs to.

**Why:** names like `spcx-audit` are opaque and tied to a project codename; `repo-audit` is self-describing and transferable.

**How to apply:** when naming CLI commands, scripts, namespaces, or assemblies, lead with the function (repo-audit, code-analysis, file-ownership) not the project prefix (spcx-, ps-core-, pwshspc-).
