# Where the reasoning lives

The plans, the ADRs, and the plan of record are **not in this repo**. They live in the durable
memory layer:

| | |
|---|---|
| **Local** | `/data1/users/antonz/pipeline/bulkiRNA-memory` (clone it as a sibling of this repo) |
| **Hub** | `/data1/users/antonz/git/bulkiRNA-memory.git` |
| **Start here** | `plans/HANDOFF.md` — live state, the pain points, the ADR summaries, what is next |

Split out on 2026-08-24 at `a852658`, immediately after `v1.0.0`, using
`git subtree split --prefix=docs/_internal`. All 78 commits went with it, so each decision still
carries the commit that made it. Nothing was archived and nothing was rewritten here.

## Why

A plan and a package have different lifetimes. This repo is versioned, tagged, gated and installed,
and a release names exactly one commit (ADR-001). The reasoning accretes and gets corrected, and is
most useful when the whole argument is readable, including the parts that turned out wrong. In one
repo, every plan correction was a commit against the released artifact.

## What stayed here, and why

`AGENTS.md` and `CONVENTIONS.md` stay: they are rules that bind code in this repo, and an agent
working here needs them without a second clone. `docs/API_REFERENCE.md` and `docs/WORKFLOWS.md` stay
for the same reason.

The split is by lifetime, not importance. The strongest statements about this package are enforced in
its own code — `bulkirna_api()`, `bulkirna_stochastic()`, `gs_validate_master()`, and the
version-identity test — not in either set of documents.
