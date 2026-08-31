# TFPT project skills

Project-local Cursor Agent Skills (`.cursor/skills/`). Loaded on demand when the task matches the skill description.

| Skill | Use when |
|-------|----------|
| [tfpt-experiment](tfpt-experiment/SKILL.md) | Work in `experiments/` — firewall, scorecard, standalone layout |
| [tfpt-empirical-search](tfpt-empirical-search/SKILL.md) | FRB/GW/pulsar/recovery preregistered searches |
| [tfpt-evidence-scorecard](tfpt-evidence-scorecard/SKILL.md) | Update `evidence_scorecard.json` rows and enums |
| [promote-to-verification](promote-to-verification/SKILL.md) | Graduate a finding to `verification/vN_*.py` + papers + ledger |
| [tfpt-deep-sync](tfpt-deep-sync/SKILL.md) | Parallel surface enumeration before/after vN integration |
| [rh-sync](rh-sync/SKILL.md) | Sync `rh/` after each RH round — INVENTORY, suite, README, paper, Lean |
| [agent-routing](agent-routing/SKILL.md) | Model tier per task, STATE handoffs, escalation signals, spend review |

## Rules (Cursor)

| Rule | Scope |
|------|-------|
| `tfpt-core.mdc` | **Always on** — invariants + skill routing |
| `agent-economics.mdc` | **Always on** — model ladder, escalation gates, context budget |
| `tfpt-workflow.mdc` | Glob: verification, TeX, changelog, README |
| `sync-maps.mdc` | Glob: verification, TeX, website, README, next.txt |
| `subagent-deep-sync.mdc` | Glob: verification, TeX, website |
| `website-sync.mdc` | Glob: website, verification, TeX |
| `rh-workspace.mdc` | Glob: rh/** — consolidated RH home rules |

## Hooks

`.cursor/hooks.json` wires three project hooks, all reading their policy from
`.cursor/hooks/model-tiers.json` and appending to `.cursor/agent-runs/routing.jsonl`:

| Hook | Event | Effect |
|------|-------|--------|
| `route-subagent.py` | `preToolUse` (`Task`) | rewrites the subagent model to the allowed tier |
| `guard-read.py` | `beforeReadFile` | blocks reads ≥ 400k chars, prices the large ones |
| `log-subagent.py` | `subagentStart` | records the model a subagent really started with |

Review with `python3 .cursor/hooks/routing_report.py`. All three fail open. Change policy in the JSON,
not in the scripts.
