# Lessons

## 1. Pass CLAUDE.md constraints to subagents explicitly

**Mistake:** Spawned a subagent to explore `inferring_emissions/` without constraining it to `paper/winter26_version/`. The subagent read code and analysis scripts, producing claims based on potentially outdated sources.

**Rule:** Before launching any subagent that touches `inferring_emissions/`, include the hard constraint from CLAUDE.md verbatim in the subagent prompt: only read files under `paper/winter26_version/`.

**Generalization:** Any project-specific constraint in CLAUDE.md must be forwarded to subagent prompts. Subagents do not inherit CLAUDE.md context automatically.
