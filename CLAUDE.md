## Project Description

See `THIS_PROJECT.md` for a full description of the project, its research questions, and data sources.

## Hardware Setup

Three desktops are used in this project: local 1, local 2, and RMD (remote desktop).

- **RMD**: Has the full NBB data (B2B, Customs, Annual Accounts, PRODCOM). Accessible only via VPN through local 2. No web browser, but connected to GitHub.
- **Local 1**: Personal desktop with Claude Code and Cursor. Has a **downsampled** version of NBB data and the **full training sample** (copied from RMD via local 2 → cloud → local 1).
- **Local 2**: Bridge machine with VPN access to RMD and regular web browser.

When copying files from RMD to local 1: RMD → local 2 → cloud (Dropbox/Claude) → local 1.

## Cross-Project Reference: `inferring_emissions/`

**HARD CONSTRAINT:** When reading anything from `inferring_emissions/`, ONLY read files under `paper/winter26_version/`. Do NOT read code, analysis scripts, or earlier paper versions — they may be outdated. This applies to all subagents and exploration tasks. Pass this constraint explicitly when delegating.

## Workflow Orchestration

### 1. Plan Node Default
- Enter plan mode for ANY non-trivial task (3+ steps or architectural decisions)
- If something goes sideways, STOP and re-plan immediately - don't keep pushing
- Use plan mode for verification steps, not just building
- Write detailed specs upfront to reduce ambiguity

### 2. Subagent Strategy
- Use subagents liberally to keep main context window clean
- Offload research, exploration, and parallel analysis to subagents
- For complex problems, throw more compute at it via subagents
- One task per subagent for focused execution

### 3. Self-Improvement Loop
- After ANY correction from the user: update `tasks/lessons.md` with the pattern
- Write rules for yourself that prevent the same mistake
- Ruthlessly iterate on these lessons until mistake rate drops
- Review lessons at session start for relevant project

### 4. Verification Before Done
- Never mark a task complete without proving it works
- Diff behavior between main and your changes when relevant
- Ask yourself: "Would a staff engineer approve this?"
- Run tests, check logs, demonstrate correctness
