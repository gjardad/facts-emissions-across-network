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

**EXCEPTION:** You MAY read `inferring_emissions/` code files (including `analysis/`, `preprocess/`, `utils/`) when the explicit purpose is to copy or directly adapt that code into `facts-emissions-across-network/`. This exception exists because this project is meant to be standalone — all necessary code from `inferring_emissions/` should be copied here so that someone can reproduce results without knowing `inferring_emissions/` exists.

## Allocation & Uncertainty Notes

- Allocation and uncertainty propagation operate at the **CRF category × year** level, in both deployment and training.
- In training, three CRF categories contain both emitters and non-emitters:
  - **Paper**: NACE 2-digit 17 and 18
  - **Refinery**: NACE 2-digit 19 only
  - **Metals**: NACE 2-digit 24 and 25


