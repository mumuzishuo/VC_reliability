# Decision log

## 2026-11-20

### Decision

Start with a frozen data inventory and analysis panel before model training or external validation.

### Reason

The reliability question depends on a stable perturbation panel, stable gene space, matched controls and traceable batch structure. Changing these after inspecting prediction results would make model comparisons and risk evaluation difficult to audit.

### Current evidence

- Four environments contain 966,728 cells in total.
- 2,053 perturbation genes are shared across all four environments.
- 6,642 expression genes are shared across all four environments.
- 1,264 perturbation tasks meet the current minimum of 30 cells per environment.
- An independent recomputation of 35 tasks found a maximum response difference of 7.34e-8.
- No low-control GEM group has been identified in the current audit.

### Pending verification

- Exact raw-file accession and checksum registry
- Complete guide-level metadata
- Per-task cell and control counts
- Batch, sample and donor confounding checks
- Reproducible audit script on the source files
