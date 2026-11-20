# Project roadmap

## Stage 0 is complete

The repository now contains version-control and reproducibility documentation. The initial commit is:

`9a216de73724b7f067a36419eb1fb25c0903c416`

The working branch for analysis development is `dev`.

## Stage 1 is the next active stage

The immediate task is to freeze the current data inventory and analysis panel before adding new models or external validation data.

### Stage 1 inputs

- Four current Perturb-seq environments
- Existing raw audit results
- Existing perturbation and gene overlap summaries
- Existing control and GEM-group checks

### Stage 1 outputs

- Current data inventory
- Analysis-panel definition
- Versioned audit configuration
- Source and checksum manifest
- Reproducible audit script
- QC decision record

### Stage 1 acceptance criteria

- Cell, gene and GEM-group counts are recorded for every environment.
- Perturbation labels and non-targeting controls are explicitly identified.
- Control matching fields are recorded.
- The 1,264-task panel is defined without using target prediction results.
- No low-control GEM group is present in the current audit.
- Raw files are never overwritten.
- Every derived table has an input manifest and generation commit.

## Stage 2 response matrix

Construct `Delta_{p,e,g}` using the frozen normalisation and matched controls. Record cell counts, control counts, GEM composition and perturbation strength for every task.

## Stage 3 baseline evaluation

Evaluate zero response, simple average, common-source average and calibrated predictions. Report task-level and aggregate metrics.

## Stage 4 reliability scoring

Freeze a target-unknown risk score using only source-side information. Candidate signals include cell count, control support, perturbation strength, source agreement, split-half reproducibility and calibration stability.

## Stage 5 external validation

First inspect metadata and processed summaries from Jiang and X-Atlas. Download full single-cell files only if coverage, controls, batch structure and task overlap pass the metadata gate.

Use an external environment as the target only after prediction, calibration, risk scoring and interpretation rules are frozen.

## Stage 6 final reproducibility package

Publish scripts, configurations, manifests, logs, result tables, figures, frozen protocol and a complete commit history.

## Analysis safeguards

- Always report the zero-response baseline.
- Preserve simple and calibrated baselines.
- Run label permutation and chance-level checks.
- Stratify by cell number and perturbation strength.
- Check GEM, batch, donor and stimulus confounding.
- Distinguish data-processing findings, statistical associations, prediction performance and causal interpretation.
