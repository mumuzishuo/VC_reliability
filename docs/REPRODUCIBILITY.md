# Reproducibility and Version Control

## Version policy

Each completed analysis stage receives a Git commit with:

- a descriptive commit message
- input file identifiers
- parameter file
- software environment
- output locations
- validation status
- SHA256 manifest when files are produced

Stable milestones receive a Git tag such as:

- `v0.1-raw-audit`
- `v0.2-response-matrix`
- `v0.3-baseline-evaluation`
- `v0.4-calibration`
- `v1.0-external-validation`

## Data leakage controls

The following must be recorded before external validation:

- prediction method
- calibration coefficients
- risk score definition
- evaluation metrics
- result interpretation rules

Target-environment ground-truth responses must not influence these frozen objects.

## Required audit trail

Each run should preserve:

- command or Slurm submission script
- parameter values
- input manifest
- software versions
- stdout and stderr logs
- output manifest
- commit SHA
- random seed where applicable

Raw datasets remain outside the repository. Their public source, version, file size and SHA256 are recorded in `data/raw/` and `manifests/`.
