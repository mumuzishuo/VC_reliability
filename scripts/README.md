# Analysis scripts

Scripts should be organised by analysis stage and should use explicit input and output paths.

Recommended stages include:

1. raw audit
2. metadata and label validation
3. gene and perturbation matching
4. response matrix construction
5. baseline evaluation
6. calibration
7. retrieval and reproducibility
8. risk scoring
9. external validation

Scripts must not overwrite raw data. New outputs belong in a new versioned results directory.
