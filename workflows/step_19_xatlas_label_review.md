# Step 19: X-Atlas label and control-design review

Inputs: the versioned outputs from step 18 and the verified Figshare+ v3
`guide_library.csv` (SHA256 `1841320684c72f4a6c5f51a98ebaab0054711367a9ee81cc2167dd3fbbd6381d`).
The guide library is kept on HPC and excluded from Git.

## Current project calculations

- Exact Ensembl target ID overlap with the 1,264-task panel: 1,231 (97.39%).
- Of these, 1,198 have identical symbols and 33 have different symbols for
  the *same* Ensembl ID. This is an annotation discrepancy to review, not
  evidence of 33 different genes or observed responses.
- Four further targets match only by symbol but have a *different* Ensembl ID;
  they are excluded from the exact-ID cohort pending independent resolution.
- The 1,026 guide-library rows with a non-Ensembl target ID all have
  `target_gene_id=Non-Targeting` and `target_gene_name=Non-Targeting`.
  This is the number of **designed guide-pair rows**, not the number of
  non-targeting *cells* in either cell line or in each GEM.
- Complete label pairs for the 33 and four appear in
  `qc/v01_xatlas_guide_design/overlap_by_gene_id.tsv` and
  `qc/v01_xatlas_guide_design/symbol_only_diagnostic.tsv`. The local review
  command/output is archived as `logs/v01_xatlas_guide_design/label_review.log`.

Verify from the archived guide-library file, without accessing expression
or target outcomes:

```bash
cd "$HOME/VC_reliability_repo"
(cd qc/v01_xatlas_guide_design && sha256sum --status -c SHA256SUMS)
python3 - <<'PY'
import csv, re
from collections import Counter
root = "qc/v01_xatlas_guide_design/"
with open(root + "overlap_by_gene_id.tsv", encoding="utf-8", newline="") as f:
    matched = list(csv.DictReader(f, delimiter="\t"))
assert len(matched) == 1231
assert sum(r["symbol_agrees"] == "0" for r in matched) == 33
with open(root + "symbol_only_diagnostic.tsv", encoding="utf-8", newline="") as f:
    assert len(list(csv.DictReader(f, delimiter="\t"))) == 4
labels = Counter()
with open("metadata/v01_external_catalog/guide_library.csv",
          encoding="utf-8-sig", newline="") as f:
    for r in csv.DictReader(f):
        gene_id = r["target_gene_id"].strip()
        if gene_id and not re.fullmatch(r"ENSG[0-9]{11}(?:\.[0-9]+)?", gene_id):
            labels[(gene_id, r["target_gene_name"].strip())] += 1
assert labels == {("Non-Targeting", "Non-Targeting"): 1026}
print("XATLAS_LABEL_REVIEW_OK")
PY
```

These are guide-*design* checks. Eligibility for independent validation remains
unresolved until observed perturbation cell counts, guide capture, control cells
per cell line/GEM, biological replicates, shared expression gene IDs and
dataset independence are checked. Freeze target selection and prediction/risk
rules before using external responses. Do not download either large H5AD
on the strength of this guide-design overlap alone.
