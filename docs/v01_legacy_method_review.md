# Legacy predictor review v01 (provisional, 2026-11-27)

This records a review of the legacy internal pilot. Step 28 archived the
source files byte for byte at commit `f999294`, with hashes in
`manifests/v01_legacy_code.sha256`. Step 29 prepares the four-source
development analysis. The external evaluation protocol remains open.

## Source identity reported by the HPC

| File | SHA256 |
| --- | --- |
| `vc_baseline.py` | `8b9f5d1f87968225c435226f45bccd1c10ca92352657bf494754e5345db9a46b` |
| `vc_internal.py` | `0d7f7e80943a78db66245f7ea690513a4ba8cab24d7c05524d1b0d24e0f94892` |
| `vc_specificity.py` | `868a9dc01fe8367e74abde0c49fd5bdff332d95ed6ffa00cdc0e56c2bebd2f5e` |

The user supplied lines 29–118 of `vc_specificity.py` and lines
148–220 of `vc_baseline.py`. Review conclusions are limited to these
excerpts; the complete functions and their callers have not yet been
audited.

## Confirmed from the excerpts

- `vc_baseline.predict()` reads source response files, aligns their
  perturbation and gene axes, computes an equal-environment source
  mean, and writes three source-only features: source heterogeneity,
  inverse minimum source cell count, and predicted response RMS. The
  shown function does not read a target response. Its target label and
  source file layout are hardcoded for the earlier HepG2 pilot.
- `scale_from_sources(a,b)` uses
  `clip(2*mean(a*b)/(mean(a*a)+mean(b*b)),0,1)` when the denominator
  is positive. This is a symmetric transfer scale estimated from two
  sources, not from the target response. The internal `freeze()`
  builds `paired=(a+b)/2` and a shared response from the source mean.
- The shown `freeze()` unpacks exactly two source arrays for each
  internal target. It therefore cannot be reused unchanged when all
  four current contexts are sources for an external HCT116 target.
  The four-source aggregation and source-only calibration must be
  specified and implemented before any external response is read.
- The retrieval function explicitly counts ties and reports chance
  references for its actual candidate count. The target-derived
  candidate mask and centering are described as evaluation-only in
  the protocol excerpt; their callers and full evaluation code still
  require inspection.

## Open issues before an external freeze

1. Archive and inspect the complete legacy source files, including
   `read_response()`, preprocessing, diagnostic code and callers.
   Verify that an external cohort cannot influence feature or task
   selection before outcomes are revealed.
2. Specify a four-source predictor, common-source predictor,
   zero-response baseline, and a calibration fit using only the
   four internal contexts. Pin code, task/feature axes, parameter
   rules and input checksums.
3. Do not promote the old 20% source-heterogeneity flag into a
   validated failure detector. In the retrospective HepG2 pilot,
   its reported AUROC for `D=MSE_pred-MSE_zero>0` was 0.208;
   candidate risk rules need an internally fixed interpretation
   before external testing.
4. Predeclare primary versus sensitivity handling of the single
   HEK293T GEM with 12 controls. The current project table shows 15
   qualifying tasks and 26 perturbation cells exposed to that GEM.
5. Check full expression/guide barcode and label membership and
   document experimental replicate provenance. A GEM capture
   batch is not by itself a biological replicate.

The internal panel was described in the excerpt as selected using all
four current environments. Its retrospective held-out-context scores
must not be presented as independent external validation. No X-Atlas
response was used for this review.

## Follow-up excerpt review

The user supplied `vc_baseline.py` lines 1–76, `vc_internal.py`
lines 1–75, and `vc_specificity.py` lines 1–28 and 118–218. These
excerpts share the SHA256 identities above. They confirm:

- The baseline's `SOURCES` are K562, RPE1 and Jurkat. Its prediction
  function and source-only risk rankings were written for the HepG2
  pilot. The internal fold's `ENVS` are the same three contexts and
  `load_source()` rejects any other environment. HepG2 is currently
  stored as `evaluation_only` and needs explicit, audited treatment
  as a *development source* in a new external-only protocol; do not
  silently relabel or overwrite its H5 file.
- The internal evaluation checks the frozen script and prediction
  digests before calling `vi.load_source(root, target, manifest)` to
  read the held-out internal response. This is a reproducible stage
  boundary in the shown code, not an access-control guarantee.
- Retrieval uses all inherited tasks for MSE, but its candidate mask
  excludes zero-norm full target responses and zero-norm responses in
  either experimental half for every repeat. The mask is shared by
  models within a representation, and its excluded tasks are written
  to `retrieval_cohort.tsv`. Thus retrieval and its `1/N` chance
  reference refer to the *actual retained candidate set*, not
  automatically all 1,264 tasks. Report the size and exclusions.
- `perturbation_centered_cosine` subtracts each matrix's mean over
  tasks, including the target matrix's mean at evaluation time. This
  is an evaluation-only diagnostic; it is not a deployable target-free
  risk feature. Raw control-relative retrieval is the primary form.
  The all-zero baseline has MSE but no defined cosine retrieval rank.
- The experimental halves are repeated splits from the same screen,
  not independent biological experiments. Reciprocal directions
  within a repeat are dependent.

No external X-Atlas response was read during the excerpt review.

## Review of the archived code

The archived preprocessing manifest uses cell and control support and exact
shared IDs from all four internal contexts to choose the inherited panel.
The normalization divides each raw gene count by the original cell UMI total
before multiplying by 10,000 and applying `log1p`. The control mean uses the
task's own proportions across GEMs. The main v01 matrix includes GEMs with
at least one control and records cells exposed to GEMs with fewer than 30.
Three HepG2 GEMs fall below 30 controls, as recorded in Step 13.

`vc_baseline.read_response()` checks the H5 role, environment, dimensions,
unique axes and finite responses. `align()` requires identical task and gene
sets and reorders them explicitly. The HepG2 baseline evaluator verifies the
code, source and prediction hashes before opening the HepG2 response H5.
`vc_internal` and `vc_specificity` evaluate K562, RPE1 and Jurkat. Their
entry points run multiple stages and can reveal internal outcomes within
the same command. Separate prediction and external evaluation entry points
are needed for the cross-study protocol.

`vc_internal` splits perturbation cells within each task and control cells
within each GEM. It reuses each control cell's assigned half across all
perturbations. It computes matched controls from each half's GEM weights
and checks reconstruction against the full response. The halves share the
same underlying screen. `vc_diagnostics` applies label permutations and
reports AUROC and average precision for positive excess MSE. Its historical
risk and split-half diagnostics concern the three internal folds. Step 29
reuses the verified loaders and metrics, fits separate paired and common
scalars, and evaluates fixed candidate scores across four nested internal
folds. See [Step 29](../workflows/step_29_four_source_calibration.md).
