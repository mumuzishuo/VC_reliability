# Legacy predictor review v01 (provisional, 2026-11-27)

This records a read-only review of excerpts supplied from the HPC. It is
not the frozen external validation protocol. The original Python files
remain at `$HOME/VC_reliability/` on the HPC and have not been copied into
this repository.

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

1. Inspect the source constants, `load_source()`, source manifest
   checks, and the full `evaluate()` functions. Verify which
   perturbation and expression axes can be selected without an
   external target response.
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
