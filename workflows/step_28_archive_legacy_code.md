# Step 28: archive exact legacy analysis scripts

The predictor/calibration review requires complete source files. This
step copies only Python scripts into a versioned archive in the GitHub
code repository. It does not touch the raw files or read an X-Atlas
response. The files are historical documentation, not a frozen method.

Run on the HPC from the repository root:

```bash
(
  set -e
  cd "$HOME/VC_reliability_repo"
  git pull --ff-only origin dev
  src="$HOME/VC_reliability"
  dst=scripts/legacy_v01
  names=(vc_preprocess.py vc_baseline.py vc_internal.py
         vc_specificity.py vc_diagnostics.py)

  for name in "${names[@]}"; do
    test ! -e "$dst/$name"
    cp -- "$src/$name" "$dst/$name"
  done

  printf '%s  %s\n' 8b9f5d1f87968225c435226f45bccd1c10ca92352657bf494754e5345db9a46b "$dst/vc_baseline.py" | sha256sum -c -
  printf '%s  %s\n' 0d7f7e80943a78db66245f7ea690513a4ba8cab24d7c05524d1b0d24e0f94892 "$dst/vc_internal.py" | sha256sum -c -
  printf '%s  %s\n' 868a9dc01fe8367e74abde0c49fd5bdff332d95ed6ffa00cdc0e56c2bebd2f5e "$dst/vc_specificity.py" | sha256sum -c -

  mkdir -p manifests
  sha256sum -- "$dst"/*.py > manifests/v01_legacy_code.sha256
  sha256sum -c manifests/v01_legacy_code.sha256
  git add -- "$dst"/*.py manifests/v01_legacy_code.sha256
  git -c core.whitespace=-blank-at-eol diff --cached --check
  git commit --only -m "archive: preserve internal pilot scripts" -- "$dst"/*.py manifests/v01_legacy_code.sha256
  git push origin dev
)
```

The three supplied SHA256 values must match exactly; stop and inspect
rather than replacing files if a check fails. The manifest also pins
`vc_preprocess.py` and `vc_diagnostics.py`. These source files are
small; this is a login-node file-copy step, not a Slurm calculation.
Do not stage the old `raw/`, `processed/` or `results/` trees.
