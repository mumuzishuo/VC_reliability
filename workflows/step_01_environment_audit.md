# Step 1. Environment and raw-file audit

The raw data directory and the Git repository are kept separate.

- Raw data directory: `$HOME/VC_reliability/raw`
- Code repository directory: `$HOME/VC_reliability_repo`

This prevents Git operations from touching the large raw files.

## 1. Clone or update the code repository

Run this on the HPC login node.

```bash
CODE_DIR="$HOME/VC_reliability_repo"

if [ -d "$CODE_DIR/.git" ]; then
    cd "$CODE_DIR"
    git fetch origin dev
    git checkout dev
    git pull --ff-only origin dev
else
    git clone --branch dev         https://github.com/mumuzishuo/VC_reliability.git         "$CODE_DIR"
    cd "$CODE_DIR"
fi

git remote -v
git branch --show-current
git log -1 --oneline
```

The raw data remains in `$HOME/VC_reliability/raw`.

## 2. Load the Python 3.9 environment

```bash
module load anaconda3/2022.10
python3 --version
```

The expected Python major and minor version is 3.9.

## 3. Check required packages without opening the count matrix

```bash
python3 - <<'PY'
import importlib.util
import sys

print("python", sys.version)
for name in ["anndata", "h5py", "pandas", "numpy"]:
    available = importlib.util.find_spec(name) is not None
    print(name, available)
PY
```

The first audit uses h5py and does not require anndata.

## 4. List the four raw files

```bash
find "$HOME/VC_reliability/raw"     -maxdepth 1     -type f     -name "*.h5ad"     -printf "%f\t%s bytes\n" | sort
```

## 5. Check storage

```bash
df -h "$HOME/VC_reliability/raw"
```

SHA256 checksums will be computed after the paths are confirmed.

## Expected output

Record:

- code repository branch and commit
- Python version
- package availability
- exact raw filenames
- file sizes
- available storage

Raw data must not be moved, renamed or overwritten during this stage.
