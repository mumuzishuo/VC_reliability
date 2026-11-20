# Step 1. Environment and raw-file audit

Run these commands on the HPC login node before any expression analysis.

## 1. Enter the project directory

```bash
cd "$HOME/VC_reliability"
pwd
```

## 2. Load the Python 3.9 environment

Use the available module if needed.

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

## 4. List the four raw files

```bash
find "$HOME/VC_reliability/raw" -maxdepth 1 -type f -name "*.h5ad" -printf "%f\t%s bytes\n" | sort
```

## 5. Check storage and file checksums

```bash
df -h "$HOME/VC_reliability/raw"
sha256sum "$HOME"/VC_reliability/raw/*.h5ad
```

The checksum command can take time for multi-gigabyte files. Save its output as a manifest after the command finishes.

## Expected output

Record:

- project directory
- Python version
- package availability
- exact raw filenames
- file sizes
- available storage
- SHA256 values

Do not run model training or download external datasets at this stage.
