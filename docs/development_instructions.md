# Compile and Modify diffpy Source Code

---

## I. Prerequisites

### 1) Install System Dependencies

Open a terminal and run the following commands in order:

```bash
# Update package list
sudo apt update

# Install build toolchain, CMake, Git, Python development headers, etc.
sudo apt install -y build-essential cmake git python3-dev

# Install SCons (used to build libdiffpy)
sudo apt install -y scons

# Install GSL (GNU Scientific Library) and ObjCryst (X-ray crystallography library)
sudo apt install -y libgsl-dev libobjcryst-dev
```

---

### 2) Create a Conda Virtual Environment

> If you haven't installed Miniconda, please download and install it first

```bash
# Create your environment (e.g., named your_env_name) with Python 3.11
conda create -n your_env_name python=3.11 -y

# Activate the environment (shell prompt will show (your_env_name))
conda activate your_env_name

# Install Python dependencies
conda install -y numpy

# Install additional scientific libraries
conda install -c conda-forge boost-cpp gsl pyobjcryst
```

All subsequent operations should be performed within the (your_env_name) environment!

---

## II. Get the Source Code

Assume our working directory is `~/workspace` (adjust to your actual setup):

```bash
# Switch to the working directory
cd ~/workspace

# Clone the two core repositories
git clone https://github.com/diffpy/libdiffpy.git diffpy.libdiffpy
git clone https://github.com/diffpy/diffpy.srreal.git diffpy.srreal
```

> The directory structure should be: `~/workspace/`
>
> ```
> ~/workspace/
> ├── diffpy.libdiffpy/   ← C++ core library
> └── diffpy.srreal/      ← Python interface layer
> ```

---

## III. Initial Build (Developer Mode)

> We use `--inplace` so after modifying code you only need to recompile, no `pip install` required.

---

### Step 1: Build libdiffpy (C++ Library)

```bash
# Enter the libdiffpy directory
cd ~/workspace/diffpy.libdiffpy

# Clean old builds (optional but recommended)
rm -rf build/

# Build the C++ library with SCons
scons
```

Success indicators:

- `build/fast-x86_64/libdiffpy.so` is generated
- No errors in terminal output; warnings can be ignored

Verify it was generated:

```bash
ls -l build/fast-x86_64/libdiffpy.so
```

> If you see `ObjCryst not found` or `gsl-config not found`, make sure you ran the `apt install` commands in Step 1.

---

### Step 2: Set Environment Variables

To allow srreal to find the newly built `libdiffpy.so`, set the following variables and replace the paths with your actual paths:

```bash
# Set the libdiffpy build path
export DIFFPY_LIBDIFFPY_BUILD="$HOME/workspace/diffpy.libdiffpy/build/fast-x86_64"

# Tell the linker and runtime where to find .so files
export LD_LIBRARY_PATH="$DIFFPY_LIBDIFFPY_BUILD:$LD_LIBRARY_PATH"
export LIBRARY_PATH="$DIFFPY_LIBDIFFPY_BUILD:$LIBRARY_PATH"

# Let Python import srreal source directly (developer mode)
export PYTHONPATH="$HOME/workspace/diffpy.srreal/src:$PYTHONPATH"
```

---

### Step 3: Build diffpy.srreal (Python Extension)

```bash
# Enter the srreal directory
cd ~/workspace/diffpy.srreal

# Clean old builds
rm -rf build/

# Build the Python extension (generate .so in place)
python setup.py build_ext --inplace
```

Success indicators:

- `diffpy/srreal/srreal_ext.cpython-*.so` is generated
- No `cannot find -ldiffpy` errors

---

### Step 4: Verify the Installation

```bash
# Test whether it can be imported
python -c "from diffpy.srreal.pdfcalculator import PDFCalculator; print('srreal Success！')"
```

If you see `srreal Success！`, everything is working!

---

## IV. Daily Development: Modify Code + Rebuild

### Scenario A: Only Modified diffpy.srreal Python or C++ Code

```bash
# 1. Activate your Conda environment
conda activate your_env_name

# 2. Enter the srreal directory
cd ~/workspace/diffpy.srreal

# 3. Clean and rebuild
rm -rf build/
python setup.py build_ext --inplace
```

No need to rebuild libdiffpy

---

### Scenario B: Modified diffpy.libdiffpy C++ Code

```bash
# 1. Activate your Conda environment
conda activate your_env_name

# 2. Rebuild libdiffpy
cd ~/workspace/diffpy.libdiffpy
rm -rf build/
scons

# 3. Update environment variables (ensure they point to the latest .so)
export DIFFPY_LIBDIFFPY_BUILD="$HOME/workspace/diffpy.libdiffpy/build/fast-x86_64"
export LD_LIBRARY_PATH="$DIFFPY_LIBDIFFPY_BUILD:$LD_LIBRARY_PATH"

# 4. Rebuild srreal (since it depends on libdiffpy)
cd ~/workspace/diffpy.srreal
rm -rf build/
python setup.py build_ext --inplace
```

Key point: if you change the underlying C++ library, you must rebuild the Python extension

---

## V. Common Issues and Fixes

### Issue 1: ModuleNotFoundError: No module named 'numpy'

→ Fix: run `conda install numpy` in your Conda environment

### Issue 2: cannot find -ldiffpy

→ Fix:

1. Confirm `libdiffpy.so` has been generated
2. Confirm `LD_LIBRARY_PATH` and `LIBRARY_PATH` are set correctly

### Issue 3: gsl-config: not found

→ Fix: run `sudo apt install libgsl-dev`

### Issue 4: Checking for C++ library ObjCryst... no

→ Fix: run `sudo apt install libobjcryst-dev`

If you encounter other issues, you can search relevant forums or contact qiaohai@tongji.edu.cn

---

## VI. Appendix: One-Click Environment Variable Setup

You can append the following to `~/.bashrc` so you don't have to set environment variables every time:

```bash
# diffpy development environment
export DIFFPY_WORKSPACE="$HOME/workspace"
export DIFFPY_LIBDIFFPY_BUILD="$DIFFPY_WORKSPACE/diffpy.libdiffpy/build/fast-x86_64"
export LD_LIBRARY_PATH="$DIFFPY_LIBDIFFPY_BUILD:$LD_LIBRARY_PATH"
export LIBRARY_PATH="$DIFFPY_LIBDIFFPY_BUILD:$LIBRARY_PATH"
export PYTHONPATH="$DIFFPY_WORKSPACE/diffpy.srreal/src:$PYTHONPATH"
```

Then run:

```bash
source ~/.bashrc
```
