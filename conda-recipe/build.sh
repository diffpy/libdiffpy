#!/bin/bash

export CPATH="${PREFIX}/include:${CPATH}"
export LIBRARY_PATH="${PREFIX}/lib:${LIBRARY_PATH}"

if [[ `uname` == Linux ]]; then
    export LD_LIBRARY_PATH="${PREFIX}/lib"
fi

MYNCPU=$(( (CPU_COUNT > 8) ? 8 : CPU_COUNT ))

# Apply sconscript.local customizations.
cp ${RECIPE_DIR}/sconscript.local ./

# Build and install the library.
scons -j $MYNCPU lib install prefix=$PREFIX

# Execute unit tests for the installed library.
scons -j $MYNCPU test prefix=$PREFIX test_installed=yes
