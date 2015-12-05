#!/bin/bash

export CPATH="${PREFIX}/include:${CPATH}"
export LIBRARY_PATH="${PREFIX}/lib:${LIBRARY_PATH}"

if [[ `uname` == Linux ]]; then
    export LD_LIBRARY_PATH="${PREFIX}/lib"
fi

if [[ `uname` == Darwin ]]; then
    export DYLD_FALLBACK_LIBRARY_PATH="${PREFIX}/lib"
fi

MYNCPU=$(( (CPU_COUNT > 8) ? 8 : CPU_COUNT ))

# Build the library and unit test program.
scons -j $MYNCPU lib alltests

# Execute unit tests.
scons test

# Install the library.
scons install prefix=$PREFIX

grep '^#define DIFFPY_VERSION_STR' < "${PREFIX}/include/diffpy/version.hpp" |
    cut -d '"' -f 2 > __conda_version__.txt
