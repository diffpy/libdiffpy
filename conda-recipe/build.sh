#!/bin/bash

export CPATH="${PREFIX}/include:${CPATH}"
export LIBRARY_PATH="$PREFIX/lib:$LIBRARY_PATH"
if [[ `uname` == Linux ]]; then
    export LDFLAGS="${LDFLAGS} "'-Wl,-rpath,\\$$ORIGIN/.'
fi
MYNCPU=$(( (CPU_COUNT > 4) ? 4 : CPU_COUNT ))

scons -j $MYNCPU lib alltests

# Execute unit tests in a modified environment so that Anaconda
# shared libraries can be found.

if [[ `uname` == Darwin ]]; then
    DYLD_FALLBACK_LIBRARY_PATH="${PREFIX}/lib" scons test
else
    LD_LIBRARY_PATH="${PREFIX}/lib" scons test
fi

scons install prefix=$PREFIX

grep '^#define DIFFPY_VERSION_STR' < "${PREFIX}/include/diffpy/version.hpp" |
    cut -d '"' -f 2 > __conda_version__.txt
