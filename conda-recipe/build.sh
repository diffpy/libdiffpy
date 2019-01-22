#!/bin/bash

MYNCPU=$(( (CPU_COUNT > 8) ? 8 : CPU_COUNT ))

# Apply sconscript.local customizations.
cp ${RECIPE_DIR}/sconscript.local ./

# Build and install the library.
scons -j $MYNCPU lib install prefix=$PREFIX
