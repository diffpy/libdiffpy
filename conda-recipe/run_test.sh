#!/bin/bash

MYNCPU=$(( (CPU_COUNT > 8) ? 8 : CPU_COUNT ))

# Remove unit tests from the build phase.
scons -C "$SRC_DIR" --clean lib alltests

# Build the unit tests program using the installed library.
scons -C "$SRC_DIR" -j $MYNCPU alltests prefix=$PREFIX test_installed=true

# Execute the unit tests.
MYALLTESTSFAST=$(ls -t ${SRC_DIR}/build/fast*/tests/alltests | head -1)
${MYALLTESTSFAST}
