# Remove cumulative coverage counts (*.gcda files) under a build directory.
#
#   cmake -DLIBDIFFPY_COVERAGE_DIR=<build directory> -P ZeroCounters.cmake

if(NOT DEFINED LIBDIFFPY_COVERAGE_DIR)
  message(FATAL_ERROR "Define LIBDIFFPY_COVERAGE_DIR.")
endif()
file(GLOB_RECURSE _gcda_files "${LIBDIFFPY_COVERAGE_DIR}/*.gcda")
if(_gcda_files)
  file(REMOVE ${_gcda_files})
endif()
list(LENGTH _gcda_files _count)
message(STATUS "Removed ${_count} coverage data files.")
