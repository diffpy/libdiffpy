# Create the source distribution libdiffpy-<version>.tar.gz from git HEAD.
#
#   cmake [-DLIBDIFFPY_SDIST_DIR=<output directory>] -P cmake/MakeSdist.cmake
#
# The archive carries the generated src/diffpy/version.hpp in place of
# version.tpl, so that it builds without git.  See .gitattributes for the
# files excluded from the archive.  Requires git 2.30 or later.

cmake_minimum_required(VERSION 3.22)

include("${CMAKE_CURRENT_LIST_DIR}/LibdiffpyVersion.cmake")
get_filename_component(_srcdir "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
if(NOT DEFINED LIBDIFFPY_SDIST_DIR)
  set(LIBDIFFPY_SDIST_DIR "${_srcdir}")
endif()

find_program(LIBDIFFPY_GIT_EXECUTABLE NAMES git)
if(LIBDIFFPY_GIT_EXECUTABLE)
  execute_process(
    COMMAND "${LIBDIFFPY_GIT_EXECUTABLE}" rev-parse --verify --quiet HEAD
    WORKING_DIRECTORY "${_srcdir}"
    RESULT_VARIABLE _rc
    OUTPUT_QUIET ERROR_QUIET)
endif()
if(NOT LIBDIFFPY_GIT_EXECUTABLE OR _rc)
  message(FATAL_ERROR "sdist must be built from a git repository.")
endif()

libdiffpy_resolve_version()
set(_name "libdiffpy-${LIBDIFFPY_VERSION_STR}")
set(_tmpdir "${LIBDIFFPY_SDIST_DIR}/${_name}.sdist-tmp")
set(_output "${LIBDIFFPY_SDIST_DIR}/${_name}.tar.gz")
file(REMOVE_RECURSE "${_tmpdir}")
file(MAKE_DIRECTORY "${_tmpdir}")

# Generate version.hpp from the committed template.
execute_process(
  COMMAND "${LIBDIFFPY_GIT_EXECUTABLE}" show HEAD:src/diffpy/version.tpl
  WORKING_DIRECTORY "${_srcdir}"
  OUTPUT_FILE "${_tmpdir}/version.tpl"
  COMMAND_ERROR_IS_FATAL ANY)
libdiffpy_write_version_header("${_tmpdir}/version.tpl" "${_tmpdir}/version.hpp")

# --add-file takes the archive path from the last preceding --prefix.
execute_process(
  COMMAND "${LIBDIFFPY_GIT_EXECUTABLE}" -c tar.umask=0022 archive
          --worktree-attributes --format=tar.gz
          "--prefix=${_name}/src/diffpy/" "--add-file=${_tmpdir}/version.hpp"
          "--prefix=${_name}/" "--output=${_output}"
          HEAD -- . ":(exclude)src/diffpy/version.tpl"
  WORKING_DIRECTORY "${_srcdir}"
  COMMAND_ERROR_IS_FATAL ANY)
file(REMOVE_RECURSE "${_tmpdir}")
message(STATUS "Created ${_output}")
