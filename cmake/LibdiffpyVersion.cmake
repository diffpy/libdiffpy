# Resolve the libdiffpy version and generate diffpy/version.hpp.
#
# The version comes from, in order of preference,
#
#   1. `git describe` of the libdiffpy repository,
#   2. the release tag in an expanded cmake/gitarchive.cfg, which
#      `git archive` fills in for a GitHub source archive,
#   3. LIBDIFFPY_FALLBACK_VERSION from cmake/FallbackVersion.cmake.
#
# Release tags have the form X.Y.Z or X.Y.Z-rc.N.  Commits after a tag
# get a ".postN" suffix, where N is the number of commits since the tag.
#
# Use from CMake:
#
#   include(LibdiffpyVersion)
#   libdiffpy_resolve_version()
#   libdiffpy_write_version_header(<version.tpl> <version.hpp>)
#
# Use as a script:
#
#   cmake -P LibdiffpyVersion.cmake
#       Print the version string.
#   cmake -DLIBDIFFPY_VERSION_HEADER=<version.hpp> -P LibdiffpyVersion.cmake
#       Write version.hpp, but only when its content changes.

include_guard(GLOBAL)

set(_LIBDIFFPY_VERSION_MODULE_DIR "${CMAKE_CURRENT_LIST_DIR}")
include("${CMAKE_CURRENT_LIST_DIR}/FallbackVersion.cmake")

set(_LIBDIFFPY_EMSG_GIT_MISSING "
Cannot determine libdiffpy version.  Compile from a git repository
or use a source archive from

    https://github.com/diffpy/libdiffpy/releases
")

# Normalize a release tag X.Y.Z or X.Y.Z-rc.N to X.Y.Z or X.Y.ZrcN.
function(_libdiffpy_release_version tag outvar)
  if(NOT tag MATCHES "^([0-9]+\\.[0-9]+\\.[0-9]+)(-rc\\.([0-9]+))?$")
    message(FATAL_ERROR
      "Invalid release tag '${tag}'; expected X.Y.Z or X.Y.Z-rc.N.")
  endif()
  set(rv "${CMAKE_MATCH_1}")
  if(NOT "${CMAKE_MATCH_3}" STREQUAL "")
    string(APPEND rv "rc${CMAKE_MATCH_3}")
  endif()
  set(${outvar} "${rv}" PARENT_SCOPE)
endfunction()

# Set _gi_commit, _gi_date and, when a release tag is found, _gi_version
# and _gi_patch in the caller.  Leave them empty outside a git repository.
function(_libdiffpy_gitinfo)
  foreach(v commit date version patch)
    set(_gi_${v} "" PARENT_SCOPE)
  endforeach()
  find_program(LIBDIFFPY_GIT_EXECUTABLE NAMES git)
  mark_as_advanced(LIBDIFFPY_GIT_EXECUTABLE)
  if(NOT LIBDIFFPY_GIT_EXECUTABLE)
    return()
  endif()
  execute_process(
    COMMAND "${LIBDIFFPY_GIT_EXECUTABLE}" log -1 "--format=%H %ci"
    WORKING_DIRECTORY "${_LIBDIFFPY_VERSION_MODULE_DIR}"
    RESULT_VARIABLE rc
    OUTPUT_VARIABLE out
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(rc OR NOT out MATCHES "^([^ ]+) (.+)$")
    return()
  endif()
  set(_gi_commit "${CMAKE_MATCH_1}" PARENT_SCOPE)
  set(_gi_date "${CMAKE_MATCH_2}" PARENT_SCOPE)
  # --long makes the commit count explicit even at a tag.  Split from
  # the right so the -rc.N suffix remains part of the release tag.
  execute_process(
    COMMAND "${LIBDIFFPY_GIT_EXECUTABLE}" describe --tags --long
            "--match=[0-9]*"
    WORKING_DIRECTORY "${_LIBDIFFPY_VERSION_MODULE_DIR}"
    RESULT_VARIABLE rc
    OUTPUT_VARIABLE out
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(rc OR NOT out MATCHES "^(.+)-([0-9]+)-[^-]+$")
    return()
  endif()
  set(count "${CMAKE_MATCH_2}")
  _libdiffpy_release_version("${CMAKE_MATCH_1}" version)
  math(EXPR count "${count}")
  if(count)
    string(APPEND version ".post${count}")
  endif()
  set(_gi_version "${version}" PARENT_SCOPE)
  set(_gi_patch "${count}" PARENT_SCOPE)
endfunction()

# Resolve the version and set in the caller
#
#   LIBDIFFPY_VERSION_STR         full version, e.g. 1.5.0rc1.post3
#   LIBDIFFPY_VERSION_MAJOR
#   LIBDIFFPY_VERSION_MINOR
#   LIBDIFFPY_VERSION_MICRO
#   LIBDIFFPY_VERSION_PATCH       number of commits since the release tag
#   LIBDIFFPY_VERSION_PRERELEASE  rcN for a release candidate, else empty
#   LIBDIFFPY_VERSION_DATE        commit date
#   LIBDIFFPY_GIT_SHA             full commit hash
function(libdiffpy_resolve_version)
  _libdiffpy_gitinfo()
  if(_gi_commit STREQUAL "")
    # Not a git repository.  Require an expanded source archive.
    file(STRINGS "${_LIBDIFFPY_VERSION_MODULE_DIR}/gitarchive.cfg" galines)
    foreach(line IN LISTS galines)
      if(line MATCHES "^([A-Za-z0-9_]+) *= *(.*[^ \t])[ \t]*$")
        set(ga_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
      endif()
    endforeach()
    if(NOT DEFINED ga_commit OR ga_commit MATCHES "\\$Format:")
      message(FATAL_ERROR "${_LIBDIFFPY_EMSG_GIT_MISSING}")
    endif()
    set(_gi_commit "${ga_commit}")
    set(_gi_date "${ga_date}")
    string(REPLACE ", " ";" refs "${ga_refnames}")
    foreach(ref IN LISTS refs)
      if(ref MATCHES "^tag: ([0-9].*)$")
        _libdiffpy_release_version("${CMAKE_MATCH_1}" _gi_version)
        set(_gi_patch 0)
        break()
      endif()
    endforeach()
  endif()
  if(_gi_version STREQUAL "")
    if(NOT LIBDIFFPY_FALLBACK_VERSION MATCHES
       "^[0-9]+\\.[0-9]+\\.[0-9]+(rc[0-9]+)?\\.post0$")
      message(FATAL_ERROR
        "LIBDIFFPY_FALLBACK_VERSION must be X.Y.Z.post0 or X.Y.ZrcN.post0.")
    endif()
    # A development checkout without a release tag, or an untagged
    # archive, has no known commit count.  Keep the fallback version.
    set(_gi_version "${LIBDIFFPY_FALLBACK_VERSION}")
    set(_gi_patch 1)
  endif()
  string(REGEX REPLACE "\\.post.*$" "" vbase "${_gi_version}")
  if(NOT vbase MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+)(rc[0-9]+)?$")
    message(FATAL_ERROR "Cannot parse libdiffpy version '${_gi_version}'.")
  endif()
  set(LIBDIFFPY_VERSION_STR "${_gi_version}" PARENT_SCOPE)
  set(LIBDIFFPY_VERSION_MAJOR "${CMAKE_MATCH_1}" PARENT_SCOPE)
  set(LIBDIFFPY_VERSION_MINOR "${CMAKE_MATCH_2}" PARENT_SCOPE)
  set(LIBDIFFPY_VERSION_MICRO "${CMAKE_MATCH_3}" PARENT_SCOPE)
  set(LIBDIFFPY_VERSION_PRERELEASE "${CMAKE_MATCH_4}" PARENT_SCOPE)
  set(LIBDIFFPY_VERSION_PATCH "${_gi_patch}" PARENT_SCOPE)
  set(LIBDIFFPY_VERSION_DATE "${_gi_date}" PARENT_SCOPE)
  set(LIBDIFFPY_GIT_SHA "${_gi_commit}" PARENT_SCOPE)
endfunction()

# Write version.hpp from version.tpl using the LIBDIFFPY_VERSION_*
# variables from libdiffpy_resolve_version.  The output file is left
# untouched when its content does not change.
function(libdiffpy_write_version_header template output)
  math(EXPR numversion
    "((${LIBDIFFPY_VERSION_MAJOR} * 1000 + ${LIBDIFFPY_VERSION_MINOR})
      * 1000 + ${LIBDIFFPY_VERSION_MICRO}) * 1000
      + ${LIBDIFFPY_VERSION_PATCH}")
  # Verify that the formulas in version.tpl work as advertised.
  if(LIBDIFFPY_VERSION_MINOR GREATER_EQUAL 1000 OR
     LIBDIFFPY_VERSION_MICRO GREATER_EQUAL 1000 OR
     LIBDIFFPY_VERSION_PATCH GREATER_EQUAL 500)
    message(FATAL_ERROR "Inconsistent value of DIFFPY_VERSION = ${numversion}")
  endif()
  set(DIFFPY_VERSION "${numversion}LL")
  if(LIBDIFFPY_VERSION_PRERELEASE)
    set(DIFFPY_VERSION "(-500 + ${DIFFPY_VERSION})")
  endif()
  set(DIFFPY_VERSION_MAJOR "${LIBDIFFPY_VERSION_MAJOR}")
  set(DIFFPY_VERSION_MINOR "${LIBDIFFPY_VERSION_MINOR}")
  set(DIFFPY_VERSION_MICRO "${LIBDIFFPY_VERSION_MICRO}")
  set(DIFFPY_VERSION_PATCH "${LIBDIFFPY_VERSION_PATCH}")
  set(DIFFPY_VERSION_STR "${LIBDIFFPY_VERSION_STR}")
  set(DIFFPY_VERSION_DATE "${LIBDIFFPY_VERSION_DATE}")
  set(DIFFPY_GIT_SHA "${LIBDIFFPY_GIT_SHA}")
  configure_file("${template}" "${output}")
endfunction()

# Set the LIBDIFFPY_VERSION_* variables from an existing version.hpp,
# for example the one bundled in a source distribution.
function(libdiffpy_read_version_header header)
  file(READ "${header}" code)
  foreach(field MAJOR MINOR MICRO PATCH)
    if(NOT code MATCHES "\n#define DIFFPY_VERSION_${field} +([0-9]+)")
      message(FATAL_ERROR "Cannot read DIFFPY_VERSION_${field} from ${header}")
    endif()
    set(LIBDIFFPY_VERSION_${field} "${CMAKE_MATCH_1}" PARENT_SCOPE)
  endforeach()
  if(NOT code MATCHES "\n#define DIFFPY_VERSION_STR +\"([^\"]*)\"")
    message(FATAL_ERROR "Cannot read DIFFPY_VERSION_STR from ${header}")
  endif()
  set(LIBDIFFPY_VERSION_STR "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction()

# Script mode -----------------------------------------------------------------

if(CMAKE_SCRIPT_MODE_FILE STREQUAL CMAKE_CURRENT_LIST_FILE)
  libdiffpy_resolve_version()
  if(DEFINED LIBDIFFPY_VERSION_HEADER)
    if(NOT DEFINED LIBDIFFPY_VERSION_TEMPLATE)
      set(LIBDIFFPY_VERSION_TEMPLATE
        "${_LIBDIFFPY_VERSION_MODULE_DIR}/../src/diffpy/version.tpl")
    endif()
    # The build configured the runtime data path for this major.minor.
    set(majorminor "${LIBDIFFPY_VERSION_MAJOR}.${LIBDIFFPY_VERSION_MINOR}")
    if(DEFINED LIBDIFFPY_CONFIGURED_MAJORMINOR AND
       NOT majorminor STREQUAL LIBDIFFPY_CONFIGURED_MAJORMINOR)
      message(FATAL_ERROR
        "libdiffpy version changed from ${LIBDIFFPY_CONFIGURED_MAJORMINOR} "
        "to ${majorminor}.  Re-run cmake to update the build.")
    endif()
    libdiffpy_write_version_header(
      "${LIBDIFFPY_VERSION_TEMPLATE}" "${LIBDIFFPY_VERSION_HEADER}")
  else()
    execute_process(COMMAND "${CMAKE_COMMAND}" -E echo "${LIBDIFFPY_VERSION_STR}")
  endif()
endif()
