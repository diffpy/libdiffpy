cmake_minimum_required(VERSION 3.22)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(LibdiffpyVersion)

# A source distribution carries a generated version.hpp.  Otherwise
# resolve the version from git or from an expanded cmake/gitarchive.cfg.
set(_libdiffpy_sdist_version_header
  "${CMAKE_CURRENT_SOURCE_DIR}/src/diffpy/version.hpp")
if(EXISTS "${_libdiffpy_sdist_version_header}")
  libdiffpy_read_version_header("${_libdiffpy_sdist_version_header}")
  set(_libdiffpy_generate_version OFF)
else()
  libdiffpy_resolve_version()
  set(_libdiffpy_generate_version ON)
endif()

project(libdiffpy
  VERSION
    ${LIBDIFFPY_VERSION_MAJOR}.${LIBDIFFPY_VERSION_MINOR}.${LIBDIFFPY_VERSION_MICRO}
  DESCRIPTION
    "C++ calculators of PDF, bond valence sum and other pair quantities"
  HOMEPAGE_URL "https://github.com/diffpy/libdiffpy"
  LANGUAGES CXX)
message(STATUS "libdiffpy version ${LIBDIFFPY_VERSION_STR}")

# Options ---------------------------------------------------------------------

if(PROJECT_IS_TOP_LEVEL)
  option(BUILD_SHARED_LIBS "Build shared libraries" ON)
elseif(NOT DEFINED BUILD_SHARED_LIBS)
  # Keep the default shared without adding a parent-wide cache option.
  set(BUILD_SHARED_LIBS ON)
endif()
option(LIBDIFFPY_BUILD_OBJECT_LIB "Build an OBJECT library for embedding" OFF)
option(LIBDIFFPY_INSTALL "Install libdiffpy and its CMake package" ON)
option(LIBDIFFPY_ENABLE_OBJCRYST "Build with the ObjCryst structure adapter" ON)
option(LIBDIFFPY_LTO
  "Link the optimized shared library with link-time optimization" ON)
option(LIBDIFFPY_BUILD_TESTS "Build the Catch2 unit tests" OFF)
option(LIBDIFFPY_TEST_INSTALLED
  "Build only the unit tests and link them with the installed libdiffpy" OFF)
set(LIBDIFFPY_TESTS "" CACHE STRING
  "Comma or space separated patterns selecting the Test*.cpp unit test files")
option(LIBDIFFPY_COVERAGE "Instrument the build for coverage analysis" OFF)
option(LIBDIFFPY_PROFILE "Build with gprof profiling information" OFF)

if(DEFINED ENV{CONDA_PREFIX} AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
  file(TO_CMAKE_PATH "$ENV{CONDA_PREFIX}" _conda_prefix)
  if(WIN32)
    string(APPEND _conda_prefix "/Library")
  endif()
  # Search explicit CMake and environment prefixes before the active conda
  # environment, and the conda environment before system installations.
  if(NOT DEFINED CMAKE_FIND_USE_CMAKE_ENVIRONMENT_PATH OR
     CMAKE_FIND_USE_CMAKE_ENVIRONMENT_PATH)
    file(TO_CMAKE_PATH "$ENV{CMAKE_PREFIX_PATH}" _env_prefixes)
    list(APPEND CMAKE_PREFIX_PATH ${_env_prefixes} "${_conda_prefix}")
  endif()
  if(PROJECT_IS_TOP_LEVEL AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${_conda_prefix}" CACHE PATH
      "Install path prefix, prepended onto install directories." FORCE)
  endif()
endif()

get_property(_libdiffpy_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(PROJECT_IS_TOP_LEVEL AND NOT _libdiffpy_multi_config AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE
    PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

include(GNUInstallDirs)

# Compiler settings -----------------------------------------------------------

if(PROJECT_IS_TOP_LEVEL)
  if(MSVC)
    # /W3 rather than /Wall, which floods the log with warnings from the
    # standard library and boost headers.
    add_compile_options(/W3)
    add_link_options(/OPT:NOREF /OPT:NOICF)
  else()
    add_compile_options(-Wall $<$<CONFIG:Release>:-ffast-math>)
  endif()
endif()

if(MSVC AND (LIBDIFFPY_COVERAGE OR LIBDIFFPY_PROFILE))
  message(FATAL_ERROR
    "LIBDIFFPY_COVERAGE and LIBDIFFPY_PROFILE require gcc or clang.")
endif()

if(LIBDIFFPY_COVERAGE)
  add_compile_options(--coverage -O0 -g)
  add_link_options(--coverage)
endif()

if(LIBDIFFPY_PROFILE)
  add_compile_options(-pg)
  add_link_options(-pg)
endif()

# Unit tests of an installed library -------------------------------------------

# Do not define any other targets, they would conflict with the installed
# files.
if(LIBDIFFPY_TEST_INSTALLED)
  find_package(libdiffpy
    ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR} CONFIG REQUIRED)
  message(STATUS "Testing the installed libdiffpy in ${libdiffpy_DIR}")
  set(LIBDIFFPY_HAS_OBJCRYST ${libdiffpy_HAS_OBJCRYST})
  enable_testing()
  add_subdirectory(src/tests)
  return()
endif()

# Dependencies ----------------------------------------------------------------

# serialization of unordered_map requires boost 1.56.0
if(NOT TARGET Boost::serialization)
  find_package(Boost 1.56 CONFIG REQUIRED COMPONENTS serialization)
endif()
if(NOT TARGET GSL::gsl)
  find_package(GSL REQUIRED)
endif()

set(LIBDIFFPY_HAS_OBJCRYST OFF)
if(LIBDIFFPY_ENABLE_OBJCRYST)
  find_package(ObjCryst 2017.2.1)
  if(NOT ObjCryst_FOUND)
    message(FATAL_ERROR
      "ObjCryst library not found.  Please install libobjcryst or configure "
      "with -DLIBDIFFPY_ENABLE_OBJCRYST=OFF.")
  endif()
  set(LIBDIFFPY_HAS_OBJCRYST ON)
endif()

# Runtime data layout ---------------------------------------------------------

set(_runtime_subdir "diffpy/libdiffpy-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")
# Empty cache values mean "compute from the current layout", so the defaults
# follow later changes to the install directories, library type or version.
set(LIBDIFFPY_RUNTIME_INSTALL_DIR "" CACHE STRING
  "Install directory for libdiffpy runtime data (empty: automatic)")
if(NOT LIBDIFFPY_RUNTIME_INSTALL_DIR)
  set(LIBDIFFPY_RUNTIME_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/${_runtime_subdir}")
endif()

# Mirror the install directories under a separate build root per configuration.
# If any destination is outside the prefix, mirror the full directory hierarchy
# so that the relative path from the library to its data stays the same.
set(_build_prefix "${PROJECT_BINARY_DIR}")
if(_libdiffpy_multi_config)
  string(APPEND _build_prefix "/$<CONFIG>")
endif()
set(_layout_root "${CMAKE_INSTALL_PREFIX}")
foreach(dir IN ITEMS BINDIR LIBDIR DATADIR)
  cmake_path(ABSOLUTE_PATH CMAKE_INSTALL_${dir}
    BASE_DIRECTORY "${CMAKE_INSTALL_PREFIX}" NORMALIZE OUTPUT_VARIABLE _install_${dir})
  cmake_path(IS_PREFIX CMAKE_INSTALL_PREFIX "${_install_${dir}}" NORMALIZE _in_prefix)
  if(NOT _in_prefix)
    set(_layout_root "")
  endif()
endforeach()
foreach(dir IN ITEMS BINDIR LIBDIR DATADIR)
  if(_layout_root)
    file(RELATIVE_PATH _relative_dir "${_layout_root}" "${_install_${dir}}")
  else()
    cmake_path(GET _install_${dir} RELATIVE_PART _relative_dir)
  endif()
  set(_build_${dir} "${_build_prefix}/${_relative_dir}")
endforeach()
if(PROJECT_IS_TOP_LEVEL)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${_build_BINDIR}")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${_build_LIBDIR}")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${_build_LIBDIR}")
endif()

# runtimepath.cpp finds data relative to the module containing its code.
# Static/OBJECT builds default to an executable in bin; an embedding
# extension can override LIBDIFFPY_RUNTIME_RELPATH for its own layout.
if(WIN32 OR NOT BUILD_SHARED_LIBS OR LIBDIFFPY_BUILD_OBJECT_LIB)
  set(_module_install_dir "${_install_BINDIR}")
  set(_module_build_dir "${_build_BINDIR}")
else()
  set(_module_install_dir "${_install_LIBDIR}")
  set(_module_build_dir "${_build_LIBDIR}")
endif()
set(LIBDIFFPY_RUNTIME_RELPATH "" CACHE STRING
  "Runtime data path relative to the module containing libdiffpy (empty: automatic)")
if(NOT LIBDIFFPY_RUNTIME_RELPATH)
  cmake_path(ABSOLUTE_PATH LIBDIFFPY_RUNTIME_INSTALL_DIR
    BASE_DIRECTORY "${CMAKE_INSTALL_PREFIX}" NORMALIZE OUTPUT_VARIABLE _runtime_install_dir)
  file(RELATIVE_PATH LIBDIFFPY_RUNTIME_RELPATH
    "${_module_install_dir}" "${_runtime_install_dir}")
endif()

set(_runtime_src_dir "${PROJECT_SOURCE_DIR}/src/runtime")
file(GLOB _runtime_data_files CONFIGURE_DEPENDS "${_runtime_src_dir}/*.dat")
list(APPEND _runtime_data_files "${_runtime_src_dir}/bvparm2011sel.cif")

set(_runtime_build_dir "${_module_build_dir}/${LIBDIFFPY_RUNTIME_RELPATH}")
if(PROJECT_IS_TOP_LEVEL AND NOT IS_ABSOLUTE "${LIBDIFFPY_RUNTIME_RELPATH}")
  set(_staged_data_files)
  foreach(f IN LISTS _runtime_data_files)
    cmake_path(GET f FILENAME fname)
    set(staged "${_runtime_build_dir}/${fname}")
    add_custom_command(
      OUTPUT "${staged}"
      COMMAND "${CMAKE_COMMAND}" -E make_directory "${_runtime_build_dir}"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${f}" "${staged}"
      DEPENDS "${f}"
      VERBATIM)
    list(APPEND _staged_data_files "${staged}")
  endforeach()
  add_custom_target(libdiffpy_runtime_data ALL DEPENDS ${_staged_data_files})
elseif(PROJECT_IS_TOP_LEVEL)
  message(WARNING
    "The library and its data install on different volumes. "
    "Set DIFFPYRUNTIME to run the build-tree library.")
endif()

# Generated headers -----------------------------------------------------------

set(_generated_include_dir "${PROJECT_BINARY_DIR}/include")

if(_libdiffpy_generate_version)
  set(_version_header "${_generated_include_dir}/diffpy/version.hpp")
  libdiffpy_write_version_header(
    "${PROJECT_SOURCE_DIR}/src/diffpy/version.tpl" "${_version_header}")
  # Refresh the git version in every build.  The header is rewritten only
  # when the version changes.
  add_custom_target(libdiffpy_version
    COMMAND "${CMAKE_COMMAND}"
            "-DLIBDIFFPY_VERSION_HEADER=${_version_header}"
            "-DLIBDIFFPY_CONFIGURED_MAJORMINOR=${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
            -P "${PROJECT_SOURCE_DIR}/cmake/LibdiffpyVersion.cmake"
    BYPRODUCTS "${_version_header}"
    COMMENT "Checking libdiffpy version"
    VERBATIM)
else()
  set(_version_header "${_libdiffpy_sdist_version_header}")
endif()

if(LIBDIFFPY_HAS_OBJCRYST)
  set(DIFFPY_HAS_OBJCRYST 1)
else()
  set(DIFFPY_HAS_OBJCRYST 0)
endif()
set(_features_header "${_generated_include_dir}/diffpy/features.hpp")
configure_file("${PROJECT_SOURCE_DIR}/src/diffpy/features.tpl"
  "${_features_header}")

# Library ---------------------------------------------------------------------

file(GLOB _diffpy_sources CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/src/diffpy/*.cpp"
  "${PROJECT_SOURCE_DIR}/src/diffpy/srreal/*.cpp")
file(GLOB _diffpy_headers CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/src/diffpy/*.hpp"
  "${PROJECT_SOURCE_DIR}/src/diffpy/*.ipp")
file(GLOB _srreal_headers CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/src/diffpy/srreal/*.hpp"
  "${PROJECT_SOURCE_DIR}/src/diffpy/srreal/*.ipp")
# The version and features headers are installed separately.
list(FILTER _diffpy_headers EXCLUDE REGEX "/(version|features)\\.hpp$")
if(NOT LIBDIFFPY_HAS_OBJCRYST)
  set(_objcryst_regex "/[^/]*[Oo][Bb][Jj][Cc][Rr][Yy][Ss][Tt][^/]*$")
  list(FILTER _diffpy_sources EXCLUDE REGEX "${_objcryst_regex}")
  list(FILTER _diffpy_headers EXCLUDE REGEX "${_objcryst_regex}")
  list(FILTER _srreal_headers EXCLUDE REGEX "${_objcryst_regex}")
endif()

if(LIBDIFFPY_BUILD_OBJECT_LIB)
  add_library(diffpy OBJECT ${_diffpy_sources})
else()
  add_library(diffpy ${_diffpy_sources})
endif()
add_library(libdiffpy::diffpy ALIAS diffpy)
if(_libdiffpy_generate_version)
  add_dependencies(diffpy libdiffpy_version)
endif()

set_target_properties(diffpy PROPERTIES
  EXPORT_NAME diffpy
  CXX_EXTENSIONS OFF
  # Export the API from the DLL, see diffpy/Export.hpp.  Code that uses
  # the library imports it instead.
  DEFINE_SYMBOL BUILD_DLL
  INSTALL_NAME_DIR "${CMAKE_INSTALL_FULL_LIBDIR}")
get_target_property(_libdiffpy_library_type diffpy TYPE)
if(NOT _libdiffpy_library_type STREQUAL "SHARED_LIBRARY")
  # Static/OBJECT consumers must not import symbols from a libdiffpy DLL.
  target_compile_definitions(diffpy PUBLIC DIFFPY_STATIC)
endif()
if(MSVC)
  # Required for the Boost.Serialization code, including embedded builds.
  target_compile_options(diffpy PRIVATE /bigobj)
else()
  # Let GCC inline calls between the exported functions of the library.
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-fno-semantic-interposition
    LIBDIFFPY_HAS_NO_SEMANTIC_INTERPOSITION)
  if(LIBDIFFPY_HAS_NO_SEMANTIC_INTERPOSITION)
    target_compile_options(diffpy PRIVATE -fno-semantic-interposition)
  endif()
  # LTO inlines the short per-pair calls between source files.  Static and
  # OBJECT libraries skip it, their objects would need the same compiler.
  if(LIBDIFFPY_LTO AND _libdiffpy_library_type STREQUAL "SHARED_LIBRARY")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _ipo_supported OUTPUT _ipo_error LANGUAGES CXX)
    if(_ipo_supported)
      set_target_properties(diffpy PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION_RELEASE ON
        INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO ON)
    else()
      message(STATUS "Link-time optimization is not supported, skipping it.")
    endif()
  endif()
endif()
target_compile_features(diffpy PUBLIC cxx_std_23)
target_include_directories(diffpy PUBLIC
  "$<BUILD_INTERFACE:${_generated_include_dir}>"
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
if(WIN32)
  # Public headers use M_PI and related constants.
  target_compile_definitions(diffpy PUBLIC _USE_MATH_DEFINES)
endif()
set_property(SOURCE "${PROJECT_SOURCE_DIR}/src/diffpy/runtimepath.cpp"
  APPEND PROPERTY COMPILE_DEFINITIONS
  "DIFFPYRUNTIMERELPATH=${LIBDIFFPY_RUNTIME_RELPATH}")

if(TARGET libdiffpy_runtime_data)
  add_dependencies(diffpy libdiffpy_runtime_data)
endif()

# The dladdr call in runtimepath.cpp requires the dl library.
target_link_libraries(diffpy
  PUBLIC Boost::serialization
  PRIVATE GSL::gsl ${CMAKE_DL_LIBS})
if(LIBDIFFPY_HAS_OBJCRYST)
  target_link_libraries(diffpy PUBLIC ObjCryst::ObjCryst)
endif()
include(CalculationKernels)

# Installation ----------------------------------------------------------------

if(LIBDIFFPY_INSTALL)
  set(LIBDIFFPY_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/libdiffpy")

  install(TARGETS diffpy EXPORT libdiffpyTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT Runtime
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT Runtime
    OBJECTS DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT Development
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT Development)
  install(FILES ${_diffpy_headers} "${_version_header}" "${_features_header}"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/diffpy"
    COMPONENT Development)
  install(FILES ${_srreal_headers}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/diffpy/srreal"
    COMPONENT Development)
  install(FILES ${_runtime_data_files}
    DESTINATION "${LIBDIFFPY_RUNTIME_INSTALL_DIR}"
    COMPONENT Runtime)

  include(CMakePackageConfigHelpers)
  configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/libdiffpyConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/libdiffpyConfig.cmake"
    INSTALL_DESTINATION "${LIBDIFFPY_CMAKE_INSTALL_DIR}"
    PATH_VARS LIBDIFFPY_RUNTIME_INSTALL_DIR)
  # The runtime data directory and the ABI change with the minor version.
  write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/libdiffpyConfigVersion.cmake"
    COMPATIBILITY SameMinorVersion)
  install(EXPORT libdiffpyTargets
    NAMESPACE libdiffpy::
    DESTINATION "${LIBDIFFPY_CMAKE_INSTALL_DIR}"
    COMPONENT Development)
  set(_config_files
    "${PROJECT_BINARY_DIR}/libdiffpyConfig.cmake"
    "${PROJECT_BINARY_DIR}/libdiffpyConfigVersion.cmake")
  if(LIBDIFFPY_HAS_OBJCRYST)
    list(APPEND _config_files "${PROJECT_SOURCE_DIR}/cmake/FindObjCryst.cmake")
  endif()
  install(FILES ${_config_files}
    DESTINATION "${LIBDIFFPY_CMAKE_INSTALL_DIR}"
    COMPONENT Development)
endif()

# Other targets ---------------------------------------------------------------

if(LIBDIFFPY_BUILD_TESTS)
  enable_testing()
  add_subdirectory(src/tests)
endif()

# sdist -- source distribution bundle libdiffpy-<version>.tar.gz
add_custom_target(sdist
  COMMAND "${CMAKE_COMMAND}" "-DLIBDIFFPY_SDIST_DIR=${PROJECT_BINARY_DIR}"
          -P "${PROJECT_SOURCE_DIR}/cmake/MakeSdist.cmake"
  COMMENT "Creating source distribution"
  VERBATIM)

# zerocounters -- remove cumulative coverage counts
if(LIBDIFFPY_COVERAGE)
  add_custom_target(zerocounters
    COMMAND "${CMAKE_COMMAND}" "-DLIBDIFFPY_COVERAGE_DIR=${PROJECT_BINARY_DIR}"
            -P "${PROJECT_SOURCE_DIR}/cmake/ZeroCounters.cmake"
    VERBATIM)
endif()
