# Unit tests with the CxxTest framework.
#
# cxxtestgen writes one source file per suite and a runner with
# the main function.  They are linked into the "alltests" program.

set(CXXTEST_ROOT "" CACHE PATH "CxxTest source checkout or install prefix")
find_path(CXXTEST_INCLUDE_DIR cxxtest/TestSuite.h HINTS "${CXXTEST_ROOT}")
find_program(CXXTEST_TESTGEN NAMES cxxtestgen cxxtestgen.py
  HINTS "${CXXTEST_ROOT}/bin")
if(NOT CXXTEST_INCLUDE_DIR OR NOT CXXTEST_TESTGEN)
  message(FATAL_ERROR
    "Cannot find CxxTest.  Install cxxtest or set CXXTEST_ROOT to a CxxTest "
    "source checkout.")
endif()

# cxxtestgen is a Python script.  Windows cannot execute it directly,
# other platforms run it with the interpreter from its #! line.
if(NOT WIN32 OR CXXTEST_TESTGEN MATCHES "\\.(exe|bat|cmd)$")
  set(_testgen "${CXXTEST_TESTGEN}")
else()
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  set(_testgen "${Python3_EXECUTABLE}" "${CXXTEST_TESTGEN}")
endif()

# Select the test suites -----------------------------------------------------

file(GLOB _suites CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/Test*.cpp")
# Windows globbing also matches lowercase helpers such as test_helpers.cpp.
list(FILTER _suites INCLUDE REGEX "/Test[^/]*\\.cpp$")
set(_objcryst_regex "/[^/]*[Oo][Bb][Jj][Cc][Rr][Yy][Ss][Tt][^/]*$")
if(NOT LIBDIFFPY_HAS_OBJCRYST)
  list(FILTER _suites EXCLUDE REGEX "${_objcryst_regex}")
endif()

if(LIBDIFFPY_TESTS)
  string(TOLOWER "${LIBDIFFPY_TESTS}" _patterns)
  string(REGEX REPLACE "[, \t]+" ";" _patterns "${_patterns}")
  list(REMOVE_ITEM _patterns "")
  set(_selected)
  foreach(suite IN LISTS _suites)
    cmake_path(GET suite FILENAME name)
    string(TOLOWER "${name}" name)
    foreach(pattern IN LISTS _patterns)
      string(FIND "${name}" "${pattern}" pos)
      if(pos GREATER_EQUAL 0)
        list(APPEND _selected "${suite}")
        break()
      endif()
    endforeach()
  endforeach()
  set(_suites ${_selected})
endif()
if(NOT _suites)
  message(FATAL_ERROR
    "Cannot find any test matching LIBDIFFPY_TESTS='${LIBDIFFPY_TESTS}'.")
endif()

# Generate the test sources --------------------------------------------------

set(_gendir "${CMAKE_CURRENT_BINARY_DIR}/cxxtest")
file(MAKE_DIRECTORY "${_gendir}")
set(_testgen_options --error-printer)

set(_generated "${_gendir}/runner.cpp")
add_custom_command(
  OUTPUT "${_gendir}/runner.cpp"
  COMMAND ${_testgen} --root ${_testgen_options} -o "${_gendir}/runner.cpp"
  DEPENDS "${CXXTEST_TESTGEN}"
  COMMENT "Generating CxxTest runner"
  VERBATIM)
foreach(suite IN LISTS _suites)
  cmake_path(GET suite STEM name)
  set(out "${_gendir}/${name}.cpp")
  add_custom_command(
    OUTPUT "${out}"
    COMMAND ${_testgen} --part ${_testgen_options} -o "${out}" "${suite}"
    DEPENDS "${suite}" "${CXXTEST_TESTGEN}"
    COMMENT "Generating CxxTest suite ${name}"
    VERBATIM)
  list(APPEND _generated "${out}")
endforeach()

# alltests -- the unit test program ------------------------------------------

set(_helpers test_helpers.cpp)
if(LIBDIFFPY_HAS_OBJCRYST)
  list(APPEND _helpers objcryst_helpers.cpp)
endif()

add_executable(alltests ${_generated} ${_helpers})
target_include_directories(alltests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(alltests SYSTEM PRIVATE "${CXXTEST_INCLUDE_DIR}")
target_link_libraries(alltests PRIVATE libdiffpy::diffpy)
# test_helpers.cpp finds the testdata directory from this macro.
set_property(SOURCE test_helpers.cpp APPEND PROPERTY COMPILE_DEFINITIONS
  "DIFFPYTESTSDIRPATH=${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME alltests COMMAND alltests)
if(WIN32 AND LIBDIFFPY_TEST_INSTALLED)
  # Windows has no rpath.  Let the test program find diffpy.dll on PATH.
  set(_path "$<TARGET_FILE_DIR:libdiffpy::diffpy>;$ENV{PATH}")
  string(REPLACE ";" "\\;" _path "${_path}")
  set_tests_properties(alltests PROPERTIES ENVIRONMENT "PATH=${_path}")
endif()

# check -- build and run the unit tests.
add_custom_target(check
  COMMAND "${CMAKE_CTEST_COMMAND}" --output-on-failure -C "$<CONFIG>"
  DEPENDS alltests
  WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
  USES_TERMINAL
  VERBATIM)
