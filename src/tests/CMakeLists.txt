# Unit tests with the Catch2 framework.
#
# CTest runs the Test*.cpp suites together in the "alltests" program.

find_package(Catch2 3.5 CONFIG REQUIRED)
message(STATUS "Using Catch2 ${Catch2_VERSION} from ${Catch2_DIR}")

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

# alltests -- the unit test program ------------------------------------------

set(_helpers test_helpers.cpp)
if(LIBDIFFPY_HAS_OBJCRYST)
  list(APPEND _helpers objcryst_helpers.cpp)
endif()

add_executable(alltests ${_suites} ${_helpers})
target_include_directories(alltests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(alltests PRIVATE libdiffpy::diffpy Catch2::Catch2WithMain)
# test_helpers.cpp finds the testdata directory from this macro.
set_property(SOURCE test_helpers.cpp APPEND PROPERTY COMPILE_DEFINITIONS
  "DIFFPYTESTSDIRPATH=${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME alltests COMMAND alltests)
if(WIN32 AND LIBDIFFPY_TEST_INSTALLED)
  set_tests_properties(alltests PROPERTIES ENVIRONMENT_MODIFICATION
    "PATH=path_list_prepend:$<TARGET_FILE_DIR:libdiffpy::diffpy>")
endif()

# check -- build and run the unit tests.
add_custom_target(check
  COMMAND "${CMAKE_CTEST_COMMAND}" --output-on-failure -C "$<CONFIG>"
  DEPENDS alltests
  WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
  USES_TERMINAL
  VERBATIM)
