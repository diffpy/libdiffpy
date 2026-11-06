# Unit tests with the Catch2 framework.
#
# The Test*.cpp files share the "alltests" program. CTest runs each tagged
# suite separately, preserving filtering without executing test discovery.

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
set_target_properties(alltests PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
  target_compile_options(alltests PRIVATE /bigobj)
endif()
target_include_directories(alltests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(alltests PRIVATE libdiffpy::diffpy Catch2::Catch2WithMain)
# test_helpers.cpp finds the testdata directory from this macro.
set_property(SOURCE test_helpers.cpp APPEND PROPERTY COMPILE_DEFINITIONS
  "DIFFPYTESTSDIRPATH=${CMAKE_CURRENT_SOURCE_DIR}")

if(NOT PROJECT_IS_TOP_LEVEL AND NOT LIBDIFFPY_TEST_INSTALLED AND
   NOT IS_ABSOLUTE "${LIBDIFFPY_RUNTIME_RELPATH}")
  # Stage data for the optional embedded tests without changing the parent's
  # output layout. Static/OBJECT code lives in alltests; shared code in diffpy.
  if(_libdiffpy_library_type STREQUAL "SHARED_LIBRARY")
    set(_test_module diffpy)
  else()
    set(_test_module alltests)
  endif()
  set(_test_runtime_dir
    "$<TARGET_FILE_DIR:${_test_module}>/${LIBDIFFPY_RUNTIME_RELPATH}")
  add_custom_command(TARGET alltests POST_BUILD
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${_test_runtime_dir}"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            ${_runtime_data_files} "${_test_runtime_dir}"
    VERBATIM)
endif()

# Some files contain several suites, e.g. the structure and bond generator
# suites in TestPeriodicStructureAdapter.cpp. Register every tag in each
# selected file so none of these cases are lost.
set(_suite_tags)
foreach(suite IN LISTS _suites)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${suite}")
  file(STRINGS "${suite}" _tag_lines REGEX "\"\\[[^]]+\\]\"")
  string(REGEX MATCHALL "\"\\[[^]]+\\]\"" _tags "${_tag_lines}")
  list(APPEND _suite_tags ${_tags})
endforeach()
list(REMOVE_DUPLICATES _suite_tags)
foreach(tag IN LISTS _suite_tags)
  string(REGEX REPLACE "^\"\\[|\\]\"$" "" name "${tag}")
  add_test(NAME "${name}" COMMAND alltests "[${name}]")
  if(WIN32 AND LIBDIFFPY_TEST_INSTALLED)
    set_tests_properties("${name}" PROPERTIES ENVIRONMENT_MODIFICATION
      "PATH=path_list_prepend:$<TARGET_FILE_DIR:libdiffpy::diffpy>")
  endif()
endforeach()

# check -- build and run the unit tests.
add_custom_target(check
  COMMAND "${CMAKE_CTEST_COMMAND}" --output-on-failure -C "$<CONFIG>"
  DEPENDS alltests
  WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
  USES_TERMINAL
  VERBATIM)
