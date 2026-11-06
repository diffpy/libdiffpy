# Find the ObjCryst++ library from libobjcryst.
#
# libobjcryst does not install a CMake package configuration, so this
# module looks up its header and library directly.
#
# Imported target:
#
#   ObjCryst::ObjCryst
#
# Result variables:
#
#   ObjCryst_FOUND
#   ObjCryst_VERSION
#   ObjCryst_INCLUDE_DIR
#   ObjCryst_LIBRARY

find_path(ObjCryst_INCLUDE_DIR NAMES ObjCryst/version.h)
# Windows builds name the import library libObjCryst.lib.
find_library(ObjCryst_LIBRARY NAMES ObjCryst libObjCryst)
mark_as_advanced(ObjCryst_INCLUDE_DIR ObjCryst_LIBRARY)

if(ObjCryst_INCLUDE_DIR)
  file(STRINGS "${ObjCryst_INCLUDE_DIR}/ObjCryst/version.h" _objcryst_vlines
       REGEX "^#define LIBOBJCRYST_VERSION_(MAJOR|MINOR|MICRO) +[0-9]+")
  set(ObjCryst_VERSION "")
  foreach(_objcryst_field MAJOR MINOR MICRO)
    foreach(_objcryst_line IN LISTS _objcryst_vlines)
      if(_objcryst_line MATCHES "_${_objcryst_field} +([0-9]+)")
        list(APPEND ObjCryst_VERSION "${CMAKE_MATCH_1}")
      endif()
    endforeach()
  endforeach()
  list(JOIN ObjCryst_VERSION "." ObjCryst_VERSION)
  unset(_objcryst_vlines)
  unset(_objcryst_line)
  unset(_objcryst_field)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ObjCryst
  REQUIRED_VARS ObjCryst_LIBRARY ObjCryst_INCLUDE_DIR
  VERSION_VAR ObjCryst_VERSION)

if(ObjCryst_FOUND AND NOT TARGET ObjCryst::ObjCryst)
  # libObjCryst.so has no SONAME.  An imported SHARED target with
  # IMPORTED_NO_SONAME links it as -lObjCryst, so that dependent libraries
  # do not record its absolute build-time path.  Windows links the import
  # library, which needs no special handling.
  if(NOT WIN32 AND ObjCryst_LIBRARY MATCHES "${CMAKE_SHARED_LIBRARY_SUFFIX}$")
    add_library(ObjCryst::ObjCryst SHARED IMPORTED)
    set_target_properties(ObjCryst::ObjCryst PROPERTIES
      IMPORTED_NO_SONAME TRUE)
  else()
    add_library(ObjCryst::ObjCryst UNKNOWN IMPORTED)
  endif()
  set_target_properties(ObjCryst::ObjCryst PROPERTIES
    IMPORTED_LOCATION "${ObjCryst_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${ObjCryst_INCLUDE_DIR}"
    # libobjcryst is built with double precision for ObjCryst's REAL type.
    # Match that ABI in all consumers of ObjCryst headers.
    INTERFACE_COMPILE_DEFINITIONS "REAL=double")
endif()
