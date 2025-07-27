#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "QuEST::QuEST" for configuration "RelWithDebInfo"
set_property(TARGET QuEST::QuEST APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(QuEST::QuEST PROPERTIES
  IMPORTED_IMPLIB_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/QuEST.lib"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/bin/QuEST.dll"
  )

list(APPEND _cmake_import_check_targets QuEST::QuEST )
list(APPEND _cmake_import_check_files_for_QuEST::QuEST "${_IMPORT_PREFIX}/lib/QuEST.lib" "${_IMPORT_PREFIX}/bin/QuEST.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
