find_path(RDKit_INCLUDE_DIR GraphMol/ROMol.h PATH_SUFFIXES rdkit)
set(_rdkit_components SmilesParse GraphMol RDGeneral)
foreach(component IN LISTS _rdkit_components)
  find_library(RDKit_${component}_LIBRARY NAMES RDKit${component})
endforeach()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(RDKit REQUIRED_VARS RDKit_INCLUDE_DIR
  RDKit_SmilesParse_LIBRARY RDKit_GraphMol_LIBRARY RDKit_RDGeneral_LIBRARY)
if(RDKit_FOUND AND NOT TARGET RDKit::RDKit)
  add_library(RDKit::RDKit INTERFACE IMPORTED)
  set_target_properties(RDKit::RDKit PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${RDKit_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${RDKit_SmilesParse_LIBRARY};${RDKit_GraphMol_LIBRARY};${RDKit_RDGeneral_LIBRARY}")
endif()
