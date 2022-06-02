# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.

# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(${PROJECT_NAME}DownloadExternal)

################################################################################
# Required libraries
################################################################################

# Eigen
if(NOT TARGET Eigen3::Eigen)
  lsc_download_eigen()
  add_library(${PROJECT_NAME}_eigen INTERFACE)
  target_include_directories(${PROJECT_NAME}_eigen SYSTEM INTERFACE
    $<BUILD_INTERFACE:${LCS_EXTERNAL}/eigen>
    $<INSTALL_INTERFACE:include>
  )
  set_property(TARGET ${PROJECT_NAME}_eigen PROPERTY EXPORT_NAME Eigen3::Eigen)
  add_library(Eigen3::Eigen ALIAS ${PROJECT_NAME}_eigen)
  # Set Eigen directory environment variable (needed for EVCTCD)
  set(ENV{EIGEN3_INCLUDE_DIR} "${LSC_EXTERNAL}/eigen/")
endif()


# function(lsc_download_openmesh)
#   execute_process(
#   COMMAND <git> [clone --recursive https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh.git]
#   WORKING_DIRECTORY <${LSC_EXTERNAL}/external>
#   )
  
# endfunction()

if(NOT TARGET OpenMeshCore)
lsc_download_openmesh()
    # Import libigl targets
    list(APPEND CMAKE_MODULE_PATH "${LSC_EXTERNAL}/openmesh/cmake")
    include(OpenMeshPackage)
  endif()



  # libigl for timing
if(NOT TARGET igl::core)
   lsc_download_libigl()
    # Import libigl targets
    list(APPEND CMAKE_MODULE_PATH "${LSC_EXTERNAL}/libigl/cmake")
    include(libigl)
  endif()

set(ANN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/external/ann.tar.gz")
set(ANN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/ann" )
if(EXISTS "${ANN_FILE}" OR EXISTS "${ANN_PATH}")
  message("LSC: ANN source file exists")
else()
  if(NOT EXISTS "${ANN_FILE}")
    message("LSC: downloading ANN")
    file(DOWNLOAD <https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz> ${ANN_FILE})
  endif()
endif()
if(NOT EXISTS "${ANN_PATH}")
    message("LSC: ANN unzipping")
    # extract_file(${CMAKE_CURRENT_SOURCE_DIR}/external/ann.zip ${CMAKE_CURRENT_SOURCE_DIR}/external/ann)
    # execute_process(
    #   COMMAND ${CMAKE_COMMAND} -E tar xzf "${ANN_FILE}"
    #   -- "${ANN_PATH}")
    file(ARCHIVE_EXTRACT INPUT ${ANN_FILE} DESTINATION ${ANN_PATH})
    message("LSC: ANN is unzipped")
else()
    message("LSC: ANN is already unzipped")
endif()
# execute_process(COMMAND ls "${CMAKE_CURRENT_SOURCE_DIR}/external/ann.zip" RESULT_VARIABLE result OUTPUT_QUIET ERROR_QUIET)


