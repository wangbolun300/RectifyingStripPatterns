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
  sparse_interp_download_eigen()
  add_library(${PROJECT_NAME}_eigen INTERFACE)
  target_include_directories(${PROJECT_NAME}_eigen SYSTEM INTERFACE
    $<BUILD_INTERFACE:${SPARSE_EXTERNAL}/eigen>
    $<INSTALL_INTERFACE:include>
  )
  set_property(TARGET ${PROJECT_NAME}_eigen PROPERTY EXPORT_NAME Eigen3::Eigen)
  add_library(Eigen3::Eigen ALIAS ${PROJECT_NAME}_eigen)
  # Set Eigen directory environment variable (needed for EVCTCD)
  set(ENV{EIGEN3_INCLUDE_DIR} "${SPARSE_EXTERNAL}/eigen/")
endif()







  # libigl for timing
if(NOT TARGET igl::core)
  sparse_interp_download_libigl()

    # Import libigl targets
    list(APPEND CMAKE_MODULE_PATH "${SPARSE_EXTERNAL}/libigl/cmake")
    include(libigl)
  endif()

 


  
    # ANN
  set(ANN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/external/ann.zip")
  set(ANN_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/ann" )
  if(EXISTS "${ANN_FILE}" OR EXISTS "${ANN_PATH}")
    message("LSC: ANN source file exists")
  else()
    if(NOT EXISTS "${ANN_FILE}")
      message("LSC: downloading ANN")
      file(DOWNLOAD https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.zip ${ANN_FILE})
    endif()
  endif()
  if(NOT EXISTS "${ANN_PATH}")
      message("LSC: ANN unzipping")
      file(ARCHIVE_EXTRACT INPUT ${ANN_FILE} DESTINATION ${ANN_PATH})
      message("LSC: ANN is unzipped")
  else()
      message("LSC: ANN is already unzipped")
  endif()
  message("LSC: ann file: ${ANN_PATH}/ann_1.1.2/")
  if (WIN32)
    set(ANN_MAKE gmake) # Windows
  elseif (APPLE)
        set(ANN_MAKE make macosx-g++)
      else ()
        set(ANN_MAKE make)#linux
  endif()

  if(NOT WIN32)
  if (NOT EXISTS "${ANN_PATH}/ann_1.1.2/lib/libANN.a") # TODO add more lib file paths for windows
    message("LSC: ANN starting installing")
    execute_process(COMMAND ${ANN_MAKE} -C ${ANN_PATH}/ann_1.1.2 )
    message("LSC: ANN get installed")
    else()
    message("LSC: ANN is already installed")
  endif()
  endif()
   
  #################openmesh
  set(OM_FILE "${CMAKE_CURRENT_SOURCE_DIR}/external/openmesh.zip")
  set(OM_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/openmesh" )
  if(EXISTS "${OM_FILE}" OR EXISTS "${OM_PATH}")
    message("LSC: OpenMesh source file exists")
  else()
    if(NOT EXISTS "${OM_FILE}")
      message("LSC: downloading OpenMesh")
      file(DOWNLOAD https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh/-/jobs/156362/artifacts/raw/OpenMesh-9.1.0.zip ${OM_FILE})
    endif()
  endif()
  if(NOT EXISTS "${OM_PATH}")
      message("LSC: OpenMesh unzipping")
      file(ARCHIVE_EXTRACT INPUT ${OM_FILE} DESTINATION ${OM_PATH})
      message("LSC: OpenMesh is unzipped")
  else()
      message("LSC: OpenMesh is already unzipped")
  endif()
  message("LSC: OpenMesh file: ${OM_PATH}/OpenMesh-9.1.0/")
  if(NOT WIN32)
  if (NOT EXISTS "${OM_PATH}/OpenMesh-9.1.0/build/Build/lib/libOpenMeshCore.a") # TODO add more lib file paths for windows
    message("LSC: OpenMesh starting installing")
    execute_process(COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/cmake/install_openmesh.py)
    message("LSC: OpenMesh get installed")
    else()
    message("LSC: OpenMesh is already installed")
  endif()
endif()
###########
