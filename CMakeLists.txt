cmake_minimum_required(VERSION 3.1)
project(lsc)



################################################################################

# project-options



# libigl
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"              OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"                     OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_CORK              "Use Cork"                     OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_LIM               "Use LIM"                      OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"                   OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"                    OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    ON)
option(LIBIGL_WITH_PNG               "Use PNG"                      OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"                   OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"            OFF)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON)

################################################################################

### Configuration
set(SPARSE_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/Catch2/contrib)

include(Warnings)
include(UseColors)
include(${PROJECT_NAME}Dependencies)
include(PrependCurrentPath)
include(${PROJECT_NAME}Utils)

include_directories(${OM_PATH}/OpenMesh-9.0.0/src)
link_directories(${OM_PATH}/OpenMesh-9.0.0/build/Build/lib)
include_directories(${OM_PATH}/OpenMesh-9.1.0/src)
link_directories(${OM_PATH}/OpenMesh-9.1.0/build/Build/lib)
############################
# libraries
include(src/SOURCE.cmake)
prepend_current_path(LSC_SOURCES)
lsc_copy_headers(${LSC_SOURCES})
lsc_set_source_group(${LSC_SOURCES})  
add_library(lsc ${LSC_SOURCES})
#target_sources(lsc PRIVATE ${LSC_SOURCES})
target_include_directories(lsc PUBLIC ${PROJECT_BINARY_DIR}/include)
target_include_directories(lsc PUBLIC ${OM_PATH}/OpenMesh-9.1.0/src)
igl_include(core)
igl_include(opengl)
igl_include(glfw)
igl_include(imgui)
target_link_libraries(lsc PUBLIC igl::core igl::opengl igl::glfw igl::imgui Eigen3::Eigen)

if(WIN32)

target_link_libraries(lsc PUBLIC OpenMeshCore.lib OpenMeshTools.lib)
else(WIN32)
target_link_libraries(lsc PUBLIC libOpenMeshCore.a libOpenMeshTools.a)
endif(WIN32)


target_compile_definitions(lsc PUBLIC
    SI_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/app/meshes/")
##############################

add_executable(${PROJECT_NAME}_bin app/main.cpp)
target_link_libraries(${PROJECT_NAME}_bin PUBLIC lsc)

add_definitions(-D_USE_MATH_DEFINES)
target_compile_definitions(${PROJECT_NAME}_bin PUBLIC
    CHECKER_BOARD_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/")
################################

add_executable(${PROJECT_NAME}_devbin app/lsc_main.cpp)
target_link_libraries(${PROJECT_NAME}_devbin PUBLIC lsc)

add_definitions(-D_USE_MATH_DEFINES)
target_compile_definitions(${PROJECT_NAME}_devbin PUBLIC
    CHECKER_BOARD_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/")


