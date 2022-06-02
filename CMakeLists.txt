cmake_minimum_required(VERSION 3.1)
project(lsc LANGUAGES CXX)# level-set curves



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
set(LSC_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#list(APPEND CMAKE_MODULE_PATH ${THIRD_PARTY_DIR}/Catch2/contrib)

include(Warnings)
include(UseColors)
include(${PROJECT_NAME}Dependencies)
include(PrependCurrentPath)
include(${PROJECT_NAME}Utils)


add_definitions(-D_USE_MATH_DEFINES)
#include_directories(C:/Users/Jiangc/Documents/libigl/external/external/stb)
#target_link_libraries(C:/Users/Jiangc/Documents/libigl/build/Release/igl_stb_image.lib)


# Add your project files
file(GLOB SRCFILES *.cpp *.h)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)
target_link_libraries(${PROJECT_NAME}_bin debug OpenMeshCored.lib OpenMeshToolsd.lib)
# target_link_libraries(${PROJECT_NAME}_bin optimized OpenMeshCore.lib OpenMeshTools.lib ANN.lib)


# ############################
# # level set library.
# # igl -> lsc
# include(src/SOURCE.cmake)
# prepend_current_path(LSC_SOURCES)
# SI_copy_headers(${LSC_SOURCES})
# SI_set_source_group(${LSC_SOURCES})  
# add_library(lsc ${LSC_SOURCES})

# target_include_directories(lsc PUBLIC ${PROJECT_BINARY_DIR}/include)
# target_link_libraries(lsc PUBLIC igl::core igl::predicates igl::opengl_glfw igl::triangle Eigen3::Eigen)
# target_compile_definitions(lsc PUBLIC
#     LSC_MESH_DIR="${CMAKE_CURRENT_SOURCE_DIR}/app/meshes/")
# SET_TARGET_PROPERTIES(lsc PROPERTIES LINKER_LANGUAGE CXX)
# ###############################

# ################################


# add_executable(${PROJECT_NAME}_bin
#     app/main.cpp 
# 	app/test.cpp
#   )
# target_include_directories(${PROJECT_NAME}_bin PUBLIC ${PROJECT_BINARY_DIR}/include)
# target_link_libraries(${PROJECT_NAME}_bin PUBLIC splines igl::core igl::opengl_glfw)
# target_link_libraries(${PROJECT_NAME}_bin PUBLIC compare)



