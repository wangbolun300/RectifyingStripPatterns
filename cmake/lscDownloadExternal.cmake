include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(LSC_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(LSC_EXTRA_OPTIONS "")
endif()

function(lsc_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${LSC_EXTERNAL}/${name}
        DOWNLOAD_DIR ${LSC_EXTERNAL}/.cache/${name}
        QUIET
        ${LSC_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################


# Eigen
function(lsc_download_eigen)
  lsc_download_project(eigen
  GIT_REPOSITORY           https://gitlab.com/libeigen/eigen.git
	GIT_TAG       3.3.7
  )
endfunction()



# libigl for timing and mesh processing
function(lsc_download_libigl)
   lsc_download_project(libigl
   GIT_REPOSITORY https://github.com/libigl/libigl.git
   GIT_TAG        aea868bd1fc64f71afecd2c51e51507a99d8e3e5
  )
endfunction()

# A modern string formatting library
function(lsc_download_fmt)
  lsc_download_project(fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG        6.2.0
  )
endfunction()

function(lsc_download_openmesh)
   lsc_download_project(openmesh
   GIT_REPOSITORY https://github.com/Lawrencemm/openmesh.git
   GIT_TAG        4e2e481f438747d64e62899021da9f469fd9daf8
  )
endfunction()


