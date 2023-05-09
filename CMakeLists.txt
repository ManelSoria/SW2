####################################################################################
# SW06_2m build system
####################################################################################

# Minimal CMake version required:
cmake_minimum_required (VERSION 2.6)

# Compiler to use (this is an OpenMPI code, so we will be using the wrapper):
SET(CMAKE_C_COMPILER /usr/bin/mpicc)
SET(CMAKE_CXX_COMPILER /usr/bin/mpicc+)

# The code uses some C99 syntax (as in for loops)
# Don't use c99 or c11, otherwise pclose will be implicit:
SET(CMAKE_C_FLAGS "-std=gnu99 -O3 ")

# Default installation prefix:
set (INSTALLDIR "$ENV{HOME}/SW/exe/")
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX ${INSTALLDIR} CACHE PATH "default install path" FORCE )
endif()

# Enable testing by default:
ENABLE_TESTING()

####################################################################################
# Options to enable / disable
####################################################################################
option (BUILD_SW
	"Build the main binary, sw" ON)
option (BUILD_POSTPROC
	"Build the postproc binary" ON)
option (BUILD_CASEGEN
	"Build the casegen binary" ON)
option (BUILD_DEMOPETITA
	"Build the demopetita binary" OFF)
option (COPY_IN_OUT_DFILES
	"During the installation process, copy all the in & out data files" ON)

# By default, the src/ directory will be used for *.h files:
include_directories("${PROJECT_BINARY_DIR}")

####################################################################################
# Required libraries
# This code looks for the required libraries before attempting to build the code
####################################################################################
# libm
find_library(LIBM m)
if (NOT LIBM)
	message (FATAL_ERROR "You need libm (-lm) to build this project")
endif()
# libopenmpi
find_library(LIBOPENMPI mpi)
if (NOT LIBOPENMPI)
	message (FATAL_ERROR "You need libopenmpi (-lmpi) to build this project.")
endif()

####################################################################################
# Common source files for all binary targets
####################################################################################s
FILE(GLOB COMMONsrc *.c)

####################################################################################
# libsppde.a library
# It is always built, otherwise the binaries won't work.
####################################################################################
include_directories ("${PROJECT_SOURCE_DIR}/sppde")
include_directories ("${PROJECT_SOURCE_DIR}")
add_subdirectory (sppde)
set (EXTRA_LIBS ${EXTRA_LIBS} sppde)

####################################################################################
# SW binary
# Only built if set to ON.
####################################################################################
if (BUILD_SW)
	include_directories ("${PROJECT_SOURCE_DIR}/main")
	add_subdirectory (main)
endif (BUILD_SW)

####################################################################################
# POSTPROC binary
# Only built if set to ON.
####################################################################################
if (BUILD_POSTPROC)
	include_directories ("${PROJECT_SOURCE_DIR}/postproc")
	add_subdirectory (postproc)
endif (BUILD_POSTPROC)

####################################################################################
# CASEGEN binary
# Only built if set to ON.
####################################################################################
if (BUILD_CASEGEN)
	include_directories ("${PROJECT_SOURCE_DIR}/casegen")
	add_subdirectory (casegen)
endif (BUILD_CASEGEN)

####################################################################################
# DEMOPETITA binary
# Only built if set to ON.
####################################################################################
if (BUILD_DEMOPETITA)
	include_directories ("${PROJECT_SOURCE_DIR}/demopetita")
	add_subdirectory (demopetita)
endif (BUILD_DEMOPETITA)
