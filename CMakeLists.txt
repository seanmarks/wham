cmake_minimum_required(VERSION 3.15)

# Include directory with custom CMake modules
# - These contain useful helper functions
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(setIfUnset)

project(wham)

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Default to 'Release' build type
setIfUnset(CMAKE_BUILD_TYPE Release)

# Verbose output
# - TODO: restrict to 'debug'?
set(CMAKE_VERBOSE_MAKEFILE ON)

# Put all binaries and libraries into subdirectories, instead of the main 'build' directory
# - Adapted from CMake tutorial (v3.15)
setIfUnset(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
setIfUnset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
setIfUnset(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")


#######################################
### Important variables and options ###
#######################################

#set(DLIB_DIR "${HOME}/programs/dlib" CACHE STRING "Location of dlib installation")

option(OPENMP_ENABLED "Enable OpenMP" OFF)

#set(DLIB_DIR   "${HOME}/programs/dlib" CACHE STRING "Location of dlib installation")
set(LIB_WHAM    "libwham"                  CACHE STRING "Name of the library to be created")
set(WHAM_DRIVER "${PROJECT_NAME}"       CACHE STRING "Name of the executable")

### dlib ###

set(dlib_name "dlib")
find_package(dlib CONFIG REQUIRED)
if(dlib_FOUND)
	message("-- Found dlib: ${dlib_INCLUDE_DIRS}")
	message("      dlib_INCLUDE_DIRS = ${dlib_INCLUDE_DIRS}")
	message("      dlib_INCLUDE_DIRS = ${dlib_LIBRARIES}")
else()
	message(FATAL_ERROR "-- dlib not fond")
endif()


#message("-- dlib location: ${DLIB_DIR}")

# Look for static library
#set(LIB_DLIB "dlib")
#find_library(LIB_DLIB_STATIC
#             NAMES lib${LIB_DLIB}.a ${LIB_DLIB}.a
#             HINTS "${DLIB_DIR}/lib"
#             )
#if(LIB_DLIB_STATIC)
#	message("-- Found static library for ${LIB_DLIB}: ${DLIB_LIB_STATIC}")
#endif()


### OpenMP ###

if(OPENMP_ENABLED)
	message("-- OpenMP enabled")
	# Add OpenMP
	find_package(OpenMP)
	if(NOT OpenMP_FOUND)
		if(NOT APPLE)
			message(FATAL_ERROR "OpenMP could not be located")
		else()
			# Mac is notoriously difficult: try to find Homebrew OpenMP semi-manually
			message("-- Will try to find Homebrew OpenMP")
			include(findHomebrewOpenMP)
			findHomebrewOpenMP()
		endif()
	endif()

	# Only set a linker flag if the user didn't set them manually
	if(NOT DEFINED OpenMP_LINKER_FLAGS)
			set(OpenMP_LINKER_FLAGS "${OpenMP_CXX_FLAGS}")
	endif()
	if(NOT DEFINED OpenMP_SHARED_LINKER_FLAGS)
			set(OpenMP_SHARED_LINKER_FLAGS "${OpenMP_CXX_FLAGS}")
	endif()

	# TODO OpenMP SIMD
	
else()
	# Unset all OpenMP flags in case OpenMP was disabled either by the user
	# or because it was only partially detected (e.g. only for C but not C++ compiler)
	unset(OpenMP_C_FLAGS CACHE)
	unset(OpenMP_CXX_FLAGS CACHE)
endif()


######################
### Compiler Flags ###
######################

# Basic flags
add_compile_options(-g -Wall)

# Optimizations
set(WHAM_OPT_DEBUG_FLAGS   -O0)
set(WHAM_OPT_RELEASE_FLAGS -O3 -ffast-math -march=native)
if(GNU)
	string(APPEND WHAM_OPT_RELEASE_FLAGS -ftree-vectorize)
endif()
add_compile_options("$<$<CONFIG:DEBUG>:${WHAM_OPT_DEBUG_FLAGS}>")
add_compile_options("$<$<CONFIG:RELEASE>:${WHAM_OPT_RELEASE_FLAGS}>")


### Build ###

add_subdirectory(src)


### Test ###

# FIXME TODO
enable_testing()
add_subdirectory(test)
