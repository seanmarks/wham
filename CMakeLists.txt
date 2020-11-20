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

# Parallelization
option(MPI_ENABLED         "Enable MPI"         OFF)
option(OPENMP_ENABLED      "Enable OpenMP"      OFF)
option(OPENMP_SIMD_ENABLED "Enable OpenMP SIMD" OFF)

#set(DLIB_DIR   "${HOME}/programs/dlib" CACHE STRING "Location of dlib installation")
set(LIB_WHAM    "libwham"                  CACHE STRING "Name of the library to be created")
set(WHAM_DRIVER "${PROJECT_NAME}"       CACHE STRING "Name of the executable")

# GPTL (for benchmarking)
option(GPTL_ENABLED "Enable GPTL" OFF)
set(GPTL_DIR "${HOME}/programs/gptl" CACHE STRING "Location of GPTL installation")


######################
### Find Libraries ###
######################

### dlib ###

set(dlib_name "dlib")
find_package(dlib CONFIG REQUIRED)
if(dlib_FOUND)
	# TODO print where it was found
	message("-- Found dlib")
	#message("-- Found dlib: ${dlib_INCLUDE_DIRS}")
	#message("      dlib_INCLUDE_DIRS = ${dlib_INCLUDE_DIRS}")
	#message("      dlib_INCLUDE_DIRS = ${dlib_LIBRARIES}")
else()
	message(FATAL_ERROR "-- dlib not fond")
endif()


### OpenMP ###

if(OPENMP_ENABLED)
	message("-- OpenMP enabled")

	if(OPENMP_SIMD_ENABLED)
		# FIXME Use a more robust method for enabling/adding OpenMP SIMD
		set(OpenMP_SIMD_FLAG "-fopenmp-simd")
		if($<COMPILE_LANG_AND_ID:CXX,Intel>)
			if(NOT WIN32)
				set(OpenMP_SIMD_FLAG "-qopenmp-simd")
			else()
				set(OpenMP_SIMD_FLAG "-Qopenmp-simd")
			endif()
		endif()
		message("-- OpenMP SIMD enabled (flag: ${OpenMP_SIMD_FLAG})")
		set(OpenMP_C_FLAGS   "${OpenMP_C_FLAGS}   ${OpenMP_SIMD_FLAG}")
		set(OpenMP_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${OpenMP_SIMD_FLAG}")
	endif()

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
else()
	# Unset all OpenMP flags in case OpenMP was disabled either by the user
	# or because it was only partially detected (e.g. only for C but not C++ compiler)
	unset(OpenMP_C_FLAGS CACHE)
	unset(OpenMP_CXX_FLAGS CACHE)
endif()


### GPTL ###

if(GPTL_ENABLED)
	# TODO: make target-wise? put in config file, that also gets included in the plumed patch?
	add_definitions(-DHAVE_GPTL)

	message("-- gptl enabled")
	message("-- gptl location: ${GPTL_DIR}")

	# Look for static library
	set(LIB_GPTL "gptl")
	find_library(LIB_GPTL_STATIC
							 NAMES lib${LIB_GPTL}.a ${LIB_GPTL}.a
							 HINTS "${GPTL_DIR}/lib"
							 )
	if(LIB_GPTL_STATIC)
		message("-- Found static library for ${LIB_GPTL}: ${LIB_GPTL_STATIC}")
	endif()
endif()


### MPI ###

if(MPI_ENABLED)
	message("-- MPI enabled")
	# Add MPI
	find_package(MPI)
	if(MPI_FOUND)
		# TODO Move to config.h?
		add_definitions(-DMPI_ENABLED)
		add_definitions(-DHAVE_LIBMPI)
	else()
    message(FATAL_ERROR "MPI could not be located")
	endif()
endif()


######################
### Compiler Flags ###
######################

# Basic flags
#add_compile_options(-g -Wall)

# Optimizations
set(OPT_DEBUG_FLAGS   -O0                           CACHE STRING "Optimizations used in Debug build configurations")
set(OPT_RELEASE_FLAGS -O3 -ffast-math -march=native CACHE STRING "Optimizations used in Release build configurations")
set(IS_RELEASE_CONFIG $<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>)
set(OPT_FLAGS
	$<$<CONFIG:Debug>:${OPT_DEBUG_FLAGS}> 
	$<${IS_RELEASE_CONFIG}:${OPT_RELEASE_FLAGS}>
	CACHE STRING "Optimizations used"
	)

### Build ###

add_subdirectory(src)


### Test ###

# FIXME TODO
enable_testing()
add_subdirectory(test)
