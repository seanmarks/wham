# Testing
#
# AUTHOR: Sean M. Marks (https://github.com/seanmarks)
#
# - All tests run the program for a set of sample input, and compare the resulting output files to 
#   "gold standard" output files stored in the repo
#   - Simplistic, but easy to implement and modify
# - TODO: Convert from bash-driven testing to cmake-driven testing
#   - Would make a number of things much cleaner
#   - Could possibly leverage more cmake/ctest features

#############
### Setup ###
#############

# Find bash
find_program(BASH_PROGRAM bash)
if(NOT BASH_PROGRAM)
	message(FATAL_ERROR "Unable to find bash for testing")
endif()

# Testing driver script (absolute path)
set(TEST_DRIVER_SCRIPT ${CMAKE_CURRENT_LIST_DIR}/run_test.sh)

set(ECHO_FAILED_DIFFS "-diff" CACHE STRING "tests: whether to print diffs for tests that fail")

# Path to program (absolute)
set(PROGRAM_PATH "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${WHAM_DRIVER}")

# Directory with testing files
# - This is inside <repo>/test, *not* the build directory
set(TEST_DIR_ROOT ${CMAKE_CURRENT_LIST_DIR})


#############
### Tests ###
#############

### Builtin Tests ###

# Root directory where testing binaries will be built
set(TEST_BIN_ROOT ${CMAKE_BINARY_DIR}/test/bin)

# Name of the test driver program
# - Tests are differentiated by their directories
set(TEST_PROGRAM test_driver)

### Basic Tests ###

# TODO List from searching file names?
set(BASIC_TESTS "")
list(APPEND BASIC_TESTS "indus/bulk_water/run_wham")

foreach(TEST_NAME IN LISTS BASIC_TESTS)
	set(TEST_DIR ${TEST_DIR_ROOT}/${TEST_NAME})
	add_test(NAME ${TEST_NAME}
					 COMMAND ${BASH_PROGRAM} ${TEST_DRIVER_SCRIPT} -e ${PROGRAM_PATH} -d ${TEST_DIR} ${ECHO_FAILED_DIFFS}
					 )
endforeach()


### OpenMP tests ###

set(OPENMP_TESTS "")
list(APPEND OPENMP_TESTS "indus/bulk_water/omp_run_wham")

foreach(TEST_NAME IN LISTS OPENMP_TESTS)
	# TODO allow this to be specified on the command line
	set(OMP_NUM_THREADS 4)

	set(TEST_DIR ${TEST_DIR_ROOT}/${TEST_NAME})
	add_test(NAME ${TEST_NAME}
					 COMMAND ${BASH_PROGRAM} ${TEST_DRIVER_SCRIPT} -e ${PROGRAM_PATH} -d ${TEST_DIR} ${ECHO_FAILED_DIFFS}
									 -nt ${OMP_NUM_THREADS}
									 #--valgrind
					 )
endforeach()
