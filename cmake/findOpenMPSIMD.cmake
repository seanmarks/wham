# Adapted from CMake module: findOpenMP

function(_OPENMP_SIMD_FLAG_CANDIDATES LANG)
  if(NOT OpenMP_${LANG}_FLAG)
    unset(OpenMP_FLAG_CANDIDATES)

    set(OMP_SIMD_GNU "-fopenmp-simd")
    set(OMP_SIMD_Clang "-fopenmp-simd"
    set(OMP_SIMD_AppleClang "-fopenmp-simd")
    #set(OMP_SIMD_HP "+Oopenmp")
    if(WIN32)
      set(OMP_SIMD_Intel "-Qopenmp-simd")
    elseif(CMAKE_${LANG}_COMPILER_ID STREQUAL "Intel" AND
           "${CMAKE_${LANG}_COMPILER_VERSION}" VERSION_LESS "15.0.0.20140528")
      set(OMP_SIMD_Intel "-openmp-simd")
    else()
      set(OMP_SIMD_Intel "-qopenmp-simd")
    endif()
    #set(OMP_SIMD_MSVC "-openmp")
    #set(OMP_SIMD_PathScale "-openmp")
    #set(OMP_SIMD_NAG "-openmp")
    #set(OMP_SIMD_Absoft "-openmp")
    #set(OMP_SIMD_PGI "-mp")
    #set(OMP_SIMD_Flang "-fopenmp")
    #set(OMP_SIMD_SunPro "-xopenmp")
    #set(OMP_SIMD_XL "-qsmp=omp")
    # Cray compiler activate OpenMP with -h omp, which is enabled by default.
    #set(OMP_SIMD_Cray " " "-h omp")

    # If we know the correct flags, use those
    if(DEFINED OMP_SIMD_${CMAKE_${LANG}_COMPILER_ID})
      set(OpenMP_FLAG_CANDIDATES "${OMP_SIMD_${CMAKE_${LANG}_COMPILER_ID}}")
    # Fall back to reasonable default tries otherwise
    else()
      set(OpenMP_FLAG_CANDIDATES "-openmp-simd" "-fopenmp-simd" "-mp-simd" " ")
    endif()
    set(OpenMP_${LANG}_FLAG_CANDIDATES "${OpenMP_FLAG_CANDIDATES}" PARENT_SCOPE)
  else()
    set(OpenMP_${LANG}_FLAG_CANDIDATES "${OpenMP_${LANG}_FLAG}" PARENT_SCOPE)
  endif()
endfunction()

# sample openmp source code to test
set(OpenMP_C_CXX_TEST_SOURCE
"
#include <omp.h>
int main(void) {
#ifdef _OPENMP
  omp_get_max_threads();
  return 0;
#else
  breaks_on_purpose
#endif
}
")

# in Fortran, an implementation may provide an omp_lib.h header
# or omp_lib module, or both (OpenMP standard, section 3.1)
# Furthmore !$ is the Fortran equivalent of #ifdef _OPENMP (OpenMP standard, 2.2.2)
# Without the conditional compilation, some compilers (e.g. PGI) might compile OpenMP code
# while not actually enabling OpenMP, building code sequentially
set(OpenMP_Fortran_TEST_SOURCE
  "
      program test
      @OpenMP_Fortran_INCLUDE_LINE@
  !$  integer :: n
      n = omp_get_num_threads()
      end program test
  "
)

function(_OPENMP_WRITE_SOURCE_FILE LANG SRC_FILE_CONTENT_VAR SRC_FILE_NAME SRC_FILE_FULLPATH)
  set(WORK_DIR ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FindOpenMP)
  if("${LANG}" STREQUAL "C")
    set(SRC_FILE "${WORK_DIR}/${SRC_FILE_NAME}.c")
    file(WRITE "${SRC_FILE}" "${OpenMP_C_CXX_${SRC_FILE_CONTENT_VAR}}")
  elseif("${LANG}" STREQUAL "CXX")
    set(SRC_FILE "${WORK_DIR}/${SRC_FILE_NAME}.cpp")
    file(WRITE "${SRC_FILE}" "${OpenMP_C_CXX_${SRC_FILE_CONTENT_VAR}}")
  elseif("${LANG}" STREQUAL "Fortran")
    set(SRC_FILE "${WORK_DIR}/${SRC_FILE_NAME}.f90")
    file(WRITE "${SRC_FILE}_in" "${OpenMP_Fortran_${SRC_FILE_CONTENT_VAR}}")
    configure_file("${SRC_FILE}_in" "${SRC_FILE}" @ONLY)
  endif()
  set(${SRC_FILE_FULLPATH} "${SRC_FILE}" PARENT_SCOPE)
endfunction()
