# Locate source files
file(GLOB_RECURSE LIB_WHAM_SRC "${CMAKE_CURRENT_LIST_DIR}/*.cpp" "${CMAKE_CURRENT_LIST_DIR}/*.h")

# Declare library
add_library(${LIB_WHAM} ${LIB_WHAM_SRC})

# dlib
target_include_directories(${LIB_WHAM} PUBLIC ${dlib_INCLUDE_DIRS})
target_link_libraries(${LIB_WHAM} PUBLIC ${dlib_LIBRARIES})

#target_include_directories(${LIB_WHAM} PUBLIC "${DLIB_DIR}/include")
# - Prefer static linkage
#if(LIB_DLIB_STATIC)
#	target_link_libraries(${LIB_WHAM} PUBLIC ${LIB_DLIB_STATIC})
#else()
#	# Fall back on dynamic linkage of DLIB
#	target_link_libraries(${LIB_WHAM} PUBLIC ${LIB_DLIB})
#	target_link_directories(${LIB_WHAM} PUBLIC "${DLIB_DIR}/lib")
#endif()

# OpenMP
if(OPENMP_ENABLED)
	target_link_libraries(${LIB_WHAM} PUBLIC OpenMP::OpenMP_CXX)
endif()

# Where the library should be installed
install(TARGETS ${LIB_WHAM} DESTINATION lib)
