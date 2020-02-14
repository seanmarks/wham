### Library ###

# This directory contains the core src files
# - Builds LIB_WHAM
add_subdirectory(wham)


### Driver ###

# Add the executable
add_executable(${WHAM_DRIVER} driver/main.cpp)
target_link_libraries(${WHAM_DRIVER} PUBLIC ${LIB_WHAM})
target_include_directories(${WHAM_DRIVER} PUBLIC ${CMAKE_SOURCE_DIR}/wham)

if(APPLE)
	# Run 'codesign' on executables built on Mac
	# - If MPI is enabled, this will prevent the endless firewall prompts about allowing
	#   incoming connections (it will only appear when the program is first run)
	add_custom_command(TARGET ${WHAM_DRIVER}
	                   POST_BUILD
	                   WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
	                   COMMAND codesign --force --deep --sign - ${WHAM_DRIVER}
	                   )
endif()

# Install rules
install(TARGETS ${WHAM_DRIVER} DESTINATION bin)
#install(TARGETS ${LIB_WHAM}    DESTINATION lib)
