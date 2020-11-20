# AUTHOR: Sean M. Marks (https://github.com/seanmarks)

# Declare library
add_library(${LIB_WHAM} "")

target_sources(${LIB_WHAM}
	PRIVATE
		Bias.cpp
		Bins.cpp
		BiasedDistributions.cpp
		Bootstrap_F_x.cpp
		BootstrapHandler.cpp
		BootstrapSubsampler.cpp
		DataSummary.cpp
		Estimator_F_x.cpp
		Estimator_F_x_y.cpp
		FileSystem.cpp
		FreeEnergyDistribution.cpp
		FreeEnergyDistributionSet.cpp
		FreeEnergyDistribution2D.cpp
		GptlWrappers.cpp
		InputParser.cpp
		ManuallyUnbiasedDistributions.cpp
		OpenMP.cpp
		OrderParameter.cpp
		OrderParameterRegistry.cpp
		Random.cpp
		RebiasedDistributions.cpp
		Simulation.cpp
		Stacktrace.cpp
		Statistics.cpp
		StringTools.cpp
		System.cpp
		TimeSeries.cpp
		Wham.cpp
		WhamEstimator.cpp
		WhamDlibWrappers.cpp
		WhamDriver.cpp
		WhamStatistics.cpp
	)

# dlib
target_include_directories(${LIB_WHAM} PUBLIC dlib::dlib)
target_link_libraries(${LIB_WHAM} PUBLIC dlib::dlib)

# OpenMP
if(OPENMP_ENABLED)
	target_link_libraries(${LIB_WHAM} PUBLIC OpenMP::OpenMP_CXX)
endif()

# GPTL
if(GPTL_ENABLED)
	target_include_directories(${LIB_WHAM} PUBLIC "${GPTL_DIR}/include")

	# Prefer static linkage
	if(LIB_GPTL_STATIC)
		target_link_libraries(${LIB_WHAM} PUBLIC ${LIB_GPTL_STATIC})
	else()
		# Fall back on dynamic linkage of GPTL
		target_link_libraries(${LIB_WHAM} PUBLIC ${LIB_GPTL})
		target_link_directories(${LIB_WHAM} PUBLIC "${GPTL_DIR}/lib")
	endif()
endif()

# Where the library should be installed
install(TARGETS ${LIB_WHAM} DESTINATION lib)
