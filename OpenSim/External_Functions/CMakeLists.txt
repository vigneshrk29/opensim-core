include_directories(${OpenSim_SOURCE_DIR}
                    ${OpenSim_SOURCE_DIR}/Vendors)
                    
option(BUILD_EXTERNAL_FUNCTIONS
    "Build external functions" ON)
					
if(BUILD_EXTERNAL_FUNCTIONS)
    add_subdirectory(Sh_GT_s0_pp_test)
endif()
