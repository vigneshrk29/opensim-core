include_directories(${OpenSim_SOURCE_DIR}
                    ${OpenSim_SOURCE_DIR}/Vendors)
                    
option(BUILD_EXTERNAL_FUNCTIONS
    "Build external functions" ON)
					
if(BUILD_EXTERNAL_FUNCTIONS)
    add_subdirectory(PredSim_v2)
    add_subdirectory(video_subject1_A_tread)
    add_subdirectory(video_subject1_tread_pp)    
    add_subdirectory(ShoulderModel)  
    add_subdirectory(ShoulderModel_consA)
    add_subdirectory(ShoulderModel_consB)
endif()
