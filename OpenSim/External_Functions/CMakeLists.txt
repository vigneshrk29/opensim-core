include_directories(${OpenSim_SOURCE_DIR}
                    ${OpenSim_SOURCE_DIR}/Vendors)
                    
option(BUILD_EXTERNAL_FUNCTIONS
    "Build external functions" ON)
					
if(BUILD_EXTERNAL_FUNCTIONS)    
    add_subdirectory(PredSim)
    add_subdirectory(PredSim_pp)
    add_subdirectory(PredSim_SSCM)
    add_subdirectory(PredSim_SSCM_pp)
    add_subdirectory(TrackSim_1)
    add_subdirectory(TrackSim_2)
    add_subdirectory(PredSimProsthesis)
    add_subdirectory(PredSimProsthesis_pp)
    add_subdirectory(PredSim_v2)
    add_subdirectory(PredSim_v2_pp)
    add_subdirectory(PredSim_2D)
    add_subdirectory(PredSim_2D_pp)
endif()
