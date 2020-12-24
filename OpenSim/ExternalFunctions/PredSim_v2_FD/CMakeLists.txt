set(TARGET_NAME PredSim_v2_FD)
add_library(${TARGET_NAME} SHARED ${TARGET_NAME}.cpp)
target_link_libraries(${TARGET_NAME} osimSimulation)
set_target_properties(${TARGET_NAME} PROPERTIES
    FOLDER "ExternalFunctions")
		