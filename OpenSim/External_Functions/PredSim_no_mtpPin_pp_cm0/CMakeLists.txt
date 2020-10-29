set(TEST_TARGET PredSim_no_mtpPin_pp_cm0)
add_executable(${TEST_TARGET} ${TEST_TARGET}.cpp)
target_link_libraries(${TEST_TARGET} osimSimulation)
set_target_properties(${TEST_TARGET} PROPERTIES
    FOLDER "External_Functions")
		