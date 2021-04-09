set(TARGET_NAME mocap_s0_pin_KA_c1)
add_executable(${TARGET_NAME} ${TARGET_NAME}.cpp)
target_link_libraries(${TARGET_NAME} osimSimulation)
set_target_properties(${TARGET_NAME} PROPERTIES
    FOLDER "External_Functions")
		