set(TARGET_NAME RError_Euler_FD)
add_library(${TARGET_NAME} SHARED ${TARGET_NAME}.cpp)
target_link_libraries(${TARGET_NAME} SimTKsimbody)
set_target_properties(${TARGET_NAME} PROPERTIES
    FOLDER "External_Functions")
		