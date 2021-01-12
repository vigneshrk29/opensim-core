set(TEST_TARGET RError_FD)
add_library(${TEST_TARGET} SHARED ${TEST_TARGET}.cpp)
target_link_libraries(${TEST_TARGET} SimTKsimbody)
set_target_properties(${TEST_TARGET} PROPERTIES
    FOLDER "External_Functions")
		