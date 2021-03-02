set(TEST_TARGET testShoulder)
set(TEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

add_executable(${TEST_TARGET} ${TEST_TARGET}.cpp)
target_link_libraries(${TEST_TARGET} ${Simbody_LIBRARIES} osimSimulation)
set_target_properties(${TEST_TARGET} PROPERTIES FOLDER "Examples")

file(GLOB TEST_FILES ThoracoscapularShoulderModel.osim) 

foreach(dataFile ${TEST_FILES})
 add_custom_command(
    TARGET ${TEST_TARGET}
    COMMAND ${CMAKE_COMMAND}
    ARGS -E copy
    ${dataFile}
    ${TEST_DIR})
endforeach(dataFile) 
