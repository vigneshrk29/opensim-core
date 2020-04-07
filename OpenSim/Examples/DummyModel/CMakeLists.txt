add_executable(DummyModel
    DummyModel.cpp)
target_link_libraries(DummyModel ${Simbody_LIBRARIES} osimSimulation)
set_target_properties(DummyModel PROPERTIES FOLDER "Examples")
