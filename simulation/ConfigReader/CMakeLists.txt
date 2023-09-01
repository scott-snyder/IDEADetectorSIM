file(GLOB ConfigReader_SOURCES "src/*.cc")
add_library(SCrd SHARED ${ConfigReader_SOURCES})
target_link_libraries(SCrd PUBLIC CLHEP::CLHEP)
