file(GLOB DRFPIC_SOURCES "src/*.cc")
include_directories(${PROJECT_SOURCE_DIR}/DRFPIC/include)
add_library(DRFPIC SHARED ${DRFPIC_SOURCES})
target_link_libraries(DRFPIC PUBLIC CLHEP::CLHEP VecGeom::vecgeom EDM4HEP::edm4hep)
