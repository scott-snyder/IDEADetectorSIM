file(GLOB MUSPC_SOURCES "src/*.cc")
include_directories(${PROJECT_SOURCE_DIR}/MUSPC/include
                    ${PROJECT_SOURCE_DIR}/SVX/include)
add_library(MUSPC SHARED ${MUSPC_SOURCES})
target_link_libraries(MUSPC PUBLIC CLHEP::CLHEP VecGeom::vecgeom)
