project(discretizeGeom)

include_directories(${CMAKE_SOURCE_DIR}/src)
find_package(MOAB REQUIRED)
include_directories(${MOAB_INCLUDE_DIR})
# get_filename_component(MOAB_LIBRARY_DIRS ${MOAB_LIBRARY} PATH)
# link_directories(${MOAB_LIBRARY_DIRS})

set(SRC
    src/discretize.cpp
    src/discretize.h
    discretize_geom.cpp)

add_executable(discretize_geom ${SRC})
target_link_libraries(discretize_geom dagmc MOAB)
