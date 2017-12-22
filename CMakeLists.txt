project(discretizeGeom LANGUAGES CXX CUDA)

include_directories(${CMAKE_SOURCE_DIR}/src)
find_package(MOAB REQUIRED)
include_directories(${MOAB_INCLUDE_DIR})


set(SRC
    src/discretize.h
    src/discretize.hu
    src/discretize.cpp
    src/discretize.cu
    discretize_geom.cpp)

add_executable(discretize_geom ${SRC})
set_target_properties(discretize_geom
                      PROPERTIES CUDA_SEPERABLE_COMPILATION ON)
target_link_libraries(discretize_geom MOAB)
