cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_COMPILER g++)

project(mesh_decimator)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
include_directories( ${PROJECT_SOURCE_DIR}/include )
set(VCG_DIR ~/manuelb_ws/src/vcglib/)
include_directories( ${VCG_DIR} )

add_library(${PROJECT_NAME} ${VCG_DIR}/wrap/ply/plylib.cpp src/mvbb_decimator.cpp)

add_executable(${PROJECT_NAME}_example  examples/mvbb_tridecimator.cpp )
target_link_libraries(${PROJECT_NAME}_example ${PROJECT_NAME})


