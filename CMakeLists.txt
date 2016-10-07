cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_COMPILER g++)

project(mvbb_tridecimator)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
include_directories( ${PROJECT_SOURCE_DIR}/include )
set(VCG_DIR ~/manuelb_ws/src/vcglib/)
include_directories( ${VCG_DIR} )


add_executable(${PROJECT_NAME}  examples/mvbb_tridecimator.cpp ${VCG_DIR}/wrap/ply/plylib.cpp src/mvbb_decimator.cpp)


