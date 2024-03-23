cmake_minimum_required(VERSION 3.10)
project(FPW)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_BUILD_TYPE Release)

include_directories(/usr/include)
include_directories(${CMAKE_SOURCE_DIR}/eigen)

add_executable(FPW ../src/fpw.cpp)