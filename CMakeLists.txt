cmake_minimum_required(VERSION 2.8)
project(Harvest)

SET(CMAKE_BUILD_TYPE "Release")

file(GLOB SOURCES src/*.cpp)

include_directories(src include)

add_executable(Harvest ${SOURCES} ${INCLUDES})
target_link_libraries(Harvest m)
