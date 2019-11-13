cmake_minimum_required(VERSION 2.8.3)

project(root_finder)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g -fPIC")

find_package(Boost REQUIRED)

find_package(Eigen3 REQUIRED)

include_directories(
    include
    ${Boost_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
)

add_executable(${PROJECT_NAME}_example src/root_finder_example.cpp src/root_finder.cpp)

target_link_libraries(${PROJECT_NAME}_example
  ${Boost_LIBRARIES}
)
