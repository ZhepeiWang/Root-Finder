cmake_minimum_required(VERSION 2.8.3)

project(root_finder)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -fPIC")
# set(CMAKE_CXX_FLAGS "-std=c++11 -pg")

find_package(Eigen3 REQUIRED)

include_directories(
    include
    ${EIGEN3_INCLUDE_DIRS}
)

add_executable(${PROJECT_NAME}_example example/root_finder_example.cpp example/JenkinsTraub/jt.cpp)