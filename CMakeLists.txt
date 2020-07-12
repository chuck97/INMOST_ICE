cmake_minimum_required(VERSION 3.10)

project(INMOST_ICE)

add_subdirectory(./Source)
include_directories(./Source/Headers)

add_executable(${PROJECT_NAME}.exe ${SOURCE})

