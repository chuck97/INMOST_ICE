cmake_minimum_required(VERSION 2.6)

project(COMPILER_TESTS)
enable_language (CXX C)

set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_C_COMPILER mpiicc)

set(CMAKE_CXX_STANDARD 20)


list(APPEND SOURCE_FILES main.cpp) 

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCE_FILES} )
