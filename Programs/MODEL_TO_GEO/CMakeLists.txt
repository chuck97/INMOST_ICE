cmake_minimum_required(VERSION 3.1)

project(MODEL_TO_GEO)
enable_language (CXX)
#set(CMAKE_CXX_COMPILER icpc)
#set(CMAKE_C_COMPILER icc)

#set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} -gcc-name=gcc-4.8)
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -gxx-name=g++-4.8)

set(CMAKE_CXX_STANDARD 11)

list(APPEND SOURCE_FILES main.cpp)

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCE_FILES} )
add_library(Test4 STATIC ${SOURCE_FILES})

target_link_libraries(${EXEC_NAME} -ldl)
target_link_libraries(${EXEC_NAME} -lgfortran)
target_link_libraries(${EXEC_NAME} -lssl)
target_link_libraries(${EXEC_NAME} -lcrypto)
