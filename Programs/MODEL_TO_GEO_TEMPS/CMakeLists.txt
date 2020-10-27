cmake_minimum_required(VERSION 3.17)

project(MODEL_TO_GEO)
enable_language (CXX)

set(CMAKE_CXX_STANDARD 11)

list(APPEND SOURCE_FILES main.cpp)

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCE_FILES} )

target_link_libraries(${EXEC_NAME} -ldl)
target_link_libraries(${EXEC_NAME} -lgfortran)
target_link_libraries(${EXEC_NAME} -lssl)
target_link_libraries(${EXEC_NAME} -lcrypto)
