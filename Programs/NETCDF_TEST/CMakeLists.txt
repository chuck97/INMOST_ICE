cmake_minimum_required(VERSION 2.6)

project(NETCDF_TEST)
enable_language (CXX C)

set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_C_COMPILER mpiicc)

set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} -gcc-name=gcc-4.8)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -gxx-name=g++-4.8)

set(CMAKE_CXX_STANDARD 11)

list(APPEND SOURCE_FILES main.cpp netcdf_INMOST.cpp)

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCE_FILES} )


target_link_libraries(${EXEC_NAME} -ldl)
target_link_libraries(${EXEC_NAME} -lgfortran)
target_link_libraries(${EXEC_NAME} -lssl)
target_link_libraries(${EXEC_NAME} -lcrypto)
target_link_libraries(${EXEC_NAME} -lnetcdf)
target_link_libraries(${EXEC_NAME} -lhdf5)
