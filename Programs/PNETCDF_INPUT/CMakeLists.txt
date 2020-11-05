cmake_minimum_required(VERSION 2.6)

project(PNETCDF_INPUT)
enable_language (CXX C)

set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_C_COMPILER mpiicc)

set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} -gcc-name=gcc-4.8)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -gxx-name=g++-4.8)

set(CMAKE_CXX_STANDARD 11)


set(PNETCDF_DIR "/data90t/geosci/spetrov/PNETCDF_BUILD/pnetcdf-intel-build")


include_directories(${PNETCDF_DIR}/include)


list(APPEND SOURCE_FILES main.cpp)

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCE_FILES} )


target_link_libraries(${EXEC_NAME} ${PNETCDF_DIR}/lib/libpnetcdf.a)
target_link_libraries(${EXEC_NAME} -ldl)
target_link_libraries(${EXEC_NAME} -lgfortran)
target_link_libraries(${EXEC_NAME} -lssl)
target_link_libraries(${EXEC_NAME} -lcrypto)
