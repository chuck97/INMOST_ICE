cmake_minimum_required(VERSION 2.6)

project(PARTITION)
enable_language (CXX C)

set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_C_COMPILER mpiicc)

set(CMAKE_CXX_STANDARD 20)

set(INMOST_BUILD_ARCH "INMOST-intel-build")
set(INMOST_SRC "/data90t/geosci/spetrov/INMOST_BUILD/INMOST-0.1")
set(INMOST_BUILD "/data90t/geosci/spetrov/INMOST_BUILD/${INMOST_BUILD_ARCH}")
set(PETSC_DIR "/data90t/geosci/spetrov/INMOST_BUILD/petsc-3.5.2")
set(PETSC_ARCH "linux-intel-debug")
set(PROJ4 "/data90t/geosci/spetrov/proj-4/proj-4-build")



include_directories(${INMOST_SRC}/Source/Headers)
include_directories(${INMOST_BUILD})
include_directories(${PROJ4}/include)


include_directories(${PETSC_DIR}/include)
include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)



list(APPEND SOURCE_FILES main.cpp INMOST_ICE_mesh.cpp)

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCE_FILES} )

target_link_libraries(${EXEC_NAME} ${INMOST_BUILD}/libinmost.a)

target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.a)

target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libparmetis.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libmetis.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2clapack.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2cblas.a)
target_link_libraries(${EXEC_NAME} ${PROJ4}/lib/libproj.a)
target_link_libraries(${EXEC_NAME} -ldl)
target_link_libraries(${EXEC_NAME} -lgfortran)
target_link_libraries(${EXEC_NAME} -lssl)
target_link_libraries(${EXEC_NAME} -lcrypto)
target_link_libraries(${EXEC_NAME} -lnetcdf)
target_link_libraries(${EXEC_NAME} -lhdf5)
