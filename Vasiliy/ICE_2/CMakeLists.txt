cmake_minimum_required(VERSION 3.5.2)

project(solver)
enable_language (CXX C)
set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_C_COMPILER mpiicc)


set(CMAKE_CXX_STANDARD 11)


set(INMOST_SRC "/data90t/geosci/spetrov/INMOST_ICE/INMOST_ROOT/INMOST-0.1")
set(INMOST_BUILD "/data90t/geosci/spetrov/INMOST_ICE/INMOST_ROOT/INMOST-intel-build")
set(PETSC_DIR "/data90t/geosci/spetrov/INMOST_ICE/INMOST_ROOT/petsc-3.5.2")
set(PETSC_ARCH "linux-intel-debug")


include_directories(${INMOST_SRC}/Source/Headers)
include_directories(${INMOST_BUILD})


include_directories(${PETSC_DIR}/include)
include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)

set(SOURCES main.cpp)

set(EXEC_NAME ${PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCES} )


target_link_libraries(${EXEC_NAME} ${INMOST_BUILD}/libinmost.a)

target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so)

target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libparmetis.so)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libmetis.so)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2clapack.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2cblas.a)
target_link_libraries(${EXEC_NAME} -ldl)
target_link_libraries(${EXEC_NAME} -lgfortran)
target_link_libraries(${EXEC_NAME} -lX11)
target_link_libraries(${EXEC_NAME} -lssl)
target_link_libraries(${EXEC_NAME} -lcrypto)
