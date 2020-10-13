cmake_minimum_required(VERSION 2.6)
project(dummy_solver)
enable_language (CXX C)
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE DEBUG)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0  -g")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0  -g")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O0  -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0  -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0  -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O0  -g")


set(INMOST_SRC "/data90t/geosci/spetrov/INMOST_BUILD/INMOST-0.1")
set(INMOST_BUILD "/data90t/geosci/spetrov/INMOST_BUILD/INMOST-gnu-build-h1")
set(PETSC_DIR "/data90t/geosci/spetrov/INMOST_BUILD/INMOST-gnu-build-h1/petsc-3.5.2")
set(PETSC_ARCH "arch-linux2-c-debug")


include_directories(${INMOST_SRC}/Source/Headers)
include_directories(${INMOST_BUILD})


include_directories(${PETSC_DIR}/include)
include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)


set(EXEC_NAME dummy_solver)

add_executable(${EXEC_NAME} main.cpp)


target_link_libraries(${EXEC_NAME} ${INMOST_BUILD}/libinmost.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.a)

target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libparmetis.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libmetis.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2clapack.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2cblas.a)
target_link_libraries(${EXEC_NAME} -ldl)
target_link_libraries(${EXEC_NAME} -lssl)
target_link_libraries(${EXEC_NAME} -lcrypto)
