cmake_minimum_required(VERSION 3.17.0)

project(TRANSPORT)
enable_language (CXX C)
set(CMAKE_CXX_COMPILER mpiicpc)
set(CMAKE_C_COMPILER mpiicc)

set(CMAKE_CXX_STANDARD 17)

set(INMOST_BUILD_ARCH "INMOST-intel-build")

set(INMOST_SRC "/data90t/geosci/spetrov/INMOST_BUILD_PVTU/INMOST")
set(INMOST_BUILD "/data90t/geosci/spetrov/INMOST_BUILD_PVTU/${INMOST_BUILD_ARCH}")
set(PETSC_DIR "/data90t/geosci/spetrov/INMOST_BUILD_PVTU/petsc-3.5.2")
set(PETSC_ARCH "linux-intel-debug")
set(JSON_DIR "/data90t/geosci/spetrov/json_lib/nlohmann_json_cmake_fetchcontent/include")

include_directories(${INMOST_SRC}/Source/Headers)
include_directories(${INMOST_BUILD})
include_directories(${JSON_DIR}/)


include_directories(${PETSC_DIR}/include)
include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)

list(APPEND HEADERS config.h INMOST_ICE_mesh.h Initialization.h Numerical.h Assembling.h Function.h Errors.h MatVec.h FCT.h Advection.h)

list(APPEND SOURCES main.cpp INMOST_ICE_mesh.cpp Initialization.cpp Numerical.cpp Assembling.cpp Function.cpp Errors.cpp MatVec.cpp FCT.cpp Advection.cpp)

set(EXEC_NAME ${CMAKE_PROJECT_NAME})

add_executable(${EXEC_NAME} ${SOURCES} ${HEADERS})

target_link_libraries(${EXEC_NAME} ${INMOST_BUILD}/libinmost.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libparmetis.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libmetis.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2clapack.a)
target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2cblas.a)
target_link_libraries(${EXEC_NAME} -ldl)
target_link_libraries(${EXEC_NAME} -lgfortran)