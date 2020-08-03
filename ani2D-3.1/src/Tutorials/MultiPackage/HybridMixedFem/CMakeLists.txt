# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the libraries once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add executable called "multi_HybridMFEM.exe", link and install
set(TEST1 multi_HybridMFEM${SUFFICS_EXE})
add_executable(${TEST1} "main.f" "forlibfem.f")
target_link_libraries(${TEST1} aft2D-${VERSION} fem2D-${VERSION} 
                               lmr2D-${VERSION} mba2D-${VERSION} 
                              view2D-${VERSION} ilu-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})


