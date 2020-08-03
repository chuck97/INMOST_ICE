# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the libraries once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 multi_DG${SUFFICS_EXE})
add_executable(${TEST1} "main.f" "forlibfem.f" "forlibaft.c")
target_link_libraries(${TEST1} aft2D-${VERSION} fem2D-${VERSION} 
                               mba2D-${VERSION} ilu-${VERSION} view2D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})


