# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the libraries once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 multi_IsoAdapHessian${SUFFICS_EXE})
add_executable(${TEST1} "main_hes.f" "forlibfem_hes.f" "forlibaft.c" "forlibmba.f")
target_link_libraries(${TEST1} aft2D-${VERSION} fem2D-${VERSION} 
                               lmr2D-${VERSION} mba2D-${VERSION} 
                                 ilu-${VERSION} view2D-${VERSION} lu-${UMFPACK_VER} ${LAPACK_LIBS}) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 multi_IsoAdapPosteriori${SUFFICS_EXE})
add_executable(${TEST2} "main_est.f" "forlibfem_est.f" "forlibaft.c" "forlibmba.f")
target_link_libraries(${TEST2} aft2D-${VERSION} fem2D-${VERSION} 
                               lmr2D-${VERSION} mba2D-${VERSION} 
                                 ilu-${VERSION} view2D-${VERSION} lu-${UMFPACK_VER} ${LAPACK_LIBS}) 
install(TARGETS ${TEST2} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

