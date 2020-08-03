# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniMBA library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniPY_mba)
add_executable(${TEST1}.exe  "mainPYmba.f")
target_link_libraries(${TEST1}.exe view2D-${VERSION} mba2D-${VERSION} ${LAPACK_LIBS} m)
install(TARGETS ${TEST1}.exe DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniPY_stokes)
add_executable(${TEST2}.exe  "mainPYstokes.f"  "forlibaft.c")
target_link_libraries(${TEST2}.exe view2D-${VERSION} mba2D-${VERSION} 
                                    aft2D-${VERSION} fem2D-${VERSION} 
                                    lmr2D-${VERSION} lu-${UMFPACK_VER} ${LAPACK_LIBS} m)
install(TARGETS ${TEST2}.exe DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})



