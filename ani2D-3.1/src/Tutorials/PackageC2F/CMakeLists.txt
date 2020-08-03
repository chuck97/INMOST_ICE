# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniC2F library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniC2F_analytic${SUFFICS_EXE})
add_executable(${TEST1} "main_analytic.c")
target_link_libraries(${TEST1} c2f2D-${VERSION} view2D-${VERSION} 
                               mba2D-${VERSION} ${LAPACK_LIBS} ${LIBRARY_G2C} m) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniC2F_nodal${SUFFICS_EXE})
add_executable(${TEST2} "main_nodal.c")
target_link_libraries(${TEST2} c2f2D-${VERSION} view2D-${VERSION} 
                               mba2D-${VERSION} ${LAPACK_LIBS} ${LIBRARY_G2C} m) 
install(TARGETS ${TEST2} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add third test
set(TEST3 aniC2F_audit${SUFFICS_EXE})
add_executable(${TEST3} "main_audit.c")
target_link_libraries(${TEST3} c2f2D-${VERSION} view2D-${VERSION} 
                               mba2D-${VERSION} ${LAPACK_LIBS} ${LIBRARY_G2C} m) 
install(TARGETS ${TEST3} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})





