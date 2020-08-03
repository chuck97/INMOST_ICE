# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniFEM library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniFEM_BC${SUFFICS_EXE})
add_executable(${TEST1} "mainBC.f")
target_link_libraries(${TEST1} fem2D-${VERSION} mba2D-${VERSION} view2D-${VERSION}) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniFEM_triangle${SUFFICS_EXE})
add_executable(${TEST2} "mainTriangle.f")
target_link_libraries(${TEST2} fem2D-${VERSION} mba2D-${VERSION}) 
install(TARGETS ${TEST2} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add third test
set(TEST3 aniFEM_template${SUFFICS_EXE})
add_executable(${TEST3} "mainTemplate.f")
target_link_libraries(${TEST3} fem2D-${VERSION} mba2D-${VERSION} view2D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST3} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})





