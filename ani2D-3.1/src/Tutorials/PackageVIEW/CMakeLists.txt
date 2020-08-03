# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniVIEW library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniVIEW_mesh${SUFFICS_EXE})
add_executable(${TEST1} "main.f")
target_link_libraries(${TEST1} view2D-${VERSION}) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniVIEW_matrix${SUFFICS_EXE})
add_executable(${TEST2} "main_matrix.f")
target_link_libraries(${TEST2} view2D-${VERSION}) 
install(TARGETS ${TEST2} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})



