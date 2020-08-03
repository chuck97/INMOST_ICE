# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniILU library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniILU_bcg0)
add_executable(${TEST1}.exe  "main_bcg_ilu0.f")
target_link_libraries(${TEST1}.exe ilu-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST1}.exe DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniILU_gmr0)
add_executable(${TEST2}.exe  "main_gmres_ilu0.f")
target_link_libraries(${TEST2}.exe ilu-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST2}.exe DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add third test
set(TEST3 aniILU_bcg2)
add_executable(${TEST3}.exe  "main_bcg_ilu2.f")
target_link_libraries(${TEST3}.exe ilu-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST3}.exe DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})






