# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniINB library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniINB_simple)
add_executable(${TEST1}.exe  "main_simple.f")
target_link_libraries(${TEST1}.exe inb-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST1}.exe DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniINB_bratu)
add_executable(${TEST2}.exe  "main_bratu.f")
target_link_libraries(${TEST2}.exe inb-${VERSION} ilu-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST2}.exe DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})






