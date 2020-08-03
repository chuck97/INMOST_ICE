# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniAFT library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniAFT_front1${SUFFICS_EXE})
add_executable(${TEST1}  "main_front1.f")
target_link_libraries(${TEST1} aft2D-${VERSION} view2D-${VERSION} m) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniAFT_front2${SUFFICS_EXE})
add_executable(${TEST2}  "main_front2.f")
target_link_libraries(${TEST2} aft2D-${VERSION} view2D-${VERSION} m) 
install(TARGETS ${TEST2} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add third test
set(TEST3 aniAFT_boundary_square${SUFFICS_EXE})
add_executable(${TEST3}  "main_boundary_square.f")
target_link_libraries(${TEST3} aft2D-${VERSION} view2D-${VERSION} m) 
install(TARGETS ${TEST3} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add fourth test
set(TEST4 aniAFT_boundary_wing${SUFFICS_EXE})
add_executable(${TEST4}  "main_boundary_wing.f" "crv_model.c")
target_link_libraries(${TEST4} aft2D-${VERSION} view2D-${VERSION} m) 
install(TARGETS ${TEST4} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})





