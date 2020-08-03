# Make sure the compiler can find include files
include_directories(${INCLUDE_OUTPUT_DIRECTORY})

# Make sure the linker can find the aniLMR library once it is built
link_directories(${LIBRARY_OUTPUT_DIRECTORY})

# Add first test
set(TEST1 aniLMR_nodal${SUFFICS_EXE})
add_executable(${TEST1} "mainNodal2Metric.f")
target_link_libraries(${TEST1} lmr2D-${VERSION} mba2D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST1} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add second test
set(TEST2 aniLMR_edgeest${SUFFICS_EXE})
add_executable(${TEST2} "mainEst2Metric.f" "esterr.f")
target_link_libraries(${TEST2} lmr2D-${VERSION} mba2D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST2} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

# Add third test
set(TEST3 aniLMR_gradient${SUFFICS_EXE})
add_executable(${TEST3} "mainFunc2GradMetric.f")
target_link_libraries(${TEST3} lmr2D-${VERSION} mba2D-${VERSION} ${LAPACK_LIBS}) 
install(TARGETS ${TEST3} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})





