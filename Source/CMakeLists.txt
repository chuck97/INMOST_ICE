add_subdirectory(Advection)
add_subdirectory(Data)
add_subdirectory(Headers)
add_subdirectory(HelloWorld)
add_subdirectory(Mesh)
add_subdirectory(Momentum)
add_subdirectory(Partitioner)
add_subdirectory(Residual)
add_subdirectory(Solver)
add_subdirectory(Thermodynamics)
#add_subdirectory(Headers)

set(
	SOURCE
	${SOURCE}
	${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  	PARENT_SCOPE
)
