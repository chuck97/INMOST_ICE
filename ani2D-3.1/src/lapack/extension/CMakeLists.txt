#
#  Ani2D 3.X
#  Extensin of Basic linear algebra package LAPACK-3.2
#

set(CMAKE_Fortran_FLAGS_RELEASE "")

add_library(lapack_ext-3.2 STATIC "dlamch.f")

install(TARGETS lapack_ext-3.2 DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




