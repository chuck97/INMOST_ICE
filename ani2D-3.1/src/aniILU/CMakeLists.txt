#
#  Ani2D 3-X
#  Incomplete LU
#

add_library(ilu-${VERSION} STATIC
            "bcg.f" "cg.f" "gmres.f" "matvecCSR.f" "dsort.f" 
            "ilu0.f" "iluoo.f"
            "../aniIO/io.f")

install(TARGETS ilu-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




