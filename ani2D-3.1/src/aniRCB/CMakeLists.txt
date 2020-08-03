#
#  Ani2D 3-X
#  Refinening/Coarsening Bisection
#

add_library(rcb2D-${VERSION} STATIC
           "auxproc.f"  "refine.f"  "coarse.f"  "error.f"
           "../aniIO/io.f")

install(TARGETS rcb2D-${VERSION} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})




